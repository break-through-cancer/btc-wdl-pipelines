/*
 * --------------------------------------------
 * Defaults / params
 * --------------------------------------------
 */

def NO_NORMAL_BAM_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAM"
def NO_NORMAL_BAI_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAI"
def NO_ALLELES_VCF_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
def NO_ALLELES_TBI_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_TBI"

new File("${workflow.projectDir}/assets").mkdirs()
new File(NO_NORMAL_BAM_PATH).createNewFile()
new File(NO_NORMAL_BAI_PATH).createNewFile()
new File(NO_ALLELES_VCF_PATH).createNewFile()
new File(NO_ALLELES_TBI_PATH).createNewFile()

if( !params.containsKey('tumor_sample') )
  params.tumor_sample = null

if( !params.containsKey('m2_extra_args') )
  params.m2_extra_args = ''

if( !params.containsKey('normal_reads') )
  params.normal_reads = null
if( !params.containsKey('normal_reads_index') )
  params.normal_reads_index = null
if( !params.containsKey('force_call_file') )
  params.force_call_file = null
if( !params.containsKey('force_call_file_index') )
  params.force_call_file_index = null

process split_intervals {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path ref_fasta
    path ref_fai
    path ref_dict
    path intervals
    val scatter_count

  output:
    path "scattered/*.intervals", emit: interval_shards

  script:
  """
  set -euo pipefail

  echo "=== split_intervals: START ==="
  echo "PWD=\$(pwd)"
  echo "ref_fasta=$ref_fasta"
  echo "intervals=$intervals"
  echo "scatter_count=$scatter_count"
  ls -lah
  mkdir -p scattered

  echo "--- Running BedToIntervalList ---"
  gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \
    -I "$intervals" \
    -SD "$ref_dict" \
    -O regions.interval_list
  echo "--- BedToIntervalList done ---"

  echo "--- Running SplitIntervals ---"
  gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \
    -R "$ref_fasta" \
    -L regions.interval_list \
    --scatter "$scatter_count" \
    -O scattered
  echo "--- SplitIntervals done ---"

  echo "Shards produced:"
  ls -lah scattered/
  echo "=== split_intervals: END ==="
  """
}


process split_bam_by_intervals {
  label 'process_high'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path tumor_bam
    path tumor_bam_index
    path interval_files

  output:
    path "shards/*.bam",       emit: bams
    path "shards/*.bam.bai",   emit: bais
    path "shards/*.intervals", emit: intervals

  script:
  """
  set -euo pipefail

  echo "=== split_bam_by_intervals: START ==="
  echo "PWD=\$(pwd)"
  echo "tumor_bam=$tumor_bam"
  echo "tumor_bam_index=$tumor_bam_index"
  echo "Interval files staged:"
  ls -lah *.intervals | head -20
  echo "Total interval files: \$(ls *.intervals | wc -l)"
  ls -lah
  mkdir -p shards

  total=\$(ls *.intervals | wc -l)
  count=0

  for interval_file in *.intervals; do
    shard_base=\$(basename "\$interval_file" .intervals)
    count=\$((count + 1))
    echo "--- Shard \${count}/\${total}: \${shard_base} ---"

    awk '!/^@/ { 
      split($1, a, /:|-/); 
      print a[1]"\t"(a[2]-1)"\t"a[3] 
    }' "$interval_file" > "${shard_base}.bed"
    echo "  BED file created: \$(wc -l < \${shard_base}.bed) regions"

    samtools view -b -L "\${shard_base}.bed" \
      -o "shards/\${shard_base}.bam" \
      "$tumor_bam"
    echo "  BAM written: \$(ls -lah shards/\${shard_base}.bam | awk '{print \$5}')"

    samtools index "shards/\${shard_base}.bam"
    echo "  BAM indexed"

    cp "\$interval_file" "shards/\${shard_base}.intervals"
    echo "  Interval copied"
  done

  echo "=== split_bam_by_intervals: DONE ==="
  echo "Final shards directory:"
  ls -lah shards/
  echo "Total BAMs: \$(ls shards/*.bam | wc -l)"
  """
}


process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"
  stageInMode 'symlink'

  input:
    tuple(
      path(interval_shard),
      path(tumor_bam),
      path(tumor_bam_index),
      path(ref_fasta),
      path(ref_fai),
      path(ref_dict),
      path(germline_resource)
    )
    path normal_bam
    path normal_bam_index
    path alleles_vcf
    path alleles_vcf_tbi
    val  tumor_sample
    val  extra_args

  output:
    path "*.vcf.gz",     emit: vcf
    path "*.vcf.gz.tbi", emit: tbi
    path "versions.yml", emit: versions

  shell:
  '''
  set -euo pipefail

  shard_base=$(basename "!{interval_shard}" .intervals)
  tumor_bam="!{tumor_bam}"
  tumor_bam_index="!{tumor_bam_index}"
  interval_shard="!{interval_shard}"
  ref_fasta="!{ref_fasta}"
  ref_fai="!{ref_fai}"
  ref_dict="!{ref_dict}"
  germline_resource="!{germline_resource}"
  normal_bam="!{normal_bam}"
  normal_bam_index="!{normal_bam_index}"
  alleles_vcf="!{alleles_vcf}"
  alleles_vcf_tbi="!{alleles_vcf_tbi}"
  tumor_sample="!{tumor_sample}"
  extra_args="!{extra_args}"

  heap_mb="!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }"

  echo "=== mutect_wrapper: START ==="
  echo "PWD=$(pwd)"
  echo "shard_base=${shard_base}"
  echo "tumor_bam=${tumor_bam}  size=$(ls -lah $tumor_bam | awk '{print $5}')"
  echo "interval_shard=${interval_shard}"
  echo "ref_fasta=${ref_fasta}"
  echo "germline_resource=${germline_resource}"
  echo "normal_bam=${normal_bam}"
  echo "alleles_vcf=${alleles_vcf}"
  echo "tumor_sample=${tumor_sample}"
  echo "heap_mb=${heap_mb}M"
  echo "extra_args='${extra_args}'"
  echo "Staged files:"
  ls -lah

  [[ -n "$tumor_sample" ]] || { echo "ERROR: tumor_sample is empty" >&2; exit 1; }

  # --- normal sample ---
  normal_args=""
  if [[ "$(basename "$normal_bam")" != "NO_NORMAL_BAM" ]]; then
    echo "Normal BAM provided: $normal_bam"
    normal_sample=$(samtools view -H "$normal_bam" \
      | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
      | sort -u)
    [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
    [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values in normal BAM" >&2; exit 1; }
    echo "normal_sample=${normal_sample}"
    normal_args="--input $normal_bam --normal-sample $normal_sample"
  else
    echo "NO_NORMAL_BAM sentinel -> tumor-only mode"
  fi

  # --- alleles ---
  alleles_args=""
  if [[ "$(basename "$alleles_vcf")" != "NO_ALLELES_VCF" ]]; then
    echo "Alleles VCF provided: $alleles_vcf"
    alleles_args="--alleles $alleles_vcf"
  else
    echo "NO_ALLELES_VCF sentinel -> no force-calling"
  fi

  # --- germline resource index ---
  echo "--- Checking germline resource index ---"
  if [[ ! -f "${germline_resource}.tbi" ]]; then
    echo "No .tbi found, creating with tabix..."
    tabix -f -p vcf "$germline_resource"
  fi
  test -s "${germline_resource}.tbi"
  echo "Germline resource index OK"

  out_prefix="out.${shard_base}"
  echo "=== Running Mutect2 (output: ${out_prefix}.vcf.gz) ==="

  gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \
    --input "$tumor_bam" \
    ${normal_args} \
    --reference "$ref_fasta" \
    --germline-resource "$germline_resource" \
    --intervals "$interval_shard" \
    --tmp-dir . \
    --tumor-sample "$tumor_sample" \
    ${alleles_args} \
    ${extra_args} \
    --output "${out_prefix}.vcf.gz"

  echo "=== Mutect2 finished ==="
  echo "Output files:"
  ls -lah ${out_prefix}*

  ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
  echo "=== mutect_wrapper: END ==="
  '''
}

process gather_vcfs {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path vcfs, arity: '1..*'

  output:
    path "started.txt"
    path "merged.vcf.gz"
    path "merged.vcf.gz.tbi"

  script:
  """
  set -euo pipefail
  echo "=== gather_vcfs: START ===" | tee started.txt
  echo "SCRIPT_STARTED \$(date)" >> started.txt
  echo "PWD=\$(pwd)"
  ls -lah

  find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list
  echo "VCFs found: \$(wc -l < vcfs.list)"
  cat vcfs.list

  sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list | sort -k1,1n | cut -f2- > vcfs.sorted.list
  awk '{print "--INPUT", \$0}' vcfs.sorted.list > gather.args
  echo "gather.args contents:"
  cat gather.args

  echo "--- Running GatherVcfs ---"
  time gatk GatherVcfs --arguments_file gather.args -O merged.vcf.gz
  echo "--- GatherVcfs done ---"

  if [ ! -s merged.vcf.gz.tbi ]; then
    echo "No index found, creating..."
    tabix -f -p vcf merged.vcf.gz || gatk IndexFeatureFile -I merged.vcf.gz
  fi
  test -s merged.vcf.gz.tbi

  echo "Final output:"
  ls -lah merged.vcf.gz merged.vcf.gz.tbi
  echo "=== gather_vcfs: END ==="
  """
}

workflow {

  log.info "=== WORKFLOW START ==="
  log.info "tumor_reads       : ${params.tumor_reads}"
  log.info "scatter_count     : ${params.scatter_count}"
  log.info "normal_reads      : ${params.normal_reads ?: 'NOT PROVIDED (tumor-only)'}"
  log.info "force_call_file   : ${params.force_call_file ?: 'NOT PROVIDED'}"
  log.info "m2_extra_args     : ${params.m2_extra_args ?: 'NONE'}"
  log.info "gatk_docker       : ${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0 (default)'}"

  // Extract tumor sample name from BAM header — no param needed
  // sample_name_res = get_tumor_sample_name(
  //   file(params.tumor_reads, checkIfExists: true),
  //   file(params.tumor_reads_index, checkIfExists: true)
  // )

  // tumor_sample_ch = sample_name_res.sample_name
  // .map { f -> 
  //   def s = f.text.trim()
  //   log.info "Tumor sample name: ${s}"
  //   s
  // }

  interval_res = split_intervals(
    file(params.ref_fasta),
    file(params.ref_fai),
    file(params.ref_dict),
    file(params.intervals),
    params.scatter_count as int
  )

  shard_res = split_bam_by_intervals(
    file(params.tumor_reads, checkIfExists: true),
    file(params.tumor_reads_index, checkIfExists: true),
    interval_res.interval_shards.collect()
  )

  // interval_res.interval_shards
  //   .flatten()
  //   .count()
  //   .view { n -> "=== split_intervals produced ${n} shards ===" }


  normal_bam_val      = params.normal_reads          ? file(params.normal_reads, checkIfExists: true)          : file(NO_NORMAL_BAM_PATH, checkIfExists: true)
  normal_bai_val      = params.normal_reads_index    ? file(params.normal_reads_index, checkIfExists: true)    : file(NO_NORMAL_BAI_PATH, checkIfExists: true)
  alleles_vcf_val     = params.force_call_file       ? file(params.force_call_file, checkIfExists: true)       : file(NO_ALLELES_VCF_PATH, checkIfExists: true)
  alleles_vcf_tbi_val = params.force_call_file_index ? file(params.force_call_file_index, checkIfExists: true) : file(NO_ALLELES_TBI_PATH, checkIfExists: true)

  log.info "normal_bam_val    : ${normal_bam_val}"
  log.info "alleles_vcf_val   : ${alleles_vcf_val}"

  if( !params.tumor_sample ) {
    error "Missing required --tumor_sample"
  }

  log.info "tumor_sample      : ${params.tumor_sample}"

  // Fan out one Mutect2 job per interval shard

  bams_ch = shard_res.bams.flatten()
  .map { bam -> tuple(bam.name.replace('.bam', ''), bam) }

  bais_ch = shard_res.bais.flatten()
    .map { bai -> tuple(bai.name.replace('.bam.bai', ''), bai) }

  intervals_ch = shard_res.intervals.flatten()
    .map { interval -> tuple(interval.name.replace('.intervals', ''), interval) }

  mutect_inputs_ch = bams_ch
    .join(bais_ch)
    .join(intervals_ch)
    .map { key, bam, bai, interval ->
      println "Preparing Mutect2 shard: ${key}"
      tuple(
        interval,
        bam,
        bai,
        file(params.ref_fasta),
        file(params.ref_fai),
        file(params.ref_dict),
        file(params.germline_resource)
      )
    }

  mutect_res = mutect_wrapper(
    mutect_inputs_ch,
    Channel.value(normal_bam_val),
    Channel.value(normal_bai_val),
    Channel.value(alleles_vcf_val),
    Channel.value(alleles_vcf_tbi_val),
    Channel.value(params.tumor_sample),
    params.m2_extra_args ?: ''
  )

  mutect_res.vcf
    .count()
    .view { n -> "=== mutect_wrapper finished: ${n} VCFs produced ===" }

  gather_vcfs(mutect_res.vcf.collect())

  log.info "=== WORKFLOW SUBMITTED ==="
}
// /*
//  * --------------------------------------------
//  * Defaults / params
//  * --------------------------------------------
//  */

// if( !params.containsKey('m2_extra_args') )
//   params.m2_extra_args = ''

// if( !params.containsKey('normal_reads') )
//   params.normal_reads = null
// if( !params.containsKey('normal_reads_index') )
//   params.normal_reads_index = null
// if( !params.containsKey('force_call_file') )
//   params.force_call_file = null
// if( !params.containsKey('force_call_file_index') )
//   params.force_call_file_index = null

// def NO_NORMAL_BAM_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAM"
// def NO_NORMAL_BAI_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAI"
// def NO_ALLELES_VCF_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
// def NO_ALLELES_TBI_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_TBI"

// process split_intervals {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path ref_fasta
//     path ref_fai
//     path ref_dict
//     path intervals
//     val scatter_count

//   output:
//     path "scattered/*.intervals", emit: interval_shards

//   script:
//   """
//   set -euo pipefail

//   echo "=== split_intervals: START ==="
//   echo "PWD=\$(pwd)"
//   echo "ref_fasta=$ref_fasta"
//   echo "intervals=$intervals"
//   echo "scatter_count=$scatter_count"
//   ls -lah
//   mkdir -p scattered

//   echo "--- Running BedToIntervalList ---"
//   gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \
//     -I "$intervals" \
//     -SD "$ref_dict" \
//     -O regions.interval_list
//   echo "--- BedToIntervalList done ---"

//   echo "--- Running SplitIntervals ---"
//   gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \
//     -R "$ref_fasta" \
//     -L regions.interval_list \
//     --scatter "$scatter_count" \
//     -O scattered
//   echo "--- SplitIntervals done ---"

//   echo "Shards produced:"
//   ls -lah scattered/
//   echo "=== split_intervals: END ==="
//   """
// }

// process split_bam_by_intervals {
//   label 'process_high'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path tumor_bam
//     path tumor_bam_index
//     path interval_files

//   output:
//     tuple path("shards/*.bam"), path("shards/*.bam.bai"), path("shards/*.intervals")

//   script:
//   """
//   set -euo pipefail

//   echo "=== split_bam_by_intervals: START ==="
//   echo "PWD=\$(pwd)"
//   echo "tumor_bam=$tumor_bam"
//   echo "tumor_bam_index=$tumor_bam_index"
//   echo "Interval files staged:"
//   ls -lah *.intervals | head -20
//   echo "Total interval files: \$(ls *.intervals | wc -l)"
//   ls -lah
//   mkdir -p shards

//   total=\$(ls *.intervals | wc -l)
//   count=0

//   for interval_file in *.intervals; do
//     shard_base=\$(basename "\$interval_file" .intervals)
//     count=\$((count + 1))
//     echo "--- Shard \${count}/\${total}: \${shard_base} ---"

//     awk 'BEGIN{OFS="\\t"} !/^@/ {print \$1, \$2-1, \$3}' "\$interval_file" > "\${shard_base}.bed"
//     echo "  BED file created: \$(wc -l < \${shard_base}.bed) regions"

//     samtools view -b -L "\${shard_base}.bed" \
//       -o "shards/\${shard_base}.bam" \
//       "$tumor_bam"
//     echo "  BAM written: \$(ls -lah shards/\${shard_base}.bam | awk '{print \$5}')"

//     samtools index "shards/\${shard_base}.bam"
//     echo "  BAM indexed"

//     cp "\$interval_file" "shards/\${shard_base}.intervals"
//     echo "  Interval copied"
//   done

//   echo "=== split_bam_by_intervals: DONE ==="
//   echo "Final shards directory:"
//   ls -lah shards/
//   echo "Total BAMs: \$(ls shards/*.bam | wc -l)"
//   """
// }

// process mutect_wrapper {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     tuple(
//       path(tumor_bam),
//       path(tumor_bam_index),
//       path(interval_shard),
//       path(ref_fasta),
//       path(ref_fai),
//       path(ref_dict),
//       path(germline_resource)
//     )

//     path normal_bam
//     path normal_bam_index
//     path alleles_vcf
//     path alleles_vcf_tbi

//     val extra_args

//   output:
//     path "*.vcf.gz",     emit: vcf
//     path "*.vcf.gz.tbi", emit: tbi
//     path "versions.yml", emit: versions

//   shell:
//   '''
//   set -euo pipefail

//   shard_base=$(basename "!{interval_shard}" .intervals)
//   tumor_bam="!{tumor_bam}"
//   tumor_bam_index="!{tumor_bam_index}"
//   interval_shard="!{interval_shard}"
//   ref_fasta="!{ref_fasta}"
//   ref_fai="!{ref_fai}"
//   ref_dict="!{ref_dict}"
//   germline_resource="!{germline_resource}"
//   normal_bam="!{normal_bam}"
//   normal_bam_index="!{normal_bam_index}"
//   alleles_vcf="!{alleles_vcf}"
//   alleles_vcf_tbi="!{alleles_vcf_tbi}"
//   extra_args="!{extra_args}"

//   avail_mem_mb="!{ task.memory ? (task.memory.mega * 0.8).intValue() : 3072 }"
//   heap_mb="!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }"

//   echo "=== mutect_wrapper: START ==="
//   echo "PWD=$(pwd)"
//   echo "shard_base=${shard_base}"
//   echo "tumor_bam=${tumor_bam}  size=$(ls -lah $tumor_bam | awk '{print $5}')"
//   echo "interval_shard=${interval_shard}"
//   echo "ref_fasta=${ref_fasta}"
//   echo "germline_resource=${germline_resource}"
//   echo "normal_bam=${normal_bam}"
//   echo "alleles_vcf=${alleles_vcf}"
//   echo "heap_mb=${heap_mb}M"
//   echo "extra_args='${extra_args}'"
//   echo "Staged files:"
//   ls -lah

//   tumor_sample="!{params.tumor_sample_name ?: ''}"
//   [[ -n "$tumor_sample" ]] || { echo "ERROR: tumor_sample_name not provided" >&2; exit 1; }
//   echo "tumor_sample=${tumor_sample}"

//   # --- normal sample ---
//   normal_args=""
//   if [[ "$(basename "$normal_bam")" != "NO_NORMAL_BAM" ]]; then
//     echo "Normal BAM provided: $normal_bam"
//     normal_sample=$(samtools view -H "$normal_bam" \
//       | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
//       | sort -u)
//     [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
//     [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values in normal BAM" >&2; exit 1; }
//     echo "normal_sample=${normal_sample}"
//     normal_args="--input $normal_bam --normal-sample $normal_sample"
//   else
//     echo "NO_NORMAL_BAM sentinel -> tumor-only mode"
//   fi

//   # --- alleles ---
//   alleles_args=""
//   if [[ "$(basename "$alleles_vcf")" != "NO_ALLELES_VCF" ]]; then
//     echo "Alleles VCF provided: $alleles_vcf"
//     alleles_args="--alleles $alleles_vcf"
//   else
//     echo "NO_ALLELES_VCF sentinel -> no force-calling"
//   fi

//   # --- germline resource index ---
//   echo "--- Checking germline resource index ---"
//   if [[ ! -f "${germline_resource}.tbi" ]]; then
//     echo "No .tbi found, creating with tabix..."
//     tabix -f -p vcf "$germline_resource"
//   fi
//   test -s "${germline_resource}.tbi"
//   echo "Germline resource index OK"

//   out_prefix="out.${shard_base}"
//   echo "=== Running Mutect2 (output: ${out_prefix}.vcf.gz) ==="

//   gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \
//     --input "$tumor_bam" \
//     ${normal_args} \
//     --reference "$ref_fasta" \
//     --germline-resource "$germline_resource" \
//     --intervals "$interval_shard" \
//     --tmp-dir . \
//     --tumor-sample "$tumor_sample" \
//     ${alleles_args} \
//     ${extra_args} \
//     --output "${out_prefix}.vcf.gz"

//   echo "=== Mutect2 finished ==="
//   echo "Output files:"
//   ls -lah ${out_prefix}*

//   ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
//   echo "=== mutect_wrapper: END ==="
//   '''
// }

// process gather_vcfs {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path vcfs, arity: '1..*'

//   output:
//     path "started.txt"
//     path "merged.vcf.gz"
//     path "merged.vcf.gz.tbi"

//   script:
//   """
//   set -euo pipefail
//   echo "=== gather_vcfs: START ===" | tee started.txt
//   echo "SCRIPT_STARTED \$(date)" >> started.txt
//   echo "PWD=\$(pwd)"
//   ls -lah

//   find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list
//   echo "VCFs found: \$(wc -l < vcfs.list)"
//   cat vcfs.list

//   sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list | sort -k1,1n | cut -f2- > vcfs.sorted.list
//   awk '{print "--INPUT", \$0}' vcfs.sorted.list > gather.args
//   echo "gather.args contents:"
//   cat gather.args

//   echo "--- Running GatherVcfs ---"
//   time gatk GatherVcfs --arguments_file gather.args -O merged.vcf.gz
//   echo "--- GatherVcfs done ---"

//   if [ ! -s merged.vcf.gz.tbi ]; then
//     echo "No index found, creating..."
//     tabix -f -p vcf merged.vcf.gz || gatk IndexFeatureFile -I merged.vcf.gz
//   fi
//   test -s merged.vcf.gz.tbi

//   echo "Final output:"
//   ls -lah merged.vcf.gz merged.vcf.gz.tbi
//   echo "=== gather_vcfs: END ==="
//   """
// }

// workflow {

//   log.info "=== WORKFLOW START ==="
//   log.info "tumor_reads       : ${params.tumor_reads}"
//   log.info "tumor_sample_name : ${params.tumor_sample_name}"
//   log.info "scatter_count     : ${params.scatter_count}"
//   log.info "normal_reads      : ${params.normal_reads ?: 'NOT PROVIDED (tumor-only)'}"
//   log.info "force_call_file   : ${params.force_call_file ?: 'NOT PROVIDED'}"
//   log.info "m2_extra_args     : ${params.m2_extra_args ?: 'NONE'}"
//   log.info "gatk_docker       : ${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0 (default)'}"

//   interval_res = split_intervals(
//     file(params.ref_fasta),
//     file(params.ref_fai),
//     file(params.ref_dict),
//     file(params.intervals),
//     params.scatter_count as int
//   )

//   // Log how many shards were produced
//   interval_res.interval_shards
//     .flatten()
//     .count()
//     .view { n -> "=== split_intervals produced ${n} shards ===" }

//   shard_res = split_bam_by_intervals(
//     file(params.tumor_reads, checkIfExists: true),
//     file(params.tumor_reads_index, checkIfExists: true),
//     interval_res.interval_shards.collect()
//   )

//   mutect_inputs_ch = shard_res
//     .transpose()
//     .map { bam, bai, interval ->
//       log.info "Queuing Mutect2 for shard: ${interval.name}"
//       tuple(
//         bam,
//         bai,
//         interval,
//         file(params.ref_fasta),
//         file(params.ref_fai),
//         file(params.ref_dict),
//         file(params.germline_resource)
//       )
//     }

//   normal_bam_val      = params.normal_reads          ? file(params.normal_reads, checkIfExists: true)          : file("NO_NORMAL_BAM")
//   normal_bai_val      = params.normal_reads_index    ? file(params.normal_reads_index, checkIfExists: true)    : file("NO_NORMAL_BAI")
//   alleles_vcf_val     = params.force_call_file       ? file(params.force_call_file, checkIfExists: true)       : file("NO_ALLELES_VCF")
//   alleles_vcf_tbi_val = params.force_call_file_index ? file(params.force_call_file_index, checkIfExists: true) : file("NO_ALLELES_TBI")

//   log.info "normal_bam_val    : ${normal_bam_val}"
//   log.info "alleles_vcf_val   : ${alleles_vcf_val}"

//   mutect_res = mutect_wrapper(
//     mutect_inputs_ch,
//     Channel.value(normal_bam_val),
//     Channel.value(normal_bai_val),
//     Channel.value(alleles_vcf_val),
//     Channel.value(alleles_vcf_tbi_val),
//     params.m2_extra_args ?: ''
//   )

//   mutect_res.vcf
//     .count()
//     .view { n -> "=== mutect_wrapper finished: ${n} VCFs produced ===" }

//   gather_vcfs(mutect_res.vcf.collect())

//   log.info "=== WORKFLOW SUBMITTED ==="
// }
