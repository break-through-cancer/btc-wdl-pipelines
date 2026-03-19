/*
 * --------------------------------------------
 * Defaults / params
 * --------------------------------------------
 */

// Avoid "Access to undefined parameter" warnings:
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

/*
 * --------------------------------------------
 * Sentinel files (distinct basenames to avoid input name collisions)
 *   - We create them at runtime so you don't have to commit assets/*
 * --------------------------------------------
 */
def NO_NORMAL_BAM_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAM"
def NO_NORMAL_BAI_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAI"
def NO_ALLELES_VCF_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
def NO_ALLELES_TBI_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_TBI"

def NO_NORMAL_BAM  = null
def NO_NORMAL_BAI  = null
def NO_ALLELES_VCF = null
def NO_ALLELES_TBI = null

/*
 * --------------------------------------------
 * split_intervals
 * --------------------------------------------
 */
// process split_intervals {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path ref_fasta
//     path ref_fai
//     path ref_dict
//     path intervals
//     val  scatter_count

//   output:
//     path "scattered/*.intervals", emit: shards

//   script:
//   """
//   set -euo pipefail
//   echo "=== split_intervals: START ==="
//   echo "PWD=\$(pwd)"
//   ls -lah

//   mkdir -p scattered

//   gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \\
//     -I "$intervals" \\
//     -SD "$ref_dict" \\
//     -O regions.interval_list

//   gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \\
//     -R "$ref_fasta" \\
//     -L regions.interval_list \\
//     --scatter "$scatter_count" \\
//     -O scattered


//   echo "=== split_intervals: END ==="
//   """
// }

/**
 * --------------------------------------------
 * do view to subset the bed file and stuff, matched up bed and bam files 
 * keep prefixes from the splitintervals bam 
* 
 * --------------------------------------------
 */

process prepare_shards_and_subset_tumor {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path tumor_bam
    path tumor_bam_index
    path ref_fasta
    path ref_fai
    path ref_dict
    path intervals
    val  scatter_count

  output:
    path "shards/*.intervals", emit: interval_shards
    path "shards/*.bam", emit: shard_bams
    path "shards/*.bam.bai", emit: shard_bais
    path "shards/manifest.tsv", emit: manifest

  script:
  """
  set -euo pipefail

  echo "=== prepare_shards_and_subset_tumor: START ==="
  echo "PWD=\$(pwd)"
  ls -lah

  mkdir -p scattered
  mkdir -p shards
  mkdir -p beds

  # Convert BED -> interval_list for GATK SplitIntervals
  gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \\
    -I "$intervals" \\
    -SD "$ref_dict" \\
    -O regions.interval_list

  # Split intervals into shard files
  gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \\
    -R "$ref_fasta" \\
    -L regions.interval_list \\
    --scatter "$scatter_count" \\
    -O scattered

  # Build one small BAM per shard and keep matching prefixes
  : > shards/manifest.tsv
  echo -e "shard_base\\tinterval\\tbam\\tbai" >> shards/manifest.tsv

  for interval_file in scattered/*.intervals; do
    shard_base=\$(basename "\$interval_file" .intervals)

    cp "\$interval_file" "shards/\${shard_base}.intervals"

    # Convert GATK interval_list -> BED for samtools
    awk 'BEGIN{OFS="\\t"} !/^@/ {print \$1, \$2-1, \$3}' "\$interval_file" > "beds/\${shard_base}.bed"

    samtools view \\
      -b \\
      -L "beds/\${shard_base}.bed" \\
      -o "shards/\${shard_base}.bam" \\
      "$tumor_bam"

    samtools index "shards/\${shard_base}.bam"

    echo -e "\${shard_base}\\tshards/\${shard_base}.intervals\\tshards/\${shard_base}.bam\\tshards/\${shard_base}.bam.bai" >> shards/manifest.tsv
  done

  echo "=== shard outputs ==="
  ls -lah shards
  echo "=== manifest ==="
  cat shards/manifest.tsv
  echo "=== prepare_shards_and_subset_tumor: END ==="
  """
}
/*
 * --------------------------------------------
 * mutect_wrapper (NO optional inputs; uses sentinel files)
 * --------------------------------------------
 */
process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple(
      val(tumor_bam),
      val(tumor_bam_index),
      path(interval_shard),
      path(ref_fasta),
      path(ref_fai),
      path(ref_dict),
      path(germline_resource)
    )

    // always provided: either real file, or sentinel file (unique per input)
    path normal_bam
    path normal_bam_index
    path alleles_vcf
    path alleles_vcf_tbi

    val extra_args

  output:
    path "*.vcf.gz",     emit: vcf
    path "*.vcf.gz.tbi", emit: tbi
    path "versions.yml", emit: versions

  shell:
  '''
  set -euo pipefail

  # Bind Nextflow inputs to bash vars (shell: does NOT auto-export them)
  shard_base="!{shard_base}"
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

  extra_args="!{extra_args}"

  # Heap sizing
  avail_mem_mb="!{ task.memory ? (task.memory.mega * 0.8).intValue() : 3072 }"
  heap_mb="!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }"

  echo "=== mutect_wrapper: START ==="
  echo "PWD=$(pwd)"
  echo "Task memory avail=${avail_mem_mb}M ; heap_mb=${heap_mb}M"
  echo "extra_args='$extra_args'"
  ls -lah

  tumor_sample="!{params.tumor_sample_name ?: ''}"
  [[ -n "$tumor_sample" ]] || { echo "ERROR: tumor_sample_name not provided" >&2; exit 1; }
  echo "Using tumor sample for Mutect2: $tumor_sample"



  # --- normal sample (if provided; sentinel name == NO_NORMAL_BAM) ---
  normal_args=""
  if [[ "$(basename "$normal_bam")" != "NO_NORMAL_BAM" ]]; then
    echo "Normal BAM provided: $normal_bam"

    normal_sample=$(samtools view -H "$normal_bam" \
      | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
      | sort -u)

    [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
    [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values found in normal BAM header:" >&2; echo "$normal_sample" >&2; exit 1; }

    normal_args="--input $normal_bam --normal-sample $normal_sample"
  else
    echo "NO_NORMAL_BAM sentinel -> tumor-only mode."
  fi

  # --- alleles (if provided; sentinel name == NO_ALLELES_VCF) ---
  alleles_args=""
  if [[ "$(basename "$alleles_vcf")" != "NO_ALLELES_VCF" ]]; then
    echo "Alleles VCF provided: $alleles_vcf"
    alleles_args="--alleles $alleles_vcf"
  else
    echo "NO_ALLELES_VCF sentinel -> no force-calling."
  fi

  # --- ensure germline resource has tabix index ---
  if [[ ! -f "${germline_resource}.tbi" ]]; then
    echo "No .tbi found for germline_resource; creating with tabix..."
    tabix -f -p vcf "$germline_resource"
  fi
  test -s "${germline_resource}.tbi"

  out_prefix="out.${shard_base}"

  echo "=== mutect_wrapper: RUN Mutect2 ==="
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

  ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
  echo "=== mutect_wrapper: END ==="
  '''
}

/*
 * --------------------------------------------
 * gather_vcfs  (no optional outputs to avoid parser weirdness)
 * --------------------------------------------
 */
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
  echo "SCRIPT_STARTED \$(date)" > started.txt
  echo "PWD=\$(pwd)"
  ls -lah

  find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list
  echo "VCFs found: \$(wc -l < vcfs.list)"

  sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list | sort -k1,1n | cut -f2- > vcfs.sorted.list
  awk '{print "--INPUT", \$0}' vcfs.sorted.list > gather.args

  time gatk GatherVcfs --arguments_file gather.args -O merged.vcf.gz

  if [ ! -s merged.vcf.gz.tbi ]; then
    tabix -f -p vcf merged.vcf.gz || gatk IndexFeatureFile -I merged.vcf.gz
  fi
  test -s merged.vcf.gz.tbi

  ls -lah merged.vcf.gz merged.vcf.gz.tbi
  """
}

/*
 * --------------------------------------------
 * workflow
 * --------------------------------------------
 */
workflow {

  /*
   * Create sentinel files in the workflow repo directory.
   * This runs on the "driver" (not in a container), before tasks are scheduled.
   */
  new File("${workflow.projectDir}/assets").mkdirs()

  def mkEmpty = { String p ->
    def f = new File(p)
    if( !f.exists() ) f.text = ""
    return file(p, checkIfExists: true)
  }

  NO_NORMAL_BAM  = mkEmpty(NO_NORMAL_BAM_PATH)
  NO_NORMAL_BAI  = mkEmpty(NO_NORMAL_BAI_PATH)
  NO_ALLELES_VCF = mkEmpty(NO_ALLELES_VCF_PATH)
  NO_ALLELES_TBI = mkEmpty(NO_ALLELES_TBI_PATH)

  interval_ch = prep.out.interval_shards
    .flatten()
    .map { f -> tuple(f.baseName, f) }

  bam_ch = prep.out.shard_bams
    .flatten()
    .map { f -> tuple(f.baseName, f) }

  bai_ch = prep.out.shard_bais
    .flatten()
    .map { f ->
      def key = f.name.replaceFirst(/\\.bam\\.bai$/, '')
      tuple(key, f)
    }

  mutect_inputs = bam_ch
    .join(bai_ch)
    .join(interval_ch)
    .map { key, bam, bai, interval ->
      tuple(
        key,
        bam,
        bai,
        interval,
        file(params.ref_fasta),
        file(params.ref_fai),
        file(params.ref_dict),
        file(params.germline_resource)
      )
    }

  normal_bam_val      = params.normal_reads          ? file(params.normal_reads)          : NO_NORMAL_BAM
  normal_bai_val      = params.normal_reads_index    ? file(params.normal_reads_index)    : NO_NORMAL_BAI
  alleles_vcf_val     = params.force_call_file       ? file(params.force_call_file)       : NO_ALLELES_VCF
  alleles_vcf_tbi_val = params.force_call_file_index ? file(params.force_call_file_index) : NO_ALLELES_TBI

  mutect_res = mutect_wrapper(
    mutect_inputs,
    Channel.value(normal_bam_val),
    Channel.value(normal_bai_val),
    Channel.value(alleles_vcf_val),
    Channel.value(alleles_vcf_tbi_val),
    params.m2_extra_args
  )

  mutect_res.vcf.view { "VCF: $it" }
  gather_vcfs(mutect_res.vcf.collect())
}



                                                                                                           // /*
//  * --------------------------------------------
//  * Defaults / params
//  * --------------------------------------------
//  */

// // Avoid "Access to undefined parameter" warnings:
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

// /*
//  * --------------------------------------------
//  * Sentinel files (distinct basenames to avoid input name collisions)
//  *   - We create them at runtime so you don't have to commit assets/*
//  * --------------------------------------------
//  */
// def NO_NORMAL_BAM_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAM"
// def NO_NORMAL_BAI_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAI"
// def NO_ALLELES_VCF_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
// def NO_ALLELES_TBI_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_TBI"

// def NO_NORMAL_BAM  = null
// def NO_NORMAL_BAI  = null
// def NO_ALLELES_VCF = null
// def NO_ALLELES_TBI = null

// /*
//  * --------------------------------------------
//  * split_intervals
//  * --------------------------------------------
//  */
// process split_intervals {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path ref_fasta
//     path ref_fai
//     path ref_dict
//     path intervals
//     val  scatter_count

//   output:
//     path "scattered/*.intervals", emit: shards

//   script:
//   """
//   set -euo pipefail
//   echo "=== split_intervals: START ==="
//   echo "PWD=\$(pwd)"
//   ls -lah

//   mkdir -p scattered

//   gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \\
//     -I "$intervals" \\
//     -SD "$ref_dict" \\
//     -O regions.interval_list

//   gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \\
//     -R "$ref_fasta" \\
//     -L regions.interval_list \\
//     --scatter "$scatter_count" \\
//     -O scattered

//   echo "=== split_intervals: END ==="
//   """
// }

// /*
//  * --------------------------------------------
//  * mutect_wrapper (NO optional inputs; uses sentinel files)
//  * --------------------------------------------
//  */
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

//     // always provided: either real file, or sentinel file (unique per input)
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

//   # Bind Nextflow inputs to bash vars (shell: does NOT auto-export them)
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

//   # Heap sizing
//   avail_mem_mb="!{ task.memory ? (task.memory.mega * 0.8).intValue() : 3072 }"
//   heap_mb="!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }"

//   echo "=== mutect_wrapper: START ==="
//   echo "PWD=$(pwd)"
//   echo "Task memory avail=${avail_mem_mb}M ; heap_mb=${heap_mb}M"
//   echo "extra_args='$extra_args'"
//   ls -lah

//   # --- tumor sample name (SM tag) ---
//   tumor_samples=$(samtools view -H "$tumor_bam" \
//     | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
//     | sort -u)

//   [[ -n "$tumor_samples" ]] || { echo "ERROR: No SM tag found in tumor BAM header" >&2; exit 1; }

//   tumor_sample_count=$(echo "$tumor_samples" | wc -l | tr -d ' ')
//   if [[ "$tumor_sample_count" -eq 1 ]]; then
//     tumor_sample="$tumor_samples"
//   else
//     echo "WARN: Multiple SM values found in tumor BAM header:" >&2
//     echo "$tumor_samples" >&2

//     # Optional override: if user supplies extra_args like: --tumor-sample-name <SM>
//     # Or you can wire a real Nextflow param (recommended) - see note below.
//     preferred="!{params.tumor_sample_name ?: ''}"

//     if [[ -n "$preferred" ]] && echo "$tumor_samples" | grep -Fxq "$preferred"; then
//       tumor_sample="$preferred"
//       echo "Using user-specified tumor_sample_name: $tumor_sample"
//     else
//       # deterministic fallback: pick first in sorted list
//       tumor_sample="$(echo "$tumor_samples" | head -n 1)"
//       echo "Using first tumor SM (fallback): $tumor_sample"
//     fi
//   fi

//   echo "Detected tumor sample used for Mutect2: $tumor_sample"


//   # --- normal sample (if provided; sentinel name == NO_NORMAL_BAM) ---
//   normal_args=""
//   if [[ "$(basename "$normal_bam")" != "NO_NORMAL_BAM" ]]; then
//     echo "Normal BAM provided: $normal_bam"

//     normal_sample=$(samtools view -H "$normal_bam" \
//       | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
//       | sort -u)

//     [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
//     [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values found in normal BAM header:" >&2; echo "$normal_sample" >&2; exit 1; }

//     normal_args="--input $normal_bam --normal-sample $normal_sample"
//   else
//     echo "NO_NORMAL_BAM sentinel -> tumor-only mode."
//   fi

//   # --- alleles (if provided; sentinel name == NO_ALLELES_VCF) ---
//   alleles_args=""
//   if [[ "$(basename "$alleles_vcf")" != "NO_ALLELES_VCF" ]]; then
//     echo "Alleles VCF provided: $alleles_vcf"
//     alleles_args="--alleles $alleles_vcf"
//   else
//     echo "NO_ALLELES_VCF sentinel -> no force-calling."
//   fi

//   # --- ensure germline resource has tabix index ---
//   if [[ ! -f "${germline_resource}.tbi" ]]; then
//     echo "No .tbi found for germline_resource; creating with tabix..."
//     tabix -f -p vcf "$germline_resource"
//   fi
//   test -s "${germline_resource}.tbi"

//   shard_base=$(basename "$interval_shard" .intervals)
//   out_prefix="out.${shard_base}"

//   echo "=== mutect_wrapper: RUN Mutect2 ==="
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

//   ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
//   echo "=== mutect_wrapper: END ==="
//   '''
// }

// /*
//  * --------------------------------------------
//  * gather_vcfs  (no optional outputs to avoid parser weirdness)
//  * --------------------------------------------
//  */
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
//   echo "SCRIPT_STARTED \$(date)" > started.txt
//   echo "PWD=\$(pwd)"
//   ls -lah

//   find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list
//   echo "VCFs found: \$(wc -l < vcfs.list)"

//   sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list | sort -k1,1n | cut -f2- > vcfs.sorted.list
//   awk '{print "--INPUT", \$0}' vcfs.sorted.list > gather.args

//   time gatk GatherVcfs --arguments_file gather.args -O merged.vcf.gz

//   if [ ! -s merged.vcf.gz.tbi ]; then
//     tabix -f -p vcf merged.vcf.gz || gatk IndexFeatureFile -I merged.vcf.gz
//   fi
//   test -s merged.vcf.gz.tbi

//   ls -lah merged.vcf.gz merged.vcf.gz.tbi
//   """
// }

// /*
//  * --------------------------------------------
//  * workflow
//  * --------------------------------------------
//  */
// workflow {

//   /*
//    * Create sentinel files in the workflow repo directory.
//    * This runs on the "driver" (not in a container), before tasks are scheduled.
//    */
//   new File("${workflow.projectDir}/assets").mkdirs()

//   def mkEmpty = { String p ->
//     def f = new File(p)
//     if( !f.exists() ) f.text = ""
//     return file(p, checkIfExists: true)
//   }

//   NO_NORMAL_BAM  = mkEmpty(NO_NORMAL_BAM_PATH)
//   NO_NORMAL_BAI  = mkEmpty(NO_NORMAL_BAI_PATH)
//   NO_ALLELES_VCF = mkEmpty(NO_ALLELES_VCF_PATH)
//   NO_ALLELES_TBI = mkEmpty(NO_ALLELES_TBI_PATH)

//   shards_ch = split_intervals(
//     params.ref_fasta,
//     params.ref_fai,
//     params.ref_dict,
//     params.intervals,
//     params.scatter_count as int
//   ).shards.flatten()

//   base_inputs = shards_ch.map { shard ->
//     tuple(
//       params.tumor_reads,
//       params.tumor_reads_index,
//       shard,
//       params.ref_fasta,
//       params.ref_fai,
//       params.ref_dict,
//       params.germline_resource
//     )
//   }

//   // Choose real file if provided; else sentinel. Wrap with file(...) so Nextflow stages it.
//   normal_bam_val      = params.normal_reads          ? file(params.normal_reads)          : NO_NORMAL_BAM
//   normal_bai_val      = params.normal_reads_index    ? file(params.normal_reads_index)    : NO_NORMAL_BAI
//   alleles_vcf_val     = params.force_call_file       ? file(params.force_call_file)       : NO_ALLELES_VCF
//   alleles_vcf_tbi_val = params.force_call_file_index ? file(params.force_call_file_index) : NO_ALLELES_TBI

//   mutect_res = mutect_wrapper(
//     base_inputs,
//     Channel.value(normal_bam_val),
//     Channel.value(normal_bai_val),
//     Channel.value(alleles_vcf_val),
//     Channel.value(alleles_vcf_tbi_val),
//     params.m2_extra_args
//   )

//   mutect_res.vcf.view { "VCF: $it" }
//   gather_vcfs(mutect_res.vcf.collect())
// }

