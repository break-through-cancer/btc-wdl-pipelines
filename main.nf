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
 * Sentinel file (for "no optional input")
 *   - We create it at runtime so you don't have to commit assets/NO_FILE
 * --------------------------------------------
 */
def NO_FILE_PATH = "${workflow.projectDir}/assets/NO_FILE"
def NO_FILE = null

/*
 * --------------------------------------------
 * split_intervals
 * --------------------------------------------
 */
process split_intervals {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path ref_fasta
    path ref_fai
    path ref_dict
    path intervals
    val  scatter_count

  output:
    path "scattered/*.intervals", emit: shards

  script:
  """
  set -euo pipefail
  echo "=== split_intervals: START ==="
  echo "PWD=\$(pwd)"
  ls -lah

  mkdir -p scattered

  gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \\
    -I "$intervals" \\
    -SD "$ref_dict" \\
    -O regions.interval_list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \\
    -R "$ref_fasta" \\
    -L regions.interval_list \\
    --scatter "$scatter_count" \\
    -O scattered

  echo "=== split_intervals: END ==="
  """
}

/*
 * --------------------------------------------
 * mutect_wrapper (NO optional inputs; uses sentinel NO_FILE)
 * --------------------------------------------
 */
process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple(
      path(tumor_bam),
      path(tumor_bam_index),
      path(interval_shard),
      path(ref_fasta),
      path(ref_fai),
      path(ref_dict),
      path(germline_resource)
    )

    // always provided: either real file, or NO_FILE sentinel
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

  avail_mem_mb=!{ task.memory ? (task.memory.mega * 0.8).intValue() : 3072 }
  heap_mb=!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }

  echo "=== mutect_wrapper: START ==="
  echo "PWD=$(pwd)"
  echo "Task memory avail=${avail_mem_mb}M ; heap_mb=${heap_mb}M"
  echo "extra_args='!{extra_args}'"
  ls -lah

  # --- tumor sample name (SM tag) ---
  tumor_sample=$(samtools view -H "$tumor_bam" \
    | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
    | sort -u)

  [[ -n "$tumor_sample" ]] || { echo "ERROR: No SM tag found in tumor BAM header" >&2; exit 1; }
  [[ $(echo "$tumor_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values found in tumor BAM header:" >&2; echo "$tumor_sample" >&2; exit 1; }
  echo "Detected tumor sample: $tumor_sample"

  # --- normal sample (if provided; sentinel name == NO_FILE) ---
  normal_args=""
  if [[ "$(basename "$normal_bam")" != "NO_FILE" ]]; then
    echo "Normal BAM provided: $normal_bam"

    normal_sample=$(samtools view -H "$normal_bam" \
      | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
      | sort -u)

    [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
    [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values found in normal BAM header:" >&2; echo "$normal_sample" >&2; exit 1; }

    normal_args="--input $normal_bam --normal-sample $normal_sample"
  else
    echo "NO_FILE sentinel for normal -> tumor-only mode."
  fi

  # --- alleles (if provided; sentinel name == NO_FILE) ---
  alleles_args=""
  if [[ "$(basename "$alleles_vcf")" != "NO_FILE" ]]; then
    echo "Alleles VCF provided: $alleles_vcf"
    alleles_args="--alleles $alleles_vcf"
  else
    echo "NO_FILE sentinel for alleles -> no force-calling."
  fi

  shard_base=$(basename "$interval_shard" .intervals)
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
    !{extra_args} \
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
    tabix -f -p vcf merged.vcf.gz || gatk IndexFeatureFile -F merged.vcf.gz
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
   * Create the sentinel file in the workflow repo directory.
   * This runs on the "driver" (not in a container), before tasks are scheduled.
   */
  new File("${workflow.projectDir}/assets").mkdirs()
  def nf = new File(NO_FILE_PATH)
  if( !nf.exists() ) {
    nf.text = ""   // create empty file
  }
  NO_FILE = file(NO_FILE_PATH, checkIfExists: true)

  shards_ch = split_intervals(
    params.ref_fasta,
    params.ref_fai,
    params.ref_dict,
    params.intervals,
    params.scatter_count as int
  ).shards.flatten()

  base_inputs = shards_ch.map { shard ->
    tuple(
      params.tumor_reads,
      params.tumor_reads_index,
      shard,
      params.ref_fasta,
      params.ref_fai,
      params.ref_dict,
      params.germline_resource
    )
  }

  // Choose real file if provided; else sentinel. Wrap with file(...) so Nextflow stages it.
  normal_bam_val       = params.normal_reads         ? file(params.normal_reads)          : NO_FILE
  normal_bai_val       = params.normal_reads_index   ? file(params.normal_reads_index)    : NO_FILE
  alleles_vcf_val      = params.force_call_file      ? file(params.force_call_file)       : NO_FILE
  alleles_vcf_tbi_val  = params.force_call_file_index? file(params.force_call_file_index) : NO_FILE

  mutect_res = mutect_wrapper(
    base_inputs,
    Channel.value(normal_bam_val),
    Channel.value(normal_bai_val),
    Channel.value(alleles_vcf_val),
    Channel.value(alleles_vcf_tbi_val),
    params.m2_extra_args
  )

  mutect_res.vcf.view { "VCF: $it" }
  gather_vcfs(mutect_res.vcf.collect())
}


// params.m2_extra_args = params.m2_extra_args ?: ''
// def do_force = params.force_call_file

// /*
//  * --------------------------------------------
//  * split_intervals
//  *   - BED -> interval_list (needs .dict)
//  *   - interval_list -> scattered/*.intervals
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
//   echo "Inputs:"
//   ls -lah
//   echo "ref_fasta=$ref_fasta"
//   echo "ref_fai=$ref_fai"
//   echo "ref_dict=$ref_dict"
//   echo "intervals=$intervals"
//   echo "scatter_count=$scatter_count"

//   mkdir -p scattered

//   echo "=== split_intervals: BedToIntervalList ==="
//   gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \\
//     -I "$intervals" \\
//     -SD "$ref_dict" \\
//     -O regions.interval_list

//   echo "Produced regions.interval_list:"
//   ls -lah regions.interval_list || true

//   echo "=== split_intervals: SplitIntervals ==="
//   gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \\
//     -R "$ref_fasta" \\
//     -L regions.interval_list \\
//     --scatter "$scatter_count" \\
//     -O scattered

//   echo "=== split_intervals: OUTPUT DIR LIST ==="
//   ls -lah scattered || true
//   echo "Count scattered files:"
//   ls -1 scattered/*.intervals 2>/dev/null | wc -l || true
//   echo "=== split_intervals: END ==="
//   """
// }


// /*
//  * --------------------------------------------
//  * mutect_wrapper
//  *   - Runs Mutect2 for ONE interval shard
//  *   - IMPORTANT: write shard-specific filenames
//  *     so outputs don't collide across shards
//  * --------------------------------------------
//  */
// process mutect_wrapper {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     tuple path(tumor_bam),
//           path(tumor_bam_index),
//           path(interval_shard),
//           path(ref_fasta),
//           path(ref_fai),
//           path(ref_dict),
//           path(germline_resource)
//     val extra_args

//   output:
//     path "*.vcf.gz",      emit: vcf
//     path "*.vcf.gz.tbi",  emit: tbi
//     path "*.stats",       optional: true, emit: stats
//     path "*.f1r2.tar.gz", optional: true, emit: f1r2
//     path "versions.yml", optional: true, emit: versions

//   script:
//   def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
//   def heap_mb   = Math.min(avail_mem, 24000)   // cap at 24 GB

//   """
//   set -euo pipefail

//   echo "=== mutect_wrapper: START ==="
//   echo "PWD=\$(pwd)"
//   echo "Task memory (avail_mem)=${avail_mem}M ; heap_mb=${heap_mb}M"
//   echo "extra_args='${extra_args}'"
//   echo "Inputs present in workdir:"
//   ls -lah

//   echo "interval_shard path: $interval_shard"
//   echo "interval_shard basename: \$(basename "$interval_shard")"
//   echo "tumor_bam: $tumor_bam"
//   echo "ref_fasta: $ref_fasta"
//   echo "germline_resource: $germline_resource"

//   echo "=== mutect_wrapper: verify interval exists ==="
//   ls -lah "$interval_shard" || { echo "ERROR: interval_shard missing" >&2; exit 2; }

//   echo "=== mutect_wrapper: detect tumor sample from BAM header ==="
//   tumor_sample=\$(samtools view -H "$tumor_bam" \\
//     | awk -F'\\t' '/^@RG/ { for (i=1;i<=NF;i++) if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i } }' \\
//     | sort -u)

//   if [ -z "\$tumor_sample" ]; then
//     echo "ERROR: No SM tag found in BAM header" >&2
//     exit 1
//   fi
//   if [ \$(echo "\$tumor_sample" | wc -l) -ne 1 ]; then
//     echo "ERROR: Multiple SM values found in BAM header:" >&2
//     echo "\$tumor_sample" >&2
//     exit 1
//   fi
//   echo "Detected tumor sample: \$tumor_sample"

//   echo "=== mutect_wrapper: ensure germline resource indexed ==="
//   if [ ! -f "${germline_resource}.tbi" ]; then
//     echo "No .tbi found; indexing germline_resource..."
//     gatk IndexFeatureFile -F "$germline_resource"
//   else
//     echo "Found existing index: ${germline_resource}.tbi"
//   fi

//   # shard-specific output prefix (prevents collisions)
//   shard_base=\$(basename "$interval_shard" .intervals)
//   out_prefix="out.\${shard_base}"

//   echo "=== mutect_wrapper: RUN Mutect2 ==="
//   echo "Output prefix: \$out_prefix"
//   echo "Command (intervals): --intervals $interval_shard"

//   gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \\
//     --input "$tumor_bam" \\
//     --reference "$ref_fasta" \\
//     --germline-resource "$germline_resource" \\
//     --intervals "$interval_shard" \\
//     --tmp-dir . \\
//     --tumor-sample "\$tumor_sample" \\
//     ${extra_args} \\
//     --output "\${out_prefix}.vcf.gz"

//   echo "=== mutect_wrapper: END-OF-TASK FILE LIST ==="
//   ls -lah
//   echo "Expected key outputs:"
//   ls -lah "\${out_prefix}.vcf.gz" "\${out_prefix}.vcf.gz.tbi" 2>/dev/null || true

//   echo "=== mutect_wrapper: versions.yml ==="
//   ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
//   cat versions.yml || true


//   echo "=== mutect_wrapper: END ==="
//   """

//   stub:
//   """
//   set -euo pipefail
//   shard_base=\$(basename "$interval_shard" .intervals)
//   out_prefix="out.\${shard_base}"

//   touch "\${out_prefix}.vcf.gz"
//   touch "\${out_prefix}.vcf.gz.tbi"
//   touch "\${out_prefix}.vcf.gz.stats"
//   touch "\${out_prefix}.f1r2.tar.gz"

//   cat <<-END_VERSIONS > versions.yml
//   "${task.process}":
//       gatk4: stub
//   END_VERSIONS
//   """
// }


// /*
//  * --------------------------------------------
//  * mutect_wrapper_force
//  *   - Runs Mutect2 for ONE interval shard with force calling
//  *   - IMPORTANT: write shard-specific filenames
//  *     so outputs don't collide across shards
//  * --------------------------------------------
//  */
// process mutect_wrapper_force {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     tuple path(tumor_bam),
//           path(tumor_bam_index),
//           path(interval_shard),
//           path(ref_fasta),
//           path(ref_fai),
//           path(ref_dict),
//           path(germline_resource),
//           path(force_call_file),
//           path(force_call_file_index)
//     val extra_args

//   output:
//     path "*.vcf.gz",      emit: vcf
//     path "*.vcf.gz.tbi",  emit: tbi
//     path "*.stats",       optional: true, emit: stats
//     path "*.f1r2.tar.gz", optional: true, emit: f1r2
//     path "versions.yml", optional: true, emit: versions

//   script:
//   def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
//   def heap_mb   = Math.min(avail_mem, 24000)   // cap at 24 GB

//   """
//   set -euo pipefail

//   echo "=== mutect_wrapper: START ==="
//   echo "PWD=\$(pwd)"
//   echo "Task memory (avail_mem)=${avail_mem}M ; heap_mb=${heap_mb}M"
//   echo "extra_args='${extra_args}'"
//   echo "Inputs present in workdir:"
//   ls -lah

//   echo "interval_shard path: $interval_shard"
//   echo "interval_shard basename: \$(basename "$interval_shard")"
//   echo "tumor_bam: $tumor_bam"
//   echo "ref_fasta: $ref_fasta"
//   echo "germline_resource: $germline_resource"

//   echo "=== mutect_wrapper: verify interval exists ==="
//   ls -lah "$interval_shard" || { echo "ERROR: interval_shard missing" >&2; exit 2; }

//   echo "=== mutect_wrapper: detect tumor sample from BAM header ==="
//   tumor_sample=\$(samtools view -H "$tumor_bam" \\
//     | awk -F'\\t' '/^@RG/ { for (i=1;i<=NF;i++) if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i } }' \\
//     | sort -u)

//   if [ -z "\$tumor_sample" ]; then
//     echo "ERROR: No SM tag found in BAM header" >&2
//     exit 1
//   fi
//   if [ \$(echo "\$tumor_sample" | wc -l) -ne 1 ]; then
//     echo "ERROR: Multiple SM values found in BAM header:" >&2
//     echo "\$tumor_sample" >&2
//     exit 1
//   fi
//   echo "Detected tumor sample: \$tumor_sample"

//   echo "=== mutect_wrapper: ensure germline resource indexed ==="
//   if [ ! -f "${germline_resource}.tbi" ]; then
//     echo "No .tbi found; indexing germline_resource..."
//     gatk IndexFeatureFile -F "$germline_resource"
//   else
//     echo "Found existing index: ${germline_resource}.tbi"
//   fi

//   # shard-specific output prefix (prevents collisions)
//   shard_base=\$(basename "$interval_shard" .intervals)
//   out_prefix="out.\${shard_base}"

//   echo "=== mutect_wrapper: RUN Mutect2 ==="
//   echo "Output prefix: \$out_prefix"
//   echo "Command (intervals): --intervals $interval_shard"

//   gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \\
//     --input "$tumor_bam" \\
//     --reference "$ref_fasta" \\
//     --germline-resource "$germline_resource" \\
//     --intervals "$interval_shard" \\
//     --tmp-dir . \\
//     --tumor-sample "\$tumor_sample" \\
//     --alleles "$force_call_file" \\
//     ${extra_args} \\
//     --output "\${out_prefix}.vcf.gz"

//   echo "=== mutect_wrapper: END-OF-TASK FILE LIST ==="
//   ls -lah
//   echo "Expected key outputs:"
//   ls -lah "\${out_prefix}.vcf.gz" "\${out_prefix}.vcf.gz.tbi" 2>/dev/null || true

//   echo "=== mutect_wrapper: versions.yml ==="
//   ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
//   cat versions.yml || true


//   echo "=== mutect_wrapper: END ==="
//   """

//   stub:
//   """
//   set -euo pipefail
//   shard_base=\$(basename "$interval_shard" .intervals)
//   out_prefix="out.\${shard_base}"

//   touch "\${out_prefix}.vcf.gz"
//   touch "\${out_prefix}.vcf.gz.tbi"
//   touch "\${out_prefix}.vcf.gz.stats"
//   touch "\${out_prefix}.f1r2.tar.gz"

//   cat <<-END_VERSIONS > versions.yml
//   "${task.process}":
//       gatk4: stub
//   END_VERSIONS
//   """
// }


// process gather_vcfs {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path vcfs, arity: '1..*'

//   output:
//     path "started.txt", optional: true
//     path "merged.vcf.gz"
//     path "merged.vcf.gz.tbi"

//   script:
//   """
//   set -euo pipefail

//   echo "SCRIPT_STARTED \$(date)" > started.txt
//   echo "PWD=\$(pwd)"
//   echo "Listing initial workdir:"
//   ls -lah
//   echo "======================================="

//   echo "Discovering staged VCF files..."
//   find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list

//   echo "Number of VCFs found:"
//   wc -l vcfs.list
//   echo "First few:"
//   head vcfs.list
//   echo "Last few:"
//   tail vcfs.list

//   if awk '{ if (\$0 !~ /out\\.[0-9]+/) { bad=1; print "BAD:", \$0 > "/dev/stderr" } } END{ exit bad }' vcfs.list; then
//     echo "All filenames match expected pattern."
//   else
//     echo "ERROR: Some VCF filenames do not match out.<num> pattern." >&2
//     exit 2
//   fi

//   echo "Sorting VCFs by shard number..."
//   sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list \
//     | sort -k1,1n \
//     | cut -f2- > vcfs.sorted.list

//   echo "Sorted list preview:"
//   head vcfs.sorted.list
//   tail vcfs.sorted.list

//   echo "Building GATK argument file..."
//   awk '{print "-I=" \$1}' vcfs.sorted.list > gather.args
//   echo "Argument preview:"
//   head gather.args
//   tail gather.args

//   echo "Running GATK GatherVcfs at \$(date)"
//   time gatk GatherVcfs --arguments_file gather.args -O merged.vcf.gz
//   echo "Gather finished at \$(date)"

//   if [ ! -s merged.vcf.gz.tbi ]; then
//     echo "Index missing, creating..."
//     gatk IndexFeatureFile -I merged.vcf.gz || tabix -p vcf merged.vcf.gz
//   fi

//   test -s merged.vcf.gz.tbi

//   echo "Final outputs:"
//   ls -lah merged.vcf.gz merged.vcf.gz.tbi
//   echo "=== gather_vcfs COMPLETE ==="
//   """
// }



// /*
//  * --------------------------------------------
//  * workflow
//  * --------------------------------------------
//  */
// workflow {

//   // Build shards channel (each item is a single scattered/*.intervals file)
//   shards_ch = split_intervals(
//     file(params.ref_fasta),
//     file(params.ref_fai),
//     file(params.ref_dict),
//     file(params.intervals),
//     params.scatter_count as int
//   ).shards.flatten()

//   def mutect_res

//   // Pair each shard with shared inputs so mutect runs once per shard
//   if(do_force) {
//     mutect_inputs = shards_ch.map { shard ->
//       tuple(
//         file(params.tumor_reads),
//         file(params.tumor_reads_index),
//         shard,
//         file(params.ref_fasta),
//         file(params.ref_fai),
//         file(params.ref_dict),
//         file(params.germline_resource),
//         file(params.force_call_file),
//         file(params.force_call_file_index)
//       )
//     }

//     mutect_res = mutect_wrapper_force(
//       mutect_inputs,
//       params.m2_extra_args
//     )
//   } else {
//     mutect_inputs = shards_ch.map { shard ->
//       tuple(
//         file(params.tumor_reads),
//         file(params.tumor_reads_index),
//         shard,
//         file(params.ref_fasta),
//         file(params.ref_fai),
//         file(params.ref_dict),
//         file(params.germline_resource)
//       )
//     }
//     mutect_res = mutect_wrapper(
//       mutect_inputs,
//       params.m2_extra_args
//     )

//   }


//   mutect_res.vcf.view { "VCF: $it" }


//   gather_vcfs(mutect_res.vcf.collect())
// }

