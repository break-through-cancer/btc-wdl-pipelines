params.m2_extra_args = params.m2_extra_args ?: ''

/*
 * --------------------------------------------
 * split_intervals
 *   - BED -> interval_list (needs .dict)
 *   - interval_list -> scattered/*.intervals
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
  echo "Inputs:"
  ls -lah
  echo "ref_fasta=$ref_fasta"
  echo "ref_fai=$ref_fai"
  echo "ref_dict=$ref_dict"
  echo "intervals=$intervals"
  echo "scatter_count=$scatter_count"

  mkdir -p scattered

  echo "=== split_intervals: BedToIntervalList ==="
  gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \\
    -I "$intervals" \\
    -SD "$ref_dict" \\
    -O regions.interval_list

  echo "Produced regions.interval_list:"
  ls -lah regions.interval_list || true

  echo "=== split_intervals: SplitIntervals ==="
  gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \\
    -R "$ref_fasta" \\
    -L regions.interval_list \\
    --scatter "$scatter_count" \\
    -O scattered

  echo "=== split_intervals: OUTPUT DIR LIST ==="
  ls -lah scattered || true
  echo "Count scattered files:"
  ls -1 scattered/*.intervals 2>/dev/null | wc -l || true
  echo "=== split_intervals: END ==="
  """
}


/*
 * --------------------------------------------
 * mutect_wrapper
 *   - Runs Mutect2 for ONE interval shard
 *   - IMPORTANT: write shard-specific filenames
 *     so outputs don't collide across shards
 * --------------------------------------------
 */
process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple path(tumor_bam),
          path(tumor_bam_index),
          path(interval_shard),
          path(ref_fasta),
          path(ref_fai),
          path(ref_dict),
          path(germline_resource)
    val extra_args

  output:
    path "*.vcf.gz",      emit: vcf
    path "*.vcf.gz.tbi",  emit: tbi
    path "*.stats",       optional: true, emit: stats
    path "*.f1r2.tar.gz", optional: true, emit: f1r2
    path "versions.yml", optional: true, emit: versions

  script:
  def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
  def heap_mb   = Math.min(avail_mem, 24000)   // cap at 24 GB

  """
  set -euo pipefail

  echo "=== mutect_wrapper: START ==="
  echo "PWD=\$(pwd)"
  echo "Task memory (avail_mem)=${avail_mem}M ; heap_mb=${heap_mb}M"
  echo "extra_args='${extra_args}'"
  echo "Inputs present in workdir:"
  ls -lah

  echo "interval_shard path: $interval_shard"
  echo "interval_shard basename: \$(basename "$interval_shard")"
  echo "tumor_bam: $tumor_bam"
  echo "ref_fasta: $ref_fasta"
  echo "germline_resource: $germline_resource"

  echo "=== mutect_wrapper: verify interval exists ==="
  ls -lah "$interval_shard" || { echo "ERROR: interval_shard missing" >&2; exit 2; }

  echo "=== mutect_wrapper: detect tumor sample from BAM header ==="
  tumor_sample=\$(samtools view -H "$tumor_bam" \\
    | awk -F'\\t' '/^@RG/ { for (i=1;i<=NF;i++) if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i } }' \\
    | sort -u)

  if [ -z "\$tumor_sample" ]; then
    echo "ERROR: No SM tag found in BAM header" >&2
    exit 1
  fi
  if [ \$(echo "\$tumor_sample" | wc -l) -ne 1 ]; then
    echo "ERROR: Multiple SM values found in BAM header:" >&2
    echo "\$tumor_sample" >&2
    exit 1
  fi
  echo "Detected tumor sample: \$tumor_sample"

  echo "=== mutect_wrapper: ensure germline resource indexed ==="
  if [ ! -f "${germline_resource}.tbi" ]; then
    echo "No .tbi found; indexing germline_resource..."
    gatk IndexFeatureFile -F "$germline_resource"
  else
    echo "Found existing index: ${germline_resource}.tbi"
  fi

  # shard-specific output prefix (prevents collisions)
  shard_base=\$(basename "$interval_shard" .intervals)
  out_prefix="out.\${shard_base}"

  echo "=== mutect_wrapper: RUN Mutect2 ==="
  echo "Output prefix: \$out_prefix"
  echo "Command (intervals): --intervals $interval_shard"

  gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \\
    --input "$tumor_bam" \\
    --reference "$ref_fasta" \\
    --germline-resource "$germline_resource" \\
    --intervals "$interval_shard" \\
    --tmp-dir . \\
    --tumor-sample "\$tumor_sample" \\
    ${extra_args} \\
    --output "\${out_prefix}.vcf.gz"

  echo "=== mutect_wrapper: END-OF-TASK FILE LIST ==="
  ls -lah
  echo "Expected key outputs:"
  ls -lah "\${out_prefix}.vcf.gz" "\${out_prefix}.vcf.gz.tbi" 2>/dev/null || true

  echo "=== mutect_wrapper: versions.yml ==="
  ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
  cat versions.yml || true


  echo "=== mutect_wrapper: END ==="
  """

  stub:
  """
  set -euo pipefail
  shard_base=\$(basename "$interval_shard" .intervals)
  out_prefix="out.\${shard_base}"

  touch "\${out_prefix}.vcf.gz"
  touch "\${out_prefix}.vcf.gz.tbi"
  touch "\${out_prefix}.vcf.gz.stats"
  touch "\${out_prefix}.f1r2.tar.gz"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gatk4: stub
  END_VERSIONS
  """
}


/*
 * --------------------------------------------
 * gather_vcfs
 *   - Collect per-shard VCFs
 * --------------------------------------------
 */
// process gather_vcfs {
//   label 'process_medium'
//   container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

//   input:
//     path vcfs

//   output:
//     path "merged.vcf.gz"
//     path "merged.vcf.gz.tbi"

//   script:
//   """
//   set -euo pipefail
//   echo "=== gather_vcfs: START ==="
//   echo "PWD=\$(pwd)"
//   echo "Inputs:"
//   ls -lah

//   echo "VCF list received by process:"
//   for f in ${vcfs}; do
//     echo " - \$f"
//   done

//   gatk GatherVcfs \\
//     ${vcfs.collect{ "-I ${it}" }.join(' ')} \\
//     -O merged.vcf.gz

//   echo "=== gather_vcfs: outputs ==="
//   ls -lah merged.vcf.gz merged.vcf.gz.tbi || true
//   echo "=== gather_vcfs: END ==="
//   """
// }
process gather_vcfs {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path vcfs, arity: '1..*'

  output:
    path "started.txt", optional: true
    path "merged.vcf.gz"
    path "merged.vcf.gz.tbi"

  script:
  """
  set -euo pipefail

  # === PROVE CONTAINER ACTUALLY STARTED ===
  echo "SCRIPT_STARTED \$(date)" > started.txt
  echo "PWD=\$(pwd)"
  echo "Listing initial workdir:"
  ls -lah
  echo "======================================="

  # === Discover staged VCFs safely ===
  echo "Discovering staged VCF files..."
  find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list

  echo "Number of VCFs found:"
  wc -l vcfs.list
  echo "First few:"
  head vcfs.list
  echo "Last few:"
  tail vcfs.list

  # === Fail fast if filenames don't match expected pattern ===
  if awk '{ if (\$0 !~ /out\\.[0-9]+/) { bad=1; print "BAD:", \$0 > "/dev/stderr" } } END{ exit bad }' vcfs.list; then
    echo "All filenames match expected pattern."
  else
    echo "ERROR: Some VCF filenames do not match out.<num> pattern." >&2
    exit 2
  fi

  # === Sort numerically by shard number ===
  echo "Sorting VCFs by shard number..."
  sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list \
    | sort -k1,1n \
    | cut -f2- > vcfs.sorted.list

  echo "Sorted list preview:"
  head vcfs.sorted.list
  tail vcfs.sorted.list

  # === Build argument file safely ===
  echo "Building GATK argument file..."
  awk '{print "-I="\\\$0}' vcfs.sorted.list > gather.args
  echo "Argument preview:"
  head gather.args
  tail gather.args

  # === Run GatherVcfs ===
  echo "Running GATK GatherVcfs at \$(date)"
  time gatk GatherVcfs --arguments_file gather.args -O merged.vcf.gz
  echo "Gather finished at \$(date)"

  # === Ensure index exists ===
  if [ ! -s merged.vcf.gz.tbi ]; then
    echo "Index missing, creating..."
    gatk IndexFeatureFile -I merged.vcf.gz || tabix -p vcf merged.vcf.gz
  fi

  test -s merged.vcf.gz.tbi

  echo "Final outputs:"
  ls -lah merged.vcf.gz merged.vcf.gz.tbi
  echo "=== gather_vcfs COMPLETE ==="
  """
}


/*
 * --------------------------------------------
 * workflow
 * --------------------------------------------
 */
workflow {

  // Build shards channel (each item is a single scattered/*.intervals file)
  shards_ch = split_intervals(
    file(params.ref_fasta),
    file(params.ref_fai),
    file(params.ref_dict),
    file(params.intervals),
    params.scatter_count as int
  ).shards.flatten()

  // Pair each shard with shared inputs so mutect runs once per shard
  mutect_inputs = shards_ch.map { shard ->
    tuple(
      file(params.tumor_reads),
      file(params.tumor_reads_index),
      shard,
      file(params.ref_fasta),
      file(params.ref_fai),
      file(params.ref_dict),
      file(params.germline_resource)
    )
  }

  mutect_res = mutect_wrapper(
    mutect_inputs,
    params.m2_extra_args
  )

  mutect_res.vcf.view { "VCF: $it" }


  gather_vcfs(mutect_res.vcf.collect())
}


// params.m2_extra_args = params.m2_extra_args ?: ''
// // JAVA_MEM=${avail_mem}
// // sample=$(basename "$tumor_bam" .bam)

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
//     path "scattered/*.intervals", emit: shards

//   script:
//   """
//   set -eo pipefail
//   mkdir -p scattered

//   gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \
//     -I "$intervals" \
//     -SD "$ref_dict" \
//     -O regions.interval_list

//   gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \
//     -R "$ref_fasta" \
//     -L regions.interval_list \
//     --scatter "$scatter_count" \
//     -O scattered

//   ls -lah scattered
//   """
// }



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
//     path "*.stats", optional: true, emit: stats
//     path "*.f1r2.tar.gz", optional: true, emit: f1r2
//     path "versions.yml",  emit: versions

//   script:
  
//   def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
//   def heap_mb   = Math.min(avail_mem, 24000)   // cap at 24 GB
//   """
//   echo "INTERVAL_SHARD=$interval_shard"
//   ls -lah

//   set -euo pipefail

//   # Get unique SM tag from BAM header (no nested quoting issues)
//   tumor_sample=\$(samtools view -H "$tumor_bam" \
//     | awk -F'\\t' '/^@RG/ { for (i=1;i<=NF;i++) if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i } }' \
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

//   # Ensure germline resource is indexed
//   if [ ! -f "${germline_resource}.tbi" ]; then
//     gatk IndexFeatureFile -F "$germline_resource"
//   fi

//   gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \\
//     --input "$tumor_bam" \\
//     --reference "$ref_fasta" \\
//     --germline-resource "$germline_resource" \\
//     --intervals "$interval_shard" \\
//     --tmp-dir . \\
//     --tumor-sample "\$tumor_sample" \\
//     ${extra_args} \\
//     --output out.vcf.gz

//   # Record versions
//   gatk --version > versions.yml 2>&1
//   """


//   stub:
//   """
//   touch ${tumor_bam.baseName}.vcf.gz
//   touch ${tumor_bam.baseName}.vcf.gz.tbi
//   touch ${tumor_bam.baseName}.vcf.gz.stats
//   touch ${tumor_bam.baseName}.f1r2.tar.gz

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
//     path vcfs

//   output:
//     path "merged.vcf.gz"
//     path "merged.vcf.gz.tbi"

//   script:
//   """
//   set -euo pipefail
//   gatk GatherVcfs \\
//     ${vcfs.collect{ "-I ${it}" }.join(' ')} \\
//     -O merged.vcf.gz
//   """
// }


// workflow {

//   shards_ch = split_intervals(
//     file(params.ref_fasta),
//     file(params.ref_fai),
//     file(params.ref_dict),
//     file(params.intervals),
//     params.scatter_count as int
//   ).shards.flatten()

//   mutect_inputs = shards_ch.map { shard ->
//     tuple(
//       file(params.tumor_reads),
//       file(params.tumor_reads_index),
//       shard,
//       file(params.ref_fasta),
//       file(params.ref_fai),
//       file(params.ref_dict),
//       file(params.germline_resource)
//     )
// }


// mutect_res = mutect_wrapper(
//   mutect_inputs,
//   params.m2_extra_args
// )


// gather_vcfs(mutect_res.vcf.collect())

// }