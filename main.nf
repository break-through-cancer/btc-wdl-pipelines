params.m2_extra_args = params.m2_extra_args ?: ''
// JAVA_MEM=${avail_mem}
// sample=$(basename "$tumor_bam" .bam)

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
    path "scattered/*.intervals", emit: shards

  script:
  """
  set -eo pipefail
  mkdir -p scattered

  gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \
    -I "$intervals" \
    -SD "$ref_dict" \
    -O regions.interval_list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \
    -R "$ref_fasta" \
    -L regions.interval_list \
    --scatter "$scatter_count" \
    -O scattered

  ls -lah scattered
  """
}



process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple path(tumor_bam), path(tumor_bam_index), path(interval_shard)
    path ref_fasta
    path ref_fai
    path ref_dict
    path germline_resource
    val  extra_args

  output:
    path "*.vcf.gz",      emit: vcf
    path "*.vcf.gz.tbi",  emit: tbi
    path "*.stats",       emit: stats
    path "*.f1r2.tar.gz", optional: true, emit: f1r2
    path "versions.yml",  emit: versions

  script:
  
  def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072

  """
  echo "INTERVAL_SHARD=$interval_shard"
  ls -lah

  set -euo pipefail

  # Get unique SM tag from BAM header (no nested quoting issues)
  tumor_sample=\$(samtools view -H "$tumor_bam" \
    | awk -F'\\t' '/^@RG/ { for (i=1;i<=NF;i++) if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i } }' \
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

  # Ensure germline resource is indexed
  if [ ! -f "${germline_resource}.tbi" ]; then
    gatk IndexFeatureFile -F "$germline_resource"
  fi

  gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" Mutect2 \\
    --input "$tumor_bam" \\
    --reference "$ref_fasta" \\
    --germline-resource "$germline_resource" \\
    --intervals "$interval_shard" \\
    --tmp-dir . \\
    --tumor-sample "\$tumor_sample" \\
    ${extra_args} \\
    --output out.vcf.gz

  # Record versions
  gatk --version > versions.yml 2>&1
  """


  stub:
  """
  touch ${tumor_bam.baseName}.vcf.gz
  touch ${tumor_bam.baseName}.vcf.gz.tbi
  touch ${tumor_bam.baseName}.vcf.gz.stats
  touch ${tumor_bam.baseName}.f1r2.tar.gz

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gatk4: stub
  END_VERSIONS
  """
}


process gather_vcfs {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path vcfs

  output:
    path "merged.vcf.gz"
    path "merged.vcf.gz.tbi"

  script:
  """
  set -euo pipefail
  gatk GatherVcfs \\
    ${vcfs.collect{ "-I ${it}" }.join(' ')} \\
    -O merged.vcf.gz
  """
}


workflow {

  shards_ch = split_intervals(
    file(params.ref_fasta),
    file(params.ref_fai),
    file(params.ref_dict),
    file(params.intervals),
    params.scatter_count as int
  ).shards

  shards_ch = split_intervals(...).shards

  tumor_triplets = shards_ch.map { shard ->
    tuple(
      file(params.tumor_reads),
      file(params.tumor_reads_index),
      shard
    )
  }

  mutect_res = mutect_wrapper(
    tumor_triplets,
    file(params.ref_fasta),
    file(params.ref_fai),
    file(params.ref_dict),
    file(params.germline_resource),
    params.m2_extra_args
  )

  gather_vcfs(mutect_res.vcf.collect())

}