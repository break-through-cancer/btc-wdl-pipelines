params.m2_extra_args = params.m2_extra_args ?: ''
// JAVA_MEM=${avail_mem}
// sample=$(basename "$tumor_bam" .bam)

process split_intervals {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    path ref_fasta
    path intervals
    val scatter_count

  output:
    path "scattered/*.interval_list", emit: shards

  script:
  """
  set -euo pipefail
  mkdir -p scattered

  gatk SplitIntervals \
    -R ${ref_fasta} \
    -L ${intervals} \
    --scatter ${scatter_count} \
    -O scattered
  """
}


script {
  def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072

  """
  set -euo pipefail

  export AVAIL_MEM=${avail_mem}
  export EXTRA_ARGS="${extra_args}"

  bash -lc '
    set -euo pipefail

    tumor_sample=$(samtools view -H "$tumor_bam" | awk -F"\\t" "
      /^@RG/ {
        for (i=1;i<=NF;i++)
          if (\\$i ~ /^SM:/) { sub(/^SM:/,\\\"\\\",\\$i); print \\$i }
      }
    " | sort -u)

    if [ -z "$tumor_sample" ]; then
      echo "ERROR: No SM tag found in BAM header" >&2
      exit 1
    fi

    if [ $(echo "$tumor_sample" | wc -l) -ne 1 ]; then
      echo "ERROR: Multiple SM values found in BAM header:" >&2
      echo "$tumor_sample" >&2
      exit 1
    fi

    echo "Detected tumor sample: $tumor_sample"

    if [ ! -f "${germline_resource}.tbi" ]; then
      echo "Index missing for germline resource; creating with IndexFeatureFile..."
      gatk IndexFeatureFile -F "$germline_resource"
    fi

    shard_id=$(basename "$interval_shard" | sed "s/\\.interval_list$//")
    sample=$(basename "$tumor_bam" .bam)

    # run mutect
    gatk --java-options "-Xmx${AVAIL_MEM}M -XX:-UsePerfData" Mutect2 \
      --input "$tumor_bam" \
      --reference "$ref_fasta" \
      --germline-resource "$germline_resource" \
      --intervals "$interval_shard" \
      --tmp-dir . \
      --tumor-sample "$tumor_sample" \
      $EXTRA_ARGS \
      --output "${sample}.${shard_id}.vcf.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: $(echo $(gatk --version 2>&1) | sed "s/^.*(GATK) v//; s/ .*\\$//")
    END_VERSIONS
  '
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
    file(params.intervals),
    params.scatter_count as int
  ).shards

  mutect_res = mutect_wrapper(
    file(params.tumor_reads),
    file(params.tumor_reads_index),
    file(params.ref_fasta),
    file(params.ref_fai),
    file(params.ref_dict),
    file(params.germline_resource),
    shards_ch,
    params.m2_extra_args
  )

  gather_vcfs(mutect_res.vcf.collect())

}