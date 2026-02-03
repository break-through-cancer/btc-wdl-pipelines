params.m2_extra_args = params.m2_extra_args ?: ''

process mutect_wrapper {
    //tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

    input:
    path tumor_bam
    path tumor_bam_index
    path ref_fasta
    path ref_fai
    path ref_dict
    path germline_resource
    path intervals
    val extra_args



    output:
    path "*.vcf.gz", emit: vcf
    path "*.vcf.gz.tbi", emit: tbi
    path "*.stats", emit: stats
    path "*.f1r2.tar.gz", optional: true, emit: f1r2
    path "versions.yml", emit: versions


    script:
    
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    


    """

    set -euo pipefail

    tumor_sample=\$(samtools view -H "$tumor_bam" | awk -F'\\t' '
      /^@RG/ { for (i=1;i<=NF;i++) if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i } }
    ' | sort -u)

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
      echo "Index missing for germline resource; creating with IndexFeatureFile..."
      gatk IndexFeatureFile -I "$germline_resource"
    fi

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" Mutect2 \
        --input $tumor_bam \
        --reference $ref_fasta \
        --germline-resource $germline_resource \
        --intervals $intervals \
        --tmp-dir . \
        --tumor-sample "\$tumor_sample" \
        $extra_args \
        --output ${tumor_bam.baseName}.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
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

workflow {
    mutect_wrapper(
    file(params.tumor_reads),
    file(params.tumor_reads_index),
    file(params.ref_fasta),
    file(params.ref_fai),
    file(params.ref_dict),
    file(params.germline_resource),
    file(params.intervals),
    params.m2_extra_args
  )
}