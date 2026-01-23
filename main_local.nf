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
    path gnomad_vcf
    path gnomad_idx
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
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" Mutect2 \
        --input $tumor_bam \
        --reference $ref_fasta \
        --germline-resource $gnomad_vcf \
        --germline-resource-index $gnomad_idx \
        
        --tmp-dir . \
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
    tumor_bam       = file("test_data/tumor.bam")
    tumor_bam_index = file("test_data/tumor.bam.bai")
    ref_fasta       = file("test_data/ref.fa")
    ref_fai         = file("test_data/ref.fa.fai")
    ref_dict        = file("test_data/ref.dict")
    gnomad_vcf      = file("test_data/gnomad.vcf.gz")
    gnomad_idx      = file("test_data/gnomad.vcf.gz.tbi")

    mutect_wrapper(
        tumor_bam,
        tumor_bam_index,
        ref_fasta,
        ref_fai,
        ref_dict,
        gnomad_vcf,
        gnomad_idx,
        params.extra_args ?: ''
    )
}







