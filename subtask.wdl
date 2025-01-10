version 1.0

workflow HapCNA {
  input {
    String participant_id
    File filtered_outputVcf
    File filtered_outputVcfIndex
    File ref_dir1000G
  }

  call iter_isec {
    input:
      refDir = ref_dir1000G,
      input_vcf = filtered_outputVcf,
      participant_id = participant_id,
      input_vcfindex = filtered_outputVcfIndex
  }

  output{
    # String message = iter_isec.message
    File final_rem = iter_isec.final_rem
    File final_remindex = iter_isec.final_remindex
    # Array[File] chr_int = iter_isec.chr_int
    # Array[File] chr_int_indices = iter_isec.chr_int_indices
  }

  meta {
    allowNestedInputs: true
  }
}

task iter_isec { #pending -O z working for wach file in the directory
  input {
    File input_vcf
    File input_vcfindex
    File refDir
    String participant_id

    # String memory = "256MiB"
    # Int timeMinutes = 1 + ceil(size(refDir, "G"))
    # String bcftools_dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
  }

  String dir_name = "hg38/1000G"
  String outDir = participant_id + "_tmp"

  command <<<
    set -e
    
    mkdir -p outDir

    tar -xf ~{refDir}
    for i in {1..2} X
    do
      echo "We are on iteration: ${i}"
      bcftools \
      isec \
      ~{input_vcf} \
      ~{dir_name}/chr${i}.1000G.genotypes.bcf \
      -w 1 \
      -O z \
      -p ~{outDir}
      mv ~{outDir}/0002.vcf.gz outDir/~{participant_id}_chr${i}_var.vcf.gz
      mv ~{outDir}/0002.vcf.gz.tbi outDir/~{participant_id}_chr${i}_var.vcf.gz.tbi
      cp -f ~{outDir}/0000.vcf.gz ~{participant_id}_isec.vcf.gz
      cp -f ~{outDir}/0000.vcf.gz.tbi ~{participant_id}_isec.vcf.gz.tbi
    done
  >>>

  runtime {
    memory: "512MiB"
    # time_minutes: timeMinutes
    docker: "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    disks: "local-disk 2000 HDD"
  }

  output {
    File final_rem = "~{participant_id}_isec.vcf.gz"
    File final_remindex = "~{participant_id}_isec.vcf.gz.tbi"
    Array[File] chr_int = glob("outDir/*_var.vcf.gz")
    Array[File] chr_int_indices = glob("outDir/*_var.vcf.gz.tbi")
  }
}