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
    File final_rem = iter_isec.final_rem
    File final_remindex = iter_isec.final_remindex
    Array[File] chr_int = iter_isec.chr_int
    Array[File] chr_int_indices = iter_isec.chr_int_indices
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
      tar -xf ~{refDir}
      for i in {1..22}
      do
        echo "We are on iteration: ${i}"
        bcftools \
        isec \
        ~{input_vcf} \
        ~{dir_name}/chr${i}.1000G.genotypes.bcf \
        -w 1 \
        -O z \
        -p ~{outDir}
        mv ~{outDir}/0002.vcf.gz ~{participant_id}"_chr${i}_var.vcf.gz"
        mv ~{outDir}/0002.vcf.gz.tbi ~{participant_id}"_chr${i}_var.vcf.gz.tbi"
        cp -f ~{outDir}/0000.vcf.gz ~{input_vcf}
        cp -f ~{outDir}/0000.vcf.gz.tbi ~{input_vcfindex}
      done
    
      echo "We are on iteration: X"
      bcftools \
      isec \
      ~{input_vcf} \
      ~{dir_name}/chrX.1000G.genotypes.bcf \
      -w 1 \
      -O z \
      -p ~{outDir}
      mv ~{outDir}/0002.vcf.gz ~{participant_id}"_chrX_var.vcf.gz"
      mv ~{outDir}/0002.vcf.gz.tbi ~{participant_id}"_chrX_var.vcf.gz.tbi"
      cp -f ~{outDir}/0000.vcf.gz ~{input_vcf}
      cp -f ~{outDir}/0000.vcf.gz.tbi ~{input_vcfindex}
    >>>

    runtime {
      memory: "256MiB"
      # time_minutes: timeMinutes
      docker: "ghcr.io/break-through-cancer/bcftools-with-stat:latest"
      disks: "local-disk 2000 HDD"
    }

    output {
      File final_rem = input_vcf
      File final_remindex = input_vcfindex
      Array[File] chr_int = glob('*_var.vcf.gz')
      Array[File] chr_int_indices = glob('*_var.vcf.gz.tbi')
    }
  }