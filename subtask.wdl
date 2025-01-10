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

  call hetvarFilter {
    input:
        chr_int = iter_isec.chr_int,
        chr_int_indices = iter_isec.chr_int_indices,
        participant_id = participant_id
  }

  output{
    # String message = iter_isec.message
    Array[File] final_rem = hetvarFilter.het_filteredandIndex
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
  String outFile = participant_id + "_isec.vcf.gz"

  command <<<
    set -e
    
    mkdir -p ~{outDir}/chr_outputs

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
      mv ~{outDir}/0002.vcf.gz ~{outDir}/chr_outputs/~{participant_id}_chr${i}_var.vcf.gz
      mv ~{outDir}/0002.vcf.gz.tbi ~{outDir}/chr_outputs/~{participant_id}_chr${i}_var.vcf.gz.tbi
      cp -f ~{outDir}/0000.vcf.gz ~{outFile}
      cp -f ~{outDir}/0000.vcf.gz.tbi ~{outFile}.tbi
    done
  >>>

  runtime {
    memory: "512MiB"
    # time_minutes: timeMinutes
    docker: "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    disks: "local-disk 2000 HDD"
  }

  output {
    File final_rem = outFile
    File final_remindex = outFile + ".tbi"
    Array[File] chr_int = glob("~{outDir}/chr_outputs/*_var.vcf.gz")
    Array[File] chr_int_indices = glob("~{outDir}/chr_outputs/*_var.vcf.gz.tbi")
  }
}

task hetvarFilter { #this task also assumes that the first sample is germline;take out the refDir extraction
    input {
        Array[File] chr_int
        Array[File] chr_int_indices
        String? include
        String? exclude
        String? softFilter
        String participant_id
        #File refDir

        String memory = "256MiB"
        Int timeMinutes = 1 + ceil(size(chr_int, "G"))
        String bcftools_dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    #String outname = basename(...)
    String outDir = participant_id + "_tmp"

    command <<<
        set -e
        for chr in ~{sep=' ' chr_int}
        do
            echo "We are on iteration: ${chr}"
            filename=`basename "${chr}" "_var.vcf.gz"`"_het.vcf.gz"
            bcftools \
            filter \
            -i 'GT[1]="het" & FORMAT/AD[1:0]>=5 & FORMAT/AD[1:1]>=5' \
            ~{"-e" + exclude} \
            ~{"-s" + softFilter} \
            "${chr}" \
            -O z \
            -o "~{outDir}/${filename}"
            bcftools index --tbi "~{outDir}/${filename}"
        done
    >>>

    output {
        Array[File] het_filteredandIndex = glob('~{outDir}/*_het.vcf.g*')
        Array[File] het_filteredonly = glob('~{outDir}/*_het.vcf.gz')
        Array[File] het_filteredIndex = glob('~{outDir}/*_het.vcf.gz.tbi')
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: bcftools_dockerImage
    }
}