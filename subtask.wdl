version 1.0

workflow HapCNA {
  input {
    String participant_id
    File filtered_outputVcf
    File filtered_outputVcfIndex
    File ref_dir1000G
    File geneticMap
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

  call eagle_phasing {
    input:
        het_var = hetvarFilter.het_filteredandIndex,
        refDir = ref_dir1000G,
        eagle_gm = geneticMap
  }

  call index {
    input:
        phased_files = eagle_phasing.eagle_phased
  }

  call query {
    input:
        chr1_het = iter_isec.final_rem,
        all_het = hetvarFilter.het_filteredonly,
        all_het_indices = hetvarFilter.het_filteredIndex,
        all_eaglephased = eagle_phasing.eagle_phased,
        all_eagle_indices = index.phased_indices,
        participant_id = participant_id
  }

  output{
    # String message = iter_isec.message
    Array[File] variant_info = query.variant_info
    Array[File] allelic_depths = query.allelic_depths
    Array[File] phased_hets = query.phased_hets
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
      mkdir -p ~{outDir}

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

task eagle_phasing{
    input {
        Array[File] het_var
        File refDir
        File eagle_gm

        String memory = "100 GB"
        Int timeMinutes = 1 + ceil(size(het_var, "G"))
        String eagle_dockerimage = "cbao/eagle:latest"

    }
        String dir_name = "hg38/1000G"
        String outDir = "outDir"


    command <<<
        set -e
        mkdir -p ~{outDir}

        tar -xf ~{refDir}
        for path in ~{sep=' ' het_var}
        do
            ext=`echo "${path##*.}"`
            if [ "${ext}" == "gz" ]
            then
              echo "We are on ${path}"
              filename=`basename "${path}" "_het.vcf.gz"`"_het_eagle2phased"
              chr=`basename "${path}" "_het.vcf.gz" | awk '{ gsub("^.*\chr","") ; print $0 }'`
              eagle \
              --vcfTarget="${path}" \
              --vcfRef=~{dir_name}/chr${chr}.1000G.genotypes.bcf \
              --geneticMapFile=~{eagle_gm} \
              --chrom=chr${chr} \
              --outPrefix="~{outDir}/${filename}" \
              --numThreads=12
            fi
        done
    >>>

    output{
        Array[File] eagle_phased = glob('~{outDir}/*_het_eagle2phased.vcf.gz')
    }

    runtime{
        memory: memory
        time_minutes: timeMinutes
        docker: eagle_dockerimage
        disks: "local-disk 2000 HDD"
    }
}

task index {
    input{
        Array[File] phased_files

        String memory = "256MiB"
        Int timeMinutes = 1 + ceil(size(phased_files, "G"))
        String bcftools_dockerImage = "ghcr.io/break-through-cancer/bcftools-with-stat:latest"
    }

    command <<<
        set -e
        mkdir -p outDir
        for path in ~{sep=' ' phased_files}
        do
           echo "We are on ${path}"
           outfilename=`basename "${path}" "_het_eagle2phased.vcf.gz"`"_het_eagle2phased.vcf.gz.tbi"
           bcftools index --tbi "${path}" -o "outDir/${outfilename}"
        done
    >>>

    output {
        Array[File] phased_indices = glob('outDir/*_het_eagle2phased.vcf.gz.tbi')
    }

    runtime{
        memory: memory
        time_minutes: timeMinutes
        docker: bcftools_dockerImage
    }
}

task query {
    input{
        File chr1_het #select it as a separate output from the task where it is created. Don't make it a workflow-level output but specifically for this purpose
        Array[File] all_het
        Array[File] all_het_indices
        Array[File] all_eaglephased
        Array[File] all_eagle_indices
        String participant_id

        String memory = "4 GB"
        Int timeMinutes = 1 + ceil(size(all_het, "G"))
        String bcftools_dockerImage = "ghcr.io/break-through-cancer/bcftools-with-stat:latest"
        
    }
    String outDir = "outDir"

    command <<<
        set -e
        mkdir -p ~{outDir}
        samples=`bcftools query -l ~{chr1_het}`
        for path in ~{sep=' ' all_het}
        do
           hetinfofile=`basename "${path}" "_het.vcf.gz"`"_HetInfo.txt"
           chr=`basename "${path}" "_het.vcf.gz" | awk '{ gsub("^.*\chr","") ; print $0 }'`
           echo "We are on chr${chr} for het querying"
           bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%INFO/DP\n' "${path}" > "~{outDir}/${hetinfofile}"
           for sample in ${samples}
           do
              bcftools query -s ${sample} -f '[%AD{0}\t%AD{1}\n]' "${path}" > ~{outDir}/${sample}.AD.chr${chr}.txt
            done
        done
        for phased in ~{sep=' ' all_eaglephased}
        do
           chr=`basename "${phased}" "_het_eagle2phased.vcf.gz" | awk '{ gsub("^.*\chr","") ; print $0 }'`
           echo "We are on chr${chr} for phased querying"
           for sample in ${samples}
           do
              bcftools query -s ${sample} -f '[%GT\n]' "${phased}" > ~{outDir}/${sample}.phasedGT.chr${chr}.1.txt
           done
       done
    >>>

    output {
        Array[File] variant_info = glob('~{outDir}/*_HetInfo.txt')
        Array[File] allelic_depths = glob('~{outDir}/*.AD.chr*')
        Array[File] phased_hets = glob("~{outDir}/*.phasedGT.chr*")
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:bcftools_dockerImage
    }
}