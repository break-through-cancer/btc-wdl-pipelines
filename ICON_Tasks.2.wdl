### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Shahab Sarmashghi, ssarmash@broadinstitute.org
### Date last updated: July 9, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute

version 1.0

# TASK DEFINITIONS

task CramToBamTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    String samtools_path
  }
    Float output_bam_size = size(input_cram, "GB") / 0.60
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20
  
  command {
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
 }
  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}

# HaplotypeCaller per-sample in GVCF mode
#HC with read-index explicit
task HaplotypeCaller {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    File interval_list
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    Int hc_scatter

    String? gcs_project_for_requester_pays

    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
  }

  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"]) 

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1
  
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float input_size_calc = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20
  Int disk_size = if input_size_calc > 500 then input_size_calc else 500

  String vcf_basename = if make_gvcf then  basename(output_filename, ".gvcf") else basename(output_filename, ".vcf")
  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      description: "a bam file",
      localization_optional: true
    }
    input_bam_index: {
      description: "an index file for the bam input",
      localization_optional: true
    }
  }
  command {
    set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --read-index ~{input_bam_index} \
      -L ~{interval_list} \
      -O ~{output_filename} \
      -contamination ~{default="0" contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
    File bamout = "~{vcf_basename}.bamout.bam"
  }
}


# Merge GVCFs generated per-interval for the same sample
task MergeVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1
  
  command {
  set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}


#Preprocess intervals
task PreprocessIntervals {
    input {
      File? intervals
      File? blacklist_intervals
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      Int? padding
      Int? bin_length
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # Determine output filename
    String filename = select_first([intervals, "wgs"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" PreprocessIntervals \
            ~{"-L " + intervals} \
            ~{"-XL " + blacklist_intervals} \
            --reference ~{ref_fasta} \
            --padding ~{default="0" padding} \
            --bin-length ~{default="25000" bin_length} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.preprocessed.interval_list
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File preprocessed_intervals = "~{base_filename}.preprocessed.interval_list"
    }
}

#CRC with read-index argumnet explicit
task CollectReadCounts {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    File intervals
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String? gcs_project_for_requester_pays

    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
  }

  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"]) 

  Int machine_mem_gb = select_first([mem_gb, 10])
  Int command_mem_gb = machine_mem_gb - 1

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float input_size_calc = ceil(size(input_bam, "GB") + ref_size) + 20
  Int disk_size = if input_size_calc > 500 then input_size_calc else 500

  parameter_meta {
    input_bam: {
      description: "a bam file",
      localization_optional: true
    }
    input_bam_index: {
      description: "an index file for the bam input",
      localization_optional: true
    }
  }
  command {
    set -e
    
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      CollectReadCounts \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --read-index ~{input_bam_index} \
      -L ~{intervals} \
      -O ~{output_filename} \
      --format TSV \
      --interval-merging-rule OVERLAPPING_ONLY \
      --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter \
      --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter \
      --read-filter PairedReadFilter --read-filter FragmentLengthReadFilter --max-fragment-length 1000 \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateDifferentStrandReadFilter \
      --read-filter NotSupplementaryAlignmentReadFilter --read-filter NotSecondaryAlignmentReadFilter \
      --read-filter FirstOfPairReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 \
      --read-filter OverclippedReadFilter --filter-too-short 50 --seconds-between-progress-updates 100 \
      ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_counts = "~{output_filename}"
  }
}


#task for GenerateSM
task GenerateSampleMapFile {
  input{
    # Command parameters
    String tumor_sample_name #changing from string array to single string and sep for tumor and normal. Optional normal
    String normal_sample_name
    String tumor_gvcf_file #Doesn't take a string but a file; sep tumor and normal GVCF and make files instead of string array
    String normal_gvcf_file
    String outfile 
    
    # Runtime parameters
    String docker = "python:latest"
    Int machine_mem_gb = 7
    Int disk_space_gb = 100
    Int preemptible_attempts = 3
  }
    command <<<
    set -oe pipefail
    
    python << CODE
    file_paths = ["~{normal_gvcf_file}","~{tumor_gvcf_file}"]
    sample_names = ["~{normal_sample_name}", "~{tumor_sample_name}"]

    if len(file_paths)!= len(sample_names):
      print("Number of File Paths does not Equal Number of File Names")
      exit(1)

    with open("~{outfile}", "w") as file_with_s3, open("local_~{outfile}", "w") as file_without_s3:
      for i in range(len(file_paths)):
        sample = sample_names[i]
        path = file_paths[i]
        # Write original paths with s3:// for when sample map lines are passed as inputs
        file_with_s3.write(sample + "\t" + path + "\n")
        # Write modified paths without s3:// for when sample map is parsed locally in tasks
        if path.startswith("s3://"):
          path = path[5:]  # Remove 's3://'
        file_without_s3.write(sample + "\t" + path + "\n")

    CODE
    >>>

    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: 3
    }

    output {
        File sample_map = outfile
        File local_sample_map = "local_" + outfile
    }
}


task Filter {
    input {
        File input_vcf
        String? include  #make this an optional input incase the user wishes to specify a different condition to include variants and specify default
        String? exclude
        String? softFilter
        String outfile
        

        String memory = "256MiB"
        Int timeMinutes = 1 + ceil(size(input_vcf, "G"))
        String bcftools_dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    command {
        set -e 
        bcftools \
        filter \
        -i 'SUM(FORMAT/AD[*:1])>=5' \
        ~{"-e " + exclude} \
        ~{"-s " + softFilter} \
        ~{input_vcf} \
        -O z \
        -o ~{outfile}
        bcftools index --tbi ~{outfile}
    }

    output {
        File filtered_outputVcf = outfile
        File filtered_outputVcfIndex = outfile + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: bcftools_dockerImage
    }

    parameter_meta {
        input_vcf: {description: "The VCF file to operate on.", category: "required"}
        include: {description: "Equivalent to the `-i` option.", category: "common"}

        bcftools_dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
    }
}

task iter_isec { #pending -O z working for wach file in the directory
    input {
        File input_vcf
        File input_vcfindex
        File refDir
        String participant_id

        String memory = "256MiB"
        Int timeMinutes = 1 + ceil(size(refDir, "G"))
        String bcftools_dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
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
        memory: memory
        time_minutes: timeMinutes
        docker: bcftools_dockerImage
        disks: "local-disk 2000 HDD"
    }

    output {
        File final_rem = input_vcf
        File final_remindex = input_vcfindex
        Array[File] chr_int = glob('*_var.vcf.gz')
        Array[File] chr_int_indices = glob('*_var.vcf.gz.tbi')
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
            -o "${filename}"
            bcftools index --tbi "${filename}"
        done
    >>>

    output {
        Array[File] het_filteredandIndex = glob('*_het.vcf.g*')
        Array[File] het_filteredonly = glob('*_het.vcf.gz')
        Array[File] het_filteredIndex = glob('*_het.vcf.gz.tbi')
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

    command <<<
        set -e
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
              --outPrefix="${filename}" \
              --numThreads=12
            fi
        done
    >>>

    output{
        Array[File] eagle_phased = glob('*_het_eagle2phased.vcf.gz')
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
        String bcftools_dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    command <<<
        set -e
        for path in ~{sep=' ' phased_files}
        do
           echo "We are on ${path}"
           outfilename=`basename "${path}" "_het_eagle2phased.vcf.gz"`"_het_eagle2phased.vcf.gz.tbi"
           bcftools index --tbi "${path}" -o "${outfilename}"
        done
    >>>

    output {
        Array[File] phased_indices = glob('*_het_eagle2phased.vcf.gz.tbi')
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
        String bcftools_dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
        
    }

    command <<<
        set -e
        samples=`bcftools query -l ~{chr1_het}`
        for path in ~{sep=' ' all_het}
        do
           hetinfofile=`basename "${path}" "_het.vcf.gz"`"_HetInfo.txt"
           chr=`basename "${path}" "_het.vcf.gz" | awk '{ gsub("^.*\chr","") ; print $0 }'`
           echo "We are on chr${chr} for het querying"
           bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%INFO/DP\n' "${path}" > "${hetinfofile}"
           for sample in ${samples}
           do
              bcftools query -s ${sample} -f '[%AD{0}\t%AD{1}\n]' "${path}" > ${sample}.AD.chr${chr}.txt
            done
        done
        for phased in ~{sep=' ' all_eaglephased}
        do
           chr=`basename "${phased}" "_het_eagle2phased.vcf.gz" | awk '{ gsub("^.*\chr","") ; print $0 }'`
           echo "We are on chr${chr} for phased querying"
           for sample in ${samples}
           do
              bcftools query -s ${sample} -f '[%GT\n]' "${phased}" > ${sample}.phasedGT.chr${chr}.1.txt
           done
       done
    >>>

    output {
        Array[File] variant_info = glob('*_HetInfo.txt')
        Array[File] allelic_depths = glob('*.AD.chr*')
        Array[File] phased_hets = glob("*.phasedGT.chr*")
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:bcftools_dockerImage
    }
}

task ac_calc_agg {
    input{
        File ac_rscript #select it as a separate output from the task where it is created. Don't make it a workflow-level output but specifically for this purpose
        File tumor_counts
        String normal_sample_id
        String tumor_sample_id
        Array[File] all_phased
        Array[File] all_allelic_depths
        Array[File] all_hetinfo
        String participant_id
        Boolean? agg_bins
        Int? desired_binsize
        

        String memory = "4 GB"
        Int timeMinutes = 1 + ceil(size(all_phased, "G"))
        String r_dockerImage = "wchukwu/r-docker:latest"
        
    }

    command <<<
        Rscript ~{ac_rscript} ~{participant_id} ~{tumor_counts} ~{agg_bins} ~{desired_binsize} ~{normal_sample_id} ~{tumor_sample_id} ~{sep=' ' all_phased} ~{sep=' ' all_allelic_depths} ~{sep=' ' all_hetinfo}
    >>>

    output {
        File read_depth = "~{participant_id}_read_depth.txt"
        File allelic_counts = "~{participant_id}_allelic_counts.txt"
        File allelicCNplots = "~{participant_id}_plots.pdf"
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:r_dockerImage
    }
}

task ac_calc_noagg {
    input{
        File ac_rscript #select it as a separate output from the task where it is created. Don't make it a workflow-level output but specifically for this purpose
        File tumor_counts
        String normal_sample_id
        String tumor_sample_id
        
        Array[File] all_phased
        Array[File] all_allelic_depths
        Array[File] all_hetinfo
        String participant_id
        Boolean? agg_bins
        

        String memory = "4 GB"
        Int timeMinutes = 1 + ceil(size(all_phased, "G"))
        String r_dockerImage = "wchukwu/r-docker:latest"
        
    }

    command <<<
        Rscript ~{ac_rscript} ~{participant_id} ~{tumor_counts} ~{agg_bins} ~{normal_sample_id} ~{tumor_sample_id} ~{sep=' ' all_phased} ~{sep=' ' all_allelic_depths} ~{sep=' ' all_hetinfo}
    >>>

    output {
        File read_depth = "~{participant_id}_read_depth.txt"
        File allelic_counts = "~{participant_id}_allelic_counts.txt"
        File allelicCNplots = "~{participant_id}_plots.pdf"
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:r_dockerImage
    }
}

task tangent_XY {
    input{
        File tangent_rscript
        File pon #panel of normals
        File tumor_counts
        File filt_probes
        File ref_clinical
        File? cohort_clinical
        String sample_sex

        #optional inputs
        String? scale_method

        

        String memory = "10 GB"
        Int timeMinutes = 1 + ceil(size(pon, "G"))
        String r_dockerImage = "wchukwu/r-docker:latest"
        
    }

    command <<<
        Rscript ~{tangent_rscript} --pon ~{pon} --tumor_cts ~{tumor_counts} --prbs ~{filt_probes} --sample_sex ~{sample_sex} --ref_cli ~{ref_clinical} --clin ~{cohort_clinical} --scale-mthd ~{scale_method}
    >>>

    output {
        File tangent_norm ='posttangent_log2.rds'
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:r_dockerImage
    }
}

task segment {
    input{
        File segment_rscript
        File allelic_counts
        File tangent_normalized
        String participant_id
        String sample_sex
        Boolean? opt_plts

        #optional arguments for segmentation
        Int? smooth_region 
        Int? outlier_SD
        Int? smooth_SD
        Float? alpha
        Int? min_width
        String? undo_split
        Float? undo_prune
        Int? undo_SD

        

        String memory = "10 GB"
        Int timeMinutes = 1 + ceil(size(allelic_counts, "G"))
        String r_dockerImage = "wchukwu/r-docker:latest"
        
    }

    command <<<
        Rscript ~{segment_rscript} --al_cts ~{allelic_counts} --tangent_cts ~{tangent_normalized} --sample_sex ~{sample_sex} --participant_id ~{participant_id} --smth_reg ~{smooth_region} --outlier_SD ~{outlier_SD} --smooth_SD ~{smooth_SD} --alpha ~{alpha} --min_width ~{min_width} --undo_split ~{undo_split} --undo_prune ~{undo_prune} --undo_SD ~{undo_SD} --add_plots ~{opt_plts}
    >>>

    output {
        #File CN_plots = "~{participant_id}_CN.Segment.2Dplots.pdf"
        File seg_file = "~{participant_id}.seg.txt"
        File ac_cts = "~{participant_id}_processed_counts.txt"
        File segment_plots = "~{participant_id}_segmentPlots.pdf"
        File segment_cts = "~{participant_id}_segment_cts.txt"
        File? outlier_plots = "~{participant_id}_outlierPlots.pdf"
    }

    runtime{
        memory:memory
        time_minutes:timeMinutes
        docker:r_dockerImage
    }
}