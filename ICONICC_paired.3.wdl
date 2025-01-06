### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Shahab Sarmashghi, ssarmash@broadinstitute.org
### Date last updated: June 24, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute

version 1.0

## Pipeline for Tumor-Normal analyses

import "https://raw.githubusercontent.com/break-through-cancer/btc-wdl-pipelines/refs/heads/iconicc/JointGeno_Tasks.1.wdl" as JGTasks
import "https://raw.githubusercontent.com/break-through-cancer/btc-wdl-pipelines/refs/heads/iconicc/ICON_Tasks.2.wdl" as ICTasks

# WORKFLOW DEFINITION 
workflow HapCNA {
  input {
    Int sample_num_threshold
    File tumor_input_bam
    File tumor_input_bam_index
    File normal_input_bam
    File normal_input_bam_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File scattered_calling_intervals_list
  
    Boolean make_gvcf = true
    Boolean make_bamout = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gatk_path = "/gatk/gatk"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String samtools_path = "samtools"

    ###################################################
    #### optional arguments for PreprocessIntervals ####
    ####################################################
    Int? padding
    Int? readcts_binlength
    Int? mem_gb_for_preprocess_intervals
    File? blacklist_intervals
    File? intervals
    File? gatk4_jar_override
    Int? preemptible_attempts

    #INPUTS FOR GENERATE SAMPLE MAP
    String tumor_sample_id #IDs for your tumor and normal samples
    String normal_sample_id 
    String participant_id
      
    #INPUTS FOR JOINT GENOTYPING
    File unpadded_intervals_file
    File dbsnp_vcf
    File dbsnp_vcf_index
    Int small_disk = 100
    Int medium_disk = 200
    Int large_disk = 1000
    Int huge_disk = 2000
    Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]
    Array[String] snp_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]
    Array[String] indel_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]
    Array[String] indel_recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]
    File haplotype_database
    File eval_interval_list
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    Float excess_het_threshold = 54.69
    Float snp_filter_level = 99.7
    Float indel_filter_level = 99.0
    Int SNP_VQSR_downsampleFactor = 10
    Int? top_level_scatter_count
    Boolean? gather_vcfs
    Int snps_variant_recalibration_threshold = 20000
    Boolean rename_gvcf_samples = true
    Float unbounded_scatter_count_scale_factor = 2.5
    Int gnarly_scatter_count = 10
    Boolean use_gnarly_genotyper = false
    Boolean use_allele_specific_annotations = false
    Boolean cross_check_fingerprints = true
    Boolean scatter_cross_check_fingerprints = false 

    #INPUTS FOR FILTERING VARIANTS AND PHASING
    File ref_dir1000G
    File geneticMap
    #optional entries for variant filtering
    String? include
    String? exclude
    String? softFilter

    #INPUTS FOR CALCULATING ALLELIC FRACTIONS
    File ac_rscript
    Boolean? agg_bins = false
    Int? desired_binsize

    #ADDITIONAL INPUTS FOR TANGENT
    File tangent_rscript
    File pon
    File filt_probes
    File ref_clinical
    File? cohort_clinical
    String? scale_method = "median"
    String sample_sex

    #ADDITIONAL INPUTS FOR SEGMENTATION
    File segment_rscript
    Boolean? opt_plts = false
    Int? smooth_region = 10 
    Int? outlier_SD = 4
    Int? smooth_SD = 2
    Float? alpha = 0.01
    Int? min_width = 2
    String? undo_split = "none"
    Float? undo_prune = 0.05
    Int? undo_SD = 2
  
  }  

    Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)
    
    Boolean is_normal_cram = sub(basename(normal_input_bam), ".*\\.", "") == "cram"
    Boolean is_tumor_cram = sub(basename(tumor_input_bam), ".*\\.", "") == "cram"
    Boolean is_normal_bam = sub(basename(normal_input_bam), ".*\\.", "") == "bam"
    Boolean is_tumor_bam = sub(basename(tumor_input_bam), ".*\\.", "") == "bam"

    String normal_sample_basename = if is_normal_cram then  basename(normal_input_bam, ".cram") else basename(normal_input_bam, ".bam")
    String tumor_sample_basename = if is_tumor_cram then  basename(tumor_input_bam, ".cram") else basename(tumor_input_bam, ".bam")
    String normal_vcf_basename = normal_sample_basename
    String tumor_vcf_basename = tumor_sample_basename
    String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String normal_output_filename = normal_sample_id + output_suffix
    String tumor_output_filename = tumor_sample_id + output_suffix

    Int potential_hc_divisor = length(scattered_calling_intervals) - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1
    

    #ADDITIONAL INPUTS FOR COLLECTREADCOUNT
    String normal_count_output_filename = normal_sample_id + "_counts.txt"
    String tumor_count_output_filename = tumor_sample_id + "_counts.txt"
    
    #INPUTS FOR PREPROCESS INTERVALS
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fasta_index, "GB"))
    Int preprocess_intervals_disk = ref_size + 20
  
  if ( is_tumor_cram ) {
    call ICTasks.CramToBamTask as TumorCramToBam{
      input:
        input_cram = tumor_input_bam,
        sample_name = tumor_sample_id,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gitc_docker,
        samtools_path = samtools_path
    }
  }

  if ( is_normal_cram ) {
    call ICTasks.CramToBamTask as NormalCramToBam{
      input:
        input_cram = normal_input_bam,
        sample_name = normal_sample_id,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gitc_docker,
        samtools_path = samtools_path
    }
  }
  # Call variants in parallel over grouped calling intervals
  
  scatter (interval_file in scattered_calling_intervals) {
    call ICTasks.HaplotypeCaller as NormalHC {
      input:
        input_bam = select_first([NormalCramToBam.output_bam, normal_input_bam]),
        input_bam_index = select_first([NormalCramToBam.output_bai, normal_input_bam_index]),
        interval_list = interval_file,
        output_filename = normal_output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_gvcf = make_gvcf,
        make_bamout = make_bamout,
        docker = gatk_docker,
        gatk_path = gatk_path
      }
    }
  
  scatter (interval_file in scattered_calling_intervals) {
    call ICTasks.HaplotypeCaller as TumorHC {
      input:
        input_bam = select_first([TumorCramToBam.output_bam, tumor_input_bam]),
        input_bam_index = select_first([TumorCramToBam.output_bai, tumor_input_bam_index]),
        interval_list = interval_file,
        output_filename = tumor_output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_gvcf = make_gvcf,
        make_bamout = make_bamout,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }
  

  call ICTasks.MergeVCFs as NormalMerge {
    input:
      input_vcfs = NormalHC.output_vcf,
      input_vcfs_indexes = NormalHC.output_vcf_index,
      output_filename = normal_output_filename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  call ICTasks.MergeVCFs as TumorMerge {
    input:
      input_vcfs = TumorHC.output_vcf,
      input_vcfs_indexes = TumorHC.output_vcf_index,
      output_filename = tumor_output_filename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  call ICTasks.GenerateSampleMapFile {
    input:
      tumor_sample_name = tumor_sample_id,
      normal_sample_name = normal_sample_id,
      tumor_gvcf_file = TumorMerge.output_vcf,
      normal_gvcf_file = NormalMerge.output_vcf,
      outfile = participant_id + ".sample_map"
  }
    Boolean allele_specific_annotations = !use_gnarly_genotyper && use_allele_specific_annotations

    Array[Array[String]] sample_name_map_lines = read_tsv(GenerateSampleMapFile.sample_map)
    Int num_gvcfs = length(sample_name_map_lines)
    Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])
    Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
    Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2

  
  call JGTasks.CheckSamplesUnique {
    input:
      sample_name_map = GenerateSampleMapFile.sample_map,
      sample_num_threshold = sample_num_threshold
  }

  call JGTasks.SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk,
      sample_names_unique_done = CheckSamplesUnique.samples_unique
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call JGTasks.ImportGVCFs {
      input:
        sample_name_map = GenerateSampleMapFile.local_sample_map,
        tumor_gvcf_file = TumorMerge.output_vcf,
        normal_gvcf_file = NormalMerge.output_vcf,
        tumor_gvcf_index = TumorMerge.output_vcf_index,
        normal_gvcf_index = NormalMerge.output_vcf_index,

        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        disk_size = medium_disk,
        batch_size = 50
    }

    if (use_gnarly_genotyper) {

      call JGTasks.SplitIntervalList as GnarlyIntervalScatterDude {
        input:
          interval_list = unpadded_intervals[idx],
          scatter_count = gnarly_scatter_count,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size = small_disk,
          sample_names_unique_done = CheckSamplesUnique.samples_unique
      }

      Array[File] gnarly_intervals = GnarlyIntervalScatterDude.output_intervals

      scatter (gnarly_idx in range(length(gnarly_intervals))) {
        call JGTasks.GnarlyGenotyper {
          input:
            workspace_tar = ImportGVCFs.output_genomicsdb,
            interval = gnarly_intervals[gnarly_idx],
            output_vcf_filename = participant_id + "." + idx + "." + gnarly_idx + ".vcf.gz",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp_vcf = dbsnp_vcf,
        }
      }

      Array[File] gnarly_gvcfs = GnarlyGenotyper.output_vcf

      call JGTasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcfs = gnarly_gvcfs,
          output_vcf_name = participant_id + "." + idx + ".gnarly.vcf.gz",
          disk_size = large_disk
      }
    }

    if (!use_gnarly_genotyper) {
      call JGTasks.GenotypeGVCFs {
        input:
          workspace_tar = ImportGVCFs.output_genomicsdb,
          interval = unpadded_intervals[idx],
          output_vcf_filename = participant_id + "." + idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          disk_size = medium_disk
      }
    }

    File genotyped_vcf = select_first([TotallyRadicalGatherVcfs.output_vcf, GenotypeGVCFs.output_vcf])
    File genotyped_vcf_index = select_first([TotallyRadicalGatherVcfs.output_vcf_index, GenotypeGVCFs.output_vcf_index])

    call JGTasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = genotyped_vcf,
        vcf_index = genotyped_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = participant_id + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = participant_id + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size = medium_disk
    }
  }

  call JGTasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      output_vcf_name = participant_id + ".sites_only.vcf.gz",
      disk_size = medium_disk
  }

  call JGTasks.IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = participant_id + ".indels.recal",
      tranches_filename = participant_id + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_vcf,
      dbsnp_resource_vcf_index = dbsnp_vcf_index,
      use_allele_specific_annotations = allele_specific_annotations,
      disk_size = small_disk
  }

  if (num_gvcfs > snps_variant_recalibration_threshold) {
    call JGTasks.SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = participant_id + ".snps.recal",
        tranches_filename = participant_id + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        downsampleFactor = SNP_VQSR_downsampleFactor,
        model_report_filename = participant_id + ".snps.model.report",
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_vcf,
        dbsnp_resource_vcf_index = dbsnp_vcf_index,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size = small_disk
    }

    scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
      call JGTasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
          sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
          recalibration_filename = participant_id + ".snps." + idx + ".recal",
          tranches_filename = participant_id + ".snps." + idx + ".tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_vcf,
          dbsnp_resource_vcf_index = dbsnp_vcf_index,
          use_allele_specific_annotations = allele_specific_annotations,
          disk_size = small_disk
        }
    }

    call JGTasks.GatherTranches as SNPGatherTranches {
      input:
        tranches = SNPsVariantRecalibratorScattered.tranches,
        output_filename = participant_id + ".snps.gathered.tranches",
        mode = "SNP",
        disk_size = small_disk
    }
  }

  if (num_gvcfs <= snps_variant_recalibration_threshold) {
    call JGTasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = participant_id + ".snps.recal",
        tranches_filename = participant_id + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_vcf,
        dbsnp_resource_vcf_index = dbsnp_vcf_index,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size = small_disk
    }
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    #for really large callsets we give to friends, just apply filters to the sites-only
    call JGTasks.ApplyRecalibration {
      input:
        recalibrated_vcf_filename = participant_id + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
        snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
        snps_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        use_allele_specific_annotations = allele_specific_annotations,
        disk_size = medium_disk
    }

    # For large callsets we need to collect metrics from the shards and gather them later.
    if (!is_small_callset) {
      call JGTasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = participant_id + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size = medium_disk
      }
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
  if (is_small_callset) {
    call JGTasks.GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs = ApplyRecalibration.recalibrated_vcf,
        output_vcf_name = participant_id + ".vcf.gz",
        disk_size = huge_disk
    }

    call JGTasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = participant_id,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size = large_disk
    }
  }

  if (!is_small_callset) {
    # For large callsets we still need to gather the sharded metrics.
    call JGTasks.GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = participant_id,
        disk_size = medium_disk
    }
  }

  # CrossCheckFingerprints takes forever on large callsets.
  # We scatter over the input GVCFs to make things faster.
  if (scatter_cross_check_fingerprints) {
    call JGTasks.GetFingerprintingIntervalIndices {
      input:
        unpadded_intervals = unpadded_intervals,
        haplotype_database = haplotype_database
    }

    Array[Int] fingerprinting_indices = GetFingerprintingIntervalIndices.indices_to_fingerprint

    scatter (idx in fingerprinting_indices) {
      File vcfs_to_fingerprint = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx]
      File vcf_indices_to_fingerprint = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx]
    }

    call JGTasks.GatherVcfs as GatherFingerprintingVcfs {
      input:
        input_vcfs = vcfs_to_fingerprint,
        output_vcf_name = participant_id + ".gathered.fingerprinting.vcf.gz",
        disk_size = medium_disk
    }

    call JGTasks.SelectFingerprintSiteVariants {
      input:
        input_vcf = GatherFingerprintingVcfs.output_vcf,
        base_output_name = participant_id + ".fingerprinting",
        haplotype_database = haplotype_database,
        disk_size = medium_disk
    }

    call JGTasks.PartitionSampleNameMap {
      input:
        sample_name_map = GenerateSampleMapFile.sample_map,
        line_limit = 1000
    }

    scatter (idx in range(length(PartitionSampleNameMap.partitions))) {

      Array[File] files_in_partition = read_lines(PartitionSampleNameMap.partitions[idx])

      call JGTasks.CrossCheckFingerprint as CrossCheckFingerprintsScattered {
        input:
          gvcf_paths = files_in_partition,
          vcf_paths = vcfs_to_fingerprint,
          vcf_indices = vcf_indices_to_fingerprint,
          sample_name_map = GenerateSampleMapFile.sample_map,
          haplotype_database = haplotype_database,
          output_base_name = participant_id + "." + idx,
          scattered = true
      }
    }

    call JGTasks.GatherPicardMetrics as GatherFingerprintingMetrics {
      input:
        metrics_files = CrossCheckFingerprintsScattered.crosscheck_metrics,
        output_file_name = participant_id + ".fingerprintcheck",
        disk_size = small_disk
    }
  }

  if (!scatter_cross_check_fingerprints) {

    scatter (line in sample_name_map_lines) {
      File gvcf_paths = line[1]
    }

    call JGTasks.CrossCheckFingerprint as CrossCheckFingerprintSolo {
      input:
        gvcf_paths = gvcf_paths,
        vcf_paths = ApplyRecalibration.recalibrated_vcf,
        vcf_indices = ApplyRecalibration.recalibrated_vcf_index,
        sample_name_map = GenerateSampleMapFile.local_sample_map,
        haplotype_database = haplotype_database,
        output_base_name = participant_id
    }
  }


  call ICTasks.PreprocessIntervals {
    input:
      intervals = intervals,
      blacklist_intervals = blacklist_intervals,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_index,
      ref_fasta_dict = ref_dict,
      padding = padding,
      bin_length = readcts_binlength,
      gatk4_jar_override = gatk4_jar_override,
      gatk_docker = gatk_docker,
      mem_gb = mem_gb_for_preprocess_intervals,
      disk_space_gb = preprocess_intervals_disk,
      preemptible_attempts = preemptible_attempts
  }

 
  #CALL COLLECT READ COUNTS
  
  call ICTasks.CollectReadCounts as NormalCount {
      input:
        input_bam = select_first([NormalCramToBam.output_bam, normal_input_bam]),
        input_bam_index = select_first([NormalCramToBam.output_bai, normal_input_bam_index]),
        intervals = PreprocessIntervals.preprocessed_intervals,
        output_filename = normal_count_output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  
  call ICTasks.CollectReadCounts as TumorCount{
    input:
      input_bam = select_first([TumorCramToBam.output_bam, tumor_input_bam]),
      input_bam_index = select_first([TumorCramToBam.output_bai, tumor_input_bam_index]),
      intervals = PreprocessIntervals.preprocessed_intervals,
      output_filename = tumor_count_output_filename,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  # Get the VCFs from either code path
  File output_vcf_file = select_first([FinalGatherVcf.output_vcf,ApplyRecalibration.recalibrated_vcf])
  File output_vcf_index_file = select_first([FinalGatherVcf.output_vcf_index,ApplyRecalibration.recalibrated_vcf_index])

    
 #CALLS FOR FILTERING AND PHASING
  call ICTasks.Filter {
    input:
        input_vcf = output_vcf_file,
        outfile = participant_id + "_HC_var.vcf.gz",
        include = include,
        exclude = exclude,
        softFilter = softFilter
  }

  call ICTasks.iter_isec {
    input:
        refDir = ref_dir1000G,
        input_vcf = Filter.filtered_outputVcf,
        participant_id = participant_id,
        input_vcfindex = Filter.filtered_outputVcfIndex
  }

  call ICTasks.hetvarFilter {
    input:
        chr_int = iter_isec.chr_int,
        chr_int_indices = iter_isec.chr_int_indices,
        participant_id = participant_id
  }

  call ICTasks.eagle_phasing {
    input:
        het_var = hetvarFilter.het_filteredandIndex,
        refDir = ref_dir1000G,
        eagle_gm = geneticMap
  }

  call ICTasks.index {
    input:
        phased_files = eagle_phasing.eagle_phased
  }

  call ICTasks.query {
    input:
        chr1_het = iter_isec.final_rem,
        all_het = hetvarFilter.het_filteredonly,
        all_het_indices = hetvarFilter.het_filteredIndex,
        all_eaglephased = eagle_phasing.eagle_phased,
        all_eagle_indices = index.phased_indices,
        participant_id = participant_id
  }
#CALLS FOR OBTAINING ALLELIC COUNTS AND FRACTIONS

  if (defined(desired_binsize)){
    call ICTasks.ac_calc_agg {
        input:
            ac_rscript = ac_rscript,
            tumor_counts = TumorCount.output_counts,
            tumor_sample_id = tumor_sample_id,
            normal_sample_id = normal_sample_id,
            participant_id = participant_id,
            all_phased = query.phased_hets,
            all_allelic_depths = query.allelic_depths,
            all_hetinfo = query.variant_info,
            agg_bins = agg_bins,
            desired_binsize = desired_binsize
      }
    }
    
  if (!(defined(desired_binsize))){
    call ICTasks.ac_calc_noagg {
        input:
            ac_rscript = ac_rscript,
            tumor_counts = TumorCount.output_counts,
            tumor_sample_id = tumor_sample_id,
            normal_sample_id = normal_sample_id,
            participant_id = participant_id,
            all_phased = query.phased_hets,
            all_allelic_depths = query.allelic_depths,
            all_hetinfo = query.variant_info,
            agg_bins = agg_bins
      }
    }
  
  call ICTasks.tangent_XY {
    	input:
          tangent_rscript = tangent_rscript,
          pon = pon,
          tumor_counts = TumorCount.output_counts,
          filt_probes = filt_probes,
          ref_clinical = ref_clinical,
          cohort_clinical = cohort_clinical,
          sample_sex = sample_sex,
            
          #optional inputs
          scale_method = scale_method
    }

  call ICTasks.segment {
    	input:
          segment_rscript = segment_rscript,
          allelic_counts = select_first([ac_calc_agg.allelic_counts,ac_calc_noagg.allelic_counts]),
          tangent_normalized = tangent_XY.tangent_norm,
          participant_id = participant_id,
          opt_plts = opt_plts,
          sample_sex = sample_sex,

            #optional inputs
          smooth_region = smooth_region, 
          outlier_SD = outlier_SD,
          smooth_SD = smooth_SD,
          alpha = alpha,
          min_width = min_width,
          undo_split = undo_split,
          undo_prune = undo_prune,
          undo_SD = undo_SD
    }

# Outputs that will be retained when execution is complete
  output{
    File read_depth = select_first([ac_calc_agg.read_depth,ac_calc_noagg.read_depth])
    File seg_file = segment.seg_file
    File ac_cts = segment.ac_cts
    File segment_plots = segment.segment_plots
    File segment_cts = segment.segment_cts
    File? outlier_plots = segment.outlier_plots

  }
  meta {
    allowNestedInputs: true
  }
}
