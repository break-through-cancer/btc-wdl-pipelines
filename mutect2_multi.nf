def NO_NORMAL_BAM_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAM"
def NO_NORMAL_BAI_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAI"
def NO_ALLELES_VCF_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
def NO_ALLELES_TBI_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_TBI"

new File("${workflow.projectDir}/assets").mkdirs()
new File(NO_NORMAL_BAM_PATH).createNewFile()
new File(NO_NORMAL_BAI_PATH).createNewFile()
new File(NO_ALLELES_VCF_PATH).createNewFile()
new File(NO_ALLELES_TBI_PATH).createNewFile()

if( !params.containsKey('tumor_sample') )
  params.tumor_sample = null

if( !params.containsKey('m2_extra_args') )
  params.m2_extra_args = ''

if( !params.containsKey('normal_reads') )
  params.normal_reads = null
if( !params.containsKey('normal_reads_index') )
  params.normal_reads_index = null
if( !params.containsKey('force_call_file') )
  params.force_call_file = null
if( !params.containsKey('force_call_file_index') )
  params.force_call_file_index = null

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
    path "scattered/*.intervals", emit: interval_shards

  script:
  """
  set -euo pipefail
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
  """
}

process subset_tumor_one_shard {
  tag "${meta.id}:${shard_id}"
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(meta),
          path(tumor_bam),
          path(tumor_bam_index),
          path(interval_file),
          val(shard_id),
          val(tumor_sample),
          path(normal_bam),
          path(normal_bam_index)

  output:
    tuple val(meta.id),
          val(shard_id),
          path("${meta.id}.${shard_id}.bam"),
          path("${meta.id}.${shard_id}.bam.bai"),
          path(interval_file),
          val(tumor_sample),
          path(normal_bam),
          path(normal_bam_index),
          emit: subsetted

  script:
  """
  set -euo pipefail

  grep -v '^@' "${interval_file}" \
    | awk 'NF>=3 {print \$1":"\$2+1"-"\$3}' \
    > regions.txt

  samtools view \
    -@ \$(( ${task.cpus} - 1 )) \
    -b \
    -o "${meta.id}.${shard_id}.bam" \
    "${tumor_bam}" \
    \$(cat regions.txt)

  samtools index "${meta.id}.${shard_id}.bam"
  """
}

process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"
  stageInMode 'symlink'

  input:
    tuple(
      val(sample_id),
      path(interval_shard),
      path(tumor_bam),
      path(tumor_bam_index),
      path(ref_fasta),
      path(ref_fai),
      path(ref_dict),
      path(germline_resource)
    )
    path normal_bam
    path normal_bam_index
    path alleles_vcf
    path alleles_vcf_tbi
    val  tumor_sample
    val  extra_args

  output:
    tuple val(sample_id), path("*.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                        emit: versions

  shell:
  '''
  set -euo pipefail

  shard_base=$(basename "!{interval_shard}" .intervals)
  tumor_bam="!{tumor_bam}"
  tumor_bam_index="!{tumor_bam_index}"
  interval_shard="!{interval_shard}"
  ref_fasta="!{ref_fasta}"
  ref_fai="!{ref_fai}"
  ref_dict="!{ref_dict}"
  germline_resource="!{germline_resource}"
  normal_bam="!{normal_bam}"
  normal_bam_index="!{normal_bam_index}"
  alleles_vcf="!{alleles_vcf}"
  alleles_vcf_tbi="!{alleles_vcf_tbi}"
  tumor_sample="!{tumor_sample}"
  extra_args="!{extra_args}"

  heap_mb="!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }"

  [[ -n "$tumor_sample" ]] || { echo "ERROR: tumor_sample is empty" >&2; exit 1; }

  normal_args=""
  if [[ "$(basename "$normal_bam")" != "NO_NORMAL_BAM" ]]; then
    normal_sample=$(samtools view -H "$normal_bam" \
      | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
      | sort -u)
    [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
    [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values in normal BAM" >&2; exit 1; }
    normal_args="--input $normal_bam --normal-sample $normal_sample"
  fi

  alleles_args=""
  if [[ "$(basename "$alleles_vcf")" != "NO_ALLELES_VCF" ]]; then
    alleles_args="--alleles $alleles_vcf"
  fi

  if [[ ! -f "${germline_resource}.tbi" ]]; then
    tabix -f -p vcf "$germline_resource"
  fi
  test -s "${germline_resource}.tbi"

  out_prefix="out.${shard_base}"

  gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \
    --input "$tumor_bam" \
    ${normal_args} \
    --reference "$ref_fasta" \
    --germline-resource "$germline_resource" \
    --intervals "$interval_shard" \
    --tmp-dir . \
    --tumor-sample "$tumor_sample" \
    ${alleles_args} \
    ${extra_args} \
    --output "${out_prefix}.vcf.gz"

  ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )
  '''
}

process gather_vcfs {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  tag "${sample_id}"

  input:
    tuple val(sample_id), path(vcfs)

  output:
    tuple val(sample_id), path("${sample_id}.merged.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.merged.vcf.gz.tbi"), emit: tbi

  script:
  """
  set -euo pipefail

  find . -maxdepth 1 -type f -name '*.vcf.gz' -print | sort > vcfs.list

  sed -E 's/.*out\\.([0-9]+).*/\\1\\t&/' vcfs.list | sort -k1,1n | cut -f2- > vcfs.sorted.list
  awk '{print "--INPUT", \$0}' vcfs.sorted.list > gather.args

  gatk GatherVcfs --arguments_file gather.args -O ${sample_id}.merged.vcf.gz

  if [ ! -s ${sample_id}.merged.vcf.gz.tbi ]; then
    tabix -f -p vcf ${sample_id}.merged.vcf.gz || gatk IndexFeatureFile -I ${sample_id}.merged.vcf.gz
  fi
  """
}

workflow {

  /*
   * 1. Split intervals → emit individual shard files
   */
  interval_res = split_intervals(
    file(params.ref_fasta,  checkIfExists: true),
    file(params.ref_fai,    checkIfExists: true),
    file(params.ref_dict,   checkIfExists: true),
    file(params.intervals,  checkIfExists: true),
    params.scatter_count as int
  )

  intervals_ch = interval_res.interval_shards
    .flatten()
    .map { f ->
      def shard_id = f.name.replaceFirst(/\.intervals$/, '')
      tuple(shard_id, f)
    }

  /*
   * 2. Build sample channel
   */
  runs_ch = Channel.fromList(params.mutect_runs)
    .map { run ->

      def nbam = run.normal_reads
        ? file(run.normal_reads)
        : file(NO_NORMAL_BAM_PATH)

      def nbai = run.normal_reads_index
        ? file(run.normal_reads_index)
        : file(NO_NORMAL_BAI_PATH)

      tuple(
        [id: run.output_prefix],
        file(run.tumor_reads),
        file(run.tumor_reads_index),
        run.tumor_sample_name,
        nbam,
        nbai
      )
    }

  /*
   * 3. CROSS product: sample × shard
   */
  sample_shard_ch = runs_ch
    .combine(intervals_ch)
    .map { meta, tbam, tbai, tsample, nbam, nbai, shard_id, interval_file ->
      tuple(meta, tbam, tbai, interval_file, shard_id, tsample, nbam, nbai)
    }

  /*
   * 4. Subset BAM per (sample, shard)
   */
  subset_res = subset_tumor_one_shard(sample_shard_ch)

  /*
   * 5. Build Mutect inputs (already shard-parallel)
   */
  mutect_main_ch = subset_res.subsetted.map {
    sid, shard_id, bam, bai, interval, tsample, nbam, nbai ->

      tuple(
        sid,
        interval,
        bam,
        bai,
        file(params.ref_fasta,         checkIfExists: true),
        file(params.ref_fai,           checkIfExists: true),
        file(params.ref_dict,          checkIfExists: true),
        file(params.germline_resource, checkIfExists: true),
        nbam,
        nbai,
        tsample
      )
  }

  /*
   * 6. Split inputs for mutect_wrapper signature
   */
  mutect_split_ch = mutect_main_ch.multiMap {
    sid, interval, bam, bai, ref, fai, dict, germ, nbam, nbai, tsample ->

      main:        tuple(sid, interval, bam, bai, ref, fai, dict, germ)
      nbam:        nbam
      nbai:        nbai
      alleles:     file(NO_ALLELES_VCF_PATH)
      alleles_tbi: file(NO_ALLELES_TBI_PATH)
      tsample:     tsample
      extra:       params.m2_extra_args ?: ''
  }

  /*
   * 7. Run Mutect fully parallel (sample × shard)
   */
  mutect_res = mutect_wrapper(
    mutect_split_ch.main,
    mutect_split_ch.nbam,
    mutect_split_ch.nbai,
    mutect_split_ch.alleles,
    mutect_split_ch.alleles_tbi,
    mutect_split_ch.tsample,
    mutect_split_ch.extra
  )

  /*
   * 8. Group per sample → gather VCFs
   */
  mutect_res.vcf
    .groupTuple(size: params.scatter_count)
    .set { grouped_vcfs_ch }

  gather_vcfs(grouped_vcfs_ch)
}