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

process subset_tumor_per_shard {
  tag "${meta.id}"
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index)
    path interval_files

  output:
    tuple val(meta.id), path("shards/*.bam"),          emit: shard_bams
    tuple val(meta.id), path("shards/*.bam.bai"),      emit: shard_bais
    tuple val(meta.id), path("shards/*.intervals"),    emit: shard_intervals
    tuple val(meta.id), path("tumor_sample_name.txt"), emit: tumor_sample

  script:
  """
  set -euo pipefail
  mkdir -p shards

  tumor_sample=\$(samtools view -H "${tumor_bam}" \\
    | awk -F'\\t' '/^@RG/ {
        for (i=1;i<=NF;i++)
          if (\$i ~ /^SM:/) { sub(/^SM:/,"",\$i); print \$i }
      }' | sort -u)

  [[ -n "\$tumor_sample" ]] \\
    || { echo "ERROR: No SM tag in tumor BAM header" >&2; exit 1; }
  [[ \$(echo "\$tumor_sample" | wc -l) -eq 1 ]] \\
    || { echo "ERROR: Multiple SM values: \$tumor_sample" >&2; exit 1; }

  echo "\$tumor_sample" > tumor_sample_name.txt

  for interval_file in *.intervals; do
    shard_base=\$(basename "\$interval_file" .intervals)

    cp "\$interval_file" "shards/\${shard_base}.intervals"

    grep -v '^@' "\$interval_file" \\
      | awk 'NF>=3 {print \$1":"\$2+1"-"\$3}' \\
      > /tmp/regions_\${shard_base}.txt

    [[ -s /tmp/regions_\${shard_base}.txt ]] \\
      || { echo "ERROR: no regions for \${shard_base}" >&2; exit 1; }

    readarray -t regions < /tmp/regions_\${shard_base}.txt

    samtools view \\
      -@ \$(( ${task.cpus} - 1 )) \\
      -b \\
      -o "shards/\${shard_base}.bam" \\
      "${tumor_bam}" \\
      "\${regions[@]}"

    samtools index "shards/\${shard_base}.bam"
  done
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

  interval_res = split_intervals(
    file(params.ref_fasta,  checkIfExists: true),
    file(params.ref_fai,    checkIfExists: true),
    file(params.ref_dict,   checkIfExists: true),
    file(params.intervals,  checkIfExists: true),
    params.scatter_count as int
  )

  println "mutect_runs size = ${params.mutect_runs?.size()}"

  Channel.fromList(params.mutect_runs)
    .view { "RUN_RAW: ${it.output_prefix} :: ${it.tumor_reads}" }

  runs_ch = Channel.fromList(params.mutect_runs)
    .map { run ->
        tuple(
            [id: run.output_prefix],
            file(run.tumor_reads),
            file(run.tumor_reads_index),
            run.normal_reads       ?: null,
            run.normal_reads_index ?: null,
            run.tumor_sample_name
        )
    }

  runs_ch.view { "RUNS_CH: $it" }

  intervals_ready = interval_res.interval_shards.collect()

  subset_res = subset_tumor_per_shard(
    runs_ch.map { meta, tbam, tbai, nbam, nbai, tsample ->
        tuple(meta, tbam, tbai)
    },
    intervals_ready
  )
    
  subset_res.shard_bams.view      { "SHARD_BAMS_RAW: $it" }
    subset_res.shard_bais.view      { "SHARD_BAIS_RAW: $it" }
    subset_res.shard_intervals.view { "SHARD_INTERVALS_RAW: $it" }

    shard_bams_ch = subset_res.shard_bams
      .transpose()
      .map { sid, f -> tuple(sid, f.name.replaceFirst(/\.bam$/, ''), f) }

    shard_bais_ch = subset_res.shard_bais
      .transpose()
      .map { sid, f -> tuple(sid, f.name.replaceFirst(/\.bam\.bai$/, ''), f) }

    shard_intervals_ch = subset_res.shard_intervals
      .transpose()
      .map { sid, f -> tuple(sid, f.name.replaceFirst(/\.intervals$/, ''), f) }

    shard_bams_ch.view      { "SHARD_BAMS: $it" }
    shard_bais_ch.view      { "SHARD_BAIS: $it" }
    shard_intervals_ch.view { "SHARD_INTERVALS: $it" }

    mutect_inputs_ch = shard_bams_ch
      .join(shard_bais_ch,      by: [0, 1])
      .join(shard_intervals_ch, by: [0, 1])
      .map { sid, base, bam, bai, interval ->
        tuple(
          sid,
          interval, bam, bai,
          file(params.ref_fasta,         checkIfExists: true),
          file(params.ref_fai,           checkIfExists: true),
          file(params.ref_dict,          checkIfExists: true),
          file(params.germline_resource, checkIfExists: true)
        )
      }

    mutect_inputs_ch.view { "MUTECT_INPUT: $it" }

    normals_ch = runs_ch.map { meta, tbam, tbai, nbam, nbai, tsample ->
      def nbam_file = nbam ? file(nbam) : file(NO_NORMAL_BAM_PATH)
      def nbai_file = nbam ? file(nbai) : file(NO_NORMAL_BAI_PATH)
      tuple(meta.id, nbam_file, nbai_file, tsample)
    }

    mutect_inputs_ch
      .join(normals_ch, by: 0)
      .map { sid, interval, bam, bai, ref, fai, dict, germ, nbam, nbai, tsample ->
        tuple(
          sid,                                   // <-- carry sid through
          tuple(interval, bam, bai, ref, fai, dict, germ),
          nbam, nbai,
          file(NO_ALLELES_VCF_PATH),
          file(NO_ALLELES_TBI_PATH),
          tsample,
          params.m2_extra_args ?: ''
        )
      }
      .multiMap { sid, main_tuple, nbam, nbai, alleles, alleles_tbi, tsample, extra ->
          main:        tuple(sid, main_tuple[0], main_tuple[1], main_tuple[2], main_tuple[3], main_tuple[4], main_tuple[5], main_tuple[6])
          nbam:        nbam
          nbai:        nbai
          alleles:     alleles
          alleles_tbi: alleles_tbi
          tsample:     tsample
          extra:       extra
      }
      .set { mutect_split_ch }

    mutect_res = mutect_wrapper(
      mutect_split_ch.main,
      mutect_split_ch.nbam,
      mutect_split_ch.nbai,
      mutect_split_ch.alleles,
      mutect_split_ch.alleles_tbi,
      mutect_split_ch.tsample,
      mutect_split_ch.extra
    )

    mutect_res.vcf
      .view { "MUTECT_VCF_OUT: $it" }
      .groupTuple(size: params.scatter_count)
      .set { grouped_vcfs_ch }

    gather_vcfs(grouped_vcfs_ch)
}