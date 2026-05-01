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
process subset_tumor_all_shards {
  tag "${meta.id}"
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
  tuple val(meta),
        path(tumor_bam),
        path(tumor_bam_index),
        val(tumor_sample),
        path(normal_bam),
        path(normal_bam_index),
        path(interval_files, stageAs: 'intervals/*')

  output:
    path("subset_bams/*.bam"),     emit: shard_bams
    path("subset_bams/*.bam.bai"), emit: shard_bais
    path("intervals/*.intervals"), emit: shard_intervals
    tuple val(meta.id),
          val(tumor_sample),
          path(normal_bam),
          path(normal_bam_index),
          emit: sample_meta

  script:
  """
  set -euo pipefail

  mkdir -p subset_bams regions

  echo "=== subset_tumor_all_shards START ==="
  echo "sample=${meta.id}"
  echo "tumor_bam=${tumor_bam}"
  ls -lh "${tumor_bam}" || true
  echo "interval count:"
  ls intervals/*.intervals | wc -l
  echo "first intervals:"
  ls intervals/*.intervals | head

  threads=\$(( ${task.cpus} > 1 ? ${task.cpus} - 1 : 1 ))
  echo "threads=\$threads"

  for interval_file in intervals/*.intervals; do
    shard_id=\$(basename "\$interval_file" .intervals)

    echo "=== START shard=\${shard_id} at \$(date) ==="

    grep -v '^@' "\$interval_file" \
      | awk 'NF>=3 {print \$1":"\$2+1"-"\$3}' \
      > "regions/\${shard_id}.regions.txt"

    echo "region count for \${shard_id}:"
    wc -l "regions/\${shard_id}.regions.txt"
    head "regions/\${shard_id}.regions.txt" || true

    samtools view \
      -@ "\$threads" \
      -b \
      -o "subset_bams/${meta.id}.\${shard_id}.bam" \
      "${tumor_bam}" \
      \$(cat "regions/\${shard_id}.regions.txt")

    samtools index "subset_bams/${meta.id}.\${shard_id}.bam"

    echo "=== DONE shard=\${shard_id} at \$(date) ==="
    ls -lh "subset_bams/${meta.id}.\${shard_id}.bam" "subset_bams/${meta.id}.\${shard_id}.bam.bai"
  done

  echo "=== FINAL COUNTS ==="
  ls subset_bams/*.bam | wc -l
  ls subset_bams/*.bam.bai | wc -l
  echo "=== subset_tumor_all_shards DONE ==="
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
   * 1. Split intervals once.
   */
  interval_res = split_intervals(
    file(params.ref_fasta,  checkIfExists: true),
    file(params.ref_fai,    checkIfExists: true),
    file(params.ref_dict,   checkIfExists: true),
    file(params.intervals,  checkIfExists: true),
    params.scatter_count as int
  )

  /*
   * 2. Batch interval shards for subsetting.
   * This follows the one-sample pattern: each item is one batch of intervals.
   */
  def batch_size = (params.extract_batch_size ?: 10) as int

  interval_batches_ch = interval_res.interval_shards
    .flatten()
    .toSortedList { a, b -> a.name <=> b.name }
    .flatMap { intervals ->
      intervals.collate(batch_size)
    }

  /*
   * 3. Create one subsetting task per sample per interval batch.
   * We do not use combine here, because combine can spread tuple/list contents.
   */
  sample_batches_ch = interval_batches_ch
    .flatMap { interval_batch ->
      params.mutect_runs.collect { run ->

        def nbam = run.normal_reads
          ? file(run.normal_reads, checkIfExists: true)
          : file(NO_NORMAL_BAM_PATH, checkIfExists: true)

        def nbai = run.normal_reads_index
          ? file(run.normal_reads_index, checkIfExists: true)
          : file(NO_NORMAL_BAI_PATH, checkIfExists: true)

        tuple(
          [id: run.output_prefix],
          file(run.tumor_reads, checkIfExists: true),
          file(run.tumor_reads_index, checkIfExists: true),
          run.tumor_sample_name,
          nbam,
          nbai,
          interval_batch
        )
      }
    }

  /*
   * 4. Subset tumor BAMs by sample and interval batch.
   * The process names outputs as: ${sample_id}.${shard_id}.bam
   */
  subset_res = subset_tumor_all_shards(sample_batches_ch)

  /*
   * 5. Build keyed BAM/BAI channels.
   * Key = [sample_id, shard_id]
   * This prevents collisions like sample1/0000-scattered vs sample2/0000-scattered.
   */
  shard_bams_ch = subset_res.shard_bams
    .flatten()
    .map { f ->
      def m = f.name =~ /^(.+)\.(\d+-scattered)\.bam$/
      if( !m.matches() )
        error "Could not parse sample/shard from BAM name: ${f.name}"

      def sid = m[0][1]
      def shard = m[0][2]

      tuple([sid, shard], f)
    }

  shard_bais_ch = subset_res.shard_bais
    .flatten()
    .map { f ->
      def m = f.name =~ /^(.+)\.(\d+-scattered)\.bam\.bai$/
      if( !m.matches() )
        error "Could not parse sample/shard from BAI name: ${f.name}"

      def sid = m[0][1]
      def shard = m[0][2]

      tuple([sid, shard], f)
    }

  /*
   * Intervals are shared across samples, so their key is only shard_id.
   */
  shard_intervals_ch = subset_res.shard_intervals
    .flatten()
    .unique { f -> f.name }
    .map { f ->
      def shard = f.name.replaceFirst(/\.intervals$/, '')
      tuple(shard, f)
    }

  /*
   * sample_meta is emitted once per sample per batch, so deduplicate by sample_id.
   */
  sample_meta_ch = subset_res.sample_meta
    .unique { sid, tsample, nbam, nbai -> sid }
    .map { sid, tsample, nbam, nbai ->
      tuple(sid, tsample, nbam, nbai)
    }

  /*
   * 6. Join BAM + BAI by [sample_id, shard_id],
   * then attach interval by shard_id,
   * then attach normal/tumor sample metadata by sample_id.
   */
  mutect_inputs_ch = shard_bams_ch
    .join(shard_bais_ch)
    .map { key, bam, bai ->
      def sid = key[0]
      def shard = key[1]
      tuple(shard, sid, bam, bai)
    }
    .join(shard_intervals_ch)
    .map { shard, sid, bam, bai, interval ->
      tuple(sid, interval, bam, bai)
    }
    .join(sample_meta_ch)
    .map { sid, interval, bam, bai, tsample, nbam, nbai ->
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
   * 7. Split inputs for mutect_wrapper signature.
   * This also passes force_call_file when provided.
   */
  mutect_split_ch = mutect_inputs_ch.multiMap {
    sid, interval, bam, bai, ref, fai, dict, germ, nbam, nbai, tsample ->

      main:        tuple(sid, interval, bam, bai, ref, fai, dict, germ)
      nbam:        nbam
      nbai:        nbai
      alleles:     params.force_call_file
                     ? file(params.force_call_file, checkIfExists: true)
                     : file(NO_ALLELES_VCF_PATH, checkIfExists: true)
      alleles_tbi: params.force_call_file_index
                     ? file(params.force_call_file_index, checkIfExists: true)
                     : file(NO_ALLELES_TBI_PATH, checkIfExists: true)
      tsample:     tsample
      extra:       params.m2_extra_args ?: ''
  }

  /*
   * 8. Run Mutect2 once per sample-shard BAM.
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
   * 9. Gather per sample.
   */
  mutect_res.vcf
    .groupTuple(size: params.scatter_count as int)
    .set { grouped_vcfs_ch }

  gather_vcfs(grouped_vcfs_ch)
}
