def NO_NORMAL_BAM_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAM"
def NO_NORMAL_BAI_PATH   = "${workflow.projectDir}/assets/NO_NORMAL_BAI"
def NO_ALLELES_VCF_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
def NO_ALLELES_TBI_PATH  = "${workflow.projectDir}/assets/NO_ALLELES_TBI"

new File("${workflow.projectDir}/assets").mkdirs()
new File(NO_NORMAL_BAM_PATH).createNewFile()
new File(NO_NORMAL_BAI_PATH).createNewFile()
new File(NO_ALLELES_VCF_PATH).createNewFile()
new File(NO_ALLELES_TBI_PATH).createNewFile()

if( !params.containsKey('merge_all_sample_vcfs') )
  params.merge_all_sample_vcfs = false
  
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
    tuple val(meta.id),
          val(tumor_sample),
          path(normal_bam),
          path(normal_bam_index),
          path("subset_bams/*.bam"),
          path("subset_bams/*.bam.bai"),
          path("out_intervals/*.intervals"),
          emit: shards

  script:
  """
  set -euo pipefail

  mkdir -p subset_bams regions out_intervals
  cp intervals/*.intervals out_intervals/

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
  echo -n "BAMs: "
  ls subset_bams/*.bam | wc -l
  echo -n "BAIs: "
  ls subset_bams/*.bam.bai | wc -l
  echo -n "Intervals: "
  ls out_intervals/*.intervals | wc -l
  echo "=== subset_tumor_all_shards DONE ==="
  """
}

process mutect_wrapper {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"
  stageInMode 'symlink'

  input:
    tuple val(sample_id),
          path(interval_shard),
          path(tumor_bam),
          path(tumor_bam_index),
          path(ref_fasta),
          path(ref_fai),
          path(ref_dict),
          path(germline_resource),
          path(normal_bam),
          path(normal_bam_index),
          path(alleles_vcf),
          path(alleles_vcf_tbi),
          val(tumor_sample),
          val(extra_args)

  output:
    tuple val(sample_id), path("*.vcf.gz"),       emit: vcf
    tuple val(sample_id), path("*.vcf.gz.tbi"),   emit: tbi
    tuple val(sample_id), path("*.vcf.gz.stats"), emit: stats
    path "versions.yml",                          emit: versions

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
  tumor_sample=$(samtools view -H "$tumor_bam" \
    | awk -F'\t' '/^@RG/ {
        for (i=1;i<=NF;i++) {
          if ($i ~ /^SM:/) {
            sub(/^SM:/,"",$i);
            print $i
          }
        }
      }' \
    | sort -u)

  [[ -n "$tumor_sample" ]] || { echo "ERROR: No SM tag found in tumor BAM header" >&2; exit 1; }
  [[ $(echo "$tumor_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple tumor SM values: $tumor_sample" >&2; exit 1; }

  echo "Using tumor sample from BAM header: $tumor_sample"

  extra_args="!{extra_args}"

  heap_mb="!{ Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000) }"

  echo "=== mutect_wrapper START ==="
  echo "sample_id=!{sample_id}"
  echo "interval_shard=$interval_shard"
  echo "tumor_bam=$tumor_bam"
  echo "tumor_bam_index=$tumor_bam_index"
  echo "tumor_sample=$tumor_sample"

  [[ -n "$tumor_sample" ]] || { echo "ERROR: tumor_sample is empty" >&2; exit 1; }

  normal_args=""
  if [[ "$(basename "$normal_bam")" != "NO_NORMAL_BAM" ]]; then
    normal_sample=$(samtools view -H "$normal_bam" \
      | awk -F'\t' '/^@RG/ { for (i=1;i<=NF;i++) if ($i ~ /^SM:/) { sub(/^SM:/,"",$i); print $i } }' \
      | sort -u)
    [[ -n "$normal_sample" ]] || { echo "ERROR: No SM tag found in normal BAM header" >&2; exit 1; }
    [[ $(echo "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: Multiple SM values in normal BAM" >&2; exit 1; }
    normal_args="--input $normal_bam --normal-sample $normal_sample"
    echo "normal_sample=$normal_sample"
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

  test -s "${out_prefix}.vcf.gz"
  test -s "${out_prefix}.vcf.gz.tbi"
  test -s "${out_prefix}.vcf.gz.stats"

  ( gatk --version > versions.yml 2>&1 || echo "gatk --version failed (non-fatal)" > versions.yml )

  echo "=== mutect_wrapper DONE ==="
  '''
}

process gather_mutect_outputs {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  tag "${sample_id}"

  input:
    tuple val(sample_id), path(vcfs), path(stats)

  output:
    tuple val(sample_id),
          path("${sample_id}.merged.vcf.gz"),
          path("${sample_id}.merged.vcf.gz.tbi"),
          path("${sample_id}.merged.vcf.gz.stats"),
          emit: calls

    tuple val(sample_id), path("${sample_id}.merged.vcf.gz"),       emit: vcf
    tuple val(sample_id), path("${sample_id}.merged.vcf.gz.tbi"),   emit: tbi
    tuple val(sample_id), path("${sample_id}.merged.vcf.gz.stats"), emit: stats

  script:
  """
  set -euo pipefail

  echo "=== GATHER MUTECT OUTPUTS: ${sample_id} ==="

  find . -maxdepth 1 -type f -name 'out.*.vcf.gz' -print \
    | sed -E 's/.*out\.([0-9]+)-scattered\.vcf\.gz/\1\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > vcfs.sorted.list

  find . -maxdepth 1 -type f -name 'out.*.vcf.gz.stats' -print \
    | sed -E 's/.*out\.([0-9]+)-scattered\.vcf\.gz\.stats/\1\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > stats.sorted.list

  vcf_count=\$(wc -l < vcfs.sorted.list)
  stats_count=\$(wc -l < stats.sorted.list)

  echo "VCF shards: \$vcf_count"
  echo "Stats shards: \$stats_count"

  if [[ "\$vcf_count" -ne "\$stats_count" ]]; then
    echo "ERROR: Number of VCF files (\$vcf_count) does not match number of stats files (\$stats_count)" >&2
    exit 1
  fi

  if [[ "\$vcf_count" -ne ${params.scatter_count as int} ]]; then
    echo "ERROR: Expected ${params.scatter_count as int} shards but found \$vcf_count" >&2
    exit 1
  fi

  awk '{print "--INPUT", \$0}' vcfs.sorted.list > gather.args

  gatk --java-options "-Xmx8g -XX:-UsePerfData" GatherVcfs \
    --arguments_file gather.args \
    -O ${sample_id}.merged.vcf.gz

  if [[ ! -s ${sample_id}.merged.vcf.gz.tbi ]]; then
    tabix -f -p vcf ${sample_id}.merged.vcf.gz \
      || gatk IndexFeatureFile -I ${sample_id}.merged.vcf.gz
  fi

  stats_args=()
  while IFS= read -r stats_file; do
    stats_args+=(-stats "\$stats_file")
  done < stats.sorted.list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" MergeMutectStats \
    "\${stats_args[@]}" \
    -O ${sample_id}.merged.vcf.gz.stats

  test -s ${sample_id}.merged.vcf.gz
  test -s ${sample_id}.merged.vcf.gz.tbi
  test -s ${sample_id}.merged.vcf.gz.stats

  echo "=== GATHERED OUTPUTS ==="
  ls -lh \
    ${sample_id}.merged.vcf.gz \
    ${sample_id}.merged.vcf.gz.tbi \
    ${sample_id}.merged.vcf.gz.stats
  """
}

process filter_mutect_calls {
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  tag "${sample_id}"

  input:
    tuple val(sample_id),
          path(vcf),
          path(vcf_tbi),
          path(stats)
    path ref_fasta
    path ref_fai
    path ref_dict

  output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"),       emit: vcf
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz.tbi"),   emit: tbi
    tuple val(sample_id), path("${sample_id}.filteringStats.tsv"),    emit: stats

  script:
  """
  set -euo pipefail

  echo "=== FILTER MUTECT CALLS: ${sample_id} ==="
  echo "VCF: ${vcf}"
  echo "Stats: ${stats}"

  gatk --java-options "-Xmx8g -XX:-UsePerfData" FilterMutectCalls \
    -R ${ref_fasta} \
    -V ${vcf} \
    -stats ${stats} \
    -O ${sample_id}.filtered.vcf.gz \
    --filtering-stats ${sample_id}.filteringStats.tsv \
    --tmp-dir .

  if [[ ! -s ${sample_id}.filtered.vcf.gz.tbi ]]; then
    tabix -f -p vcf ${sample_id}.filtered.vcf.gz \
      || gatk IndexFeatureFile -I ${sample_id}.filtered.vcf.gz
  fi

  test -s ${sample_id}.filtered.vcf.gz
  test -s ${sample_id}.filtered.vcf.gz.tbi
  test -s ${sample_id}.filteringStats.tsv

  echo "=== FILTERED OUTPUTS ==="
  ls -lh \
    ${sample_id}.filtered.vcf.gz \
    ${sample_id}.filtered.vcf.gz.tbi \
    ${sample_id}.filteringStats.tsv
  """
}

process merge_all_sample_vcfs {
  label 'process_medium'
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

  input:
    path vcfs
    path tbis

  output:
    path "all_samples.merged.vcf.gz",     emit: vcf
    path "all_samples.merged.vcf.gz.tbi", emit: tbi

  script:
  """
  set -euo pipefail

  echo "=== INPUT SAMPLE-LEVEL VCFS ==="
  ls -lh *.merged.vcf.gz

  rm -f tumor_vcfs.list
  mkdir -p tumor_only_vcfs

  echo "=== FILTERING OUT PBMC/NORMAL SAMPLES ==="

  for f in \$(ls -1 *.merged.vcf.gz | sort); do
    echo "Input: \$f"

    tumor_samples=\$(bcftools query -l "\$f" | grep -viE 'PBMC|NORMAL|BLOOD' || true)

    if [ -z "\$tumor_samples" ]; then
      echo "ERROR: No tumor samples found in \$f after removing PBMC/NORMAL/BLOOD"
      exit 1
    fi

    echo "Keeping tumor samples:"
    echo "\$tumor_samples"

    out="tumor_only_vcfs/\$(basename "\$f" .vcf.gz).tumor_only.vcf.gz"

    echo "\$tumor_samples" > keep_samples.txt

    bcftools view \
      -S keep_samples.txt \
      -Oz \
      -o "\$out" \
      "\$f"

    tabix -f -p vcf "\$out"

    echo "\$out" >> tumor_vcfs.list
  done

  echo "=== TUMOR-ONLY VCFS TO MERGE ==="
  cat tumor_vcfs.list

  bcftools merge \
    --force-samples \
    -Oz \
    -o all_samples.merged.vcf.gz \
    --file-list tumor_vcfs.list

  tabix -f -p vcf all_samples.merged.vcf.gz

  echo "=== FINAL SAMPLE NAMES ==="
  bcftools query -l all_samples.merged.vcf.gz

  echo "=== FINAL COHORT VCF ==="
  ls -lh all_samples.merged.vcf.gz all_samples.merged.vcf.gz.tbi
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
   */
  def batch_size = (params.extract_batch_size ?: 10) as int

  interval_batches_ch = interval_res.interval_shards
    .flatten()
    .toSortedList { a, b -> a.name <=> b.name }
    .flatMap { intervals ->
      intervals.collate(batch_size)
    }

  /*
   * 3. One subsetting task per sample per interval batch.
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
   * 4. Subset tumor BAMs.
   */
  subset_res = subset_tumor_all_shards(sample_batches_ch)

  /*
   * 5. One emitted BAM shard -> one Mutect task.
   */
  mutect_inputs_ch = subset_res.shards
    .flatMap { sid, tsample, nbam, nbai, bams, bais, intervals ->

      def bam_list = bams instanceof List ? bams : [bams]
      def bai_list = bais instanceof List ? bais : [bais]
      def int_list = intervals instanceof List ? intervals : [intervals]

      bam_list = bam_list.sort { it.name }
      bai_list = bai_list.sort { it.name }
      int_list = int_list.sort { it.name }

      bam_list.collect { bam ->

        def m = bam.name =~ /^(.+)\.(\d+-scattered)\.bam$/
        if( !m.matches() )
          error "Could not parse sample/shard from BAM name: ${bam.name}"

        def sample_from_bam = m[0][1]
        def shard_id = m[0][2]

        def bai = bai_list.find { it.name == "${sample_from_bam}.${shard_id}.bam.bai" }
        if( bai == null )
          error "Could not find BAI for BAM: ${bam.name}"

        def interval = int_list.find { it.name == "${shard_id}.intervals" }
        if( interval == null )
          error "Could not find interval for shard: ${shard_id}"

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
          params.force_call_file
            ? file(params.force_call_file, checkIfExists: true)
            : file(NO_ALLELES_VCF_PATH, checkIfExists: true),
          params.force_call_file_index
            ? file(params.force_call_file_index, checkIfExists: true)
            : file(NO_ALLELES_TBI_PATH, checkIfExists: true),
          tsample,
          params.m2_extra_args ?: ''
        )
      }
    }

  mutect_inputs_ch.view { x ->
    "MUTECT FANOUT: sample=${x[0]}, interval=${x[1].name}, bam=${x[2].name}"
  }

  mutect_inputs_ch
    .map { 1 }
    .reduce { a, b -> a + b }
    .view { n -> "TOTAL MUTECT TASKS EXPECTED: ${n}" }

  /*
   * 6. Run Mutect2 once per BAM shard.
   */
  mutect_res = mutect_wrapper(mutect_inputs_ch)

     /*
   * 7. Gather per-shard Mutect VCFs into one merged VCF per sample.
   *    This is sample-level stitching, not cross-sample merging.
   */
  grouped_vcfs_ch = mutect_res.vcf
    .groupTuple(size: params.scatter_count as int)

  grouped_stats_ch = mutect_res.stats
    .groupTuple(size: params.scatter_count as int)

  grouped_mutect_outputs_ch = grouped_vcfs_ch
    .join(grouped_stats_ch)
    .map { sample_id, vcfs, stats ->
      tuple(sample_id, vcfs, stats)
    }

  grouped_mutect_outputs_ch.view { sample_id, vcfs, stats ->
    "GATHER INPUT: sample=${sample_id}, vcfs=${vcfs.size()}, stats=${stats.size()}"
  }

  gathered_mutect_res = gather_mutect_outputs(grouped_mutect_outputs_ch)

  gathered_mutect_res.vcf.view { sid, vcf ->
    "MERGED MUTECT VCF PER SAMPLE: sample=${sid}, vcf=${vcf}"
  }

  /*
   * 8. Filter each gathered per-sample Mutect2 VCF using its merged stats.
   */
  filtered_mutect_res = filter_mutect_calls(
    gathered_mutect_res.calls,
    file(params.ref_fasta, checkIfExists: true),
    file(params.ref_fai,   checkIfExists: true),
    file(params.ref_dict,  checkIfExists: true)
  )

  filtered_mutect_res.vcf.view { sid, vcf ->
    "FILTERED MUTECT VCF PER SAMPLE: sample=${sid}, vcf=${vcf}"
  }

  /*
   * 9. Optionally merge all filtered per-sample VCFs into one cohort-level VCF.
   */
  if( params.merge_all_sample_vcfs ) {

    filtered_mutect_res.vcf
      .map { sid, vcf -> vcf }
      .collect()
      .set { all_sample_vcfs_ch }

    filtered_mutect_res.tbi
      .map { sid, tbi -> tbi }
      .collect()
      .set { all_sample_tbis_ch }

    merge_all_sample_vcfs(all_sample_vcfs_ch, all_sample_tbis_ch)

  } else {
    log.info "Skipping cross-sample VCF merge because params.merge_all_sample_vcfs=false"
  }

}