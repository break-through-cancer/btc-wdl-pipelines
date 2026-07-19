def NO_ALLELES_VCF_PATH = "${workflow.projectDir}/assets/NO_ALLELES_VCF"
def NO_ALLELES_TBI_PATH = "${workflow.projectDir}/assets/NO_ALLELES_TBI"
def NO_PON_VCF_PATH = "${workflow.projectDir}/assets/NO_PON_VCF"
def NO_PON_TBI_PATH = "${workflow.projectDir}/assets/NO_PON_TBI"

new File("${workflow.projectDir}/assets").mkdirs()
new File(NO_ALLELES_VCF_PATH).createNewFile()
new File(NO_ALLELES_TBI_PATH).createNewFile()
new File(NO_PON_VCF_PATH).createNewFile()
new File(NO_PON_TBI_PATH).createNewFile()

if( !params.containsKey('merge_all_sample_vcfs') )
  params.merge_all_sample_vcfs = false

if( !params.containsKey('m2_extra_args') )
  params.m2_extra_args = ''



/*
 * Split the requested calling regions into balanced GATK interval-list shards.
 */
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
    path "scattered/*.interval_list", emit: interval_shards

  script:
  """
  set -euo pipefail

  mkdir -p scattered

  gatk --java-options "-Xmx8g -XX:-UsePerfData" BedToIntervalList \
    -I "${intervals}" \
    -SD "${ref_dict}" \
    -O regions.interval_list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" SplitIntervals \
    -R "${ref_fasta}" \
    -L regions.interval_list \
    --scatter "${scatter_count}" \
    -O scattered
  """
}


/*
 * Create tumor BAM shards for one tumor sample.
 */
process subset_tumor_all_shards {
  tag "${meta.id}"
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(meta),
          path(tumor_bam),
          path(tumor_bam_index),
          path(interval_files, stageAs: 'intervals/*')

  output:
    tuple val(meta.id),
          path("subset_bams/*.tumor.*.bam"),
          path("subset_bams/*.tumor.*.bam.bai"),
          path("out_intervals/*.interval_list"),
          emit: shards

  script:
  """
  set -euo pipefail

  mkdir -p subset_bams regions out_intervals
  cp intervals/*.interval_list out_intervals/

  threads=\$(( ${task.cpus} > 1 ? ${task.cpus} - 1 : 1 ))

  for interval_file in intervals/*.interval_list; do
    shard_id=\$(basename "\${interval_file}" .interval_list)
    regions_file="regions/\${shard_id}.regions.txt"

    grep -v '^@' "\${interval_file}" \
      | awk 'NF >= 3 { print \$1 ":" \$2 "-" \$3 }' \
      > "\${regions_file}"

    tumor_out="subset_bams/${meta.id}.tumor.\${shard_id}.bam"

    samtools view \
      -@ "\${threads}" \
      -b \
      -o "\${tumor_out}" \
      "${tumor_bam}" \
      \$(cat "\${regions_file}")

    samtools index \
      -@ "\${threads}" \
      "\${tumor_out}"

    test -s "\${tumor_out}"
    test -s "\${tumor_out}.bai"
  done

  tumor_count=\$(find subset_bams -maxdepth 1 -name '*.tumor.*.bam' | wc -l)
  interval_count=\$(find out_intervals -maxdepth 1 -name '*.interval_list' | wc -l)

  if [[ "\${tumor_count}" -ne "\${interval_count}" ]]; then
    echo "ERROR: Tumor BAM shard count does not equal interval count" >&2
    exit 1
  fi
  """
}


/*
 * Create shards for the one shared matched-normal BAM.
 */
process subset_shared_normal_all_shards {
  tag "shared_normal"
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple path(normal_bam),
          path(normal_bam_index),
          path(interval_files, stageAs: 'intervals/*')

  output:
    tuple path("subset_bams/shared_normal.normal.*.bam"),
          path("subset_bams/shared_normal.normal.*.bam.bai"),
          emit: shards

  script:
  """
  set -euo pipefail

  mkdir -p subset_bams regions

  threads=\$(( ${task.cpus} > 1 ? ${task.cpus} - 1 : 1 ))

  for interval_file in intervals/*.interval_list; do
    shard_id=\$(basename "\${interval_file}" .interval_list)
    regions_file="regions/\${shard_id}.regions.txt"

    grep -v '^@' "\${interval_file}" \
      | awk 'NF >= 3 { print \$1 ":" \$2 "-" \$3 }' \
      > "\${regions_file}"

    normal_out="subset_bams/shared_normal.normal.\${shard_id}.bam"

    samtools view \
      -@ "\${threads}" \
      -b \
      -o "\${normal_out}" \
      "${normal_bam}" \
      \$(cat "\${regions_file}")

    samtools index \
      -@ "\${threads}" \
      "\${normal_out}"

    test -s "\${normal_out}"
    test -s "\${normal_out}.bai"
  done

  normal_count=\$(find subset_bams -maxdepth 1 -name 'shared_normal.normal.*.bam' | wc -l)
  interval_count=\$(find intervals -maxdepth 1 -name '*.interval_list' | wc -l)

  if [[ "\${normal_count}" -ne "\${interval_count}" ]]; then
    echo "ERROR: Normal BAM shard count does not equal interval count" >&2
    exit 1
  fi
  """
}


/*
 * Run Mutect2, emit F1R2 counts, and calculate tumor/normal pileup summaries
 * while the interval-specific BAM shards are already staged in the task.
 */
process mutect_wrapper {
  tag "${sample_id}.${interval_shard.name}"
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"
  stageInMode 'symlink'

  input:
    tuple val(sample_id),
          path(interval_shard),
          path(tumor_bam),
          path(tumor_bam_index),
          path(normal_bam),
          path(normal_bam_index),
          path(ref_fasta),
          path(ref_fai),
          path(ref_dict),
          path(germline_resource),
          path(germline_resource_index),
          path(panel_of_normals),
          path(panel_of_normals_index),
          path(alleles_vcf),
          path(alleles_vcf_tbi),
          val(extra_args)

  output:
    tuple val(sample_id), path("out.*-scattered.vcf.gz"),                emit: vcf
    tuple val(sample_id), path("out.*-scattered.vcf.gz.tbi"),            emit: tbi
    tuple val(sample_id), path("out.*-scattered.vcf.gz.stats"),          emit: stats
    tuple val(sample_id), path("out.*-scattered.f1r2.tar.gz"),           emit: f1r2
    tuple val(sample_id), path("out.*-scattered.tumor.pileups.table"),   emit: tumor_pileup
    tuple val(sample_id), path("out.*-scattered.normal.pileups.table"),  emit: normal_pileup
    path "versions.yml",                                                  emit: versions

  script:
  def heap_mb = Math.min(task.memory ? (task.memory.mega * 0.8).intValue() : 3072, 24000)

  """
  set -euo pipefail

  shard_base=\$(basename "${interval_shard}" .interval_list)
  out_prefix="out.\${shard_base}"

  tumor_sample=\$(samtools view -H "${tumor_bam}" \
    | awk -F'\\t' '/^@RG/ {
        for (i=1; i<=NF; i++) {
          if (\$i ~ /^SM:/) {
            sub(/^SM:/, "", \$i)
            print \$i
          }
        }
      }' \
    | sort -u)

  normal_sample=\$(samtools view -H "${normal_bam}" \
    | awk -F'\\t' '/^@RG/ {
        for (i=1; i<=NF; i++) {
          if (\$i ~ /^SM:/) {
            sub(/^SM:/, "", \$i)
            print \$i
          }
        }
      }' \
    | sort -u)

  [[ -n "\${tumor_sample}" ]] || {
    echo "ERROR: No SM tag found in tumor BAM header" >&2
    exit 1
  }

  [[ -n "\${normal_sample}" ]] || {
    echo "ERROR: No SM tag found in normal BAM header" >&2
    exit 1
  }

  [[ \$(printf '%s\\n' "\${tumor_sample}" | wc -l) -eq 1 ]] || {
    echo "ERROR: Multiple tumor SM values: \${tumor_sample}" >&2
    exit 1
  }

  [[ \$(printf '%s\\n' "\${normal_sample}" | wc -l) -eq 1 ]] || {
    echo "ERROR: Multiple normal SM values: \${normal_sample}" >&2
    exit 1
  }

  if [[ "\${tumor_sample}" == "\${normal_sample}" ]]; then
    echo "ERROR: Tumor and normal have the same SM value: \${tumor_sample}" >&2
    exit 1
  fi

  pon_args=()
  if [[ "\$(basename "${panel_of_normals}")" != "NO_PON_VCF" ]]; then
    pon_args+=(--panel-of-normals "${panel_of_normals}")
  fi

  alleles_args=()
  if [[ "\$(basename "${alleles_vcf}")" != "NO_ALLELES_VCF" ]]; then
    alleles_args+=(--alleles "${alleles_vcf}")
  fi

  echo "=== MUTECT2: ${sample_id} / \${shard_base} ==="
  echo "tumor_sample=\${tumor_sample}"
  echo "normal_sample=\${normal_sample}"

  gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" Mutect2 \
    --input "${tumor_bam}" \
    --input "${normal_bam}" \
    --normal-sample "\${normal_sample}" \
    --tumor-sample "\${tumor_sample}" \
    --reference "${ref_fasta}" \
    --germline-resource "${germline_resource}" \
    "\${pon_args[@]}" \
    --intervals "${interval_shard}" \
    "\${alleles_args[@]}" \
    --f1r2-tar-gz "\${out_prefix}.f1r2.tar.gz" \
    --tmp-dir . \
    ${extra_args} \
    --output "\${out_prefix}.vcf.gz"

  # Use the same common germline resource as nf-core's pileup branch.
  # It must be a sites-only/biallelic SNP resource with INFO/AF values.
  gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" GetPileupSummaries \
    --input "${tumor_bam}" \
    --variant "${germline_resource}" \
    --intervals "${interval_shard}" \
    --reference "${ref_fasta}" \
    --output "\${out_prefix}.tumor.pileups.table"

  gatk --java-options "-Xmx${heap_mb}M -XX:-UsePerfData" GetPileupSummaries \
    --input "${normal_bam}" \
    --variant "${germline_resource}" \
    --intervals "${interval_shard}" \
    --reference "${ref_fasta}" \
    --output "\${out_prefix}.normal.pileups.table"

  test -s "\${out_prefix}.vcf.gz"
  test -s "\${out_prefix}.vcf.gz.tbi"
  test -s "\${out_prefix}.vcf.gz.stats"
  test -s "\${out_prefix}.f1r2.tar.gz"
  test -s "\${out_prefix}.tumor.pileups.table"
  test -s "\${out_prefix}.normal.pileups.table"

  (
    gatk --version > versions.yml 2>&1 \
      || echo "gatk --version failed (non-fatal)" > versions.yml
  )
  """
}


/*
 * Gather per-shard VCFs and merge the companion Mutect stats.
 */
process gather_mutect_outputs {
  tag "${sample_id}"
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(sample_id), path(vcfs), path(stats)
    val expected_shards

  output:
    tuple val(sample_id),
          path("${sample_id}.merged.vcf.gz"),
          path("${sample_id}.merged.vcf.gz.tbi"),
          path("${sample_id}.merged.vcf.gz.stats"),
          emit: calls

  script:
  """
  set -euo pipefail

  find . -maxdepth 1 -type f -name 'out.*-scattered.vcf.gz' -print \
    | sed -E 's/.*out\\.([0-9]+)-scattered\\.vcf\\.gz/\\1\\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > vcfs.sorted.list

  find . -maxdepth 1 -type f -name 'out.*-scattered.vcf.gz.stats' -print \
    | sed -E 's/.*out\\.([0-9]+)-scattered\\.vcf\\.gz\\.stats/\\1\\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > stats.sorted.list

  vcf_count=\$(wc -l < vcfs.sorted.list)
  stats_count=\$(wc -l < stats.sorted.list)

  if [[ "\${vcf_count}" -ne "${expected_shards}" ]]; then
    echo "ERROR: Expected ${expected_shards} VCF shards but found \${vcf_count}" >&2
    exit 1
  fi

  if [[ "\${stats_count}" -ne "${expected_shards}" ]]; then
    echo "ERROR: Expected ${expected_shards} stats shards but found \${stats_count}" >&2
    exit 1
  fi

  awk '{ print "--INPUT", \$0 }' vcfs.sorted.list > gather.args

  gatk --java-options "-Xmx8g -XX:-UsePerfData" GatherVcfs \
    --arguments_file gather.args \
    --OUTPUT "${sample_id}.merged.vcf.gz"

  if [[ ! -s "${sample_id}.merged.vcf.gz.tbi" ]]; then
    tabix -f -p vcf "${sample_id}.merged.vcf.gz" \
      || gatk IndexFeatureFile \
        --input "${sample_id}.merged.vcf.gz"
  fi

  stats_args=()
  while IFS= read -r stats_file; do
    stats_args+=(-stats "\${stats_file}")
  done < stats.sorted.list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" MergeMutectStats \
    "\${stats_args[@]}" \
    --output "${sample_id}.merged.vcf.gz.stats"

  test -s "${sample_id}.merged.vcf.gz"
  test -s "${sample_id}.merged.vcf.gz.tbi"
  test -s "${sample_id}.merged.vcf.gz.stats"
  """
}


/*
 * Learn one read-orientation artifact model per tumor sample from all Mutect2
 * F1R2 shard archives.
 */
process learn_read_orientation_model {
  tag "${sample_id}"
  label 'process_low'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(sample_id), path(f1r2_files)
    val expected_shards

  output:
    tuple val(sample_id),
          path("${sample_id}.read-orientation-model.tar.gz"),
          emit: artifact_prior

  script:
  """
  set -euo pipefail

  find . -maxdepth 1 -type f -name 'out.*-scattered.f1r2.tar.gz' -print \
    | sed -E 's/.*out\\.([0-9]+)-scattered\\.f1r2\\.tar\\.gz/\\1\\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > f1r2.sorted.list

  f1r2_count=\$(wc -l < f1r2.sorted.list)

  if [[ "\${f1r2_count}" -ne "${expected_shards}" ]]; then
    echo "ERROR: Expected ${expected_shards} F1R2 shards but found \${f1r2_count}" >&2
    exit 1
  fi

  f1r2_args=()
  while IFS= read -r f1r2_file; do
    f1r2_args+=(-I "\${f1r2_file}")
  done < f1r2.sorted.list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" LearnReadOrientationModel \
    "\${f1r2_args[@]}" \
    --output "${sample_id}.read-orientation-model.tar.gz"

  test -s "${sample_id}.read-orientation-model.tar.gz"
  """
}


/*
 * Gather scattered tumor and normal GetPileupSummaries tables.
 */
process gather_pileup_summaries {
  tag "${sample_id}"
  label 'process_low'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(sample_id), path(tumor_tables), path(normal_tables)
    path ref_dict
    val expected_shards

  output:
    tuple val(sample_id),
          path("${sample_id}.tumor.pileups.table"),
          path("${sample_id}.normal.pileups.table"),
          emit: tables

  script:
  """
  set -euo pipefail

  find . -maxdepth 1 -type f -name 'out.*-scattered.tumor.pileups.table' -print \
    | sed -E 's/.*out\\.([0-9]+)-scattered\\.tumor\\.pileups\\.table/\\1\\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > tumor.sorted.list

  find . -maxdepth 1 -type f -name 'out.*-scattered.normal.pileups.table' -print \
    | sed -E 's/.*out\\.([0-9]+)-scattered\\.normal\\.pileups\\.table/\\1\\t&/' \
    | sort -k1,1n \
    | cut -f2- \
    > normal.sorted.list

  tumor_count=\$(wc -l < tumor.sorted.list)
  normal_count=\$(wc -l < normal.sorted.list)

  if [[ "\${tumor_count}" -ne "${expected_shards}" ]]; then
    echo "ERROR: Expected ${expected_shards} tumor pileup shards but found \${tumor_count}" >&2
    exit 1
  fi

  if [[ "\${normal_count}" -ne "${expected_shards}" ]]; then
    echo "ERROR: Expected ${expected_shards} normal pileup shards but found \${normal_count}" >&2
    exit 1
  fi

    tumor_args=()
  while IFS= read -r table_file; do
    tumor_args+=(-I "\${table_file}")
  done < tumor.sorted.list

  normal_args=()
  while IFS= read -r table_file; do
    normal_args+=(-I "\${table_file}")
  done < normal.sorted.list

  gatk --java-options "-Xmx8g -XX:-UsePerfData" GatherPileupSummaries \
    "\${tumor_args[@]}" \
    --sequence-dictionary "${ref_dict}" \
    -O "${sample_id}.tumor.pileups.table"

  gatk --java-options "-Xmx8g -XX:-UsePerfData" GatherPileupSummaries \
    "\${normal_args[@]}" \
    --sequence-dictionary "${ref_dict}" \
    -O "${sample_id}.normal.pileups.table"

  test -s "${sample_id}.tumor.pileups.table"
  test -s "${sample_id}.normal.pileups.table"
  """
}


/*
 * Estimate cross-sample contamination and produce tumor segmentation.
 */
process calculate_contamination {
  tag "${sample_id}"
  label 'process_low'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(sample_id),
          path(tumor_table),
          path(normal_table)

  output:
    tuple val(sample_id),
          path("${sample_id}.contamination.table"),
          path("${sample_id}.segments.table"),
          emit: results

  script:
  """
  set -euo pipefail

  gatk --java-options "-Xmx8g -XX:-UsePerfData" CalculateContamination \
    --input "${tumor_table}" \
    --matched-normal "${normal_table}" \
    --output "${sample_id}.contamination.table" \
    --tumor-segmentation "${sample_id}.segments.table"

  test -s "${sample_id}.contamination.table"
  test -s "${sample_id}.segments.table"
  """
}


/*
 * Apply Mutect stats, read-orientation, contamination, and segmentation filters.
 */
process filter_mutect_calls {
  tag "${sample_id}"
  label 'process_medium'
  container "${params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'}"

  input:
    tuple val(sample_id),
          path(vcf),
          path(vcf_tbi),
          path(stats),
          path(orientation_model),
          path(contamination_table),
          path(segmentation_table)

    path ref_fasta
    path ref_fai
    path ref_dict

  output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"),        emit: vcf
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz.tbi"),    emit: tbi
    tuple val(sample_id), path("${sample_id}.filteringStats.tsv"),     emit: stats

  script:
  """
  set -euo pipefail

  gatk --java-options "-Xmx8g -XX:-UsePerfData" FilterMutectCalls \
    --reference "${ref_fasta}" \
    --variant "${vcf}" \
    --stats "${stats}" \
    --orientation-bias-artifact-priors "${orientation_model}" \
    --contamination-table "${contamination_table}" \
    --tumor-segmentation "${segmentation_table}" \
    --output "${sample_id}.filtered.vcf.gz" \
    --filtering-stats "${sample_id}.filteringStats.tsv" \
    --tmp-dir .

  if [[ ! -s "${sample_id}.filtered.vcf.gz.tbi" ]]; then
    tabix -f -p vcf "${sample_id}.filtered.vcf.gz" \
      || gatk IndexFeatureFile \
        --input "${sample_id}.filtered.vcf.gz"
  fi

  test -s "${sample_id}.filtered.vcf.gz"
  test -s "${sample_id}.filtered.vcf.gz.tbi"
  test -s "${sample_id}.filteringStats.tsv"
  """
}


/*
 * Optional post hoc cohort VCF merge.
 * This is not joint somatic calling.
 */
process merge_all_sample_vcfs {
  label 'process_medium'
  container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

  input:
    path vcfs
    path tbis

  output:
    path "all_samples.merged.vcf.gz",      emit: vcf
    path "all_samples.merged.vcf.gz.tbi",  emit: tbi

  script:
  """
  set -euo pipefail

  rm -f tumor_vcfs.list
  mkdir -p tumor_only_vcfs

  for f in \$(find . -maxdepth 1 -type f -name '*.filtered.vcf.gz' -print | sort); do
    tumor_samples=\$(bcftools query -l "\${f}" \
      | grep -viE 'PBMC|NORMAL|BLOOD' \
      || true)

    if [[ -z "\${tumor_samples}" ]]; then
      echo "ERROR: No tumor sample found in \${f}" >&2
      exit 1
    fi

    keep_file="\$(basename "\${f}" .vcf.gz).keep_samples.txt"
    out="tumor_only_vcfs/\$(basename "\${f}" .vcf.gz).tumor_only.vcf.gz"

    printf '%s\\n' "\${tumor_samples}" > "\${keep_file}"

    bcftools view \
      --samples-file "\${keep_file}" \
      --output-type z \
      --output "\${out}" \
      "\${f}"

    tabix -f -p vcf "\${out}"
    echo "\${out}" >> tumor_vcfs.list
  done

  bcftools merge \
    --force-samples \
    --output-type z \
    --output all_samples.merged.vcf.gz \
    --file-list tumor_vcfs.list

  tabix -f -p vcf all_samples.merged.vcf.gz

  bcftools query -l all_samples.merged.vcf.gz
  """
}


workflow {

  /*
   * This implementation follows the paired tumor-normal nf-core logic.
   * A matched normal is required for each run.
   */
  params.mutect_runs.each { run ->
    if( !run.normal_reads || !run.normal_reads_index ) {
      error "Run ${run.output_prefix} is missing normal_reads or normal_reads_index. This full filtering workflow requires a matched normal."
    }
  }

  ref_fasta_ch = file(params.ref_fasta, checkIfExists: true)
  ref_fai_ch = file(params.ref_fai, checkIfExists: true)
  ref_dict_ch = file(params.ref_dict, checkIfExists: true)

  germline_resource_ch = file(
    params.germline_resource,
    checkIfExists: true
  )

  germline_resource_index_ch = file(
    params.germline_resource_index ?: "${params.germline_resource}.tbi",
    checkIfExists: true
  )

  panel_of_normals_ch = (
    params.panel_of_normals
      ? file(params.panel_of_normals, checkIfExists: true)
      : file(NO_PON_VCF_PATH, checkIfExists: true)
  )

  panel_of_normals_index_ch = (
    params.panel_of_normals_index
      ? file(params.panel_of_normals_index, checkIfExists: true)
      : file(NO_PON_TBI_PATH, checkIfExists: true)
  )

  alleles_vcf_ch = (
    params.force_call_file
      ? file(params.force_call_file, checkIfExists: true)
      : file(NO_ALLELES_VCF_PATH, checkIfExists: true)
  )

  alleles_vcf_tbi_ch = (
    params.force_call_file_index
      ? file(params.force_call_file_index, checkIfExists: true)
      : file(NO_ALLELES_TBI_PATH, checkIfExists: true)
  )

  expected_shards = params.scatter_count as int
  batch_size = (params.extract_batch_size ?: 10) as int

  /*
   * 1. Split calling intervals.
   */
  interval_res = split_intervals(
    ref_fasta_ch,
    ref_fai_ch,
    ref_dict_ch,
    file(params.intervals, checkIfExists: true),
    expected_shards
  )

  /*
   * 2. Batch interval shards.
   */
  interval_batches_ch = interval_res.interval_shards
    .flatten()
    .toSortedList { a, b -> a.name <=> b.name }
    .flatMap { intervals ->
      intervals.collate(batch_size)
    }

  /*
   * 3. Validate that every tumor run points to the same normal BAM.
   */
  shared_normal_paths = params.mutect_runs
    .collect { run -> run.normal_reads }
    .unique()

  shared_normal_index_paths = params.mutect_runs
    .collect { run -> run.normal_reads_index }
    .unique()

  if( shared_normal_paths.size() != 1 ) {
    error "Expected exactly one shared normal BAM, but found: ${shared_normal_paths}"
  }

  if( shared_normal_index_paths.size() != 1 ) {
    error "Expected exactly one shared normal BAM index, but found: ${shared_normal_index_paths}"
  }

  /*
   * 4. Subset every tumor sample independently.
   */
  tumor_batches_ch = interval_batches_ch
    .flatMap { interval_batch ->
      params.mutect_runs.collect { run ->
        tuple(
          [id: run.output_prefix],
          file(run.tumor_reads, checkIfExists: true),
          file(run.tumor_reads_index, checkIfExists: true),
          interval_batch
        )
      }
    }

  tumor_subset_res = subset_tumor_all_shards(tumor_batches_ch)

  /*
   * 5. Subset the single shared normal once.
   */
  shared_normal_batches_ch = interval_batches_ch
    .map { interval_batch ->
      tuple(
        file(shared_normal_paths[0], checkIfExists: true),
        file(shared_normal_index_paths[0], checkIfExists: true),
        interval_batch
      )
    }

  normal_subset_res = subset_shared_normal_all_shards(
    shared_normal_batches_ch
  )

  /*
   * Flatten tumor output into one tuple per sample and shard.
   */
  tumor_shards_ch = tumor_subset_res.shards
    .flatMap { sample_id, tumor_bams, tumor_bais, intervals ->
      def bam_list = tumor_bams instanceof List ? tumor_bams : [tumor_bams]
      def bai_list = tumor_bais instanceof List ? tumor_bais : [tumor_bais]
      def interval_list = intervals instanceof List ? intervals : [intervals]

      bam_list.collect { tumor_bam ->
        def matcher = tumor_bam.name =~ /^(.+)\.tumor\.(\d+-scattered)\.bam$/

        if( !matcher.matches() )
          error "Could not parse tumor BAM shard name: ${tumor_bam.name}"

        def shard_id = matcher[0][2]

        def tumor_bai = bai_list.find {
          it.name == "${sample_id}.tumor.${shard_id}.bam.bai"
        }

        def interval = interval_list.find {
          it.name == "${shard_id}.interval_list"
        }

        if( tumor_bai == null )
          error "Missing tumor BAI for ${sample_id}, shard ${shard_id}"

        if( interval == null )
          error "Missing interval list for ${sample_id}, shard ${shard_id}"

        tuple(shard_id, sample_id, interval, tumor_bam, tumor_bai)
      }
    }

  /*
   * Flatten shared-normal output into one tuple per shard.
   */
  normal_shards_ch = normal_subset_res.shards
    .flatMap { normal_bams, normal_bais ->
      def bam_list = normal_bams instanceof List ? normal_bams : [normal_bams]
      def bai_list = normal_bais instanceof List ? normal_bais : [normal_bais]

      bam_list.collect { normal_bam ->
        def matcher = normal_bam.name =~ /^shared_normal\.normal\.(\d+-scattered)\.bam$/

        if( !matcher.matches() )
          error "Could not parse normal BAM shard name: ${normal_bam.name}"

        def shard_id = matcher[0][1]

        def normal_bai = bai_list.find {
          it.name == "shared_normal.normal.${shard_id}.bam.bai"
        }

        if( normal_bai == null )
          error "Missing shared-normal BAI for shard ${shard_id}"

        tuple(shard_id, normal_bam, normal_bai)
      }
    }

  /*
   * Reuse each shared-normal shard for every tumor with the same shard ID.
   */
  mutect_inputs_ch = tumor_shards_ch
    .combine(normal_shards_ch, by: 0)
    .map {
      shard_id,
      sample_id,
      interval,
      tumor_bam,
      tumor_bai,
      normal_bam,
      normal_bai ->

      tuple(
        sample_id,
        interval,
        tumor_bam,
        tumor_bai,
        normal_bam,
        normal_bai,
        ref_fasta_ch,
        ref_fai_ch,
        ref_dict_ch,
        germline_resource_ch,
        germline_resource_index_ch,
        panel_of_normals_ch,
        panel_of_normals_index_ch,
        alleles_vcf_ch,
        alleles_vcf_tbi_ch,
        params.m2_extra_args ?: ''
      )
    }

  mutect_inputs_ch.view { x ->
    "MUTECT FANOUT: sample=${x[0]}, interval=${x[1].name}, tumor=${x[2].name}, normal=${x[4].name}"
  }

  mutect_res = mutect_wrapper(mutect_inputs_ch)

  /*
   * 5. Gather VCFs/stats and learn the orientation model.
   */
  grouped_vcfs_ch = mutect_res.vcf
    .groupTuple(size: expected_shards)

  grouped_stats_ch = mutect_res.stats
    .groupTuple(size: expected_shards)

  grouped_mutect_outputs_ch = grouped_vcfs_ch
    .join(
      grouped_stats_ch,
      failOnDuplicate: true,
      failOnMismatch: true
    )
    .map { sample_id, vcfs, stats ->
      tuple(sample_id, vcfs, stats)
    }

  gathered_mutect_res = gather_mutect_outputs(
    grouped_mutect_outputs_ch,
    expected_shards
  )

  grouped_f1r2_ch = mutect_res.f1r2
    .groupTuple(size: expected_shards)

  orientation_res = learn_read_orientation_model(
    grouped_f1r2_ch,
    expected_shards
  )

  /*
   * 6. Gather tumor/normal pileup summaries and calculate contamination.
   */
  grouped_tumor_pileups_ch = mutect_res.tumor_pileup
    .groupTuple(size: expected_shards)

  grouped_normal_pileups_ch = mutect_res.normal_pileup
    .groupTuple(size: expected_shards)

  grouped_pileups_ch = grouped_tumor_pileups_ch
    .join(
      grouped_normal_pileups_ch,
      failOnDuplicate: true,
      failOnMismatch: true
    )
    .map { sample_id, tumor_tables, normal_tables ->
      tuple(sample_id, tumor_tables, normal_tables)
    }

  gathered_pileup_res = gather_pileup_summaries(
    grouped_pileups_ch,
    ref_dict_ch,
    expected_shards
  )

  contamination_res = calculate_contamination(
    gathered_pileup_res.tables
  )

  /*
   * 7. Join all per-sample filtering inputs.
   */
  filter_inputs_ch = gathered_mutect_res.calls
    .join(
      orientation_res.artifact_prior,
      failOnDuplicate: true,
      failOnMismatch: true
    )
    .join(
      contamination_res.results,
      failOnDuplicate: true,
      failOnMismatch: true
    )
    .map {
      sample_id,
      vcf,
      vcf_tbi,
      stats,
      orientation_model,
      contamination_table,
      segmentation_table ->

      tuple(
        sample_id,
        vcf,
        vcf_tbi,
        stats,
        orientation_model,
        contamination_table,
        segmentation_table
      )
    }

  filtered_mutect_res = filter_mutect_calls(
    filter_inputs_ch,
    ref_fasta_ch,
    ref_fai_ch,
    ref_dict_ch
  )

  filtered_mutect_res.vcf.view { sample_id, vcf ->
    "FILTERED MUTECT VCF: sample=${sample_id}, vcf=${vcf}"
  }

  /*
   * 8. Optional post hoc merge across tumor samples.
   */
  if( params.merge_all_sample_vcfs ) {
    all_sample_vcfs_ch = filtered_mutect_res.vcf
      .map { sample_id, vcf -> vcf }
      .collect()

    all_sample_tbis_ch = filtered_mutect_res.tbi
      .map { sample_id, tbi -> tbi }
      .collect()

    merge_all_sample_vcfs(
      all_sample_vcfs_ch,
      all_sample_tbis_ch
    )
  }
  else {
    log.info "Skipping cross-sample VCF merge because params.merge_all_sample_vcfs=false"
  }
}
