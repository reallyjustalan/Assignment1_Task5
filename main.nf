nextflow.enable.dsl=2

process FASTQ_QC {
    tag "$sample"
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/qc/raw", mode: 'symlink'

    input:
    tuple val(sample), val(protocol), path(reads)
    path qc_script

    output:
    tuple val(sample), path("${sample}.raw.json")

    script:
    """
    python3 '${qc_script}' --sample '${sample}' \
        --fastq '${reads}' --output '${sample}.raw.json'
    """
}

process REFERENCE_INDEX {
    cpus 1
    memory '2 GB'

    input:
    path genome

    output:
    path "${genome}.fai"

    script:
    """
    samtools faidx '${genome}'
    """
}

process ALIGN_SORT {
    tag "$sample ($protocol)"
    cpus 6
    memory '24 GB'
    maxForks 1
    publishDir "${params.outdir}/bam", mode: 'symlink'

    input:
    tuple val(sample), val(protocol), path(reads)
    path genome

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    path "${sample}.minimap2.log", emit: log
    path "${sample}.versions.txt", emit: versions

    script:
    def flags = protocol == 'directRNA' ? '-ax splice -uf -k14' : '-ax splice'
    """
    minimap2 ${flags} -t 4 '${genome}' '${reads}' \
        2> '${sample}.minimap2.log' \
        | samtools sort -@ 1 -m 1G -T sort_tmp -o '${sample}.bam' -
    samtools index '${sample}.bam'
    samtools quickcheck -v '${sample}.bam'
    { minimap2 --version; samtools --version; } > '${sample}.versions.txt'
    """
}

process ALIGNMENT_QC {
    tag "$sample"
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/qc/alignment", mode: 'symlink'

    input:
    tuple val(sample), path(bam), path(bai), path(raw_stats)
    path qc_script

    output:
    tuple val(sample), path("${sample}.qc.tsv"), emit: metrics
    path "${sample}.flagstat.json", emit: flagstat
    path "${sample}.mapq20.txt", emit: mapq

    script:
    """
    samtools quickcheck -v '${bam}'
    samtools flagstat -O json '${bam}' > '${sample}.flagstat.json'
    samtools view -c -F 2308 -q 20 '${bam}' > '${sample}.mapq20.txt'
    python3 '${qc_script}' \
        --raw '${raw_stats}' --flagstat '${sample}.flagstat.json' \
        --mapq '${sample}.mapq20.txt' --min-reads ${params.min_reads} \
        --min-mapped-pct ${params.min_mapped_pct} --output '${sample}.qc.tsv'
    """
}

process COLLECT_QC {
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/qc", mode: 'symlink'

    input:
    path tables
    path merge_script

    output:
    path 'qc_summary.tsv'

    script:
    """
    python3 '${merge_script}' ${tables.join(' ')} > qc_summary.tsv
    """
}

process BAMBU {
    tag "$mode"
    cpus 4
    memory '40 GB'
    publishDir "${params.outdir}/bambu/${mode}", mode: 'symlink'

    input:
    tuple path(bams), path(indices)
    path genome
    path genome_index
    path annotation
    val mode
    path qc_summary
    path bambu_script

    output:
    path 'bambu_output'

    script:
    def annotation_arg = mode == 'annotated' ? annotation : '-'
    """
    Rscript '${bambu_script}' '${mode}' '${genome}' \
        '${annotation_arg}' ${task.cpus} ${bams.join(' ')}
    """
}

workflow {
    if (!(params.mode in ['annotated', 'unannotated']))
        error "--mode must be annotated or unannotated"
    if (!params.samples || !params.genome)
        error 'Provide --samples and --genome'
    if (params.mode == 'annotated' && !params.gtf)
        error 'Annotated mode requires --gtf'
    if (params.min_reads < 0 || params.min_mapped_pct < 0 || params.min_mapped_pct > 100)
        error 'Invalid QC thresholds'

    samples = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !(row.sample ==~ /[A-Za-z0-9_-]+/))
                error 'Sample IDs must contain only letters, numbers, underscores or hyphens'
            if (!(row.protocol in ['directRNA', 'cDNA']))
                error "Unknown protocol for ${row.sample}: ${row.protocol}"
            if (!row.fastq || !row.fastq.startsWith('/'))
                error 'Use absolute FASTQ paths in the sample sheet'
            tuple(row.sample, row.protocol, file(row.fastq, checkIfExists: true))
        }
        .ifEmpty { error 'Empty sample sheet' }
        .collect(flat: false)
        .flatMap { rows ->
            if (!rows) error 'Empty sample sheet'
            if (rows.collect { it[0] }.unique(false).size() != rows.size())
                error 'Duplicate sample IDs'
            rows
        }

    genome = Channel.value(file(params.genome, checkIfExists: true))
    annotation = Channel.value(params.mode == 'annotated' ? file(params.gtf, checkIfExists: true) : [])

    FASTQ_QC(samples, file("${projectDir}/bin/fastq_qc.py"))
    REFERENCE_INDEX(genome)
    ALIGN_SORT(samples, genome)
    ALIGNMENT_QC(ALIGN_SORT.out.bam.join(FASTQ_QC.out), file("${projectDir}/bin/summarise_qc.py"))
    COLLECT_QC(ALIGNMENT_QC.out.metrics.map { id, table -> table }.collect().map { it.sort { a, b -> a.name <=> b.name } }, file("${projectDir}/bin/merge_qc.py"))

    // Stable sample order makes Bambu inputs reproducible and cache-friendly.
    all_bams = ALIGN_SORT.out.bam.collect(flat: false).map { rows ->
        def ordered = rows.sort { a, b -> a[0] <=> b[0] }
        tuple(ordered.collect { it[1] }, ordered.collect { it[2] })
    }
    BAMBU(all_bams, genome, REFERENCE_INDEX.out, annotation, params.mode, COLLECT_QC.out, file("${projectDir}/bin/run_bambu.R"))
}

workflow.onComplete {
    log.info "Workflow ${workflow.runName}: ${workflow.success ? 'SUCCESS' : 'FAILED'} | mode=${params.mode} | elapsed=${workflow.duration}"
}
