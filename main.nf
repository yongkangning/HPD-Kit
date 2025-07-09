#!/usr/bin/env nextflow
// nextflow run main.nf --metadataFile=samples_infor.csv --outdir=result -profile singularity

process QC_AND_DEHOST {
    // debug true
    container = params.base_img
    containerOptions = params.container_options
    maxForks params.qc_and_dehost_max_run
    input:
        tuple val(sampleId), path(reads), val(host_genome)
    output:
        tuple val(sampleId), path("${sampleId}.dehost*.fq.gz")

    script:
    paired_end = reads.size() == 2 ? "yes" : "no"
    """
        # Check whether the gzip file is complete
        pigz -t \$(readlink -f ${reads})

        # read length filter
        read_length=\$(seqkit stats ${reads[0]} | tail -n+2 | awk '{print \$8}' | sed 's/,//g')
        read_length_threshold=\$(python -c "print(max(${params.min_read_length}, round(\${read_length} * ${params.read_coverage_threshold})))")

        echo "read_length: \${read_length}  read_length_threshold: \${read_length_threshold}"

        # fastp quality control
        if [ "${paired_end}" == "yes" ]; then
            # paired-end
            fastp -w 8 \
            -i ${reads[0]} -o ${sampleId}.clean.R1.fq.gz \
            -I ${reads[1]} -O ${sampleId}.clean.R2.fq.gz \
            -j ${sampleId}.json \
            -h ${sampleId}.html \
            -q 20 -u 40 -n 10 -l \${read_length_threshold}

            # Dehost
            if [ "${host_genome}" == "no" ]; then
                ln -s ${sampleId}.clean.R1.fq.gz ${sampleId}.dehost.R1.fq.gz
                ln -s ${sampleId}.clean.R2.fq.gz ${sampleId}.dehost.R2.fq.gz
            else
                bowtie2 --very-sensitive -p ${params.bwt2_dehost_thread} -x ${host_genome} -1 ${sampleId}.clean.R1.fq.gz -2 ${sampleId}.clean.R2.fq.gz --un-conc-gz ${sampleId}.dehost.R%.fq.gz > /dev/null
            fi

        else
            # single-end
            fastp -w 8 \
            -i ${reads[0]} -o ${sampleId}.clean.fq.gz \
            -j ${sampleId}.json \
            -h ${sampleId}.html \
            -q 20 -u 40 -n 10 -l \${read_length_threshold}

            # Dehost
            if [ "${host_genome}" == "no" ]; then
                ln -s ${sampleId}.clean.fq.gz ${sampleId}.dehost.fq.gz
            else
                bowtie2 --very-sensitive -p ${params.bwt2_dehost_thread} -x ${host_genome} -U ${sampleId}.clean.fq.gz --un-gz ${sampleId}.dehost.fq.gz > /dev/null
            fi
        fi
    """
}

process KRAKEN2_IDENTIFICATION {
    // debug true
    container = params.base_img
    containerOptions = params.container_options
    maxForks params.kraken2_identification_max_run
    input:
        tuple val(sampleId), path(reads)
        each pathogen_type
    output:
        tuple val(sampleId), val(pathogen_type), path("${sampleId}.${pathogen_type}.k2.report.txt"), path(reads)

    script:
    paired_end = reads.size() == 2 ? "yes" : "no"
    """
        # kraken2 identification
        if [ "${paired_end}" == "yes" ]; then
            # paired-end
            kraken2 \
              --db ${params.pathogen_db_path}/${pathogen_type}/k2 \
              --threads ${params.k2_classify_thread} \
              --minimum-hit-groups ${params.k2_minimum_hit_groups} \
              --report-minimizer-data \
              --paired ${sampleId}.dehost.R1.fq.gz ${sampleId}.dehost.R2.fq.gz \
              --report ${sampleId}.${pathogen_type}.k2.report.txt \
              --output ${sampleId}.${pathogen_type}.k2.output.txt
        else
            # single-end
            kraken2 \
              --db ${params.pathogen_db_path}/${pathogen_type}/k2 \
              --threads ${params.k2_classify_thread} \
              --minimum-hit-groups ${params.k2_minimum_hit_groups} \
              --report-minimizer-data \
              ${sampleId}.dehost.fq.gz \
              --report ${sampleId}.${pathogen_type}.k2.report.txt \
              --output ${sampleId}.${pathogen_type}.k2.output.txt
        fi
    """
}


process KRAKEN2_RESULT_FILTER {
    // debug true
    container = params.base_img
    containerOptions = params.container_options
    maxForks params.kraken2_result_filter_max_run
    publishDir "${params.outDir}/kraken2", mode:'copy', pattern: '*.{txt,xls}'
    input:
        tuple val(sampleId), val(pathogen_type), path("${sampleId}.${pathogen_type}.k2.report.txt"), path(reads)
    output:
        tuple val(sampleId), val(pathogen_type), env(passed_filter_count), path("${sampleId}.k2.${pathogen_type}.raw.xls"), path("${sampleId}.k2.${pathogen_type}.filtered.xls"), path("${sampleId}.k2.${pathogen_type}.unpassed.xls"), path(reads)

    script:
    """
        # filter kraken2 result
        kraken2_report_processed.R \
          ${sampleId} \
          ${sampleId}.${pathogen_type}.k2.report.txt \
          ${params.min_k2_reads} \
          ${params.min_relative_abundance} \
          ${params.min_unique_kmers} \
          ${params.min_unique_kmers_fold} \
          ${params.min_unique_kmers_rate} \
          ${pathogen_type} \
          ${params.pathogen_db_path}/${pathogen_type}/${pathogen_type}.genome.infor.xls

        passed_filter_count=`tail -n+2 "${sampleId}.k2.${pathogen_type}.filtered.xls" | wc -l`
    """
}

process BWT2_BLAST_VERIFICATION {
    // debug true
    container = params.base_img
    containerOptions = params.container_options
    maxForks params.bwt2_blast_verification_max_run
    publishDir "${params.outDir}/coverage_plot", mode:'copy', pattern: '*_cov.png'
    publishDir "${params.outDir}/identified_pathogens", mode:'copy', pattern: '*.xls'
    publishDir "${params.outDir}/bwt2", mode:'copy', pattern: '*.markdup.bam*'
    input:
        tuple val(sampleId), val(pathogen_type), path(k2_filter_result), path(dehost_read)
    output:
        tuple val(sampleId), path("${sampleId}.${pathogen_type}.final.xls"), emit: final_xls
        path "*.markdup.bam*"
        // Optional output file
        path "*.png" optional true

    script:
    paired_end = dehost_read.size() == 2 ? "yes" : "no"
    """

        if [ "${paired_end}" == "yes" ]; then
            # paired-end
            tail -n+2  ${k2_filter_result} | rush  -k -j ${params.rush_parallel_number} \
              '''
                bwt2_blast_verification_paired.sh \
                ${sampleId} \
                {2} \
                ${pathogen_type} \
                ${params.pathogen_db_path} \
                ${dehost_read[0]} \
                ${dehost_read[1]} \
                ${params.min_map_quality} \
                ${params.plot_coverage} \
                $baseDir/bin \
                ${params.bwt2_pathogen_thread} \
                ${params.blast_pathogen_thread}
              '''
        else
            # single-end
            tail -n+2  ${k2_filter_result} | rush  -k -j ${params.rush_parallel_number} \
              '''
                bwt2_blast_verification_single.sh \
                ${sampleId} \
                {2} \
                ${pathogen_type} \
                ${params.pathogen_db_path} \
                ${dehost_read[0]} \
                ${params.min_map_quality} \
                ${params.plot_coverage} \
                $baseDir/bin \
                ${params.bwt2_pathogen_thread} \
                ${params.blast_pathogen_thread}
              '''
        fi

        # merge the result
        cat ${sampleId}.*.bwt2.txt > ${sampleId}.${pathogen_type}.bwt2.xls
        sed -i '1i\\sample_id\\ttaxid\\taccession\\tmapped_chrom_count\\tbwt2_reads\\tgenome_base_count\\tmapped_base\\tcoverage_positions\\tavg_depth\\tavg_coverage' ${sampleId}.${pathogen_type}.bwt2.xls
        cat ${sampleId}.*.blast.res.txt | grep -v "sample_id" > ${sampleId}.${pathogen_type}.blast.xls
        sed -i '1i\\sample_id\\ttaxid\\tblast_reads\\tavg_pident\\tmin_pident\\tmax_pident' ${sampleId}.${pathogen_type}.blast.xls

        # sample total reads
        sample_reads=\$(seqkit stats ${dehost_read[0]} | tail -n+2 | awk '{print \$4}' | sed 's/,//g')

        merge_result.R \
          ${sampleId} \
          ${pathogen_type} \
          ${sampleId}.k2.${pathogen_type}.filtered.xls \
          ${sampleId}.${pathogen_type}.bwt2.xls \
          ${sampleId}.${pathogen_type}.blast.xls \
          ${params.pathogen_db_path}/${pathogen_type}/${pathogen_type}.genome.infor.xls \
          ${params.min_unique_reads} \
          ${params.min_normalized_pathogen_abundance} \
          ${params.min_base_coverage} \
          ${params.min_sequence_coverage} \
          \${sample_reads}
    """
}

process MERGE_PATHOGENS_RESULT {
    // debug true
    container = params.base_img
    containerOptions = params.container_options
    maxForks params.merge_pathogens_max_run
    publishDir "${params.outDir}/identified_pathogens", mode:'copy'
    input:
        tuple val(sampleId), path(xls_result)
    output:
        path "${sampleId}.pathogens.final.xls"
    // when:
    //    params.merge_pathogens_type == 1
    script:
    """
        head -1 `find ./ -name "${sampleId}*.final.xls" | head -1` > ${sampleId}.final.tmp
        cat `find ./ -name "${sampleId}*.final.xls"` | { grep -v "^taxid" || true; } >> ${sampleId}.final.tmp
        sort_result.R ${sampleId}.final.tmp ${sampleId}.pathogens.final.xls NPAS reverse
    """
}


process CONTRAST_PATHOGENS_RESULT {
    // debug true
    container = params.base_img
    containerOptions = params.container_options
    publishDir "${params.outDir}/identified_pathogens", mode:'copy'
    input:
        path xls_result
    output:
        path "*_vs_control*.xls" optional true
    when:
        params.comparisonFile != ""
    shell:
    '''
        cat !{params.comparisonFile} | while read -r treat control
        do
            cat /dev/null > ${treat}_vs_control.tsv
            for sp in `echo ${control} | tr ',' ' '`
            do
                if [ -f "${sp}.pathogens.final.xls" ]; then
                    tail -n+2 "${sp}.pathogens.final.xls" | cut -f1,13 >> ${treat}_vs_control.tsv
                fi
            done
            sed -i '1i taxid\tNPA_control' ${treat}_vs_control.tsv
            rows_count=`cat ${treat}_vs_control.tsv | grep -v taxid | wc -l`
            if [ ${rows_count} -ne 0 ]; then
                adjusted_NPAS.R ${treat}.pathogens.final.xls ${treat}_vs_control.tsv !{params.control_samples_npa_method}
            fi
        done
    '''
}

workflow {
    // metaChanel = Channel.fromPath(params.metadataFile).splitCsv(sep:',', header: true).map{row -> ["${row.SampleId}", "${row.Read1}", "${row.Read2}", "${row.Host}"]}
    fastq_list = params.fastq1 ? (params.fastq2 ? [params.sample_id, [params.fastq1, params.fastq2], params.host_fasta_index] : [params.sample_id, [params.fastq1], params.host_fasta_index]) : []
    // fastqChanel = Channel.of(fastq_list)
    // fastqChanel.view()
    metaChanel = params.metadataFile ? Channel.fromPath(params.metadataFile).splitCsv(sep:',', header: true).map{row -> ["${row.SampleId}", "${row.Read2}" ? ["${row.Read1}","${row.Read2}"] : ["${row.Read1}"], "${row.Host}"]} : Channel.of(fastq_list)
    // metaChanel.view()
    // metaChanel = Channel.fromPath(params.metadataFile).splitCsv(sep:',', header: true).map{row -> ["${row.SampleId}", "${row.Read2}" ? ["${row.Read1}","${row.Read2}"] : ["${row.Read1}"], "${row.Host}"]}
    pathogenType =  Channel.of(params.pathogenTypes.split(','))
    ch_qc_dehost = QC_AND_DEHOST(metaChanel)
    ch_k2 = KRAKEN2_IDENTIFICATION(ch_qc_dehost, pathogenType) | KRAKEN2_RESULT_FILTER
    // only k2 filter result greater 0 will be verification
    ch_k2_verification = ch_k2.filter{ it[2].toInteger() > 0 }.map{ [it[0], it[1], it[4], it[6]] }
    BWT2_BLAST_VERIFICATION(ch_k2_verification)
    ch_merge = BWT2_BLAST_VERIFICATION.out.final_xls.groupTuple()
    ch_final = MERGE_PATHOGENS_RESULT(ch_merge)
    ch_final.collect() | CONTRAST_PATHOGENS_RESULT
}