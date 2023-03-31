process BRACKEN {

  tag { "abundance estimation in ${sample_id}" }

  if (params.output_dir != "") { publishDir(params.output_dir, mode:'copy') }

  input:
  tuple val(sample_id), path(kraken_report)
  path(brackendb)
  val(taxon)


  output:
  path("*_bracken"),                                        emit: bracken_out_ch
  tuple val(sample_id), path("*bracken.kraken.report.txt"), emit: bracken_kraken_report_ch
  tuple val(sample_id), path("*bracken.log.file.txt"),      emit: bracken_log_ch
  path "versions.yml",                                      emit: versions


  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${sample_id}"
  def VERSION = '2.8'
  """

  bracken ${args} -d ${brackendb} -i ${kraken_report} -o ${prefix}.${taxon}_bracken -r 100 -l ${taxon} -t ${task.cpus} -w ${prefix}.${taxon}.bracken.kraken.report.txt > ${prefix}.bracken.log.file.txt

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${VERSION}
  END_VERSIONS

  """

  }


process EXTRACT_TAXON_SPECIFIC_INFO {
    //publishDir "${params.output_dir}", mode:'copy'
    tag "extract $taxon from ${sample_id}.braken.kraken.report.txt"


    input:
    tuple val(sample_id), path(bracken_kraken_report)
    val(taxon)


    output:
    path("*.${taxon}.bracken.kraken.txt"), emit: taxon_bracken_kraken_ch


    script:
    """

    echo "percentage\tcladeReads\ttaxonReads\ttaxRank\ttaxID\tspecies" > ${sample_id}.${taxon}.bracken.kraken.txt

    grep "\\s${taxon}\\s" ${bracken_kraken_report} >> ${sample_id}.${taxon}.bracken.kraken.txt

    """
}

process COMBINE_BRACKEN_KRAKEN_REPORT_FROM_TAXA {
    publishDir "${params.output_dir}", mode:'copy'
    tag "combine bracken.kraken.report.txt for $taxon"


    input:
    path(bracken_kraken_taxon_report_files)
    val(taxon)
    val(date)


    output:
    path("combined_bracken_kraken_report_${taxon}_${date}.txt"), emit: bracken_kraken_comb_ch


    script:
    """
    BRACKEN_KRAKEN_TAXON_REPORT_FILES=(${bracken_kraken_taxon_report_files})

    for index in \${!BRACKEN_KRAKEN_TAXON_REPORT_FILES[@]}; do
    BRACKEN_KRAKEN_TAXON_REPORT_FILE=\${BRACKEN_KRAKEN_TAXON_REPORT_FILES[\$index]}
    sample_id=\${BRACKEN_KRAKEN_TAXON_REPORT_FILE%.${taxon}.bracken.kraken.txt}

    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "samplename\t\$(head -1 \${BRACKEN_KRAKEN_TAXON_REPORT_FILE})" >> combined_bracken_kraken_report_${taxon}_${date}.txt
    fi
    awk -F '\\t' 'FNR>=2 { print FILENAME, \$0 }' \${BRACKEN_KRAKEN_TAXON_REPORT_FILE} |  sed 's/\\.${taxon}\\.bracken\\.kraken\\.txt//g' >> combined_bracken_kraken_report_${taxon}_${date}.txt
    done

    """
}


process COMBINE_BRACKEN_FILES_FROM_TAXA {
    publishDir "${params.output_dir}", mode:'copy'
    tag "combine bracken files for $taxon"


    input:
    path(bracken_taxon_files)
    val(taxon)
    val(date)


    output:
    path("combined_bracken_${taxon}_${date}.txt"), emit: bracken_comb_ch


    script:
    """
    BRACKEN_TAXON_FILES=(${bracken_taxon_files})

    for index in \${!BRACKEN_TAXON_FILES[@]}; do
    BRACKEN_TAXON_FILE=\${BRACKEN_TAXON_FILES[\$index]}
    sample_id=\${BRACKEN_TAXON_FILE%.${taxon}_bracken}

    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "Filename\t\$(head -1 \${BRACKEN_TAXON_FILE})" >> combined_bracken_${taxon}_${date}.txt
    fi
    awk -F '\\t' 'FNR>=2 { print FILENAME, \$0 }' \${BRACKEN_TAXON_FILE} |  sed 's/\\.${taxon}\\_bracken//g' >> combined_bracken_${taxon}_${date}.txt
    done


    """
}


process EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG {
    //publishDir "${params.output_dir}", mode:'copy'
    tag "extract from ${sample_id}.bracken.log.file.txt"


    input:
    tuple val(sample_id), path(bracken_log_file)

    output:
    path("*.bracken.unclassified.log.txt"), emit: bracken_uncl_ch


    script:
    """

    echo "Filename\tUnclassifiedReads" > ${sample_id}.bracken.unclassified.log.txt

    FILENAME=\$(cat $bracken_log_file | grep 'BRACKEN OUTPUT PRODUCED' | cut -d':' -f2 | sed 's/\\.S\\_bracken//g')
    UNCLASSIFIEDREADS=\$(cat $bracken_log_file | grep 'Unclassified' | cut -d':' -f2)

    paste <(printf %s \${FILENAME}) <(printf %s \${UNCLASSIFIEDREADS}) >> ${sample_id}.bracken.unclassified.log.txt

    """
}


process COMBINE_UNCLASSIFIED_LOGS {
    publishDir "${params.output_dir}", mode:'copy'


    input:
    path(bracken_unclassified_files)
    val(date)


    output:
    path("Number_of_Unclassified_Reads_bracken_${date}.txt"), emit: bracken_comb_unclassified_ch


    script:
    """
    BRACKEN_UNCLASSIFIED_FILES=(${bracken_unclassified_files})

    for index in \${!BRACKEN_UNCLASSIFIED_FILES[@]}; do
    BRACKEN_UNCLASSIFIED_FILE=\${BRACKEN_UNCLASSIFIED_FILES[\$index]}

    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "\$(head -1 \${BRACKEN_UNCLASSIFIED_FILE})" >> Number_of_Unclassified_Reads_bracken_${date}.txt
    fi
    awk 'FNR>=2 { print }' \${BRACKEN_UNCLASSIFIED_FILE}  >> Number_of_Unclassified_Reads_bracken_${date}.txt
    done

    """
}
