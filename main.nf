#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()

// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)
final_params = check_params(merged_params)

// starting pipeline
pipeline_start_message(version, final_params)

// include processes
include { BRACKEN; EXTRACT_TAXON_SPECIFIC_INFO; COMBINE_BRACKEN_KRAKEN_REPORT_FROM_TAXA; COMBINE_BRACKEN_FILES_FROM_TAXA;  EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG; COMBINE_UNCLASSIFIED_LOGS } from './modules/processes.nf' addParams(final_params)


workflow  {
          kraken_ch = channel
                          .fromPath( final_params.kraken_report )
                          .map { file -> tuple(file.simpleName, file) }
			  .ifEmpty { error "Cannot find any reads matching: ${final_params.reads}" }


          database_ch = channel.fromPath( final_params.database )


          BRACKEN(kraken_ch, database_ch.collect(), final_params.taxon)

          //for bracken kraken reports
          EXTRACT_TAXON_SPECIFIC_INFO(BRACKEN.out.bracken_kraken_report_ch, final_params.taxon)

          collected_bracken_kraken_taxon_ch = EXTRACT_TAXON_SPECIFIC_INFO.out.taxon_bracken_kraken_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

          COMBINE_BRACKEN_KRAKEN_REPORT_FROM_TAXA(collected_bracken_kraken_taxon_ch, final_params.taxon, final_params.sequencing_date)

          //for bracken reports
          collected_bracken_taxon_ch = BRACKEN.out.bracken_out_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

          COMBINE_BRACKEN_FILES_FROM_TAXA(collected_bracken_taxon_ch, final_params.taxon, final_params.sequencing_date)

          //for extracting unclassified read info from bracken log

          EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG(BRACKEN.out.bracken_log_ch)

          collected_unclassified_log_ch = EXTRACT_UNCLASSIFIED_INFO_FROM_BRACKEN_LOG.out.bracken_uncl_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )
          COMBINE_UNCLASSIFIED_LOGS(collected_unclassified_log_ch, final_params.sequencing_date)
}
