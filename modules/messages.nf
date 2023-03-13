def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" --database "PathToDB" --sequencing_date "GYYMMDD" --taxon "taxon_symbol" --kraken_report "*.kraken.report.txt"

        Mandatory arguments:
         --reads                        Query fastqz file of sequences you wish to supply as input (e.g., "/MIGE/01_DATA/01_FASTQ/T055-8-*.fastq.gz")
         --database                	KRAKEN database directory (full path required, e.g., "/KRAKEN_DB")
         --output_dir                   Output directory to place final combined kraken output (e.g., "/MIGE/01_DATA/03_ASSEMBLY")
         --sequencing_date              Sequencing Date (for TAPIR, must start with G e.g., "G230223")
	 --kraken_report                kraken report file generated from kraken analysis

        Optional arguments:
         --taxon                        taxonomic rank symbol (e.g., "D", "P", "C", "O", "F", "G", "S", "S1", Default: "S")
         --help                         This usage statement.
         --version                      Version statement
        """
}

version = '1.0dev'

def Version() {
      println(
            """
            ===============================================
             BRACKEN TAPIR Pipeline version ${version}
            ===============================================
            """.stripIndent()
        )

}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Show version
if (params.version) {
    Version()
    exit 0
}


