## Workflow to estimate microbial taxon abundance using bracken.
### Usage

```

===============================================
 BRACKEN TAPIR Pipeline version 1.0dev
===============================================
 The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" --database "PathToDB" --sequencing_date "GYYMMDD" --taxon "taxon_symbol" --kraken_report "*.kraken.report.txt"

        Mandatory arguments:
         --database                     KRAKEN database directory (full path required, e.g., "/KRAKEN_DB")
         --output_dir                   Output directory to place final combined kraken output (e.g., "/MIGE/01_DATA/11_QC_REPORTS")
         --sequencing_date              Sequencing Date (for TAPIR, must start with G e.g., "G230223")
         --kraken_report                kraken report file generated from kraken analysis

        Optional arguments:
         --taxon                        taxonomic rank symbol (e.g., "D", "P", "C", "O", "F", "G", "G1", "S", "S1", Default: "S")
         --help                         This usage statement.
         --version                      Version statement

```


## Introduction
This pipeline estimates microbial abundance at any taxonomic level derived from kraken report files. This Nextflow pipeline was adapted from the Bracken [github page](https://github.com/jenniferlu717/Bracken), and the NF Core Bracken Module [github page](https://github.com/nf-core/modules/tree/master/modules/nf-core/bracken/bracken).  


## Sample command
An example of a command to run this pipeline is:

```
nextflow run main.nf  --output_dir "test2" --database "/KRAKEN_DB" --sequencing_date "G230202" --taxon "S" --kraken_report "test2/*.kraken.report.txt"
```

## Word of Note
This is an ongoing project at the Microbial Genome Analysis Group, Institute for Infection Prevention and Hospital Epidemiology, Üniversitätsklinikum, Freiburg. The project is funded by BMBF, Germany, and is led by [Dr. Sandra Reuter](https://www.uniklinik-freiburg.de/institute-for-infection-prevention-and-control/microbial-genome-analysis.html).


## Authors and acknowledgment
The TAPIR (Track Acquisition of Pathogens In Real-time) team.
