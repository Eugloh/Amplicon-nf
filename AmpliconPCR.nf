#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null
params.input_folder = "FASTQ/"
params.output_folder = "."
params.ext = "fastq.gz"
params.multiqc_config = 'NO_FILE'
params.cpu = 2
params.mem = 10
params.mem_QC       = 2
params.prefix      = "pool"
params.suffix1      = "_1"
params.suffix2      = "_2"

params.fastq_dir = null
params.fastq = null
params.working_dir = null
params.dbnt = "nt"
params.info = null
params.identity = 98



log.info ""
log.info "--------------------------------------------------------"
log.info " AmpliconPCR - November 2020 - Amplicon analysis         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run Eugloh/AmpliconPCR [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_folder         FOLDER               Folder containing fastq files"
    log.info "--output_folder        FOLDER                 Output directory for html and zip files (default=fastqc_ouptut)"
    log.info "--prefix        STRING                 prefix (default=pool)"
    log.info ""
    log.info "Optional arguments:"
    log.info '--ext                  STRING               Extension of files (default : fastq.gz)'
    log.info '--multiqc_config       PATH                 config yaml file for multiqc (default : none)'
    log.info '--cpu                  INTEGER              Number of cpu used by fastqc (default: 2).'
    log.info '--mem                  INTEGER              Size of memory used for mapping (in GB) (default: 10).'
    log.info '--mem_QC         INTEGER                Size of memory used for QC and cutadapt (in GB) (default: 32).'
    log.info '--suffix1        STRING                 Suffix of fastq files 1 (default : _1)'
    log.info '--suffix2        STRING                 Suffix of fastq files 2 (default : _2)'
    log.info "--fastq_dir         PATH"
    log.info "--fastq             STRING"
    log.info "--working_dir       PATH "
    log.info "--dbnt              STRING  (default:nt) "
    log.info "--info              PATH"
    log.info "--identity          FLOAT"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

assert (params.input_folder != null) : "please provide the --input_folder option"

infiles = Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .into  { read_files_fastqc; read_files_trimG }

//file data from Channel.fromPath(${params.fastq_dir}+'*').collect()

//multiqc config file
ch_config_for_multiqc = file(params.multiqc_config)


log.info'##########################################\n##\tFastQC control of the raw reads\t##\n##########################################'

process fastqc_fastq {
  //tag "$pair_id"
  publishDir "results/fastq/fastqc/", mode: 'copy', overwrite: false

  input:
  set pair_id, file(reads) from read_files_fastqc

  output:
    file "*.{zip,html}" into fastqc_report

  script:
  """
  fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
  ${reads[0]} ${reads[1]}
  """
}



process multiqc {
  //tag "$report[0].baseName"

  publishDir "results/fastq/multiqc/", mode: 'copy', overwrite: false
  cpus = 1

  input:
    file report from fastqc_report.collect()

  output:
    file "*multiqc_*" into multiqc_report

  script:
  """
  multiqc -f .
  """
}

log.info'##########################################\n##\tRemove adapter sequence\t\t##\n##########################################'
// Trim Galore & FastQC


process trim_galore {
    //tag "$name"
    publishDir "results/fastq/trim_galore", mode: 'copy', overwrite: false,
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else null}

    input:
    set val(name), file(reads) from read_files_trimG

    output:
    file("*_val_*.fq.gz") into fastqc_postpairs    
    file("*_unpaired_*.fq.gz") into unpaired
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports


    script:
        """
        trim_galore --paired --fastqc --length 30 --retain_unpaired --gzip $reads 
  
        """
    }


process multiqc_trim {
  //tag "$report[0].baseName"

  publishDir "results/fastq/trim_galore/multiqc/", mode: 'copy' , overwrite: false
  cpus = 1

  input:
    file report from trimgalore_fastqc_reports.collect()

  output:
    file "*multiqc_*" into multiqc_trimgalore_report

  script:
  """
  multiqc -f .
  """
}
