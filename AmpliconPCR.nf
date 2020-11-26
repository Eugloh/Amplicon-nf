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
params.ext = "fastq.gz"
params.cpus = 2
params.mem = 10
params.prefix      = "pool"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.dbnt = "nt"
params.identity = 98
params.rawReads="/home/lohmanne/AmpliconPCR-nf/FASTQ"

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
    log.info ""
    log.info "Optional arguments:"
    log.info "--prefix        STRING                 prefix (default=pool)"
    log.info '--ext                  STRING               Extension of files (default : fastq.gz)'
    log.info '--cpus                  INTEGER              Number of cpu used by fastqc (default: 2).'
    log.info '--mem                  INTEGER              Size of memory used for mapping (in GB) (default: 10).'
    log.info '--suffix1        STRING                 Suffix of fastq files 1 (default : _1)'
    log.info '--suffix2        STRING                 Suffix of fastq files 2 (default : _2)'
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

raw_reads = params.rawReads

read_pair = Channel
  .fromFilePairs( "${raw_reads}/*_R[1,2]*.fastq", type: 'file')
  .ifEmpty { error "Cannot find any fastq files matching: ${params.input_folder}" }
  .into  { read_files_fastqc; read_files_trimG }



//file data from Channel.fromPath(${params.fastq_dir}+'*').collect()


log.info'##########################################\n##\tFastQC control of the raw reads\t##\n##########################################'

process fastqc_fastq {
  publishDir 'out/fastqc/' , mode: 'copy', overwrite: false

  input:
  set pair_id, file(reads) from read_files_fastqc

  output:
    file "*.{zip,html}" into fastqc_report
    stdout into name

  script:
  """
  fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
  ${reads[0]} ${reads[1]}
  """
}



process multiqc {

  publishDir 'out/multiqc/', mode: 'copy', overwrite: false
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


process trim_galore {
    publishDir 'out/trimgalore/', mode: 'copy', overwrite: false,
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (filename.indexOf("fq.gz") > 0) "out/$filename"
            else null}

    input:
    set val(name), file(reads) from read_files_trimG

    output:
    set val(name), file("*_val_*.fq.gz") into fastq_files_trim    
    file("*_unpaired_*.fq.gz") into unpaired
    file("*trimming_report.txt") into trimgalore_results
    file("*_fastqc.{zip,html}") into trimgalore_fastqc_reports


    script:
        """
        echo ${name};
        echo ${reads};

        trim_galore --paired --fastqc --length 30 --retain_unpaired --gzip ${reads[0]} ${reads[1]} 
  
        """
    }

// fastq_files_trim.view { "value: $it" }



// pool30_S30_L001_R1_001_val_1.fq.gz
// pool30_S30_L001_R1_001_val_2.fq.gz

process multiqc_trim {

  publishDir 'out/trimgalore/multiqc/', mode: 'copy' , overwrite: false
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




log.info'##########################################\n##\tClustering step : VSEARCH\t##\n##########################################'

/*usearch -threads 1 -fastq_mergepairs ${sample}_forward_renamed.fastq \
        -reverse ${sample}_reverse_renamed.fastq \
        -fastqout ${sample}_merged.fastq \
        -fastq_maxdiffs ${params.fastqMaxdiffs}
// */


process Merging {

    tag { "${params.prefix}.${sample}" }
    label 'vsearch'
    memory { 4.GB * task.attempt }
    publishDir 'out/vsearch/log/' , mode: 'copy', overwrite: false

    input:  
    set val(name), file(reads) from fastq_files_trim


    output:
    set val(name), file("*_merged.fasta") into merged_read_pair
    file("*_UnMFwd.fasta") into unmergedForward
    file("*_UnMRev.fasta") into unmergedReverse

    """

    vsearch --quiet --fastq_mergepairs ${reads[0]}  \
    --reverse  ${reads[1]} \
    --threads ${params.cpus} \
    --fastq_allowmergestagger --label_suffix ${params.prefix} \
    --fastaout_notmerged_fwd ${name}_UnMFwd.fasta \
    --fastaout_notmerged_rev ${name}_UnMRev.fasta \
    --fastaout ${name}_merged.fasta

    
    """
}

log.info'##########################################\n##\t FUSION DES UNMERGED \t##\n##########################################'

process Fusion {

    tag { "${params.prefix}.${sample}" }
    label 'fusion'
    memory { 4.GB * task.attempt }
    publishDir 'out/vsearch/log/' , mode: 'copy', overwrite: false

    input:  
    set val(name), file(reads) from merged_read_pair
    file(readsF) from unmergedForward
    file(readsR) from unmergedReverse
  
    output:
    set val(name), file("*_paired.fasta") into fusionfasta


    """
cat ${readsF} ${reads}  >  ${name}_paired.fasta ;

cat ${readsR} >> ${name}_paired.fasta 
    
    """
}

log.info'##########################################\n##Dereplication step : VSEARCH\t##\n##########################################'

process  vparseDerepWorkAround {
    tag { "${params.prefix}.${sample}" }
    label 'derep'
    memory { 4.GB * task.attempt }
    publishDir "out/vsearch/log/", mode: 'copy', overwrite: false,
      saveAs: {filename ->
          if (filename.indexOf("lin_der") > 0) "$filename"
          else if (filename.indexOf("log") > 0) "log_files/$filename"
          else null}


    input:
	    set val(name), file(read) from fusionfasta

    output:
      set val(name), file('*lin_der.fas') into derep_fasta

    """
   vsearch --quiet --derep_fulllength ${read} \
    --sizeout --threads ${params.cpus} --relabel_sha1 --fasta_width 0 \
    --strand both \
    --minuniquesize 1 --output ${name}_lin_der.fas \
    --log ${name}_lin_der.log


    """
}

log.info'##########################################\n##Chimerics removal step : VSEARCH\t##\n##########################################'

process  ChimericRemov {
    tag { "${params.prefix}.${sample}" }
    label 'derep'
    memory { 4.GB * task.attempt }
    publishDir "out/vsearch/log/", mode: 'copy', overwrite: false,
      saveAs: {filename ->
          if (filename.indexOf("no_chim") > 0) "$filename"
          else if (filename.indexOf("log") > 0) "log_files/$filename"
          else null}

    input:
	   set val(name), file(read) from derep_fasta

    output:
    set val(name), file("*no_chim.fasta") into no_chimfasta

    """
vsearch --quiet --uchime_denovo ${read}  \
    --sizein --threads ${params.cpus} --relabel ${params.prefix} \
    --sizeout --xsize --nonchimeras ${name}_no_chim.fasta \
    --log ${name}_chimeria.log 

    """
}


log.info'##########################################\n##Clustering step : VSEARCH\t##\n##########################################'


process  Clustering {
    tag { "${params.prefix}.${sample}" }
    label 'derep'
    memory { 4.GB * task.attempt }
    publishDir "out/vsearch/log/", mode: 'copy', overwrite: false,
      saveAs: {filename ->
          if (filename.indexOf("clustered") > 0) "$filename"
          else if (filename.indexOf("log") > 0) "log_files/$filename"
          else null}

    input:
	   set val(name), file(read) from no_chimfasta

    output:
    set val(name), file("*_clustered.fasta") into clusteredfasta



"""
vsearch --quiet --cluster_size ${read} \
    --id 0.97 --threads  ${params.cpus} \
    --sizein --clusterout_id --clusterout_sort \
    --sizeout --xsize --relabel ${params.prefix} \
    --centroids ${name}_clustered.fasta \
    --log ${name}_clustered.log
        """
}

log.info'##########################################\n##Blast\t##\n##########################################'



process  Blast {
    tag { "${params.prefix}.${sample}" }
    label 'derep'
    memory { 30.GB * task.attempt }
    publishDir "out/blast_result/", mode: 'copy', overwrite: false

    input:
	  set val(name), file(read) from clusteredfasta

    output:
	   set val(name), file("*blast") into blast


"""
 blastn -task megablast -db nt -query ${read}  -out ${name}.blast \
     -evalue 1e-05 -max_target_seqs 1 -num_threads ${params.cpus}  \
     -outfmt "6 qseqid sseqid evalue bitscore length pident frames staxids sskingdoms sscinames scomnames sblastnames stitle qseq qstart qend" ;

     sed -i '1iQueryID\tSubjectID\tevalue\tbitscore\tlength query\tperc id\tframes\ttaxid\tkingdom\tscientifique name\tcommon name\tblast name\ttitle\tseq query\tstartq\tstopq'  ${name}.blast

"""
}