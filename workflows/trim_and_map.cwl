cwlVersion: v1.0
class: Workflow

doc: |
  This workflow trims adapters from single or paired-end sequencing reads and alings them to a reference genome.
  For trimming, Trim Galore v0.4.4 is used. Alignmnent is performed using Bowtie2 v2.2.6-2.

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  fastq1:
    doc: |
      FASTQ file containing raw reads. 
      For paired end data, please only provide the first reads of read pair.
    type: File
  fastq2: 
    doc: |
      FASTQ file containing the second raw reads of a pair. 
      Only relevant for paired end data.
    type: File?
  genome:
    doc: |
      Path to reference genome in fasta format. 
      Bowtie2 index files (".1.bt2", ".2.bt2", ...) as well as a samtools index (".fai") 
      has to be located in the same directory.\n
      All of these files can be downloaded for the most common genome builds at  
      https://support.illumina.com/sequencing/sequencing_software/igenome.html. 
      Alternatively, you can use "bowtie2-build" or "samtools index" to create them yourself.
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
  adapters: 
    doc: |
      Adapter to be trimmed from raw reads. Cab be one of the following: \n
      - "nextera" for the Nextera adapter (CTGTCTCTTATA)\n
      - "illumina" for the Illumina universal adapter (AGATCGGAAGAGC)\n
      - "small_rna" for the Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)\n
      - "auto" to automatically detect the write setting
    type:
      type: enum
      symbols:
        - nextera
        - illumina
        - small_rna
        - auto
  is_paired_end:
    doc: |
      Specify whether the data is paired or single ended.
    type: boolean
  max_mapping_insert_length:
    doc: |
      Maximum insert length between two reads of a pair. In case of ATACseq, 
      very long insert sizes are possible. So it is recommended to use at least 
      a value of 1500. However, please note that alignment will take significantly 
      longer for higher insert sizes. The default is 2500.
    type: long?
    default: 500
        
steps:
  qc_raw:
    doc: fastqc - quality control for trimmed fastq
    run: "../tools/fastqc.cwl"
    in:
      fastq1:
        source: fastq1
      fastq2:
        source: fastq2
    out:
      - fastqc_zip
      - fastqc_html

  adaptor_trimming_and_qc_trimmed:
    doc: trim galore - adapter trimming using trim_galore
    run: "../tools/trim_galore.cwl"
    in:
      fastq1:
        source: fastq1
      fastq2:
        source: fastq2
      adapter1:
        source: adapters
      adapter2:
        source: adapters   
    out:
      - fastq1_trimmed
      - fastq2_trimmed
      - fastq1_trimmed_unpaired
      - fastq2_trimmed_unpaired
      - trim_galore_log
      - trimmed_fastqc_html
      - trimmed_fastqc_zip

  mapping:
    doc: bowite2 - mapper, produces sam file
    run: "../tools/bowtie2.cwl"
    in:
      fastq1:
        source: adaptor_trimming_and_qc_trimmed/fastq1_trimmed
      fastq2:
        source: adaptor_trimming_and_qc_trimmed/fastq2_trimmed
      genome_index:
        source: genome
      is_paired_end:
        source: is_paired_end
      max_mapping_insert_length:
        source: max_mapping_insert_length
    out:
      - sam

  sam2bam:
    doc: samtools view - convert sam to bam
    run: "../tools/samtools_view_sam2bam.cwl"
    in:
      sam:
        source: mapping/sam
    out:
      - bam_unsorted
      
  sort_bam: 
    doc: samtools sort - sorts unsorted bam file by coordinates.
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: sam2bam/bam_unsorted
    out:
      - bam_sorted
      
outputs:
  raw_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: qc_raw/fastqc_zip
  raw_fastqc_html:
    type:
      type: array
      items: File
    outputSource: qc_raw/fastqc_html
  fastq1_trimmed:
    type: File
    outputSource: adaptor_trimming_and_qc_trimmed/fastq1_trimmed
  fastq2_trimmed:
    type: File?
    outputSource: adaptor_trimming_and_qc_trimmed/fastq2_trimmed
  trim_galore_log:
    type:
      type: array
      items: File
    outputSource: adaptor_trimming_and_qc_trimmed/trim_galore_log
  trimmed_fastqc_html:
    type:
      type: array
      items: File
    outputSource: adaptor_trimming_and_qc_trimmed/trimmed_fastqc_html
  trimmed_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: adaptor_trimming_and_qc_trimmed/trimmed_fastqc_zip
  bam:
    type: File
    outputSource: sort_bam/bam_sorted
    