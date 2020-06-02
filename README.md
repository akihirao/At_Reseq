# Mutation identification in whole genome sequencing data of <i>Arabidopsis thaliana</i>
Bioinformatics pipeline for idintifying mutations in whole genome sequences of <i>Arabidopsis thaliana</i> 
 

## Requirement

* BEDOPS: the fast, highly scalable and easily-parallelizable genome analysis toolkit (https://bedops.readthedocs.io)
* bcftools: Tools for manipulating VCF and BCF files (http://samtools.github.io/bcftools/bcftools.html)
* BioAlcidaeJdk: java-based version of awk for bioinformatics (http://lindenb.github.io/jvarkit/BioAlcidaeJdk.html)
* BWA: Burrow-Wheeler Aligner (http://bio-bwa.sourceforge.net) 
* GATK: Genome Analysis Toolkit (https://gatk.broadinstitute.org)
* samtools: Tools for manipulating NGS data (https://github.com/samtools/samtools)
* vcftools: A set of tools for working with VCF files (https://github.com/vcftools/vcftools)


The environment under CentOS 7.5 is tested. The versions of the tools used are documented in a series of shell scripts.

## <i>Arabidopsis thaliana</i> genomic resources
TAIR10 genomic sequences were downloaded from the TAIR FTP site.


## Flowchart
This workflow is referred to the practice "Germline short variant discovery (SNPs + INDELs)" in GATK.

<p align="left">
  <img src="https://github.com/akihirao/At_Reseq/blob/master/images/AtReseq.workflow.jpeg"/>
</p>


## Usage
Run a series of the shell scripts in the order listed after changing paths according to your environemt:
 
```bash
Pipe.01.Map.sh
Pipe.02.MarkduplicatesSpark.sh
Pipe.03.BaseRecalibrator.sh
...
Pipe.10.MutationIdentification.sh
```



## Note
The whole-genome sequencing data analyzed were deposited into the DNA Data Bank of Japan Sequence Read Archive (https://ddbj.nig.ac.jp/dra) with the accession numbers DRA009784.  
  
This project is currently under development. Thank you!

