# Mutation identification in whole genome sequencing data of <i>Arabidopsis thaliana</i>
Bioinformatics pipeline for idintifying mutations in whole genome sequences of <i>Arabidopsis thaliana</i>   
  

## Requirement
###
* BEDOPS: the fast, highly scalable and easily-parallelizable genome analysis toolkit https://bedops.readthedocs.io
* bedtools: a powerful toolset for genome arithmetic https://bedtools.readthedocs.io
* bcftools: Tools for manipulating VCF and BCF files http://samtools.github.io/bcftools/bcftools.html
* BioAlcidaeJdk: java-based version of awk for bioinformatics http://lindenb.github.io/jvarkit/BioAlcidaeJdk.html
* BWA: Burrow-Wheeler Aligner http://bio-bwa.sourceforge.net 
* GATK: Genome Analysis Toolkit https://gatk.broadinstitute.org
* samtools: Tools for manipulating NGS data https://github.com/samtools/samtools
* snpEff: Genomic variant annotations and functional effect prediction toolbox http://snpeff.sourceforge.net
* Perl: https://www.perl.org
* Pindel: A indel calling tool based on the split-read approach http://gmt.genome.wustl.edu/packages/pindel
* vcftools: A set of tools for working with VCF files https://github.com/vcftools/vcftools 

The environment under CentOS 7.5 is tested. The versions of the tools used are documented in a series of shell scripts.

## <i>Arabidopsis thaliana</i> genomic resources
TAIR10 genomic sequences were downloaded from the TAIR FTP site. The whole-genome sequencing data analyzed were deposited into the DNA Data Bank of Japan Sequence Read Archive with the accession number DRA009784.


## Flowchart
Variant call procedure was adapted from the germline short variant discovery workflow in GATK.

<p align="left">
  <img src="https://github.com/akihirao/At_Reseq/blob/master/images/AtReseq.workflow.pindel.jpeg"/>
</p>


## Usage
Run a series of the shell scripts in the order listed after changing paths according to your environemt:
 
```bash
Pipe.01.Map.sh
Pipe.02.MarkduplicatesSpark.sh
Pipe.03.BaseRecalibrator.sh
...
Pipe.14.Annotation.sh
```
A full list of all identified mutations is output to the current folder as "M2.mutations.full.list.csv"

## Filtering parameters
The called raw variants by gatk VariantFiltration in Pipe10   

###
* Common for filtering out: Depth ＜ 10x, Depth > 200x, GenotypeQuality < 20 
* SNPs for flitering out: QualByDepth < 2.0, FisherStrand > 60.0, RMSMappingQuality < 40.0, MQRankSum < -12.5, ReadPosRankUsm < -8.0, StrandOddsRatio > 4.0, and ExcessHet > 13.0  
* INDELs for flitering out: QualByDepth < 2.0, FisherStrand > 200.0, RMSMappingQuality < 20.0, StrandOddsRatio > 10.0, and ExcessHet > 13.0    
  

The mutation identification in Pipe11

* --mendelian-violation-qual-threshold: 30  
This setting in gatk SelectVariants will select only variants that correspond to a mendelian violation as determined on the basis of family structure as <i>P</i> < 0.01.
* --max-nocall-fraction: 0  
This setting in gatk SelectVariants will select only variants having genotyping rate among samples of 100%.
* The candidate mutation sites having allele frequencies (AF; proportions of mutant reads at a site) of 25% or less were excluded.
* The candidate mutation sites with the other neighbor sites within 150bp range in a single sample were disregarded because a high proportion of false positives of neighboring mutation sites are known to be cause by mismapping (Keightley <i>et al</i>. 2014)


## Statistical analyses
* [Statistical modeling of the number of each type of mutation and radiation dose](https://github.com/akihirao/AT_Reseq/R_work/blob/master/Plot.fig2.AT.mutation.vs.dose.md)

## Note
This project is currently under development. Thank you!

