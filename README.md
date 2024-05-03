# Mutation identification in whole genome sequencing data of <i>Arabidopsis thaliana</i>
Bioinformatics pipeline for identifying mutations in whole genome sequences of <i>Arabidopsis thaliana</i>   
  

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
  <img src="https://github.com/akihirao/At_Reseq/blob/main/images/AtReseq.workflow.pindel.jpeg"/>
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
A list of all mutations identified is provided as [M2.mutations.full.list.csv](https://github.com/akihirao/At_Reseq/blob/main/M2.mutations.full.list.csv)

## Filtering parameters
The called raw variants by gatk VariantFiltration in Pipe 10   

###
* Common for filtering out: Depth ＜ 10x, Depth > 200x, GenotypeQuality < 20 
* SNPs for flitering out: QualByDepth < 2.0, FisherStrand > 60.0, RMSMappingQuality < 40.0, MQRankSum < -12.5, ReadPosRankUsm < -8.0, StrandOddsRatio > 4.0, and ExcessHet > 13.0  
* INDELs for flitering out: QualByDepth < 2.0, FisherStrand > 200.0, RMSMappingQuality < 20.0, StrandOddsRatio > 10.0, and ExcessHet > 13.0    
  

The mutation identification in Pipe 11

* --mendelian-violation-qual-threshold: 30  
This setting in gatk SelectVariants will select only variants that correspond to a mendelian violation as determined on the basis of family structure as <i>P</i> < 0.01.
* --max-nocall-fraction: 0  
This setting in gatk SelectVariants will select only variants having genotyping rate among samples of 100%.
* The candidate mutation sites having allele frequencies (AF; proportions of mutant reads at a site) of 25% or less were excluded.
* The candidate mutation sites with the other neighbor sites within 150bp range in a single sample were disregarded because a high proportion of false positives of neighboring mutation sites are known to be cause by mismapping (Keightley <i>et al</i>. 2014)


## R codes for Statistical analyses
* [Statistical modeling of the effect of radiation on the number of each type of mutation & plotting fig.2](https://github.com/akihirao/At_Reseq/blob/main/R_work/Plot.fig2.AT.mutation.vs.dose.md)
* [Statistical test for the ratio of heterozygous/homozygous mutations & plotting fig.3](https://github.com/akihirao/At_Reseq/blob/main/R_work/Plot.fig3.AT.homo2hetero.md)
* [Statistical modeling of the effect of radiation on size of mutation & plotting fig.S2](https://github.com/akihirao/At_Reseq/blob/main/R_work/Plot.figS2.INDEL.size.dose.md)
* [Statistical modeling of the effect of radiation on gene function $ plotting fig.S3](https://github.com/akihirao/At_Reseq/blob/main/R_work/Plot.figS3.AT.Snpeff.md)

## Reference
Hirao AS, Watanabe Y, Hasegawa Y, Takagi T, Ueno S, Kaneko S (2022) Mutational effects of chronic gamma radiation throughout the life cycle of <i>Arabidopsis thaliana</i>: insight into radiosensitivity in the reproductive stage. Science of the Total Environment. https://doi.org/10.1016/j.scitotenv.2022.156224

