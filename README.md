# Mutation identification in whole genome sequencing data of <i>Arabidopsis thaliana</i>
Bioinformatics pipeline for idintifying mutations in whole genome sequences of <i>Arabidopsis thaliana</i>   
  
  

## Requirement
###
* BEDOPS: the fast, highly scalable and easily-parallelizable genome analysis toolkit (https://bedops.readthedocs.io)
* bedtools: a powerful toolset for genome arithmetic (https://bedtools.readthedocs.io)
* bcftools: Tools for manipulating VCF and BCF files (http://samtools.github.io/bcftools/bcftools.html)
* BioAlcidaeJdk: java-based version of awk for bioinformatics (http://lindenb.github.io/jvarkit/BioAlcidaeJdk.html)
* BWA: Burrow-Wheeler Aligner (http://bio-bwa.sourceforge.net) 
* GATK: Genome Analysis Toolkit (https://gatk.broadinstitute.org)
* samtools: Tools for manipulating NGS data (https://github.com/samtools/samtools)
* vcftools: A set of tools for working with VCF files (https://github.com/vcftools/vcftools)  
* Perl: (https://www.perl.org)  

The environment under CentOS 7.5 is tested. The versions of the tools used are documented in a series of shell scripts.

## <i>Arabidopsis thaliana</i> genomic resources
TAIR10 genomic sequences were downloaded from the TAIR FTP site. The whole-genome sequencing data analyzed were deposited into the DNA Data Bank of Japan Sequence Read Archive with the accession numbers DRA009784.


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

## Filtering parameters
The called raw variants by gatk VariantFiltration in Pipe.09   

###
* Common for filtering out: Depth ＜ 10x, Depth > 200x, GenotypeQuality < 20 
* SNPs for flitering out: QualByDepth < 2.0, FisherStrand > 60.0, RMSMappingQuality < 40.0, MQRankSum < -12.5, ReadPosRankUsm < -8.0, StrandOddsRatio > 4.0, and ExcessHet > 13.0  
* INDELs for flitering out: QualByDepth < 2.0, FisherStrand > 200.0, RMSMappingQuality < 20.0, StrandOddsRatio > 10.0, and ExcessHet > 13.0    
  

The mutation identification in Pipe.10

* mendelian-violation-qual-threshold: 30  
This setting in gatk SelectVariants wiil select only variants that correspond to a mendelian violation as determined on the basis of family structure as <i>P</i> < 0.01.
* Candidate muatation sites having allele frequencies (e.g., proportions of mutant reads at a site) of more than 25% && leth than 80%.
* Final-identified mutation sites having allele frequencies (proportions of mutant reads at a site) of more than 25% && leth than 80% with GQ = 99.


## Note
This project is currently under development. Thank you!

