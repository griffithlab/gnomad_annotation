# gnomad_annotation
Python package to add gnomAD allele frequencies to MGI annotated files
- gnomAD exome and genome information are provided separately

### Things to note about gnomAD 
- Does not provide Y chromosome information
- Allele counts, numbers and frequencies provided by this script are filtered high quality genotypes, not raw counts (GQ >= 20, DP >= 10, allele balance > 0.2 for heterozygote genotypes), previously referred to as adjusted allele frequencies for ExAC
- [Read more](https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad/)

### Uses Ensembl VCF files as initial source. 
- Source: ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/ GRCh[version]/variation_genotype/gnomad.[version].vcf.gz
- Provides smaller files than gnomAD directly
- Provides GRCH38 liftover coordinates
- Imported to an MGI allocation using https://github.com/genome/genome/blob/master/lib/perl/Genome/Db/Ensembl/Gnomad.pm

### VCF files are first processed and stored in trie structures for fast annotation
- vcf_to_pickle.py
- Initial creation takes ~1 day for genomes file

# Code is currently actively developed, use at your own risk
- There is no lower limit on the number of alleles contributing to the allele frequency provided
- Alleles not found and alleles with no information (no high quality genotype data) are annotated with NA
