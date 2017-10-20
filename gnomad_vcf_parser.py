## Parses the imported gnomAD VCF file
## demultiplex alleles
## extract the allele frequencies
## left align the reference
## Export a hash where key is chr,start,ref,var and value is allele freq (allele counts?)

import vcf

class GnomadVcfParser:
    def __init__(self, path_g_vcf):
        self.path_g_vcf = path_g_vcf
    def parse_vcf(self):
        vcf_reader = vcf.Reader(open(self.path_g_vcf))
        for line in vcf_reader:
            