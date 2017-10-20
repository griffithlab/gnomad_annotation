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
        vcf_reader = vcf.Reader(open(self.path_g_vcf, "rb"))
        g_hash = {}
        for line in vcf_reader:
            for alt, af, ac in zip(line.ALT, line.INFO['AF'], line.INFO['AC']):
                new_pos, new_ref, new_alt = GnomadVcfParser.get_minimal_representation(line.POS,line.REF,str(alt))
                key = "_".join([str(line.CHROM),str(new_pos),new_ref,new_alt])
                g_hash[key] = {}
                g_hash[key]['af'] = af
                g_hash[key]['ac'] = ac
                g_hash[key]['an'] = line.INFO['AN']
        return(g_hash)

    @classmethod
    def get_minimal_representation(cls, pos, ref, alt):
        # If it's a simple SNV, don't remap anything
        if len(ref) == 1 and len(alt) == 1:
            return pos, ref, alt
        else:
            # strip off identical suffixes
            while (alt[-1] == ref[-1] and min(len(alt), len(ref)) > 1):
                alt = alt[:-1]
                ref = ref[:-1]
            # strip off identical prefixes and increment position
            while (alt[0] == ref[0] and min(len(alt), len(ref)) > 1):
                alt = alt[1:]
                ref = ref[1:]
                pos += 1
            # print 'returning: ', pos, ref, alt
            return pos, ref, alt