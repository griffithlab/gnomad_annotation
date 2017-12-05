import vcf
import sys

class GnomadVcfParser:
    def __init__(self, path_g_vcf):
        self.path_g_vcf = path_g_vcf

    def parse_vcf(self):
        # create a reader object
        line_count=0
        vcf_reader = vcf.Reader(filename=self.path_g_vcf, compressed=True)
        prog = '1'
        print("\nReading in gnomAD chromosome ", prog)
        # reads in the vcf line by line
        for line in vcf_reader:
            line_count += 1
            if str(line.CHROM) != prog:
                print("Reading in gnomAD chromosome ", str(line.CHROM))
                prog = str(line.CHROM)
            # if line_count % 1000000 == 0:
            #     print("Processing line {}".format(line_count))
            #     sys.stdout.flush()
            # print(line)
            # Simultaneously iterates over the three lists (alt allele, alt allele freq, alt allele count)
            # For each iteration i (allele), gets the values at position i (all values for that allele)
            for alt, af, ac in zip(line.ALT, line.INFO['AF'], line.INFO['AC']):
                # Uses class method to left align (normalize) the position, ref, alt alleles for single allele representation
                new_pos, new_ref, new_alt = GnomadVcfParser.get_minimal_representation(line.POS,line.REF,str(alt))
                # Uses the chr, start, ref, alt as a hash key to provide INFO field information
                key = "_".join([str(line.CHROM),str(new_pos),new_ref,new_alt])
                error_rate = 0
                if af is not None:
                    if new_ref == new_alt:
                        error_rate += 1
                        if error_rate % 1000 == 0:
                            print('Error rate = {}'.format(error_rate))
                        print(key)

    # Removes extra bases added to allow for collapsing of nearby indels into multi-allele VCF representation
    # Provides new start, ref, alt
    # Adapted from http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
    @classmethod
    def get_minimal_representation(cls, pos, ref, alt):
        # If it's a simple SNV, don't remap anything
        if len(ref) == 1 and len(alt) == 1:
            return pos, ref, alt
        # Process indels
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
            # convert to mgi format for insertions/deletions
            if (len(alt) > len(ref)):
                alt = alt[1:]
                ref = "-"
            elif (len(alt) < len(ref)):
                alt = "-"
                ref = ref[1:]
                # Removing first base of the reference allele so increment position
                pos += 1
            return pos, ref, alt
