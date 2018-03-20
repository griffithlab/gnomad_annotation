import vcf
import sys


class GnomadVcfParserExome:
    def __init__(self, path_g_vcf):
        self.path_g_vcf = path_g_vcf

    def parse_vcf(self):
        # create a reader object
        line_count = 0
        vcf_reader = vcf.Reader(filename=self.path_g_vcf, compressed=True)
        prog = '1'
        print("\nReading in gnomAD chromosome ", prog)
        # reads in the vcf line by line
        for line in vcf_reader:
            line_count += 1
            if str(line.CHROM) != prog:
                print("Reading in gnomAD chromosome ", str(line.CHROM))
                prog = str(line.CHROM)
            if line_count % 1000000 == 0:
                print("Processing line {}".format(line_count))
                sys.stdout.flush()
            # print(line)
            # Simultaneously iterates over the three lists (alt allele, alt allele freq, alt allele count)
            # For each iteration i (allele), gets the values at position i (all values for that allele)
            for alt, af, ac, af_afr, ac_afr, af_amr, ac_amr, af_asj, ac_asj, af_eas, ac_eas, af_fin, ac_fin, af_nfe, \
                ac_nfe, af_sas, ac_sas, af_oth, ac_oth in zip(line.ALT, line.INFO['AF'], line.INFO['AC'],
                                                              line.INFO['AF_AFR'], line.INFO['AC_AFR'],
                                                              line.INFO['AF_AMR'], line.INFO['AC_AMR'],
                                                              line.INFO['AF_ASJ'], line.INFO['AC_ASJ'],
                                                              line.INFO['AF_EAS'], line.INFO['AC_EAS'],
                                                              line.INFO['AF_FIN'], line.INFO['AC_FIN'],
                                                              line.INFO['AF_NFE'], line.INFO['AC_NFE'],
                                                              line.INFO['AF_SAS'], line.INFO['AC_SAS'],
                                                              line.INFO['AF_OTH'], line.INFO['AC_OTH']):
                # Uses class method to left align (normalize) the position, ref, alt alleles
                #  for single allele representation
                new_pos, new_ref, new_alt = GnomadVcfParserExome.get_minimal_representation(line.POS, line.REF, str(alt))
                # Uses the chr, start, ref, alt as a hash key to provide INFO field information
                key = "_".join([str(line.CHROM), str(new_pos), new_ref, new_alt])
                if af is not None:
                    yield key, (float(af), int(ac), int(line.INFO['AN']), float(af_afr) if af_afr is not None else 0,
                                int(ac_afr), int(line.INFO['AN_AFR']), float(af_amr) if af_amr is not None else 0,
                                int(ac_amr), int(line.INFO['AN_AMR']), float(af_asj) if af_asj is not None else 0,
                                int(ac_asj), int(line.INFO['AN_ASJ']), float(af_eas) if af_eas is not None else 0,
                                int(ac_eas), int(line.INFO['AN_EAS']), float(af_fin) if af_fin is not None else 0,
                                int(ac_fin), int(line.INFO['AN_FIN']), float(af_nfe) if af_nfe is not None else 0,
                                int(ac_nfe), int(line.INFO['AN_NFE']), float(af_sas) if af_sas is not None else 0,
                                int(ac_sas), int(line.INFO['AN_SAS']), float(af_oth) if af_oth is not None else 0,
                                int(ac_oth), int(line.INFO['AN_OTH']))

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


class GnomadVcfParserGenome:
    def __init__(self, path_g_vcf):
        self.path_g_vcf = path_g_vcf

    def parse_vcf(self):
        # create a reader object
        line_count = 0
        vcf_reader = vcf.Reader(filename=self.path_g_vcf, compressed=True)
        prog = '1'
        print("\nReading in gnomAD chromosome ", prog)
        # reads in the vcf line by line
        for line in vcf_reader:
            line_count += 1
            if str(line.CHROM) != prog:
                print("Reading in gnomAD chromosome ", str(line.CHROM))
                prog = str(line.CHROM)
            if line_count % 1000000 == 0:
                print("Processing line {}".format(line_count))
                sys.stdout.flush()
            # print(line)
            # Simultaneously iterates over the three lists (alt allele, alt allele freq, alt allele count)
            # For each iteration i (allele), gets the values at position i (all values for that allele)
            for alt, af, ac, af_afr, ac_afr, af_amr, ac_amr, af_asj, ac_asj, af_eas, ac_eas, af_fin, ac_fin, af_nfe,\
                ac_nfe, af_oth, ac_oth in zip(line.ALT, line.INFO['AF'], line.INFO['AC'],
                                              line.INFO['AF_AFR'], line.INFO['AC_AFR'], line.INFO['AF_AMR'],
                                              line.INFO['AC_AMR'], line.INFO['AF_ASJ'], line.INFO['AC_ASJ'],
                                              line.INFO['AF_EAS'], line.INFO['AC_EAS'], line.INFO['AF_FIN'],
                                              line.INFO['AC_FIN'], line.INFO['AF_NFE'], line.INFO['AC_NFE'],
                                              line.INFO['AF_OTH'], line.INFO['AC_OTH']):
                # Uses class method to left align (normalize) the position, ref, alt alleles
                # for single allele representation
                new_pos, new_ref, new_alt = GnomadVcfParserGenome.get_minimal_representation(line.POS, line.REF, str(alt))
                # Uses the chr, start, ref, alt as a hash key to provide INFO field information
                key = "_".join([str(line.CHROM), str(new_pos), new_ref, new_alt])
                if af is not None:
                    yield key, (float(af), int(ac), int(line.INFO['AN']), float(af_afr) if af_afr is not None else 0,
                                int(ac_afr), int(line.INFO['AN_AFR']), float(af_amr) if af_amr is not None else 0,
                                int(ac_amr), int(line.INFO['AN_AMR']), float(af_asj) if af_asj is not None else 0,
                                int(ac_asj), int(line.INFO['AN_ASJ']), float(af_eas) if af_eas is not None else 0,
                                int(ac_eas), int(line.INFO['AN_EAS']), float(af_fin) if af_fin is not None else 0,
                                int(ac_fin), int(line.INFO['AN_FIN']), float(af_nfe) if af_nfe is not None else 0,
                                int(ac_nfe), int(line.INFO['AN_NFE']), float(af_oth) if af_oth is not None else 0,
                                int(ac_oth), int(line.INFO['AN_OTH']))

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
