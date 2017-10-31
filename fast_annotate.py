import marisa_trie
import gzip
import csv


def dict_reader_zipped(dict_filename):
    cur_line = 0
    with gzip.open(dict_filename, 'r') as zipfile:
        for line in zipfile:
            if len(line) == 0:
                break
            fields = line.strip().split('\t')
            cur_line += 1
            if cur_line % 1000000 == 0:
                print('Processing line {}'.format(cur_line))
            yield (fields[0]).decode('utf-8'), (float(fields[1]), int(fields[2]), int(fields[3]))
        print('Finished processing at line {}'.format(cur_line))


def dict_reader_unzipped(dict_filename):
    cur_line = 0
    with open(dict_filename, 'r') as dictfile:
        for line in dictfile:
            if len(line) == 0:
                break
            fields = line.strip().split('\t')
            cur_line += 1
            if cur_line % 1000000 == 0:
                print('Processing line {}'.format(cur_line))
            yield (fields[0]), (float(fields[1]), int(fields[2]), int(fields[3]))
        print('Finished processing at line {}'.format(cur_line))


def annotate(mutation_filename, gnomad_annotations):
    annotated_filename = mutation_filename + ".out"
    with open(mutation_filename, "r") as mgi_tsv, open(annotated_filename, "w") as outfile:
        print('Beginning comparison')
        mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
        header = mgi_tsv_reader.fieldnames
        header_new = header + ["gnomAD_AC", "gnomAD_AN", "gnomAD_AF"]
        mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t")
        mgi_tsv_writer.writeheader()
        line_count = 0
        for line in mgi_tsv_reader:
            line_count += 1
            if line_count % 50000 == 0:
                print('Processing line {} from inputfile'.format(line_count))
            mgi_key = "_".join([str(line["chromosome_name"]), str(line["start"]), line["reference"], line["variant"]])
            if mgi_key not in gnomad_annotations:
                mgi_key = mgi_key.replace("-", "0")
            if mgi_key in gnomad_annotations:
                gnomad_record = gnomad_annotations[mgi_key]
                print(gnomad_record)
                new_line = line.copy()
                new_line["gnomAD_AC"] = gnomad_record[0][1]
                new_line["gnomAD_AN"] = gnomad_record[0][2]
                new_line["gnomAD_AF"] = gnomad_record[0][0]
                mgi_tsv_writer.writerow(new_line)
            else:
                new_line = line.copy()
                new_line["gnomAD_AC"] = "NA"
                new_line["gnomAD_AN"] = "NA"
                new_line["gnomAD_AF"] = "NA"
                mgi_tsv_writer.writerow(new_line)

reader = dict_reader_unzipped('gnomad_genome_37_dict.tsv')
records = marisa_trie.RecordTrie('<fII', reader)
records.save('my_trie.marisa')
print('Finished reading records')
annotate('/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/filter_variants/756d8642be0a488ca8dcb424d5769ae2/all_variants_whitelist.756d8642be0a488ca8dcb424d5769ae2.chr.coverage.trv.noerror.tsv', records)
