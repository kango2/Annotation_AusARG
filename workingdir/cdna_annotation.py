import argparse
import csv

def parse_gff(gff_path):
    with open(gff_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        entries = [line for line in reader if not line[0].startswith('#')]
    return entries

def parse_cds_map(cds_map_path):
    cds_map = {}
    with open(cds_map_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            # Assuming the file format: mRNA_ID, CDS_start, CDS_end
            cds_map[line[0]] = (int(line[1]), int(line[2]))
    return cds_map

def add_features(entries, cds_map):
    new_entries = []
    for entry in entries:
        if entry[2] == 'mRNA':
            chrom, _, _, start, end, _, strand, _, attributes = entry
            start = int(start)
            end = int(end)
            attr_dict = dict(item.split('=') for item in attributes.split(';') if item)
            mrna_id = attr_dict.get('ID')

            if mrna_id in cds_map:
                cds_start, cds_end = cds_map[mrna_id]

                if strand == '+':
                    utr5_start = start
                    utr5_end = start + cds_start - 1
                    utr3_start = start + cds_end
                    utr3_end = end
                else:
                    utr5_start = start + cds_end
                    utr5_end = end
                    utr3_start = start
                    utr3_end = start + cds_start - 1
                
                if utr5_start < utr5_end:
                    new_entries.append([chrom, 'GFF', 'five_prime_UTR', utr5_start, utr5_end, '.', strand, '.', attributes])
                if utr3_start < utr3_end:
                    new_entries.append([chrom, 'GFF', 'three_prime_UTR', utr3_start, utr3_end, '.', strand, '.', attributes])
                
                new_entries.append([chrom, 'GFF', 'CDS', start + cds_start, start + cds_end, '.', strand, '0', attributes])

                # Assuming exons are sorted and no overlapping occurs
                previous_end = start
                for exon in sorted((e for e in entries if e[2] == 'exon' and e[8] == attributes), key=lambda x: int(x[3])):
                    exon_start = int(exon[3])
                    if previous_end < exon_start - 1:
                        new_entries.append([chrom, 'GFF', 'intron', previous_end + 1, exon_start - 1, '.', strand, '.', attributes])
                    previous_end = int(exon[4])

    return new_entries

def main(input_file, output_file, cds_map_file):
    entries = parse_gff(input_file)
    cds_map = parse_cds_map(cds_map_file)
    new_entries = add_features(entries, cds_map)
    
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(new_entries)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map cDNA coordinates to genome annotations using a CDS to cDNA map.')
    parser.add_argument('--input', type=str, help='Input GFF file path')
    parser.add_argument('--output', type=str, help='Output GFF file path')
    parser.add_argument('--cds_map', type=str, help='CDS to cDNA map file path')
    
    args = parser.parse_args()
    main(args.input, args.output, args.cds_map)
