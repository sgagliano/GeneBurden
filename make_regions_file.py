import sys
import argparse
from intervaltree import IntervalTree

argparser = argparse.ArgumentParser(description = 'Transforms list of LoF variants into regions file for savvy tool.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_lof', required = True, help = 'Input file with LoF variants.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = True, type = str, help = 'Output file.')


if __name__ == '__main__':
    args = argparser.parse_args()

    genes_by_chrom = {}

    with open(args.in_lof, 'rt') as ifile:
        for line in ifile:
            fields = line.rstrip().split()
            if len(fields) < 2:
                continue
            variant_name = fields[0]
            variant_gene = fields[1]
            fields = variant_name.split(':') # variant name is "chr:pos:ref:alt"
            if len(fields) != 4:
                raise Exception('Variant name {} is not correct.'.format(variant_name))
            chrom = fields[0]
            pos = int(fields[1])
            genes = genes_by_chrom.setdefault(chrom, {})
            genes.setdefault(variant_gene, []).append(pos)

    regions_by_chrom = {}
    for chrom, genes in genes_by_chrom.items():
        regions = regions_by_chrom.setdefault(chrom, IntervalTree())
        for gene, positions in genes.items():
            regions.addi(min(positions) - 100, max(positions) + 100) # +/-100 for padding, just to be conservative while extracting regions

    for chrom, regions in regions_by_chrom.items(): # merge overlapping/touching intervals
        regions.merge_overlaps(strict = False)

    with open(args.out_file, 'wt') as ofile:
        for chrom, regions in regions_by_chrom.items():
            for region in sorted(regions):
                ofile.write('{}\t{}\t{}\n'.format(chrom, region.begin, region.end))

