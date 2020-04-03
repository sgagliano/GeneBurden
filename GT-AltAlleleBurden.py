##CMC gene-based test using genotypes, e.g. if no allele of interest in gene assign person as 0/0, if at least one allele of interest in the gene then 0/1

import sys
import gzip
import argparse
import pysam
from contextlib import closing
import itertools
from collections import OrderedDict

argparser = argparse.ArgumentParser(description = 'Create new VCF for burden test with omni dosage')
argparser.add_argument('--inGZVCF', metavar = 'filename', dest = 'inGZVCF', required = True, help = 'Input input VCF.gz')
argparser.add_argument('--outGZAnno', metavar = 'filename', dest = 'outGZAnno', required = True, help = 'Output output Anno.gz')
#argparser.add_argument('--minR2', metavar='number', dest = 'minR2', type = float, required = True, help = 'Imputatation R2 > this value')
argparser.add_argument('--maxAF', metavar='number', dest = 'maxAF', type = float, required = True, help = 'AC < this value')
argparser.add_argument('--minAF', metavar = 'number', dest = 'minAF', type = float, default = 0, help = 'Minimal AF. Default is 0.')
argparser.add_argument('--BED', metavar='filename', dest = 'Map', required = True, help = 'Tab-delimited file: chr:pos:ref:alt, gene name for variants of interest')

NEW_META_LINES = [
    '##ALT=<ID=ALT, Description="No alternate alleles at gene level test.">',
    '##INFO=<ID=N,Number=1,Type=Float,Description="Number of variants in set.">'
]

REF= "N"
ALT= "<ALT>"
QUAL= "."
FILTER= "PASS"
INFO= "."
FORMAT= "GT" #0/1 if person has at least 1 singleton in gene

def MakeDummy(inGZVCF, outGZAnno, maxAF, Map):
    genes_all = OrderedDict() # stores key-value pairs. keys = genes, values = list of variant positions to use in analysis.

    # BEGIN: This part initializes the gene sets
    with open(Map, 'r') as file:
        for line in file:
            mappings = line.rstrip().split('\t')
            variant_id = mappings[0]
	    gene_name = mappings[1]
            # add variant name to approporiate gene name (the "key"). dictionary keys store unique set of gene names
            gene_variants = genes_all.get(gene_name, None)
            if gene_variants is None:
                genes_all[gene_name] = set()
            genes_all[gene_name].add(variant_id)
    # END

    n_variants = [len(variants) for variants in genes_all.values()]
    print(f'{len(genes_all)} gene(s) and {sum(n_variants)} variant(s) in BED input file.')
    print(f'Avg. no. of variant(s) per gene = {sum(n_variants)/ float(len(genes_all))}')
    print(f'Max. no. of variant(s) in gene = {max(n_variants)} ({float(len(genes_all))} gene(s))')	  
    print(f'Min. no. of variant(s) in gene = {min(n_variants)} ({n_variants.count(min(n_variants))} gene(s))')

    # BEGIN: This part aggragates dosages
    with closing(pysam.VariantFile(inGZVCF, 'r')) as vfile, gzip.GzipFile(outGZAnno, 'w') as oz:
        # write meta and header to the output file
        for line in vfile.header.records:
            oz.write('{}'.format(line))
        for line in NEW_META_LINES:
            oz.write('{}\n'.format(line))

        IDList = list(vfile.header.samples)
        print(f'{len(IDList)} individual(s) in VCF input file.')

        oz.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(IDList)))

        for gene_name, variants in genes_all.iteritems():
            # create array of dosages for each individuals
            omni_dosage = [0] * len(IDList)
            assert sum(omni_dosage) == 0
	    omni_genotype = [0] * len(IDList)

            # get all variants positions and find min/max
            positionlist = [int(x.split(':')[1]) for x in variants]
            start = min(positionlist)
            end = max(positionlist)

            # get all variants chromosomes: must be only one unique chromosome
            chromset = set([x.split(':')[0] for x in variants])
            if len(chromset) != 1:
                raise Exception('{} gene is annotated with variants from different chromosomes!'.format(gene_name))
            chrom = chromset.pop()

            # extract all variants in that region
	    print(f'Processing {len(variants)} variant(s) {gene_name} in gene in {chrom,}:{start}-{end} region.')
            n_variants_processed = 0
            n_variants_total = 0
            for vcfrow in vfile.fetch(chrom, start, end):
                n_variants_total += 1
                if vcfrow.id is None:
		    myid = '{}:{}:{}:{}'.format(vcfrow.chrom, vcfrow.pos, vcfrow.ref, vcfrow.alts[0])
		else:
		    myid = vcfrow.id
		if myid not in variants:
		    continue
                AFvalue = vcfrow.info['AF'] if 'AF' in vcfrow.info else vcfrow.info['AC'][0]/vcfrow.info['AN'] # We assume all entries are bi-allelic
                if AFvalue >= maxAF:
                    continue
                # Standard dosage sum
                for i, sample_name in enumerate(IDList):
		    omni_dosage[i] += sum(vcfrow.samples[sample_name]['GT'])
                n_variants_processed += 1
            if n_variants_processed > 0:
                oz.write('{}\t{}\t{}\tN\t<ALT>\t.\tPASS\tN={}\tGT\t{}\n'.format(chrom, start, gene_name, n_variants_processed, '\t'.join(['0/0' if ds == 0 else '0/1' for ds in omni_dosage])))

		oz.flush()
            print 'Finished processing {} gene in {}:{}-{} region. Fetched {} variant(s) from VCF, included {} in burden.'.format(gene_name, chrom, start, end, n_variants_total, n_variants_processed)


if __name__ == '__main__':
   args = argparser.parse_args()
   MakeDummy(args.inGZVCF, args.outGZAnno, args.maxAF, args.Map)
