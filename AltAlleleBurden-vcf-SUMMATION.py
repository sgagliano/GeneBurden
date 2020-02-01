##Updates:
##Use AF rather than MAF, b/c it is the alt allele that is the LoF
##Apply minumum default AF, needs to be greater than >0.00001 b/c low allele counts are unstable --> removed this b/c found that didn't make a difference in the association test
##Match by variant ID (chr:pos:ref:alt) in LoF variant list and VCF

# In[98]:

import sys
sys.path.insert(0, '/net/snowwhite/home/sarahgag/.local/lib/python2.7/site-packages/')
import gzip
import argparse
import pysam
from contextlib import closing
import itertools
from collections import OrderedDict

argparser = argparse.ArgumentParser(description = 'Create new VCF for burden test with omni dosage')
argparser.add_argument('--inGZVCF', metavar = 'filename', dest = 'inGZVCF', required = True, help = 'Input input VCF.gz')
argparser.add_argument('--outGZAnno', metavar = 'filename', dest = 'outGZAnno', required = True, help = 'Output output Anno.gz')
argparser.add_argument('--minR2', metavar='number', dest = 'minR2', type = float, required = True, help = 'Imputatation R2 > this value')
argparser.add_argument('--maxAF', metavar='number', dest = 'maxAF', type = float, required = True, help = 'AF < this value')
#argparser.add_argument('--minAF', metavar = 'number', dest = 'minAF', type = float, default = 0, help = 'Minimal AF. Default is 0.')
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
FORMAT= "DS"

def MakeDummy(inGZVCF, outGZAnno, minR2, maxAF, Map):
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
    print '{} gene(s) and {} variant(s) in BED input file.'.format(len(genes_all), sum(n_variants))
    print 'Avg. no. of variant(s) per gene = {}'.format(sum(n_variants)/ float(len(genes_all)))
    print 'Max. no. of variant(s) in gene = {} ({} gene(s))'.format(max(n_variants), n_variants.count(max(n_variants)))
    print 'Min. no. of variant(s) in gene = {} ({} gene(s))'.format(min(n_variants), n_variants.count(min(n_variants)))

    # BEGIN: This part aggragates dosages
    with closing(pysam.VariantFile(inGZVCF, 'r')) as vfile, gzip.GzipFile(outGZAnno, 'w') as oz:
        # write meta and header to the output file
        for line in vfile.header.records:
            oz.write('{}'.format(line))
        for line in NEW_META_LINES:
            oz.write('{}\n'.format(line))

        IDList = list(vfile.header.samples)
        print '{} individual(s) in VCF input file.'.format(len(IDList))

        oz.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(IDList)))

        for gene_name, variants in genes_all.iteritems():
            # create array of dosages for each individuals
            omni_dosage = [0] * len(IDList)
            assert sum(omni_dosage) == 0

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
	    print 'Processing {} variant(s) in {} gene in {}:{}-{} region.'.format(len(variants), gene_name, chrom, start, end)
            n_variants_processed = 0
            n_variants_total = 0
            for vcfrow in vfile.fetch(chrom, start, end):
                n_variants_total += 1
                if vcfrow.id is None:
                    raise Exception('Input VCF/BCF file must have non-empty variant ID column.')
                if vcfrow.id not in variants:
                    continue
		#create own id instead of vcfrow.id (e.g. to match id in --BED):
		#if vcfrow.id is None:
                #    myid = '{}:{}:{}:{}'.format(vcfrow.chrom, vcfrow.pos, vcfrow.ref, vcfrow.alts[0])
                #else:
                #    myid = vcfrow.id
                #if myid not in variants:
                #    continue
                AFvalue = vcfrow.info['AF'] # We assume all entries are bi-allelic
                R2value = float(vcfrow.info['R2'])
                if R2value <= minR2 or AFvalue >= maxAF:
                    continue
                # Standard dosage sum
                for i, sample_name in enumerate(IDList):
                    omni_dosage[i] += vcfrow.samples[sample_name]['DS']
                n_variants_processed += 1
            if n_variants_processed > 0:
                oz.write('{}\t{}\t{}\tN\t<ALT>\t.\tPASS\tN={}\tDS\t{}\n'.format(chrom, start, gene_name, n_variants_processed, '\t'.join([str(ds) for ds in omni_dosage])))

		oz.flush()
            print 'Finished processing {} gene in {}:{}-{} region. Fetched {} variant(s) from VCF, included {} in burden.'.format(gene_name, chrom, start, end, n_variants_total, n_variants_processed)


if __name__ == '__main__':
   args = argparser.parse_args()
   MakeDummy(args.inGZVCF, args.outGZAnno, args.minR2, args.maxAF, args.Map)
