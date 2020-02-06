import sys
import argparse
import pysam

argparser = argparse.ArgumentParser(description = 'Extracts LoF variants from VEP annotated VCF/BCF.')
argparser.add_argument('-i', '--input', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input BCF/VCF with CSQ INFO field.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = True, type = str, help = 'Output file.')


if __name__ == '__main__':
    args = argparser.parse_args()
    genes = dict()
    ac = dict()
    with pysam.VariantFile(args.in_VCF) as ifile:
        csq_meta = ifile.header.info.get('CSQ', None)
        if csq_meta is None:
            raise Exception('No meta-information entry about CSQ INFO field found!')
        csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')
        if 'LoF' not in csq_header:
            raise Exception('No LoF sub-field found in CSQ INFO field!')


        for record in ifile.fetch():
            variant_id = '{}:{}:{}:{}'.format(record.contig, record.pos, record.ref, record.alts[0]) # assume only bi-allelic variants
            ac[variant_id] = record.info['AC'][0]
            for csq in record.info['CSQ']:
                csq = dict(zip(csq_header, csq.split('|')))
                consequences = set(csq['Consequence'].split('&'))
                if csq['BIOTYPE'] != 'protein_coding':
                    continue
                gene_id = csq['Gene']
                gene_name = csq['SYMBOL']
                lof = csq['LoF']
                if lof is not None and lof == 'HC':
                    gene = genes.setdefault(gene_id, { 'name': gene_name, 'variants': dict() })
                    gene['variants'].setdefault(variant_id, set()).update(consequences)

    with open(args.out_file, 'wt') as ofile:
        for gene_id, gene in genes.items():
            for variant, consequences in gene['variants'].items():
                ofile.write(f'{variant}\t{gene["name"]}\t{gene_id}\t{",".join(consequences)}\t{ac[variant]}\n')

