import sys
import pandas as pd

sys.path.insert(0, '../../preprocessing')
from pylib import load_list

def load_tags(in_list, in_tag):
    '''
    Prioritize the loading of in_list.
    If in_list is None, try loading in_tag.
    '''
    if in_list is not None:
        o = load_list(in_list)
    else:
        if in_tag is not None:
            o = [ in_tag ]
        else:
            raise ValueError
    return o

def load_gwas(ff, cols):
    tmp = pd.read_csv(ff, compression='gzip', sep='\t')
    tmp = tmp[ cols ].copy()
    tmp.columns = [ 'variant_id', 'pval', 'chromosome', 'position' ]
    return tmp

def write_output(df, filename):
    df.to_csv(filename, index=False, compression='gzip', sep='\t')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='construct_count_table.py', description='''
        Obtain 2-by-2 table for a list of GWASs and biosamples.
    ''')
    parser.add_argument('--gwas_pattern', help='''
        The GWAS file pattern (use {gwas_tag} as wildcard).
    ''')
    parser.add_argument('--bed_pattern', help='''
        The BED file pattern (use {bed_tag} as wildcard).
    ''')
    parser.add_argument('--gwas_tag', default=None, help='''
        The tag name of GWAS file.
    ''')
    parser.add_argument('--bed_tag', default=None, help='''
        The tag name of BED file.
    ''')
    parser.add_argument('--gwas_list', default=None, help='''
        A list of GWAS tags.
    ''')
    parser.add_argument('--bed_list', default=None, help='''
        A list of BED tags.
    ''')
    parser.add_argument('--pval_cutoff_list', default=None, type=float, nargs='+', help='''
        A list of pval cutoffs. (default: 5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
    ''')
    parser.add_argument('--snp_list', default=None, help='''
        A list of SNPs for output.
    ''')
    parser.add_argument('--col_meta', nargs='+', help='''
        Column name in GWAS file.
        List the names of: SNP ID, p-value, chromosome, position,
        in the order.
    ''')
    parser.add_argument('--output', help='''
        A TSV.GZ table containing the 2-by-2 table for all GWAS / BED pairs.
    ''')
    parser.add_argument('--bedtools_pylib', help='''
        Path to bedtools pylib script.
    ''')
    parser.add_argument('--liftover_gwas', default=None, nargs='+', help='''
        If want to liftover GWAS result, set the chain file and 
        path to liftover python lib.
    ''')
    args = parser.parse_args()
 
    import logging, time, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    from tqdm import tqdm
    
    sys.path.insert(0, args.bedtools_pylib)
    from bedtools_pylib import intersect_with_bed
   
    # setup liftover
    if args.liftover_gwas is not None:
        sys.path.insert(0, args.liftover_gwas[1])
        from lib import liftover
        
    # load pval cutoffs
    if args.pval_cutoff_list is None:
        pval_cutoffs = [5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    else:
        pval_cutoffs = args.pval_cutoff_list 
    
    # load tags
    gwas_tags = load_tags(args.gwas_list, args.gwas_tag)
    bed_tags = load_tags(args.bed_list, args.bed_tag)
    logging.info('Working on {} GWASs and {} BED files.'.format(len(gwas_tags), len(bed_tags)))
    
    # load SNP list
    if args.snp_list is not None:
        snp_list = load_list(args.snp_list)
    else:
        snp_list = None
    
    # build the table
    collector = []
    for gwas_tag in gwas_tags:
        
        logging.info('Working on GWAS {}'.format(gwas_tag))
        # load gwas and limit to variants in snp_list
        gwas_file = args.gwas_pattern.format(gwas_tag=gwas_tag)
        df_gwas = load_gwas(gwas_file, args.col_meta)
        df_gwas = df_gwas[ ~ df_gwas.pval.isna() ].reset_index(drop=True)
        if snp_list is not None:
            df_gwas = df_gwas[ df_gwas.variant_id.isin(snp_list) ].reset_index(drop=True)
        if args.liftover_gwas is not None:
            tmp = liftover(df_gwas.chromosome, df_gwas.position, args.liftover_gwas[0])
            df_gwas.chromosome = tmp.liftover_chr
            df_gwas.position = tmp.liftover_pos
            df_gwas = df_gwas[ df_gwas.position > 0 ].reset_index(drop=True)
        
        for bed_tag in tqdm(bed_tags):
            
            bed_file = args.bed_pattern.format(bed_tag=bed_tag)
            df_gwas_x_bed = intersect_with_bed(
                df_gwas, bed_file, 
                inplace=False,
                tmp_prefix='{}_x_{}_x_{}'.format(args.output, gwas_tag, bed_tag)
            )
            n_total = df_gwas_x_bed.shape[0]
            n_in = df_gwas_x_bed.annot.sum(axis=0)
            n_cutoff_in_bed = []
            n_cutoff_total = []
            for p_cutoff in pval_cutoffs:
                tmp = df_gwas_x_bed[ df_gwas_x_bed.pval < p_cutoff ]
                n_cutoff_total.append(tmp.shape[0])
                n_cutoff_in_bed.append(tmp.annot.sum(axis=0))
            df_counts = pd.DataFrame(
                [[gwas_tag, bed_tag] + [ n_total, n_in ] + n_cutoff_total + n_cutoff_in_bed], 
                columns=['GWAS', 'BED', 'N_total', 'N_total_in_bed'] + [ 'N_total_at_pval_lt_{}'.format(i) for i in pval_cutoffs ] + [ 'N_in_bed_at_pval_lt_{}'.format(i) for i in pval_cutoffs ]
            )
            collector.append(df_counts)
    
    df_out = pd.concat(collector, axis=0)
    write_output(df_out, args.output)
    
    
