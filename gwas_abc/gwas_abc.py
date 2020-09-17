import pandas as pd

def read_gwas_lead(ff, cols):
    cols_to_use = [ 'chromosome', 'position', 'snpid', 'trait' ]
    tmp = pd.read_csv(ff, compression='gzip', sep='\t')
    tmp = tmp[cols].copy()
    tmp.columns = cols_to_use
    return tmp

def read_ldblock(ff):
    return pd.read_csv(ff, sep='\t')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='gwas_abc.py', description='''
        Annotate GWAS leading variants with genes based on ABC.
    ''')
    parser.add_argument('--gwas_leading_variants', nargs='+', help='''
        A TSV.GZ table containing the GWAS leading variant.
        If should have at least columns for: chromosome, position, snpid, trait.
        List the column names of these columns in order.
    ''')
    parser.add_argument('--ld_block', help='''
        We assume that LD block file and GWAS use the same genome build.
        TSV file with region_name, chromosome, start, end columns.
    ''')
    parser.add_argument('--liftover', default=None, nargs='+', help='''
        In the case when GWAS and LD block files do not match 
        the assembly version of ABC bed, 
        the GWAS and LD block will be lifted over.
        Specify the chain file and the code path to
        misc-tools/liftover_snp/ 
    ''')
    parser.add_argument('--abc_pattern', nargs='+', help='''
        The ABC BED file pattern. Use {biosample} as wildcard.
        And add the column index of the ABC score, mapped gene in the order.
    ''')
    parser.add_argument('--abc_sample_yaml', help='''
        A YAML file listing the ABC biosamples to use for each trait.
        The format: key:value pairs with trait name as key and the list of 
        biosamples as value.
    ''')
    parser.add_argument('--bedtools_pylib', help='''
        Path to the bedtools pylib misc-tools/bedtools_pylib
    ''')
    parser.add_argument('--output', help='''
        Output filename TSV.GZ containing snpid, trait, gene, abc_max.
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    
    from tqdm import tqdm
    
    sys.path.insert(0, args.bedtools_pylib)
    from bedtools_pylib import annotate_region_with_bed, annotate_region_with_df
    
    sys.path.insert(0, '../preprocessing')
    from pylib import get_intersect, read_yaml
    
    # setup liftover
    if args.liftover_gwas is not None:
        sys.path.insert(0, args.liftover_gwas[1])
        from lib import liftover
        
    # load gwas lead variants
    logging.info('Loading GWAS leading variants.')
    df_gwas = read_gwas_lead(args.gwas_leading_variants[0], args.gwas_leading_variants[1:])
    df_gwas['start'] = df_gwas['position'] # we do 1-based and include both
    df_gwas['end'] = df_gwas['position']
    
    # annotate with LD block
    logging.info('Annotating GWAS variants with LD block.')
    df_ldblock = read_ldblock(args.ld_block)
    tmp_output = '{}_ldblock'.format(args.output)
    region_cols = ['chromosome', 'start', 'end']
    df_gwas = annotate_region_with_df(
        df_gwas[region_cols], df_ldblock[region_cols + ['region_name']], 
        suffix='_ldblock', tmp_prefix=tmp_output
    )
    
    # loop over traits
    trait_meta = read_yaml(args.abc_sample_yaml)
    trait_list = df_gwas.trait.unique().tolist()
    trait_list = get_intersect(trait_list, list(trait_meta.keys()))
    logging.info('There are {} traits to work with'.format(len(trait_list)))
    
    biosample_pattern, score_idxs = args.abc_pattern
    
    gwas_region_cols = ['chromosome_ldblock', 'start_ldblock', 'end_ldblock']
    collector = []
    for trait in tqdm(trait_list):
        
        biosamples = trait_meta[trait]
        df_gwas_sub = df_gwas[ df_gwas.trait == trait ].reset_index(drop=True)
        
        for bio in biosamples:
            
            tmp_output = '{}_x_{}_x_{}'.format(args.output, trait, bio)
            tmp = annotate_region_with_bed(
                df_gwas_sub[gwas_region_cols], 
                biosample_pattern.format(biosample=bio), 
                tmp_prefix=tmp_output,
                score_idxs
            )
            tmp.columns[-2] = 'ABC_score'
            tmp.columns[-1] = 'Mapped_gene'
            collector.append(tmp)
            
    df_abc = pd.concat(collector, axis=0)
    
    # save results
    logging.info('Saving results.')
    df_abc.to_csv(args.output, compression='gzip', sep='\t')
