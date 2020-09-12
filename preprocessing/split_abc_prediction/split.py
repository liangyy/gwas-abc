import pandas as pd

def read_abc_file(filename):
    return pd.read_csv(filename, compression='gzip', sep='\t')

def update_region(df, start, end, shrink_size):
    if shrink_size is None:
        pass
    else:
        # start pos
        df[start] = df[start] + shrink_size
        # end pos
        df[end] = df[end] - shrink_size
        df = df[ df[start] < df[end] ].reset_index(drop=True)
    return df

def filter_by_abc(df, abc, cutoff):
    return df[ df[abc] >= cutoff ].reset_index(drop=True)        
    
def write_to_file(df, filename):
    df.to_csv(filename, compression='gzip', sep='\t', index=False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='split.py', description='''
        Split ABC prediction file by biosample.
        For each biosample, output the original score and a clean-up version.
        In the clean-up version, it allows to apply range shrinkage and filter 
        by ABC score.
        It requires reading in the whole ABC prediction file into the memory.
    ''')
    parser.add_argument('--input', nargs='+', help='''
        ABC prediction file.
        Assume it is TSV.GZ.
        And show the column names of chromosome, start, end, cell type, ABC score 
        in the order.
    ''')
    parser.add_argument('--output_prefix', help='''
        The prefix of the output.
    ''')
    parser.add_argument('--shrink_region_by', default=None, type=int, help='''
        Shrink the region by []. 
        It matters only for the clean-up file.
    ''')
    parser.add_argument('--remove_abc_lt', default=None, type=float, help='''
        Remove ABC score that is less than []. 
        It matters only for the clean-up file.
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
    
    logging.info('Load ABC score file.')
    infile, chrm_col, start_col, end_col, celltype_col, abc_col = args.input
    df_abc = read_abc_file(infile)
    celltype_list = df_abc[celltype_col].unique().tolist()
    logging.info('There are {} cell types'.format(len(celltype_list)))
    
    for celltype in tqdm(celltype_list):
        original_out = '{}.{}.original.tsv.gz'.format(args.output_prefix, celltype)
        cleanup_out = '{}.{}.cleanup.tsv.gz'.format(args.output_prefix, celltype)
        df_sub = df_abc[ df_abc[celltype_col] == celltype ].reset_index(drop=True)
        
        # write the original file
        write_to_file(df_sub, original_out)
        
        # work on the cleanup file
        df_sub = update_region(df_sub, start_col, end_col, args.shrink_region_by)
        df_sub = filter_by_abc(df_sub, abc_col, args.remove_abc_lt)
        write_to_file(df_sub, cleanup_out) 
    
