
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='extract_relevant_cols.py', description='''
        Extract columns from TSV.GZ.
    ''')
    parser.add_argument('--input', help='''
        Assume it is TSV.GZ.
    ''')
    parser.add_argument('--col_yaml', help='''
        Column names in a YAML
    ''')
    parser.add_argument('--output', help='''
        Output filename
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
    
    import sys, os, gzip
    sys.path.insert(0, '../')
    from pylib import read_yaml
    
    col_dict = read_yaml(args.col_yaml)
    cols_to_extract = list(col_dict.values())
    
    col_idxs = []
    with gzip.open(args.input, 'rt') as f:
        for i in f:
            i = i.strip().split('\t') 
            for kk in cols_to_extract:
                if kk not in i:
                    raise ValueError(f'{kk} not in header. Exit.')
                else:
                    col_idxs.append(str(i.index(kk) + 1))
            break
    
    cmd = "zcat {} | cut -f {} | gzip > {}".format(args.input, ','.join(col_idxs), args.output)
    os.system(cmd)
    
