import pathlib, gzip, os
import yaml

def get_intersect(l1, l2):
    s1 = set(l1)
    return list(s1.intersection(set(l2)))

def get_arg(str_, config):
    if str_ not in config:
        return ''
    else:
        return '--{} {}'.format(str_, config[str_])

def load_list(ff):
    o = []
    with open(ff, 'r') as f:
        for i in f:
            i = i.strip()
            o.append(i)
    return o

def write_list(ofile, abc, cellcol):
    with gzip.open(abc, 'rt') as f:
        for i in f:
            i = i.strip().split('\t')
            if cellcol not in i:
                raise ValueError(f'{cellcol} not in header! Exit')
            cellidx = i.index(cellcol) + 1    
            break
    cmd = f"zcat {abc} | tail -n +2 | awk '{{print ${cellidx}}}' | sort | uniq > {ofile}"
    os.system(cmd)

def read_yaml(ff):
    with open(ff, 'r') as f:
        o = yaml.safe_load(f)
    return o
    
def get_prefix(ss):
    return '.'.join(ss.split('.')[:-1])

# for split_abc_prediction
def get_all(abc_file, celltype_col):
    mylist = abc_file + '.all_celltype_list'
    if pathlib.Path(mylist).is_file():
        pass
    else:
        write_list(mylist, abc_file, celltype_col)
    cell_list = load_list(mylist)
    outlist = [ '{outdir}/{name_tag}.' + cell + '.original.tsv.gz' for cell in cell_list ]
    outlist += [ '{outdir}/{name_tag}.' + cell + '.cleanup.tsv.gz' for cell in cell_list ]   
    return outlist

    
