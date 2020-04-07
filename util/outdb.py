import numpy as np
import pandas as pd
import exdir
import collections
from pathlib import Path

ignored_keys = ['config.xenon', 'config.helium', 'config.neutrals', 'config.physical', 'config.project']

def flatten(d, parent_key='', sep='.', ignored=[]): 
        items = [] 
        for k, v in d.items(): 
            new_key = parent_key + sep + k if parent_key else k
            if new_key in ignored: continue
            if isinstance(v, collections.abc.MutableMapping): 
                items.extend(flatten(v, new_key, sep=sep, ignored=ignored).items()) 
            else: 
                items.append((new_key, v)) 
        return dict(items) 

def fetch_data(search_path='.'):
    
    search_dir = Path(search_path)
    output_files = sorted(search_dir.glob('**/hy*.exdir')) 

    metadata = []

    for f in output_files:
        exdir_file = exdir.File(f)
        attrs = exdir_file.attrs.to_dict()
        attrs['metadata']['group'] = f.relative_to(search_path).parent
        metadata_dict = flatten(attrs, ignored=ignored_keys)
        metadata.append(metadata_dict)
    
    return metadata

    
def save_to_xls(metadata, path='outdb.html'):
    df = pd.DataFrame(metadata)
    df.to_html(path)


if __name__ == "__main__":
    save_to_xls(fetch_data())

