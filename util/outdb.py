import exdir
import pandas as pd
import collections, argparse
from pathlib import Path

ignored_keys = ['config.xenon', 'config.helium', 'config.neutrals', 'config.physical', 'config.project',
                'config.particles.n_factor', 'config.particles.n_factor_dsmc', 'config.particles.n_max_particles',
                'config.thruster.eta_p', 'config.thruster.mfc', 'config.time.n_steps_dsmc', 'config.time.k_sub_dsmc',
                'config.time.k_sub', 'config.time.dt', 'config.thruster.mfc', 'config.ions.q_i', 'config.ions.gas_name',
                'config.electrons.sqr_duty_cycle', 'config.electrons.q_el', 'config.electrons.m_el', 'config.boundaries.volt_1',
                'config.boundaries.volt_0', 'config.time.n_average_dsmc', 'simulation.verbosity']

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

def sort_dict(d):
    return {key: d[key] for key in sorted(d.keys(), reverse=True)}

def fetch_data(search_path='.'):
    
    search_dir = Path(search_path)
    output_files = sorted(search_dir.glob('**/hy*.exdir')) 

    metadata = []

    for f in output_files:
        exdir_file = exdir.File(f)
        attrs = exdir_file.attrs.to_dict()

        attrs['file'] = {}
        attrs['file']['parent'] = f.relative_to(search_path).parent
        
        metadata_dict = sort_dict(flatten(attrs, ignored=ignored_keys))
        metadata.append(metadata_dict)
    
    return metadata

    
def save_to_xls(metadata, path='out.xlsx'):
    df = pd.DataFrame(metadata)
    df.to_excel(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate xls file with list of outputs from Hyplex.')
    parser.add_argument("-s", "--search", help="directory to be searched.", default='..')
    parser.add_argument("-o", "--output", help="output xls file.", default='out.xlsx')
    args = parser.parse_args()

    save_to_xls(fetch_data(args.search), args.output)

