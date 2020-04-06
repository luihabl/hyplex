import strictyaml as yaml
import collections
import copy
from pathlib import Path

def ismap_yaml(d):
    return d.is_mapping()

def ismap_collections(d):
    return isinstance(d, collections.abc.Mapping)

def update(d, u, ismap=ismap_yaml):
    if not ismap(d):
        return u
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v, ismap)
        else:
            d[k] = v
    return d

def load_yaml(path):
    return yaml.load(Path(path).read_text())

def save_yaml(y, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(y.as_yaml())

def create_cases(cases, input_file, output_path, modtype='numbering'):
    
    output_dir = Path(output_path)
    original_config = load_yaml(input_file)

    for i, (k, v) in enumerate(cases.items()):
        config = copy.deepcopy(original_config)
        config = update(config, update(v, {'simulation':{'job_name': k}}, 
                                                    ismap=ismap_collections))
        original_filename = Path(input_file).stem
        original_file_extension = Path(input_file).suffix

        if modtype == 'numbering':
            mod = '_' + str(i)
        
        if modtype == 'key':
            mod = '_' + k

        new_filename = original_filename + mod + original_file_extension
        save_yaml(config, output_dir / new_filename)
        print(output_dir / new_filename)


if __name__ == "__main__":
    
    cases = {
        'a2p5_25mhz': {'thruster': {'v_sb': '1003.938', 'freq': '25e6'}}
        # 'a2p5_20mhz': {'thruster': {'v_sb': '1003.938', 'freq': '20e6'}},
        # 'a2p5_15mhz': {'thruster': {'v_sb': '1003.938', 'freq': '15e6'}},
        # 'a2p5_10mhz': {'thruster': {'v_sb': '1003.938', 'freq': '10e6'}}
    }

    create_cases(cases, '../input/config/config.yaml', './config', 'numbering')






