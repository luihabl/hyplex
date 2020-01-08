import os, shutil, configparser


config_file_path = 'input/config/config.ini'
output_path = 'output'
run_file_path = 'run'

def change_val(file_path, group, key, new_value):
    
    file = configparser.ConfigParser()
    file.optionxform = str
    file.read(file_path)
    
    old_value_line = file[group][key]
    
    old_value_list = old_value_line.split(' ')
    old_value_list[0] = str(new_value)
    new_value_line = ' '.join(old_value_list)
    
    file[group][key] = new_value_line
    
    with open(file_path, 'w') as configfile:
        file.write(configfile)

def modify_config_file(file_path, config_dict):
    for group, group_dict in config_dict.items():
        for config_key, config_val in group_dict.items():
            change_val(file_path, group, config_key, config_val)
            
def rename_output(output_folder, suffix):
    try: 
        os.mkdir(os.path.join(output_folder, suffix))
    except:
        pass
    for file in os.listdir(output_folder):
        if file.endswith(".h5"):
            new_file = file
            # new_file = file.split('.')
            # new_file[-1] = '_' + suffix + '.h5'
            # new_file = ''.join(new_file)
            # os.rename(os.path.join(output_folder, file), os.path.join(output_folder, new_file))
            destination_folder = os.path.join(output_folder, suffix)
            shutil.move(os.path.join(output_folder, new_file), os.path.join(destination_folder, new_file))

def run_cases(cases, run_file, config_file, output_folder):

    for case_name, case_dict in cases.items():
        
        modify_config_file(config_file, case_dict)
        print("running: " + case_name)
        os.system('mpirun -np 1 ' + run_file)
        os.system('sleep 1')
        rename_output(output_folder, case_name)

if __name__ == '__main__':
    
    cases = {
        'a2p5_25mhz': {'thruster': {'V_SB': 1003.938, 'FREQ': 25e6}},
        'a2p5_20mhz': {'thruster': {'V_SB': 1003.938, 'FREQ': 20e6}},
        'a2p5_15mhz': {'thruster': {'V_SB': 1003.938, 'FREQ': 15e6}},
        'a2p5_10mhz': {'thruster': {'V_SB': 1003.938, 'FREQ': 10e6}},
        'a3_30mhz': {'thruster': {'V_SB': 1003.027, 'FREQ': 30e6}},
        'a3_25mhz': {'thruster': {'V_SB': 1003.027, 'FREQ': 25e6}},
        'a3_20mhz': {'thruster': {'V_SB': 1003.027, 'FREQ': 20e6}},
        'a3_15mhz': {'thruster': {'V_SB': 1003.027, 'FREQ': 15e6}},
        'a3_10mhz': {'thruster': {'V_SB': 1003.027, 'FREQ': 10e6}},
    }

    run_cases(cases, run_file_path, config_file_path, output_path)






