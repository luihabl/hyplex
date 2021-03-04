import yaml
from pathlib import Path
import sys
from ctxt import * 

def load_yaml(path):
    return yaml.load(Path(path).read_text(),Loader=yaml.FullLoader)



# print(CColor.F_Green + k + " as key not in d2", "\n" + CBase.END)
def find_dict_diff2(dict_1, dict_2, dict_1_name='d1', dict_2_name='d2', path=""):

    for k in dict_1.keys():
        if k not in dict_2:
            print(CFormatting.Bold + f"{path} :" + CBase.END)
            print(f"    {CColor.F_Yellow}{k}{CBase.END} key not in {dict_1_name}")
        else:
            if isinstance(dict_1[k], dict) and isinstance(dict_2[k], dict):
                find_dict_diff2(dict_1[k],dict_2[k],'d1','d2', path + '/' + k)
            else:
                if dict_1[k] != dict_2[k]:
                    print(CFormatting.Bold + f"{path} :" + CBase.END)
                    print(f"    {CColor.F_Red}(d1) {k} : {dict_1[k]}{CBase.END}")
                    print(f"    {CColor.F_Green}(d2) {k} : {dict_2[k]}{CBase.END}")

    for k in dict_2.keys():
        if k not in dict_1:
            print(CFormatting.Bold + f"{path} :" + CBase.END)
            print(f"    {CColor.F_Yellow}{k}{CBase.END} key not in {dict_1_name}")


if __name__ == "__main__":
    
    if len(sys.argv) < 3: 
        print('Not enough arguments')
        exit()
    
    dict1 = load_yaml(sys.argv[1])
    dict2 = load_yaml(sys.argv[2])

    find_dict_diff2(dict1, dict2)
