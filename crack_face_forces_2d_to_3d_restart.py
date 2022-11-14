#!/usr/bin/env python
import sys
import os
import glob
import matplotlib.pyplot as plt
import numpy as np

def load_ASCII(file_name):
    sep =' '
    with open(file_name,'r') as f:
        read_data = f.read()
    read_data = read_data.split(sep)
    return [float(r) for r in read_data]

def write_ASCII(file_name,field):
    sep =' '
    with open(file_name,'w') as f:
        f.write(sep.join(['{:1.12e}'.format(fld) for fld in field]))    
    

def get_nbz(bname):
    with open(bname+'.in','r') as f:
        for line in f:
            if line[0]=='#':
                continue
            if 'nb_z_elements' in line:
                line = line.replace(' ','')
                print(line)
                nbz = int(line.split('nb_z_elements=')[1])
                return nbz
            
def to_3d(bname_3d,field):
    nbx = len(field)
    nbz = get_nbz(bname_3d)
    ff = np.array(list(field)*nbz).reshape((nbz,nbx)).T
    return ff.flatten()
    

###################
input_dir = 'input_data'
bname_2d = sys.argv[1]
bname_3d = bname_2d.replace('2d','3d')

if '3d' in bname_2d:
    raise RuntimeError('{} already 3d'.format(bname_2d))
if '3d' not in bname_3d:
    raise RuntimeError('{} not 3d'.format(bname_3d))

path_to_restart_2d = os.path.join(input_dir,bname_2d+'-restart')
path_to_restart_3d = os.path.join(input_dir,bname_3d+'_from2d-restart')
if os.path.exists(path_to_restart_3d):
    print('{} exists'.format(path_to_restart_3d))
else:
    os.mkdir(path_to_restart_3d)

field_paths = glob.glob(path_to_restart_2d+'/*proc0*.out')

for field_path in field_paths:
    f = load_ASCII(field_path)
    ff = to_3d(bname_3d,f)


    write_ASCII(field_path.replace(path_to_restart_2d,
                                   path_to_restart_3d),
                ff)
    
    if '_1.proc' in field_path:
        print('add dir 2 field')
        ff=np.array(ff)*0.0
        field_path = field_path.replace('_1.proc','_2.proc')
        write_ASCII(field_path.replace(path_to_restart_2d,
                                       path_to_restart_3d),
                    ff)
        

    
    #break
