import os
import sys
import subprocess
from aflow_sym import Symmetry
from siman.bands import plot_bands
from siman.header import _update_configuration, db
from siman.calc_manage import smart_structure_read, add_loop, res_loop
from siman.set_functions import read_vasp_sets
from siman.database import read_database
import project_sets
from siman.dos_functions import plot_dos
from my_band_plot import plot_bands
import time
from siman.calc_manage import get_structure_from_matproj
# from mp_api.client import MPRester
from pymatgen.core import Structure
from siman import header
from siman.core.calculation import Calculation
import matplotlib.pyplot as plt
import re
import subprocess
from project_sets import band_pack


api_key = "your_Materials_Project_key" #fill
_update_configuration('simanrc.py')
read_database()


def create_KPOINTS_mass(k_folder, k_path):
    with open(k_folder, 'w', newline='') as f:
        f.write('Explicit k-point list\n')
        f.write('{:} ! full number of points\n'.format(k_path[0]))
        f.write('r\n')  # now only reciprocal are supported
        for point in k_path[1:]:
            f.write('{:6.3f} {:6.3f} {:6.3f} {:d}\n'.format(point[1], point[2], point[3], point[0]))


def create_KPOINTS(structure_name, st, method="AFLOW"):
    if method == "AFLOW": #using AFLOW software
        input_POSCAR = structure_name + ".POSCAR"
        subprocess.run(f"aflow --kpath < {input_POSCAR} > KPOINTS", shell=True, env=env)
        k_path = band_set_from_KPOINTS()
        os.remove('KPOINTS')
    elif method == 'reference': #using your own KPOINTS file
        k_path = band_set_from_KPOINTS()
    else:
        print("Неверно задан метод расчета k-path")
        return 0

    band_pack['k_band_structure'] = k_path
    band_pack['kpoints_file'] = 1
    read_vasp_sets([('local_band', 'sp', band_pack, 'override')])


def band_set_from_KPOINTS():
    k_path = []
    pattern = r'(?i)Gamma'
    with open('KPOINTS', 'r') as file:
        for line in file:
            if "!" in line:
                data = line.replace("\\", "").split()
                data = [re.sub(pattern, 'G', elem) for elem in data]
                if 'grids' in data:
                    k_path.append(int(data[0])) # number of points per line
                else:
                    point_name = data[-1]
                    if k_path and (len(k_path) == 1 or (len(k_path) >= 1 and k_path[-1][0] != point_name)):
                        k_path.append((point_name, float(data[0]), float(data[1]), float(data[2])))
    return k_path


def charge_calk(st_name, st=None): #send task to the cluster
    file_name = st_name + '.rho'
    it_folder = st_name + '/charge'

    add_loop(file_name, 'sp', 1, input_st=st, it_folder=it_folder, run=2) #chardge_density
    while res_loop(file_name, 'sp', 1, up='x') == ({}, []):
       time.sleep(60)


def band_calk(st_name): #send task to the cluster
    file_name = st_name + '.rho'
    input_geo_file = st_name + '/' + st_name + '.POSCAR'

    add_loop(file_name , 'sp', 1, ise_new='local_band', inherit_option='full', savefile='ocxe', input_geo_file = input_geo_file, it_folder = st_name, override=1, run=2)
    while res_loop(file_name + '.if', 'local_band', 1, up='x')== ({}, []):
        time.sleep(40)

    res_loop(file_name + '.if', 'local_band', 1, up='x')
    db[file_name + '.if', 'local_band', 1].get_file('1.vasprun.xml')
    db[file_name + '.if', 'local_band', 1].get_file('1.EIGENVAL')


def mass_calk(st_name, num_of_calc): #send task to the cluster
    file_name = st_name + '.rho'
    input_geo_file = st_name + '/' + st_name + '.POSCAR'

    add_loop(file_name, 'sp', 1, ise_new=f'mass_{num_of_calc}', inherit_option='full', savefile='ocxe',
             input_geo_file = input_geo_file, it_folder = st_name, override=1, run=2)

    while res_loop(file_name + '.if', f'mass_{num_of_calc}', 1, up='x')== ({}, []):
        time.sleep(20)

    db[file_name + '.if', f'mass_{num_of_calc}', 1].get_file('1.vasprun.xml')


def create_band_plot(st_name, mode='total', **kwargs):
    band_file_name = st_name + '.rho.if'

    with open(f'{st_name}/{band_file_name}.local_band/KPOINTS') as file:
        for line in file.readlines():
            print(line)

    while res_loop(st_name + '.rho' + '.if', 'local_band', 1, up='x') == ({}, []):
        time.sleep(3)

    plot_bands(st_name, vasprun_bands=f'{st_name}/{band_file_name}.local_band/1.vasprun.xml',
               kpoints=f'{st_name}/{band_file_name}.local_band/KPOINTS',
               folder='out/',
               ylim=(-12, 8), mode=mode, **kwargs)















