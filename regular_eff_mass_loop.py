from effective_mass import find_extr, parse_EIGENVAL_VASP, calc_eff_mass_st3, generate_kpoints
from pymatgen.io.vasp import BSVasprun
from siman.set_functions import read_vasp_sets
import os
import project_sets
from project_sets import band_pack
from band_structure import  create_KPOINTS_mass, mass_calk


def eff_mass_tensor(st_name='Si', st='st3', band=4, k_extr=(0, 0, 0), name_KPOINTS='/KPOINTS_mass', num_of_calc=0, direct=None):
    varset = read_vasp_sets(project_sets.user_vasp_sets, override_global=0)
    script_dir = os.getcwd()

    folder_name = st_name + '.rho.if'
    folder_name_mass = folder_name + f'.mass_{num_of_calc}'
    folder_name_k_mass = st_name + name_KPOINTS
    folder_name_band = folder_name + '.local_band'

    vasprun_mass_path = f"{st_name}/{folder_name_mass}/1.vasprun.xml"
    kpoints_band_path = f"{st_name}/{folder_name_band}/KPOINTS"

    kpoints = generate_kpoints(k_extr, h=0.01)

    create_KPOINTS_mass(folder_name_k_mass, kpoints)
    band_pack['kpoints_file'] = folder_name_k_mass
    read_vasp_sets([(f'mass_{num_of_calc}', 'sp', band_pack, 'override')])

    mass_calk(st_name, num_of_calc)

    if os.path.exists(vasprun_mass_path):
        os.remove(vasprun_mass_path)
        print(f"Файл {vasprun_mass_path} удален.")
    else:
        print(f"Файл {vasprun_mass_path} не существует.")

    
    v = BSVasprun(vasprun_mass_path)
    bs_mass = v.get_band_structure(kpoints_filename=kpoints_band_path)

    return calc_eff_mass_st3(bs_mass, st, band)


def loop_exst_mass(st_name='Si', file_name_mass='eff_mass_tensor.txt', band=4):
    folder_name = st_name + '.rho.if'
    folder_name_band = st_name + '/' + folder_name + '.local_band'

    script_dir = os.getcwd()
    eigenval_path = os.path.join(script_dir, f"{folder_name_band}/1.EIGENVAL")

    E, all_points = parse_EIGENVAL_VASP(eigenval_path, band)
    extr_k, extr_E = find_extr(E, all_points)

    with open(file_name_mass, 'w') as file:
        file.write(f'Effective mass tensors for {st_name}')
        for i in range(len(extr_k)):
            print("\nCalculations in the point:")
            k = extr_k[i]
            print(i, k)
            file.write('\n\n\nCalculations in the point:')
            file.write(f'\n\n{k}\n\n')
            mass = eff_mass_tensor(st_name='Si', st='st3', band=band, k_extr=k, name_KPOINTS='KPOINTS' + "_" + str(i), num_of_calc=i)
            for row in mass:
                file.write(" ".join(f"{num:.3f}  " for num in row) + '\n')

if __name__ == '__main__':
    loop_exst_mass(st_name='Si', band=4)