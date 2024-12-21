import os
from multiprocessing import Process, Lock, Manager
from band_structure import create_KPOINTS_mass
from siman.header import _update_configuration, db
from siman.calc_manage import smart_structure_read, add_loop, res_loop
from siman.set_functions import read_vasp_sets
import time
from project_sets import band_pack
from effective_mass import calc_eff_mass_st3, generate_kpoints, parse_EIGENVAL_VASP, find_extr
from pymatgen.io.vasp import BSVasprun
import project_sets
from project_sets import band_pack
from band_structure import create_KPOINTS_mass


def mass_calk(st_name, num_of_calc, lock):
    file_name = st_name + '.rho'
    input_geo_file = st_name + '/' + st_name + '.POSCAR'

    with lock:
        add_loop(file_name, 'sp', 1, ise_new=f'mass_{num_of_calc}', inherit_option='full', savefile='ocxe',
                 input_geo_file = input_geo_file, it_folder = st_name, override=1, run=2)

    while res_loop(file_name + '.if', f'mass_{num_of_calc}', 1, up='x')== ({}, []):
        time.sleep(20)


def eff_mass_tensor(st_name='Si', st='st3', band=4, k_extr=(0, 0, 0), name_KPOINTS='/KPOINTS_mass', num_of_calc=0, lock=None):
    folder_name = st_name + '.rho.if'
    folder_name_mass = folder_name + f'.mass_{num_of_calc}'
    folder_name_k_mass = st_name + name_KPOINTS
    folder_name_band = folder_name + '.local_band'

    vasprun_mass_path = f"{st_name}/{folder_name_mass}/1.vasprun.xml"
    kpoints_band_path = f"{st_name}/{folder_name_band}/KPOINTS"

    kpoints = generate_kpoints(k_extr, 'st3', 0.01)
    create_KPOINTS_mass(folder_name_k_mass, kpoints)

    with lock:
        band_pack['kpoints_file'] = folder_name_k_mass
        read_vasp_sets([(f'mass_{num_of_calc}', 'sp', band_pack, 'override')])

    mass_calk(st_name, num_of_calc, lock)

    if os.path.exists(vasprun_mass_path):
        os.remove(vasprun_mass_path)
        print(f"Файл {vasprun_mass_path} удален.")
    else:
        print(f"Файл {vasprun_mass_path} не существует.")

    db[folder_name, f'mass_{num_of_calc}', 1].get_file('1.vasprun.xml')

    v = BSVasprun(vasprun_mass_path)
    bs_mass = v.get_band_structure(kpoints_filename=kpoints_band_path)

    return calc_eff_mass_st3(bs_mass, st, band)


def process_point(st_name, scheme, band, k, i, file_name_mass, lock):
    import cProfile
    import pstats
    def worker():
        try:
            print(f"Calculations at point {i}: {k}")
            mass = eff_mass_tensor(st_name=st_name, st=scheme, band=band, k_extr=k, name_KPOINTS=f"KPOINTS_{i}", num_of_calc=i, lock=lock)

            with lock:
                with open(file_name_mass, 'a') as file:
                    file.write(f'\n\n\nCalculations at point {i}:\n{k}\n\n')
                    if mass is not None:
                        for row in mass:
                            file.write(" ".join(f"{num:.3f}" for num in row) + '\n')
                    else:
                        file.write("Failed to compute mass tensor.\n")
        except Exception as e:
            print(f"Error processing point {i}: {e}")

    profiler = cProfile.Profile()
    profiler.enable()
    worker()
    profiler.disable()

    profile_file = f"profile_process_{i}.prof"
    with open(profile_file, "w") as f:
        ps = pstats.Stats(profiler, stream=f)
        ps.dump_stats(profile_file)

    print(f"Profile saved for process {i} to {profile_file}")


def loop_exst_mass_multiproc(st_name='Si', file_name_mass='eff_mass_tensor_multi.txt', band=4, num_processes=4):
    try:
        folder_name = st_name + '.rho.if'
        folder_name_band = st_name + '/' + folder_name + '.local_band'

        script_dir = os.getcwd()
        eigenval_path = os.path.join(script_dir, f"{folder_name_band}/1.EIGENVAL")

        E, all_points = parse_EIGENVAL_VASP(eigenval_path, band)
        extr_k, extr_E = find_extr(E, all_points)

        with open(file_name_mass, 'w') as file:
            file.write(f'Effective mass tensors for {st_name}\n')

        lock = Lock()

        processes = []

        for i, k in enumerate(extr_k):
            process = Process(target=process_point, args=(st_name, 'st3', band, k, i, file_name_mass, lock))
            processes.append(process)
            process.start()

            if len(processes) >= num_processes:
                for p in processes:
                    p.join()
                processes = []

        for p in processes:
            p.join()

    except Exception as e:
        print(f"Error in loop_exst_mass: {e}")


if __name__ == '__main__':
    loop_exst_mass_multiproc(st_name='Si', file_name_mass='eff_mass_tensor_multi_new.txt', band=4, num_processes=10)


