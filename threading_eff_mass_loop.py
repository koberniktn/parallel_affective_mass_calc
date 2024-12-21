import os
import threading
from regular_eff_mass_loop import eff_mass_tensor
from effective_mass import parse_EIGENVAL_VASP, find_extr
import project_sets
from project_sets import band_pack
from band_structure import create_KPOINTS_mass, mass_calk


def process_point(st_name, scheme, band, k, i, file_name_mass, lock):
    try:
        print(f"Calculations at point {i}: {k}")
        mass = eff_mass_tensor(st_name=st_name, st=scheme, band=band, k_extr=k, name_KPOINTS=f"KPOINTS_{i}", num_of_calc=i)

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


def loop_exst_mass_thread(st_name='Si', file_name_mass='eff_mass_tensor.txt', band=4, num_threads=4):
    try:
        folder_name = st_name + '.rho.if'
        folder_name_band = st_name + '/' + folder_name + '.local_band'

        script_dir = os.getcwd()
        eigenval_path = os.path.join(script_dir, f"{folder_name_band}/1.EIGENVAL")

        E, all_points = parse_EIGENVAL_VASP(eigenval_path, band)
        extr_k, extr_E = find_extr(E, all_points)

        with open(file_name_mass, 'w') as file:
            file.write(f'Effective mass tensors for {st_name}\n')

        lock = threading.Lock()

        threads = []
        for i, k in enumerate(extr_k):
            thread = threading.Thread(target=process_point, args=(st_name, 'st3', band, k, i, file_name_mass, lock))
            threads.append(thread)
            thread.start()

            if len(threads) >= num_threads:
                for t in threads:
                    t.join()
                threads = []

        for t in threads:
            t.join()

    except Exception as e:
        print(f"Error in loop_exst_mass: {e}")


if __name__ == '__main__':
    loop_exst_mass_thread(st_name='Si', file_name_mass='eff_mass_tensor_treads_check.txt', band=4, num_threads=4)



