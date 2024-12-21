from band_structure import create_KPOINTS, charge_calk, band_calk, create_band_plot
import psutil
import time
import cProfile
from threading_eff_mass_loop import loop_exst_mass_thread
from regular_eff_mass_loop import loop_exst_mass
from multiproc_eff_mass_loop import loop_exst_mass_multiproc
from mp_api.client import MPRester
from input_file_maker import maker
from siman.header import _update_configuration
from siman.calc_manage import smart_structure_read
from siman.database import read_database
import os
import re

api_key = "your_Materials_Project_key" #fill
_update_configuration('simanrc.py')
read_database()

def initial_calculations(mp_id, debug=False):
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            material_ids=silicon
        )
    for doc in docs:
        st = doc.structure
        st_name = re.sub(r'([()\\[\]])', r'', st.reduced_formula)
        path = maker(st, st_name) #создаем POSCAR
        siman_st = smart_structure_read(f'{st_name}.POSCAR')

        create_KPOINTS(st_name, st)
        os.chdir(f'{path}')
        charge_calk(st_name, siman_st)
        band_calk(st_name)

        if debug:
            create_band_plot(st_name, mode='total') 


def check_resorces_regular_script():
    def track_resources():
        cpu_usage = psutil.cpu_percent(interval=1)
        memory_info = psutil.virtual_memory()
        return cpu_usage, memory_info.percent

    start_time = time.time()
    cpu_before, memory_before = track_resources()

    loop_exst_mass(st_name='Si', band=4)

    cpu_after, memory_after = track_resources()
    end_time = time.time()

    print(f"Время выполнения: {end_time - start_time} секунд")
    print(f"CPU до: {cpu_before}%, CPU после: {cpu_after}%")
    print(f"Память до: {memory_before}%, Память после: {memory_after}%")


def check_resorces_threading_script():    
    def track_resources():
        cpu_usage = psutil.cpu_percent(interval=1)
        memory_info = psutil.virtual_memory()
        return cpu_usage, memory_info.percent
    
    start_time = time.time()
    cpu_before, memory_before = track_resources()
    
    loop_exst_mass_thread(st_name='Si', file_name_mass='eff_mass_tensor_treads.txt', band=4, num_threads=10)
    
    cpu_after, memory_after = track_resources()
    end_time = time.time()
    
    print(f"Время выполнения: {end_time - start_time} секунд")
    print(f"CPU до: {cpu_before}%, CPU после: {cpu_after}%")
    print(f"Память до: {memory_before}%, Память после: {memory_after}%")


def check_resorces_multiproc_script():
    def track_resources():
        # Мониторим загрузку процессора и памяти
        cpu_usage = psutil.cpu_percent(interval=1)
        memory_info = psutil.virtual_memory()
        return cpu_usage, memory_info.percent
    
    start_time = time.time()
    cpu_before, memory_before = track_resources()
    
    loop_exst_mass_multiproc(st_name='Si', file_name_mass='eff_mass_tensor_multi_new.txt', band=4, num_processes=10)
    
    cpu_after, memory_after = track_resources()
    end_time = time.time()
    
    print(f"Время выполнения: {end_time - start_time} секунд")
    print(f"CPU до: {cpu_before}%, CPU после: {cpu_after}%")
    print(f"Память до: {memory_before}%, Память после: {memory_after}%")


if __name__ == '__main__':
    silicon = ['mp-149']
    initial_calculations(silicon) #preparatory calculations to obtain the band structure, the results are stored in the folder Si

    print("Profiling a regular function:")
    cProfile.run("loop_exst_mass(st_name='Si', band=4)")

    print("Profiling a treading function:")
    cProfile.run("loop_exst_mass_thread(st_name='Si', file_name_mass='eff_mass_tensor_treads.txt', band=4, num_threads=10)")

    print("Profiling a multiproc function:")
    cProfile.run("loop_exst_mass_multiproc(st_name='Si', file_name_mass='eff_mass_tensor_multi_new.txt', band=4, num_processes=10)")

    check_resorces_regular_script()
    check_resorces_threading_script()
    check_resorces_multiproc_script()

        

