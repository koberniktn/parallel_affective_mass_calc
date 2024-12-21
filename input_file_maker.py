import os

def maker(st, st_name):
    folder_path = st_name
    os.makedirs(folder_path, exist_ok=True)
    path = os.getcwd()
    os.chdir(f'{path}/{folder_path}')
    st.to(f"{st_name}.POSCAR", folder_path)
    return path
