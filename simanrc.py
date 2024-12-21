"""
User-related settings for siman
"""

local_path = '/Users/irinabrodyagina/PycharmProjects/find_k-path/'
PATH2POTENTIALS = local_path+'POT_GGA_PAW_PBE'
AUTO_UPDATE_DB = True
pmgkey = "your_Materials_Project_key" #fill


"""Cluster settings"""
DEFAULT_CLUSTER = 'ccmm'
user = "user_name" #fill
host = "your_host" #fill

from siman.header import CLUSTERS


CLUSTERS['ccmm'] = {
'address':user+'@'+host,
'vasp_com':'vasp_std',
'homepath':'your home local path', #fill
'schedule':'SLURM',
'corenum':1,
'modules':
'\nulimit -s unlimited\n\
'
}


