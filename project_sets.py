# -*- coding: utf-8 -*-
import copy

from siman import header
from siman.set_functions import InputSet, inherit_iset, make_sets_for_conv, init_default_sets

"""List of user VASP sets obtained on inheritance principle """
"""Syntax:  ("set_new", "set_old", {"param1":value1, "param2":value2, ...})
        - set_new - name of new set
        - set_old - name of base set used for creating new set
        - {} -      dictionary of parameters to be updated in set_new
"""


static_run_packet = {'NSW':0, 'EDIFF'     : 6e-06, 'NELM':50}

band_pack = {'LORBIT': 11, 'KSPACING': 0.10}

user_vasp_sets = [
('sp', 'static', {'KSPACING': 0.35, 'NPAR':None,'ISTART':None,'NELM':None,'PREC':None,'ALGO':None,'KGAMMA':None
                         ,'ENCUT':600,'SIGMA':None,'LPLANE':None, 'LREAL':'.FALSE.'}) ]

mass_pack = {'LORBIT': 11, 'KSPACING': 0.10}
