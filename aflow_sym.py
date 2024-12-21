import json
import subprocess
import os

class Symmetry:

    def __init__(self, aflow_executable='aflow'):
        self.aflow_executable = aflow_executable


    def _check_command(self, command):
        try:
            return subprocess.check_output(
                self.aflow_executable + command,
                shell=True
            )
        except subprocess.CalledProcessError:
            print ("Error aflow executable not found at: " + self.aflow_executable)


    def _aflow_command(self, input_file, command, tol=None, magmoms=None):
        fpath = os.path.realpath(input_file.name)
        output = ''

        if tol:
            command += '=' + str(tol)
        if magmoms:
            command += ' --magmom=' + magmoms

        output = self._check_command(
            command + ' --print=json --screen only'
            + ' < ' + fpath
        )

        return output


    def get_symmetry(self, input_file, tol=None, magmoms=None):
        resss = []
        #res = self._aflow_command(input_file, 'aflow --aflowSYM', tol, magmoms)
        #for line in res:
            #resss.append(json.loads(line))
        #return resss


    def get_edata(self, input_file, tol=None, magmoms=None):
        return json.loads(self._aflow_command(input_file, 'aflow --edata', tol, magmoms))


    def get_sgdata(self, input_file, tol=None, magmoms=None):
        return json.loads(self._aflow_command(input_file, 'aflow --sgdata', tol, magmoms))


    def get_k_path(self, input_file, tol=None, magmoms=None):
        return json.loads(self._aflow_command(input_file, 'aflow --kpath', tol, magmoms))