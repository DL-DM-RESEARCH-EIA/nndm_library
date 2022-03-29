from asyncore import read
import sys
import os
from pathlib import Path
import numpy as np
import glob
import pandas as pd
import pickle
import uproot
import tqdm
import xml.etree.ElementTree as ET
import numpy as np
import tqdm


class ReadFileBase:
    def __init__(self, path):
        self.path = path

    def extract_params_from_path(self):
        """
        path: path of the file, which is assumed to contain the values of the eps and the mass to be extracted from it
        format example is as follows: eta_decay_events_mk_0.38_eps2_5.404557191441203e-07.lhe

        return:
            particle type, mass, eps2
        """
        splitted = self.path.split("_")
        mk = float(splitted[-3])
        eps2 = float(splitted[-1].split(".lhe")[0])
        if "eta" in self.path:
            type_particle = 0
        elif "pion" in self.path:
            type_particle = 1
        return type_particle, mk, eps2


class ReadLhe(ReadFileBase):
    def __init__(self, path, variables_interest=["energy", "ps"], catch=1, verbose=1):
        """
        ReadLhe: Class
        -------------------------------------------------------------
        * Parameters:
        - variables_interest: List. Variable to extract from the lhe.
        egs. ["energy","angle"], ["energy","px","py"] .... By default: ["energy","angle"]
            KeyArgs:
            - path: String. Path to search for the txt with the path for lhe's
            - catch: Int. Identificate the particle. 1 is the outgoing electron (measured) . Default 1
        * Output
        - DataFrame with all readed events.
        """
        ReadFileBase.__init__(self, path) 

        self.doc = ET.parse(path)
        # catch is used to select the electrons that go away
        self.catch = catch
        self.cartessian = ["px", "py", "pz"]

        self.root = self.doc.getroot()

        self.dic_var = {}
        if "ps" in variables_interest:
            variables_interest = variables_interest + self.cartessian
        for var in [v for v in variables_interest if v != "ps"]:
            self.dic_var[var] = np.array([])

        self.verbose = verbose
        self._read()

    def _getAngle(self, momentums, axis=2):
        """
        Calculate the angle of the particles starting from a list of the form [px, py, pz].
        This with respect to the "axis" element.
        """
        momentumTotal = np.linalg.norm(momentums)
        angule = np.arccos(momentums[axis] / momentumTotal)
        return angule

    def _read(self):
        """
        procced with certain LHE reading
        """
        if self.verbose:
            apply_function = tqdm
        else:
            apply_function = lambda x: x

        for bloque in apply_function(self.root):
            Evento = []
            if bloque.tag == "event":  # Evaluate if data comes from valid event
                datos = bloque.text
                for x in datos.splitlines():
                    a = x.split()
                    if len(a) == 13:
                        Evento.append(a)
                energia_electrones = [float(Evento[1][9]), float(Evento[3][9])]
                momentum_electrones = [
                    [float(Evento[1][i]) for i in [6, 7, 8]],
                    [float(Evento[3][i]) for i in [6, 7, 8]],
                ]  # [[px,py,pz]]
                for var in self.dic_var.keys():
                    if var in self.cartessian:  # Append all momentums
                        # self.dic_var[var] = np.append(self.dic_var[var],momentum_electrones[self.catch])
                        dif = self.cartessian.index(var)
                        self.dic_var[self.cartessian[dif]] = np.append(
                            self.dic_var[self.cartessian[dif]],
                            momentum_electrones[self.catch][dif],
                        )
                    elif "angle" == var:  # Append the angle
                        self.dic_var[var] = np.append(
                            self.dic_var[var],
                            self._getAngle(momentum_electrones[self.catch], axis=2),
                        )
                    elif "energy" == var:  # Append the energy :: 'energy'
                        self.dic_var[var] = np.append(
                            self.dic_var[var], energia_electrones[self.catch]
                        )
                    elif "P" == var:  # Get magnitude of Four Momentum :: P
                        self.dic_var[var] = np.append(
                            self.dic_var[var],
                            np.linalg.norm(momentum_electrones[self.catch]),
                        )

    def get(self):
        return self.dic_var

    def get_momenta(self):
        """
        name: name of the file, which is assumed to contain the values of the eps and the mass to be extracted from it
        format example is as follows: eta_decay_events_mk_0.38_eps2_5.404557191441203e-07.lhe

        return: momentum and energy
            px, py, pz, energy
        """
        return self.dic_var["energy"], self.dic_var["px"], self.dic_var["py"], self.dic_var["pz"]

# This class read the output data relating to the electron scaterings
class ReadRoot:
    def __init__(
        self, path: str, out_channel="v_e_scattering/onantinuelepton10125.root"
    ):
        """
        DataReader: Class to read the labeled data coming in ROOT format
        -----------------------------------------------------------------

        Input:
        - path: string. The direction to the root file containing the event information
        - out_channel: string. Path to the onantinuelepton file to the tree output
        """

        self.tree = uproot.open(path)
        # Read only what it it needed for time purposes, note that this does not works
        #   when there are several branches
        branches = ["out.t", "out.x", "out.y", "out.z", "out._mass"]

        if path == out_channel:
            self.df = self.tree["treeout;2"]["e/out"].arrays(branches, library="pd")
        else:
            self.df = self.tree["treeout;3"]["e/out"].arrays(branches, library="pd")



class FilesManipulator:
    def __init__(self, path, verbose=0):
        self.path = path
        self.verbose = verbose
        # Create structure to save the data in the .lhe files
        # Here we have pictorical description of the data structure
        #   First, a list of ints is [int, int, ...] == [(int)]
        #   So a list of a list of floats is: [[(float)], [(float)], ...] == [( [(float)] )]
        #   {id: [(int)], typ: [(str)], mk: [(float)], eps2: [(float)], px: [[(float)], [(float)], ...], py: [( [(float)] )], pz: [( [(float)] )] }
        # Note that despite momentum and energy being list, they will contain a list of list
        #   to be more specific, a list of arrays, where each array correspons to a param point
        self.scan = {
            "id": [],
            "typ": [],
            "mk": [],
            "eps2": [],
            "name": [],
            "px": [],
            "py": [],
            "pz": [],
            "energy": [],
        }

    def fill_up_scan(self):
        file_list = glob.glob(self.path)

        # initialize scan dictionary as described above going over all the names
        for id, name in enumerate(file_list):
            read_file = ReadLhe(name, verbose=self.verbose)
            typ, mk, eps2 = read_file.extract_params_from_path()
            self.scan["id"].append(id), self.scan["typ"].append(typ), self.scan["mk"].append(mk)
            self.scan["eps2"].append(eps2), self.scan["name"].append(os.path.basename(name))
            energy, px, py, pz = read_file.get_momenta()
            for var in ["px", "py", "pz", "energy"]:
                self.scan[var].append(eval(var))

    def save_scan(self, save_name="complete_task.pickle"):
        with open(save_name, "wb") as f:
            pickle.dump(self.scan, f, protocol=pickle.HIGHEST_PROTOCOL)
