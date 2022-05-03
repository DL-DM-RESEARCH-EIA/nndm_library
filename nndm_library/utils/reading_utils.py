from asyncore import read
from posixpath import split
from reprlib import recursive_repr
import sys
import os
from pathlib import Path
from matplotlib.container import BarContainer
import numpy as np
import glob
import pandas as pd
import pickle
import uproot
import tqdm
import numpy as np
import tqdm
import functools
import pylhe
import re 

def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

class Constants:
    # Chosen from standard
    ELECTRON_ID = 11

    # TO DO: find standar name for particles
    #   chosen in arbitrary way to the moment
    ETA_ID = 1
    PION_ID = 2

    # Definig map dictionary
    name_to_id = {'electron' : ELECTRON_ID, 'eta' : ETA_ID, 'pion' : PION_ID}

class ReadFileBase:
    def __init__(self, path):
        self.path = path

    def extract_params_from_path(self):
        """
        path: path of the file, which is assumed to contain the values of the eps and the mass to be extracted from it
        format example is as follows: eta_decay_events_mk_0.38_eps2_5.404557191441203e-07.lhe. In more general term, it i
            {particle_name}_{param1}_{value1}_{param2}_{value2}_{param3}_{value3}*.lhe

        return:
            dictionary with all extracted data
        """
        splitted = os.path.splitext(os.path.basename(self.path))[0].split("_")
        res_dict = {}

        res_dict['particle_type'] = Constants.name_to_id[splitted[0]]

        # going over splitted values and creating dictionary
        for i in range(len(splitted)):
            print(splitted[i], splitted[i].isnumeric())
            if isfloat(splitted[i]) and (splitted[i - 1] != splitted[0]):
                res_dict[splitted[i - 1]] = splitted[i]

        return res_dict


# This class read the output data relating to the electron scaterings
class ReadRoot(ReadFileBase):
    
    def __init__(
        self, path: str, output_base_tree="treeout", pattern_output="first",
        output_base_middle_branch = "/e/out",
        leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"], recursive=False,
        files_dir=None, relabel_events=True 
    ):
        """
        Class to read the labeled data coming in ROOT format

        Parameters: path: string. 
                        The direction to the root file containing the event information
                    output_base_name: string. 
                        Name bas of the first node of the tree that has the data. For instance, if the base name 
                        is treeout, there options could be treeout1, treeout2, ...., treeoutN. 
                    pattern_output: string 
                        The idea is this parameter define a methodology to choose from the possible first nodes 
                        that have a given output_base_name. As an example, first would choose treeoout1 in the 
                        example before.
                    output_base_middle_branch: string 
                        middle branch that goes after the selected first node chosen by  the output pattern. If this 
                        variable is "e/out", following the example the tree to consult at the moment would be 
                        treeout1/e/out/
                    leafs: list 
                        what are the leafs to exaplore in the actual branch. If out.a is the ouput name for the a momenta, 
                        giving a list [out.x, out.y] will give the data to consult. That is, treeout1/e/out/out.x and
                        treeout1/e/out/out.y
                    data: dataframe that contains all the branches specified by the parameters below as name for the rows, and the
                        data associated with such names goes in the columns.
                    files_dir: dictionary
                        it saves all the files from which data was readed to construct the data frame and it shows the integer
                        associeted with it that identifies it in the dataframe constructed.
                    relabel_events: bool 
                        there is an id for each possible event. For instance a collision have an id for it and two sub ids for the
                            particle that interact in it. When relabel_events is True, the values of id are associated unequivocally
                            with each event.  

        Methods
        ------------
            No public methods.

        """
        ReadFileBase.__init__(self, path) 

        self.output_base_tree = output_base_tree
        self.pattern_output = pattern_output
        self.output_base_middle_branch = output_base_middle_branch
        self.leafs = leafs

        self.recursive = recursive
        self.relabel_events = relabel_events

        self.files_dir = files_dir

        self._read()

    def _sanity_check(self):
        """
        This function cheks that the parameters entered, or the combination of them are acutally valid.
        """
        pass

    def _read(self):
        """
        General method to read independently from the initialization from the class.
        """
        self._read_single_file_safe()        
        self._read_recursive_fliles()

    def _read_recursive_fliles(self):
        ##########
        # Read when self.path is a directory
        if os.path.isdir(self.path):
            # files_dir is no mor None
            self.files_dir = {}

            file_list = self._file_list()
            if not len(file_list):
                print("Please check the existense of .root files in the specified directory %s or use recursive=True"
                      % self.path)
                exit(0)
            else:
                self._append_df_of_file_list(file_list)

    def _file_list(self):
        # Recursive reading means that will find all the .root files inside in 
        #  any of the subsequen directories
        if (self.recursive):
            file_list = [file for sub_dir in os.walk(self.path) for file in glob.glob(os.path.join(sub_dir[0], '*.root'))]
        # When non-recursive it will simply try to read the .root files in such a directory
        else:
            file_list = glob.glob(self.path + "*.root")

        return file_list

    def _read_single_file_safe(self):
        #########
        # Read when self.path is a file
        if os.path.isfile(self.path):
            _, ext = os.path.splitext(self.path)
            if (ext == ".root"):
                if (self.pattern_output == "first"):
                    # create list of number termination
                    self.data = self._read_single_file()
            else:
                # TO DO: create standard errors
                print("please enter a valid .root file")
                exit(0)



    def _read_single_file(self, path=""):
        """
        Return dataframe associated with a single file.
        """

        # Check is path is no specified. It is specified when working with directories to go over
        #   each file.

        # opening correct file for each case (directory and file)
        if path:
            self.tree = uproot.open(path)
        else:
            self.tree = uproot.open(self.path)

        # choose keys with the correct base name
        keys = self.tree.keys()
        filter_by_base_name = [key for key in keys if self.output_base_tree in key]
        
        if (self.pattern_output == "first"):
            # chose the first one
            numeric_values = [int(filter_by_base_name[i].split(";", 1)[1]) for i in range(len(filter_by_base_name))]
            index_min = numeric_values.index(min(numeric_values))
        
            # get the data we are interested
            final_branch = filter_by_base_name[index_min] + self.output_base_middle_branch
            data_frame = self.tree[final_branch].arrays(self.leafs, library="pd")

        return data_frame

    def _append_df_of_file_list(self, file_list):
        self.data = pd.DataFrame()
        last_index = 0
        for i, file_path in enumerate(file_list):
            df_to_append = self._read_single_file(path=file_path)
            df_to_append["path"] = i
            self.files_dir[i] = file_path

            # Make first row is a unique identifier
            if (self.relabel_events):
                df_to_append = df_to_append.rename(functools.partial(index_mapper, last_index), 
                                                axis=0, level=df_to_append.index.names[0])
                last_index = df_to_append.index.tolist()[-1]
                if type(last_index) != int:
                    last_index = last_index[0]

        self.data = pd.concat([self.data, df_to_append], axis = 0)

def index_mapper(index, last_index):
    """
    When appending data frames, the first index of the new dataframe must be shifted
        to end up with unique ids. 
    """

    if type(index) == int:
        return index + last_index
    else:
        index_values = list(index)
        index_values[0] = index_values[0] + last_index
        return tuple(index_values)

class ReadLhe(ReadFileBase):
    def __init__(self, path, particle_ids=None, var_of_interest=None, outgoing=False, recursive=False,
        files_dir=None, relabel_events=True, verbose=1):
        """
        Class to read the labeled data coming in lhe format
        Parameters: variables_interest: list
                        names of the variables to extract from the lhe. eg. ["energy","angle"], 
                        ["energy","px","py"] .... By default: ["energy","angle"]
                    path: string. 
                        The direction to the file containing the event information
        """

        ReadFileBase.__init__(self, path) 
        self.var_of_interest = var_of_interest

        self.verbose = verbose
        self.outgoing = outgoing
        self.particle_ids = particle_ids

        self.data = self._init_data()
        self._read_single_file()
        self.data = pd.DataFrame.from_dict(self.data)

    def _init_data(self, path=''):
        data = {}
        if self.var_of_interest is not None:
            for var in [v for v in self.var_of_interest]:
                data[var] = np.array([])
            return data

        if self.var_of_interest is None and os.path.isfile(self.path):
            self._check_right_path_is_open(path)

            for obj in self.file:
                for particle in obj.particles:
                    if isinstance(particle, pylhe.LHEParticle):
                        for name in particle.fieldnames:
                            data[name] = np.array([])
                        return data

    def _check_right_path_is_open(self, path=""):
        # Check is path is no specified. It is specified when working with directories to go over
        #   each file.

        # opening correct file for each case (directory and file)
        if path:
            self.file = pylhe.readLHE(path)
        else:
            self.file = pylhe.readLHE(self.path)

    def _filtrate_outgoing(self, particle):
        """ Function returns the if particle is outgoing ot not when outgoing option is on, 
            otherwise filter does not do anything."""
        if self.outgoing:
            return particle.status == 1
        else:
            return True

    def _filtrate_by_id(self, particle):
        """
        When particle_ids are given filtering is evaluated
        """
        if self.particle_ids is not None:
            return particle.id in self.particle_ids
        else:
            return True

    def _add_particle_data(self, particle):
        """
        Add particle data to the dict_var output
        """
        for name in self.data.keys():
            self.data[name] = np.append(self.data[name], getattr(particle, name))

    def _read_single_file(self, path=""):
        """
        Return dataframe associated with a single file.
        """

        self._check_right_path_is_open()

        for obj in self.file:
            for particle in obj.particles:
                if self._filtrate_outgoing(particle):
                    # filtrating by type
                    if self._filtrate_by_id(particle):
                        # saving vars of interest
                        self._add_particle_data(particle)

    def get(self):
        return self.data

    def get_momenta(self):
        """
        name: name of the file, which is assumed to contain the values of the eps and the mass to be extracted from it
        format example is as follows: eta_decay_events_mk_0.38_eps2_5.404557191441203e-07.lhe

        return: momentum and energy
            px, py, pz, energy
        """
        return self.data["energy"], self.data["px"], self.data["py"], self.data["pz"]

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
