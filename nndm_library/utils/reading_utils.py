from asyncore import read
from importlib.resources import path
from posixpath import split
from reprlib import recursive_repr
import os
import numpy as np
import glob
import pandas as pd
import pickle
import uproot
import numpy as np
# import tqdm #TO DO
import functools
import pylhe
import json

from nndm_library.utils.utils import ColumnFunctionsMixin

def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def index_mapper(index, last_index):
    """
    When appending data frames, the first index of the new dataframe must be shifted
        to end up with unique ids. 
    """

    if type(index) == int:
        return index + last_index + 1
    else:
        index_values = list(index)
        index_values[0] = index_values[0] + last_index + 1
        return tuple(index_values)

class Constants:
    # Chosen from standard
    ELECTRON_ID = 11
    ETA_ID = 221
    PION_ID = 111
    DM_ID = 50

    # Definig map dictionary
    name_to_id = {'electron' : ELECTRON_ID, 'eta' : ETA_ID, 'pion' : PION_ID, 'dm' : DM_ID}

class ReadFileBase(ColumnFunctionsMixin):
    """
    TODO: end in case no files is foound in recursive

    Class to read the labeled data coming in a format like the following:\n
    data1 data2 data3\n
    v11   v12   v13\n
    .     .     .  \n
    .     .     .  \n
    .     .     .  \n
    vn1   vn2   vn3\n

    where data1, data2, ... respresent names and vij, a value in the given i row and column dataj. 

    :param path: the direction to the file containing all the events information.
    :type path: string

    :param recursive: read all the .lhe files found in all paths inside a given files_dir
    :type recursive: bool

    :param ext: extension of the files to read
    :type ext: str

    :type relabel_events: bool 
    :param relabel_events: there is an id for each possible event. For instance a collision have an id for it and two sub ids for the
        particle that interact in it. When relabel_events is True, the values of id are associated unequivocally with each event.

    :var data: dataframe with the read events
    :type data: dataframe

    :var files_dir: directory with the name of the files read and its id
    :type data: dict
    """
    
    def __init__(self, path, recursive=False, ext='.txt', relabel_events=True):
        self.path = path
        self.recursive = recursive
        self.ext = ext
        self.relabel_events = relabel_events
        
        self._read()

    def extract_params_from_path(self):
        """
        Format is as follows: {particle_name}_{param1}_{value1}_{param2}_{value2}_{param3}_{value3}*.lhe
        An example would be eta_decay_events_mk_0.38_eps2_5.404557191441203e-07.lhe.

        :return: dictionary with all extracted data
        """
        res_dict = {}

        for name in self.files_dir.values():
            splitted = os.path.splitext(os.path.basename(name))[0].split("_")

            if res_dict.get('particle_type'):
                if type(res_dict['particle_type']) == int:
                    for key in res_dict.keys():
                        res_dict[key] = [res_dict[key]]

                res_dict['particle_type'].append(Constants.name_to_id[splitted[0]])

                # going over splitted values and creating dictionary
                for i in range(len(splitted)):
                    if isfloat(splitted[i]) and (splitted[i - 1] != splitted[0]):
                        res_dict[splitted[i - 1]].append(float(splitted[i]))
            else:
                res_dict['particle_type'] = Constants.name_to_id[splitted[0]]

                # going over splitted values and creating dictionary
                for i in range(len(splitted)):
                    if isfloat(splitted[i]) and (splitted[i - 1] != splitted[0]):
                        res_dict[splitted[i - 1]] = float(splitted[i])

        return res_dict

    def assign_process_weights(self, file_weights):
        weights_dir = json.load(open(file_weights))

        self.data['weight'] = 0

        # Go over files for which we know the weight and modify
        #   and assign correct values for weitgh in dataframe if
        #   they share such procedence
        total = np.array(list(weights_dir.values())).sum()
        for file_name, weigth in weights_dir.items():
            for index_path, file_path in self.files_dir.items():
                if file_name in file_path:
                    # normalize by proportion
                    self.data.loc[self.data.path == index_path, 'weight'] = weigth / total


    def _fill_files_dir(self):
        self.files_dir = {}
        file_list = self._file_list()

        for i, file_path in enumerate(file_list):
            self.files_dir[i] = file_path

    def _read_func(self, path):
        return path
    
    def _check_right_path_is_open(self, path=""):
        """
        Check if path is no specified to return self.path. It is specified when working with directories 
            to go over each file, thenit would retuen path.
        """

        # opening correct file for each case (directory and file)
        if path:
            self.file = path
        else:
            self.file = self.path

        self.read_object = self._read_func(self.file)

    def _read_single_file_safe(self, path=''):
        self._check_right_path_is_open(path)

        if os.path.isfile(self.file):
            _, ext = os.path.splitext(self.file)
            if ext == self.ext:
                self._read_single_file()
            else:
                # TO DO: create standard errors
                print("please enter a valid %s file" % (self.ext))
                exit(0)

    def _read_single_file(self):
        self.data = pd.read_csv(self.file, sep='\s+')

    def _file_list(self):
        # Recursive reading means that will find all the .root files inside in 
        #  any of the subsequen directories
        if self.recursive:
            file_list = [file for sub_dir in os.walk(self.path) for file in glob.glob(os.path.join(sub_dir[0], '*' + self.ext))]
        # When non-recursive it will simply try to read the .root files in such a directory
        elif os.path.isdir(self.path):
            file_list = glob.glob(self.path + '*' + self.ext)
        elif os.path.isfile(self.path):
            file_list = [self.path]
        else:
            print("please enter a valid path to a file or directory; %s was not found" % self.path)
            exit(0)

        return file_list

    def _append_df_of_file_list(self):
        res_dataframe = pd.DataFrame()
        last_index = 0
        for i, file_path in self.files_dir.items():
            self._read_single_file_safe(path=file_path)
            df_to_append = self.data 
            df_to_append["path"] = i

            # Make first row is a unique identifier
            if self.relabel_events:
                df_to_append = df_to_append.rename(functools.partial(index_mapper, last_index), 
                                                axis=0, level=df_to_append.index.names[0])
                last_index = df_to_append.index.tolist()[-1]
                if type(last_index) != int:
                    last_index = last_index[0]

            res_dataframe = pd.concat([res_dataframe, df_to_append], axis = 0)
        self.data = res_dataframe

    def _read_recursive_files(self):
        """
        Read when self.path is a directory, either all datain the given path or all
            the files found recursively from the given path.
        """
        if os.path.isdir(self.path):
            # files_dir is no mor None

            if len(self.files_dir) == 0:
                print("Please check the existense of %s files in the specified directory %s or use recursive=True"
                      % (self.ext, self.path) )
                exit(0)
            else:
                self._append_df_of_file_list()

    def _read(self):
        """
        General method to read independently from the initialization from the class.
        """
        self._fill_files_dir()
        self._read_single_file_safe()
        self._read_recursive_files()


class ReadLhe(ReadFileBase):
    """
    Class to read the data coming in lhe format. By default it will read all the particles. Filters used apply to the such default data.

    :param path: the direction to the file containing all the events information.
    :type path: string

    :param partcile_ids: ids of the particles to extract from the file according to the pdg, 
        By default: None, which means exctract all the particles.
    :type partcile_ids: list of integers

    :param var_of_interest: names of the variables to extract from the lhe. eg. ["e","angle"], ["e","px","py"] 
        .... By default: None, which means exctract all the variables.
    :type var_of_interest: list of strings

    :param outgoing: filtrate to obtain all the outgoing particles
    :type outgoing: bool

    :param files_dir: directory where the files are to be found
    :type files_dir: string

    :param recursive: read all the .lhe files found in all paths inside a given files_dir
    :type recursive: bool

    :param verbose: show progress reading all the .lhe files
    :type verbose: bool
    
    :var data: dataframe with the read events
    :type data: dataframe

    :var files_dir: directory with the name of the files read and its id
    :type data: dict
    """
    def __init__(self, path, particle_ids=None, var_of_interest=None, outgoing=False, 
        recursive=False, relabel_events=True, verbose=1):
        
        self.path = path 
        self.ext = ".lhe"
        self.verbose = verbose
        self.outgoing = outgoing
        self.particle_ids = particle_ids
        self.var_of_interest = var_of_interest

        self.recursive = recursive
        self.relabel_events = relabel_events

        self.data = self._init_data()

        ReadFileBase.__init__(self, path, ext=self.ext, recursive=self.recursive, relabel_events=self.relabel_events) 
        self.data = pd.DataFrame.from_dict(self.data)

    def _init_data(self, path=''):
        data = {}
        # Initialize when var_of_interest interest is passed as
        #   parameter
        if self.var_of_interest is not None:
            for var in [v for v in self.var_of_interest]:
                data[var] = np.array([])
            return data

        # Initialize when var_of_interest is None and there is a valid file
        #  here is assumed that the same info is given for each event
        elif self.var_of_interest is None and os.path.isfile(self.path):
            self._check_right_path_is_open(path)

            for obj in self.read_object:
                for particle in obj.particles:
                    if isinstance(particle, pylhe.LHEParticle):
                        for name in particle.fieldnames:
                            data[name] = np.array([])
                        return data

        elif path and os.path.isfile(path):
            self._check_right_path_is_open(path)

            for obj in self.read_object:
                for particle in obj.particles:
                    if isinstance(particle, pylhe.LHEParticle):
                        for name in particle.fieldnames:
                            data[name] = np.array([])
                        return data

    def _read_func(self, path):
        return pylhe.readLHE(path)

    def _filtrate_outgoing(self, particle):
        """ 
        Function returns the if particle is outgoing ot not when outgoing option is on, 
            otherwise filter does not do anything.
        """
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

    def _read_single_file(self):
        """
        Return dataframe associated with a single file.
        """
        self.data = self._init_data(self.file)

        for obj in self.read_object:
            for particle in obj.particles:
                if self._filtrate_outgoing(particle):
                    # filtrating by type
                    if self._filtrate_by_id(particle):
                        # saving vars of interest
                        self._add_particle_data(particle)
    
        self.data = pd.DataFrame.from_dict(self.data)

# This class read the output data relating to the electron scaterings
class ReadRoot(ReadFileBase):
    """
    Class to read the labeled data coming in ROOT format. By default it assumes values for output_base_tree, pattern_output,
    output_base_middle_branch, and leafs. This is for a fast reading.

    :param path: the direction to the root file(s)
    :type path: str 

    :param output_base_name: Name bas of the first node of the tree that has the data. For instance, if the base name 
        is treeout, there options could be treeout1, treeout2, ...., treeoutN. 
    :type output_base_name: str 

    :param pattern_output: The idea is this parameter define a methodology to choose from the possible first nodes 
        that have a given output_base_name. As an example, first would choose treeoout1 in the 
        example before.
    :type pattern_output: str 

    :type output_base_middle_branch: str 
    :param output_base_middle_branch: middle branch that goes after the selected first node chosen by  the output pattern. If this 
        variable is "e/out", following the example the tree to consult at the moment would be treeout1/e/out/.

    :type leafs: list of strings
    :param leafs: what are the leafs to exaplore in the actual branch. If out.a is the ouput name for the a momenta, 
        giving a list [out.x, out.y] will give the data to consult. That is, treeout1/e/out/out.x and
        treeout1/e/out/out.y
    
    :type files_dir: list of strings
    :param files_dir: it saves all the files from which data was readed to construct the data frame and it shows 
        the integer associeted with it that identifies it in the dataframe constructed.
    :type relabel_events: bool 
    :param relabel_events: there is an id for each possible event. For instance a collision have an id for it and two sub ids for the
        particle that interact in it. When relabel_events is True, the values of id are associated unequivocally with each event.

    """

    def __init__(
        self, path: str, output_base_tree="treeout", pattern_output="first",
        output_base_middle_branch = "/e/out",
        leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"], recursive=False,
        files_dir=None, relabel_events=True 
    ):
        self.output_base_tree = output_base_tree
        self.pattern_output = pattern_output
        self.output_base_middle_branch = output_base_middle_branch
        self.leafs = leafs

        self.recursive = recursive
        self.relabel_events = relabel_events

        self.ext = ".root"

        self.files_dir = files_dir

        ReadFileBase.__init__(self, path, ext=self.ext, relabel_events=self.relabel_events, recursive=self.recursive)

    def _read_func(self, path):
        if os.path.isfile(path):
            return uproot.open(path)

    def _read_single_file(self):
        """
        Return dataframe associated with a single file.
        """

        if os.path.isfile(self.file):
            # choose keys with the correct base name
            keys = self.read_object.keys()
            filter_by_base_name = [key for key in keys if self.output_base_tree in key]
            
            if self.pattern_output == "first":
                # chose the first one
                numeric_values = [int(filter_by_base_name[i].split(";", 1)[1]) for i in range(len(filter_by_base_name))]
                index_min = numeric_values.index(min(numeric_values))
            
                # get the data we are interested
                final_branch = filter_by_base_name[index_min] + self.output_base_middle_branch
                data_frame = self.read_object[final_branch].arrays(self.leafs, library="pd")

            self.data = data_frame


class FilesManipulator:
    """
    Create structure to save the data in the .lhe files in the path given
    Here we have pictorical description of the data structure:  
        First, a list of ints is [int, int, ...] == [(int)].  So a list of a list of floats is: [[(float)], [(float)], ...] == [( [(float)] )]
        {id: [(int)], typ: [(str)], mk: [(float)], eps2: [(float)], px: [[(float)], [(float)], ...], py: [( [(float)] )], pz: [( [(float)] )] }
        Note that momentum and energy are a list of arrays, where each array correspons to a param point
    
    :param path: the direction to the file containing all the events information.
    :type path: string
    """

    def __init__(self, path, particle_ids=None, var_of_interest=None, outgoing=False, verbose=0):
        self.path = path
        self.verbose = verbose
        self.particle_ids = particle_ids
        self.var_of_interest = var_of_interest
        self.outgoing = outgoing

        self.scan = {}

    def fill_up_scan(self):
        file_list = glob.glob(self.path)

        # initialize scan dictionary as described above going over all the names
        for i, name in enumerate(file_list):
            read_file = ReadLhe(name, var_of_interest=self.var_of_interest, particle_ids=self.particle_ids, 
                                outgoing=self.outgoing, verbose=self.verbose)
            params = read_file.extract_params_from_path()
            if i == 0:            
                for param, value in params.items():
                    self.scan[param] = []
                self.scan["file_name"] = [] 
                for variable in list(read_file.data.columns):
                    self.scan[variable] = []

            for param, value in params.items():
                self.scan[param].append(value)
            self.scan["file_name"].append(os.path.basename(name))
            for variable in list(read_file.data.columns):
                self.scan[variable].append(read_file.data[variable].to_numpy())

    def save_scan(self, save_name="complete_task.pickle"):
        with open(save_name, "wb") as f:
            pickle.dump(self.scan, f, protocol=pickle.HIGHEST_PROTOCOL)
