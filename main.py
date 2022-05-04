import numpy as np
from nndm_library import ReadFileBase
from nndm_library import ReadLhe
from nndm_library import FilesManipulator
from nndm_library import ReadRoot
from nndm_library import Constants

# Read standard basic file
print("Read single basic file")
print("*" * 50)
file = ReadFileBase("tests/data_tests/basic_type/some_particle_names/quarks.txt")
print(file.data)

print("Read basic files in a directory")
print("*" * 50)
file = ReadFileBase("tests/data_tests/basic_type/some_particle_names/")
print(file.data)

print("Read basic files recursively")
print("*" * 50)
file = ReadFileBase("tests/data_tests/basic_type/", recursive=True, relabel_events=True)
print(file.data)

# # Returns the parameters in the name
# print("Read parameters")
# print("*" * 50)
# read_file = ReadFileBase("eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe")
# res_dict = read_file.extract_params_from_path()
# print(res_dict)

# # Return the quadri-momenta values of the file in question 
# print("\n" * 2)
# print("*" * 50)
# print("Read momenta")
# read_file = ReadLhe("tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe", verbose=0)
# print(read_file.path)
# print(read_file.get_df())

# # Reading files from signal
# print("\n" * 2)
# print("*" * 50)
# print("Read momenta of set of files")
# # Scan data
# file_manipulator = FilesManipulator("tests/data_tests/signal/*.lhe")
# file_manipulator.fill_up_scan()
# print(file_manipulator.scan)

# print("\n" * 2)
# print("*" * 50)
# print("Save scan of files")
# from os.path import exists
# save_file = "tests/data_tests/signal/signal.pickle"
# file_manipulator.save_scan(save_file)
# print("file " + save_file + " exists?", "yes" if exists(save_file) else "no")

# # Testing reading of root files
# print("\n" * 2)
# print("*" * 50)
# print("Read Root file")
# file = ReadRoot("tests/data_tests/background/v_e_scattering/onantinuelepton10125.root")
# print(file.get_df())

# print("\n" * 2)
# print("*" * 50)
# print("Read background from all the roots in a given directory")
# file = ReadRoot("tests/data_tests/background/v_e_scattering/")
# print(file.get_df())
# print(file.files_dir)

# last_index = file.get_df().index.tolist()[-1]
# if type(last_index) != int:
#     last_index = last_index[0]

# print( (last_index + 1) * 2, len(file.get_df()))

# print("\n" * 2)
# print("*" * 50)
# print("Read background all the root found walking through all the directories in the directory specified")
# file = ReadRoot("tests/data_tests/background/", recursive=True)
# print(file)
# print(file.get_df())

# last_index = file.get_df().index.tolist()[-1]
# if type(last_index) != int:
#     last_index = last_index[0]

# print( (last_index + 1) * 2, len(file.get_df()))

# Testing reading of lhe files
# filename = "tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.8832099149961855e-09.lhe"

# file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], var_of_interest=['px'], outgoing=True)
# print(file.get_df())

# file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], outgoing=True)
# print(file.get_df())

# file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID])
# print(file.get_df())
# print(np.unique(file.get_df()["status"]))

# file = ReadLhe(path=filename)
# print(file.get_df())
# print(np.unique(file.get_df()["id"]))