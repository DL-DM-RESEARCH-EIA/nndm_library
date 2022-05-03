# nndm_library
A library design to apply the machine learning and traditional methods to the discovery of new particles.

### Installation
```
pip install nndm_library
```

### Get started
First, for the following examples to work, the files in the directory tests/data_tests/ must be downloaded. In short

Reading lhe files

we will read all data associated with a file in a given data frame by default. In addition to data, there are several
options available to filtrate what is intender to be readed exactly. 

```Python
from nndm_library import ReadLhe

filename = "tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.8832099149961855e-09.lhe"

# Read only the px data associated with the electrons that are going out.
file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], var_of_interest=['px'], outgoing=True)
print(file.data)

# Read all data associated with the electrons that are going out.
file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], outgoing=True)
print(file.data)

# Read all data associated with electrons.
file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID])
print(file.data)
print(np.unique(file.data["status"]))

# Read all data in given lhe file
file = ReadLhe(path=filename)
print(file.data)
print(np.unique(file.data["id"]))
```

Reading root files

Again, we will read all the data in the files and save it in a dataframe by default.

```Python
from nndm_library import ReadRoot

# The root file will be readed to the leafs indicated. The pattern_output first indicates that if the a same root file
#   has names like treeout;2, treeout;3, treeout;4, treeout;5,etc. Then, the one with the lest nomber will be chosen.
#   In the and the complete path to the leafst is like this treeout;{lest_number}/output_base_middle_branch/leaf
file = ReadRoot("tests/data_tests/background/v_e_scattering/onantinuelepton10125.root",
                output_base_tree="treeout", pattern_output="first", output_base_middle_branch = "/e/out",
                leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"])
print(file.df)

# The reading can be done to all the .root files that are in a given. In the below example they will be
#   somehting like tests/data_tests/background/v_e_scattering/file1.root,
#                  tests/data_tests/background/v_e_scattering/file2.root.

file = ReadRoot("tests/data_tests/background/v_e_scattering/", output_base_tree="treeout", 
                pattern_output="first", output_base_middle_branch = "/e/out",
                leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"])
print(file.df)

# Finally, the reading can also be done to all the files in a given directory recursively. In the example below, it is
#   something like tests/data_tests/background/dir1/file1.root, tests/data_tests/background/dir2/file1.root,
#   tests/data_tests/background/file1.root
file = ReadRoot("tests/data_tests/background/", output_base_tree="treeout", 
                pattern_output="first", output_base_middle_branch = "/e/out",
                leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"], recursive=True)
print(file.df)
```

Othe utils. 

```Python
# Extract param names and values, along with particl type (check Constants class for the last)

read_file = ReadFileBase("eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe")
res_dict = read_file.extract_params_from_path()

# If you need to read different files and classify them according to their parameters.

file_manipulator = FilesManipulator("tests/data_tests/signal/*.lhe")
file_manipulator.fill_up_scan()
print(file_manipulator.self_scan())

# saving your scan results
from os.path import exists
save_file = "tests/data_tests/signal/signal.pickle"
file_manipulator.save_scan(save_file)
print("file " + save_file + " exists?", "yes" if exists(save_file) else "no")
```
