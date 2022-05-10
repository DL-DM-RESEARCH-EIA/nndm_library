# nndm_library
A library design to apply the machine learning and traditional methods to the discovery of new particles.

### Installation
```
pip install nndm_library
```

### Get started
First, for the following examples to work, the files in the directory tests/data_tests/ must be downloaded. In short:


Reading column files o the type
data1 data2 data3
v11   v12   v13
.     .     .
.     .     .
.     .     .
vn1   vn2   vn3

```Python
print("Read single basic file that contains two rows: particle id and particle name")
file = "tests/data_tests/basic_type/some_particle_names/quarks.txt"
print("Reading %s" % file)
print("*" * 50)
file = ReadFileBase("tests/data_tests/basic_type/some_particle_names/quarks.txt")
print(file.data)

print("\n\nRead all basic files like the one before, but in a given directory")
file = "tests/data_tests/basic_type/some_particle_names/"
print("Reading directory %s" % file)
print("*" * 50)
file = ReadFileBase("tests/data_tests/basic_type/some_particle_names/")
print(file.data)

print("\n\nRead all basic files like the one before, but all found recursively in a given directory")
file = "tests/data_tests/basic_type/"
print("Reading all files found recursively in directory %s" % file)
print("*" * 50)
file = ReadFileBase("tests/data_tests/basic_type/", recursive=True, relabel_events=True)
print(file.data)
```

Reading lhe files

we will read all data associated with a file in a given data frame by default. In addition to data, there are several
options available to filtrate what is intender to be readed exactly. 

```Python
# Single lhe file
filename = "tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.8832099149961855e-09.lhe"
print("Reading %s filetring by particle_ids, var_of_interest, and outgoing particles" % filename)
print("*" * 50)
file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], var_of_interest=['px'], outgoing=True)
print(file.data)

print("Reading %s filtering by particle_ids, and outgoing particles (all available variables are read" % filename)
print("*" * 50)
file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], outgoing=True)
print(file.data)

print("Reading %s filtering by particle_ids (all available variables for particles going in and out are read" % filename)
print("*" * 50)
file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID])
print(file.data)
print("in out out states")
print(np.unique(file.data["status"]))

print("Reading from %s all particles for all available variables for particles going in and out" % filename)
print("*" * 50)
file = ReadLhe(path=filename)
print(file.data)
print("particle ids")
print(np.unique(file.data["id"]))

# files in a directory
print("\n\nRead all files like the one before, but in a given directory")
filename = "tests/data_tests/signal/"
print("Reading directory %s" % filename)
print("*" * 50)
file = ReadLhe(path=filename)
print(file.data)
print(file.files_dir)
res_dict = file.extract_params_from_path()
print(res_dict)

# files in a directory recursively
print("\n\nRead all basic files like the one before, but all found recursively in a given directory")
filename = "tests/data_tests/signal/"
print("Reading all files found recursively in directory %s" % filename)
print("*" * 50)
file = ReadLhe(filename, recursive=True, relabel_events=True)
print(file.data)
print(file.files_dir)

res_dict = file.extract_params_from_path()
print(res_dict)
```

Reading root files

Again, we will read all the data in the files and save it in a dataframe by default.

```Python
from nndm_library import ReadRoot

# The root file will be readed to the leafs indicated. The pattern_output first indicates that if the a same root file
#   has names like treeout;2, treeout;3, treeout;4, treeout;5,etc. Then, the one with the lest nomber will be chosen.
#   In the and the complete path to the leafst is like this treeout;{lest_number}/output_base_middle_branch/leaf
print("\n" * 2)
print("*" * 50)
print("Read Root file")
file = ReadRoot("tests/data_tests/background/v_e_scattering/onantinuelepton10125.root",
                output_base_tree="treeout", pattern_output="first", output_base_middle_branch = "/e/out",
                leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"])
print(file.data)

# The reading can be done to all the .root files that are in a given. In the below example they will be
#   somehting like tests/data_tests/background/v_e_scattering/file1.root,
#                  tests/data_tests/background/v_e_scattering/file2.root.

print("\n" * 2)
print("*" * 50)
print("Read background from all the roots in a given directory")
file = ReadRoot("tests/data_tests/background/v_e_scattering/", output_base_tree="treeout", 
                pattern_output="first", output_base_middle_branch = "/e/out",
                leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"])
print(file.data)
print(file.files_dir)


# Finally, the reading can also be done to all the files in a given directory recursively. In the example below, it is
#   something like tests/data_tests/background/dir1/file1.root, tests/data_tests/background/dir2/file1.root,
#   tests/data_tests/background/file1.root
print("\n" * 2)
print("*" * 50)
print("Read background all the root found walking through all the directories in the directory specified")
file = ReadRoot("tests/data_tests/background/", output_base_tree="treeout", 
                pattern_output="first", output_base_middle_branch = "/e/out",
                leafs = ["out.t", "out.x", "out.y", "out.z", "out._mass"], recursive=True)
print(file.data)
```

Othe utils. 

```Python
# Extract param names and values, along with particl type (check Constants class for the last)

name = "tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe"
read_file = ReadLhe(name)
res_dict = read_file.extract_params_from_path()

print("\n" * 2)
print("*" * 50)
print("Read momenta of set of files")
# Scan data
file_manipulator = FilesManipulator("tests/data_tests/signal/*.lhe", var_of_interest=['e', 'px', 'py', 'pz'])
file_manipulator.fill_up_scan()
print(file_manipulator.scan)

print("\n" * 2)
print("*" * 50)
print("Save scan of files")
from os.path import exists
save_file = "tests/data_tests/signal/signal.pickle"
file_manipulator.save_scan(save_file)
print("file " + save_file + " exists?", "yes" if exists(save_file) else "no")
```
