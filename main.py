from numpy import save
from nndm_library import ReadFileBase
from nndm_library import ReadLhe
from nndm_library import FilesManipulator

# Returns the parameters in the name
print("Read parameters")
print("*" * 50)
read_file = ReadFileBase("eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe")
type_particle, mk, eps2 = read_file.extract_params_from_path()
print(mk, type(mk))


# Return the quadri-momenta values of the file in question 
print("\n" * 2)
print("*" * 50)
print("Read momenta")
read_file = ReadLhe("tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe", verbose=0)
print(read_file.path)
print(read_file.get())


print("\n" * 2)
print("*" * 50)
print("Read momenta of set of files")
# Scan data
file_manipulator = FilesManipulator("tests/data_tests/signal/*.lhe")
file_manipulator.fill_up_scan()
print(file_manipulator.scan)


print("\n" * 2)
print("*" * 50)
print("Save scan of files")
from os.path import exists
save_file = "tests/data_tests/signal/signal.pickle"
file_manipulator.save_scan(save_file)
print("file " + save_file + " exists?", "yes" if exists(save_file) else "no")

print("\n" * 2)
print("*" * 50)
print("Read Root file")