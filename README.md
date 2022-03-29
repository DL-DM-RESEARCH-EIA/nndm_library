# nndm_library
A library design to apply the machine learning and traditional methods to the discovery of new particles.

### Installation
```
pip install nndm_library
```

### Get started
Reading of lhe or root files for signal and background and applying some methods to identify the signal.

```Python
from nndm_library import ReadLhe
from nndm_library import ReadRoot

# read lhe
read_file = ReadLhe("tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe", verbose=0)

# read root file
file = ReadRoot("tests/data_tests/background/v_e_scattering/onantinuelepton10125.root")
print(file.df)
```
