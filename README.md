Insall python 

https://www.python.org/downloads/

If you have python3 as default python set in your environment do

# How to run

```
$ test.py 20
```

Otherwise if you have python3 installed but not as default you can run
as:

```
$ python3 test.py 20
```

where 20 is the mesh size.

# Files in Repo

## Library files

```
analysis_utils.py
heat_problem.py
roc_model.py
spice_gen.py
```

## Test files

```
current\_split.py
left\_src\_others\_sink.py
test.py
```

# Versions and Dependencies

I use the following software/libs

```
Linux 4.13.4
Python 3.6.2
  numpy
  matplotlib (Probably wont hit this import unless you plot something)

  (following should come with the python installation)
  itertools
  re
  sys
  subprocess
ngspice 27
```

# Using other SPICE simulators

If/when we try to use Windows, we need to move to a different SPICE
simulator with a command line interface. 

To change the simulator to be run `SpiceGenerator.run` needs to be
modified.

We may also need to change the way the output is parsed depending on
the simulator, which can be done in `SpiceGenerator.generate_results_dict`
