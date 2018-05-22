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

# Tuples in the test.py

https://github.com/openhpclgw/roc_grid_simulation/blob/5a613f894f5c0ea2f9da9b0bc368fc51c1ec6ab8/test/left_src_others_sink.py#L10

This and the following line sets the locations of sinks and sources.
They are 4-tuples, that define the bounding box for a boundary
condition. The tuple is in form (left_pos, top_pos, width, height).

Every source and sink must be a rectangle and must fit into a mesh node
when discretized. if you want to do complex shapes, as the sink in this
example you can pass a list of tuples.

# Environment Setup after Cloning
```
git clone https://github.com/openhpclgw/roc_grid_simulation.git
cd roc_grid_simulation/
unset PYTHONPATH
source setenv.sh
cd test/
python3 test.py 20
sudo apt-get install python3-numpy
python3 test.py 20
ls
mkdir tmp
python3 test.py 20
sudo apt-get install ngspice
python3 test.py 20
sudo apt-get install python3-matplotlib
python3 test.py 20
```

# Questions about the problem Setup

```
#Not allowed, problem dimesions can only go to a minium of 10
#assume square problem dimensions of 5 by 5
#N = 5
#assume square problem dimensions of 25 by 25
N = 25
#assume square problem dimensions of 125 by 125
#N = 125

#Grid Layout

#0
#1
#2
#3
#.
#.
#.
#   0   1   2   ...

# The tuple is in form (left_pos, top_pos, width, height)

source = (0, 0, 25, 1)
sink = (0, 5, 1, 5)

# sink = [(11, 2, 1, 10), (5, 15, 2, 2)]
# source = (4, 4, 4, 12)

cond_exp = -3
conductance = 10**cond_exp

#mesh dimension entered by user
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.0)

```
