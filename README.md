**Contact DEM Evaluator** is a high-performance OOP C++ tool designed for analyzing contact time and force distributions in Discrete Element Method (DEM) simulations, specifically optimized for LIGGGHTS output.

# Build:
For build you can use attached CMakeList.txt.
Or you can use custom build.
```
# For Windows
g++ main.cpp -o contactDem.exe
```
```
# For Linux
g++ main.cpp -o contactDem
```
Note. that it requires g++ compiler. Also note that for intel architecture is better to use icpx build.
# LIGGGHTS implementation:
For LIGGGHTS input script its necessary to include following commands.

**1) Compute force-chains:**
```
compute         forcechainspair all pair/gran/local id pos vel force
compute         forcechainswall all wall/gran/local id pos vel force
```

**2) For dumps:**
```
dump    dumpoutforcepair all local ${dumpfreq} post/pair_forcechains_*.dump c_forcechainspair[1] c_forcechainspair[2] c_forcechainspair[3] c_forcechainspair[4] c_forcechainspair[5] c_forcechainspair[6] c_forcechainspair[7] c_forcechainspair[8] c_forcechainspair[9] c_forcechainspair[10] c_forcechainspair[11] c_forcechainspair[12] c_forcechainspair[13] c_forcechainspair[14] c_forcechainspair[15] c_forcechainspair[16] c_forcechainspair[17] c_forcechainspair[18]
dump    dumpoutforcewall all local ${dumpfreq} post/wall_forcechains_*.dump c_forcechainswall[1] c_forcechainswall[2] c_forcechainswall[3] c_forcechainswall[4] c_forcechainswall[5] c_forcechainswall[6] c_forcechainswall[7] c_forcechainswall[8] c_forcechainswall[9] c_forcechainswall[10] c_forcechainswall[11] c_forcechainswall[12] c_forcechainswall[13] c_forcechainswall[14] c_forcechainswall[15] c_forcechainswall[16] c_forcechainswall[17] c_forcechainswall[18]  
```
# Usage:
