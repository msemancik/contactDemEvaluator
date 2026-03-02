**Contact DEM Evaluator** is a high-performance OOP C++ tool designed for analyzing contact time and force distributions in Discrete Element Method (DEM) simulations, specifically optimized for LIGGGHTS output.

# Build:
For build you can use attached CMakeList.txt.
```bash
sudo apt update
cd ~
git clone https://github.com/msemancik/contactDemEvaluator.git
cd contactDemEvaluator
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
make install -j
```
Or optionally you can use make.sh file for common C++ file compilation.
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
After you run a simulation just open terminal and type:
```
./contactDem -dry folderPath initialStep dumpStep endStep
```
Example might look like this:
```
cd /home/user/Documents/mySimulation/
./contactDem -dry post 0 4000 12000
```
# Potential Troubleshooting
**1) Finding initial dump:**
To set initial parameters for calculation it is necessary to load initial dump.
From attached code it is clear that initial dump has to be in this format:
a) file starts with dump prefix, b) file has in its name rev_, c) file ends with .atom extension. For your environment it is necessary to change this part of code.
```cpp
bool findInitDump(const std::string& name) {
    bool starts = (name.size() >= 5 && name.compare(0, 5, "dump_") == 0);

    bool ends = (name.size() >= 5 && name.compare(name.size() - 5, 5, ".atom") == 0);

    return starts && (name.find("rev_") != std::string::npos) && ends;
}
```
**2) Timestep determination:**
Unfortunately from initial dump there is no possibility to load timestep info. For default settings timestep is determinated via radius. For your different timestep you have to change following setter.
```cpp
void setDeltaT(double radius) {
  if (radius < 5e-4) {
    deltaT = 1e-5;
          }
  else {
    deltaT = 2e-5;
          }
}
```
