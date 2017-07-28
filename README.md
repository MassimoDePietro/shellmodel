# Description
This scientific code solves the equations of motion for the "SABRA" 
shell models (see https://doi.org/10.1103/PhysRevE.58.1811).  


# Dependencies
Cmake and a c++11 compiler.  


# How to compile
Use cmake. E.g. in linux:
```bash
mkdir build
cd build
cmake -g "Unix Makefiles" .. -DCMAKE_BUILD_TYPE=Debug/Release
make
```


# How to run a simulation
Run the executable with no arguments. It should write on screen all the 
parameters it needs for initializing the simulation, with a 
description. Some parameters (K0, Lambda, the coefficients A, B and C 
of the model) are hardcoded in the source file, so they cannot be 
passed as command line arguments.

The simulation, if not continuing from an old velocity field, 
randomizes the phases of the initial condition. You may want to control 
the seed used for the random number generator. One of the command line 
arguments will take care of this.

There are command line arguments for specifying the name of the output 
files. The output files will be created in the execution directory. 


# License
This code is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the 
Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

This code is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have 
received a copy of the GNU General Public License along with this code. 
If not, see http://www.gnu.org/licenses/.
