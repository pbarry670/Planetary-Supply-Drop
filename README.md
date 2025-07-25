# Planetary-Supply-Drop
This repository serves as a toolset for the concept of a "planetary supply drop," or the notion that some payload in orbit can be delivered to any specified location on Earth. This idea is advantageous because the time it takes an object to reach sea level from a 500 km altitude orbit around Earth is typically between 15-25 minutes, depending on the entry, descent, and landing system employed and associated parameters. For such an idea to be brought to fruition, it is necessary to draw from a body of knowledge consisting of launch vehicle dynamics, orbital mechanics, control theory, ballistic entry physics, reference frames, and trajectory optimization - as well as implementing  all of these big ideas in code! 
The repository is divided into a handful of folders. First, the MATLAB folder contains some algorithm development code used to quickly prototype ideas and implementations. This is done before translation into C++, the folder for which contains the primary executable, main.exe, as well as header and source files for all the associated code. The data collected during simulation is written into .csv files, which is why the Python folder exists: Python's vpython library is leveraged to create animations to link a visual representation to the math being carried out. The Python folder also contains code that carries out the Guidance for Fuel-Optimal Large Diverts (G-FOLD) algorithm, which requires solving a convex optimization problem. This, however, is called by and returns data to the C++ code, keeping all execution as simple as running main.exe.
## Running The Code
In order to run the code, there are a handful of requirements:
- MATLAB
  - CVX
- C++
  - g++ compiler
  - Eigen
- Python
  - numpy
  - cvxpy
  - pandas
  - vpython
  - matplotlib
Then, in a terminal, change your current directory to the C++ folder. Run the command 
'''
g++ -std=c++11 -I path/to/eigen-X.X.X -I"C:\PythonXX\include" -L"C:\PythonXX\libs" src/edl.cpp src/orbits.cpp src/main.cpp -o main -lpythonXX
'''
In  the above, eigen-X.X.X refers to the version of Eigen installed (anything 3.3 or newer should be fine), and PythonXX refers to the Python version installed. I have Python 3.7, so this would be Python37 for me. Additionally, I used the default install location for the Python path to its C/C++ header files, so if your install is different, this path may be different than simply 'C:\PythonXX'.
This generates a 
## Simulation Summary
Below is a summary of the flow of simulation from beginning to end, including an explanation of the algorithms at play in this code package. Feel free to take a look at the code and if you have any suggestions or ideas for additions, email pbarry670@gmail.com.
