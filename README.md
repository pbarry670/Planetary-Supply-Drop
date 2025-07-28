# Planetary-Supply-Drop
This repository serves as a toolset for the concept of a "planetary supply drop," or the notion that some payload in orbit can be delivered to any specified location on Earth. This idea is advantageous because the time it takes an object to reach sea level from a 500 km altitude orbit around Earth is typically between 15-25 minutes, depending on the entry, descent, and landing system employed and its associated parameters. For such an idea to be brought to fruition, it is necessary to draw from a body of knowledge consisting of launch vehicle dynamics, orbital mechanics, control theory, ballistic entry physics, reference frames, and trajectory optimization - as well as implementing  all of these big ideas in code! 

The repository is divided into a handful of folders. First, the MATLAB folder contains some algorithm development code used to quickly prototype ideas and implementations. This is done before translation into C++, the folder for which contains the primary executable, main.exe, as well as header and source files for all associated C++ code. The data collected during simulation is written into .csv files, which is why the Python folder exists: Python's vpython library is leveraged to create animations to link a visual representation to the math being carried out. The Python folder also contains code that carries out the Guidance for Fuel-Optimal Large Diverts (G-FOLD) algorithm, which requires solving a convex optimization problem. This, however, is called by and returns data to the C++ code, keeping all execution within one file: main.exe.
<p align="center">
<img width="361" height="308" alt="image" src="https://github.com/user-attachments/assets/4b243c4c-d25e-4388-a8a1-23cfb6da29f3" />
<img width="361" height="308" alt="image" src="https://github.com/user-attachments/assets/0583c675-a954-4b27-81af-c7b07a659cbc" />
</p>


<p align="center">
<img width="900" height="299" alt="image" src="https://github.com/user-attachments/assets/9cf76346-92cb-46cc-a360-fb90508b46c9" />
</p>
<p align="center">
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/34679b89-ab01-4c02-911d-661399765ad5" />
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/199407d7-e83e-4c72-8769-706742f01a06" />
</p>


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
```
g++ -std=c++11 -I path/to/eigen-X.X.X -I"C:\PythonXX\include" -L"C:\PythonXX\libs" src/edl.cpp src/orbits.cpp src/main.cpp -o main -lpythonXX
```
In  the above, `eigen-X.X.X` refers to the version of Eigen installed (anything 3.3 or newer should be fine), and PythonXX refers to the Python version installed. I have Python 3.7, so this would be Python37 for me. Additionally, I used the default install location for the Python path to its C/C++ header files, so if your install is different, this path may be different than simply `C:\PythonXX`.

This generates an executable, main.exe, within the C++ folder. To run the .exe file, simply type in the terminal `main.exe`. The code contains several helpful print statements aimed at helping the user understand what is happening during simulation.

The executable loads data into an assortment of .csv files; these .csv files may be read into Python for visualization with vpython. To run the visualization, ensure the .csv files have been produced, navigate to the Python folder, and run the command `python masterViz.py`. This will open up a browser window that displays all visualizations in sequence.

## Simulation Summary
Below is a summary of the flow of simulation from beginning to end, including an explanation of the algorithms at play in this code package. Feel free to take a look at the code and if you have any suggestions or ideas for additions, email pbarry670@gmail.com.

### Precursors
The idea behind this simulation is that it is an end-to-end planetary supply drop. While it is not implemented yet, the idea is that a launch vehicle takes some spacecraft to orbit. This spacecraft contains its own propulsion system, to perform orbital correction and trajectory adjustment maneuvers, but primarily contains a capsule that carries the payload for the supply drop. The capsule is deployed from the satellite for atmospheric entry; at a specified altitude, the heatshield jettisons and a parachute deploys, slowing the capsule more, and at an even lower specified altitude, the capsule separates from a lander module that contains thrusters to guide it softly to the desired drop location. All dynamical propagation is performed with 4th-order Runge-Kutta, except for powered descent, which performs the trajectory optimization with Forward Euler. An exponential atmospheric model is used, and gravity is modeled as decreasing with altitude.

###  ⚠️ Under Construction ⚠️  Launch To Orbit

### Orbit Propagation
The initial state for the spacecraft orbit is a circular, equatorial LEO orbit. The spacecraft is parked in this orbit for some number of orbital periods before changing its inclination to the latitude of the desired drop point. While in orbit, there are a couple things happening behind the scenes. First, the spacecraft is being subjected to several forces other than just Earth's gravity. These include:
- non-spherical gravity (J2 effect)
- atmospheric drag
- solar radiation pressure
#### Orbit Control and Trajectory Adjustment
If the spacecraft experienced all of these perturbing forces with no corrective maneuvers, its trajectory would slowly deviate from the desired trajectory. As such, an orbit controller is utilized to bring the actual trajectory into alignment with the reference trajectory. This orbit controller was derived using the Clohessy-Wiltshire (CW) equations for relative orbital motion:
<p align="center">
<img src="https://tse3.mm.bing.net/th/id/OIP.S5tzGvMvLkwzoutPf3v2qAHaEg?r=0&rs=1&pid=ImgDetMain&o=7&rm=3" alt="Description of image" width="600" />
</p>
In this formulation, the reference trajectory is said to be the "chief," with the actual satellite under perturbations as the "deputy." Because the relative orbit equations of motion are linear, the controller selected to perform orbit control is Linear Quadratic Regulator (LQR). The CW dynamics are placed in the dynamics (state) matrix, while the controller is framed as forces which affect the vehicle acceleration in the CW frame. A 6x6 Q matrix can be formed to denote the emphasis on controller performance while a 3x3 R matrix influences the level of control effort allowable by the system. Using MATLAB's lqr(), a gain matrix K can be selected to find the control inputs according to u = -Kx, where x is the state of the actual satellite relative to the reference trajectory in the CW frame. All navigation and control in-orbit is performed in the Earth-centered inertial (ECI) frame, so the control inputs are rotated from the CW frame into the ECI frame. In a practical application, for spacecraft thruster control, these could be further rotated into the vehicle's body frame from the ECI frame.
<br>
<br>
In order for the spacecraft to drop the capsule so it lands on desired drop point, the spacecraft must have its orbit pass over the drop point. This is simple enough; starting from an equatorial orbit, the spacecraft must perform an inclination change equal to the latitude of the desired drop point. Put more generally, the desired inclination of the spacecraft's orbit is equal to the latitude of the desired drop point. Subsequently, the spacecraft must execute a Delta-V manuever to change its inclination. Because the spacecraft occupies an equatorial orbit, the Delta-V could be equal to 2*sin(i/2), where i is the inclination change. However, to keep this as general as possible, the Delta-V magnitude and direction is applied in several steps: <br>
1. Save the current velocity of the spacecraft. <br>
2. Obtain the orbital elements of the spacecraft's orbit. <br>
3. Find the velocity of the spacecraft if it were to contain the same orbital elements, only changing the desired inclination to the drop point's latitude. <br>
4. To find the Delta-V, subtract the current velocity from the desired velocity. <br>
<p align="center">
<img width="455" height="194" alt="image" src="https://github.com/user-attachments/assets/41abff8b-bb81-4a4a-8eb5-5715d82ac730" />
</p>

#### ⚠️ Under Construction ⚠️ Orbit Determination
#### Capsule Drop Point Calculation
With the spacecraft in the desired orbit so that it passes over the drop point, it is now relatively straightforward to compute when the capsule should be dropped. Ahead of orbital simulation, the ballistic descent parameters are simulated to determine the total angle over Earth that is covered from the time the capsule drops to the time ballistic descent ends (at which point it is assumed the vehicle is traveling with the rotation of the Earth, an approximation). Knowing this total angle over Earth covered by ballistic descent, $\theta$, informs us that the spacecraft should travel over the angle (360&#176; - $\theta$) after passing over the drop point latitude for the final time. How do we know when it is the final time the spacecraft passes over the drop point latitude? I use the following logic:
- The first time the spacecraft passes over the drop point latitude, the spacecraft longitude is saved.
- The second time the spacecraft passes over the drop point latitude, the spacecraft longitude is again saved.
- By knowing both of these longitude values, one can determine the longitudinal phase shift of the point of Earth at the drop point latitude that the spacecraft passes over.
  - Ex. Let's say the drop point is located at 30&#176; latitude. The first time the spacecraft passes over this latitude, its longitude is 5&#176;. The second time, its longitude is 10&#176;. If the drop point is located at 20&#176; longitude, we know the capsule should deploy after the third orbit but before completing the fourth orbit so it intercepts the drop point upon falling.
- Knowing the longitudinal phase shift, mark the capsule as being ready to deploy when the spacecraft's longitude at the drop point latitude is equal to the negative of the longitudinal phase shift.
- Deploy the capsule after the spacecraft has traveled (360&#176; - $\theta$) through its orbital plane from the drop point latitude, where $\theta$ is the angle over Earth spanned by ballistic descent. 
<p align="center">
<img width="263" height="275" alt="image" src="https://github.com/user-attachments/assets/e4997999-105b-400c-9784-e2e6516c8165" />
</p>

### Entry, Descent, and Landing
After the capsule is dropped from the satellite, it follows an entry, descent, and landing sequence to reach Earth's surface. The architecture followed for this is a three-stage EDL sequence: <br>
1. Perform an uncontrolled ballistic descent to the chute deploy height.
2. At the chute deploy height, jettison the capsule's heatshield and deploy the parachute, entering chute descent.
3. At the powered descent height, cut the chute lines and have the lander module activate its thrusters to guide itself to the drop point, landing with zero velocity.
#### Ballistic Descent
Ballistic descent is performed by propagating three states: speed, flight path angle, and altitude. The flight path angle is defined as the angle between the local horizontal and the velocity vector. At handoff from ballistic descent, the altitude is set as the orbit altitude, the speed is set as the norm of the spacecraft's velocity at the capsule deployment point, and the flight path angle is initialized at 5&#176;. Note that the ballistic descent is considered to be planar, i.e., all motion occurs in the orbital plane from which the capsule was deployed. Ballistic dynamics are propagated from orbit altitude until the altitude is equal to the chute deploy height. Over the course of ballistic descent, the flight path angle slowly increases to 90&#176; due to atmospheric drag effects, leading the capsule to be falling straight down.
<p align="center">
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/ce62680b-a4b9-40e6-953b-df0fbdc16fbb" />
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/0a2e4b21-b2ff-496b-9aa9-65940f3ed24c" />
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/b4324dc2-22dd-4741-b127-1360734d9882" />
</p>
#### Chute Descent
At the beginning of chute descent, the mass of the capsule has the mass of the heatshield subtracted from it. A simple drag model is utilized for the parachute, where the upward force provided is proportional to air density, drag coefficient, and speed squared. Wind and other perturbations are modeled as perturbing cross-range accelerations, represented as random variables centered around some nonzero mean. From the capsule deployment calculation in the orbit, the ballistic descent is assumed to end with the capsule directly above the drop point. These perturbing accelerations introduce cross-range error, resulting in the need for powered descent to correct these errors to land exactly at the drop point.
<p align="center">
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/4eb7774b-9df2-494a-95f6-4eb196f348f3" />
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/092c8357-9317-428e-85e2-1b600e6c2483" />
</p>

#### Powered Descent
The powered descent phase carries the lander from some starting position and velocity to a zero-velocity state with minimum positional landing error and minimal fuel usage. This is performed through a convexification of the planetary soft landing problem, which is nicknamed as Guidance for Fuel-Optimal Large Diverts (G-FOLD). G-FOLD takes into account many parameters of the landing problem:
- Lander mass and fuel usage constraints
- Maximum and minimum thrust constraints
- Thrust pointing direction with constraint
- Lander ground approach angle constraints
- Maximum velocity constraints
- Mass and translational dynamics <br>
<br>
Attitude dynamics are ignored because it is assumed that attitude corrections may be performed with much smaller thrust values at a much higher frequency than trajectory corrections. Aerodynamic effects are also ignored, as gravity and thrust are the only modeled forces on the lander; this assumption is nominally valid when the lander is traveling at low speeds.
<br>
<br>
My implementation of the G-FOLD algorithm is based off of a paper by Acikmese, Carson, and Blackmore, titled "Lossless Convexification of Nonconvex Control Bound and Pointing Constraints of the Soft Landing Optimal Control Problem." The planetary soft landing problem constraints are first placed in the form of convex constraints, and the trajectory optimization is framed as two subsequent convex optimization problems. First, the landing error is minimized (given dynamics, mass, thrust, glideslope, and boundary condition constraints); the landing error is then saved as a maximum bound on landing error in a subsequent problem that minimizes the fuel usage to reach a final condition with that landing error or better. Because of the complexities in implementing a convex optimization problem in C++, the C++ code calls Python code that utilizes the convex optimization library "cvxpy" to carry out the G-FOLD algorithm. The solved trajectory with mass and thrust control history is then returned to C++ for data logging. In the event of the optimization failing, the user is prompted to relax their powered descent constraints so the solver can converge to a feasible solution. If the G-FOLD solver does find a solution, the final velocity of the vehicle is zero, with minimum landing error relative to the drop point and minimum fuel usage from the lander to reach that point. The payload aboard the lander has been delivered to the drop point, all the way from launch and orbit. <br>
<p align="center">
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/52be09ad-4046-43fc-8c5a-eecb547e7932" />
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/e02ba098-1141-4bef-b526-e5f3b0dbcf49" />
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/3c874d2d-e3c0-4308-82b6-c5faf1ee6537" />
</p>














