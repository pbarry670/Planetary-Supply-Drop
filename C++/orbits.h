#ifndef ORBITS_H
#define ORBITS_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

const float R_EARTH = 6378000; // m
const float ORB_ALT = 500000; // m
const float PI = 3.14159;
const float W_E = (2*PI)/86164; // rad/s
const float MU_E = 3.986e14; // m^3/s^2
const float R_ORBIT = R_EARTH + ORB_ALT;
const float SPEED_ORBIT = sqrt(MU_E/R_ORBIT); // m/s
const float T_ORBIT =  2*PI*sqrt(pow(R_ORBIT,3)/MU_E); // s
const float ALPHA_0 = 0; // rad
const float ORBIT_DT = 1; // s





#endif