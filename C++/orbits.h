#ifndef ORBITS_H
#define ORBITS_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

const float R_EARTH = 6378000; // m
const float ORB_ALT = 500000; // m
const float PI = 3.14159;
const float W_E = (2*PI)/86164; // rad/s
const float MU_E = 3.986*pow(10,14); // m^3/s^2
const float MU_S = 1.327*pow(10,20); // m^3/s^2
const float R_ORBIT = R_EARTH + ORB_ALT;
const float SPEED_ORBIT = sqrt(MU_E/R_ORBIT); // m/s
const float T_ORBIT =  2*PI*sqrt(pow(R_ORBIT,3)/MU_E); // s
<<<<<<< HEAD
const float ORBIT_DT = 0.5; // s
=======
const float ORBIT_DT = 1; // s
>>>>>>> parent of ca89aaa... No runtime errors on orbit propagation, debugging math/algos now
const int NUM_ORBITS_BEFORE_INCL_CHANGE = 2;

const float J2 = 0.0010827; // J2 constant
const float P_SMF = 4.56*pow(10,-6); // solar momentum flux, N/m^2
const float CR = 1; // reflectivity coefficient

struct ActualSatellite {

    float m;
    float CD;
    float A;

    Eigen::Matrix<float,6,1> x_ECI;
    Eigen::Matrix<float,6,1> x_ECEF;
    Eigen::Matrix<float,6,1> x_LVLH;
    Eigen::Vector3d x_LLA;

    Eigen::Vector3d u_ECI;
    Eigen::Vector3d u_LVLH;

    ActualSatellite(float m, float CD, float A, Eigen::Matrix<float,6,1> x_ECI);

};

struct ReferenceSatellite {

    float m;
    Eigen::Matrix<float,6,1> x_ECI;

    ReferenceSatellite(float m, Eigen::Matrix<float,6,1> x_ECI);


};

struct Earth {

    Eigen::Matrix<float,6,1> x;
    Eigen::Vector3d r_S2E;

    Earth(Eigen::Matrix<float,6,1> x);

};

struct OrbitParams {

    float alpha0;
    float alpha;

    Eigen::Matrix<float,3,6> K;
    Eigen::Vector3d desiredDropLocation;

    Eigen::Matrix<float,3,3> R_LVLH_2_ECI;
    Eigen::Matrix<float,3,3> R_ECI_2_LVLH;
    Eigen::Matrix<float,3,3> R_ECEF_2_ECI;
    Eigen::Matrix<float,3,3> R_ECI_2_ECEF;

    bool hasOrbitedOnce = false;
    bool hasOrbitedTwice = false;
    bool readyToDrop = false;
    bool phaseShiftDetermined = false;

    float lonAtDesiredLatFirst;
    float lonAtDesiredLatSecond;
    float lonPhaseShiftPerOrbit;

    float latTolerance;
    float lonTolerance;

    float timeOfFirstPass;
    float timeOfSecondPass;
    float timeOfFinalPass;

    float DTheta;
    float DThetaTolerance;
    float timeToDropFromFinalPass;

    float timeOfDeployment;
    bool targetLocAtChuteDeployLogged = false;

    bool hasDeployed = false;

    OrbitParams(float alpha0, Eigen::Matrix<float,3,6> K, Eigen::Vector3d ddl, float DTheta, float DThetaTolerance, float ttdffp, float latTol, float lonTol);


};

Eigen::Matrix<float,6,1> satelliteActualOrbitDynamics(ActualSatellite sat, Eigen::Matrix<float,6,1> x, Eigen::Vector3d u, Eigen::Vector3d r_S2E);
Eigen::Matrix<float,6,1> satelliteReferenceOrbitDynamics(ReferenceSatellite sat, Eigen::Matrix<float,6,1> x);
Eigen::Matrix<float,6,1> earthOrbitDynamics(Eigen::Matrix<float,6,1> x);

Eigen::Matrix<float,6,1> propagateActualSatellite(ActualSatellite sat, Eigen::Matrix<float,6,1> x, Eigen::Vector3d u, Eigen::Vector3d r_S2E);
Eigen::Matrix<float,6,1> propagateReferenceSatellite(ReferenceSatellite sat, Eigen::Matrix<float,6,1> x);
Eigen::Matrix<float,6,1> propagateEarthState(Eigen::Matrix<float,6,1> x);

Eigen::Matrix<float,6,1> elements2RV(float a, float e, float i, float O, float w, float f);
Eigen::Matrix<float,5,1> RV2elements(Eigen::Matrix<float,6,1> x);

Eigen::Matrix<float,3,3> computeR_ECI_2_LVLH(float i, float O, float T);
Eigen::Matrix<float,3,3> computeR_ECEF_2_ECI(float alpha);
Eigen::Matrix<float,3,3> computeR_ECI_2_ECEF(float alpha);

Eigen::Vector3d ECEF2LLA(Eigen::Matrix<float,6,1> x_ECEF);
Eigen::Vector3d LLA2ECEF(Eigen::Vector3d x_LLA);

#endif