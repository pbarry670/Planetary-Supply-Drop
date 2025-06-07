#include "edl.h"
#include "orbits.h"
#include <random>
#include <iostream>
using namespace std;

ActualSatellite::ActualSatellite(float m, float CD, float A, Eigen::Matrix<float,6,1> x_ECI) :
m(m), CD(CD), A(A), x_ECI(x_ECI) { 
    u_ECI << 0,0,0; 
    // Parameters not set after constructor: x_ECEF, x_LVLH, x_LLA, u_LVLH
}

ReferenceSatellite::ReferenceSatellite(float m, Eigen::Matrix<float,6,1> x_ECI) :
m(m), x_ECI(x_ECI) {
    // Parameters not set after constructor: x_ECEF, x_LVLH, x_LLA
}

Earth::Earth(Eigen::Matrix<float,6,1> x) : x(x) {
    r_S2E(0) = x(0);
    r_S2E(1) = x(1);
    r_S2E(2) = x(2);
}

OrbitParams::OrbitParams(float alpha0, Eigen::Matrix<float,3,6> K, Eigen::Vector3d ddl, float DTheta, float DThetaTolerance, float ttdffp, float latTol, float lonTol) :
alpha0(alpha0), K(K), desired_drop_location(ddl), DTheta(DTheta), DThetaTolerance(DThetaTolerance), timeToDropFromFinalPass(ttdffp), latTolerance(latTol), lonTolerance(lonTol) {
    alpha = alpha0;

    Eigen::Matrix<float,3,3> R_LVLH_2_ECI;
    Eigen::Matrix<float,3,3> R_ECI_2_LVLH;
    Eigen::Matrix<float,3,3> R_ECEF_2_ECI;
    Eigen::Matrix<float,3,3> R_ECI_2_ECEF;

    R_LVLH_2_ECI << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;

    R_ECI_2_LVLH << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;

    R_ECEF_2_ECI << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;

    R_ECI_2_ECEF << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;

    // Parameters not set after constructor: lonAtDesiredLatFirst; lonPhaseShiftPerOrbit, timeOfFirstPass, timeOfSecondPass, timeOfFinalPass

}


Eigen::Matrix<float,6,1> satelliteActualOrbitDynamics(ActualSatellite sat, Eigen::Matrix<float,6,1> x, Eigen::Vector3d u, Eigen::Vector3d r_S2E){

    float r = sqrt(pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2));
    float V = sqrt(pow(x(3), 2) + pow(x(4), 2) + pow(x(5), 2));

    Eigen::Vector3d J2AccelVec(  (1 - 5*pow((x(2)/r),2))*(x(0)/r),  (1 - 5*pow((x(2)/r),2))*(x(1)/r), (3 - 5*pow((x(2)/r),2))*(x(2)/r) );

    Eigen::Vector3d f_J2;
    f_J2(0) = (-3/2)*J2*(MU_E/pow(r,2))*(pow(R_EARTH/r, 2))*J2AccelVec(0);
    f_J2(1) = (-3/2)*J2*(MU_E/pow(r,2))*(pow(R_EARTH/r, 2))*J2AccelVec(1);
    f_J2(2) = (-3/2)*J2*(MU_E/pow(r,2))*(pow(R_EARTH/r, 2))*J2AccelVec(2);

    float rho = pow(SEA_LEVEL_RHO, (-1*(r-R_EARTH)/SCALE_HEIGHT));

    Eigen::Vector3d f_drag;
    f_drag(0) = -0.5*rho*((sat.CD*sat.A)/sat.m)*V*x(3);
    f_drag(1) = -0.5*rho*((sat.CD*sat.A)/sat.m)*V*x(4);
    f_drag(2) = -0.5*rho*((sat.CD*sat.A)/sat.m)*V*x(5);

    float cosTheta = (r_S2E(0)*x(0) + r_S2E(1)*x(1) + r_S2E(2)*x(2))/(r*sqrt(pow(r_S2E(0), 2) + pow(r_S2E(1), 2) + pow(r_S2E(2), 2)));
    float theta = acos(cosTheta);
    float v; // light indicator, =0 if in darkness
    if (theta < PI/2 - acos(R_EARTH/r)) {
        v = 0;
    } else {
        v = 1;
    }

    Eigen::Vector3d r_S2Sat( r_S2E(0)+x(0), r_S2E(1)+x(1), r_S2E(2)+x(2) );
    Eigen::Vector3d unit_r_S2Sat = r_S2Sat/r_S2Sat.norm();

    Eigen::Vector3d f_SRP;
    f_SRP(0) = -P_SMF*((v*sat.A)/sat.m)*CR*unit_r_S2Sat(0);
    f_SRP(1) = -P_SMF*((v*sat.A)/sat.m)*CR*unit_r_S2Sat(1);
    f_SRP(2) = -P_SMF*((v*sat.A)/sat.m)*CR*unit_r_S2Sat(2);

    Eigen::Vector3d f_grav;
    f_grav(0) = (-MU_E/pow(r,3))*x(0);
    f_grav(1) = (-MU_E/pow(r,3))*x(1);
    f_grav(2) = (-MU_E/pow(r,3))*x(2);

    Eigen::Matrix<float,6,1> xdot;
    xdot(0) = x(3);
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = f_grav(0) + f_J2(0) + f_drag(0) + f_SRP(0) + u(0)/sat.m;
    xdot(4) = f_grav(1) + f_J2(1) + f_drag(1) + f_SRP(1) + u(1)/sat.m;
    xdot(5) = f_grav(2) + f_J2(2) + f_drag(2) + f_SRP(2) + u(2)/sat.m;

    return xdot;

}

Eigen::Matrix<float,6,1> satelliteReferenceOrbitDynamics(ReferenceSatellite sat, Eigen::Matrix<float,6,1> x) {

    float r = sqrt(pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2));

    Eigen::Vector3d f_grav;
    f_grav(0) = (-MU_E/pow(r,3))*x(0);
    f_grav(1) = (-MU_E/pow(r,3))*x(1);
    f_grav(2) = (-MU_E/pow(r,3))*x(2);

    Eigen::Matrix<float,6,1> xdot;
    xdot(0) = x(3);
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = f_grav(0);
    xdot(4) = f_grav(1);
    xdot(5) = f_grav(2);

    return xdot;

}

Eigen::Matrix<float,6,1> earthOrbitDynamics(Eigen::Matrix<float,6,1> x) {

    float r = sqrt(pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2));

    Eigen::Vector3d f_grav;
    f_grav(0) = (-MU_S/pow(r,3))*x(0);
    f_grav(1) = (-MU_S/pow(r,3))*x(1);
    f_grav(2) = (-MU_S/pow(r,3))*x(2);

    Eigen::Matrix<float,6,1> xdot;
    xdot(0) = x(3);
    xdot(1) = x(4);
    xdot(2) = x(5);
    xdot(3) = f_grav(0);
    xdot(4) = f_grav(1);
    xdot(5) = f_grav(2);

    return xdot;

}
