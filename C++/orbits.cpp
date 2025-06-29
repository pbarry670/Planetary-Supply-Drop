#include "edl.h"
#include "orbits.h"
#include <random>
#include <iostream>
using namespace std;

ActualSatellite::ActualSatellite(float m, float CD, float A, Eigen::Matrix<float,6,1> x_ECI) :
m(m), CD(CD), A(A), x_ECI(x_ECI) { 
    u_ECI << 0,0,0; 
    u_LVLH << 0,0,0;
    x_LVLH << 0,0,0,0,0,0;
    // Parameters not set after constructor: x_ECEF, x_LLA
}

ReferenceSatellite::ReferenceSatellite(float m, Eigen::Matrix<float,6,1> x_ECI) :
m(m), x_ECI(x_ECI) {}

Earth::Earth(Eigen::Matrix<float,6,1> x) : x(x) {
    r_S2E(0) = x(0);
    r_S2E(1) = x(1);
    r_S2E(2) = x(2);
}

OrbitParams::OrbitParams(float alpha0, Eigen::Matrix<float,3,6> K, Eigen::Vector3d ddl, float DTheta, float DThetaTolerance, float ttdffp, float latTol, float lonTol) :
alpha0(alpha0), K(K), desiredDropLocation(ddl), DTheta(DTheta), DThetaTolerance(DThetaTolerance), timeToDropFromFinalPass(ttdffp), latTolerance(latTol), lonTolerance(lonTol) {
    alpha = alpha0;

    Eigen::Matrix3d R_LVLH_2_ECI;
    Eigen::Matrix3d R_ECI_2_LVLH;
    Eigen::Matrix3d R_ECEF_2_ECI;
    Eigen::Matrix3d R_ECI_2_ECEF;

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

    // Parameters not set after constructor: lonAtDesiredLatFirst; lonPhaseShiftPerOrbit, timeOfFirstPass, timeOfSecondPass, timeOfFinalPass, timeOfDeployment

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

Eigen::Matrix<float,6,1> propagateActualSatellite(ActualSatellite sat, Eigen::Matrix<float,6,1> x, Eigen::Vector3d u, Eigen::Vector3d r_S2E) {

    //Eigen::Matrix<float,6,1> k1 = satelliteActualOrbitDynamics(sat, x, u, r_S2E);
    //Eigen::Matrix<float,6,1> k2 = satelliteActualOrbitDynamics(sat, x + k1*(ORBIT_DT/2), u, r_S2E);
    //Eigen::Matrix<float,6,1> k3 = satelliteActualOrbitDynamics(sat, x + k2*(ORBIT_DT/2), u, r_S2E);
    //Eigen::Matrix<float,6,1> k4 = satelliteActualOrbitDynamics(sat, x + k3*ORBIT_DT, u, r_S2E);
    //Eigen::Matrix<float,6,1> x_next = x + ORBIT_DT*( (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4 );

    Eigen::Matrix<float,6,1> x_dot = satelliteActualOrbitDynamics(sat, x, u, r_S2E);
    Eigen::Matrix<float,6,1> x_next = x + ORBIT_DT*x_dot;

    return x_next;

}

Eigen::Matrix<float,6,1> propagateReferenceSatellite(ReferenceSatellite sat, Eigen::Matrix<float,6,1> x) {

    //Eigen::Matrix<float,6,1> k1 = satelliteReferenceOrbitDynamics(sat, x);
    //Eigen::Matrix<float,6,1> k2 = satelliteReferenceOrbitDynamics(sat, x + k1*(ORBIT_DT/2));
    //Eigen::Matrix<float,6,1> k3 = satelliteReferenceOrbitDynamics(sat, x + k2*(ORBIT_DT/2));
    //Eigen::Matrix<float,6,1> k4 = satelliteReferenceOrbitDynamics(sat, x + k3*ORBIT_DT);
    //Eigen::Matrix<float,6,1> x_next = x + ORBIT_DT*( (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4 );

    Eigen::Matrix<float,6,1> x_dot = satelliteReferenceOrbitDynamics(sat, x);
    Eigen::Matrix<float,6,1> x_next = x + ORBIT_DT*x_dot;

    return x_next;
}

Eigen::Matrix<float,6,1> propagateEarthState(Eigen::Matrix<float,6,1> x) {

    //Eigen::Matrix<float,6,1> k1 = earthOrbitDynamics(x);
    //Eigen::Matrix<float,6,1> k2 = earthOrbitDynamics(x + k1*(ORBIT_DT/2));
    //Eigen::Matrix<float,6,1> k3 = earthOrbitDynamics(x + k2*(ORBIT_DT/2));
    //Eigen::Matrix<float,6,1> k4 = earthOrbitDynamics(x + k3*ORBIT_DT);
    //Eigen::Matrix<float,6,1> x_next = x + ORBIT_DT*( (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4 );

    Eigen::Matrix<float,6,1> x_dot = earthOrbitDynamics(x);
    Eigen::Matrix<float,6,1> x_next = x + ORBIT_DT*x_dot;

    return x_next;

}


Eigen::Matrix<float,6,1> elements2RV(float a, float e, float i, float O, float w, float f) {

    float rNorm = (a*(1 - pow(e,2)))/(1 + e*cos(f));
    Eigen::Vector3d rPeri;
    rPeri << rNorm*cos(f), rNorm*sin(f), 0;

    float p = a*(1 - pow(e,2));
    Eigen::Vector3d vPeri;
    vPeri << sqrt(MU_E/p)*-1*sin(f), sqrt(MU_E/p)*(e+cos(f)), 0;

    Eigen::Matrix3d R; // rotation matrix from perifocal frame to ECI frame

    R(0,0) = cos(O)*cos(w) - sin(O)*sin(w)*cos(i);
    R(0,1) = -cos(O)*sin(w) - sin(O)*cos(w)*cos(i);
    R(0,2) = sin(O)*sin(i);

    R(1,0) = sin(O)*cos(w) + cos(O)*sin(w)*cos(i);
    R(1,1) = -sin(O)*sin(w) + cos(O)*cos(w)*cos(i);
    R(1,2) = -cos(O)*sin(i);

    R(2,0) = sin(w)*sin(i);
    R(2,1) = cos(w)*sin(i);
    R(2,2) = cos(i);

    Eigen::Vector3d r = R*rPeri;
    Eigen::Vector3d v = R*vPeri;

    Eigen::Matrix<float,6,1> x;
    x << r(0), r(1), r(2), v(0), v(1), v(2);
    return x;

}

Eigen::Matrix<float,5,1> RV2elements(Eigen::Matrix<float,6,1> x) {

    // Returns elements a, e, i, O, w, and T=Theta. Theta is w+f because this package deals with circular orbits.

    Eigen::Vector3d r;
    r << x(0), x(1), x(2);
    Eigen::Vector3d v;
    v << x(3), x(4), x(5);
    Eigen::Vector3d h = r.cross(v);

    float eps = pow(v.norm(), 2)/2 - MU_E/r.norm();

    float a = -MU_E/(2*eps);
    float e;
    if (1 + 2*eps*h.norm()*h.norm()/(pow(MU_E,2)) <= 0) {
        e = 0.00001;
    } else {
        e = sqrt(1 + 2*eps*h.norm()*h.norm()/(pow(MU_E,2)));
    }

    float p = pow(h.norm(), 2)/MU_E;

    float cosI = h(2)/h.norm();
    float i = acos(cosI);

    float O;
    if (h(1) == 0 && h(0) == 0) {
        O = 0;
    } else {
        float hxy = sqrt(pow(h(0), 2) + pow(h(1), 2));
        float sinO = h(0)/hxy;
        float cosO = -h(1)/hxy;
        O = atan2(sinO, cosO);
    }

    float T;
    Eigen::Vector3d n;
    n << cos(O), sin(O), 0;
    float cosT = n.dot(r)/r.norm();
    T = acos(cosT);
    if (r.dot(v) < 0) {
        T = 2*PI - T;
    }

    Eigen::Matrix<float,5,1> elems;
    elems << a, e, i, O, T;
    return elems;

}

Eigen::Matrix3d computeR_ECI_2_LVLH(float i, float O, float T) {

    Eigen::Matrix3d R;
    R(0,0) = cos(O)*cos(T) - sin(O)*sin(T)*cos(i);
    R(0,1) = sin(O)*cos(T) + cos(O)*sin(T)*cos(i);
    R(0,2) = sin(T)*sin(i);
    R(1,0) = -cos(O)*sin(T) - sin(O)*cos(T)*cos(i);
    R(1,1) = -sin(O)*sin(T) + cos(O)*cos(T)*cos(i);
    R(1,2) = cos(T)*sin(i);
    R(2,0) = sin(O)*sin(i);
    R(2,1) = -cos(O)*sin(i);
    R(2,2) = cos(i);

    return R;

}

Eigen::Matrix3d computeR_ECEF_2_ECI(float alpha) {

    Eigen::Matrix3d R;
    R << cos(alpha), -sin(alpha), 0,
         sin(alpha), cos(alpha), 0,
         0, 0, 1;

    return R;

}

Eigen::Matrix3d computeR_ECI_2_ECEF(float alpha) {

    Eigen::Matrix3d R;
    R << cos(alpha), sin(alpha), 0,
         sin(alpha), cos(alpha), 0,
         0, 0, 1;
    return R;

}

Eigen::Vector3d ECEF2LLA(Eigen::Matrix<float,6,1> x_ECEF) {

    Eigen::Vector3d x_LLA;
    Eigen::Vector3d r_ECEF;
    r_ECEF << x_ECEF(0), x_ECEF(1), x_ECEF(2);

    float lat = asin(r_ECEF(2)/r_ECEF.norm());

    float cosLon = r_ECEF(0)/(cos(lat)*r_ECEF.norm());
    float sinLon = r_ECEF(1)/(cos(lat)*r_ECEF.norm());

    float lon = atan2(sinLon, cosLon);

    float alt = r_ECEF.norm() - R_EARTH;

    x_LLA << lat, lon, alt;
    return x_LLA;

}

Eigen::Vector3d LLA2ECEF(Eigen::Vector3d x_LLA) {

    Eigen::Vector3d r_ECEF;

    float x = (R_EARTH + x_LLA(2))*cos(x_LLA(0))*cos(x_LLA(1));
    float y = (R_EARTH + x_LLA(2))*cos(x_LLA(0))*sin(x_LLA(1));
    float z = (R_EARTH + x_LLA(2))*sin(x_LLA(0));

    r_ECEF << x, y, z;
    return r_ECEF;

}
