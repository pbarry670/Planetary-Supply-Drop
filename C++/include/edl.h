#ifndef EDL_H
#define EDL_H

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "orbits.h"

const float EDL_DT = 0.05; // s
const float GFOLD_DT = 0.5; // s
const float SCALE_HEIGHT = 8420; // m
const float SEA_LEVEL_RHO = 1.225; // kg/m^3
const float SEA_LEVEL_G = 9.807; // m/s^2


struct Parachute {
        float d; // parachute diameter, m
        float A; // parachute area, m^2
        float CD; // parachute drag coefficient
        float m_parachute; //parachute mass, kg

        Parachute(float d, float CD, float m_parachute);
    };


struct Capsule {
    float m;
    float m_heatshield;
    float m_fuel;


    float CD; // drag coefficient of capsule during descent
    float A; // cross-sectional area of capsule facing descent side
    float beta; //ballistic coefficient of capsule

    Eigen::Vector3d x_ballistic; // ballistic state vector: V (m/s), gamma (rads), h (m)
    Eigen::Matrix<float, 6, 1> x_chute; // chute state vector: [vx, x, vy, y, vz, z] (m/s, m)

    Parachute chute;

    Capsule(float m, float m_heatshield, float m_fuel, float CD, float A, Parachute chute, Eigen::Vector3d x_ballistic);
    void setBallisticState(float V, float gamma, float h);
    void setChuteState(float vx, float x, float vy, float y, float vz, float z);

};

struct Lander {

    float m_fuel;
    float m_wet;

    float Isp; // specific impulse, s
    float Tmin; // minimum thrust, N
    float Tmax; // maximum thrust, N
    float alpha; // fuel consumption parameter

    float n_Tpoint[3];
    Eigen::Matrix<float,6,1> x_lander;
    
    Lander(float m_wet, float m_fuel, float Isp, float Tmax, float n_Tpoint[3], Eigen::Matrix<float,6,1> x_lander);

};

Eigen::Vector3d ballisticDynamics(Capsule capsule, Eigen::Vector3d x_ballistic);
Eigen::Matrix<float,6,1> chuteDynamics(Capsule capsule, Eigen::Matrix<float,6,1> x_chute, float vy_perturb, float vz_perturb);
Eigen::Vector3d propagateBallistic(Capsule capsule);
Eigen::Matrix<float,6,1> propagateChute(Capsule capsule, float vy_perturb, float vz_perturb);

#endif // EDL_H