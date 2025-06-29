#include "edl.h"
#include "orbits.h"
#include <random>
#include <iostream>
using namespace std;


Parachute::Parachute(float d, float CD, float m_parachute) : 
d(d), CD(CD), m_parachute(m_parachute) { A = PI*(d/2)*(d/2); }

Capsule::Capsule(float m, float m_heatshield, float m_fuel, float CD, float A, Parachute chute, Eigen::Vector3d x_ballistic) : 
m(m), m_heatshield(m_heatshield), m_fuel(m_fuel), chute(chute), x_ballistic(x_ballistic) { beta = m/(CD*A); }

Lander::Lander(float m_wet, float m_fuel, float Isp, float Tmax, float theta_alt, float alpha, Eigen::Matrix<float,6,1> x_lander) :
m_wet(m_wet), m_fuel(m_fuel), Isp(Isp), Tmax(Tmax), theta_alt(theta_alt), alpha(alpha), x_lander(x_lander) { Tmin = 0.3*Tmax; }


void Capsule::setBallisticState(float V, float gamma, float h) {
        Capsule::x_ballistic(0) = V;
        Capsule::x_ballistic(1) = gamma;
        Capsule::x_ballistic(2) = h;
    }

void Capsule::setChuteState(float vx, float x, float vy, float y, float vz, float z) {
        x_chute(0) = vx;
        x_chute(1) = x;
        x_chute(2) = vy;
        x_chute(3) = y;
        x_chute(4) = vz;
        x_chute(5) = z;
    }




Eigen::Vector3d ballisticDynamics(Capsule capsule, Eigen::Vector3d x_ballistic){

    float rho = pow(SEA_LEVEL_RHO, -1*x_ballistic(2)/SCALE_HEIGHT);

    float g = SEA_LEVEL_G*pow((R_EARTH/(R_EARTH+x_ballistic(2))), 2);

    Eigen::Vector3d xdot; // state derivative vector: Vdot, gammadot, hdot

    xdot(0) = (-0.5/capsule.beta)*(rho*pow(x_ballistic(0),2)) + g*sin(x_ballistic(1)); //Vdot

    xdot(0) = (-1*rho*x_ballistic(0)*x_ballistic(0))/(2*capsule.beta) + g*sin(x_ballistic(1)); // Vdot
    xdot(1) = (-1*x_ballistic(0)*x_ballistic(0)*cos(x_ballistic(1)))/(x_ballistic(0)*(R_EARTH+x_ballistic(2))) + g*cos(x_ballistic(1))/x_ballistic(0); // gammadot
    xdot(2) = -1*x_ballistic(0)*sin(x_ballistic(1)); // hdot

    return xdot;
}


Eigen::Matrix<float,6,1> chuteDynamics(Capsule capsule, Eigen::Matrix<float,6,1> x, float vy_perturb, float vz_perturb){
    // vy_perturb and vz_perturb should be < 0.1

    float rho = pow(SEA_LEVEL_RHO, -1*x(1)/SCALE_HEIGHT);

    float F_D = 0.5*rho*x(0)*x(0)*capsule.chute.CD*capsule.chute.A; // parachute drag force, N

    Eigen::Matrix<float,6,1> xdot;
    xdot(0) = -1*F_D/capsule.m + SEA_LEVEL_G;
    xdot(1) = -1*x(0);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<> dist_vy(vy_perturb, 1.0); // args = mean, standard deviation
    xdot(2) = dist_vy(gen);
    xdot(3) = x(2) + EDL_DT*xdot(2);

    std::normal_distribution<> dist_vz(vz_perturb, 1.0); // args = mean, standard deviation
    xdot(4) = dist_vz(gen);
    xdot(5) = x(4) + EDL_DT*vz_perturb;

    return xdot;
}

Eigen::Vector3d propagateBallistic(Capsule capsule){

    Eigen::Vector3d k1 = ballisticDynamics(capsule, capsule.x_ballistic);
    Eigen::Vector3d k2 = ballisticDynamics(capsule, capsule.x_ballistic + k1*(EDL_DT/2));
    Eigen::Vector3d k3 = ballisticDynamics(capsule, capsule.x_ballistic + k2*(EDL_DT/2));
    Eigen::Vector3d k4 = ballisticDynamics(capsule, capsule.x_ballistic + k3*EDL_DT);
    Eigen::Vector3d x_next = capsule.x_ballistic + EDL_DT*(k1/6 + k2/3 + k3/3 + k4/6);

    return x_next;
}

Eigen::Matrix<float,6,1> propagateChute(Capsule capsule, float vy_perturb, float vz_perturb){

    Eigen::Matrix<float,6,1> k1 = chuteDynamics(capsule, capsule.x_chute, vy_perturb, vz_perturb);
    Eigen::Matrix<float,6,1> k2 = chuteDynamics(capsule, capsule.x_chute + k1*(EDL_DT/2), vy_perturb, vz_perturb);
    Eigen::Matrix<float,6,1> k3 = chuteDynamics(capsule, capsule.x_chute + k2*(EDL_DT/2), vy_perturb, vz_perturb);
    Eigen::Matrix<float,6,1> k4 = chuteDynamics(capsule, capsule.x_chute + k3*EDL_DT, vy_perturb, vz_perturb);
    Eigen::Matrix<float,6,1> x_next = capsule.x_chute + EDL_DT*(k1/6 + k2/3 + k3/3 + k4/6);

    return x_next;
}

