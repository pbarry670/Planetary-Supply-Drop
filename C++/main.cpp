#include <iostream>
using namespace std;
#include <Eigen/Dense>
#include "edl.h"
#include "orbits.h"

int main(){

    
    Eigen::Matrix<double, 3, 6> K;
    K << 3.1696, -0.1525, 0, 137.9086, 0.0058, 0,
         0.1525, 3.1586, 0, 0.0058, 137.6687, 0,
         0, 0, 3.1586, 0, 0, 137.6685; // gain matrix for orbit LQR control

    float d_chute = 7.5; // m
    float CD_chute = 0.8; // drag coefficient of capsule + parachute
    float m_parachute = 30; // kg

    Parachute chute(d_chute, CD_chute, m_parachute); // Initialize capsule's parachute parameters

    float m_capsule_init = 2530; // starting capsule mass, kg. Includes capsule/payload/lander, heatshield, parachute, and fuel.
    float m_heatshield = 500; // heatshield mass, kg
    float CD_capsule = 1.1; // capsule drag coefficient
    float A_capsule = 2.5; // capsule cross sectional area on descent side, m^2

    float V_0_ballistic = sqrt(MU_E/(R_EARTH+ORB_ALT)); // initial velocity, m/s
    float gamma_0_ballistic = 5*(PI/180); // initial flight path angle, rads
    float h_0_ballistic = ORB_ALT; // initial altitude, m

    Eigen::Vector3d x_ballistic_init(V_0_ballistic, gamma_0_ballistic, h_0_ballistic); // initial ballistic descent state

    Capsule capsule(m_capsule_init, m_heatshield, CD_capsule, A_capsule, chute, x_ballistic_init);

    // Dynamical perturbations during chute descent
    float vy_perturb = 0.02; // these represent acceleration biases during chute descent, roughly analogous to wind
    float vz_perturb = -0.04; //both vy_perturb and vz_perturb should have absolute value < 0.05

    Eigen::Vector3d x_ballistic_next;
    float total_dtheta = 0; // change in angle relative to chute descent point durign ballistic descent

    while (capsule.x_ballistic(2) >= 5000) { // propagate ballistic dynamics and count dtheta for orbit optimization
        total_dtheta = total_dtheta + (capsule.x_ballistic(0)*cos(capsule.x_ballistic(1))*EDL_DT)/(R_EARTH + capsule.x_ballistic(2));
        // TODO: Log x_ballistic_current or capsule.x_ballistic, either one. Log total dtheta.
        x_ballistic_next = propagateBallistic(capsule);

        //cout << capsule.x_ballistic;
        //cout << x_ballistic_next;

        capsule.x_ballistic = x_ballistic_next;
    }
    // TODO: Log last x_ballistic_current, final total_dtheta
    cout << "State at End of Ballistic Descent: ";
    cout << capsule.x_ballistic;
    cout << "Total DTheta: ";
    cout << total_dtheta*(180/PI);

    // Initialize capsule parameters and state. Initialize vy and vz perturb.
    // Propagate ballistic dynamics, adding up dtheta values
    // At certain altitude, halt ballistic propagation, transfer state over to chute state
    // Propagate chute state
    // At certain altitude, halt chute propagation.

    // Make sure to data log during this whole process to a .csv.


    return 0;
}