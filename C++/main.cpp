#include <iostream>
using namespace std;
#include <Eigen/Dense>
#include <fstream>
#include "edl.h"
#include "orbits.h"

int main(){

    std::ofstream ballistic_descent_file;
    std::ofstream chute_descent_file;

    float global_time = 0; // seconds

    ballistic_descent_file.open("bd.csv"); // time, V, gamma, h, dtheta_total
    chute_descent_file.open("cd.csv"); // time, vx, px, vy, py, vz, pz

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
    float m_fuel = 900; // fuel mass, kg
    float CD_capsule = 1.1; // capsule drag coefficient
    float A_capsule = 2.5; // capsule cross sectional area on descent side, m^2

    float V_0_ballistic = sqrt(MU_E/(R_EARTH+ORB_ALT)); // initial velocity, m/s
    float gamma_0_ballistic = 5*(PI/180); // initial flight path angle, rads
    float h_0_ballistic = ORB_ALT; // initial altitude, m

    Eigen::Vector3d x_ballistic_init(V_0_ballistic, gamma_0_ballistic, h_0_ballistic); // initial ballistic descent state

    Capsule capsule(m_capsule_init, m_heatshield, m_fuel, CD_capsule, A_capsule, chute, x_ballistic_init);

    // Dynamical perturbations during chute descent
    float vy_perturb = 0.02; // these represent acceleration biases during chute descent, roughly analogous to wind
    float vz_perturb = -0.04; //both vy_perturb and vz_perturb should have absolute value < 0.05

    Eigen::Vector3d x_ballistic_next;
    float total_dtheta = 0; // change in angle relative to chute descent point durign ballistic descent

    ballistic_descent_file << global_time << ","
                            << capsule.x_ballistic(0) << ","
                            << capsule.x_ballistic(1) << ","
                            << capsule.x_ballistic(2) << ","
                            << total_dtheta
                            << "\n";

    while (capsule.x_ballistic(2) >= 5000) { // propagate ballistic dynamics and count dtheta for orbit optimization
        total_dtheta = total_dtheta + (capsule.x_ballistic(0)*cos(capsule.x_ballistic(1))*EDL_DT)/(R_EARTH + capsule.x_ballistic(2));
        x_ballistic_next = propagateBallistic(capsule);
        capsule.x_ballistic = x_ballistic_next;
        global_time = global_time + EDL_DT;

        ballistic_descent_file << global_time << ","
                               << capsule.x_ballistic(0) << ","
                               << capsule.x_ballistic(1) << ","
                               << capsule.x_ballistic(2) << ","
                               << total_dtheta
                               << "\n";
    }

    cout << "State at End of Ballistic Descent: ";
    cout << capsule.x_ballistic;
    cout << "Total DTheta: ";
    cout << total_dtheta*(180/PI);

    // BALLISTIC DESCENT COMPLETE. MOVE TO CHUTE DESCENT.

    capsule.x_chute << capsule.x_ballistic(0), capsule.x_ballistic(2), 0, 0, 0, 0; // V (vx), h (px), vy, py, vz, pz
    capsule.m = capsule.m - capsule.m_heatshield;
    Eigen::Matrix<float,6,1> x_chute_next;

    chute_descent_file << global_time << ","
                       << capsule.x_chute(0) << ","
                       << capsule.x_chute(1) << ","
                       << capsule.x_chute(2) << ","
                       << capsule.x_chute(3) << ","
                       << capsule.x_chute(4) << ","
                       << capsule.x_chute(5) << ","
                       << "\n";


    while (capsule.x_chute(1) >= 1000) {

        x_chute_next = propagateChute(capsule, vy_perturb, vz_perturb);
        capsule.x_chute = x_chute_next;
        global_time = global_time + EDL_DT;

        chute_descent_file << global_time << ","
                       << capsule.x_chute(0) << ","
                       << capsule.x_chute(1) << ","
                       << capsule.x_chute(2) << ","
                       << capsule.x_chute(3) << ","
                       << capsule.x_chute(4) << ","
                       << capsule.x_chute(5) << ","
                       << "\n";
    }

    cout << "State at End of Chute Descent: ";
    cout << capsule.x_chute;

    // CHUTE DESCENT COMPLETE, MOVE TO POWERED DESCENT

    capsule.m = capsule.m - capsule.chute.m_parachute;

    float Isp = 325; // s
    float Tmax = 5000; // N
    float theta_alt = 87*(PI/180); // radians. Maximum glide angle
    float alpha = 1/(Isp*SEA_LEVEL_G);
    Eigen::Matrix<float,6,1> x_lander_init;
    x_lander_init << capsule.x_chute(1), capsule.x_chute(3), capsule.x_chute(5), -1*capsule.x_chute(0), capsule.x_chute(2), capsule.x_chute(4);

    Lander lander(capsule.m, capsule.m_fuel, Isp, Tmax, theta_alt, alpha, x_lander_init);
    // Convex Optimization G-FOLD algorithm for powered descent


    ballistic_descent_file.close();
    chute_descent_file.close();

    return 0;
}