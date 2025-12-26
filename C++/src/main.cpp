#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "../include/edl.h"
#include "../include/orbits.h"
#include "../include/orbitDet.h"
#include <thread>
#include <chrono>
#include <Python.h>

using namespace std;

// To compile: g++ -std=c++11 -I path/to/eigen-3.4.0 src/edl.cpp src/orbits.cpp src/orbitDet.cpp src/main.cpp -I/usr/include/python3.8 -lpython3.8 -o main

// To execute: "./main" or "main.exe"

// My path to eigen: ../../../extensions/eigen-3.4.0
// Python version used: 3.8.0

int main(){

    // Initialize orbit parameters
    float desiredDropLocationLat = 38.836587*(PI/180.0); // desired drop latitude, rad
    float desiredDropLocationLon = -77.196691*(PI/180.0); // desired drop longitude, rad
    float desiredDropLocationAlt = 0; // m, assumes a spherical Earth
    Eigen::Vector3d desiredDropLocation (desiredDropLocationLat, desiredDropLocationLon, desiredDropLocationAlt);

    float latTol = 0.01*(PI/180.0); // rad. Denotes the tolerance of determining when the satellite is considered to be passing over the latitude of the target location
    float lonTol = 0.05*(PI/180.0); // rad. Denotes the tolerance of the +/- of longitude error when satellite is passing over target location
    float DThetaTol = 0.05*(PI/180.0); // rad. Denotes the tolerance of when the satellite is DTheta phase angle from the desired drop point in the orbital plane.
    float timeToDropFromFinalPass; // s. Assumes a circular orbit - represents the elapsed time from making the final pass over the drop point latitude to the drop time
    float alpha0 = 0; // rad. Denotes the angle between the ECI and ECEF frames at start of orbit propagation

    float globalTime = 0; // s. This variable tracks the total elapsed time of the simulation. At orbit initialization, t = 0.

    // Initialize actual and reference satellites
    float mSat = 3000; // satellite mass, kg
    float CD_Sat = 2.3; // satellite drag coefficient
    float A_Sat = 4; // satellite ram-side cross-sectional area, m^2

    float a0_Sat = ORB_ALT; // m
    float e0_Sat = 0; // Must be circular or essentially circular for orbit control assumptions to hold valid
    float i0_Sat = 0; // radians
    float O0_Sat = 0; // radians
    float w0_Sat = 0; // radians
    float f0_Sat = 0; // radians
    Eigen::Matrix<float,6,1> x0_Sat = elements2RV(a0_Sat, e0_Sat, i0_Sat, O0_Sat, w0_Sat, f0_Sat);
    ActualSatellite sat(mSat, CD_Sat, A_Sat, x0_Sat); 
    ReferenceSatellite refSat(mSat, x0_Sat);
    Eigen::Matrix<float,5,1> orbElems;
    float numIterationsPerOrbit = round(T_ORBIT/ORBIT_DT); // number of timestep propagations per orbital period

    // Initialize Earth state
    Eigen::Matrix<float,6,1> x0_Earth;
    x0_Earth << 1000*-27590008.8805254, 1000*143923575.304665, 1000*19210.6298233792, 1000*-29.7767856127157, 1000*-5.53623265806239, 1000*-0.000194504010588403; 
    // The above is the state of Earth in the ecliptic J2000 frame at JAN 01 2025 00:00:00.00, obtained using SPICE
    Earth earth(x0_Earth);

    // Initialize array of ground station structs based on input .csv file
    std::string groundStationInputFilename = "../GroundStations.csv";
    std::vector<GroundStation> GroundStations = readInputGroundStationFile(groundStationInputFilename, alpha0);
    int numObservingGroundStations; 

    // Initialize Kalman filter parameters
    Eigen::Matrix<float,6,1> x0Estimate = x0_Sat;
    Eigen::Matrix<float,6,6> P0Estimate;
    P0Estimate << 5.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 5.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f;
    KalmanEstimate kalmanEstimate = KalmanEstimate(x0Estimate, P0Estimate);
    KalmanEstimate kalmanEstimatePrior = kalmanEstimate;

    Eigen::Matrix<float,6,6> Q;
    Q << 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.5f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.5f;
    Eigen::Matrix<float, Eigen::Dynamic, 1> Y;
    float observationChancePerSecond = 0.2; // between 0 and 1; if 1, an observation comes in every time step, and if 0, an observation never comes

    // Initialize orbit controller
    Eigen::Matrix<float, 3, 6> K; // gain matrix for orbit LQR control (constraining actual orbit to reference orbit)
    K << 3.1696f, -0.1525f, 0.0f, 137.9086f, 0.0058f, 0.0f,
        0.1525f, 3.1586f, 0.0f, 0.0058f, 137.6687f, 0.0f,
        0.0f, 0.0f, 3.1586f, 0.0f, 0.0f, 137.6685f;
    // The above was derived using CW frame dynamics and MATLAB's lqr().

    // Variables needed for drop point
    float timeAtChuteDeployment;
    Eigen::Vector3d r_targetECI;
    Eigen::Vector3d r_DropPointECI;

    // Initialize capsule parameters
    float m_capsule_init = 2530; // starting capsule mass, kg. Includes capsule/payload/lander, heatshield, parachute, and fuel.
    float m_heatshield = 500; // heatshield mass, kg
    float m_fuel = 1200; // fuel mass, kg
    float CD_capsule = 1.1; // drag coefficient of capsule alone (during ballistic descent)
    float A_capsule = 2.5; // capsule cross sectional area on descent side, m^2

    // Initialize capsule state
    float V_0_ballistic = sqrt(MU_E/ORB_ALT); // initial velocity, m/s
    float gamma_0_ballistic = 5*(PI/180.0); // initial flight path angle, rads
    float h_0_ballistic = ORB_ALT - R_EARTH; // initial altitude, m
    Eigen::Vector3d V_0_ballistic_vec; // will hold velocity vector at orbit drop point
    Eigen::Vector3d r_0_ballistic_vec; // will hold position vector at orbit drop point

    // Initialize parachute parameters
    float d_chute = 7.5; // diameter of parachute, m
    float CD_chute = 0.8; // drag coefficient of capsule + parachute
    float m_parachute = 30; // mass of parachute, kg
    Parachute chute(d_chute, CD_chute, m_parachute); // Initialize capsule's parachute parameters
    float chute_deploy_height = 5000; // height of chute deployment, meters
    Eigen::Matrix<float,6,1> x_chute_next; // Chute state, [vx, px, vy, py, vz, pz]

    // Dynamical perturbations during chute descent
    float vy_perturb = 0.04; // these represent acceleration biases during chute descent, roughly analogous to wind
    float vz_perturb = -0.05; //both vy_perturb and vz_perturb should have absolute value < 0.07 (for realistic wind impact)

    // Initialize parameters for powered descent
    float powered_descent_height = 1000; // height at which chute descent terminates and powered descent begins, m
    float Isp = 300; // s, denotes fuel efficiency
    float Tmax = 20000; // N, maximum thrust of lander's thrusters. Tmin is 0.3*Tmax
    Eigen::Matrix<float,6,1> x_lander_init; // Lander state, [px py pz vx vy vz]
    float n_Tpoint[] = {-1.0, 0.0, 0.0}; // unit vector for nominal thruster pointing direction. [-1.0, 0.0, 0.0] is straight down.
    float gamma_gs = 20*(PI/180); // radians. Minimum glide angle
    float theta = 75*(PI/180); //radians. Maximum deviation of thrusters from n_Tpoint
    float Vmax = 150; // m/s. Maximum speed allowable for lander
    float tf = 75; // seconds. Time for powered descent trajectory to minimize distance to target with minimum fuel use and zero final velocity.

    // Initialize .csv files for data logging
    std::ofstream ballistic_descent_info_file; // file for precursor ballistic descent simulation information to determine drop point
    ballistic_descent_info_file.open("results/bdinfo.csv"); // time, V, gamma, h, dtheta_total

    std::ofstream orbitsFile; // file for logging parameters during orbital propagation
    orbitsFile.open("results/orbs.csv");

    std::ofstream filterFile; // file for logging state estimate and uncertainties
    filterFile.open("results/filter.csv");

    ofstream oneTimeOrbitParamsFile; // file for logging one-time orbital variables of interest
    oneTimeOrbitParamsFile.open("results/orbitVariablesOfInterest.txt");

    std::ofstream ballisticDescentFile; // file for logging actual ballistic descent data
    ballisticDescentFile.open("results/bd.csv"); // time, V, gamma, h, dtheta_total

    std::ofstream chuteDescentFile; // file for logging chute descent data
    chuteDescentFile.open("results/cd.csv"); // time, vx, px, vy, py, vz, pz

    std::ofstream poweredDescentFile;
    poweredDescentFile.open("results/pd.csv"); // time, rx, ry, rz, vx, vy, vz, m, Tx, Ty, Tz

    // Begin by simulating a ballistic descent sequence, to determine the angle spanned in the orbital plane 
    // and time taken by the capsule's travel during ballistic descent. 
    // This informs the capsule drop point during satellite orbit.
    float bd_time = 0; // tracks ballistic descent precursor simulation, s
    float total_dtheta = 0; // change in angle relative to chute descent point durign ballistic descent, radians

    Eigen::Vector3d x_ballistic_init_sim(V_0_ballistic, gamma_0_ballistic, h_0_ballistic); // initial ballistic descent state
    Capsule capsule_sim(m_capsule_init, m_heatshield, m_fuel, CD_capsule, A_capsule, chute, x_ballistic_init_sim);
    Eigen::Vector3d x_ballistic_next;
    ballistic_descent_info_file << bd_time << ","
                            << capsule_sim.x_ballistic(0) << ","
                            << capsule_sim.x_ballistic(1) << ","
                            << capsule_sim.x_ballistic(2) << ","
                            << total_dtheta
                            << "\n";

    // Propagate ballistic dynamics and count total dtheta to input values to orbit propagation for drop point determination
    while (capsule_sim.x_ballistic(2) >= chute_deploy_height) { 
        bd_time = bd_time + EDL_DT; // update ballistic sim time
        total_dtheta = total_dtheta + (capsule_sim.x_ballistic(0)*cos(capsule_sim.x_ballistic(1))*EDL_DT)/(R_EARTH + capsule_sim.x_ballistic(2)); // cumulatively sum dtheta
        x_ballistic_next = propagateBallistic(capsule_sim); // propagate dynamics by a timestep
        capsule_sim.x_ballistic = x_ballistic_next;
        ballistic_descent_info_file << bd_time << ","
                               << capsule_sim.x_ballistic(0) << ","
                               << capsule_sim.x_ballistic(1) << ","
                               << capsule_sim.x_ballistic(2) << ","
                               << total_dtheta
                               << "\n"; // Log data
    }

    timeToDropFromFinalPass = ((2*PI-total_dtheta)/(2*PI))*T_ORBIT; // Determine amount of time after final pass over drop point latitude to drop

    // Initialize orbit parameters
    OrbitParams orbitParams(alpha0, K, desiredDropLocation, total_dtheta, DThetaTol, timeToDropFromFinalPass, latTol, lonTol);
    orbitParams.R_ECI_2_LVLH = computeR_ECI_2_LVLH(i0_Sat, O0_Sat, w0_Sat);
    orbitParams.R_LVLH_2_ECI = orbitParams.R_ECI_2_LVLH.transpose();
    orbitParams.R_ECI_2_ECEF = computeR_ECI_2_ECEF(orbitParams.alpha);
    orbitParams.R_ECEF_2_ECI = computeR_ECEF_2_ECI(orbitParams.alpha);

    orbitsFile << globalTime << ","
                << earth.x(0) << "," << earth.x(1) << "," << earth.x(2) << "," << earth.x(3) << "," << earth.x(4) << "," << earth.x(5) << ","
                << sat.x_ECI(0) << "," << sat.x_ECI(1) << "," << sat.x_ECI(2) << "," << sat.x_ECI(3) << "," << sat.x_ECI(4) << "," << sat.x_ECI(5) << ","
                << refSat.x_ECI(0) << "," << refSat.x_ECI(1) << "," << refSat.x_ECI(2) << "," << refSat.x_ECI(3) << "," << refSat.x_ECI(4) << "," << refSat.x_ECI(5) << ","
                << sat.u_ECI(0) << "," << sat.u_ECI(1) << "," << sat.u_ECI(2) << ","
                << sat.u_LVLH(0) << "," << sat.u_LVLH(1) << "," << sat.u_LVLH(2) << ","
                << sat.x_LVLH(0) << "," << sat.x_LVLH(1) << "," << sat.x_LVLH(2) << "," << sat.x_LVLH(3) << "," << sat.x_LVLH(4) << "," << sat.x_LVLH(5) << ","
                << sat.x_LLA(0) << "," << sat.x_LLA(1) << "," << sat.x_LLA(2) << ","
                << orbitParams.alpha
                << "\n";

    // Begin orbital simulation by orbiting in initial orbit for NUM_ORBITS_BEFORE_INCL_CHANGE
    cout << "Starting satellite state: " << endl;
    cout << "Position ECI, x (m): " << sat.x_ECI(0) << endl;
    cout << "Position ECI, y (m): " << sat.x_ECI(1) << endl;
    cout << "Position ECI, z (m): " << sat.x_ECI(2) << endl;
    cout << "Velocity ECI, x (m/s): " << sat.x_ECI(3) << endl;
    cout << "Velocity ECI, y (m/s): " << sat.x_ECI(4) << endl;
    cout << "Velocity ECI, z (m/s): " << sat.x_ECI(5) << endl;

    cout << "Satellite orbital parameters before equatorial orbits: " << endl;
    cout << "Semimajor axis, a (m): " << a0_Sat << endl;
    cout << "Eccentricity, e: " << e0_Sat << endl;
    cout << "Inclination, i (deg): " << (180/PI)*i0_Sat << endl;
    cout << "Right Ascension of Ascending Node, Omega (deg): " << (180/PI)*O0_Sat << endl;
    cout << "Argument of Latitude at Epoch, omega + f (deg): " << (180/PI)*(w0_Sat+f0_Sat) << endl;

    for ( int i = 1; i <= NUM_ORBITS_BEFORE_INCL_CHANGE*numIterationsPerOrbit; i++) {

        globalTime = globalTime + ORBIT_DT; // Update global time
        orbitParams.alpha = W_E*globalTime + orbitParams.alpha0; // Update angle from ECI to ECEF frames

        kalmanEstimatePrior = timeUpdate(kalmanEstimate, Q, MU_E, ORBIT_DT); // Perform kalman filter time update

        Eigen::Matrix<float,6,1> xSatNext = propagateActualSatellite(sat, sat.x_ECI, sat.u_ECI, earth.r_S2E); // propagate actual satellite traj
        Eigen::Matrix<float,6,1> xRefSatNext = propagateReferenceSatellite(refSat.x_ECI); // propagate reference traj
        Eigen::Matrix<float,6,1> xEarthNext = propagateEarthState(earth.x); // propagate Earth state
        sat.x_ECI = xSatNext;
        refSat.x_ECI = xRefSatNext;
        earth.x = xEarthNext;
        earth.r_S2E << xEarthNext(0), xEarthNext(1), xEarthNext(2); // update vector from Sun to Earth

        // Update ground station parameters at the current time step
        GroundStations = updateGroundStationLocations(GroundStations, orbitParams.alpha);
        GroundStations = updateGroundStationObservability(GroundStations, sat.x_ECI);
        numObservingGroundStations = getNumberObservingGroundStations(GroundStations);

        // Perform Kalman filter measurement update (or don't, if a measurement is not available)
        float randomNum = static_cast<float>(std::rand()) / RAND_MAX; // Generate random number between 0 and 1 to see if measurement is taken
        if (numObservingGroundStations > 0 & randomNum < observationChancePerSecond){
            Y = getMeasurement(sat.x_ECI, GroundStations, W_E, numObservingGroundStations);
            kalmanEstimate = measurementUpdate(kalmanEstimatePrior, Y, GroundStations, W_E, numObservingGroundStations);
        } else {
            kalmanEstimate = kalmanEstimatePrior;
        }

        orbElems = RV2elements(refSat.x_ECI); // determine orbital elements
        orbitParams.R_ECI_2_LVLH = computeR_ECI_2_LVLH(orbElems(2), orbElems(3), orbElems(4)); // rotation from ECI frame to LVLH frame
        Eigen::Matrix<float,6,1> dx_ECI = sat.x_ECI - refSat.x_ECI;
        Eigen::Vector3d dr_ECI(dx_ECI(0), dx_ECI(1), dx_ECI(2));
        Eigen::Vector3d drdot_ECI(dx_ECI(3), dx_ECI(4), dx_ECI(5));
        Eigen::Vector3d r_LVLH = orbitParams.R_ECI_2_LVLH*dr_ECI; // Find state in LVLH frame
        Eigen::Vector3d v_LVLH = orbitParams.R_ECI_2_LVLH*drdot_ECI;
        Eigen::Matrix<float,6,1> state_LVLH;
        state_LVLH << r_LVLH(0), r_LVLH(1), r_LVLH(2), v_LVLH(0), v_LVLH(1), v_LVLH(2);
        sat.x_LVLH = state_LVLH;

        // Compute control input in LVLH frame
        sat.u_LVLH(0) = -(K(0,0)*sat.x_LVLH(0) + K(0,1)*sat.x_LVLH(1) + K(0,2)*sat.x_LVLH(2) + K(0,3)*sat.x_LVLH(3) + K(0,4)*sat.x_LVLH(4) + K(0,5)*sat.x_LVLH(5));
        sat.u_LVLH(1) = -(K(1,0)*sat.x_LVLH(0) + K(1,1)*sat.x_LVLH(1) + K(1,2)*sat.x_LVLH(2) + K(1,3)*sat.x_LVLH(3) + K(1,4)*sat.x_LVLH(4) + K(1,5)*sat.x_LVLH(5));
        sat.u_LVLH(2) = -(K(2,0)*sat.x_LVLH(0) + K(2,1)*sat.x_LVLH(1) + K(2,2)*sat.x_LVLH(2) + K(2,3)*sat.x_LVLH(3) + K(2,4)*sat.x_LVLH(4) + K(2,5)*sat.x_LVLH(5));

        orbitParams.R_LVLH_2_ECI = orbitParams.R_ECI_2_LVLH.transpose();
        sat.u_ECI = orbitParams.R_LVLH_2_ECI*sat.u_LVLH; // Compute control input in ECI frame - all control/nav is in ECI frame

        orbitParams.R_ECI_2_ECEF = computeR_ECI_2_ECEF(orbitParams.alpha);
        Eigen::Vector3d r_ECI(sat.x_ECI(0), sat.x_ECI(1), sat.x_ECI(2));
        Eigen::Vector3d rdot_ECI(sat.x_ECI(3), sat.x_ECI(4), sat.x_ECI(5));
        Eigen::Vector3d r_ECEF = orbitParams.R_ECI_2_ECEF*r_ECI;
        Eigen::Vector3d v_ECEF = orbitParams.R_ECI_2_ECEF*rdot_ECI;
        Eigen::Matrix<float,6,1> state_ECEF;
        state_ECEF << r_ECEF(0), r_ECEF(1), r_ECEF(2), v_ECEF(0), v_ECEF(1), v_ECEF(2); // Determine ECEF state for logging
        sat.x_ECEF = state_ECEF;
        sat.x_LLA = ECEF2LLA(sat.x_ECEF);

        orbitsFile << globalTime << ","
                << earth.x(0) << "," << earth.x(1) << "," << earth.x(2) << "," << earth.x(3) << "," << earth.x(4) << "," << earth.x(5) << ","
                << sat.x_ECI(0) << "," << sat.x_ECI(1) << "," << sat.x_ECI(2) << "," << sat.x_ECI(3) << "," << sat.x_ECI(4) << "," << sat.x_ECI(5) << ","
                << refSat.x_ECI(0) << "," << refSat.x_ECI(1) << "," << refSat.x_ECI(2) << "," << refSat.x_ECI(3) << "," << refSat.x_ECI(4) << "," << refSat.x_ECI(5) << ","
                << sat.u_ECI(0) << "," << sat.u_ECI(1) << "," << sat.u_ECI(2) << ","
                << sat.u_LVLH(0) << "," << sat.u_LVLH(1) << "," << sat.u_LVLH(2) << ","
                << sat.x_LVLH(0) << "," << sat.x_LVLH(1) << "," << sat.x_LVLH(2) << "," << sat.x_LVLH(3) << "," << sat.x_LVLH(4) << "," << sat.x_LVLH(5) << ","
                << sat.x_LLA(0) << "," << sat.x_LLA(1) << "," << sat.x_LLA(2) << ","
                << orbitParams.alpha
                << "\n"; // log orbit data

        filterFile << globalTime << ","
            << kalmanEstimate.x(0) << "," << kalmanEstimate.x(1) << "," << kalmanEstimate.x(2) << "," 
            << kalmanEstimate.x(3) << "," << kalmanEstimate.x(4) << "," << kalmanEstimate.x(5) << ","
            << 3*sqrt(kalmanEstimate.P(0,0)) << "," << 3*sqrt(kalmanEstimate.P(1,1)) << "," << 3*sqrt(kalmanEstimate.P(2,2)) << ","
            << 3*sqrt(kalmanEstimate.P(3,3)) << "," << 3*sqrt(kalmanEstimate.P(4,4)) << "," << 3*sqrt(kalmanEstimate.P(5,5))
            << "\n"; // Log filter data
    }

    cout << endl;
    cout << "Satellite state after equatorial orbits: " << endl;
    cout << "Position ECI, x (m): " << sat.x_ECI(0) << endl;
    cout << "Position ECI, y (m): " << sat.x_ECI(1) << endl;
    cout << "Position ECI, z (m): " << sat.x_ECI(2) << endl;
    cout << "Velocity ECI, x (m/s): " << sat.x_ECI(3) << endl;
    cout << "Velocity ECI, y (m/s): " << sat.x_ECI(4) << endl;
    cout << "Velocity ECI, z (m/s): " << sat.x_ECI(5) << endl;

    cout << endl;

    cout << "Satellite orbital parameters after equatorial orbits: " << endl;
    cout << "Semimajor axis, a (m): " << orbElems(0) << endl;
    cout << "Eccentricity, e: " << orbElems(1) << endl;
    cout << "Inclination, i (deg): " << (180/PI)*orbElems(2) << endl;
    cout << "Right Ascension of Ascending Node, Omega (deg): " << (180/PI)*orbElems(3) << endl;
    cout << "Argument of Latitude at Epoch, omega + f (deg): " << (180/PI)*orbElems(4) << endl;

    cout << endl;

    cout << "Satellite state estimate after equatorial orbits: " << endl;
    cout << "Position ECI, x (m): " << kalmanEstimate.x(0) << endl;
    cout << "Position ECI, y (m): " << kalmanEstimate.x(1) << endl;
    cout << "Position ECI, z (m): " << kalmanEstimate.x(2) << endl;
    cout << "Velocity ECI, x (m/s): " << kalmanEstimate.x(3) << endl;
    cout << "Velocity ECI, y (m/s): " << kalmanEstimate.x(4) << endl;
    cout << "Velocity ECI, z (m/s): " << kalmanEstimate.x(5) << endl;

    cout << endl;

    // Now, compute and execute the Delta-V necessary for the satellite to reach a desired orbit necessary to pass over the drop point.
    float iDesired = orbitParams.desiredDropLocation(0); // desired orbital inclination is latitude of drop point

    // Find Delta-V for actual satellite
    Eigen::Vector3d vCurrentAct;
    vCurrentAct << sat.x_ECI(3), sat.x_ECI(4), sat.x_ECI(5);
    Eigen::Matrix<float, 5, 1> orbElemsAct = RV2elements(sat.x_ECI);
    Eigen::Matrix<float, 6, 1> xDesiredAct = elements2RV(orbElemsAct(0), orbElemsAct(1), iDesired, orbElemsAct(3), orbElemsAct(4), 0);
    Eigen::Vector3d vDesiredAct(xDesiredAct(3), xDesiredAct(4), xDesiredAct(5));
    Eigen::Vector3d DeltaVAct = vDesiredAct - vCurrentAct;

    // Find Delta-V for reference satellite
    Eigen::Vector3d vCurrentRef;
    vCurrentRef << refSat.x_ECI(3), refSat.x_ECI(4), refSat.x_ECI(5);
    Eigen::Matrix<float, 5, 1> orbElemsRef = RV2elements(refSat.x_ECI);
    Eigen::Matrix<float, 6, 1> xDesiredRef = elements2RV(orbElemsAct(0), orbElemsAct(1), iDesired, orbElemsAct(3), orbElemsAct(4), 0);
    Eigen::Vector3d vDesiredRef(xDesiredRef(3), xDesiredRef(4), xDesiredRef(5));
    Eigen::Vector3d DeltaVRef = vDesiredRef - vCurrentRef;

    // Apply Delta-V
    sat.x_ECI(3) = sat.x_ECI(3) + DeltaVAct(0);
    sat.x_ECI(4) = sat.x_ECI(4) + DeltaVAct(1);
    sat.x_ECI(5) = sat.x_ECI(5) + DeltaVAct(2);
    refSat.x_ECI(3) = refSat.x_ECI(3) + DeltaVRef(0);
    refSat.x_ECI(4) = refSat.x_ECI(4) + DeltaVRef(1);
    refSat.x_ECI(5) = refSat.x_ECI(5) + DeltaVRef(2);

    // Change state estimate based on applied Delta-V
    kalmanEstimate.x(3) = kalmanEstimate.x(3) + DeltaVAct(0);
    kalmanEstimate.x(4) = kalmanEstimate.x(4) + DeltaVAct(1);
    kalmanEstimate.x(5) = kalmanEstimate.x(5) + DeltaVAct(2);

    cout << "Actual Satellite Delta-V applied: " << endl;
    cout << "DVx (m/s): " << DeltaVAct(0) << endl;
    cout << "DVy (m/s): " << DeltaVAct(1) << endl;
    cout << "DVz (m/s): " << DeltaVAct(2) << endl;
    cout << endl;
    cout << "Reference Delta-V applied: " << endl;
    cout << "DVx Ref (m/s): " << DeltaVRef(0) << endl;
    cout << "DVy Ref (m/s): " << DeltaVRef(1) << endl;
    cout << "DVz Ref (m/s): " << DeltaVRef(2) << endl;
    cout << endl;
  
    // Continue propagating the orbit that passes over the drop point at its zenith until the capsule drops
    while (!orbitParams.hasDeployed) {

        globalTime = globalTime + ORBIT_DT;
        orbitParams.alpha = W_E*globalTime + orbitParams.alpha0;
        kalmanEstimatePrior = timeUpdate(kalmanEstimate, Q, MU_E, ORBIT_DT); // Perform kalman filter time update

        Eigen::Matrix<float,6,1> xSatNext = propagateActualSatellite(sat, sat.x_ECI, sat.u_ECI, earth.r_S2E); // propagate actual satellite traj
        Eigen::Matrix<float,6,1> xRefSatNext = propagateReferenceSatellite(refSat.x_ECI); // propagate reference satellite traj
        Eigen::Matrix<float,6,1> xEarthNext = propagateEarthState(earth.x); // propagate Earth state
        sat.x_ECI = xSatNext;
        refSat.x_ECI = xRefSatNext;
        earth.x = xEarthNext;
        earth.r_S2E << xEarthNext(0), xEarthNext(1), xEarthNext(2); // update vector from Sun to Earth

        // Update ground station parameters at the current time step
        groundStations = updateGroundStationLocations(GroundStations, orbitParams.alpha);
        GroundStations = updateGroundStationObservability(GroundStations, sat.x_ECI);
        numObservingGroundStations = getNumberObservingGroundStations(GroundStations);

        // Perform Kalman filter measurement update (or don't, if a measurement is not available)
        float randomNum = static_cast<float>(std::rand()) / RAND_MAX; // Generate random number between 0 and 1 to see if measurement is taken
        if (numObservingGroundStations > 0 & randomNum < observationChancePerSecond) {
            Y = getMeasurement(sat.x_ECI, GroundStations, W_E, numObservingGroundStations);
            kalmanEstimate = measurementUpdate(kalmanEstimatePrior, Y, GroundStations, W_E, numObservingGroundStations);
        } else {
            kalmanEstimate = kalmanEstimatePrior;
        }

        Eigen::Matrix<float,5,1> orbElems = RV2elements(refSat.x_ECI); // determine orbital elements
        Eigen::Matrix<float,6,1> dx_ECI = sat.x_ECI - refSat.x_ECI;
        Eigen::Vector3d dr_ECI(dx_ECI(0), dx_ECI(1), dx_ECI(2));
        Eigen::Vector3d drdot_ECI(dx_ECI(3), dx_ECI(4), dx_ECI(5));
        Eigen::Vector3d r_LVLH = orbitParams.R_ECI_2_LVLH*dr_ECI;
        Eigen::Vector3d v_LVLH = orbitParams.R_ECI_2_LVLH*drdot_ECI;
        Eigen::Matrix<float,6,1> state_LVLH;
        state_LVLH << r_LVLH(0), r_LVLH(1), r_LVLH(2), v_LVLH(0), v_LVLH(1), v_LVLH(2); // determine satellite LVLH state
        sat.x_LVLH = state_LVLH;

        // Compute control input in LVLH frame
        sat.u_LVLH(0) = -(K(0,0)*sat.x_LVLH(0) + K(0,1)*sat.x_LVLH(1) + K(0,2)*sat.x_LVLH(2) + K(0,3)*sat.x_LVLH(3) + K(0,4)*sat.x_LVLH(4) + K(0,5)*sat.x_LVLH(5));
        sat.u_LVLH(1) = -(K(1,0)*sat.x_LVLH(0) + K(1,1)*sat.x_LVLH(1) + K(1,2)*sat.x_LVLH(2) + K(1,3)*sat.x_LVLH(3) + K(1,4)*sat.x_LVLH(4) + K(1,5)*sat.x_LVLH(5));
        sat.u_LVLH(2) = -(K(2,0)*sat.x_LVLH(0) + K(2,1)*sat.x_LVLH(1) + K(2,2)*sat.x_LVLH(2) + K(2,3)*sat.x_LVLH(3) + K(2,4)*sat.x_LVLH(4) + K(2,5)*sat.x_LVLH(5));

        orbitParams.R_LVLH_2_ECI = orbitParams.R_ECI_2_LVLH.transpose();
        sat.u_ECI = orbitParams.R_LVLH_2_ECI*sat.u_LVLH; // compute control input in ECI frame

        orbitParams.R_ECI_2_ECEF = computeR_ECI_2_ECEF(orbitParams.alpha);
        Eigen::Vector3d r_ECI(sat.x_ECI(0), sat.x_ECI(1), sat.x_ECI(2));
        Eigen::Vector3d rdot_ECI(sat.x_ECI(3), sat.x_ECI(4), sat.x_ECI(5));
        Eigen::Vector3d r_ECEF = orbitParams.R_ECI_2_ECEF*r_ECI;
        Eigen::Vector3d v_ECEF = orbitParams.R_ECI_2_ECEF*rdot_ECI;
        Eigen::Matrix<float,6,1> state_ECEF;
        state_ECEF << r_ECEF(0), r_ECEF(1), r_ECEF(2), v_ECEF(0), v_ECEF(1), v_ECEF(2); // determine satellite ECEF state
        sat.x_ECEF = state_ECEF;
        sat.x_LLA = ECEF2LLA(sat.x_ECEF);

        orbitsFile << globalTime << ","
                << earth.x(0) << "," << earth.x(1) << "," << earth.x(2) << "," << earth.x(3) << "," << earth.x(4) << "," << earth.x(5) << ","
                << sat.x_ECI(0) << "," << sat.x_ECI(1) << "," << sat.x_ECI(2) << "," << sat.x_ECI(3) << "," << sat.x_ECI(4) << "," << sat.x_ECI(5) << ","
                << refSat.x_ECI(0) << "," << refSat.x_ECI(1) << "," << refSat.x_ECI(2) << "," << refSat.x_ECI(3) << "," << refSat.x_ECI(4) << "," << refSat.x_ECI(5) << ","
                << sat.u_ECI(0) << "," << sat.u_ECI(1) << "," << sat.u_ECI(2) << ","
                << sat.u_LVLH(0) << "," << sat.u_LVLH(1) << "," << sat.u_LVLH(2) << ","
                << sat.x_LVLH(0) << "," << sat.x_LVLH(1) << "," << sat.x_LVLH(2) << "," << sat.x_LVLH(3) << "," << sat.x_LVLH(4) << "," << sat.x_LVLH(5) << ","
                << sat.x_LLA(0) << "," << sat.x_LLA(1) << "," << sat.x_LLA(2) << ","
                << orbitParams.alpha
                << "\n"; // log orbit data

        filterFile << globalTime << ","
            << kalmanEstimate.x(0) << "," << kalmanEstimate.x(1) << "," << kalmanEstimate.x(2) << "," 
            << kalmanEstimate.x(3) << "," << kalmanEstimate.x(4) << "," << kalmanEstimate.x(5) << ","
            << 3*sqrt(kalmanEstimate.P(0,0)) << "," << 3*sqrt(kalmanEstimate.P(1,1)) << "," << 3*sqrt(kalmanEstimate.P(2,2)) << ","
            << 3*sqrt(kalmanEstimate.P(3,3)) << "," << 3*sqrt(kalmanEstimate.P(4,4)) << "," << 3*sqrt(kalmanEstimate.P(5,5))
            << "\n"; // Log filter data

        // Log time of passage over desired drop point latitude and mark down longitudes at each
        if (abs(orbitParams.desiredDropLocation(0) - sat.x_LLA(0)) <= orbitParams.latTolerance) {
            if (!orbitParams.hasOrbitedOnce) {
                orbitParams.lonAtDesiredLatFirst = sat.x_LLA(1);
                orbitParams.hasOrbitedOnce = true;
                orbitParams.timeOfFirstPass = globalTime;
            } else if (!orbitParams.hasOrbitedTwice && ((globalTime - orbitParams.timeOfFirstPass) >= 0.5*T_ORBIT)) {
                orbitParams.lonAtDesiredLatSecond = sat.x_LLA(1);
                orbitParams.hasOrbitedTwice = true;
                orbitParams.timeOfSecondPass = globalTime;
            }
        }

        // Determine how much the longitude shifts each time the satellite passes over the desired latitude
        if (orbitParams.hasOrbitedOnce && orbitParams.hasOrbitedTwice && !orbitParams.phaseShiftDetermined) {
            orbitParams.lonPhaseShiftPerOrbit = orbitParams.lonAtDesiredLatSecond - orbitParams.lonAtDesiredLatFirst;
            orbitParams.phaseShiftDetermined = true;
        }

        // If satellite is one orbit's longitudinal phase shift away from passing over desired point, mark it as ready to drop on next orbit
        if (orbitParams.phaseShiftDetermined && abs(orbitParams.desiredDropLocation(0) - sat.x_LLA(0))) {
            if (abs(sat.x_LLA(1) + orbitParams.lonPhaseShiftPerOrbit) <= orbitParams.lonTolerance) {
                orbitParams.readyToDrop = true;
                orbitParams.timeOfFinalPass = globalTime;
            }
        }

        // If timeToDropFromFinalPass has elapsed, drop the capsule.
        if (orbitParams.readyToDrop && (globalTime - orbitParams.timeOfFinalPass >= orbitParams.timeToDropFromFinalPass) && !orbitParams.hasDeployed) {
            orbitParams.hasDeployed = true;
            r_DropPointECI << sat.x_ECI(0), sat.x_ECI(1), sat.x_ECI(2);
            orbitParams.timeOfDeployment = globalTime;
            cout << "Dropping capsule!" << endl;
        }

        // Log chute deployment time
        if (orbitParams.hasDeployed && (globalTime >= orbitParams.timeOfDeployment+bd_time) && !orbitParams.targetLocAtChuteDeployLogged) {
            Eigen::Vector3d r_targetECEF = LLA2ECEF(orbitParams.desiredDropLocation);
            orbitParams.R_ECEF_2_ECI = computeR_ECEF_2_ECI(orbitParams.alpha);
            r_targetECI = orbitParams.R_ECEF_2_ECI*r_targetECEF;
            timeAtChuteDeployment = globalTime;
            orbitParams.targetLocAtChuteDeployLogged = true;
        }
    }

    // Log one-time variables of interest to their file.
    oneTimeOrbitParamsFile << "Longitude at first pass over desired latitude (deg): " << orbitParams.lonAtDesiredLatFirst*(180/PI) << "\n";
    oneTimeOrbitParamsFile << "Time of first pass over desired latitude (s): " << orbitParams.timeOfFirstPass << "\n";
    oneTimeOrbitParamsFile << "Longitude at second pass over desired latitude (deg): " << orbitParams.lonAtDesiredLatSecond*(180/PI) << "\n";
    oneTimeOrbitParamsFile << "Time of second pass over desired latitude (s): " << orbitParams.timeOfSecondPass << "\n";
    oneTimeOrbitParamsFile << "Longitude phase shift per orbit (deg): " << orbitParams.lonPhaseShiftPerOrbit*(180/PI) << "\n";
    oneTimeOrbitParamsFile << "Time of final pass over desired latitude (s): " << orbitParams.timeOfFinalPass << "\n";
    oneTimeOrbitParamsFile << "Time Between Final Pass and Capsule Deployment (s): " << orbitParams.timeToDropFromFinalPass << "\n";
    oneTimeOrbitParamsFile << "Ballistic Descent DTheta (deg): " << orbitParams.DTheta*(180/PI) << "\n";
    oneTimeOrbitParamsFile << "Time of Capsule Deployment (s): " << orbitParams.timeOfDeployment << "\n";
    oneTimeOrbitParamsFile << "Time that Ballistic Descent Takes (s): " << bd_time << "\n";
    oneTimeOrbitParamsFile << "Time at Chute Descent Phase Beginning (s): " << timeAtChuteDeployment << "\n";
    oneTimeOrbitParamsFile << "ECI Location of Satellite at Capsule Deployment Point (m):" << r_DropPointECI << "\n";
    oneTimeOrbitParamsFile << "ECI Location of Target Point at Ballistic Descent End (m):" << r_targetECI << "\n";

    cout << endl;
    cout << "Satellite state at point of capsule deployment: " << endl;
    cout << "Position ECI, x: " << sat.x_ECI(0) << endl;
    cout << "Position ECI, y: " << sat.x_ECI(1) << endl;
    cout << "Position ECI, z: " << sat.x_ECI(2) << endl;
    cout << "Velocity ECI, x: " << sat.x_ECI(3) << endl;
    cout << "Velocity ECI, y: " << sat.x_ECI(4) << endl;
    cout << "Velocity ECI, z: " << sat.x_ECI(5) << endl;

    cout << "Satellite orbital parameters after equatorial orbits: " << endl;
    cout << "Semimajor axis, a (m): " << orbElems(0) << endl;
    cout << "Eccentricity, e: " << orbElems(1) << endl;
    cout << "Inclination, i (deg): " << (180/PI)*orbElems(2) << endl;
    cout << "Right Ascension of Ascending Node, Omega (deg): " << (180/PI)*orbElems(3) << endl;
    cout << "Argument of Latitude at Epoch, omega + f (deg): " << (180/PI)*orbElems(4) << endl;
    cout << endl;

    // Now, link end of orbits to beginning of EDL ballistic phase.
    V_0_ballistic_vec << sat.x_ECI(3), sat.x_ECI(4), sat.x_ECI(5);
    r_0_ballistic_vec << sat.x_ECI(0), sat.x_ECI(1), sat.x_ECI(2);

    V_0_ballistic = V_0_ballistic_vec.norm(); // initial velocity, m/s
    h_0_ballistic =  r_0_ballistic_vec.norm() - R_EARTH; // initial altitude, m
    // gamma_0_ballistic determined by capsule separation mechanism --> set at top of main()

    cout << "Capsule state at beginning of ballistic descent: " << endl;
    cout << "Speed (m/s): " << V_0_ballistic << endl;
    cout << "Flight Path Angle (deg): " << gamma_0_ballistic*(180/PI) << endl;
    cout << "Altitude (m): " << h_0_ballistic << endl;
    cout << endl;

    // Initialize actual ballistic descent parameters
    Eigen::Vector3d x_ballistic_init(V_0_ballistic, gamma_0_ballistic, h_0_ballistic); // initial ballistic descent state
    Capsule capsule(m_capsule_init, m_heatshield, m_fuel, CD_capsule, A_capsule, chute, x_ballistic_init);

    total_dtheta = 0; // change in angle relative to chute descent point durign ballistic descent
    ballisticDescentFile << globalTime << ","
                            << capsule.x_ballistic(0) << ","
                            << capsule.x_ballistic(1) << ","
                            << capsule.x_ballistic(2) << ","
                            << total_dtheta
                            << "\n";

    // Propagate ballistic dynamics and count total dtheta
    while (capsule.x_ballistic(2) >= chute_deploy_height) { 
        globalTime = globalTime + EDL_DT; // update global time
        total_dtheta = total_dtheta + (capsule.x_ballistic(0)*cos(capsule.x_ballistic(1))*EDL_DT)/(R_EARTH + capsule.x_ballistic(2)); // cumulatively add dtheta
        x_ballistic_next = propagateBallistic(capsule); // propagate ballistic dynamics by a timestep
        capsule.x_ballistic = x_ballistic_next;
        ballisticDescentFile << globalTime << ","
                               << capsule.x_ballistic(0) << ","
                               << capsule.x_ballistic(1) << ","
                               << capsule.x_ballistic(2) << ","
                               << total_dtheta
                               << "\n";  // log ballistic data
    }
    cout << "Capsule state at end of ballistic descent: " << endl;
    cout << "Speed (m/s): " << capsule.x_ballistic(0) << endl;
    cout << "Flight Path Angle (deg): " << (180/PI)*capsule.x_ballistic(1) << endl;
    cout << "Altitude (m): " << capsule.x_ballistic(2) << endl;
    cout << "Total DTheta over the course of ballistic descent (deg): " << total_dtheta*(180/PI) << endl;
    cout << endl;

    // Ballistic descent phase complete. Move to chute descent.

    capsule.x_chute << capsule.x_ballistic(0), capsule.x_ballistic(2), 0, 0, 0, 0; // V (vx), h (px), vy, py, vz, pz
    capsule.m = capsule.m - capsule.m_heatshield; // heatshield separation

    cout << "Capsule state at beginning of chute descent in local-level frame: " << endl;
    cout << "Position, x (m): " << capsule.x_chute(1) << endl;
    cout << "Position, y (m): " << capsule.x_chute(3) << endl;
    cout << "Position, z (m): " << capsule.x_chute(5) << endl;
    cout << "Velocity, x (m): " << capsule.x_chute(0) << endl;
    cout << "Velocity, y (m): " << capsule.x_chute(2) << endl;
    cout << "Velocity, z (m): " << capsule.x_chute(4) << endl;
    cout << endl;

    chuteDescentFile << globalTime << ","
                       << capsule.x_chute(0) << ","
                       << capsule.x_chute(1) << ","
                       << capsule.x_chute(2) << ","
                       << capsule.x_chute(3) << ","
                       << capsule.x_chute(4) << ","
                       << capsule.x_chute(5) << ","
                       << "\n";

    while (capsule.x_chute(1) >= powered_descent_height) {
        globalTime = globalTime + EDL_DT; // update global time
        x_chute_next = propagateChute(capsule, vy_perturb, vz_perturb); // propagate chute dynamics
        capsule.x_chute = x_chute_next;
        chuteDescentFile << globalTime << ","
                       << capsule.x_chute(0) << ","
                       << capsule.x_chute(1) << ","
                       << capsule.x_chute(2) << ","
                       << capsule.x_chute(3) << ","
                       << capsule.x_chute(4) << ","
                       << capsule.x_chute(5) << ","
                       << "\n"; // log chute descent data
    }
    cout << "Capsule state at end of chute descent/beginning of powered descent in local-level frame: " << endl;
    cout << "Position, x (m): " << capsule.x_chute(1) << endl;
    cout << "Position, y (m): " << capsule.x_chute(3) << endl;
    cout << "Position, z (m): " << capsule.x_chute(5) << endl;
    cout << "Velocity, x (m): " << capsule.x_chute(0) << endl;
    cout << "Velocity, y (m): " << capsule.x_chute(2) << endl;
    cout << "Velocity, z (m): " << capsule.x_chute(4) << endl;
    cout << endl;

    // Chute descent phase complete. Move to powered descent and landing.

    capsule.m = capsule.m - capsule.chute.m_parachute; // Parachute separates to leave just lander, its payload, and fuel.
    x_lander_init << capsule.x_chute(1), capsule.x_chute(3), capsule.x_chute(5), -1*capsule.x_chute(0), capsule.x_chute(2), capsule.x_chute(4); // convert chute state to lander state
    // The above state is px, py, pz, vx, vy, vz. px points up from the origin, which is the target location.
    // It is assumed that chute descent begins perfectly over the target location based on the drop point calculations performed.
    Lander lander(capsule.m, capsule.m_fuel, Isp, Tmax, n_Tpoint, x_lander_init); // initialize lander

    // Values to pass to Python convex optimization function
    float passedVals[20] = {lander.x_lander(0), lander.x_lander(1), lander.x_lander(2), lander.x_lander(3), lander.x_lander(4), lander.x_lander(5),
                                lander.n_Tpoint[0], lander.n_Tpoint[1], lander.n_Tpoint[2],
                                lander.m_wet, lander.m_fuel, lander.Tmin, lander.Tmax, lander.alpha,
                                gamma_gs, theta, Vmax, SEA_LEVEL_G, tf, GFOLD_DT}; // synthesize variables for G-FOLD algorithm
    int passedValsSize = sizeof(passedVals)/sizeof(passedVals[0]);

    // Convex Optimization G-FOLD algorithm for powered descent
    Py_Initialize();
    PyObject* pyList = PyList_New(passedValsSize); // allocate space in memory for values to be passed to convex optimizer
    for (int i = 0; i < passedValsSize; ++i) {
        PyList_SetItem(pyList, i, PyFloat_FromDouble(passedVals[i])); // set values for passed args
    }

    PyRun_SimpleString("import sys; sys.path.append('../Python/');");  // Import path to custom Python module
    PyObject *numpy_module = PyImport_ImportModule("numpy"); // import numpy
    PyObject *cvxpy_module = PyImport_ImportModule("cvxpy"); // import cvxpy
    PyObject* pModule = PyImport_ImportModule("gfoldSolver"); // import custom Python module

    if (!numpy_module || !cvxpy_module || !pModule) {
        PyErr_Print();
        throw std::runtime_error("Import failure");
    }

    PyObject* pFunc = PyObject_GetAttrString(pModule, "solveGfoldOptim"); // import function within Python module

    if (!pFunc || !PyCallable_Check(pFunc)) {
        PyErr_Print();
        throw std::runtime_error("Function error");
    }

    PyObject* pArgs = PyTuple_Pack(1, pyList); // pack passed arguments
    PyObject* returnedList = PyObject_CallObject(pFunc, pArgs); // execute Python code and receive result

    if (!returnedList || !PyList_Check(returnedList)) {
        PyErr_Print();
        throw std::runtime_error("Bad return value");
    }

    // Convert returned PyObject* to a std::vector<float>
    Py_ssize_t size = PyList_Size(returnedList);
    std::vector<float> gfoldDataVec;
    gfoldDataVec.reserve(size);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyList_GetItem(returnedList, i);
        gfoldDataVec.push_back(PyFloat_AsDouble(item));
    }
   
    // Clean up Python interpreter
    Py_DECREF(pArgs);
    Py_DECREF(returnedList);
    Py_DECREF(pFunc);
    Py_DECREF(pyList);
    Py_DECREF(pModule);
    Py_DECREF(numpy_module);
    Py_DECREF(cvxpy_module);

    Py_Finalize();
    // End of G-FOLD optimization
    
    // Copy data into arrays for each variable
    std::size_t numGfoldDatapts = gfoldDataVec.size();
    int numTrajPts = (numGfoldDatapts-1)/10;
    float status_flag = gfoldDataVec[0];
    float rx[numTrajPts];
    float ry[numTrajPts];
    float rz[numTrajPts];
    float vx[numTrajPts];
    float vy[numTrajPts];
    float vz[numTrajPts];
    float m[numTrajPts];
    float Tx[numTrajPts];
    float Ty[numTrajPts];
    float Tz[numTrajPts];
    float pd_time[numTrajPts];
    for (int i = 0; i < numTrajPts; ++i) {
        rx[i] = gfoldDataVec[i+1];
        ry[i] = gfoldDataVec[i+1 + 1*numTrajPts];
        rz[i] = gfoldDataVec[i+1 + 2*numTrajPts];
        vx[i] = gfoldDataVec[i+1 + 3*numTrajPts];
        vy[i] = gfoldDataVec[i+1 + 4*numTrajPts];
        vz[i] = gfoldDataVec[i+1 + 5*numTrajPts];
        m[i] = exp(gfoldDataVec[i+1 + 6*numTrajPts]);
        Tx[i] = gfoldDataVec[i+1 + 7*numTrajPts]*m[i];
        Ty[i] = gfoldDataVec[i+1 + 8*numTrajPts]*m[i];
        Tz[i] = gfoldDataVec[i+1 + 9*numTrajPts]*m[i];

        if (i >= 1) {
            globalTime = globalTime + GFOLD_DT; // update global time
        }
        pd_time[i] = globalTime;
    }

    // Log the data to the .csv file
    for (int i = 0; i < numTrajPts; ++i) {
        poweredDescentFile << pd_time[i] << ","
                       << rx[i] << ","
                       << ry[i] << ","
                       << rz[i] << ","
                       << vx[i] << ","
                       << vy[i] << ","
                       << vz[i] << ","
                       << m[i] << ","
                       << Tx[i] << ","
                       << Ty[i] << ","
                       << Tz[i] << ","
                       << "\n";

        if (i == numTrajPts-1) {
                cout << endl;
                cout << "Capsule state at end powered descent in local-level frame: " << endl;
                cout << "Position, x (m): " << rx[i] << endl;
                cout << "Position, y (m): " << ry[i] << endl;
                cout << "Position, z (m): " << rz[i] << endl;
                cout << "Velocity, x (m/s): " << vx[i] << endl;
                cout << "Velocity, y (m/s): " << vy[i] << endl;
                cout << "Velocity, z (m/s): " << vz[i] << endl;
                cout << endl;
        }
    }
    // Close out data logging files
    ballistic_descent_info_file.close();
    oneTimeOrbitParamsFile.close();
    orbitsFile.close();
    filterFile.close();
    ballisticDescentFile.close();
    chuteDescentFile.close();
    poweredDescentFile.close();

    cout << "G-FOLD Status Flag (1 if Optimal, 0 if Not): " << status_flag << endl;
    if ((int)status_flag == 1) { // if G-FOLD solver finds optimal solution
        cout << "Landed!" << endl;
    }

    return 0;
}
