#include <iostream>
using namespace std;
#include <Eigen/Dense>
#include <fstream>
#include "edl.h"
#include "orbits.h"

int main(){

    // Begin by simulating a ballistic descent sequence, to determine the angle spanned in the orbital plane and time taken by the capsule's
    // travel during ballistic descent. This informs the capsule drop point during satellite orbit.
    std::ofstream ballistic_descent_info_file;
    ballistic_descent_info_file.open("bdinfo.csv");

    float bd_time = 0; // seconds
    float h_chute_deploy = 5000; // meters

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

    ballistic_descent_info_file << bd_time << ","
                            << capsule.x_ballistic(0) << ","
                            << capsule.x_ballistic(1) << ","
                            << capsule.x_ballistic(2) << ","
                            << total_dtheta
                            << "\n";

    while (capsule.x_ballistic(2) >= h_chute_deploy) { // propagate ballistic dynamics and count dtheta for orbit optimization
        bd_time = bd_time + EDL_DT;
        
        total_dtheta = total_dtheta + (capsule.x_ballistic(0)*cos(capsule.x_ballistic(1))*EDL_DT)/(R_EARTH + capsule.x_ballistic(2));
        x_ballistic_next = propagateBallistic(capsule);
        capsule.x_ballistic = x_ballistic_next;

        ballistic_descent_info_file << bd_time << ","
                               << capsule.x_ballistic(0) << ","
                               << capsule.x_ballistic(1) << ","
                               << capsule.x_ballistic(2) << ","
                               << total_dtheta
                               << "\n";
    }
    ballistic_descent_info_file.close();
    // Important parameters for determining orbit stuff: total_dtheta, bd_time

    // Satellite Orbit Initialization Steps
    float globalTime = 0; // s. This variable tracks the total elapsed time of the simulation. At orbit initialization, t = 0.

    // Initialize actual and reference satellites
    float mSat = 3000; // kg
    float CD_Sat = 2.3; // satellite drag coefficient
    float A_Sat = 4; // satellite ram-side cross-sectional area, m^2

    float a0_Sat = ORB_ALT;
    float e0_Sat = 0; // Must be circular or essentially circular for orbit control assumptions to hold valid
    float i0_Sat = 0;
    float O0_Sat = 0;
    float w0_Sat = 0;
    float f0_Sat = 0;
    Eigen::Matrix<float,6,1> x0_Sat = elements2RV(a0_Sat, e0_Sat, i0_Sat, O0_Sat, w0_Sat, f0_Sat);
    
    ActualSatellite sat(mSat, CD_Sat, A_Sat, x0_Sat); 
    ReferenceSatellite refSat(mSat, x0_Sat);

    // Initialize Earth state
    Eigen::Matrix<float,6,1> x0_Earth;
    x0_Earth << 1000*-27590008.8805254, 1000*143923575.304665, 1000*19210.6298233792, 1000*-29.7767856127157, 1000*-5.53623265806239, 1000*-0.000194504010588403; 
    // The above is the state of Earth in the ecliptic J2000 frame at JAN 01 2025 00:00:00.00, obtained using SPICE
    Earth earth(x0_Earth);

    // Initialize orbit parameters
    float desiredDropLocationLat = 38.836587*(PI/180); // rad
    float desiredDropLocationLon = -77.196691*(PI/180); // rad
    float desiredDropLocationAlt = 0; // m, assumes a spherical Earth
    Eigen::Vector3d desiredDropLocation (desiredDropLocationLat, desiredDropLocationLon, desiredDropLocationAlt);

    float latTol = 0.01*(PI/180); // rad. Denotes the tolerance of determining when the satellite is considered to be passing over the latitude of the target location
    float lonTol = 0.05*(PI/180); // rad. Denotes the tolerance of the +/- of longitude error when satellite is passing over target location
    float DThetaTol = 0.05*(PI/180); // rad. Denotes the tolerance of when the satellite is DTheta phase angle from the desired drop point in the orbital plane.
    float timeToDropFromFinalPass = ((2*PI-total_dtheta)/(2*PI))*T_ORBIT;

    float alpha0 = 0; // rad. Denotes the angle between the ECI and ECEF frames at start of orbit propagation
    
    Eigen::Matrix<double, 3, 6> K;
    K << 3.1696, -0.1525, 0, 137.9086, 0.0058, 0,
         0.1525, 3.1586, 0, 0.0058, 137.6687, 0,
         0, 0, 3.1586, 0, 0, 137.6685; // gain matrix for orbit LQR control (constraining actual orbit to reference orbit)
         // The above was derived using the CW frame dynamics and MATLAB's lqr()
    OrbitParams orbitParams(alpha0, K, desiredDropLocation, total_dtheta, DThetaTol, timeToDropFromFinalPass, latTol, lonTol);

    orbitParams.R_ECI_2_LVLH = computeR_ECI_2_LVLH(i0_Sat, O0_Sat, w0_Sat);
    orbitParams.R_LVLH_2_ECI = orbitParams.R_ECI_2_LVLH.transpose();
    orbitParams.R_ECI_2_ECEF = computeR_ECI_2_ECEF(orbitParams.alpha);
    orbitParams.R_ECEF_2_ECI = computeR_ECEF_2_ECI(orbitParams.alpha);

    float timeAtChuteDeployment;
    Eigen::Vector3d r_targetECI;
    Eigen::Vector3d r_DropPointECI;

    // Open up orbital data logging
    std::ofstream orbitsFile;
    orbitsFile.open("orbs.csv");

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


    // Now that all things necessary for orbit propagation have been initialized, begin orbital simulation by orbiting in initial orbit for NUM_ORBITS_BEFORE_INCL_CHANGE
    float numIterationsPerOrbit = round(T_ORBIT/ORBIT_DT);

    for ( int i = 1; i <= NUM_ORBITS_BEFORE_INCL_CHANGE*numIterationsPerOrbit; i++) {

        globalTime = globalTime + ORBIT_DT;
        orbitParams.alpha = W_E*globalTime + orbitParams.alpha0;

        Eigen::Matrix<float,6,1> xSatNext = propagateActualSatellite(sat, sat.x_ECI, sat.u_ECI, earth.r_S2E);
        Eigen::Matrix<float,6,1> xRefSatNext = propagateReferenceSatellite(refSat, refSat.x_ECI);
        Eigen::Matrix<float,6,1> xEarthNext = propagateEarthState(earth.x);
        sat.x_ECI = xSatNext;
        refSat.x_ECI = xRefSatNext;
        earth.x = xEarthNext;
        earth.r_S2E = xEarthNext(Eigen::seq(0,2));

        Eigen::Matrix<float,5,1> orbElems = RV2elements(refSat.x_ECI);
        orbitParams.R_ECI_2_LVLH = computeR_ECI_2_LVLH(orbElems(2), orbElems(3), orbElems(4));
        Eigen::Matrix<float,6,1> dx_ECI = sat.x_ECI - refSat.x_ECI;
        Eigen::Vector3d r_LVLH = orbitParams.R_ECI_2_LVLH*dx_ECI(Eigen::seq(0,2));
        Eigen::Vector3d v_LVLH = orbitParams.R_ECI_2_LVLH*dx_ECI(Eigen::seq(3,5));
        Eigen::Matrix<float,6,1> state_LVLH;
        state_LVLH << r_LVLH(0), r_LVLH(1), r_LVLH(2), v_LVLH(0), v_LVLH(1), v_LVLH(2);
        sat.x_LVLH = state_LVLH;
        sat.u_LVLH = -K*sat.x_LVLH;

        orbitParams.R_LVLH_2_ECI = orbitParams.R_ECI_2_LVLH.transpose();
        sat.u_ECI = orbitParams.R_LVLH_2_ECI*sat.u_LVLH;

        orbitParams.R_ECI_2_ECEF = computeR_ECI_2_ECEF(orbitParams.alpha);
        Eigen::Vector3d r_ECEF = orbitParams.R_ECI_2_ECEF*sat.x_ECI(Eigen::seq(0,2));
        Eigen::Vector3d v_ECEF = orbitParams.R_ECI_2_ECEF*sat.x_ECI(Eigen::seq(3,5));
        Eigen::Matrix<float,6,1> state_ECEF;
        state_ECEF << r_ECEF(0), r_ECEF(1), r_ECEF(2), v_ECEF(0), v_ECEF(1), v_ECEF(2);
        sat.x_ECEF = state_ECEF;

        sat.x_LLA = ECEF2LLA(sat.x_ECEF);
        cout << "Latitude: " << sat.x_LLA(0)*(180/PI);
        cout << "\n";

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
    }

    // Now, compute and execute the Delta-V necessary for the satellite to reach a desired orbit necessary to pass over the drop point.

    float iDesired = orbitParams.desiredDropLocation(0);
    
    // Find Delta-V for actual satellite
<<<<<<< HEAD
    Eigen::Vector3d vCurrentAct;
    vCurrentAct << sat.x_ECI(3), sat.x_ECI(4), sat.x_ECI(5);
    float vAct = vCurrentAct.norm();
    float DVAct = 2*vAct*sin(iDesired/2);
    sat.x_ECI(5) = sat.x_ECI(5) + DVAct;

    // Find Delta-V for reference satellite
    Eigen::Vector3d vCurrentRef;
    vCurrentRef << refSat.x_ECI(3), refSat.x_ECI(4), refSat.x_ECI(5);
    float vRef = vCurrentRef.norm();
    float DVRef = 2*vRef*sin(iDesired/2);
    refSat.x_ECI(5) = sat.x_ECI(5) + DVRef;

    // Find Delta-V for actual satellite
    //Eigen::Vector3d vCurrentAct;
    //vCurrentAct << sat.x_ECI(3), sat.x_ECI(4), sat.x_ECI(5);
    //Eigen::Matrix<float, 5, 1> orbElemsAct = RV2elements(sat.x_ECI);
    //Eigen::Matrix<float, 6, 1> xDesiredAct = elements2RV(orbElemsAct(0), orbElemsAct(1), iDesired, orbElemsAct(3), orbElemsAct(4), 0);
    //Eigen::Vector3d vDesiredAct(xDesiredAct(3), xDesiredAct(4), xDesiredAct(5));
    //Eigen::Vector3d DeltaVAct = vDesiredAct - vCurrentAct;

    // Find Delta-V for reference satellite
    //Eigen::Vector3d vCurrentRef;
    //vCurrentRef << refSat.x_ECI(3), refSat.x_ECI(4), refSat.x_ECI(5);
    //Eigen::Matrix<float, 5, 1> orbElemsRef = RV2elements(refSat.x_ECI);
    //Eigen::Matrix<float, 6, 1> xDesiredRef = elements2RV(orbElemsAct(0), orbElemsAct(1), iDesired, orbElemsAct(3), orbElemsAct(4), 0);
    //Eigen::Vector3d vDesiredRef(xDesiredRef(3), xDesiredRef(4), xDesiredRef(5));
    //Eigen::Vector3d DeltaVRef = vDesiredRef - vCurrentRef;

    // Apply Delta-V
    //sat.x_ECI(3) = sat.x_ECI(3) + DeltaVAct(0);
    //sat.x_ECI(4) = sat.x_ECI(4) + DeltaVAct(1);
    //sat.x_ECI(5) = sat.x_ECI(5) + DeltaVAct(2);
    //refSat.x_ECI(3) = refSat.x_ECI(3) + DeltaVRef(0);
    //refSat.x_ECI(4) = refSat.x_ECI(4) + DeltaVRef(1);
    //refSat.x_ECI(5) = refSat.x_ECI(5) + DeltaVRef(2);
  
=======
    Eigen::Vector3d vCurrentAct = sat.x_ECI(Eigen::seq(3,5));
    Eigen::Matrix<float,5,1> orbElemsAct = RV2elements(sat.x_ECI);
    Eigen::Matrix<float,6,1> xDesiredAct = elements2RV(orbElemsAct(0), orbElemsAct(1), iDesired, orbElemsAct(3), orbElemsAct(4), 0);
    Eigen::Vector3d vDesiredAct = xDesiredAct(Eigen::seq(3,5));
    Eigen::Vector3d DeltaVAct = vDesiredAct - vCurrentAct;

    // Find Delta-V for reference satellite
    Eigen::Vector3d vCurrentRef = refSat.x_ECI(Eigen::seq(3,5));
    Eigen::Matrix<float,5,1> orbElemsRef = RV2elements(refSat.x_ECI);
    Eigen::Matrix<float,6,1> xDesiredRef = elements2RV(orbElemsAct(0), orbElemsAct(1), iDesired, orbElemsAct(3), orbElemsAct(4), 0);
    Eigen::Vector3d vDesiredRef = xDesiredRef(Eigen::seq(3,5));
    Eigen::Vector3d DeltaVRef = vDesiredRef - vCurrentRef;

    // Apply Delta-V
    sat.x_ECI(3) = sat.x_ECI(3) + DeltaVAct(0);
    sat.x_ECI(4) = sat.x_ECI(4) + DeltaVAct(1);
    sat.x_ECI(5) = sat.x_ECI(5) + DeltaVAct(2);
    refSat.x_ECI(3) = refSat.x_ECI(3) + DeltaVRef(0);
    refSat.x_ECI(4) = refSat.x_ECI(4) + DeltaVRef(1);
    refSat.x_ECI(5) = refSat.x_ECI(5) + DeltaVRef(2);

>>>>>>> parent of ca89aaa... No runtime errors on orbit propagation, debugging math/algos now
    // Continue propagating the orbit that passes over the drop point at its zenith until the capsule drops
    while (!orbitParams.hasDeployed) {

        globalTime = globalTime + ORBIT_DT;
        orbitParams.alpha = W_E*globalTime + orbitParams.alpha0;

        Eigen::Matrix<float,6,1> xSatNext = propagateActualSatellite(sat, sat.x_ECI, sat.u_ECI, earth.r_S2E);
        Eigen::Matrix<float,6,1> xRefSatNext = propagateReferenceSatellite(refSat, refSat.x_ECI);
        Eigen::Matrix<float,6,1> xEarthNext = propagateEarthState(earth.x);
        sat.x_ECI = xSatNext;
        refSat.x_ECI = xRefSatNext;
        earth.x = xEarthNext;
        earth.r_S2E = xEarthNext(Eigen::seq(0,2));

        Eigen::Matrix<float,5,1> orbElems = RV2elements(refSat.x_ECI);
        orbitParams.R_ECI_2_LVLH = computeR_ECI_2_LVLH(orbElems(2), orbElems(3), orbElems(4));
        Eigen::Matrix<float,6,1> dx_ECI = sat.x_ECI - refSat.x_ECI;
        Eigen::Vector3d r_LVLH = orbitParams.R_ECI_2_LVLH*dx_ECI(Eigen::seq(0,2));
        Eigen::Vector3d v_LVLH = orbitParams.R_ECI_2_LVLH*dx_ECI(Eigen::seq(3,5));
        Eigen::Matrix<float,6,1> state_LVLH;
        state_LVLH << r_LVLH(0), r_LVLH(1), r_LVLH(2), v_LVLH(0), v_LVLH(1), v_LVLH(2);
        sat.x_LVLH = state_LVLH;
        sat.u_LVLH = -K*sat.x_LVLH;

        orbitParams.R_LVLH_2_ECI = orbitParams.R_ECI_2_LVLH.transpose();
        sat.u_ECI = orbitParams.R_LVLH_2_ECI*sat.u_LVLH;

        orbitParams.R_ECI_2_ECEF = computeR_ECI_2_ECEF(orbitParams.alpha);
        Eigen::Vector3d r_ECEF = orbitParams.R_ECI_2_ECEF*sat.x_ECI(Eigen::seq(0,2));
        Eigen::Vector3d v_ECEF = orbitParams.R_ECI_2_ECEF*sat.x_ECI(Eigen::seq(3,5));
        Eigen::Matrix<float,6,1> state_ECEF;
        state_ECEF << r_ECEF(0), r_ECEF(1), r_ECEF(2), v_ECEF(0), v_ECEF(1), v_ECEF(2);
        sat.x_ECEF = state_ECEF;

        sat.x_LLA = ECEF2LLA(sat.x_ECEF);
<<<<<<< HEAD
        cout << "Latitude: " << sat.x_LLA(0)*(180/PI);
        cout << "\n";
=======
>>>>>>> parent of ca89aaa... No runtime errors on orbit propagation, debugging math/algos now

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

        // If satellite is one orbits' longitudinal phase shift away from passing over desired point, mark it as ready to drop on next orbit
        if (orbitParams.phaseShiftDetermined && abs(orbitParams.desiredDropLocation(0) - sat.x_LLA(0))) {
            if (abs(sat.x_LLA(1) + orbitParams.lonPhaseShiftPerOrbit) <= orbitParams.lonTolerance) {
                orbitParams.readyToDrop = true;
                orbitParams.timeOfFinalPass = globalTime;
            }
        }

        if (orbitParams.readyToDrop && (globalTime - orbitParams.timeOfFinalPass >= orbitParams.timeToDropFromFinalPass) && !orbitParams.hasDeployed) {
            orbitParams.hasDeployed = true;
            r_DropPointECI = sat.x_ECI(Eigen::seq(0,2));
            orbitParams.timeOfDeployment = globalTime;
        }

        if (orbitParams.hasDeployed && (globalTime >= orbitParams.timeOfDeployment+bd_time) && !orbitParams.targetLocAtChuteDeployLogged) {
            Eigen::Vector3d r_targetECEF = LLA2ECEF(orbitParams.desiredDropLocation);
            orbitParams.R_ECEF_2_ECI = computeR_ECEF_2_ECI(orbitParams.alpha);
            r_targetECI = orbitParams.R_ECEF_2_ECI*r_targetECEF;
            timeAtChuteDeployment = globalTime;
            orbitParams.targetLocAtChuteDeployLogged = true;
        }

    }

    // Log one-time variables of interest to a file.
    ofstream oneTimeOrbitParamsFile;
    oneTimeOrbitParamsFile.open("orbitVariablesOfInterest.txt");
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
    oneTimeOrbitParamsFile.close();

    orbitsFile.close();

    // Now, link end of orbits to beginning of EDL ballistic phase.
    std::ofstream ballisticDescentFile;
    std::ofstream chuteDescentFile;
    ballisticDescentFile.open("bd.csv"); // time, V, gamma, h, dtheta_total
    chuteDescentFile.open("cd.csv"); // time, vx, px, vy, py, vz, pz

    float V_0_ballistic = sat.x_ECI(Eigen::seq(3,5)).norm(); // initial velocity, m/s
    float gamma_0_ballistic = 5*(PI/180); // initial flight path angle, rads
    float h_0_ballistic =  sat.x_ECI(Eigen::seq(0,2)).norm() - R_EARTH; // initial altitude, m
    Eigen::Vector3d x_ballistic_init(V_0_ballistic, gamma_0_ballistic, h_0_ballistic); // initial ballistic descent state
    Capsule capsule(m_capsule_init, m_heatshield, m_fuel, CD_capsule, A_capsule, chute, x_ballistic_init);

    // Dynamical perturbations during chute descent
    float vy_perturb = 0.02; // these represent acceleration biases during chute descent, roughly analogous to wind
    float vz_perturb = -0.04; //both vy_perturb and vz_perturb should have absolute value < 0.05
    float chute_deploy_height = 5000; // m
    float powered_descent_height = 1000; // m

    Eigen::Vector3d x_ballistic_next;
    float total_dtheta = 0; // change in angle relative to chute descent point durign ballistic descent
    ballisticDescentFile << globalTime << ","
                            << capsule.x_ballistic(0) << ","
                            << capsule.x_ballistic(1) << ","
                            << capsule.x_ballistic(2) << ","
                            << total_dtheta
                            << "\n";

    while (capsule.x_ballistic(2) >= chute_deploy_height) { // propagate ballistic dynamics and count dtheta for orbit optimization
        total_dtheta = total_dtheta + (capsule.x_ballistic(0)*cos(capsule.x_ballistic(1))*EDL_DT)/(R_EARTH + capsule.x_ballistic(2));
        x_ballistic_next = propagateBallistic(capsule);
        capsule.x_ballistic = x_ballistic_next;
        globalTime = globalTime + EDL_DT;

        ballisticDescentFile << globalTime << ","
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

    ballisticDescentFile.close();

    // Ballistic descent phase complete. Move to chute descent.

    capsule.x_chute << capsule.x_ballistic(0), capsule.x_ballistic(2), 0, 0, 0, 0; // V (vx), h (px), vy, py, vz, pz
    capsule.m = capsule.m - capsule.m_heatshield;
    Eigen::Matrix<float,6,1> x_chute_next;

    chuteDescentFile << globalTime << ","
                       << capsule.x_chute(0) << ","
                       << capsule.x_chute(1) << ","
                       << capsule.x_chute(2) << ","
                       << capsule.x_chute(3) << ","
                       << capsule.x_chute(4) << ","
                       << capsule.x_chute(5) << ","
                       << "\n";


    while (capsule.x_chute(1) >= powered_descent_height) {

        x_chute_next = propagateChute(capsule, vy_perturb, vz_perturb);
        capsule.x_chute = x_chute_next;
        globalTime = globalTime + EDL_DT;

        chuteDescentFile << globalTime << ","
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

    chuteDescentFile.close();

    // Chute descent phase complete. Move to powered descent and landing.
    capsule.m = capsule.m - capsule.chute.m_parachute;

    float Isp = 325; // s
    float Tmax = 5000; // N
    float theta_alt = 87*(PI/180); // radians. Maximum glide angle
    float alpha = 1/(Isp*SEA_LEVEL_G);
    Eigen::Matrix<float,6,1> x_lander_init;
    x_lander_init << capsule.x_chute(1), capsule.x_chute(3), capsule.x_chute(5), -1*capsule.x_chute(0), capsule.x_chute(2), capsule.x_chute(4);
    // The above state is px, py, pz, vx, vy, vz. px points up from the origin, which is the target location.
    // It is assumed that chute descen begins perfectly over the target location
    Lander lander(capsule.m, capsule.m_fuel, Isp, Tmax, theta_alt, alpha, x_lander_init);
    
    // Convex Optimization G-FOLD algorithm for powered descent


    return 0;
}
