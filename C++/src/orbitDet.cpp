#include "../include/edl.h"
#include "../include/orbits.h"
#include "../include/orbitDet.h"
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

const float VISIBILITY_ANGLE = 100.0*(PI/180.0);

GroundStation::GroundStation(std::string locationName, float lat, float lon, float alt, float rangeAcc, float rangeRateAcc, float alpha) : 
    locationName(locationName), lat(lat), lon(lon), alt(alt), rangeAcc(rangeAcc), rangeRateAcc(rangeRateAcc) {

        Eigen::Vector3d x_LLA = {lat, lon, alt};
        r_ECEF = LLA2ECEF(x_LLA);
        Eigen::Matrix3d R_ECEF_2_ECI = computeR_ECEF_2_ECI(alpha);
        r_ECI = R_ECEF_2_ECI*r_ECEF;

}

KalmanEstimate::KalmanEstimate(Eigen::Matrix<float,6,1> x, Eigen::Matrix<float,6,6> P) : x(x), P(P) {}


std::vector<GroundStation> readInputGroundStationFile(std::string groundStationInputFilename, float alpha){

    std::vector<GroundStation> groundStations;
    std::ifstream file(groundStationInputFilename);

    std::string line;


    if (!file.is_open()) {
        throw std::runtime_error("Could not open CSV file: " + groundStationInputFilename);
    }

    while (std::getline(file, line)) {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::string token;

        // Station name
        std::string locationName;
        std::getline(ss, locationName, ',');

        // Latitude
        std::getline(ss, token, ',');
        float lat = (PI/180.0)*std::stof(token);

        // Longitude
        std::getline(ss, token, ',');
        float lon = (PI/180.0)*std::stof(token);

        // Altitude
        std::getline(ss, token, ',');
        float alt = std::stod(token);

        // Range accuracy
        std::getline(ss, token, ',');
        float rangeAcc = std::stof(token);

        // Range rate accuracy
        std::getline(ss, token, ',');
        float rangeRateAcc = std::stof(token);

        GroundStation groundStation = GroundStation(locationName, lat, lon, alt, rangeAcc, rangeRateAcc, alpha);
        groundStations.push_back(groundStation);

    }
    return groundStations;

}

KalmanEstimate timeUpdate(KalmanEstimate posterior, Eigen::Matrix<float,6,6> Q, float mu, float dt){

    float r = sqrt(pow(posterior.x(0), 2) + pow(posterior.x(1), 2) + pow(posterior.x(2), 2));

    Eigen::Matrix<float,6,6> A;
    A << 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
        -mu/pow(r,3) + (3*mu*pow(posterior.x(0),2))/(pow(r,5)), (3*mu*posterior.x(0)*posterior.x(1))/(pow(r,5)), (3*mu*posterior.x(0)*posterior.x(2))/(pow(r,5)), 0.0f, 0.0f, 0.0f,
        (3*mu*posterior.x(0)*posterior.x(1))/(pow(r,5)), -mu/pow(r,3) + (3*mu*pow(posterior.x(1),2))/(pow(r,5)), (3*mu*posterior.x(1)*posterior.x(2))/(pow(r,5)), 0.0f, 0.0f, 0.0f,
        (3*mu*posterior.x(0)*posterior.x(2))/pow(r,5), (3*mu*posterior.x(1)*posterior.x(2))/pow(r,5), -mu/pow(r,3) + (3*mu*pow(posterior.x(2),2))/(pow(r,5)), 0.0f, 0.0f, 0.0f;

    Eigen::Matrix<float,6,6> Phi = Eigen::Matrix<float,6,6>::Identity() + A*dt + (1.0/2.0)*A*A*(pow(dt,2)) + (1.0/6.0)*A*A*A*(pow(dt,3)) + (1.0/24.0)*A*A*A*A*(pow(dt,4));

    Eigen::Matrix<float,6,1> xPrior = propagateReferenceSatellite(posterior.x);
    Eigen::Matrix<float,6,6> PPrior = Phi*posterior.P*Phi.transpose() + Q;
    KalmanEstimate prior = KalmanEstimate(xPrior, PPrior);
    return prior;

}

std::vector<GroundStation> updateGroundStationLocations(std::vector<GroundStation> GroundStations, float alpha){

    Eigen::Matrix3d R_ECEF_2_ECI = computeR_ECEF_2_ECI(alpha);
    for (size_t i = 0; i < GroundStations.size(); i++){
        GroundStations[i].r_ECI = R_ECEF_2_ECI*GroundStations[i].r_ECEF;
    }
    return GroundStations;

}

std::vector<GroundStation> updateGroundStationObservability(std::vector<GroundStation> GroundStations, Eigen::Matrix<float,6,1> xSatECI){

    Eigen::Vector3d rSatECI = {xSatECI(0), xSatECI(1), xSatECI(2)};

    for (size_t i = 0; i < GroundStations.size(); i++){
        Eigen::Vector3d rGroundStationECI = GroundStations[i].r_ECI;
        Eigen::Vector3d rGS2SatECI = rSatECI - rGroundStationECI;

        float rGroundStationECIMag = sqrt(pow(rGroundStationECI(0), 2) + pow(rGroundStationECI(1), 2) + pow(rGroundStationECI(2), 2));
        float rGS2SatECIMag = sqrt(pow(rGS2SatECI(0), 2) + pow(rGS2SatECI(1), 2) + pow(rGS2SatECI(2), 2));

        float thetaObs = acos((rGS2SatECI(0)*rGroundStationECI(0) + rGS2SatECI(1)*rGroundStationECI(1) + rGS2SatECI(2)*rGroundStationECI(2))/(rGS2SatECIMag*rGroundStationECIMag));
        if (thetaObs <= VISIBILITY_ANGLE) {
            GroundStations[i].isObserving = true;
        } else {
            GroundStations[i].isObserving = false;
        }
    }
    return GroundStations;
}

int getNumberObservingGroundStations(std::vector<GroundStation> GroundStations){
    
    int NObs = 0;
    for (size_t i = 0; i < GroundStations.size(); i++){
        if (GroundStations[i].isObserving) {
            NObs = NObs + 1;
        }
    }
    return NObs;

}

Eigen::Matrix<float,Eigen::Dynamic,1> getMeasurement(Eigen::Matrix<float,6,1> xSatECI, std::vector<GroundStation> GroundStations, float wE, int NObs){

    Eigen::Matrix<float, Eigen::Dynamic, 1> Y;
    Y.resize(2*NObs, 1);
    int index = 0;

    for (size_t i = 0; i < GroundStations.size(); i++){
        if (GroundStations[i].isObserving) {

            Eigen::Vector3d rStation = GroundStations[i].r_ECI;
            Eigen::Vector3d vStation = {-wE*rStation(1), wE*rStation(0), 0.0};

            float rho_i = sqrt( pow((xSatECI(0) - rStation(0)),2) + pow((xSatECI(1) - rStation(1)),2) + pow((xSatECI(2) - rStation(2)),2) );
            float rhodot_i = ((xSatECI(0) - rStation(0))*(xSatECI(3) - vStation(0)) + (xSatECI(1) - rStation(1))*(xSatECI(4) - vStation(1)) + (xSatECI(2) - rStation(2))*(xSatECI(5) - vStation(2)))/rho_i;


            std::random_device rd;   // Non-deterministic seed
            std::mt19937 gen(rd());  // Mersenne Twister PRNG
            std::normal_distribution<float> rangeDist(rho_i, GroundStations[i].rangeAcc/2.0);
            std::normal_distribution<float> rangeRateDist(rhodot_i, GroundStations[i].rangeRateAcc/2.0);

            Y(index) = rangeDist(gen);
            Y(index+1) = rangeRateDist(gen);
            index = index+2;
        }
    }

    return Y;

}

Eigen::Matrix<float,Eigen::Dynamic,1> measurementModel(Eigen::Matrix<float,6,1> xSatECI, std::vector<GroundStation> GroundStations, float wE, int NObs){

    Eigen::Matrix<float, Eigen::Dynamic, 1> Y;
    Y.resize(2*NObs, 1);
    int index = 0;

    for (size_t i = 0; i < GroundStations.size(); i++){
        if (GroundStations[i].isObserving) {

            Eigen::Vector3d rStation = GroundStations[i].r_ECI;
            Eigen::Vector3d vStation = {-wE*rStation(1), wE*rStation(0), 0.0};

            float rho_i = sqrt( pow((xSatECI(0) - rStation(0)),2) + pow((xSatECI(1) - rStation(1)),2) + pow((xSatECI(2) - rStation(2)),2) );
            float rhodot_i = ((xSatECI(0) - rStation(0))*(xSatECI(3) - vStation(0)) + (xSatECI(1) - rStation(1))*(xSatECI(4) - vStation(1)) + (xSatECI(2) - rStation(2))*(xSatECI(5) - vStation(2)))/rho_i;

            Y(index) = rho_i;
            Y(index+1) = rhodot_i;
            index = index+2;
        }
    }
    return Y;
}

KalmanEstimate measurementUpdate(KalmanEstimate prior, Eigen::Matrix<float, Eigen::Dynamic, 1> Y, std::vector<GroundStation> GroundStations, float wE, int NObs){

    Eigen::Matrix<float, Eigen::Dynamic, 1> y;
    y.resize(2*NObs, 1);
    y = Y - measurementModel(prior.x, GroundStations, wE, NObs);

    Eigen::Matrix<float, Eigen::Dynamic, 6> H;
    H.resize(2*NObs, 6);
    H.setZero();

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> R;
    R.resize(2*NObs, 2*NObs);
    R.setZero();

    int rowIndex = 0;
    for (size_t i = 0; i < GroundStations.size(); i++){
        if (GroundStations[i].isObserving) {
            Eigen::Vector3d rStation = GroundStations[i].r_ECI;
            Eigen::Vector3d vStation = {-wE*rStation(1), wE*rStation(0), 0.0};

            H(rowIndex, 0) = -(0.5*(2.0*rStation(0) - 2.0*prior.x(0)))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 1) = -(0.5*(2.0*rStation(1) - 2.0*prior.x(1)))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 2) = -(0.5*(2.0*rStation(2) - 2.0*prior.x(2)))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            R(rowIndex, rowIndex) = GroundStations[i].rangeAcc/2.0;

            rowIndex = rowIndex + 1;

            H(rowIndex, 0) = (0.5*(2.0*rStation(0) - 2.0*prior.x(0))*((rStation(0) - prior.x(0))*(vStation(0) - prior.x(3)) + (rStation(1) - prior.x(1))*(vStation(1) - prior.x(4)) + (rStation(2) - prior.x(2))*(vStation(2) - prior.x(5))))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),1.5) - (vStation(0) - prior.x(3))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 1) = (0.5*(2.0*rStation(1) - 2.0*prior.x(1))*((rStation(0) - prior.x(0))*(vStation(0) - prior.x(3)) + (rStation(1) - prior.x(1))*(vStation(1) - prior.x(4)) + (rStation(2) - prior.x(2))*(vStation(2) - prior.x(5))))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),1.5) - (vStation(1) - prior.x(4))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 2) = (0.5*(2.0*rStation(2) - 2.0*prior.x(2))*((rStation(0) - prior.x(0))*(vStation(0) - prior.x(3)) + (rStation(1) - prior.x(1))*(vStation(1) - prior.x(4)) + (rStation(2) - prior.x(2))*(vStation(2) - prior.x(5))))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),1.5) - (vStation(2) - prior.x(5))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 3) = -(rStation(0) - prior.x(0))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 4) = -(rStation(1) - prior.x(1))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 5) = -(rStation(2) - prior.x(2))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            R(rowIndex, rowIndex) = GroundStations[i].rangeRateAcc/2.0;

            rowIndex = rowIndex + 1;

        }
    }

    Eigen::Matrix<float,6,Eigen::Dynamic> K;
    K.resize(6, 2*NObs);

    K = prior.P*H.transpose()*(H*prior.P*H.transpose() + R).inverse();
    Eigen::Matrix<float,6,1> xPosterior = prior.x + K*y;
    Eigen::Matrix<float,6,6> PPosterior = (Eigen::Matrix<float,6,6>::Identity() - K*H)*prior.P*(Eigen::Matrix<float,6,6>::Identity() - K*H).transpose() + K*R*K.transpose();

    return KalmanEstimate(xPosterior, PPosterior);

}