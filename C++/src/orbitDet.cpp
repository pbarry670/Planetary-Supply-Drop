#include "../include/edl.h"
#include "../include/orbits.h"
#include "../include/orbitDet.h"
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

const float groundStationVisibilityAngle = 100.0;

GroundStation::GroundStation(string locationName, float lat, float lon, float alt, float rangeAcc, float rangeRateAcc, float alpha) : 
    locationName(locationName), lat(lat), lon(lon), alt(alt), rangeAcc(rangeAcc), rangeRateAcc(rangeRateAcc) {

        Eigen::Vector3d x_LLA = {lat, lon, alt};
        r_ECEF = LLA2ECEF(x_LLA);
        Eigen::Matrix3d R_ECEF_2_ECI = computeR_ECEF_2_ECI(alpha);
        r_ECI = R_ECEF_2_ECI*r_ECEF;

}

KalmanEstimate::KalmanEstimate(Eigen::Matrix<float,6,1> x, Eigen::Matrix<float,6,6> P) : x(x), P(P) {}


std::vector<GroundStation> readInputGroundStationFile(string groundStationInputFilename, float alpha0){

    fstream fin;
    fin.open(groundStationInputFilename, ios::in);

    std::vector<string> row;
    std::vector<GroundStation> GroundStations;
    GroundStations.clear();

    string line, word, temp;
    while (fin >> temp){
        row.clear();

        getline(fin, line); //Read an entire row of the .csv and store in a string variable "line"
        stringstream s(line); // Break the elements of the .csv apart

        while (getline(s, word, ',')){ // Read every column of a row and store it in "row"
            row.push_back(word);
        }

        string locName = row[0];
        float lat = stof(row[1]);
        float lon = stof(row[2]);
        float alt = stof(row[3]);
        float rangeAcc = stof(row[4]);
        float rangeRateAcc = stof(row[5]);
        GroundStation groundStation_i = GroundStation(locName, lat, lon, alt, rangeAcc, rangeRateAcc, alpha0);
        GroundStations.push_back(groundStation_i);

    }
    fin.close();
    return GroundStations;

}

KalmanEstimate timeUpdate(kalmanEstimate posterior, Eigen::Matrix<float,6,6> Q, float mu, float dt){

    float r = sqrt(pow(posterior.x(0), 2) + pow(posterior.x(1), 2) + pow(posterior.x(2), 2));

    Eigen::Matrix<float,6,6> A = [0, 0, 0, 1, 0, 0;
                                0, 0, 0, 0, 1, 0;
                                0, 0, 0, 0, 0, 1;
                                -mu/pow(r,3) + (3*mu*pow(posterior.x(0),2))/(pow(r,5)), (3*mu*posterior.x(0)*posterior.x(1))/(pow(r,5)), (3*mu*posterior.x(0)*posterior.x(2))/(pow(r,5)), 0, 0, 0;
                                (3*mu*posterior.x(0)*posterior.x(1))/(pow(r,5)), -mu/r^3 + (3*mu*pow(posterior.x(1),2))/(pow(r,5)), (3*mu*posterior.x(1)*posterior.x(2))/(pow(r,5)), 0, 0, 0;
                                (3*mu*posterior.x(0)*posterior.x(2))/(r^5), (3*mu*posterior.x(1)*posterior.x(2))/(r^5), -mu/r^3 + (3*mu*pow(posterior.x(2),2))/(pow(r,5)), 0, 0, 0];

    Phi = Eigen::Matrix<float,6,6>::Identity() + A*dt + (1.0/2.0)*A*A*(pow(dt,2)) + (1.0/6.0)*A*A*A*(pow(dt,3)) + (1.0/24.0)*A*A*A*A*(pow(dt,4));

    Eigen::Matrix<float,6,1> xPrior = propagateReferenceSatellite(posterior.x);
    Eigen::Matrix<float,6,6> PPrior = Phi*posterior.P*Phi.transpose() + Q;
    KalmanEstimate prior = KalmanEstimate(xPrior, PPrior)
    return prior;

}

std::vector<GroundStation> updateGroundStationLocations(std::vector<GroundStation> GroundStations, float alpha){

    Eigen::Matrix3d R_ECEF_2_ECI = compute_R_ECEF_2_ECI(alpha);
    for (size_t i = 0; i < GroundStations.size(); i++){
        GroundStations[i].r_ECI = R_ECEF_2_ECI*GroundStations[i].r_ECI;
    }
    return GroundStations;

}

std::vector<GroundStation> updateGroundStationObservability(std::vector<GroundStation> GroundStations, Eigen::Matrix<float,6,1> xSatECI){

    Eigen::Vector3d rSatECI = {xSatECI(0), xSatECI(1), xSatECI(2)};
    float rSatECIMag = sqrt(pow(rSatECI(0), 2) + pow(rSatECI(1), 2) + pow(rSatECI(2), 2));
    for (size_t i = 0; i < GroundStations.size(); i++){
        Eigen::Vector3d rGroundStationECI = GroundStations[i].r_ECI;
        float rGroundStationECIMag = sqrt(pow(rGroundStationECI(0), 2) + pow(rGroundStationECI(1), 2) + pow(rGroundStationECI(2), 2));
        float thetaObs = acos((rSatECI(0)*rGroundStationECI(0) + rSatECI(1)*rGroundStationECI(1) + rSatECI(2)*rGroundStationECI(2))/(rSatECIMag*rGroundStationECIMag));
        if thetaObs <= groundStationVisibilityAngle{
            GroundStations[i].isObserving = true;
        }
    }
    return GroundStations;
}

int getNumberObservingGroundStations(std::vector<GroundStation> GroundStations){
    
    int NObs = 0;
    for (size_t i = 0; i < GroundStations.size(); i++){
        if GroundStations[i].isObserving{
            NObs = NObs + 1;
        }
    }
    return NObs;

}

Eigen::Matrix<float,Dynamic,1> getMeasurement(Eigen::Matrix<float,6,1> xSatECI, std::vector<GroundStation> GroundStations, float wE, int NObs){

    Eigen::Matrix<float, 2*NObs, 1> Y;
    int index = 0;
    for (size_t i = 0; i < GroundStations.size(); i++){
        if GroundStations[i].isObserving{

            Eigen::Vector3d rStation = GroundStations[i].r_ECI;
            Eigen::Vector3d vStation = {-wE*rStation(1), wE*rStation(0), 0.0};

            float rho_i = sqrt( pow((xSatECI(0) - rStation(0)),2) + pow((xSat(1) - rStation(1)),2) + pow((xSatECI(2) - rStation(2)),2) );
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

Eigen::Matrix<float,Dynamic,1> measurementModel(Eigen::Matrix<float,6,1> xSatECI, std::vector<GroundStation> GroundStations, float wE, int NObs){

    Eigen::Matrix<float, 2*NObs, 1> Y;
    int index = 0;
    for (size_t i = 0; i < GroundStations.size(); i++){
        if GroundStations[i].isObserving{

            Eigen::Vector3d rStation = GroundStations[i].r_ECI;
            Eigen::Vector3d vStation = {-wE*rStation(1), wE*rStation(0), 0.0};

            float rho_i = sqrt( pow((xSatECI(0) - rStation(0)),2) + pow((xSat(1) - rStation(1)),2) + pow((xSatECI(2) - rStation(2)),2) );
            float rhodot_i = ((xSatECI(0) - rStation(0))*(xSatECI(3) - vStation(0)) + (xSatECI(1) - rStation(1))*(xSatECI(4) - vStation(1)) + (xSatECI(2) - rStation(2))*(xSatECI(5) - vStation(2)))/rho_i;

            Y(index) = rho_i;
            Y(index+1) = rhodot_i;
            index = index+2;
        }
    }
    return Y;
}

KalmanEstimate measurementUpdate(kalmanEstimate prior, Eigen::Matrix<float, Dynamic, 1> Y, std::vector<GroundStation> GroundStations, float wE, int NObs){

    Eigen::Matrix<float, 2*NObs, 1> y = Y - measurementModel(prior.x, GroundStations, wE, NObs);

    Eigen::Matrix<float, 2*NObs, 6> H;
    H.setZero();

    Eigen::Matrix<float, 2*NObs, 2*NObs> R;
    R.setZero();

    int rowIndex = 0;
    for (size_t i = 0; i < GroundStations.size(); i++){
        if GroundStations[i].isObserving{
            Eigen::Vector3d rStation = GroundStations[i].r_ECI;
            Eigen::Vector3d vStation = {-wE*rStation(1), wE*rStation(0), 0.0};

            H(rowIndex, 0) = -(0.5*(2.0*rStation(0) - 2.0*prior.x(0)))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 1) = -(0.5*(2.0*rStation(1) - 2.0*prior.x(1)))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 2) = -(0.5*(2.0*rStation(2) - 2.0*prior.x(2)))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            R(rowIndex, rowIndex) = GroundStations[i].RangeAcc/2.0;

            rowIndex = rowIndex + 1;

            H(rowIndex, 0) = (0.5*(2.0*rStation(0) - 2.0*prior.x(0))*((rStation(0) - prior.x(0))*(vStation(0) - prior.x(3)) + (rStation(1) - prior.x(1))*(vStation(1) - prior.x(4)) + (rStation(2) - prior.x(2))*(vStation(2) - prior.x(5))))/pow(((rStation(0) - prior.x(0))^2 + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),1.5) - (vStation(0) - prior.x(3))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 1) = (0.5*(2.0*rStation(1) - 2.0*prior.x(1))*((rStation(0) - prior.x(0))*(vStation(0) - prior.x(3)) + (rStation(1) - prior.x(1))*(vStation(1) - prior.x(4)) + (rStation(2) - prior.x(2))*(vStation(2) - prior.x(5))))/pow(((rStation(0) - prior.x(0))^2 + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),1.5) - (vStation(1) - prior.x(4))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 2) = (0.5*(2.0*rStation(2) - 2.0*prior.x(2))*((rStation(0) - prior.x(0))*(vStation(0) - prior.x(3)) + (rStation(1) - prior.x(1))*(vStation(1) - prior.x(4)) + (rStation(2) - prior.x(2))*(vStation(2) - prior.x(5))))/pow(((rStation(0) - prior.x(0))^2 + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),1.5) - (vStation(2) - prior.x(5))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 3) = -(rStation(0) - prior.x(0))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 4) = -(rStation(1) - prior.x(1))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            H(rowIndex, 5) = -(rStation(2) - prior.x(2))/pow((pow((rStation(0) - prior.x(0)),2) + pow((rStation(1) - prior.x(1)),2) + pow((rStation(2) - prior.x(2)),2)),0.5);
            R(rowIndex, rowIndex) = GroundStations[i].RangeRateAcc/2.0;

            rowIndex = rowIndex + 1;

        }
    }

    Eigen::Matrix<float,6,2*NObs> K = prior.P*H.transpose()*(H*prior.P*H.transpose() + R).inverse();
    Eigen::Matrix<float,6,1> xPosterior = prior.x + K*y;
    Eigen::Matrix<float,6,6> PPosterior = (Eigen::Matrix<float,6,6>::Identity() - K*H)*prior.P*(Eigen::Matrix<float,6,6>::Identity() - K*H).transpose() + K*R*K.transpose();
    return KalmanEstimate(xPosterior, PPosterior);

}