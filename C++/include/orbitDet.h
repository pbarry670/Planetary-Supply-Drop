#ifndef ORBITDET_H
#define ORBITDET_H

#include "orbits.h"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <vector>


const float visibilityConeAngle = 100.0;

struct GroundStation {

    std::string locationName;

    float lat;
    float lon;
    float alt;

    float rangeAcc;
    float rangeRateAcc;

    Eigen::Vector3d r_ECEF;
    Eigen::Vector3d r_ECI;

    bool isObserving = false;
    GroundStation(std::string locationName, float lat, float lon, float alt, float rangeAcc, float rangeRateAcc, float alpha);

};

struct KalmanEstimate {

    Eigen::Matrix<float,6,1> x;
    Eigen::Matrix<float,6,6> P;

    KalmanEstimate(Eigen::Matrix<float,6,1> x, Eigen::Matrix<float,6,6> P);

};

std::vector<GroundStation> readInputGroundStationFile(std::string groundStationInputFilename, float alpha);
KalmanEstimate timeUpdate(KalmanEstimate posterior, Eigen::Matrix<float,6,6> Q, float mu, float dt);
std::vector<GroundStation> updateGroundStationLocations(std::vector<GroundStation> GroundStations, float alpha);
std::vector<GroundStation> updateGroundStationObservability(std::vector<GroundStation> GroundStations, Eigen::Matrix<float,6,1> xSatECI);
int getNumberObservingGroundStations(std::vector<GroundStation> GroundStations);
Eigen::Matrix<float,Eigen::Dynamic,1> getMeasurement(Eigen::Matrix<float,6,1> xSatECI, std::vector<GroundStation> GroundStations, float wE, int NObs);
Eigen::Matrix<float,Eigen::Dynamic,1> measurementModel(Eigen::Matrix<float,6,1> xSatECI, std::vector<GroundStation> GroundStations, float wE, int NObs);
KalmanEstimate measurementUpdate(KalmanEstimate prior, Eigen::Matrix<float, Eigen::Dynamic, 1> Y, std::vector<GroundStation> GroundStations, float wE, int NObs);

#endif
