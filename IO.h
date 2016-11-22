/*
 * File:   IO.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:06 PM
 */

#ifndef IO_H
#define	IO_H


#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <ctime>
#include <csignal>
#include <iomanip>
// only to create a directory
#include <sys/stat.h>
#include <sys/types.h>
//
#include "LB.h"
#include "DEM.h"
#include "vector.h"
#include "elmt.h"
#include "node.h"

using namespace std;
extern ProblemName problemName;

// define precision outputs for statistic file
#define wSt int(15)

class IO{
public:
    // problem name
    string problemNameString;
    // time at invocation
    struct tm *now;
    // config file
    string lbmCfgName;
    // time integration
    unsigned int currentTimeStep, maximumTimeSteps;
    // maximum time variable
    double maxTime;
    // intervals for output
    double screenExpTime, stateExpTime, fluidExpTime, partExpTime, objectExpTime, recycleExpTime, outputExpTime;
    // support variables for output
    unsigned int lastScreenExp, lastFluidExp, lastPartExp, lastObjectExp, lastOutputExp, lastRecycleExp;
    // formats for output file names
    string fluidFileFormat, partFileFormat, objectFileFormat, recycleFileFormat;
    // support variable to store file names
    char filePathBuffer [1024];
    // boolean for exit control
    bool problem;
    // switchers for solvers
    bool demSolve, lbmSolve;
    // simulation time in real units
    double realTime;
    // locations
    string workDirectory, partDirectory, fluidDirectory;
    // output file names
    string exportFileName, maxSpeedFileName, plasticityFileName;
    string forceFileName, viscFileName;
    // generic file output streams
    ofstream exportFile, maxSpeedFile, plasticityFile;
    ofstream forceFile, viscFile;
    // drum files
    string centerOfMassFileName, positionFileName;//, fluidPositionFileName, particlePositionFileName; // particleDensityFileName
    ofstream centerOfMassFile, positionFile;//, fluidPositionFile, particlePositionFile; // particleDensityFile
    // new drum files
    string profileAngleFileName, fluidLevelFileName, particleLevelFileName, pressureFileName, fluidSpeedFileName, particleSpeedFileName, particleSurfaceFileName, fluidSurfaceFileName;
    ofstream profileAngleFile, fluidLevelFile, particleLevelFile, pressureFile, fluidSpeedFile, particleSpeedFile, particleSurfaceFile, fluidSurfaceFile;
    string particleAngleFileName, particleNumberFileName, collisionNumberFileName, hydraulicNumberFileName, collisionForceFileName, hydraulicForceFileName;
    ofstream particleAngleFile, particleNumberFile, collisionNumberFile, hydraulicNumberFile, collisionForceFile, hydraulicForceFile;
    // avalanche files
    string obstacleFileName;
    ofstream obstacleFile;

    // time-averaged quantites
    double msdSum, msdSum2;
    unsigned int msdSumN;
    tVect msdVecSum, msdVecSum2;
    //
    doubleList thetaList, collisionForceAngle, hydraulicForceAngle;
    intList particleNumberAngle, collisionNumberAngle, hydraulicNumberAngle;
    int timeStepCount;
    
    // new numbering for recycle particle (e.g. recyclePart_1, recyclePart_2 to ....) (Devis)
    unsigned int recNum;        

    // stuff for energy file(Devis)
    string energyFileName;
    ofstream energyFile;
	
    // variable for force statistic (.fstat) (Devis)
    string dataDirectory; 
    string statDirectory;      
    string restartDirectory;
    string dataFileName;
    string statisticFileName;
    string restartFileName;
    ofstream restartFile;
    ofstream dataFile;
    ofstream statParticleFile;
	
    int saveCount;
    
public:
    void initialize();
    void outputStep(const double& demTime, const LB& lb, const DEM& dem);
    void outputFinal();
    void initViscFile();
    void writeViscFile(const LB& lb, const wallList& walls, const elmtList& elmts);
    void exportRecycleParticles(const elmtList& elmts, const pbcList& pbcs);
private:
    // paraview
    void createParaviewFiles(const LB& lb, const DEM& dem);
    void exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile);
    void exportParaviewObjects(const objectList& particles, const string& objectFile);
    void exportParaviewBox(const double& time);
    void exportParaviewFluidOld(const LB& lb, const string& fluidFile);
    void exportParaviewFluidNew(const LB& lb);
    // print export stuff
    void exportMaxSpeedFluid(const LB& lb);
    void exportMaxSpeedParticles(const DEM& dem);
    void exportTotalMass(const LB& lb);
    void exportPlasticity(const LB& lb);
    void exportShearCell(const LB& lb, const DEM& dem);
    void exportForces(const DEM& dem);
    void exportForceObstacle(const DEM& dem);
    void exportDrum(const LB& lb, const DEM& dem);
    void exportEnergy(const DEM& dem);
    // data elaboration
    double totPlastic(const LB& lb) const;
    double totFluidMass(const LB& lb) const;
    double totParticleMass(const elmtList& elmts) const;
    tVect  centerOfMassFluid(const LB& lb) const;
    tVect centerOfMassParticles(const elmtList& elmts) const;
    tVect totForceObject(const objectList& objects) const;
    void  apparentViscosity(const LB& lb, const wallList& walls, double& externalShear, double& wallStress, double& appVisc) const;
    double hydraulicForceTot(const elmtList& elmts) const;
    double collisionForceTot(const elmtList& elmts) const;
    double totVolumeParticles(const LB& lb) const;
    // drum specific
    void drumMassCenter(const LB& lb, const elmtList& elmts, const double& particleDensity, const cylinder& drumCylinder);
    void drumPositions(const LB& lb, const elmtList& elmts, const cylinder& drumCylinder);
//    void drumSurfaceProperties(const LB& lb, const elmtList& elmts, const cylinder& drum);
//    void drumBaseProperties(const LB& lb, const elmtList& elmts, const cylinder& drum);
    void drumProfileProperties(const LB& lb, const elmtList& elmts, const cylinder& drum);
    void drumParticleProperties(const DEM& dem);
    void updateDistributions(const DEM& dem, const cylinder& drum);
    // energy (devis)
    void exportEnergySystem(const DEM& dem);
};


#endif	/* IO_H */

