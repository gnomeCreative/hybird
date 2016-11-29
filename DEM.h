/*
 * File:   DEM.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:07 PM
 */

#ifndef DEM_H
#define	DEM_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <fstream>
//
#include "vector.h"
#include "elmt.h"
#include "utils.h"
#include "LB.h"

using namespace std;
extern ProblemName problemName;

// forward declaration of IO
class IO;

class DEM{
    
    //IO *io;
    
private:
    // domain size: DEM grid is orthogonal
    doubleList demSize;
    // neighbor list (here COMMENTS are needed)
    // neighbor parameters
    double maxDisp;
    double nebrRange;
    double nebrShell;
    double cellWidth[3];
    unsigned int nCells[3];
    // neighbor tables
    intList cellTable;
    unsIntList neighborTable;
    //unsIntList nearWallTable;
    //unsIntList nearObjectTable;
    //unsIntList nearCylinderTable;
    unsIntList newNearCylinderTable;
    // capping distance to avoid infinite lubrication forces
    double minDistLub;
    // numerical viscosity to avoid energy problems
    double numVisc;
    
    // prototype shape collection
    std::vector <vecList> prototypes;
    
    typedef std::vector<unsIntList> matrix;
    matrix nearWallTable;
    matrix nearObjectTable;
    matrix nearCylinderTable;
public:
    // energy
    energy particleEnergy;
    // number of dem time steps before activating LB
    double demInitialRepeat;
    // time step
    double deltat;
    // force field (if not constant need to be calculated inside cycle) // should not be here!!!!!
    tVect demF;
    // material for elements // should not be here!!!
    material sphereMat;
    // time
    double demTime;
    // relative difference between time step in DEM and time step in LBM
    unsigned int multiStep;
    // ratio between time step and estimated duration of contacts
    double criticalRatio;
    // solid walls
    wallList walls;
    // solid walls
    cylinderList cylinders;
    // periodic boundary conditions
    pbcList pbcs;
    // list with elements
    elmtList elmts;
    // list with particles
    particleList particles;
    // list with objects
    objectList objects;
    //couple elongation list
    map<string, Elongation> elongTable;
    // particles contains both standard particles and ghost particles. The total number of standard particles needs to be saved
    unsigned int stdParticles;
    // list with ghost elements
    ghostList ghosts;
    // true if indices for the LB need to be regenerated
    bool newNeighborList;
    // stuff for drum
    double drumSpeed;
    // Total number of standard objects (i.e. not ghosts) which needs to be saved (Devis))
    unsigned int stdObjects;
    
public:
    DEM(){
        demInitialRepeat=0.0;
        newNeighborList=false;
        demSize.resize(3);
        demSize[0]=demSize[1]=demSize[2]=1.0;
        demTime=0.0;
        deltat=1.0;
        demF.reset();
        walls.clear();
        pbcs.clear();
        elmts.clear();
        particles.clear();
        objects.clear();
        stdParticles=0;
        ghosts.clear();
//        ghosts.clear();
        cellTable.clear();
        neighborTable.clear();
        nearWallTable.clear();
        nearWallTable.clear();
        nearObjectTable.clear();
        nearCylinderTable.clear();
        newNearCylinderTable.clear();
        // neighbor list variables
        maxDisp=0.0;
        nebrRange=0.0;
        cellWidth[0]=cellWidth[0]=cellWidth[0]=0.0;
        nCells[0]=nCells[0]=nCells[0]=0;
        prototypes.clear();
//        energy.reset();
        criticalRatio=0.1;
        multiStep=1;
        numVisc=0.0;
        //
        drumSpeed=0.0;
        
}
    void discreteElementStep(IO& io);
    void discreteElementGet(GetPot& lbmCfgFile, GetPot& command_line);
    void discreteElementInit(unsIntList& boundary, unsIntList& lbSize, measureUnits& unit, tVect& lbF);
    void evolveBoundaries();
    // devis
    void exportDataFile(IO& io);
    void exportRestartFile(IO& io);
    void printHeaderStatFile(IO& io);
private:
    // initialization functions
    void compositeProperties();
    void initializeWalls(unsIntList& boundary, unsIntList& lbSize, measureUnits& unit);
    void initializeCylinders();
    void initializePbcs(unsIntList& boundary, unsIntList& lbSize, measureUnits& unit);
    // integration functions
    void predictor();
    void corrector();
    void evaluateForces(IO& io);
    void updateParticlesPredicted();
    void updateParticlesCorrected();
    double criticalTimeStep() const;
    // neighbor list functions
    void evalNeighborTable(); // should be private: correct!
    void initNeighborParameters();
    void evalMaxDisp();
    void evalCellTable();
    void evalNearWallTable();
    void evalNearObjectTable();
    void evalNearCylinderTable();
    // periodicity functions
    void createGhosts();
    void pbcShift();
    void periodicObjects();
    // force computation functions
    void particleParticleContacts(IO& io);
    void wallParticleContacts(IO& io);
    void cylinderParticelContacts();
    void objectParticleContacts(IO& io);
    
    inline void particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance, IO& io, Elongation& elongation);
    inline void wallParticleCollision(wall *walli, const particle *partJ, const double& overlap, IO& io, Elongation& elongation);
    inline void cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap, Elongation& elongation);
    inline void objectParticleCollision(object *iObject, const particle *partJ, const tVect& vectorDistance, IO& io, Elongation& elongation);
    
    
    //inline void particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance, IO& io);
    //inline void wallParticleCollision(wall *walli, const particle *partJ, const double& overlap, IO& io);
    //inline void cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap);
    //inline void objectParticleCollision(object *iObject, const particle *partJ, const tVect& vectorDistance, IO& io);
    double normalContact(const double& overlap, const double& vrelnnorm, const double& rEff, const double& massEff) const;
    double tangentialContact(const double& vreltNorm, const double& fn, const double& effRad, const double& effMass, const double& friction) const;
    
    tVect FRtangentialContact(const tVect& tangRelVelContact, const double& fn, 
             const double& effMass, Elongation& elongation,const double& friction, const double& tangStiff, const double& viscTang);
    
    void lubrication(int i, int j, double ri, double rj, tVect x0ij);
    // energy functions
    void updateEnergy();
    /// this are the (new) ones that save the contacts (Devis))
    void saveStatPPcollision(IO& io,const particle *partI, const particle *partJ, const double& overlap, const double& normNormalForce, const tVect& en, const double normTangForce, const tVect et) const;
    void saveStatWPcollision(IO& io,const wall *wallI, const particle *partJ, const double& overlap, const double& normNormalForce, const tVect& en, const double normTangForce, const tVect et) const ;
    void saveStatOPcollision(IO& io, const object *ObjectI, const particle *partJ, const double& overlap, const double& normNormalForce, const tVect& en, const double normTangForce, const tVect et) const;
    
    void deleteParticleElongation();
    void deleteWallElongation();
    tVect rollingContact(const tVect& wI,const tVect& wJ,const double effRad,const double& fn, 
             const double& rolling);
};

#endif	/* DEM_H */

