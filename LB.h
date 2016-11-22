/* 
 * File:   LB.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:06 PM
 */

#ifndef LB_H
#define LB_H

//class IO;

#include <stdio.h>
//#include <iostream>
//#include <vector>
#include <stdlib.h>
#include <math.h>
//#include <algorithm>
#include <string.h>
//
#include "macros.h"
#include "vector.h"
#include "lattice.h"
#include "node.h"
#include "utils.h"
#include "elmt.h"
#include "getpot.h"

using namespace std;
extern ProblemName problemName;

class LB{
    // Characteristics ///////////////////////////////
    // lattice-Boltzmann equation is solved in a dimensionless space
    // with lattice spacing and time step both unitary
public: //private

    // initial viscosity
    double initDynVisc;
    // parameters for non-Newtonian fluid
    double yieldStress, plasticVisc;
    // parameter for Smagorinsky turbulence model
    double turbConst;
    // slip coefficient
    double slipCoefficient;
    // initial density of the fluid
    double initDensity;
    // initial velocity for the fluid
    tVect initVelocity;
    // force field (if not constant need to be calculated inside cycle)
    tVect lbF;
    // total number of nodes
    unsigned int totNodes;
    // total mass (initial)
    double totalMass;
    // standard neighbors shifting
    intList ne;
    unsIntList shift;
    intList domain;
    // node types (boundary are located halfway between nodes)
    typeList types;
    vecList positions;
    // curved info container
    curveList curves;
    // neighbor list for every node
    neighborList neighbors;
    // boundary conditions
    unsIntList boundary;
    // list with nodes
    nodeList nodes;
    // indexes of active nodes
    unsIntList activeNodes;
    // indexes of fluid nodes
    unsIntList fluidNodes;
    // indexes of solid nodes
    unsIntList particleNodes;
    // indexes of interface nodes
    unsIntList interfaceNodes;
public:
    // switchers for force field, non-Newtonian and everything
    bool freeSurface;
    bool forceField;
    bool nonNewtonian;
    bool turbulenceOn;
    // number of LBM step before activating the free surface
    unsigned int lbmInitialRepeat;
    // absolute time
    unsigned int time;
    unsIntList lbSize;
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    measureUnits unit;
    // problem-specific stuff ///////////////////////////////
    // stuff for shear cell
    double maxVisc, maxPlasticVisc, maxYieldStress, maxShearRate;
    unsigned int viscSteps, shearRateSteps;
    // stuff for drum
    double fluidMass;
    // stuff for avalanches (net-like)
    double avalanchePosit;
    // standard constructor
public:
    LB(){
        //
        freeSurface=false;
        nonNewtonian=false;
        forceField=false;
        turbulenceOn=false;
        //
        time=0;
        lbSize.resize(lbmDim);
        lbSize[0]=lbSize[1]=lbSize[2]=1;
        //
        yieldStress=1.0;
        plasticVisc=1.0;
        turbConst=0.0;
        slipCoefficient=0.0;
        //
        initDynVisc=(1.0-0.5)/3.0/lbmDt;
        initDensity=0.0;
        initVelocity.reset();
        lbF.reset();
        //
        totNodes=1;
        //
        shift.resize(lbmDim);
        shift[0]=shift[1]=shift[2]=1;
        domain.resize(lbmDim);
        domain[0]=domain[1]=domain[2]=1;
        //
        boundary.resize(2*lbmDim);
        ne.resize(lbmDirec);
        // initialize variable containers
        nodes.clear();
        types.clear();
        neighbors.clear();
        // initialize lists
        activeNodes.clear();
        particleNodes.clear();
        interfaceNodes.clear();
        //
        lbmInitialRepeat=0;
        // stuff for shear cell
        maxVisc=0.0;
        maxPlasticVisc=0.0;
        maxYieldStress=0.0;
        maxShearRate=0.0;
        viscSteps=0;
        shearRateSteps=0;
        // stuff for drum
        fluidMass=0.0;
        }
    // showing Lattice characteristics
    void LBShow() const;
    void latticeDefinition();
    void latticeBoltzmannGet(GetPot& lbmCfgFile, GetPot& command_line);
    void latticeBolzmannInit(cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects);
    void latticeBolzmannStep(elmtList& elmts, particleList& particles, wallList& walls);
    void latticeBoltzmannCouplingStep(bool& newNeighborList, elmtList& eltms, particleList& particles);
    void latticeBoltzmannFreeSurfaceStep();
    void evolveBoundary();
    void accelerateDrum(double& drumSpeed, cylinderList& cylinders, wallList& walls);
private:
    // initialize neighbors of a node coherently with the lattice
    void initializeNodeNeighbors(node& dummyNode);
    // initialization functions
    void initializeTypes(particleList& particles, wallList& walls, cylinderList& cylinders, objectList& objects);
    void initializeNodes();
    void initializeLatticeBoundaries();
    void initializeParticleBoundaries(particleList& particles);
    void initializeWallBoundaries(wallList& walls);
    void initializeObjectBoundaries(objectList& objects);
    void initializeCylinderBoundaries(cylinderList& cylinders);
    void initializeCurved(cylinderList& cylinders);
    void initializeInterface(double totParticles);
    void initializeLists();
    void resetLists();
    void initializeVariables();
    void initializeWalls(wallList& walls, cylinderList& cylinders, objectList& objects);
    // integration functions
    void reconstruction();
    void collision();
    void streaming(elmtList& elmts, particleList& particles, wallList& walls);
    // error dumping functions
    void cleanLists();
    // coupling functions
    void updateLists(unsIntList& newPopUpNodes, unsIntList& newSolidNodes);
    void computeHydroForces(elmtList& elmts, particleList& particles);
    void findNewActive(unsIntList& newPopUpNodes, elmtList& elmts, particleList& particles);
    void findNewSolid(unsIntList& newSolidNodes, elmtList& elmts, particleList& particles);
    void activeToSolid(unsIntList& newSolidNodes, elmtList& elmts, double& massSurplus);
    void solidToActive(unsIntList& newPopUpNodes, elmtList& elmts, double& massSurplus);
    void updateBoundary(unsIntList& newPopUpNodes, unsIntList& newSolidNodes, elmtList& elmts);
    void updateIndices(elmtList& elmts, particleList& particles);
//    void updateBoundaryOld(intList& newPopUpNodes, intList& newSolidNodes, particle& dummyParticle);
    // interface functions
    void updateMass();
    void updateInterface();
    void findInterfaceMutants(unsIntList& filledNodes, unsIntList& emptiedNodes);
    void smoothenInterface(unsIntList& newFluidNodes, unsIntList& newGasNodes, unsIntList& newInterfaceNodes, double& massSurplus);
    void updateMutants(unsIntList& filledNodes, unsIntList& emptiedNodes, unsIntList& newInterfaceNodes, double& massSurplus);
    void removeIsolated(double& massSurplus);
    void redistributeMass(const double& massSurplus);
    void enforceMassConservation();
    //boundary evolution function
    void gassification(unsIntList& gassificationNodes);
public:
    // function for linearized index management
    tVect setPosition(const unsigned int& index) const;
    unsigned int getX(const unsigned int& index) const;
    unsigned int getY(const unsigned int& index) const;
    unsigned int getZ(const unsigned int& index) const;

};

#endif	/* LB_H */