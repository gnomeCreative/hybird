
#include "LB.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PUBLIC FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

void LB::LBShow() const {
    cout<<"LATTICE CHARACTERISTICS (lattice units):\n";
    cout<<"directions="<<lbmDirec<<";\n";
    cout<<"internal time step"<<lbmDt<<";\n";
    cout<<"time unit = "<<unit.Time<<"; length unit = "<<unit.Length<<"; density unit = "<<unit.Density<<";\n";
    cout<<"dynVisc unit = "<<unit.DynVisc<<"; kinVisc unit = "<<unit.KinVisc<<"; speed unit = "<<unit.Speed<<";\n";
    cout<<"xdim ="<<lbSize[0]<<"; ydim= "<<lbSize[1]<<"; zdim= "<<lbSize[2]<<";\n";
    cout<<"Init tau = "<<(0.5+3.0*initDynVisc)<<";init visc ="<<initDynVisc<<"; init density = "<<initDensity<<";\n";
    cout<<"Min tau = "<<minTau<<"; Max tau = "<<maxTau<<";\n";
    if (nonNewtonian) {
        cout<<"Plastic viscosity="<<plasticVisc<<"(tau = "<<(0.5+3.0*plasticVisc)<<" ); yield stress="<<yieldStress<<";\n";
    }
    cout<<"Initial Velocity=";
    initVelocity.show();
    cout<<"; F=";
    lbF.show();
    cout<<";\n";
    cout<<"X(-)="<<boundary[0]<<"; X(+)="<<boundary[1]<<"; Y(-)="<<boundary[2]<<"; Y(+)="<<boundary[3]<<"; Z(-)="<<boundary[4]<<"; Z(+)="<<boundary[5]<<";\n";
    cout<<"PHISICAL CHARACTERISTICS (physical units)\n";
    cout<<"xdim="<<(lbSize[0]-2)*unit.Length<<"; ydim= "<<(lbSize[1]-2)*unit.Length<<"; zdim= "<<(lbSize[2]-2)*unit.Length<<";\n";
    cout<<"init kin viscosity="<<initDynVisc*unit.KinVisc<<", init dyn viscosity="<<initDynVisc*unit.DynVisc<<"; init density="<<initDensity*unit.Density<<";\n";
    if (nonNewtonian) {
        cout<<"plastic kin viscosity="<<plasticVisc*unit.KinVisc<<", plastic dyn viscosity="<<plasticVisc*unit.DynVisc<<"; yield stress="<<yieldStress*unit.Stress<<";\n";
    }
    cout<<"Min dyn viscosity = "<<(minTau-0.5)/3.0*unit.DynVisc<<"; Max dym viscosity = "<<(maxTau-0.5)/3.0*unit.DynVisc<<";\n";
    cout<<"Initial Velocity=";
    tVect initVelocityPhysic=initVelocity*unit.Speed;
    initVelocityPhysic.show();
    cout<<"; F=";
    tVect FPhysic=lbF*unit.Accel;
    FPhysic.show();
    cout<<";"<<endl;
}

void LB::latticeDefinition() {
    // LATTICE PARAMETERS  ////////////////
    //size-dependent; the others are in the lattice.h
    
    // domain movement variables
    shift[2]=lbSize[0]*lbSize[1];
    shift[1]=lbSize[0];
    shift[0]=1;
    //
    domain[2]=lbSize[2]*lbSize[1]*lbSize[0]-2*shift[2];
    domain[1]=lbSize[1]*lbSize[0]-2*shift[1];
    domain[0]=lbSize[0]-2*shift[0];
    
    // standard neighbors shifting
    // order is O,x,y,z,xy,yz,zx.
    ne[0]=0;
    //
    ne[1]=shift[0];
    ne[2]=-shift[0];
    //
    ne[3]=shift[1];
    ne[4]=-shift[1];
    //
    ne[5]=shift[2];
    ne[6]=-shift[2];
    //
    ne[7]=shift[0]+shift[1];
    ne[8]=-shift[0]-shift[1];
    ne[9]=-shift[0]+shift[1];
    ne[10]=shift[0]-shift[1];
    //
    ne[11]=shift[1]+shift[2];
    ne[12]=-shift[1]-shift[2];
    ne[13]=-shift[1]+shift[2];
    ne[14]=shift[1]-shift[2];
    //
    ne[15]=shift[2]+shift[0];
    ne[16]=-shift[2]-shift[0];
    ne[17]=-shift[2]+shift[0];
    ne[18]=shift[2]-shift[0];
    
}

void LB::latticeBoltzmannGet(GetPot& lbmCfgFile, GetPot& command_line) {
    
    // conversion units //////////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // measure units for DEM solver are in the international system
    
    // primary
    PARSE_CLASS_MEMBER(lbmCfgFile, unit.Length, "unitLength",1.0);
    ASSERT(unit.Length > 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, unit.Time, "unitTime",1.0);
    //ASSERT(unit.Time > 0);
    PARSE_CLASS_MEMBER(lbmCfgFile, unit.Density, "unitDensity",1.0);
    ASSERT(unit.Density > 0);
    
    // fixing length: domain size
    double lbSizeX, lbSizeY, lbSizeZ;
    lbSize.resize(lbmDim);
    PARSE_CLASS_MEMBER(lbmCfgFile, lbSizeX, "lbSizeX",0.0);
    ASSERT(lbSizeX > 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, lbSizeY, "lbSizeY",0.0);
    ASSERT(lbSizeY > 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, lbSizeZ, "lbSizeZ",0.0);
    ASSERT(lbSizeZ > 0.0);
    // scaling and adapting length unit to get exact domain discretization
    lbSize[0]=int(floor(lbSizeX/unit.Length+0.5));
    lbSize[1]=int(floor(lbSizeY/unit.Length+0.5));
    lbSize[2]=int(floor(lbSizeZ/unit.Length+0.5));
    
    // getting viscosity
    PARSE_CLASS_MEMBER(lbmCfgFile, plasticVisc, "plasticVisc",1.0);
    ASSERT(plasticVisc >= 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, yieldStress, "yieldStress",1.0);
    ASSERT(yieldStress >= 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, initDynVisc, "initVisc",1.0);
    ASSERT(initDynVisc > 0.0);
    
    // if time step is chosen as optimal, compute it
    // see thesis, §4.1.1
    if (unit.Time==0.0) {
        const double tauOpt=1.0;
        if (nonNewtonian) {
            unit.Time=unit.Density*unit.Length*unit.Length/plasticVisc*(1.0/3.0)*(tauOpt-0.5);
        }
        else {
            unit.Time=unit.Density*unit.Length*unit.Length/initDynVisc*(1.0/3.0)*(tauOpt-0.5);
        }
        cout<<"Time step automatically chosen: "<<unit.Time<<" s"<<endl;
    }
    
    
    // compute all non-primary conversion units
    unit.setComposite();
    
    // scaling viscosity
    initDynVisc/=unit.DynVisc;
    plasticVisc/=unit.DynVisc;
    yieldStress/=unit.Stress;
    
    //    PARSE_CLASS_MEMBER(lbmCfgFile, initDensity, "initDensity",1.0);
    //    initDensity/=unit.Density;
    // to avoid errors we set this to be just 1
    initDensity=1.0;
    ASSERT(initDensity > 0.0);
    
    double initVelocityX, initVelocityY, initVelocityZ;
    PARSE_CLASS_MEMBER(lbmCfgFile, initVelocityX, "initVelocityX",0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, initVelocityY, "initVelocityY",0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, initVelocityZ, "initVelocityZ",0.0);
    initVelocity=tVect(initVelocityX, initVelocityY, initVelocityZ);
    initVelocity/=unit.Speed;
    
    double lbFX, lbFY, lbFZ;
    PARSE_CLASS_MEMBER(lbmCfgFile, lbFX, "lbFX",0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, lbFY, "lbFY",0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, lbFZ, "lbFZ",0.0);
    lbF=tVect(lbFX, lbFY, lbFZ);
    switch (problemName) {
        case demChute: {
            double chuteInclination=0.0;
            PARSE_CLASS_MEMBER(lbmCfgFile, chuteInclination, "chuteInclination",0.0);
            const double  gravity=9.086;
            lbFX = -1.0*gravity*sin(chuteInclination*(M_PI/180));
            lbFY = -1.0*0.0;
            lbFZ = -1.0*gravity*cos(chuteInclination*(M_PI/180));
            lbF=tVect(lbFX, lbFY, lbFZ);
        }
    }
    lbF/=unit.Accel;   
    
    boundary.resize(lbmDim);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[0], "boundary0",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[1], "boundary1",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[2], "boundary2",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[3], "boundary3",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[4], "boundary4",1);
    PARSE_CLASS_MEMBER(lbmCfgFile, boundary[5], "boundary5",1);
    
    PARSE_CLASS_MEMBER(lbmCfgFile, slipCoefficient, "slipCoefficient",0.0);
    ASSERT(slipCoefficient >= 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, turbConst, "turbConst",0.0);
    ASSERT(turbConst >= 0.0);
    
}

void LB::latticeBolzmannInit(cylinderList& cylinders, wallList& walls, particleList& particles, objectList& objects) {
    //  Lattice Boltzmann initialization steps
    
    // first comes the initialization of the data structures
    initializeNodes();
    // then the initial node type must be identified for every node
    initializeTypes(particles, walls, cylinders,objects);
    // initialize node properties and lists
    initializeLists();
    initializeVariables();
    initializeWalls(walls, cylinders,objects);
    //create active list and checks for errors in the list
    cleanLists();
    // calculate totalMass
    totalMass=0.0;
    switch (problemName) {
        case DRUM: {
            totalMass=fluidMass/unit.Mass;
            break;
        }
        default: {
            for (int it=0; it<nodes.size(); ++it){
                if (types[it].isActive() && !types[it].isInsideParticle()) {
                    totalMass+=nodes[it]->mass;
                }
            }
            break;
        }
    }
}

void LB::latticeBolzmannStep(elmtList& elmts, particleList& particles, wallList& walls) {
    // Lattice Boltzmann core steps
    //create active list and checks for errors in the list
    cleanLists();
    // reconstruction of macroscopic variables
    reconstruction();
    // compute interaction forces
    computeHydroForces(elmts, particles);
    //collision operator
    collision();
    //streaming operator
    streaming(elmts, particles, walls);
}

void LB::latticeBoltzmannFreeSurfaceStep() {
    // mass and free surface update
    updateMass();
    updateInterface();
    switch (problemName) {
        case DRUM: {
            enforceMassConservation();
            break;
        }
    }
}

void LB::latticeBoltzmannCouplingStep(bool& newNeighborList, elmtList& elmts, particleList& particles) {
    // identifies which nodes need to have an update due to particle movement
    // the complexity arises from trying to keep the scaling as close to linear as possible
    // maybe the best idea is to do this in two step:
    // 1) the first is to check for new active nodes and initialize them
    // 2) the second is checking for new solid nodes.
    // this automatically check also for the hideous case of particle to particle double transition
    
    // first we check if a new neighbor table has been defined. In that case, the indexing needs to be reinitialized
    if (newNeighborList) {
        updateIndices(elmts, particles);
        newNeighborList=false;
    }
    // declaring list for mutant nodes
    unsIntList newPopUpNodes, newSolidNodes;
    // emptying list of mutant nodes
    newPopUpNodes.clear();
    newSolidNodes.clear();
    // double massSurplus=0.0;
    
    // SOLID TO ACTIVE CHECK
    findNewActive(newPopUpNodes, elmts, particles);
    
    //    solidToActive(newPopUpNodes, elmts, massSurplus);
    
    // ACTIVE TO SOLID CHECK
    findNewSolid(newSolidNodes, elmts, particles);
    
    //    activeToSolid(newSolidNodes, elmts, massSurplus);
    
    // redistributing extra mass due to popping of nodes
    // redistributeMass(massSurplus);
    
}

void LB::accelerateDrum(double& drumSpeed, cylinderList& cylinders, wallList& walls) {
    
    static double accelerationTime=0.5;
    
    double currentTime=double(time)*unit.Time;
    if ((currentTime>=0.0)&&(currentTime<=accelerationTime)) {
        
        double currentSpeed=std::max(drumSpeed*currentTime/accelerationTime,drumSpeed);
        cout<<"Current speed = "<<currentSpeed<<"\n";
        
        // updating the geometrical elements (for DEM and coupling)
        // walls
        walls[2].omega=tVect(0.0,currentSpeed,0.0);
        walls[3].omega=tVect(0.0,currentSpeed,0.0);
        // drum
        cylinders[0].omega=tVect(0.0,currentSpeed,0.0);
        
        // updating solid nodes velocity
        for (int n=0; n<totNodes; ++n) {
            if (types[n].isDynWall()||types[n].isSlipDynWall()||types[n].isCurvedWall()) {
                // updating velocity for cylinder nodes
                if (!positions[n].insideCylinder(cylinders[0].p1/unit.Length, cylinders[0].naxes, cylinders[0].R/unit.Length)) {
                    nodes[n]->u=cylinders[0].getSpeed(positions[n]*unit.Length)/unit.Speed;
                }
                // updating velocity for walls nodes
                else if (positions[n].insidePlane(walls[2].p/unit.Length, walls[2].n)) {
                    nodes[n]->u=walls[2].getSpeed(positions[n]*unit.Length)/unit.Speed;
                }
                else if (positions[n].insidePlane(walls[3].p/unit.Length, walls[3].n)) {
                    nodes[n]->u=walls[3].getSpeed(positions[n]*unit.Length)/unit.Speed;
                }
            }
        }
    }
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PRIVATE FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

// initialization functions

void LB::initializeTypes(particleList& particles, wallList& walls, cylinderList& cylinders, objectList& objects) {
    
    // application of lattice boundaries
    initializeLatticeBoundaries();
    // application of particle initial position
    initializeParticleBoundaries(particles);
    // application of solid walls
    initializeWallBoundaries(walls);
    // application of solid cylinders
    initializeCylinderBoundaries(cylinders);
    // application of objects
    initializeObjectBoundaries(objects);
    // initializing curved properties
    initializeCurved(cylinders);
    // initialize interface
    initializeInterface(particles.size());
}

void LB::initializeNodes() {
    
    cout<<"Initializing nodes"<<endl;
    
    // total number of nodes
    totNodes=lbSize[0]*lbSize[1]*lbSize[2];
    // vector with node pointers
    nodes.resize(totNodes);
    for (int it=0; it<totNodes; ++it) {
        nodes[it]=0;
    }
    // vector with types
    types.resize(totNodes);
    for (int it=0; it<totNodes; ++it) {
        types[it].setFluid();
    }
    // vector with positions
    positions.resize(totNodes);
    for (int it=0; it<totNodes; ++it) {
        positions[it]=setPosition(it);
    }
    // vector with curves
    curves.resize(totNodes);
    for (int it=0; it<totNodes; ++it) {
        curves[it]=0;
    }
}

void LB::initializeLatticeBoundaries() {
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)
    
    // NEIGHBORING ////////////////////////////////////////////////
    // assigning neighboring conditions (this needs to be done after applying boundary conditions)
    cout<<"Initializing neighbors"<<endl;
    neighbors.resize(totNodes);
#pragma omp parallel for
    for (int it=0; it<totNodes; ++it) {
        // runs through free cells and identifies neighboring cells. If neighbor cell is
        // a special cell (periodic) then the proper neighboring condition is applied
        
        // resetting neighbors and vectors
        for (int j=1; j<lbmDirec; ++j) {
            neighbors[it].d[j]=it+ne[j];
        }
    }
    
    // BOUNDARY CONDITIONS ///////////////////////////
    // solid boundary wins over all in corners, where more than 1 bc is defined
    // nodes on the boundary have no neighbors
    cout<<"Initializing boundaries"<<endl;
#pragma omp parallel for
    for (int it=0; it<totNodes; ++it) {
        if (getX(it)==0 && !types[it].isWall()) {
            types[it].setType(boundary[0]);
            for (int j=0; j<lbmDirec; ++j) {
                neighbors[it].d[j]=it;
            }
        }
        else if (getX(it)==lbSize[0]-1 && !types[it].isWall()) {
            types[it].setType(boundary[1]);
            for (int j=0; j<lbmDirec; ++j) {
                neighbors[it].d[j]=it;
            }
        }
        if (getY(it)==0 && !types[it].isWall()) {
            types[it].setType(boundary[2]);
            for (int j=0; j<lbmDirec; ++j) {
                neighbors[it].d[j]=it;
            }
        }
        else if (getY(it)==lbSize[1]-1 && !types[it].isWall()) {
            types[it].setType(boundary[3]);
            for (int j=0; j<lbmDirec; ++j) {
                neighbors[it].d[j]=it;
            }
            
        }
        if (getZ(it)==0 && !types[it].isWall()) {
            types[it].setType(boundary[4]);
            for (int j=0; j<lbmDirec; ++j) {
                neighbors[it].d[j]=it;
            }
        }
        else if (getZ(it)==lbSize[2]-1 && !types[it].isWall()) {
            types[it].setType(boundary[5]);
            for (int j=0; j<lbmDirec; ++j) {
                neighbors[it].d[j]=it;
            }
        }
    }
    
    // PERIODICITY ////////////////////////////////////////////////
    // assigning periodicity conditions (this needs to be done after applying boundary conditions)
    cout<<"Initializing periodicity"<<endl;
#pragma omp parallel for
    for (int it=0; it<totNodes; ++it) {
        // runs through free cells and identifies neighboring cells. If neighbor cell is
        // a special cell (periodic) then the proper neighboring condition is applied
        
        // neighboring and periodicity vector for boundary update
        int pbc[lbmDirec];
        for (int j=1; j<lbmDirec; ++j) {
            pbc[j]=0;
        }
        
        // calculates the effect of periodicity
        if (types[it].isFluid()) {
            tVect vx(1.0,0.0,0.0);
            tVect vy(0.0,1.0,0.0);
            tVect vz(0.0,0.0,1.0);
            for (int j=1; j<lbmDirec; ++j) {
                if (types[neighbors[it].d[1]].isPeriodic() && vx.dot(v[j])>0.0)
                    pbc[j]-=domain[0];
                if (types[neighbors[it].d[2]].isPeriodic() && vx.dot(v[j])<0.0)
                    pbc[j]+=domain[0];
                if (types[neighbors[it].d[3]].isPeriodic() && vy.dot(v[j])>0.0)
                    pbc[j]-=domain[1];
                if (types[neighbors[it].d[4]].isPeriodic() && vy.dot(v[j])<0.0)
                    pbc[j]+=domain[1];
                if (types[neighbors[it].d[5]].isPeriodic() && vz.dot(v[j])>0.0)
                    pbc[j]-=domain[2];
                if (types[neighbors[it].d[6]].isPeriodic() && vz.dot(v[j])<0.0)
                    pbc[j]+=domain[2];
            }
            // apply periodicity
            for (int j=1; j<lbmDirec; ++j) {
                neighbors[it].d[j]+=pbc[j];
            }
        }
    }
}

void LB::initializeParticleBoundaries(particleList& particles) {
    
    // INITIAL PARTICLE POSITION ////////////////////////
    cout<<"Initializing particle nodes"<<endl;
#pragma omp parallel for
    for (int it=0; it<totNodes; ++it) {
        // reset old list
        types[it].setOutsideParticle();
        // checking if we are not in a boundary
        // whose properties have already been assigned and must be invariant in time
        if (types[it].isActive()) {
            // checking if node is inside a particle
            for (int n=0; n<particles.size(); ++n) {
                if (positions[it].insideSphere(particles[n].x0/unit.Length,particles[n].r/unit.Length)) { //-0.5?
                    types[it].setInsideParticle();
                    types[it].setSolidIndex(particles[n].particleIndex);
                }
            }
        }
    }
}

void LB::initializeWallBoundaries(wallList& walls) {
    
    // SOLID WALLS ////////////////////////
    cout<<"Initializing solid walls"<<endl;
    for (int iw=0; iw<walls.size(); ++iw) {
        wall dummyWall;
        dummyWall=walls[iw];
        if (!dummyWall.flag) {
            for (int it=0; it<totNodes; ++it) {
                // creating solid cells
                if (positions[it].insidePlane(dummyWall.p/unit.Length, dummyWall.n)) {
                    // setting solidIndex
                    types[it].setSolidIndex(dummyWall.index);
                    // setting type: 5-6=slip, 7-8=no-slip
                    if (dummyWall.slip) {
                        // setting type for slip: 5=static, 6=moving
                        if (dummyWall.moving) {
                            types[it].setSlipDynWall();
                        }
                        else {
                            types[it].setSlipStatWall();
                        }
                    }
                    else{
                        // setting type for no-slip: 7=static, 8=moving
                        if (dummyWall.moving) {
                            types[it].setDynWall();
                        }
                        else {
                            types[it].setStatWall();
                        }
                    }
                }
            }
        }
    }
}

void LB::initializeObjectBoundaries(objectList& objects) {
    
    // SOLID WALLS ////////////////////////
    cout<<"Initializing objects"<<endl;
    for (int io=0; io<objects.size(); ++io) {
        object dummyObject;
        dummyObject=objects[io];
        for (int it=0; it<totNodes; ++it)  {
            const tVect nodePosition=positions[it];
            if (nodePosition.insideSphere(dummyObject.x0/unit.Length,dummyObject.r/unit.Length)) {// setting solidIndex
                types[it].setSolidIndex(dummyObject.index);
                types[it].setStatWall();
            }
        }
    }
}

void LB::initializeCylinderBoundaries(cylinderList& cylinders) {
    
    // SOLID CYLINDERS ////////////////////////
    cout<<"Initializing solid cylinders"<<endl;
    for (int ic=0; ic<cylinders.size(); ++ic) {
        cylinder dummyCylinder;
        dummyCylinder=cylinders[ic];
        for (int it=0; it<totNodes; ++it) {
            // creating solid cells
            if (!positions[it].insideCylinder(dummyCylinder.p1/unit.Length, dummyCylinder.naxes, dummyCylinder.R/unit.Length)) {
                // setting type: 5-6=slip, 7-8=no-slip
                if (dummyCylinder.slip) {
                    // setting type for slip: 5=static, 6=moving
                    if (dummyCylinder.moving) {
                        types[it].setSlipDynWall();
                    }
                    else {
                        types[it].setSlipStatWall();
                    }
                }
                else{
                    // setting type for no-slip: 7=static, 8=moving
                    if (dummyCylinder.moving) {
                        types[it].setDynWall();
                    }
                    else {
                        types[it].setStatWall();
                    }
                }
                // any case, initialize curved <- TO BE CHANGED!!!
                types[it].setCurvedWall();
                curves[it]=new curve;
            }
        }
    }
}

void LB::initializeCurved(cylinderList& cylinders) {
    cout<<"Initializing curved boundaries"<<endl;
    for (int it=0; it<totNodes; ++it) {
        if (types[it].isCurvedWall()) {
            for (int j=1; j<lbmDirec; ++j) {
                const unsigned int link=neighbors[it].d[j];
                if (!types[link].isWall()){
                    const tVect nodePos=unit.Length*positions[it];
                    curves[it]->delta[j]=1.0-cylinders[0].segmentIntercept(nodePos,unit.Length*v[j]);
                }
                curves[it]->computeCoefficients();
            }
        }
    }
}

void LB::initializeInterface(double totParticles) {
    // creates an interface electing interface cells from active cells
    // defining a surface
    
    cout<<"Initializing interface"<<endl;
    
    switch (problemName) {
        case demChute: {
            // PLANE
            tVect p,n;
            p=tVect(0.0, 0.0, 0.025/unit.Length+0.5);
            n=tVect(0.0, 0.0, 1.0);
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    tVect dist=positions[it]-p;
                    if (dist.dot(n)>0.0) {
                        types[it].setGas();
                    }
                }
            }
            break;
        }
        case SPLASH: {
            // PLANE
            tVect p,n;
            p=tVect(0.0, 0.0, 20.0);
            n=tVect(0.0, 0.0, 1.0);
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    tVect dist=positions[it]-p;
                    if (dist.dot(n)>0.0) {
                        types[it].setGas();
                    }
                }
            }
            break;
        }
        case AVALANCHE: {
            // avalanche
            double xMin=-10.0;
            double xMax=250.0;
            double yMin=-100.0;
            double yMax=10000.0;
            double zMin=-10.0;
            double zMax=20.0;
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    if ((getX(it)<xMin)||(getX(it)>xMax)||(getY(it)<yMin)||(getY(it)>yMax)||(getZ(it)<zMin)||(getZ(it)>zMax)) {
                        types[it].setGas();
                    }
                }
            }
            break;
        }
        case NET:
        case BARRIER:{
            // BOX BEFORE THE NET
            unsigned int xMin=int((avalanchePosit/unit.Length));
            unsigned int xMax=lbSize[0];
            unsigned int yMin=0;
            unsigned int yMax=lbSize[1];
            unsigned int zMin=0;
            unsigned int zMax=int((4.0/unit.Length));
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    if ((getX(it)<xMin)||(getX(it)>xMax)||(getY(it)<yMin)||(getY(it)>yMax)||(getZ(it)<zMin)||(getZ(it)>zMax)) {
                        types[it].setGas();
                    }
                }
            }
            // PLANE
            tVect p,n;
            p=tVect(int((avalanchePosit/unit.Length)), 0.0, 0.0);
            n=tVect(-1.0, 0.0, 1.0);
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    tVect dist=positions[it]-p;
                    if (dist.dot(n)>0.0) {
                        types[it].setGas();
                    }
                }
            }
            break;
        }
        case BOX: {
            double xMin=50.0;
            double xMax=125.0;
            double yMin=88.0;
            double yMax=162.0;
            double zMax=70.0;
            double zMin=-10.0;
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    if ((getX(it)<xMin)||(getX(it)>xMax)||(getY(it)<yMin)||(getY(it)>yMax)||(getZ(it)<zMin)||(getZ(it)>zMax)) {
                        types[it].setGas();
                    }
                }
            }
            cout<<"\nnoooo"<<endl;
            // SPECIAL FOX BOX: checking for boundary between gas and fluid and assigning interface properties
            double boxSpeed=(0.15/1.5)/unit.Speed;
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    for (int j=1; j<lbmDirec; ++j) {
                        int link=0;
                        link=neighbors[it].d[j];
                        if (types[link].isGas()) {
                            // BOX
                            if ((getX(it)==xMin)||(getX(it)==xMax)||(getY(it)==yMin)||(getY(it)==yMax)) {
                                types[it].setDynWall();
                                nodes[it]->u=tVect(0.0,0.0,boxSpeed);
                            }
                            // INTERFACE
                            else {
                                types[it].setInterface();
                            }
                            break;
                        }
                    }
                }
            }
            break;
        }
        case DRUM: {
            // curved shape
            cylinder dummyCylinder;
            dummyCylinder.p1=tVect(0.45,0.0,-1.0+0.0+0.04); //tVect(0.7+0.35,0.0,-1.0+0.0);
            dummyCylinder.p2=tVect(0.45,1.0,-1.0+0.0+0.04); //tVect(0.7+0.35,0.0,-1.0+0.0);
            const double fluidRadius=1.01+0.00035298*sqrt(7891604-5713*(1389.9-fluidMass));
            // const double fluidRadius=1.01+0.00035298*sqrt(7891604-5666*(1389.9-0.0045*totParticles-fluidMass));
            cout<<"Fluid radius = "<<fluidRadius<<"\n";
            dummyCylinder.R=fluidRadius;
            dummyCylinder.initAxes();
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    if (!positions[it].insideCylinder(dummyCylinder.p1/unit.Length, dummyCylinder.naxes, dummyCylinder.R/unit.Length)) {
                        types[it].setGas();
                    }
                }
            }
            //            // liquid veil
            //            cylinder dummyCylinder;
            //            dummyCylinder.p1=tVect(0.0 + 0.4,0.0,1.26);
            //            dummyCylinder.p2=tVect(0.0 + 0.4,1.0,1.26);
            //            dummyCylinder.R=1.229;
            //            dummyCylinder.initAxes();
            //            for (int it=0; it<totNodes; ++it) {
            //                if (types[it].isFluid()) {
            //                    // creating gas cells
            //                    if (getPosition(it).insideCylinder(dummyCylinder.p1/unit.Length, dummyCylinder.naxes, dummyCylinder.R/unit.Length)) {
            //                        types[it].setGas();
            //                    }
            //                }
            //            }
            break;
        }
        case DIFF: {
            // PLANE
            tVect p,n;
            p=tVect(0.0, 0.0, 40.0);
            n=tVect(0.0, 0.0, 1.0);
            for (int it=0; it<totNodes; ++it) {
                if (types[it].isFluid()) {
                    // creating gas cells
                    tVect dist=positions[it]-p;
                    if (dist.dot(n)>0.0) {
                        types[it].setGas();
                    }
                }
            }
            break;
        }
    }
    
    
    // checking for boundary between gas and fluid and assigning interface properties
    for (int it=0; it<totNodes; ++it) {
        if (types[it].isFluid()) {
            for (int j=1; j<lbmDirec; ++j) {
                unsigned int link=neighbors[it].d[j];
                if (types[link].isGas()) {
                    types[it].setInterface();
                }
            }
        }
    }
    
    // eliminating isolated interface nodes
    for (int it=0; it<totNodes; ++it) {
        if (types[it].isInterface()) {
            bool closeToFluid=false;
            for (int j=1; j<lbmDirec; ++j) {
                unsigned int link=neighbors[it].d[j];
                if (types[link].isFluid()) {
                    closeToFluid=true;
                    break;
                }
            }
            if (!closeToFluid) {
                types[it].setGas();
            }
        }
    }
}

void LB::initializeLists() {
    
    cout<<"Initializing lists"<<endl;
    
    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    fluidNodes.clear();
    particleNodes.clear();
    interfaceNodes.clear();
    
    /*    // position in space with pressure equal to import initial value
     //    tVect p0((double)lbSize[0], (double)lbSize[1], (double)lbSize[2]);
     //    #pragma omp parallel for ordered
     //    for (int z=0; z<lbSize[2]; ++z) {
     //        for (int y=0; y<lbSize[1]; ++y) {
     //            for (int x=0; x<lbSize[0]; ++x) {
     //                // vectorized position
     //                tVect p((double)(x), (double)(y), (double)(z));
     //                // distance from p0
     //                tVect deltah=p-p0;
     //                // depth increase for hydrostatic pressure computation
     //                double density=initDensity;//+3.0*dt*dt*initDensity*deltah.dot(lbF);
     //                double mass=density;
     //                double viscosity=initVisc;
     //                // velocity profile
     //                tVect velocity=initVelocity;
     //                // position variables
     //                node dummyNode;
     //
     //
     //                //          initializeNode(dummyNode, density, velocity, mass, viscosity);
     //                #pragma omp ordered
     //                nodes[index]=dummyNode;
     //                //          nodes.push_back(dummyNode);
     //            }
     //        }
     //    } */
    // creating list and initialize macroscopic variables for all nodes except walls
    for (int it=0; it<totNodes; ++it) {
        // FLUID NODES ////
        if (types[it].isFluid()) {
            // creating node
            nodes[it]=new node;
            // adding the free node indexes to the list
            fluidNodes.push_back(it);
        }
        // INTERFACE NODES ////
        else if (types[it].isInterface()) {
            // creating node
            nodes[it]=new node;
            // interface cells list update
            interfaceNodes.push_back(it);
        }
        // PARTICLE NODES ////
        if (types[it].isInsideParticle()) {
            // adding the particle node to the list
            particleNodes.push_back(it);
        }
    }
    
}

void LB::resetLists() {
    
    cout<<"Resetting lists"<<endl;
    
    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    fluidNodes.clear();
    particleNodes.clear();
    interfaceNodes.clear();
    
    
    // creating list and initialize macroscopic variables for all nodes except walls
    for (int it=0; it<totNodes; ++it) {
        // FLUID NODES ////
        if (types[it].isFluid()) {
            // adding the free node indexes to the list
            fluidNodes.push_back(it);
        }
        // INTERFACE NODES ////
        else if (types[it].isInterface()) {
            // interface cells list update
            interfaceNodes.push_back(it);
        }
        // PARTICLE NODES ////
        if (types[it].isInsideParticle()) {
            // adding the particle node to the list
            particleNodes.push_back(it);
        }
    }
    
}

void LB::initializeVariables() {
    
    cout<<"Initializing variables"<<endl;
    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    //fluidNodes.clear();
    //particleNodes.clear();
    //interfaceNodes.clear();
    // calculate maximum height of the fluid
    unsigned int maxX(0), maxY(0), maxZ(0);
    for (int it=0; it<totNodes; ++it) {
        if (types[it].isActive()) {
            maxX=std::max(maxX,getX(it));
            maxY=std::max(maxY,getY(it));
            maxZ=std::max(maxZ,getZ(it));
        }
    }
    tVect maxP(maxX,maxY,maxZ);
    // creating list and initialize macroscopic variables for all nodes except walls
    for (int it=0; it<totNodes; ++it) {
        // FLUID NODES ////
        if (types[it].isFluid()) {
            
            // setting macroscopic variables
            // density is calculated using hydrostatic profile
            tVect deltah=positions[it]-maxP;
            // density=initDensity; velocity=initVelocity, mass=initDensity; viscosity=initVisc; force=lbF
            nodes[it]->initialize(initDensity+3.0*lbmDt*lbmDt*initDensity*deltah.dot(lbF), initVelocity, initDensity, initDynVisc, lbF);
        }
        // INTERFACE NODES ////
        else if (types[it].isInterface()) {
            // setting macroscopic variables
            // density=initDensity; velocity=initVelocity, mass=0.5*initDensity; viscosity=initVisc; force=lbF
            nodes[it]->initialize(initDensity, initVelocity, 0.5*initDensity, initDynVisc, lbF);
        }
    }
}

void LB::initializeWalls(wallList& walls, cylinderList& cylinders, objectList& objects) {
    cout<<"Initializing wall nodes"<<endl;
    double zero=0.0;
    
    // initializing wall nodes
    // notice that, in the hypothesis that these walls are not evolving, only nodes at the interface need creation
    for (int it=0; it<totNodes; ++it) {
        if (!types[it].isWall()) {
            unsigned int link=0;
            for (int j=0; j<lbmDirec; ++j) {
                link=neighbors[it].d[j];
                // check if the node needs to be initialized
                if (nodes[link]==0) {
                    // initialize node
                    // STATIC WALL NODES ////
                    if ( types[link].isStatWall() || types[link].isSlipStatWall() ) {
                        // creating node
                        nodes[link]=new node;
                        // reset velocity and mass (useful for plotting)
                        // density=0.0; velocity=(0.0,0.0,0.0), mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                        nodes[link]->initialize(initDensity, Zero, zero, zero, Zero);
                    }
                    // DYNAMIC WALL NODES ////
                    else if ( types[link].isDynWall() || types[link].isSlipDynWall() || types[link].isCurvedWall() ) {
                        // creating node
                        nodes[link]=new node;
                        // need to define velocity. It could be part of a cylinder or wall, we check both
                        tVect solidVelocity;
                        tVect nodePosition=positions[link];
                        unsigned int solidIndex=types[link].getSolidIndex();
                        // wall
                        if (nodePosition.insidePlane(walls[solidIndex].p/unit.Length, walls[solidIndex].n)) {
                            solidVelocity=walls[solidIndex].getSpeed(nodePosition*unit.Length)/unit.Speed;
                        }
                        // cylinder
                        else if (!nodePosition.insideCylinder(cylinders[solidIndex].p1/unit.Length, cylinders[solidIndex].naxes, cylinders[solidIndex].R/unit.Length)) {
                            solidVelocity=cylinders[solidIndex].getSpeed(nodePosition*unit.Length)/unit.Speed;
                        }
                        // objects
                        else if (nodePosition.insideSphere(objects[solidIndex].x0/unit.Length,objects[solidIndex].r/unit.Length)) {
                            solidVelocity=objects[solidIndex].x1/unit.Speed;
                        }
                        // reset velocity and mass (useful for plotting)
                        // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
                        nodes[link]->initialize(initDensity, solidVelocity, zero, zero, Zero);
                    }
                }
            }
        }
    }
}

// integration functions

void LB::cleanLists() {
    // sorts the active-node list and removes duplicates
    
    std::sort (fluidNodes.begin(), fluidNodes.end());
    std::sort (interfaceNodes.begin(), interfaceNodes.end());
    std::sort (particleNodes.begin(), particleNodes.end());
    
    for (int ind=fluidNodes.size()-1; ind>0; --ind) {
        const unsigned int i=fluidNodes[ind];
        const unsigned int j=fluidNodes[ind-1];
        if (i==j) {
            cout<<"duplicate-fluid!"<<endl;
            fluidNodes.erase(fluidNodes.begin()+ind);
        }
    }
    for (int ind=interfaceNodes.size()-1; ind>0; --ind) {
        const unsigned int i=interfaceNodes[ind];
        const unsigned int j=interfaceNodes[ind-1];
        if (i==j) {
            cout<<"duplicate-interface!"<<endl;
            interfaceNodes.erase(interfaceNodes.begin()+ind);
        }
    }
    for (int ind=particleNodes.size()-1; ind>0; --ind) {
        const unsigned int i=particleNodes[ind];
        const unsigned int j=particleNodes[ind-1];
        if (i==j) {
            cout<<"duplicate-interface!"<<endl;
            particleNodes.erase(particleNodes.begin()+ind);
        }
    }
    
    // list with active nodes i.e. nodes where collision and streaming are solved
    // solid nodes, particle nodes and gas nodes are excluded
    activeNodes.clear();
    activeNodes.reserve(fluidNodes.size()+interfaceNodes.size());
    activeNodes.insert(activeNodes.end(), fluidNodes.begin(), fluidNodes.end());
    activeNodes.insert(activeNodes.end(), interfaceNodes.begin(), interfaceNodes.end());
    
    std::sort (activeNodes.begin(), activeNodes.end());
    
    for (int ind=activeNodes.size()-1; ind>0; --ind) {
        unsigned int i=activeNodes[ind];
        unsigned int j=activeNodes[ind-1];
        if (types[i].isGas()) {
            activeNodes.erase(activeNodes.begin()+ind);
            cout<<"wrong-gas!"<<endl;
        }
        else if (i==j) {
            activeNodes.erase(activeNodes.begin()+ind);
            cout<<"duplicate-active!"<<endl;
        }
    }
    
    for (int it=0; it<activeNodes.size(); ++it) {
        if (types[activeNodes[it]].isGas()) {
            cout<<"GAS ERROR"<<endl;
            exit(0);
        }
    }
}

void LB::reconstruction() {
    // reconstruction of macroscopic variables from microscopic distribution
    // this step is necessary to proceed to the collision step
    
#pragma omp parallel for
    for (int it=0; it<activeNodes.size(); ++it) {
        nodes[activeNodes[it]]->reconstruct();
    }
}

void LB::collision() {
    
    if (!forceField) {
        lbF.reset();
    }
    
    //    const double amplitude=lbF.norm();
    
#pragma omp parallel for
    for (int it=0; it<activeNodes.size(); ++it) {
        //        tVect forceField=lbF;
        //        switch (problemName) {
        //            case DRUM: {
        //                if (getX(it)>0.9*lbSize[0]) {
        ////                    forceField=lbF+Xm*2.0*amplitude;
        ////                    cout<<endl;
        ////                    lbF.show();
        ////                    forceField.show();
        //                }
        //                break;
        //            }
        //        }
        const unsigned int index=activeNodes[it];
        // equilibrium distributions
        double feq[lbmDirec];
        // shift velocity field to F/2
        nodes[index]->shiftVelocity(lbF);
        // compute equilibrium distributions
        nodes[index]->computeEquilibrium(feq);
        //        if (types[index].isInsideParticle()) {
        //            nodes[index]->resetViscosity(initDynVisc);
        //        }
        if (nonNewtonian || turbulenceOn) {
            // compute shear rate tensor, find invariant calculate viscosity (Bingham)
            nodes[index]->computeShearRate(feq, turbConst, plasticVisc, yieldStress,nonNewtonian,turbulenceOn);
        }
        // compute new distributions
        nodes[index]->solveCollision(feq);
        // add force term to new distributions
        nodes[index]->addForce(lbF);
    }
    
    
    //    if (!forceField && !nonNewtonian) {
    //        #pragma omp parallel for
    //        for (int it=0; it<activeNodes.size(); ++it) {
    //            nodes[activeNodes[it]]->collide(initDynVisc);
    //        }
    //    }
    //
    //    if (forceField && !nonNewtonian) {
    ////        #pragma omp parallel for
    //        for (int it=0; it<activeNodes.size(); ++it) {
    //            nodes[activeNodes[it]]->collideF(lbF, initDynVisc);
    //        }
    //    }
    //
    //    if (!forceField && nonNewtonian) {
    //        #pragma omp parallel for
    //        for (int it=0; it<activeNodes.size(); ++it) {
    //            nodes[activeNodes[it]]->collideNN(turbConst, plasticVisc, yieldStress);
    //        }
    //    }
    //
    //    if (forceField && nonNewtonian) {
    //        #pragma omp parallel for
    //        for (int it=0; it<activeNodes.size(); ++it) {
    //            nodes[activeNodes[it]]->collideFNN(turbConst, lbF, plasticVisc, yieldStress);
    //        }
    //    }
}

void LB::streaming(elmtList& elmts, particleList& particles, wallList& walls) {
    // STREAMING STEP
    // coefficients for equilibrium
    static double C1=3.0*lbmDt*lbmDt;
    static double C2=4.5*lbmDt*lbmDt*lbmDt*lbmDt;
    static double C3=1.5*lbmDt*lbmDt;
    // coefficient for free-surface
    static const double C2x2=9.0*lbmDt*lbmDt*lbmDt*lbmDt;
    static const double C3x2=3.0*lbmDt*lbmDt;
    // coefficient for slip conditions
    static const double S1=slipCoefficient;
    static const double S2=(1.0-slipCoefficient);
    // creating list for collision function
    double staticPres[lbmDirec];
    for (int j=0; j<lbmDirec; j++) {
        staticPres[j]=initDensity*coeff[j];
    }
    // coefficient for bounce-back
    static const double BBCoeff=2.0*3.0*lbmDt*lbmDt;
    // extra mass due to bounce-back and moving walls
    double extraMass=0.0;
    
    //    // initializing the elements forces (lattice units)
    //    #pragma omp parallel for
    //    for (int n=0; n<elmts.size(); ++n) {
    //        // storing for stabilization
    //        if (elmts[n].stabNumber>=2) {
    //        for (int s=elmts[n].stabNumber-1; s>=1; --s) {
    //
    //                elmts[n].FHydroOld[s]=elmts[n].FHydroOld[s-1];
    //                elmts[n].MHydroOld[s]=elmts[n].MHydroOld[s-1];
    //            }
    //        }
    //        if (elmts[n].stabNumber>=1) {
    //            elmts[n].FHydroOld[0]=elmts[n].FHydro;
    //            elmts[n].MHydroOld[0]=elmts[n].MHydro;
    //        }
    //        //initializing this time step hydrodynamic force
    //        elmts[n].FHydro.reset();
    //        elmts[n].MHydro.reset();
    //    }
    
    // initializing wall forces
    for (int iw=0; iw<walls.size(); ++iw) {
        walls[iw].FHydro.reset();
    }
    
    //  Saving in support variables f->fs
#pragma omp parallel for
    for (int it=0; it<activeNodes.size(); ++it) {
        nodes[activeNodes[it]]->store();
    }
    
    
    
    //    for (int it=0; it<activeNodes.size(); ++it) { // this part is to be DELETED!!!
    //        nodes[activeNodes[it]].flag=false;
    //        nodes[activeNodes[it]].force.reset();
    //        for (int j=1; j<direc; ++j) {
    //            if (nodes[nodes[activeNodes[it]].d[j]].isWall()) {
    //                nodes[activeNodes[it]].flag=true;
    //            }
    //        }
    //    }
    
    //    for (int n=0; n<nodes.size(); ++n) {
    //        nodes[n].u.reset();//+dt*lbF/2.0/n;
    //        nodes[n].setEquilibrium(vZero);
    //    }
    
    
    //  Streaming
#pragma omp parallel for ordered reduction(+:extraMass)
    // cycling through active nodes
    for (int in=0; in<activeNodes.size(); ++in) {
        
        // index of active node
        const unsigned int it=activeNodes[in];
        // saving in a support node variable
        //        node dummyNode=nodes[activeNodes[it]];
        //cycling through neighbors
        for (int j=1; j<lbmDirec; ++j) {
            // getting neighbor index
            // variable for linked node index
            const unsigned int link=neighbors[it].d[j];
            
            // if neighbor is normal fluid cell what follows is true
            if (types[link].isActive()) {
                // streaming is solved normally
                //                dummyNode.stream(opp[j], nodes[link], opp[j], one, zero);
                nodes[it]->f[opp[j]]=nodes[link]->fs[opp[j]];
            }
            
            else if (types[link].isGas()) {
                // additional variables for equilibrium f computation
                const double usq=nodes[it]->u.norm2();
                const double vuj=nodes[it]->u.dot(v[j]);
                //streaming with constant pressure interface
                // comes from dummyNode.f[opp[j]]=coeff[j]*minPres*(1.0+3.0*vuj+4.5*vuj*vuj-1.5*usq)+coeff[opp[j]]*minPres*(1.0+3.0*vujOpp+4.5*vujOpp*vujOpp-1.5*usq)-dummyNode.fs[j];
                //                dummyNode.stream(opp[j], dummyNode, j, -1.0, coeff[j]*minPres*(2.0+3.0*dt*dt*(vuj+vujOpp)+4.5*dt*dt*dt*dt*(vuj*vuj+vujOpp*vujOpp)-3.0*dt*dt*usq));
                //                dummyNode.f[opp[j]]=-dummyNode.fs[j]+coeff[j]*initDensity*(2.0+3.0*dt*dt*(vuj+vujOpp)+4.5*dt*dt*dt*dt*(vuj*vuj+vujOpp*vujOpp)-3.0*dt*dt*usq);
                nodes[it]->f[opp[j]]=-nodes[it]->fs[j]+coeff[j]*initDensity*(2.0+C2x2*(vuj*vuj)-C3x2*usq);
            }
            /*
             //            else if (types[link].isInsideParticle()) {
             //                // getting the index of the particle to compute force in the right object
             //                unsigned int particleIndex=types[link].getSolidIndex();
             //                unsigned int clusterIndex=particles[particleIndex].clusterIndex;
             //                // calculating velocity of the solid boundary at the node (due to rotation of particles)
             //                // vectorized radius (lattice units)
             //                tVect radius=getPosition(link)-particles[particleIndex].x0/unit.Length+particles[particleIndex].radiusVec/unit.Length;
             //                // update velocity of the particle node (u=v_center+omega x radius) (lattice units)
             //                tVect localVel=elmts[clusterIndex].x1/unit.Speed+(elmts[clusterIndex].w.cross(radius))/unit.AngVel;
             //                // variation in Bounce-Back due to moving object (lattice units)
             //                double BBi=BBCoeff*nodes[it]->n*coeff[j]*localVel.dot(v[j]); //CHECK THIS LINE!!!! (1.0-initDensity/dummyNode.n)? *dummyNode.mass
             //                // force due to momentum exchange. (lattice units)
             //                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
             //                tVect BBforce=nodes[it]->bounceBackForce(j, staticPres, BBi);
             //                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
             ////                if (!dummyNode.flag) {
             //                #pragma omp critical
             //                {
             //                    elmts[clusterIndex].FHydro+=BBforce;
             //                    elmts[clusterIndex].MHydro+=radius.cross(BBforce);
             ////                    dummyNode.force+=BBforce;
             //                }
             ////                }
             //                // adding the extra mass to the surplus
             //                extraMass+=nodes[it]->mass*BBi;
             //                // streaming in case of Bounce-Back
             ////                dummyNode.stream(opp[j], dummyNode, j, 1.0, -BBi);
             //                nodes[it]->f[opp[j]]=nodes[it]->fs[j]-BBi;
             //            }    */
            // for curved walls there is the rule of Mei-Luo-Shyy
            else if (types[link].isCurvedWall()) {
                //                    // getting the index of the wall to compute force in the right object
                //                    const int solidIndex=types[link].getSolidIndex();
                // velocity of the wall
                const tVect vel=nodes[link]->u;
                // variation in Bounce-Back due to moving object
                const double BBi=BBCoeff*nodes[it]->n*coeff[j]*vel.dot(v[j]);// mass!!!!!
                
                //                    // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                //                    const tVect BBforce=nodes[it]->bounceBackForce(j, staticPres, BBi);
                //                    // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
                //                    #pragma omp ordered
                //                    {
                //                        walls[solidIndex].FHydro+=BBforce;
                //                    }
                // get value for chi coefficient
                const double chi=curves[link]->getChi(opp[j], 0.5+3.0*nodes[it]->visc);
                // get fluid reference velocity
                tVect ubf(0.0,0.0,0.0);
                if (curves[link]->delta[opp[j]]>=0.5) {
                    ubf=curves[link]->m1[opp[j]]*nodes[it]->u+curves[link]->m2[opp[j]]*vel;
                }
                else {
                    const unsigned int linkBack=neighbors[it].d[opp[j]];
                    if (types[linkBack].isFluid() && nodes[linkBack]!=0) {
                        ubf=nodes[linkBack]->u;
                    }
                    else {
                        ubf=vel;
                    }
                }
                // calculate missing distribution
                const double usq=nodes[it]->u.norm2();
                const double vu=nodes[it]->u.dot(v[j]);
                const double fStar=nodes[it]->n*coeff[j]*(1+C1*v[j].dot(ubf)+C2*vu*vu+C3*usq);
                
                // stream
                //                    dummyNode.stream(opp[j], dummyNode, j, 1.0, -BBi);
                nodes[it]->f[opp[j]]=(1.0-chi)*nodes[it]->fs[j]+chi*fStar-BBi;
                // adding the extra mass to the surplus
                extraMass+=nodes[it]->mass*(chi*nodes[it]->fs[j]-chi*fStar+BBi);
            }
            // for moving walls there is simple bounce-back with velocity correction
            else if (types[link].isDynWall()) {
                //                    // getting the index of the wall to compute force in the right object
                const int solidIndex=types[link].getSolidIndex();
                // velocity of the wall
                const tVect vel=nodes[link]->u;
                // variation in Bounce-Back due to moving object
                const double BBi=BBCoeff*nodes[it]->n*coeff[j]*vel.dot(v[j]);// mass!!!!!
                
                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce=nodes[it]->bounceBackForce(j, staticPres, BBi);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                {
                    walls[solidIndex].FHydro+=BBforce;
                }
                // stream
                //                    dummyNode.stream(opp[j], dummyNode, j, 1.0, -BBi);
                nodes[it]->f[opp[j]]=nodes[it]->fs[j]-BBi;
                // adding the extra mass to the surplus
                extraMass+=BBi*nodes[it]->mass;
            }
            
            // for walls there is simple bounce-back
            else if (types[link].isStatWall()) {
                //                // getting the index of the wall to compute force in the right object
                //                int solidIndex=nodes[link].solidIndex;
                //                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                //                tVect BBforce=dummyNode.bounceBackForce(j, staticPres, zero);
                //                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
                //                #pragma omp ordered
                //                {
                //                    walls[solidIndex].force+=BBforce;
                //                }
                //                dummyNode.stream(opp[j], dummyNode, j, 1.0, 0.0);
                nodes[it]->f[opp[j]]=nodes[it]->fs[j];
            }
            
            else if (types[link].isSlipStatWall()) {
                //                // getting the index of the wall to compute force in the right object
                //                int solidIndex=nodes[link].solidIndex;
                //                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                //                tVect BBforce=dummyNode.bounceBackForce(j, staticPres, zero);
                //                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
                //                #pragma omp ordered
                //                {
                //                    walls[solidIndex].force+=BBforce;
                //                }
                //                dummyNode.stream(opp[j], dummyNode, j, 1.0, 0.0);
                if (j>6) {
                    bool active1=false;
                    bool active2=false;
                    const unsigned int nodeCheck1=neighbors[it].d[slip1Check[j]];
                    const unsigned int nodeCheck2=neighbors[it].d[slip2Check[j]];
                    // check for the environment
                    if (types[nodeCheck1].isActive()) {
                        active1=true;
                    }
                    if (types[nodeCheck2].isActive()) {
                        active2=true;
                    }
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes[it]->f[opp[j]]=S1*nodes[nodeCheck1]->fs[slip1[j]]+S2*nodes[it]->fs[j];
                        
                    }
                    else if (!active1 && active2) {
                        // second
                        nodes[it]->f[opp[j]]=S1*nodes[nodeCheck2]->fs[slip2[j]]+S2*nodes[it]->fs[j];
                    }
                    else {
                        // standard BB
                        nodes[it]->f[opp[j]]=nodes[it]->fs[j];
                    }
                }
                else {
                    // standard BB
                    nodes[it]->f[opp[j]]=nodes[it]->fs[j];
                }
            }
            
            else if (types[link].isSlipDynWall()) {
                //                // getting the index of the wall to compute force in the right object
                //                int solidIndex=nodes[link].solidIndex;
                //                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                //                tVect BBforce=dummyNode.bounceBackForce(j, staticPres, zero);
                //                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
                //                #pragma omp ordered
                //                {
                //                    walls[solidIndex].force+=BBforce;
                //                }
                //                dummyNode.stream(opp[j], dummyNode, j, 1.0, 0.0);
                // velocity of the wall
                const tVect vel=nodes[link]->u;
                // variation in Bounce-Back due to moving object
                const double BBi=BBCoeff*nodes[it]->n*coeff[j]*vel.dot(v[j]);
                
                if (j>6) {
                    bool active1=false;
                    bool active2=false;
                    const unsigned int nodeCheck1=neighbors[it].d[slip1Check[j]];
                    const unsigned int nodeCheck2=neighbors[it].d[slip2Check[j]];
                    // check for the environment
                    if (types[nodeCheck1].isActive()) {
                        active1=true;
                    }
                    if (types[nodeCheck2].isActive()) {
                        active2=true;
                    }
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes[it]->f[opp[j]]=S1*nodes[nodeCheck1]->fs[slip1[j]]+S2*(nodes[it]->fs[j]-BBi);
                        // adding the extra mass to the surplus
                        extraMass+=S2*nodes[it]->mass*BBi;
                    }
                    else if (!active1 && active2) {
                        // second
                        nodes[it]->f[opp[j]]=S1*nodes[nodeCheck2]->fs[slip2[j]]+S2*(nodes[it]->fs[j]-BBi);
                        // adding the extra mass to the surplus
                        extraMass+=S2*nodes[it]->mass*BBi;
                    }
                    else {
                        // standard BB
                        nodes[it]->f[opp[j]]=nodes[it]->fs[j]-BBi;
                        // adding the extra mass to the surplus
                        extraMass+=nodes[it]->mass*BBi;
                    }
                }
                else {
                    // standard BB
                    nodes[it]->f[opp[j]]=nodes[it]->fs[j]-BBi;
                    // adding the extra mass to the surplus
                    extraMass+=nodes[it]->mass*BBi;
                }
            }
            
            else {
                cout<<link<<" "<<getX(link)<<" "<<getY(link)<<" "<<getZ(link)<<" "<<types[link].getType()<<" TYPE ERROR"<<endl;
                exit(0);
            }
        }
    }
    
    //    unsigned int targetNode=0;
    //    double maxForce=0.0;
    //    for (int it=0; it<activeNodes.size(); ++it) { // this part is to be DELETED!!!
    //        if (nodes[activeNodes[it]].force.norm()>maxForce) {
    //            maxForce=nodes[activeNodes[it]].force.norm();
    //            targetNode=activeNodes[it];
    //        }
    //    }
    //    cout<<"Target="<<targetNode<<"\t";
    
    
    // redistributing extra mass due to bounce back to interface cells
    redistributeMass(extraMass);
    
    //    // shifting elements forces and torques to physical units
    //    for (int n=0; n<elmts.size(); ++n) {
    //        elmts[n].FHydro*=unit.Force;
    //        elmts[n].MHydro*=unit.Torque;
    //    }
    
    for (int n=0; n<walls.size(); ++n) {
        walls[n].FHydro*=unit.Force;
    }
}

// free surface functions

void LB::updateMass() {
    // refer to the article of Miller or of Svec et al. for this
    
#pragma omp parallel for
    for (int it=0; it<interfaceNodes.size(); ++it) {
        nodes[interfaceNodes[it]]->newMass=nodes[interfaceNodes[it]]->mass;
    }
    
    // mass for interface nodes is regulated by the evolution equation
#pragma omp parallel for
    for (int it=0; it<interfaceNodes.size(); ++it) {
        const int index=interfaceNodes[it];
        // additional mass streaming to/from interface
        double deltaMass=0.0;
        //cycling through neighbors
        for (int j=1; j<lbmDirec; ++j) {
            // getting neighbor index
            const unsigned int link=neighbors[index].d[j];
            // average liquid fraction
            double averageMass=0.0;
            if (types[link].isInterface()) {
                // average liquid fraction
                averageMass=0.5*(nodes[link]->mass+nodes[index]->mass);
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            }
            else if (types[link].isFluid()) {
                averageMass=1.0;
                //                dummyNode.f[opp[j]]=nodes[link].fs[opp[j]];
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            }
            //            else if (types[link].isParticle()) {
            //                averageMass=1.0*nodes[index]->mass;
            ////                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            //            }
            else if (types[link].isGas()) {
                averageMass=0.0;
                //                averageMass=1.0;//*nodes[index].mass;
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            }
            else if (types[link].isDynWall()) {
                averageMass=1.0*nodes[index]->mass;//1.0*nodes[index].mass;//0.5*nodes[index].u.norm()/(nodes[link].u.norm()+nodes[index].u.norm());
                //              double BBi=2.0*3.0*dummyNode.n*coeff[j]*vel.dot(v[j]);
                //              dummyNode.f[opp[j]]=dummyNode.fs[j]-BBi;
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
                //                cout<<"delta mass ="<<averageMass<<","<<deltaMass<<"\n";
            }
            else if (types[link].isCurvedWall()) {
                averageMass=1.0*nodes[index]->mass;//1.0*nodes[index].mass;//0.5*nodes[index].u.norm()/(nodes[link].u.norm()+nodes[index].u.norm());
                //              double BBi=2.0*3.0*dummyNode.n*coeff[j]*vel.dot(v[j]);
                //              dummyNode.f[opp[j]]=dummyNode.fs[j]-BBi;
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
                //                cout<<"delta mass ="<<averageMass<<","<<deltaMass<<"\n";
            }
            else if (types[link].isSlipDynWall()) {
                if (j>6) {
                    bool active1=false;
                    bool active2=false;
                    const unsigned int nodeCheck1=neighbors[index].d[slip1Check[j]];
                    const unsigned int nodeCheck2=neighbors[index].d[slip2Check[j]];
                    // check for the environment
                    if (types[nodeCheck1].isActive()) {
                        active1=true;
                    }
                    if (types[nodeCheck2].isActive()) {
                        active2=true;
                    }
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // adding the extra mass to the surplus
                        averageMass+=1.0*(1.0-slipCoefficient)*nodes[index]->mass;
                    }
                    else if (!active1 && active2) {
                        // adding the extra mass to the surplus
                        averageMass+=1.0*(1.0-slipCoefficient)*nodes[index]->mass;
                    }
                    else {
                        // adding the extra mass to the surplus
                        averageMass+=1.0*nodes[index]->mass;
                    }
                }
                else {
                    // adding the extra mass to the surplus
                    averageMass+=1.0*nodes[index]->mass;
                }
            }
            deltaMass+=lbmDt*averageMass*nodes[index]->massStream(j);
        }
        nodes[index]->newMass+=deltaMass;
    }
    // mass for fluid nodes is equal to density
#pragma omp parallel for
    for (int it=0; it<fluidNodes.size(); ++it) {
        nodes[fluidNodes[it]]->mass=nodes[fluidNodes[it]]->n;
    }
#pragma omp parallel for
    for (int it=0; it<interfaceNodes.size(); ++it) {
        nodes[interfaceNodes[it]]->mass=nodes[interfaceNodes[it]]->newMass;
    }
}

void LB::updateInterface() {
    // updates the location of the interface, creating and deleting interface nodes
    
    // variable for storage of mass surplus
    double massSurplus=0.0;
    
    // lists for "mutant" nodes
    unsIntList filledNodes, emptiedNodes;    
    // interface nodes created due to interface evolution
    unsIntList newInterfaceNodes;
    
    // filling lists of mutant nodes and changing their type
    findInterfaceMutants(filledNodes, emptiedNodes);
    
    // fixing the interface (always one interface between fluid and gas)
    smoothenInterface(filledNodes, emptiedNodes, newInterfaceNodes, massSurplus);
    
    // updating characteristics of mutant nodes
    updateMutants(filledNodes, emptiedNodes, newInterfaceNodes, massSurplus);
    
    // remove isolated interface cells (both surrounded by gas and by fluid)
    removeIsolated(massSurplus);
    
    // distributing surplus to interface cells
    redistributeMass(massSurplus);
    
}

void LB::findInterfaceMutants(unsIntList& filledNodes, unsIntList& emptiedNodes) {
    // checks for nodes going through the following transitions
    // interface -> fluid (filled node)
    // interface -> gas (emptied node)
    filledNodes.clear();
    emptiedNodes.clear();
    
    // CHECKING FOR MUTANT NODES (old backward style for cycle...)
    for (int ind=interfaceNodes.size()-1; ind>=0; --ind) {
        const unsigned int i=interfaceNodes[ind];
        // CHECKING FOR NEW FLUID NODES from filling
        if (nodes[i]->mass>nodes[i]->n) {
            // updating mutant list
            filledNodes.push_back(i);
            // updating type
            types[i].setFluid();
            // updating lists (interface -> fluid)
            fluidNodes.push_back(i);
            interfaceNodes.erase(interfaceNodes.begin()+ind);
        }
        // CHECKING FOR NEW GAS NODES from emptying
        else if (nodes[i]->mass<0.0) {
            // updating mutant list
            emptiedNodes.push_back(i);
            // updating type
            types[i].setGas();
            // updating lists (interface ->noList)
            interfaceNodes.erase(interfaceNodes.begin()+ind);
        }
    }
}

void LB::smoothenInterface(unsIntList& filledNodes, unsIntList& emptiedNodes, unsIntList& newInterfaceNodes, double& massSurplus) {
    
    newInterfaceNodes.clear();
    
    // CHECKING FOR NEW INTERFACE NODES from neighboring a new fluid node
    for (int it=0; it<filledNodes.size(); ++it) {
        const unsigned int index=filledNodes[it];
        // cycling through neighbors
        for (int j=1; j<lbmDirec; ++j) {
            // index of neighbor node
            const unsigned int link=neighbors[index].d[j];
            // checking if node is gas (so to be tranformed into interface)
            if (types[link].isGas()) {
                types[link].setInterface();
                newInterfaceNodes.push_back(link);
                // allocating memory
                nodes[link]=new node;
                // node is becoming active and needs to be initialized
                // same density and velocity; 1% of the mass
                nodes[link]->initialize(initDensity, nodes[index]->u, 0.01*initDensity, nodes[index]->visc, nodes[index]->hydroForce+lbF);
                // the 1% of the mass is taken form the surplus
                massSurplus-=0.01*initDensity;
                //                activeNodes.push_back(nodes[link].index);
            }
        }
    }
    
    // CHECKING FOR NEW INTERFACE NODES from neighboring a new gas node
    for (int it=0; it<emptiedNodes.size(); ++it) {
        const unsigned int index=emptiedNodes[it];
        for (int j=0; j<lbmDirec; ++j) {
            const unsigned int link=neighbors[index].d[j];
            if (types[link].isFluid()) {
                types[link].setInterface();
                newInterfaceNodes.push_back(link);
                // characteristics are inherited by previous fluid cell. Only mass must be updated to 99% of initial mass
                nodes[link]->mass=0.99*nodes[link]->n;
                // the remaining 1% of the mass is added to the surplus
                massSurplus+=0.01*nodes[link]->n;
                // removing from fluid nodes
                for (int ind=fluidNodes.size()-1; ind>=0; --ind) {
                    const unsigned int i=fluidNodes[ind];
                    if (i==link) {
                        fluidNodes.erase(fluidNodes.begin()+ind);
                    }
                }
            }
        }
    }
}

void LB::updateMutants(unsIntList& filledNodes, unsIntList& emptiedNodes, unsIntList& newInterfaceNodes, double& massSurplus) {
    
    // resetting new gas macroscopic quantities
    for (int it=0; it<emptiedNodes.size(); ++it) {
        bool isNotInterface=true;
        for (int it1=0; it1<newInterfaceNodes.size(); ++it1) {
            if (emptiedNodes[it]==newInterfaceNodes[it1]) {
                isNotInterface=false;
            }
        }
        if (isNotInterface) {
            // updating mass surplus
            massSurplus+=nodes[emptiedNodes[it]]->mass; // maybe problems
            // deleting node
            delete nodes[emptiedNodes[it]];
            nodes[emptiedNodes[it]]=0;
        }
    }
    
    // resetting new fluid macroscopic quantities
    for (int it=0; it<filledNodes.size(); ++it) {
        bool isNotInterface=true;
        for (int it1=0; it1<newInterfaceNodes.size(); ++it1) {
            if (filledNodes[it]==newInterfaceNodes[it1]) {
                isNotInterface=false;
            }
        }
        if (isNotInterface) {
            // updating mass surplus
            massSurplus+=nodes[filledNodes[it]]->mass-nodes[filledNodes[it]]->n;
            // setting liquid fraction for new fluid cell (other macroscopic characteristics stay the same)
            nodes[filledNodes[it]]->mass=nodes[filledNodes[it]]->n;
        }
    }
    // resetting neighbors for newly created interface cells
    for (int it=0; it<newInterfaceNodes.size(); ++it) {
        interfaceNodes.push_back(newInterfaceNodes[it]);
    }
    
}

void LB::removeIsolated(double& massSurplus) {
    // remove isolated interface cells (surrounded either by only fluid or only solid cells)
    
    // checking if it is surrounded by fluid (in that case is converted to fluid). Solid is an exception
    // reverse cycle is needed because of deletion function
    for (int ind=interfaceNodes.size()-1; ind>=0; --ind) {
        const unsigned int i=interfaceNodes[ind];
        bool surroundedFluid=true;
        for (int j=1; j<lbmDirec; ++j) {
            const unsigned int link=neighbors[i].d[j];
            if  (types[link].isGas()) {
                surroundedFluid=false;
                break;
            }
        }
        if (surroundedFluid) {
            //            cout<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].z<<" NEW FLUID NODE from surrounding\n";
            // update mass storage for balance
            massSurplus+=nodes[i]->mass-nodes[i]->n;
            // update characteristics (inherited from the gas node)
            nodes[i]->mass=nodes[i]->n;
            types[i].setFluid();
            interfaceNodes.erase(interfaceNodes.begin()+ind);
            fluidNodes.push_back(i);
        }
    }
    
    // checking if it is surrounded by gas (in that case is converted to gas)
    // or, better, if it is not connected to fluid (could be connected to walls or particles)
    for (int ind=interfaceNodes.size()-1; ind>=0; --ind) {
        const unsigned int i=interfaceNodes[ind];
        bool surroundedGas=true;
        for (int j=1; j<lbmDirec; ++j) {
            const unsigned int link=neighbors[i].d[j];
            if  (types[link].isFluid()) {
                surroundedGas=false;
                break;
            }
        }
        // updating mass surplus
        if (surroundedGas) {
            //            cout<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].z<<" NEW GAS NODE from surrounding\n";
            // update mass
            massSurplus+=nodes[i]->mass;
            interfaceNodes.erase(interfaceNodes.begin()+ind);
            types[i].setGas();
            delete nodes[i];
            nodes[i]=0;
        }
    }
}

void LB::redistributeMass(const double& massSurplus) {
    // redistribute the mass surplus among interface cells
    const double addMass=massSurplus/interfaceNodes.size();
    
#pragma omp parallel for
    for (int it=0; it<interfaceNodes.size(); ++it) {
        nodes[interfaceNodes[it]]->mass+=addMass;
    }
}

void LB::enforceMassConservation(){
    // calculate total mass
    double thisMass=0.0;
    for (int it=0; it<nodes.size(); ++it){
        if (types[it].isActive() && !types[it].isInsideParticle()) {
            thisMass+=nodes[it]->mass;
        }
    }
    
    // mass deficit
    const double massDeficit=(thisMass-totalMass);
    // cout<<endl<<"This="<<thisMass<<" tot="<<totalMass;
    
    // fix it
    redistributeMass(-0.01*massDeficit);
    
}

// particle coupling functions

void LB::updateIndices(elmtList& elmts, particleList& particles) {
    // checks if solid indices need to be changed due to periodic shifting
    // this is a compulsory step when working with ghost particles
    
    
    initializeParticleBoundaries(particles);
    resetLists();
    
    
    ////    cout<<"Updating indices\n";
    ////    #pragma omp parallel for
    //    for (int ip=0; ip<particleNodes.size(); ++ip) {
    //        const unsigned int index=particleNodes[ip];
    //        const tVect nodePosition=positions[index];
    //        for (int p=0; p<particles.size();++p) {
    //            if (nodePosition.insideSphere(particles[p].x0/unit.Length,particles[p].r/unit.Length)) { //-0.5?
    //                // if the node is inside the element, then set the indices again and leave
    //                types[index].setSolidIndex(particles[p].particleIndex);
    //                break;
    //            }
    //        }
    //    }
    
}

void LB::computeHydroForces(elmtList& elmts, particleList& particles) {
    
    // initializing the elements forces (lattice units)
#pragma omp parallel for
    for (int n=0; n<elmts.size(); ++n) {
        //initializing this time step hydrodynamic force
        elmts[n].FHydro=tVect(0.0,0.0,0.0);
        elmts[n].MHydro=tVect(0.0,0.0,0.0);
        // initializing the fluid mass for buoyancy
        elmts[n].fluidVolume=0.0;
    }
    
#pragma omp parallel for
    // cycling through active nodes
    for (int ip=0; ip<activeNodes.size(); ++ip) {
        unsigned int index=activeNodes[ip];
        // resetting hydrodynamic forces on nodes
        nodes[index]->hydroForce.reset();
        if (types[index].isInsideParticle()) { // && types[index].isFluid()
            // getting the index of the particle to compute force in the right object
            const unsigned int particleIndex=types[index].getSolidIndex();
            const unsigned int clusterIndex=particles[particleIndex].clusterIndex;
            // calculating velocity of the solid boundary at the node (due to rotation of particles)
            // vectorized radius (real units)
            const tVect radius=positions[index]-particles[particleIndex].x0/unit.Length+particles[particleIndex].radiusVec/unit.Length;
            // update velocity of the particle node (u=v_center+omega x radius) (real units)
            const tVect localVel=elmts[clusterIndex].x1/unit.Speed+(elmts[clusterIndex].wGlobal.cross(radius))/unit.AngVel;
            
            // calculate differential velocity
            const tVect diffVel=nodes[index]->liquidFraction()*(nodes[index]->u-localVel);
            
            // force on fluid
            nodes[index]->hydroForce+=-1.0*diffVel; // -1.0*nodes[index]->liquidFraction()*diffVel;
            
            //             // force on fluid
            //            # pragma omp critical
            //            {
            //                for (int j=0; j<lbmDirec; ++j) {
            //                    const unsigned int link=neighbors[index].d[j];
            //                    if (types[link].isActive()) {
            //                        nodes[link]->hydroForce+=-1.0/nodes[index]->liquidFraction()*coeff[j]*diffVel;
            //                    }
            //                }
            //            }
            
            // force on particle
#pragma omp critical
            {
                elmts[clusterIndex].fluidVolume+=nodes[index]->mass;
                elmts[clusterIndex].FHydro+=1.0*diffVel;
                elmts[clusterIndex].MHydro+=1.0*radius.cross(diffVel);
            }
        }
    }
    
    
    // shifting elements forces and torques to physical units
#pragma omp parallel for
    for (int n=0; n<elmts.size(); ++n) {
        
        // adding buoyancy
        const tVect buoyancy=(-1.0)*elmts[n].fluidVolume*lbF;
        
        //        elmts[n].FHydro+=buoyancy;
        elmts[n].FHydro*=unit.Force;
        elmts[n].MHydro*=unit.Torque;
        elmts[n].fluidVolume*=unit.Volume;
    }
}

void LB::findNewActive(unsIntList& newPopUpNodes, elmtList& elmts, particleList& particles) {
    
    // SOLID TO ACTIVE CHECK
    // cycling through particle nodes
    //    #pragma omp parallel for ordered
    for (int ip=0; ip<particleNodes.size(); ++ip) {
        const unsigned int index=particleNodes[ip];
        const tVect nodePosition=positions[index];
        // solid index to identify cluster
        const unsigned int particleIndex=types[index].getSolidIndex();
        const unsigned int clusterIndex=particles[particleIndex].clusterIndex;
        // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
        bool newActive=true;
        // cycling through component particles
        unsigned int componentIndex=0;
        for (int j=0; j<elmts[clusterIndex].components.size(); ++j) {
            // getting indexes from particle composing the cluster
            componentIndex=elmts[clusterIndex].components[j];
            // checking if it has been uncovered in component j of the cluster
            // radius need to be increased by half a lattice unit
            // this is because solid boundaries are located halfway between soli and fluid nodes
            if (nodePosition.insideSphere(particles[componentIndex].x0/unit.Length,particles[componentIndex].r/unit.Length)) { //-0.5?
                // if the node is still inside the element, the hypothesis of new active is not true anymore
                newActive=false;
                // and we can get out of the cycle
                break;
            }
        }
        if (newActive) {
            // turning up the cell
            //            #pragma omp ordered
            newPopUpNodes.push_back(index);
            types[index].setOutsideParticle();
            //            cout<<"new active\n";
        }
    }
    
    for (int it=0; it<newPopUpNodes.size(); ++it) {
        //        cout<<"New active\n";
        for (int ind=particleNodes.size()-1; ind>=0; --ind) {
            if (newPopUpNodes[it]==particleNodes[ind]) {
                // deleting the cell from particle list
                particleNodes.erase(particleNodes.begin()+ind);
            }
        }
    }
}

void LB::findNewSolid(unsIntList& newSolidNodes, elmtList& elmts, particleList& particles) {
    // ACTIVE TO SOLID CHECK
    // we check only first order neighbors of particle nodes. This is done in order to avoid cycling through all active cells
    //    #pragma omp parallel for ordered
    for (int it=0; it<particleNodes.size(); ++it) {
        const unsigned int index=particleNodes[it];
        //    for (intlist::iterator it=particleNodes.begin(); it<particleNodes.end(); ++it) {
        // solid index to identify cluster
        const unsigned int particleIndex=types[index].getSolidIndex();
        const unsigned int clusterIndex=particles[particleIndex].clusterIndex;
        // cycle through first neighbors
        for (int k=1; k<lbmMainDirec; ++k) {
            const unsigned int link=neighbors[index].d[k];
            // checking if solid particle is close to an active one -> we have an active node to check
            if (!types[link].isInsideParticle()) {
                const tVect linkPosition=positions[link];
                // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                bool newSolid=false;
                unsigned int componentIndex=0;
                // cycling through all components of the cluster
                for (int j=0; j<elmts[clusterIndex].components.size(); ++j) {
                    //                    cout<<j<<"\n";
                    // getting component particle index
                    componentIndex=elmts[clusterIndex].components[j];
                    // check if it getting inside
                    // radius need to be increased by half a lattice unit
                    // this is because solid boundaries are located halfway between soli and fluid nodes
                    if  (linkPosition.insideSphere(particles[componentIndex].x0/unit.Length,particles[componentIndex].r/unit.Length)) { //-0.5?
                        // if so, then the false hypothesis does not hold true anymore
                        newSolid=true;
                        // and we exit the cycle
                        break;
                    }
                }
                // an active cell can become part of multiple particle at the same time
                // we need to check if it is already on the list
                if (newSolid) {
                    // check if we are creating a duplicate - we start with a false hypothesis
                    bool alreadyInside=false;
                    for (int it1=0; it1<newSolidNodes.size(); ++it1) {
                        //                    for (intlist::iterator it1=newSolidNodes.begin(); it1!=newSolidNodes.end(); ++it1) {
                        if (link==newSolidNodes[it1]) {
                            // if it is already in the list than the true hypothesis does not hold true anymore
                            alreadyInside=true;
                            // and we exit the cycle
                            break;
                        }
                    }
                    // at this point, if it is both newSolid and not alreadyInside we can add it to particle nodes
                    if (!alreadyInside) {
                        // solid index is assigned here to avoid sending this information to the switching functions
                        //                        # pragma omp ordered
                        {
                            types[link].setSolidIndex(componentIndex);
                            types[link].setInsideParticle();
                            // sending node index for updating
                            newSolidNodes.push_back(link);
                            particleNodes.push_back(link);
                        }
                    }
                }
            }
        }
    }
}

void LB::solidToActive(unsIntList& newPopUpNodes, elmtList& elmts, double& massSurplus) {
    
    //    for (int it=0; it<newPopUpNodes.size(); ++it) {
    ////        cout<<"New active\n";
    //        for (int ind=particleNodes.size()-1; ind>=0; --ind) {
    //            if (newPopUpNodes[it]==particleNodes[ind]) {
    //                // deleting the cell from particle list
    //                particleNodes.erase(particleNodes.begin()+ind);
    //            }
    //        }
    //    }
    //
    //    // PRELIMINARY STEPS  ////////////////////////////////////////////////////////////
    //    // Initializing types for new pop-up
    //    for (int it=0; it<newPopUpNodes.size(); ++it) {
    //        // index of the popping up node
    //        unsigned int index=newPopUpNodes[it];
    //        // routine to check if it should be interface, fluid or gas
    //        bool outOfFluid=true;
    //        bool closeToGas=false;
    //        for (int j=1; j<lbmDirec; ++j) {
    //            unsigned int link=neighbors[index].d[j];
    //            if (types[link].isFluid()) {
    //                outOfFluid=false;
    //            }
    //            else if(types[link].isGas()) {
    //                closeToGas=true;
    //            }
    //        }
    //        if (outOfFluid) {
    //            //then node is fluid
    //            types[index].setGas();
    //        }
    //        if ((!outOfFluid)&&(closeToGas)) {
    //            //then node is interface
    //            types[index].setInterface();
    //        }
    //        if ((!outOfFluid)&&(!closeToGas)) {
    //            // then node is fluid
    //            types[index].setFluid();
    //        }
    //    }
    //
    //    // NEW POP-UP NODES - Initialization  ////////////////////////////////////////////////////////////
    //   for (int it=0; it<newPopUpNodes.size(); ++it) {
    //       unsigned int index=newPopUpNodes[it];
    //       // check if cell is active (it could be gas) - the type has already been assigned in the lines over here
    //        if (types[index].isActive()) {
    ////            cout<<nodes[index].x<<" "<<nodes[index].y<<" "<<nodes[index].z<<" NEW ACTIVE NODE from pop-up\n";
    //            //  macroscopic velocity and density for the new cell (average of neighbors)
    //            tVect uAverage;
    //            double nAverage=0.0;
    //            double viscAverage=0.0;
    //            double mass=0.0;
    //
    //            // Routine for average properties of new nodes. CHECK!!
    //            // number of neighbor active cells
    //            unsigned int ku=0;
    //            unsigned int kn=0;
    //            uAverage.reset();
    //            // calculating average
    //            for (int j=0; j<lbmDirec; j++) {
    //                unsigned int link=neighbors[index].d[j];
    //                // checking if neighbor is active
    //                if (types[link].isActive()) {
    //                    //we need to check if this cell is a new cell, in that case its properties are not good to use
    //                    bool isNew=false;
    //                    for (int it1=0; it1<newPopUpNodes.size(); ++it1) {
    ////                    for (intlist::iterator it1=newPopUpNodes.begin(); it1<newPopUpNodes.end(); ++it1) {
    //                        if (link==newPopUpNodes[it1]) {
    //                            isNew=true;
    //                            break;
    //                        }
    //                    }
    //                    if (!isNew) {
    //                         // incrementing velocity average
    //                        uAverage=uAverage+nodes[link]->u;
    //                        nAverage+=nodes[link]->n;
    //                        viscAverage+=nodes[link]->visc;
    //                        // incrementing number active neighbors
    //                        ++ku;
    //                        ++kn;
    //                    }
    //                }
    ////                if (nodes[link].isParticle()) {
    ////                    // incrementing velocity average
    ////                    // getting local velocity
    ////                    tVect radius=tVect(((double)nodes[link].x), ((double)nodes[link].y), ((double)nodes[link].z))-elmts[nodes[link].solidIndex].x0/unit.Length;
    ////                    // local velocity (v_local=v_center+omega x radius) (lattice units)
    ////                    tVect localVel=(elmts[nodes[link].solidIndex].x1/unit.Speed+(elmts[nodes[link].solidIndex].w.cross(radius))/unit.AngVel);
    ////                    uAverage=uAverage+localVel;
    ////                    // incrementing number of active neighbors
    ////                    ++ku;
    ////                }
    //            }
    //            if (ku==0) {
    //                uAverage.reset();
    //            }
    //            else {
    //                // dividing to get true average
    //                uAverage=uAverage/ku;
    //            }
    //            if (kn==0) {
    //                viscAverage=initDynVisc;
    //                nAverage=1.0;
    //            }
    //            else {
    //                // dividing to get true average
    //                nAverage/=kn;
    //                viscAverage/=kn;
    //            }
    //
    //            if (types[index].isInterface()) {
    ////                activeNodes.push_back(nodes[index].index);
    //                interfaceNodes.push_back(index);
    //                mass=0.01*nAverage;
    //            }
    //            else if (types[index].isFluid()) {
    //                fluidNodes.push_back(index);
    ////                activeNodes.push_back(nodes[index].index);
    //                mass=nAverage;
    //            }
    //            nodes[index]->initialize(nAverage, uAverage, mass, viscAverage, lbF);
    //            massSurplus-=mass;
    ////            initializeNode(nodes[index], nAverage, uAverage, mass, viscAverage);
    //            unsigned int zero=0;
    //            types[index].setSolidIndex(zero);
    //        }
    //        else {
    //            // remove node
    //            delete nodes[index];
    //            nodes[index]=0;
    //        }
    //    }
    
}

void LB::activeToSolid(unsIntList& newSolidNodes, elmtList& elmts, double& massSurplus) {
    
    //    for (int it=0; it<newSolidNodes.size(); ++it) {
    ////    for (intlist::iterator it=newSolidNodes.begin(); it<newSolidNodes.end(); ++it) {
    ////        cout<<"New solid\n";
    //        //adding cell to solid list
    //        particleNodes.push_back(newSolidNodes[it]);
    //        // deleting the cell from active lists (fluid or interface, depending on the type)
    //        if (types[newSolidNodes[it]].isFluid()) {
    //            for (int ind=fluidNodes.size()-1; ind>=0; --ind) {
    //                if (newSolidNodes[it]==fluidNodes[ind]) {
    //                    fluidNodes.erase(fluidNodes.begin()+ind);
    //                    break;
    //                }
    //            }
    //        }
    //        else if (types[newSolidNodes[it]].isInterface()) {
    //            for (int ind=interfaceNodes.size()-1; ind>=0; --ind) {
    //                if (newSolidNodes[it]==interfaceNodes[ind]) {
    //                    interfaceNodes.erase(interfaceNodes.begin()+ind);
    //                    break;
    //                }
    //            }
    //        }
    //    }
    //    // initializing type for new solid nodes and redistributing mass in case of interface
    //    for (int it=0; it<newSolidNodes.size(); ++it) {
    //        // index of new solid node
    //        unsigned int index=newSolidNodes[it];
    //        massSurplus+=nodes[index]->mass;
    //        // initializing type
    //        types[index].setParticle();
    //    }
    //
    //    // NEW SOLID NODES - Initialization ///////////////////////////////////////////////////////////
    //    for (int it=0; it<newSolidNodes.size(); ++it) {
    //        unsigned int index=newSolidNodes[it];
    ////        cout<<nodes[index].x<<" "<<nodes[index].y<<" "<<nodes[index].z<<" NEW SOLID NODE from absorption\n";
    //        // reset velocity and mass (useful for plotting)
    //        nodes[index]->n=0.0;
    //        nodes[index]->visc=0.0;
    //        nodes[index]->u.reset();
    //        nodes[index]->mass=0.0;
    //        nodes[index]->shearRate=0.0;
    ////        nodes[index].force.reset();
    //
    //        // solid index is added in the movement function
    //    }
    
}

// other functions

void LB::evolveBoundary() {
    // makes a boundary change, by switching gas and solid nodes
    
    unsIntList gassificationNodes;
    gassificationNodes.clear();
    if (time>0.0) {
        static double vel=(0.15/1.5)/unit.Speed;
        double height=1.0+vel*double(time);
        //    cout<<"height = "<<height<<"\t";
        
        for (int it=0; it<nodes.size(); ++it) {
            if (types[it].isWall()) {
                if ((getX(it)!=lbSize[2]-1)&&(getZ(it)!=0)&&
                        (getY(it)!=lbSize[1]-1)&&(getY(it)!=0)&&
                        (getX(it)!=lbSize[0]-1)&&(getX(it)!=0)) {
                    if (getZ(it)<height) {
                        gassificationNodes.push_back(it);
                    }
                }
            }
        }
        gassification(gassificationNodes);
    }
}

void LB::gassification(unsIntList& gassificationNodes) {
    unsIntList newInterfaceNodes;
    newInterfaceNodes.clear();
    //    cout<<"gasify "<<gassificationNodes.size()<<"nodes\n";
    
    for (int it =0; it<gassificationNodes.size(); ++it) {
        types[gassificationNodes[it]].setGas();
    }
    
    // NEW GAS NODES - Initialization  ////////////////////////////////////////////////////////////
    for (int it=0; it<gassificationNodes.size(); ++it) {
        unsigned int index=gassificationNodes[it];
        for (int j=0; j<lbmDirec; ++j) {
            unsigned int link=neighbors[index].d[j];
            if (types[link].isFluid()) {
                //                cout<<nodes[link].x<<" "<<nodes[link].y<<" "<<nodes[link].z<<" NEW INTERFACE NODE from fluid\n";
                types[link].setInterface();
                newInterfaceNodes.push_back(link);
                // characteristics are inherited by previous fluid cell. Only mass must be updated
                nodes[link]->mass=0.99*nodes[link]->n;
                // let's forget for the moment the mass update
                //                massSurplus+=0.01*nodes[link].n;
                for (int ind=fluidNodes.size()-1; ind>=0; --ind) {
                    if (link==fluidNodes[ind]) {
                        fluidNodes.erase(fluidNodes.begin()+ind);
                        break;
                    }
                }
            }
        }
    }
    // initializing new interface nodes type
    for (int it=0; it<newInterfaceNodes.size(); ++it) {
        // isn't this already done? check!
        types[newInterfaceNodes[it]].setInterface();
    }
    
    // resetting new gas macroscopic quantities and updating neighbors for gas cells
    for (int it=0; it<gassificationNodes.size(); ++it) {
        bool isNotInterface=true;
        for (int it1=0; it1<newInterfaceNodes.size(); ++it1) {
            if (gassificationNodes[it]==newInterfaceNodes[it1]) {
                isNotInterface=false;
            }
        }
        if (isNotInterface) {
            // delete memory
            delete nodes[gassificationNodes[it]];
            nodes[gassificationNodes[it]]=0;
        }
    }
    // resetting neighbors for newly created interface cells
    for (int it=0; it<newInterfaceNodes.size(); ++it) {
        interfaceNodes.push_back(newInterfaceNodes[it]);
    }
    //    // distributing surplus to interface cells
    //    double addMass=massSurplus/interfaceNodes.size();
    //    for (int it=0; it<newInterfaceNodes.size(); ++it) {
    //        nodes[newInterfaceNodes[it]].massFunction+=addMass;
    //    }
}

// functions for index management

tVect LB::setPosition(const unsigned int& index) const {
    unsigned int x,y,z;
    
    // index is calculated in this fashion:
    // index = x + y*X + z*X*Y
    // where X and Y are sizes of the lattice in x and y direction
    
    // from this stems that
    // x + y*X = index MOD X*Y
    // x = x + y*X MOD X
    // y = x + y*X DIV X
    // z = index DIV X*Y
    
    // see online documentation for class div_t (stdlib.h)
    div_t firstDiv, secondDiv;
    
    firstDiv=div(int(index), int(lbSize[0]*lbSize[1]));
    secondDiv=div(firstDiv.rem, int(lbSize[0]));
    
    x=secondDiv.rem;
    y=secondDiv.quot;
    z=firstDiv.quot;
    
    return tVect(double(x), double(y), double(z));
}

unsigned int LB::getX(const unsigned int& index) const {
    
    // see function getPosition for documentation
    div_t firstDiv, secondDiv;
    
    firstDiv=div(int(index), int(lbSize[0]*lbSize[1]));
    secondDiv=div(firstDiv.rem, int(lbSize[0]));
    
    return secondDiv.rem;
}

unsigned int LB::getY(const unsigned int& index) const {
    
    // see function getPosition for documentation
    div_t firstDiv, secondDiv;
    
    firstDiv=div(int(index), int(lbSize[0]*lbSize[1]));
    secondDiv=div(firstDiv.rem, int(lbSize[0]));
    
    return secondDiv.quot;
}

unsigned int LB::getZ(const unsigned int& index) const {
    
    // see function getPosition for documentation
    div_t firstDiv;
    
    firstDiv=div(int(index), int(lbSize[0]*lbSize[1]));
    
    return firstDiv.quot;
}