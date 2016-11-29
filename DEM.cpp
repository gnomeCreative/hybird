
#include "DEM.h"
#include "IO.h"

#define USE_MATH_DEFINES

using namespace std;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

void DEM::discreteElementGet(GetPot& lbmCfgFile, GetPot& command_line) {
    // getting material properties
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.density, "density", 1.0);
    ASSERT(sphereMat.density > 0.0);

    string contactModelString;
    PARSE_CLASS_MEMBER(lbmCfgFile, contactModelString, "contactModel", "none");
    if (contactModelString == "HERTZIAN") sphereMat.contactModel = HERTZIAN;
    else if (contactModelString == "LINEAR") sphereMat.contactModel = LINEAR;


    // Hertzian contact model /////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.youngMod, "youngMod", 1.0);
    ASSERT(sphereMat.youngMod > 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.poisson, "poisson", 1.0);
    ASSERT(sphereMat.poisson >= 0.0);
    // normal stiffness (constant part, calculated here to avoid repetition)
    sphereMat.knConst = 2.0 / 3.0 * sphereMat.youngMod / (1 - sphereMat.poisson * sphereMat.poisson);
    // there was a bracket mistake here
    sphereMat.ksConst = 2.0 * sphereMat.youngMod / (2.0 - sphereMat.poisson) / (1.0 + sphereMat.poisson);

    // linear contact model /////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.linearStiff, "linearStiff", 1.0);
    ASSERT(sphereMat.linearStiff >= 0.0);

    // normal damping ///////////////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.restitution, "restitution", 1.0);
    ASSERT(sphereMat.restitution > 0.0);
    ASSERT(sphereMat.restitution <= 1.0);
    // calculating coefficient for normal damping
    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            // see "On the Determination of the Damping Coefficient of Non-linear Spring-dashpot System to Model Hertz Contact for Simulation by Discrete Element Method"
            // Hu, Hu, Jian, Liu, Wan, Journal of Computers, 6 (2011) OR BETTER Antypov & Elliott
            sphereMat.dampCoeff = -1.0 * sqrt(5) * log(sphereMat.restitution) / sqrt((log(sphereMat.restitution) * log(sphereMat.restitution) + M_PI * M_PI));
        }
        case LINEAR:
        {
            sphereMat.dampCoeff = -1.0 * log(sphereMat.restitution) / sqrt((log(sphereMat.restitution) * log(sphereMat.restitution) + M_PI * M_PI));
        }
    }

    // tangential model //////////////////////
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.viscTang, "viscTang", 1.0);
    ASSERT(sphereMat.viscTang >= 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.frictionCoefPart, "frictionCoefPart", 1.0);
    ASSERT(sphereMat.frictionCoefPart >= 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, sphereMat.frictionCoefWall, "frictionCoefWall", 1.0);
    ASSERT(sphereMat.frictionCoefWall >= 0.0);

    // particle initial state //////////////////////
    string particleFile;
    PARSE_CLASS_MEMBER(lbmCfgFile, particleFile, "particleFile", "particles.dat");
    double translateX(0.0), translateY(0.0), translateZ(0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, translateX, "translateX", 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, translateY, "translateY", 0.0);
    PARSE_CLASS_MEMBER(lbmCfgFile, translateZ, "translateZ", 0.0);
    tVect translate(translateX, translateY, translateZ);
    double scale = 1.0;
    PARSE_CLASS_MEMBER(lbmCfgFile, scale, "scale", 1.0);

    ifstream particleFileID;
    particleFileID.open(particleFile.c_str(), ios::in);
    ASSERT(particleFileID.is_open());

    unsigned int totElmt;
    particleFileID>>totElmt;

    for (int n = 0; n < totElmt; ++n) {
        elmt dummyElmt;

        // import variables
        particleFileID >> dummyElmt.index;
        particleFileID >> dummyElmt.size;
        particleFileID >> dummyElmt.radius;
        dummyElmt.radius = dummyElmt.radius*scale;
        // position
        double x0, y0, z0;
        particleFileID>>x0;
        particleFileID>>y0;
        particleFileID>>z0;
        dummyElmt.x0 = tVect(x0*scale, y0*scale, z0 * scale) + translate;
        // translational velocity
        double x1, y1, z1;
        particleFileID>>x1;
        particleFileID>>y1;
        particleFileID>>z1;
        dummyElmt.x1 = tVect(x1, y1, z1); //.reset();
        // rotational velocity
        double w01, w02, w03;
        particleFileID>>w01;
        particleFileID>>w02;
        particleFileID>>w03;
        dummyElmt.w0 = tVect(w01, w02, w03); //.reset();
        //dummyElmt.x1.reset();
        // orientation
        double p0, q0, r0, s0;
        particleFileID>>p0;
        particleFileID>>q0;
        particleFileID>>r0;
        particleFileID>>s0;
        dummyElmt.q0 = tQuat(p0, q0, r0, s0);
        //dummyElmt.q0.resetSoft();
        // translational velocity (in quaternion rates))
        double p1, q1, r1, s1;
        particleFileID>>p1;
        particleFileID>>q1;
        particleFileID>>r1;
        particleFileID>>s1;
        dummyElmt.q1 = tQuat(p1, q1, r1, s1);
        //dummyElmt.q1.resetHard();
        // add to list
        elmts.push_back(dummyElmt);

    }

    // objects initial state //////////////////////
    string objectFile;
    PARSE_CLASS_MEMBER(lbmCfgFile, objectFile, "objectFile", "objects.dat");
    ifstream objectFileID;
    objectFileID.open(objectFile.c_str(), ios::in);
    ASSERT(objectFileID.is_open());

    unsigned int totObjects;
    objectFileID>>totObjects;

    for (int n = 0; n < totObjects; ++n) {
        object dummyObject;

        double size; //dummy, just to use same style

        // import variables
        objectFileID >> dummyObject.index;
        objectFileID>>size; // must be one
        objectFileID >> dummyObject.r;
        double x0, y0, z0;
        objectFileID>>x0;
        objectFileID>>y0;
        objectFileID>>z0;
        dummyObject.x0 = tVect(x0, y0, z0);
        double x1, y1, z1;
        objectFileID>>x1;
        objectFileID>>y1;
        objectFileID>>z1;
        dummyObject.x1 = tVect(x1, y1, z1);
        // the next eight values are for rotation, and are not used
        double trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;
        objectFileID>>trash;

        objects.push_back(dummyObject);
    }


    // numerical viscosity for stability
    PARSE_CLASS_MEMBER(lbmCfgFile, numVisc, "numVisc", 0.0);
    // set multiplication number (sets the number of DEM steps between two LBM steps)
    PARSE_CLASS_MEMBER(lbmCfgFile, multiStep, "multiStep", 1);
    // set ratio between time step and estimated duration of contacts (only if multiStep=0)
    PARSE_CLASS_MEMBER(lbmCfgFile, criticalRatio, "criticalRatio", 0.1);

}

void DEM::discreteElementInit(unsIntList& boundary, unsIntList& lbSize, measureUnits& unit, tVect& lbF) {

    // initializing DEM parameters from lattice parameters

    // domain size is the same of LBM. This is expressed in lattice coordinates
    demSize[0] = double(lbSize[0]) * unit.Length;
    demSize[1] = double(lbSize[1]) * unit.Length;
    demSize[2] = double(lbSize[2]) * unit.Length;
    // acceleration field
    demF = lbF * unit.Accel;
    // capping distance for lubrication forces
    minDistLub = 0.01 * unit.Length;

    // DEM time step
    // if multistep is 0, it should be calculated by the program here
    if (elmts.size()) {
        if (multiStep == 0) {
            // find critical deltaT
            const double crit = criticalRatio * criticalTimeStep();
            // if the critical time is bigger thatn the LBM time step, then just use the LBM time step
            if (crit >= unit.Time) {
                multiStep = 1;
                deltat = unit.Time;
            }                // if it is lower, than calculate the number of substeps
            else {
                const double ratio = unit.Time / crit;
                ASSERT(ratio >= 1);
                multiStep = std::floor(ratio) + 1;
                deltat = unit.Time / (double(multiStep));
            }
        }            // multistep can also be imposed by the user
        else {
            deltat = unit.Time / (double(multiStep));
        }
    } else {
        deltat = unit.Time;
        multiStep = 1;
    }


    // initializing particles
    const double partDensity = sphereMat.density;

    // initializing composite particle properties
    compositeProperties();
    // clear particle list
    particles.clear();

    unsigned int globalIndex = 0;
    for (int n = 0; n < elmts.size(); ++n) {
        // initialize element
        elmts[n].initialize(partDensity, prototypes, demF);
        // generate particles
        elmts[n].generateParticles(globalIndex, particles, prototypes);
    }
    // the number of standard particles (=not ghosts) is now fixed, and WILL NOT BE CHANGED
    stdParticles = particles.size();

    // initializing wall for DEM
    initializeWalls(boundary, lbSize, unit);
    // initializing cylinders for DEM
    initializeCylinders();
    // initializing periodic boundary conditions for DEM
    initializePbcs(boundary, lbSize, unit);

    // used to get tot number of objects for CG (coarse graining)
    stdObjects = objects.size();

    // initialize neighbor list parameters (also sets ghosts)
    if (elmts.size()) {
        initNeighborParameters();
        periodicObjects();
        evalNeighborTable();
    }

    double totMass = 0.0;
    for (int n = 0; n < elmts.size(); ++n) {
        // calculate mass
        totMass += elmts[n].m;
    }

    cout << "DEM parameters\n";
    cout << "domain size: xdim =" << demSize[0] << "; ydim= " << demSize[1] << "; zdim= " << demSize[2] << ";\n";
    cout << "Tot elements: " << elmts.size() << ";\t";
    cout << "Tot objects: " << objects.size() << ";\t";
    cout << "Tot standard particles: " << stdParticles << endl;
    if (multiStep > 0) {
        cout << "Deltat =" << deltat << ", therefore " << multiStep << " multiple steps, as a " << criticalRatio << ":1 ratio to estimated collision time" << endl;
    } else {
        cout << "Deltat =" << deltat << ", by imposing " << multiStep << " substeps" << endl;
    }

    cout << "Total mass =" << totMass << endl;
    switch (sphereMat.contactModel) {
        case LINEAR:
        {
            cout << "Contact model: linear dashpot" << endl;
            cout << "Normal stiffness = " << sphereMat.linearStiff << endl;
            break;
        }
        case HERTZIAN:
        {
            cout << "Contact model: damped Hertzian contact" << endl;
            break;
        }
    }
    cout << "Damping ratio = " << sphereMat.dampCoeff << ", equivalent to a coefficient of restitution of c=" << sphereMat.restitution << endl;
    cout << "Tangential viscosity = " << sphereMat.viscTang << endl;
    cout << "Particle-particle friction = " << sphereMat.frictionCoefPart << ", wall-particle friction = " << sphereMat.frictionCoefWall << endl;
    cout << "Numerical viscosity =" << numVisc << endl;
}

void DEM::discreteElementStep(IO& io) {

    // set trigger for new neighbor list
    double neighListTrigger = 0.25 * nebrRange;

    for (int demIter = 0; demIter < multiStep; ++demIter) {

        // neighbor management
        evalMaxDisp();
        if (maxDisp > neighListTrigger) { // maxDisp>0.25*(nebrRange-2.0*cutOff)
            maxDisp = 0.0;
            //cout<<"new neighbor list"<<endl;
            evalNeighborTable();
        }

        // predictor step
        predictor();

        // particles generation
        updateParticlesPredicted();

        // the statement below makes sure that .restart and .data are saved every multiple of saveCount
        if (io.saveCount != 0 && io.currentTimeStep % io.saveCount == 0) {
            // export .data and .restart file
            exportDataFile(io);
            exportRestartFile(io);
        }

        // force evaluation
        evaluateForces(io);

        // force evaluation (old version))
        //evalForces();

        // corrector step
        corrector();

        // particles re-generation
        updateParticlesCorrected();

        // updating energy
        updateEnergy();

        //delete particle particle elongation to solve non collision particles
        deleteParticleElongation();

        //delete wall particle elongation to solve non collision
        deleteWallElongation();

    }

}

void DEM::evolveBoundaries() {

    double frac = 1.0 / demInitialRepeat / ((double) multiStep);
    //    cout<<"frac="<<frac<<"\n";

    for (int n = 0; n < walls.size(); ++n) {

        if (walls[n].translating) {
            walls[n].p = walls[n].p + frac * walls[n].trans;

            if (n == 1) {
                //            walls[n].p.show();
                //            cout<<"\n";
            }
        }
    }

    evalNeighborTable();
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// initialization functions

void DEM::compositeProperties() {

    vecList prototype1, prototype2, prototype3, prototype4;


    // prototypes for shapes
    // every vector defines the position of a particle in the object reference frame
    // unit is radius

    prototypes.resize(5);
    prototype1.resize(1);
    prototype1[0].reset();
    prototypes[1] = prototype1;
    prototype2.resize(2);
    prototype2[0] = tVect(0.5, 0.0, 0.0);
    prototype2[1] = tVect(-0.5, 0.0, 0.0);
    prototypes[2] = prototype2;
    prototype3.resize(3);
    prototype3[0] = tVect(0.0, 1.0, 0.0);
    prototype3[1] = tVect(-sqrt(3) / 2, -1 / 2, 0.0);
    prototype3[2] = tVect(sqrt(3) / 2, -1 / 2, 0.0);
    prototypes[3] = prototype3;
    prototype4.resize(4);
    prototype4[0] = tVect(0.0, 0.0, 1.0);
    prototype4[1] = tVect(0.0, 2.0 * sqrt(2) / 3.0, -1.0 / 3.0);
    prototype4[2] = tVect(2.0 * sqrt(6) / 6.0, -2.0 * sqrt(2) / 6.0, -1.0 / 3.0);
    prototype4[3] = tVect(-2.0 * sqrt(6) / 6.0, -2.0 * sqrt(2) / 6.0, -1.0 / 3.0);
    prototypes[4] = prototype4;

}

void DEM::initializeWalls(unsIntList& boundary, unsIntList& lbSize, measureUnits& unit) {
    walls.clear();
    wall dummyWall;
    unsigned int index = 0;
    double boxVel = 0.15 / 1.5;

    // wall #0
    if ((boundary[0] == 5) || (boundary[0] == 6) || (boundary[0] == 7) || (boundary[0] == 8)) {
        dummyWall.p = tVect(0.5 * unit.Length, 0.0, 0.0);
        dummyWall.n = tVect(1.0, 0.0, 0.0);
        dummyWall.index = index;
        dummyWall.flag = false;
        dummyWall.translating = false;
        dummyWall.trans.reset();
        if (boundary[0] == 5) {
            dummyWall.moving = false;
            dummyWall.slip = true;
        } else if (boundary[0] == 6) {
            dummyWall.moving = true;
            dummyWall.slip = true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        } else if (boundary[0] == 7) {
            dummyWall.moving = false;
            dummyWall.slip = false;
        } else if (boundary[0] == 8) {
            dummyWall.moving = true;
            dummyWall.slip = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel = tVect(0.00 * unit.Speed, 0.0, 0.0); // HERE
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #1
    if ((boundary[1] == 5) || (boundary[1] == 6) || (boundary[1] == 7) || (boundary[1] == 8)) {
        dummyWall.p = tVect((double(lbSize[0]) - 1.5) * unit.Length, 0.0, 0.0);
        dummyWall.n = tVect(-1.0, 0.0, 0.0);
        dummyWall.index = index;
        dummyWall.flag = false;
        dummyWall.translating = false;
        dummyWall.trans.reset();
        if (boundary[1] == 5) {
            dummyWall.moving = false;
            dummyWall.slip = true;
        } else if (boundary[1] == 6) {
            dummyWall.moving = true;
            dummyWall.slip = true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        } else if (boundary[1] == 7) {
            dummyWall.moving = false;
            dummyWall.slip = false;
        } else if (boundary[1] == 8) {
            dummyWall.moving = true;
            dummyWall.slip = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel = tVect(0.00 * unit.Speed, 0.0, 0.0); // HERE
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #2
    if ((boundary[2] == 5) || (boundary[2] == 6) || (boundary[2] == 7) || (boundary[2] == 8)) {
        dummyWall.p = tVect(0.0, 0.5 * unit.Length, 0.0);
        dummyWall.n = tVect(0.0, 1.0, 0.0);
        dummyWall.index = index;
        dummyWall.flag = false;
        dummyWall.translating = false;
        dummyWall.trans.reset();
        if (boundary[2] == 5) {
            dummyWall.moving = false;
            dummyWall.slip = true;
        } else if (boundary[2] == 6) {
            dummyWall.moving = true;
            dummyWall.slip = true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        } else if (boundary[2] == 7) {
            dummyWall.moving = false;
            dummyWall.slip = false;
        } else if (boundary[2] == 8) {
            dummyWall.moving = true;
            dummyWall.slip = false;
            dummyWall.rotCenter = tVect(0.0, 0.5 * unit.Length + 0.4, 1.260);
            dummyWall.omega = tVect(0.0, 0.0, 0.0); //0.942
            dummyWall.vel = tVect(-0.0, 0.0, 0.0);
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #3
    if ((boundary[3] == 5) || (boundary[3] == 6) || (boundary[3] == 7) || (boundary[3] == 8)) {
        dummyWall.p = tVect(0.0, (double(lbSize[1]) - 1.5) * unit.Length, 0.0);
        dummyWall.n = tVect(0.0, -1.0, 0.0);
        dummyWall.index = index;
        dummyWall.flag = false;
        dummyWall.translating = false;
        dummyWall.trans.reset();
        if (boundary[3] == 5) {
            dummyWall.moving = false;
            dummyWall.slip = true;
        } else if (boundary[3] == 6) {
            dummyWall.moving = true;
            dummyWall.slip = true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        } else if (boundary[3] == 7) {
            dummyWall.moving = false;
            dummyWall.slip = false;
        } else if (boundary[3] == 8) {
            dummyWall.moving = true;
            dummyWall.slip = false;
            dummyWall.rotCenter = tVect(0.0, (lbSize[1] - 1.5) * unit.Length + 0.4, 1.260);
            dummyWall.omega = tVect(0.0, 0.0, 0.0); //0.942
            dummyWall.vel = tVect(-0.0, 0.0, 0.0);
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #4
    if ((boundary[4] == 5) || (boundary[4] == 6) || (boundary[4] == 7) || (boundary[4] == 8)) {
        dummyWall.p = tVect(0.0, 0.0, 0.5 * unit.Length);
        dummyWall.n = tVect(0.0, 0.0, 1.0);
        dummyWall.index = index;
        dummyWall.flag = false;
        dummyWall.translating = false;
        dummyWall.trans.reset();
        if (boundary[4] == 5) {
            dummyWall.moving = false;
            dummyWall.slip = true;
        } else if (boundary[4] == 6) {
            dummyWall.moving = true;
            dummyWall.slip = true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        } else if (boundary[4] == 7) {
            dummyWall.moving = false;
            dummyWall.slip = false;
        } else if (boundary[4] == 8) {
            dummyWall.moving = true;
            dummyWall.slip = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel = tVect(0.0, 0.0, 0.0); // HERE (-1.0*)wallVel
        }
        ++index;
        walls.push_back(dummyWall);
    }
    // wall #5
    if ((boundary[5] == 5) || (boundary[5] == 6) || (boundary[5] == 7) || (boundary[5] == 8)) {
        dummyWall.p = tVect(0.0, 0.0, (double(lbSize[2]) - 1.5) * unit.Length);
        dummyWall.n = tVect(0.0, 0.0, -1.0);
        dummyWall.index = index;
        dummyWall.flag = false;
        dummyWall.translating = true; // CHANGE HERE
        dummyWall.trans = tVect(0.0, 0.0, -150.0 * unit.Length); // CHANGE HERE
        if (boundary[5] == 5) {
            dummyWall.moving = false;
            dummyWall.slip = true;
        } else if (boundary[5] == 6) {
            dummyWall.moving = true;
            dummyWall.slip = true;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
        } else if (boundary[5] == 7) {
            dummyWall.moving = false;
            dummyWall.slip = false;
        } else if (boundary[5] == 8) {
            dummyWall.moving = true;
            dummyWall.slip = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel = tVect(0.0, 0.0, 0.0); // HERE  wallVel
        }
        ++index;
        walls.push_back(dummyWall);
    }

    // additional walls
    switch (problemName) {
        case AVALANCHE:
        {
            dummyWall.p = tVect(60, 15, 15 - 14.86);
            dummyWall.n = tVect(0.42262, 0.0, 0.9063);
            dummyWall.index = index;
            dummyWall.flag = false;
            dummyWall.moving = false;
            dummyWall.rotCenter.reset();
            dummyWall.omega.reset();
            dummyWall.vel.reset();
            dummyWall.translating = false;
            ++index;
            walls.push_back(dummyWall);

            //            dummyWall.p=tVect((double(lbSize[0])-1.5)*unit.Length,0.0,35.0*unit.Length);
            //            dummyWall.n=tVect(-0.087156,0.0,0.996195);
            //            dummyWall.index=index;
            //            dummyWall.flag=false;
            //            dummyWall.moving=false;
            //            dummyWall.rotCenter.reset();
            //            dummyWall.omega.reset();
            //            dummyWall.vel.reset();
            //            dummyWall.translating=false;
            //            ++index;
            //            walls.push_back(dummyWall);
            break;
        }
        case BOX:
        {
            //            // xMin
            //            dummyWall.p=Tvect((30.0+0.5)*unit.Length,0.0,0.0);
            //            dummyWall.n=Tvect(1.0,0.0,0.0);
            //            dummyWall.index=index;
            //            dummyWall.flag=true;
            //            dummyWall.moving=true;
            //            dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //            dummyWall.omega=Tvect(0.0,0.0,0.0);
            //            dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //            dummyWall.translating=true;
            //            dummyWall.trans=Tvect((20.0)*unit.Length,0.0,0.0);
            //            ++index;
            //            walls.push_back(dummyWall);
            //
            //        //    dummyWall.p=Tvect((50.0-1.5)*unit.Length,0.0,0.0);
            //        //    dummyWall.n=Tvect(-1.0,0.0,0.0);
            //        //    dummyWall.index=index;
            //        //    dummyWall.flag=true;
            //        //    dummyWall.moving=true;
            //        //    dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.omega=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //        //    ++index;
            //        //    walls.push_back(dummyWall);
            //
            //            // xMax
            //        //    dummyWall.p=Tvect((125.0-0.5)*unit.Length,0.0,0.0);
            //        //    dummyWall.n=Tvect(1.0,0.0,0.0);
            //        //    dummyWall.index=index;
            //        //    dummyWall.flag=true;
            //        //    dummyWall.moving=true;
            //        //    dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.omega=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //        //    ++index;
            //        //    walls.push_back(dummyWall);
            //
            //            dummyWall.p=Tvect((145.0-0.5)*unit.Length,0.0,0.0);
            //            dummyWall.n=Tvect(-1.0,0.0,0.0);
            //            dummyWall.index=index;
            //            dummyWall.flag=true;
            //            dummyWall.moving=true;
            //            dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //            dummyWall.omega=Tvect(0.0,0.0,0.0);
            //            dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //            dummyWall.translating=true;
            //            dummyWall.trans=Tvect((-20.0)*unit.Length,0.0,0.0);
            //            ++index;
            //            walls.push_back(dummyWall);
            //
            //            // yMin
            //            dummyWall.p=Tvect(0.0,(68.0+0.5)*unit.Length,0.0);
            //            dummyWall.n=Tvect(0.0,1.0,0.0);
            //            dummyWall.index=index;
            //            dummyWall.flag=true;
            //            dummyWall.moving=true;
            //            dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //            dummyWall.omega=Tvect(0.0,0.0,0.0);
            //            dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //            dummyWall.translating=true;
            //            dummyWall.trans=Tvect(0.0,(20.0)*unit.Length,0.0);
            //            ++index;
            //            walls.push_back(dummyWall);
            //
            //        //    dummyWall.p=Tvect(0.0,(88.0-1.5)*unit.Length,0.0);
            //        //    dummyWall.n=Tvect(0.0,-1.0,0.0);
            //        //    dummyWall.index=index;
            //        //    dummyWall.flag=true;
            //        //    dummyWall.moving=true;
            //        //    dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.omega=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //        //    ++index;
            //        //    walls.push_back(dummyWall);
            //
            //            // yMax
            //        //    dummyWall.p=Tvect(0.0,(162.0-0.5)*unit.Length,0.0);
            //        //    dummyWall.n=Tvect(0.0,1.0,0.0);
            //        //    dummyWall.index=index;
            //        //    dummyWall.flag=true;
            //        //    dummyWall.moving=true;
            //        //    dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.omega=Tvect(0.0,0.0,0.0);
            //        //    dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //        //    ++index;
            //        //    walls.push_back(dummyWall);
            //
            //            dummyWall.p=Tvect(0.0,(182.0-0.5)*unit.Length,0.0);
            //            dummyWall.n=Tvect(0.0,-1.0,0.0);
            //            dummyWall.index=index;
            //            dummyWall.flag=true;
            //            dummyWall.moving=true;
            //            dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //            dummyWall.omega=Tvect(0.0,0.0,0.0);
            //            dummyWall.vel=Tvect(0.0,0.0,boxVel);
            //            dummyWall.translating=true;
            //            dummyWall.trans=Tvect(0.0,(-20.0)*unit.Length,0.0);
            //            ++index;
            //            walls.push_back(dummyWall);
            //
            //            dummyWall.p=Tvect(0.0,0.0,(70.0-0.5)*unit.Length);
            //            dummyWall.n=Tvect(0.0,0.0,-1.0);
            //            dummyWall.index=index;
            //            dummyWall.flag=true;
            //            dummyWall.moving=true;
            //            dummyWall.rotCenter=Tvect(0.0,0.0,0.0);
            //            dummyWall.omega=Tvect(0.0,0.0,0.0);
            //            dummyWall.vel=Tvect(0.0,0.0,0.0);
            //            dummyWall.translating=true;
            //            dummyWall.trans=Tvect(0.0,0.0,0.0);
            //            ++index;
            //            walls.push_back(dummyWall);
            break;
        }
        case DRUM:
        {
            // left wall
            if ((boundary[2] == 5) || (boundary[2] == 6) || (boundary[2] == 7) || (boundary[2] == 8)) {
                dummyWall.p = tVect(0.0 + 0.6 + 0.35, 0.50001 * unit.Length, 0.0);
                dummyWall.n = tVect(0.0, 1.0, 0.0);
                dummyWall.index = index;
                dummyWall.flag = false;
                dummyWall.translating = false;
                dummyWall.trans.reset();
                if (boundary[2] == 5) {
                    dummyWall.moving = false;
                    dummyWall.slip = true;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                } else if (boundary[2] == 6) {
                    dummyWall.moving = true;
                    dummyWall.slip = true;
                    dummyWall.rotCenter = tVect(0.0 + 0.6 + 0.35, 0.50001 * unit.Length, 1.260);
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                    dummyWall.vel = tVect(-0.0, 0.0, 0.0);
                } else if (boundary[2] == 7) {
                    dummyWall.moving = false;
                    dummyWall.slip = false;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                } else if (boundary[2] == 8) {
                    dummyWall.moving = true;
                    dummyWall.slip = false;
                    dummyWall.rotCenter = tVect(0.0 + 0.6 + 0.35, 0.50001 * unit.Length, 1.260);
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                    dummyWall.vel = tVect(-0.0, 0.0, 0.0);
                }
                ++index;
                walls.push_back(dummyWall);
            }
            // right wall
            if ((boundary[3] == 5) || (boundary[3] == 6) || (boundary[3] == 7) || (boundary[3] == 8)) {
                dummyWall.p = tVect(0.0 + 0.6 + 0.35, (double(lbSize[1]) - 1.50001) * unit.Length, 0.0);
                dummyWall.n = tVect(0.0, -1.0, 0.0);
                dummyWall.index = index;
                dummyWall.flag = false;
                dummyWall.translating = false;
                dummyWall.trans.reset();
                if (boundary[3] == 5) {
                    dummyWall.moving = false;
                    dummyWall.slip = true;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                } else if (boundary[3] == 6) {
                    dummyWall.moving = true;
                    dummyWall.slip = true;
                    dummyWall.rotCenter = tVect(0.0 + 0.6 + 0.35, (lbSize[1] - 1.50001) * unit.Length, 1.260);
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                    dummyWall.vel = tVect(-0.0, 0.0, 0.0);
                } else if (boundary[3] == 7) {
                    dummyWall.moving = false;
                    dummyWall.slip = false;
                    dummyWall.rotCenter.reset();
                    dummyWall.omega.reset();
                    dummyWall.vel.reset();
                } else if (boundary[3] == 8) {
                    dummyWall.moving = true;
                    dummyWall.slip = false;
                    dummyWall.rotCenter = tVect(0.0 + 0.6 + 0.35, (lbSize[1] - 1.50001) * unit.Length, 1.260);
                    dummyWall.omega = tVect(0.0, drumSpeed, 0.0);
                    dummyWall.vel = tVect(-0.0, 0.0, 0.0);
                }
                ++index;
                walls.push_back(dummyWall);
            }
            break;
        }
    }

    for (int n = 0; n < walls.size(); ++n) {
        walls[n].wallShow();
    }
}

void DEM::initializeCylinders() {
    cylinders.clear();
    unsigned int index = 0;

    switch (problemName) {
        case DRUM:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0.0 + 0.6 + 0.35, 0.0, 1.26);
            dummyCylinder.p2 = tVect(0.0 + 0.6 + 0.35, 1.0, 1.26);
            dummyCylinder.R = 1.243;
            dummyCylinder.omega = tVect(0.0, drumSpeed, 0.0);
            dummyCylinder.initAxes();
            dummyCylinder.moving = true;
            dummyCylinder.slip = false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
        case AVALANCHE:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0.0, 15, 15); //tVect(0.0,14.86,14.86);
            dummyCylinder.p2 = tVect(1.0, 15, 15);
            dummyCylinder.R = 14.86;
            dummyCylinder.omega.reset();
            dummyCylinder.initAxes();
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
            //            cylinder dummyCylinder;
            //            dummyCylinder.index=index;
            //            dummyCylinder.p1=tVect(0.0,10.0,14.86);
            //            dummyCylinder.p2=tVect(1.0,10.0,14.86);
            //            dummyCylinder.R=14.86;
            //            dummyCylinder.omega.reset();
            //            dummyCylinder.initAxes();
            //            dummyCylinder.moving=false;
            //            dummyCylinder.slip=false;
            //            cylinders.push_back(dummyCylinder);
            //            ++index;
            //            break;
        }
        case NET:
        {
            cylinder dummyCylinder;
            dummyCylinder.index = index;
            dummyCylinder.p1 = tVect(0.0, 4.0, 4.0);
            dummyCylinder.p2 = tVect(1.0, 4.0, 4.0);
            dummyCylinder.R = 4.0;
            dummyCylinder.omega.reset();
            dummyCylinder.initAxes();
            dummyCylinder.moving = false;
            dummyCylinder.slip = false;
            cylinders.push_back(dummyCylinder);
            ++index;
            break;
        }
    }

    for (int n = 0; n < cylinders.size(); ++n) {
        cylinders[n].cylinderShow();
    }
}

void DEM::initializePbcs(unsIntList& boundary, unsIntList& lbSize, measureUnits& unit) {

    pbcs.clear();
    pbc dummyPbc;
    unsigned int index = 0;

    if (boundary[0] == 4) {
        if (boundary[1] != 4) {
            cout << "Periodic boundaries error!" << endl;
            exit(0);
        } else {
            dummyPbc.p = tVect(0.5 * unit.Length, 0.0, 0.0);
            dummyPbc.v = tVect((double(lbSize[0]) - 2.0) * unit.Length, 0.0, 0.0);
            dummyPbc.setPlanes();
            dummyPbc.index = index;
            ++index;
            pbcs.push_back(dummyPbc);
        }
    }
    if (boundary[2] == 4) {
        if (boundary[3] != 4) {
            cout << "Periodic boundaries error!" << endl;
            exit(0);
        } else {
            dummyPbc.p = tVect(0.0, 0.5 * unit.Length, 0.0);
            dummyPbc.v = tVect(0.0, (double(lbSize[1]) - 2.0) * unit.Length, 0.0);
            dummyPbc.setPlanes();
            dummyPbc.index = index;
            ++index;
            pbcs.push_back(dummyPbc);
        }
    }
    if (boundary[4] == 4) {
        if (boundary[5] != 4) {
            cout << "Periodic boundaries error!" << endl;
            exit(0);
        } else {
            dummyPbc.p = tVect(0.0, 0.0, 0.5 * unit.Length);
            dummyPbc.v = tVect(0.0, 0.0, (double(lbSize[2]) - 2.0) * unit.Length);
            dummyPbc.setPlanes();
            dummyPbc.index = index;
            ++index;
            pbcs.push_back(dummyPbc);
        }
    }
    for (int b = 0; b < pbcs.size(); ++b) {
        pbcs[b].pbcShow();
    }
}

void DEM::periodicObjects() {

    // copy objects
    unsigned int originalObjects = objects.size();
    unsigned int objectIndex = objects.size() - 1;
    for (int b = 0; b < pbcs.size(); ++b) {
        const tVect pbcVector = pbcs[b].v;
        // cycle through elements
        for (int o = 0; o < objects.size(); ++o) {
            const double leftDist = pbcs[b].pl1.dist(objects[o].x0);
            const double rightDist = pbcs[b].pl2.dist(objects[o].x0);
            if (leftDist < 0.0) {
                objectIndex++;
                object dummyObject = objects[o];
                dummyObject.x0 = dummyObject.x0 + pbcVector;
                dummyObject.index = objectIndex;
                objects.push_back(dummyObject);
            } else if (leftDist < 2.0 * nebrRange && o < originalObjects) {
                objectIndex++;
                object dummyObject = objects[o];
                dummyObject.x0 = dummyObject.x0 + pbcVector;
                dummyObject.index = objectIndex;
                objects.push_back(dummyObject);
            }
            // second plane of couple (right)
            if (rightDist < 0.0) {
                objectIndex++;
                object dummyObject = objects[o];
                dummyObject.x0 = dummyObject.x0 - pbcVector;
                dummyObject.index = objectIndex;
                objects.push_back(dummyObject);
            } else if (rightDist < 2.0 * nebrRange && o < originalObjects) {
                objectIndex++;
                object dummyObject = objects[o];
                dummyObject.x0 = dummyObject.x0 - pbcVector;
                dummyObject.index = objectIndex;
                objects.push_back(dummyObject);
            }
        }
    }
}

// integration functions

void DEM::predictor() {
    static const double c1[5] = {deltat, deltat * deltat / 2.0, deltat * deltat * deltat / 6.0, deltat * deltat * deltat * deltat / 24.0, deltat * deltat * deltat * deltat * deltat / 120.0};
    static const double c2[5] = {deltat, deltat * deltat / 2.0, deltat * deltat * deltat / 6.0, deltat * deltat * deltat * deltat / 24.0, deltat * deltat * deltat * deltat * deltat / 120.0};
    //    c[0] = deltat;
    //    c[1] = c[0] * deltat / 2.0;
    //    c[2] = c[1] * deltat / 3.0;
    //    c[3] = c[2] * deltat / 4.0;
    //    c[4] = c[3] * deltat / 5.0;

#pragma omp parallel for
    for (int n = 0; n < elmts.size(); ++n) {
        elmts[n].predict(c1, c2);
    }
}

void DEM::corrector() {
    //static double gear[6] = {3.0/16.0, 251.0/360.0, 1.0, 11.0/18.0, 1.0/6.0, 1.0/60.0};
    //static const double c[5]={deltat, deltat*deltat/2.0, deltat*deltat*deltat/6.0, deltat*deltat*deltat*deltat/24.0, deltat*deltat*deltat*deltat*deltat/120.0};
    //static const double coeff[6]={gear[0]*c[1], gear[1]*c[1]/c[0], gear[2]*c[1]/c[1], gear[3]*c[1]/c[2], gear[4]*c[1]/c[3], gear[5]*c[1]/c[4]};


    // see: Computer Simulation of Liquids By M. P. Allen, D. J. Tildesle 
    static const double gear1ord[6] = {95.0 / 288.0, 1.0, 25.0 / 24.0, 35.0 / 72.0, 5.0 / 48.0, 1.0 / 120.0}; // coefficients for a first-order equation 
    static const double gear2ord[6] = {3.0 / 16.0, 251.0 / 360.0, 1.0, 11.0 / 18.0, 1.0 / 6.0, 1.0 / 60.0}; // coefficients for a second-order equation

    static const double c[5] = {deltat, deltat * deltat / 2.0, deltat * deltat * deltat / 6.0, deltat * deltat * deltat * deltat / 24.0, deltat * deltat * deltat * deltat * deltat / 120.0};
    // Careful the c's must be changed for eq of 1st order
    static const double coeff1ord[6] = {gear1ord[0] * c[0], gear1ord[1] * c[0] / c[0], gear1ord[2] * c[0] / c[1], gear1ord[3] * c[0] / c[2], gear1ord[4] * c[0] / c[3], gear1ord[5] * c[0] / c[4]};
    // coeff2ord has different c's sequence 
    static const double coeff2ord[6] = {gear2ord[0] * c[1], gear2ord[1] * c[1] / c[0], gear2ord[2] * c[1] / c[1], gear2ord[3] * c[1] / c[2], gear2ord[4] * c[1] / c[3], gear2ord[5] * c[1] / c[4]};


    //    doubleList coeff;
    //    coeff.resize(6);
    //    coeff[0]=gear[0]*c[1];
    //    coeff[1]=gear[1]*c[1]/c[0];
    //    coeff[2]=gear[2]*c[1]/c[1];
    //    coeff[3]=gear[3]*c[1]/c[2];
    //    coeff[4]=gear[4]*c[1]/c[3];
    //    coeff[5]=gear[5]*c[1]/c[4];

#pragma omp parallel for
    for (int n = 0; n < elmts.size(); n++) {
        elmts[n].correct(coeff1ord, coeff2ord);
    }
}

void DEM::evaluateForces(IO& io) {

    for (int n = 0; n < elmts.size(); ++n) {
        elmts[n].FParticle.reset();
        elmts[n].FWall.reset();
        elmts[n].MParticle.reset();
        elmts[n].MWall.reset();
        elmts[n].FSpringP.reset();
        elmts[n].FSpringW.reset();
        elmts[n].MRolling.reset();
    }

    for (int w = 0; w < walls.size(); ++w) {
        walls[w].FParticle.reset();
    }

    for (int o = 0; o < objects.size(); ++o) {
        objects[o].FParticle.reset();
    }

    if (io.saveCount != 0 && io.currentTimeStep % io.saveCount == 0) {
        // Open statistic file and print header
        printHeaderStatFile(io);
    }

    // forces due to particle overlap and lubrication
    particleParticleContacts(io);
    // forces due to contact with plane walls and lubrication
    wallParticleContacts(io);
    // forces due to contact with stationary spheres (objects)
    objectParticleContacts(io);
    // forces due to contact with cylinders
    cylinderParticelContacts();

    // close file here to be able to append the data from next saving
    io.statParticleFile.close();

    //  Newton equations solution
    for (int n = 0; n < elmts.size(); ++n) {

        // numerical viscosity for stability
        // see "Viscous torque on a sphere under arbitrary rotation" by Lei,  Yang, and Wu, Applied Physics Letters 89, 181908 (2006)
        const tVect FVisc = -6.0 * M_PI * numVisc * elmts[n].radius * elmts[n].xp1;
        const tVect MVisc = -8.0 * M_PI * numVisc * elmts[n].radius * elmts[n].radius * elmts[n].radius * elmts[n].wpGlobal;


        // translational motion
        // acceleration
        elmts[n].x2 = (FVisc + elmts[n].FHydro + elmts[n].FParticle + elmts[n].FWall) / elmts[n].m + demF;

        // rotational motion
        // adjoint of orientation quaternion
        //const tQuat q0adj=elmts[n].qp0.adjoint();
        // rotational velocity (body-fixed reference frame)
        //const tVect wBf=2.0*quat2vec( q0adj.multiply( elmts[n].qp1 ) );
        // moment in global reference frame
        const tVect moment = MVisc + elmts[n].MHydro + elmts[n].MParticle + elmts[n].MWall + elmts[n].MRolling;

        // moment in body-fixed reference frame
        //if (elmts[n].size
        const tVect momentBf = project(moment, elmts[n].qp0.adjoint());
        // rotational acceleration (body-fixed reference frame) (Newton equation for principal system)
        const tVect waBf = newtonAcc(momentBf, elmts[n].I, elmts[n].wpLocal);
        // rotational acceleration (vector)
        elmts[n].w1 = project(waBf, elmts[n].qp0);
        // rotational acceleration (quaternion)
        if (elmts[n].size > 1) {
            const tQuat waQuat = quatAcc(waBf, elmts[n].qp1);
            elmts[n].q2 = 0.5 * elmts[n].qp0.multiply(waQuat);
        }
    }
}

void DEM::updateParticlesPredicted() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

    //    #pragma omp parallel for
    for (int p = 0; p < stdParticles; ++p) {
        //getting belonging element index
        const unsigned int clusterIndex = particles[p].clusterIndex;
        particles[p].updatePredicted(elmts[clusterIndex], prototypes);
    }
    if (ghosts.size() != 0) {
        // updating ghost particles
        for (int g = 0; g < ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex = ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex = stdParticles + g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
}

void DEM::updateParticlesCorrected() {
    // this subroutine assumes that there is no change in the structure of the particles
    // any creation, deletion or reorganization of particles requires a different subroutine

    //    #pragma omp parallel for
    for (int p = 0; p < stdParticles; ++p) {
        //getting belonging element index
        const unsigned int clusterIndex = particles[p].clusterIndex;
        particles[p].updateCorrected(elmts[clusterIndex], prototypes);
    }

    if (ghosts.size() != 0) {
        // updating ghost particles
        for (int g = 0; g < ghosts.size(); ++g) {
            //getting original particle index
            const unsigned int originParticleIndex = ghosts[g].ghostIndex;
            const unsigned int ghostParticleIndex = stdParticles + g;
            // reconstructing ghost particle
            particles[ghostParticleIndex].ghostUpdate(particles[originParticleIndex], ghosts[g].pbcVector);
        }
    }
}

double DEM::criticalTimeStep() const {
    // determines the critical time step based on the stiffness and mass of the elements
    // we use YADE documentation, see https://yade-dem.org/doc/formulation.html

    double minRad = elmts[0].radius;
    double minMass = elmts[0].m;
    for (int n = 0; n < elmts.size(); ++n) {
        minRad = std::min(minRad, elmts[n].radius);
        minMass = std::min(minMass, elmts[n].m);
    }
    for (int n = 0; n < elmts.size(); ++n) {
        minRad = std::min(minRad, elmts[n].radius);
        minMass = std::min(minMass, elmts[n].m);
    }


    // double const k=8.0/15.0*sphereMat.youngMod/(1-sphereMat.poisson*sphereMat.poisson)*sqrt(minRad);

    double deltaTCrit = 0.0;
    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            // maximum length scale
            const double maxDist = std::max(demSize[0], std::max(demSize[1], demSize[2]));
            // modulus of acceleration
            const double maxAccel = demF.norm();
            // estimate maximum velocity (assuming governed by acceleration field))
            const double maxVel = std::sqrt(2.0 * maxAccel * maxDist);
            // see Landau & Lifshitz, or better Antypov & Elliott
            deltaTCrit = 2.214 * (2.0 * minRad) * std::pow(sphereMat.density / sphereMat.youngMod, 2.0 / 5.0) * std::pow(maxVel, 1.0 / 5.0);
            // deltaTCrit=minRad*sqrt(sphereMat.density/sphereMat.youngMod);
            break;
        }
        case LINEAR:
        {
            // sphereMat.dampCoeff is probably worng, check again
            deltaTCrit = M_PI / sqrt(sphereMat.linearStiff / minMass - sphereMat.dampCoeff * sphereMat.dampCoeff * sphereMat.linearStiff / minMass);
            break;
        }
    }
    return deltaTCrit;
}

// neighbor list functions

void DEM::initNeighborParameters() {
    // initializes all parameters useful for neighbor list algorithm

    cout << "Initialize neighbor list\n";
    // maximum radius
    double maxRad = 0.0;
    for (int i = 0; i < elmts.size(); ++i) {
        if (maxRad < elmts[i].radius) {
            maxRad = elmts[i].radius;
        }
    }

    // if there are no particles, just to avoid nonsense numbers, we use an only cell
    // therefore the cell width is the same as the simulation box (which comes from the LB settings)
    if (elmts.size() == 0) {
        cellWidth[0] = demSize[0];
        cellWidth[1] = demSize[1];
        cellWidth[2] = demSize[2];
    }        // otherwise the cell size is determined by the size of the particles
    else {
        cellWidth[0] = std::min(maxRad * 5.0, demSize[0]);
        cellWidth[1] = std::min(maxRad * 5.0, demSize[1]);
        cellWidth[2] = std::min(maxRad * 5.0, demSize[2]);
    }

    for (int k = 0; k < 3; ++k) {
        // ncells = numbers of cells for the linked cell algorithm
        nCells[k] = int( ceil(demSize[k] / cellWidth[k]));
        // width of the cells (actual)
        cellWidth[k] = demSize[k] / double(nCells[k]);
        // increase by two to give a container for ghost cells, in case of periodicity
        nCells[k] += 2;
    }

    // may need a revision
    nebrRange = std::max(maxRad * 3.0, 0.5 * std::min(cellWidth[0], std::min(cellWidth[1], cellWidth[2])));
    maxDisp = 0.5 * nebrRange;
    cout << "Neighbor list parameters\n";
    cout << "Number of Cells " << nCells[0] << " " << nCells[1] << " " << nCells[2] << "\n";
    cout << "Cell width " << cellWidth[0] << " " << cellWidth[1] << " " << cellWidth[2] << "\n";
    cout << "Range " << nebrRange << "\n";
}

void DEM::evalMaxDisp() {

    double maxVel = 0.0;
    for (int n = 0; n < elmts.size(); ++n) {
        if (elmts[n].x1.norm2() > maxVel) {
            maxVel = elmts[n].x1.norm2();
        }
    }
    maxVel = sqrt(maxVel);
    maxDisp += maxVel*deltat;
}

void DEM::evalCellTable() {
    // updates cellTable, a table which associates every particle to the belonging cell

    // delete precedent cell table
    cellTable.clear();
    // create a new cell table of size: number of particles + number of cells
    cellTable.resize(particles.size() + nCells[0] * nCells[1] * nCells[2]);

    // assigns value -1 to all elements of cellTable (the reason becomes clear when looking at the last two lines of this function)
    for (int n = 0; n < cellTable.size(); ++n) {
        // -1 is the default value for unassigned node
        cellTable[n] = -1;
    }

    // cycle through n = cellTable size (number of particles + number of cells)
    for (int n = 0; n < particles.size(); ++n) {

        // c is a identifier for the belonging cell of a particle. If particle A belongs to cell (2,5,7) in a cell
        // system with 4*6*9 cells, then it will be c = particleSize + 2 + (5*4 + 7)*6 = particleSize + 164 <- bullshit :D
        // basically is a system to store the cell coordinate of a particle in a line vector
        // Here the code: floor(x/cellWidth[0]) + nCells[0]*( floor(y/cellWidth[1]) + (nCells[1]*floor(z/cellWidth[2])) );
        const int c = particles[n].x0.linearizePosition(cellWidth, nCells) + particles.size();
        particles[n].tableCell = c;
        // if the coordinate exceeds the borders of the box, a message is displayed
        if (c > cellTable.size() || c < 0) { // just put control over ghost, here !!!!!!!!!!!!!!!!!!!!!!
            cout << "#neighborList, " << c << " initCellTable: particle " << n << " outside box, ignoring for force calculation." << endl;
            elmts[n].resetVelocity();
            elmts[n].radius = 0.0;
            elmts[n].x0 = Zero;
            //            exit(0);
            //            continue;
        }

        // cellTable is a structure to contain data in an efficient ways
        // every element of the array contains the pointer to another particle of the cell, until -1 is found
        // this explains why it needs to have particle.size() + cell.size() elements:
        // every particle points to another particle, and we need cell.size() elements with value -1
        cellTable[n] = cellTable[c];
        cellTable[c] = n;
    }

}

void DEM::evalNeighborTable() {
    //    cout<<"NEW NEIGHBOR TABLE"<<endl;
    //    cout<<"Neighbor table evaluation"<<endl;
    // prototype neighbors
    const unsigned int neighborCellNumber = 14;

    const static int shiftCell[neighborCellNumber][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
        {-1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1},
        {-1, 1, 1},
        {-1, 0, 1},
        {-1, -1, 1},
        {0, -1, 1},
        {1, -1, 1}
    };

    const double nebrRange2 = nebrRange*nebrRange;

    //updateParticles();
    // if periodicity is active: shift particles, identify ghosts and set the update signal for the LBM
    if (pbcs.size() != 0) {
        // cleaning up
        ghosts.clear();
        // resizing containers
        particles.resize(stdParticles);
        //resizing components of elements
        for (int n = 0; n < elmts.size(); ++n) {
            elmts[n].components.resize(elmts[n].size);
        }
        // periodicity shift
        pbcShift();
        // regenerate particles
        updateParticlesCorrected();
        // identify ghosts
        createGhosts();
        // update particle listing to get ready for neighborList
        updateParticlesCorrected();
        // update signal for LBM
        newNeighborList = true;

    }

    // every element in its cell
    evalCellTable();
    // delete precedent neighbor table
    neighborTable.clear();
    // elements size
    unsigned int particleSize = particles.size();

    // cycle through all cells ->i0
    // variables for cell coordinate
    unsigned int i0[3];
    for (i0[0] = 0; i0[0] < nCells[0]; ++i0[0]) {
        for (i0[1] = 0; i0[1] < nCells[1]; ++i0[1]) {
            for (i0[2] = 0; i0[2] < nCells[2]; ++i0[2]) {

                // cycle through neighbors vectors -> [s][k]
                for (int s = 0; s < neighborCellNumber; s++) {
                    // variables for cell coordinate
                    unsigned int i1[3];
                    for (int k = 0; k < 3; k++) {
                        // determines the neighbor cells starting from the considered cell and the
                        // prototype for neighbors (shifCell)
                        i1[k] = i0[k] + shiftCell[s][k];
                        // if neighbor cell is in upper border
                        if (i1[k] == nCells[k]) {
                            i1[k] = 0;
                        }
                        // if neighbor cell is in lower border
                        if (i1[k] < 0) {
                            i1[k] = nCells[k] - 1;
                        }
                    }

                    // linearized coordinate of cell0 (first cell of couples)
                    const unsigned int cell0 = (i0[2] * nCells[1] + i0[1]) * nCells[0] + i0[0] + particleSize;
                    // linearized coordinate of cell0 (first cell of couples)
                    const unsigned int cell1 = (i1[2] * nCells[1] + i1[1]) * nCells[0] + i1[0] + particleSize;

                    // this cycles through the elements of a cell, checking all neighbors
                    // the storage system is efficient but not trivial
                    int n0 = cellTable[cell0];
                    while (n0>-1) {
                        int n1 = cellTable[cell1];
                        while (n1>-1) {
                            // we save only the pair n1<n0 (otherwise twice the number of pairs), the pair n0<n1 (which is just the same) is excluded 
                            if (cell0 != cell1 || n1 < n0) {
                                // we exclude particles belonging to the same cluster
                                if (particles[n0].clusterIndex != particles[n1].clusterIndex) {
                                    // we also exclude too far away pairs
                                    const tVect r10 = particles[n1].x0 - particles[n0].x0;
                                    if (r10.norm2() < nebrRange2) {
                                        neighborTable.push_back(n0);
                                        neighborTable.push_back(n1);
                                    }
                                }
                            }
                            // to the next element (cell1)
                            n1 = cellTable[n1];
                        }
                        // to the next element (cell0)
                        n0 = cellTable[n0];
                    }
                }
            }
        }
    }

    if (walls.size() != 0) {
        evalNearWallTable();
    }
    if (cylinders.size() != 0) {
        evalNearCylinderTable();
    }
    if (objects.size() != 0) {
        evalNearObjectTable();
    }

}

void DEM::evalNearWallTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearWallTable.clear();
    unsIntList part;

    for (wallList::iterator ip = walls.begin(); ip != walls.end(); ++ip) {
        part.clear();

        for (int n = 0; n < stdParticles; n++) {
            if (ip->dist(particles[n].x0) < nebrRange) {
                part.push_back(n);
            }

        }

        nearWallTable.push_back(part); // to be reviewed

    }
}

void DEM::evalNearObjectTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearObjectTable.clear();
    unsIntList part;

    const double nebrRange2 = nebrRange*nebrRange;

    for (int o = 0; o < objects.size(); ++o) {
        part.clear();
        for (int n = 0; n < stdParticles; ++n) {
            tVect x0ij = particles[n].x0 - objects[o].x0;

            if (x0ij.norm2() < nebrRange2) {
                part.push_back(n);
            }
        }
        nearObjectTable.push_back(part);
    }
}

void DEM::evalNearCylinderTable() {
    // evaluate the distance of all particles to the cylinders. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearCylinderTable.clear();
    unsIntList part;

    for (cylinderList::iterator ip = cylinders.begin(); ip != cylinders.end(); ++ip) {
        part.clear();
        for (int n = 0; n < stdParticles; ++n) {
            if (ip->dist(particles[n].x0) < nebrRange) {
                part.push_back(n);

            }
        }
        nearCylinderTable.push_back(part);
    }
}

/*void DEM::evalNearWallTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearWallTable.clear();

    for (int n=0; n<stdParticles; ++n) {
        for (int w=0; w<walls.size(); ++w) {
            if (walls[w].dist(particles[n].x0)<nebrRange) {
        // for (wallList::iterator ip=walls.begin(); ip!=walls.end(); ++ip) {
            //if (ip->dist(particles[n].x0)<nebrRange) {
                nearWallTable.push_back(n);
                nearWallTable.push_back(w);
                break;
            }
        }
    }
}*/

/*void DEM::evalNearObjectTable() {
    // evaluate the distance of all particles to the walls. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearObjectTable.clear();

    const double nebrRange2=nebrRange*nebrRange;

    for (int n=0; n<stdParticles; ++n) {
        for (int o=0; o<objects.size(); ++o) {
            tVect x0ij = particles[n].x0 - objects[o].x0;
            if (x0ij.norm2()<nebrRange2) {
                nearObjectTable.push_back(n);
                nearObjectTable.push_back(o);
                //break;
            }
        }
    }
}*/

/*void DEM::evalNearCylinderTable() {
    // evaluate the distance of all particles to the cylinders. Creates a list of all particles
    // inside a certain range (nebrRange)

    nearCylinderTable.clear();

    for (int n=0; n<stdParticles; ++n) {
        for (cylinderList::iterator ip=cylinders.begin(); ip!=cylinders.end(); ++ip) {
            if (ip->dist(particles[n].x0) < nebrRange) {
                nearCylinderTable.push_back(n);
                break;
            }
        }
    }
}*/

// periodicity functions

void DEM::pbcShift() {
    //shifts elements and particle according to boundary conditions

    // this subroutine works on elements: particles need be generated/updated elsewhere
    for (int b = 0; b < pbcs.size(); ++b) {
        const tVect pbcVector = pbcs[b].v;
        // cycle through elements
        for (int n = 0; n < elmts.size(); ++n) {
            const double leftDist = pbcs[b].pl1.dist(elmts[n].x0);
            const double rightDist = pbcs[b].pl2.dist(elmts[n].x0);
            // first plane of couple (left)
            if (leftDist < 0.0) {
                // cout<<"Shift element number"<<elmts[n].index<<"\n";
                elmts[n].translate(pbcVector);
                // shift particles (not the ghosts)
                //                for (int j=0; j<elmts[n].size; ++j) {
                //                    particles[elmts[n].components[j]].x0+=pbcVector;
                //                }
            }
            // second plane of couple (right)
            if (rightDist < 0.0) {
                //                cout<<"Shift element number"<<elmts[n].index<<"\n";
                elmts[n].translate(-1.0 * pbcVector);
                // shift particles (not the ghosts)
                //                for (int j=0; j<elmts[n].size; ++j) {
                //                    particles[elmts[n].components[j]].x0-=pbcVector;
                //                }
            }
        }
    }
}

void DEM::createGhosts() {

    ghosts.clear();

    // this subroutine implies that elements (and their non-ghost particles) have already been shifted
    for (int b = 0; b < pbcs.size(); ++b) {
        const tVect pbcVector = pbcs[b].v;
        // cycle through elements
        for (int n = 0; n < elmts.size(); ++n) {
            // cycle through standard particles
            for (int j = 0; j < elmts[n].size; ++j) {
                const unsigned int p = elmts[n].components[j];
                // distances from the periodic walls
                const double leftDist = pbcs[b].pl1.dist(particles[p].x0);
                const double rightDist = pbcs[b].pl2.dist(particles[p].x0);
                //            ASSERT(leftDist>0);
                //            ASSERT(rightDist>0)
                // first plane of couple (left)
                if (leftDist < nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = particles[p].particleIndex;
                    dummyGhost.pbcVector = pbcVector;
                    ghosts.push_back(dummyGhost);
                }                    // second plane of couple (right), we copy only one time
                else if (rightDist < nebrRange) {
                    ghost dummyGhost;
                    dummyGhost.ghostIndex = particles[p].particleIndex;
                    dummyGhost.pbcVector = -1.0 * pbcVector;
                    ghosts.push_back(dummyGhost);
                }
            }
        }
        //        for (int p=0; p<stdParticles; ++p) {
        //            // distances from the periodic walls
        //            const double leftDist=pbcs[b].pl1.dist(particles[p].x0);
        //            const double rightDist=pbcs[b].pl2.dist(particles[p].x0);
        ////            ASSERT(leftDist>0);
        ////            ASSERT(rightDist>0)
        //            // first plane of couple (left)
        //            if (leftDist<nebrRange) {
        //                ghost dummyGhost;
        //                dummyGhost.ghostIndex=particles[p].particleIndex;
        //                dummyGhost.pbcVector=pbcVector;
        //                ghosts.push_back(dummyGhost);
        //            }
        //            // second plane of couple (right), we copy only one time
        //            else if (rightDist<nebrRange) {
        //                ghost dummyGhost;
        //                dummyGhost.ghostIndex=particles[p].particleIndex;
        //                dummyGhost.pbcVector=-1.0*pbcVector;
        //                ghosts.push_back(dummyGhost);
        //            }
        //        }
    }
    // now we need to check for particles at corners
    // these are particles have already been shifted by two periodic walls IMPORTANT: (the case of three periodic walls is missing)
    ghostList cornerGhosts;
    cornerGhosts.clear();
    for (int g1 = 0; g1 < ghosts.size(); ++g1) {
        for (int g2 = g1 + 1; g2 < ghosts.size(); ++g2) {
            if (ghosts[g1].ghostIndex == ghosts[g2].ghostIndex) {
                // we need to create a corner ghost
                ghost dummyGhost;
                dummyGhost.ghostIndex = ghosts[g1].ghostIndex;
                dummyGhost.pbcVector = ghosts[g1].pbcVector + ghosts[g2].pbcVector;
                cornerGhosts.push_back(dummyGhost);
            }
        }
    }
    //    cout<<cornerGhosts.size()<<endl;
    ghosts.insert(ghosts.end(), cornerGhosts.begin(), cornerGhosts.end());

    // resizing containers
    particles.resize(stdParticles);
    //resizing components of elements
    for (int t = 0; t < elmts.size(); ++t) {
        elmts[t].components.resize(elmts[t].size);
    }

    // creating particles from ghosts
    // index for ghost particles
    unsigned int ind = stdParticles;
    // adding ghost particles to particle list
    for (int g = 0; g < ghosts.size(); ++g) {
        //getting original particle index
        const unsigned int originParticleIndex = ghosts[g].ghostIndex;
        // reconstructing ghost particle
        particle dummyParticle = particles[originParticleIndex];
        //        dummyParticle.x0+=ghosts[g].pbcVector;
        dummyParticle.particleIndex = ind;
        // adding ghost particle to components of the mother element
        elmts[dummyParticle.clusterIndex].components.push_back(dummyParticle.particleIndex);
        // updating particle list
        particles.push_back(dummyParticle);
        // updating particle index
        ind++;
    }
}

// force computation functions

void DEM::particleParticleContacts(IO& io) {

    for (unsIntList::iterator ipi = neighborTable.begin(); ipi != neighborTable.end(); ipi = ipi + 2) {
        // couple of contact candidates
        unsIntList::iterator ipj = ipi + 1;
        // pointers to particles


        particle *parti = &particles[*ipi];
        particle *partj = &particles[*ipj];


        string s;
        if (parti->particleIndex > partj->particleIndex) {
            s = "p" + to_string(parti->particleIndex) + ":p" + to_string(partj->particleIndex);
        } else {
            s = "p" + to_string(parti->particleIndex) + ":p" + to_string(partj->particleIndex);
            partj = &particles[*ipi];
            parti = &particles[*ipj];
        }



        // checking for overlap
        const double ri = parti->r;
        const double rj = partj->r;
        const double sigij = ri + rj;
        const double sigij2 = sigij*sigij;
        // distance between centers
        const tVect vectorDistance = partj->x0 - parti->x0;
        const double distance2 = vectorDistance.norm2();

        tVect spring = parti->x0 + vectorDistance / 2;

        Elongation elongation_here = elongTable[s];
        // check for lubrication contact
        // check for contact
        if (distance2 < sigij2) {
            //particleParticleCollision(parti,partj,vectorDistance,io);

            //cout<<parti->particleIndex<<" "<<partj->particleIndex<<" "<<s<<" ";


            particleParticleCollision(parti, partj, vectorDistance, io, elongation_here);
            elongation_here.p = spring;
            elongTable[s] = elongation_here;


            //elongation_here.e.show();
            //cout<<elongation_here.sliping<<endl;

        } else {
            elongTable.erase(s);
        }
    }
}

void DEM::wallParticleContacts(IO& io) {
    // to keep conventions, the index i refers to the wall, and j the particle

    for (wallList::iterator ip = walls.begin(); ip != walls.end(); ip++) {
        wall *wallI = &*ip;
        //cout<<"muro n "<<wallI->index<<endl;
        // cycling through particle in wall neighbor list
        for (unsIntList::iterator it = nearWallTable[wallI->index].begin(); it != nearWallTable[wallI->index].end(); it++) {

            // particle
            const particle *partJ = &particles[*it];

            // radius
            const double rj = partJ->r;
            // distance from wall (norm)
            const double distance = wallI->dist(partJ->x0);
            bool wallIsOk = true;

            tVect spring = partJ->x0 - wallI->p;

            const string s = "w" + to_string(wallI->index) + ":p" + to_string(partJ->particleIndex);

            Elongation elongation_here = elongTable[s];
            if (wallIsOk) {
                // distance before contact
                const double overlap = rj - distance;
                if (overlap > 0.0) {
                    wallParticleCollision(wallI, partJ, overlap, io, elongation_here);
                    elongation_here.p = spring;
                    elongTable[s] = elongation_here;
                } else {
                    elongTable.erase(s);
                }
            }
        }
    }
}

void DEM::cylinderParticelContacts() {

    // to keep conventions, the index i refers to the cylinder, and j the particle
    // cycling through the cylinders
    for (cylinderList::iterator ip = cylinders.begin(); ip != cylinders.end(); ip++) {
        cylinder *cylinderI = &*ip;
        // cycling through the cylinder neighbor particles
        for (unsIntList::iterator it = nearCylinderTable[cylinderI->index].begin(); it != nearCylinderTable[cylinderI->index].end(); it++) {

            // particle
            const particle *partJ = &particles[*it];

            // radius
            const double rj = partJ->r;
            // distance from wall (norm)
            const double distance = cylinderI->dist(partJ->x0);
            // distance before contact
            const double overlap = rj - distance;

            // distance to point 1 of axis
            const tVect p1dist = partJ->x0 - cylinderI->p1;
            // same but projected on the axis
            const tVect p1distax = (p1dist.dot(cylinderI->naxes)) * cylinderI->naxes;
            // distance of point from cylinder axis
            const tVect p1distcylinder = p1distax - p1dist;
            tVect spring = cylinderI->p1 - p1distcylinder;
            const string s = "c" + to_string(cylinderI->index) + ":p" + to_string(partJ->particleIndex);

            Elongation elongation_here = elongTable[s];

            // check for contact
            if (overlap > 0.0) {
                cylinderParticleCollision(cylinderI, partJ, overlap, elongation_here);
                elongation_here.p = spring;
                elongTable[s] = elongation_here;
            } else {
                elongTable.erase(s);
            }
        }
    }
}

void DEM::objectParticleContacts(IO& io) {

    // to keep conventions, the index i refers to the object, and j the particle
    for (objectList::iterator iob = objects.begin(); iob != objects.end(); iob++) {
        object *objectI = &*iob;

        for (unsIntList::iterator it = nearObjectTable[objectI->index].begin(); it != nearObjectTable[objectI->index].end(); it++) {
            // particle
            const particle *partJ = &particles[*it];

            // radius
            const double rj = partJ->r;
            // distance from object (vector)
            const tVect vectorDistance = partJ->x0 - objectI->x0;
            // distance from object (norm)
            const double distance = vectorDistance.norm();
            // distance before contact
            const double overlap = rj + objectI->r - distance;

            const string s = "o" + to_string(objectI->index) + ":p" + to_string(partJ->particleIndex);
            tVect spring = objectI->x0 - vectorDistance / 2;

            Elongation elongation_here = elongTable[s];
            elongation_here.p = spring;

            if (overlap > 0.0) {
                objectParticleCollision(objectI, partJ, vectorDistance, io, elongation_here);
                elongation_here.p = spring;
                elongTable[s] = elongation_here;
            } else {
                elongTable.erase(s);
            }
        }
    }
}

//void DEM::objectParticleContacts(IO& io) {
//
//    // to keep conventions, the index i refers to the object, and j the particle
//
//    for (objectList::iterator ib=objects.begin(); ib!=objects.end(); ib++) {
//        object *objectI=&*ib;
//        for (unsIntList::iterator it=nearObjectTable.begin(); it!=nearObjectTable.end(); it=it++) {
//            // particle
//            const particle *partJ=&particles[*it];
//            // radius
//            const double rj=partJ->r;
//            // distance from object (vector)
//            const tVect vectorDistance=partJ->x0-objectI->x0;
//            // distance from object (norm)
//            const double distance=vectorDistance.norm();
//            // distance before contact
//            double overlap=rj+objectI->r-distance;
//            if (overlap>-cutOff){
//                if (overlap>0.0){
//                    objectParticleCollision(objectI,partJ,vectorDistance,io);
//                }
//                else {
//                    // lubrication
//                }
//            }
//        }
//    }
//}

inline void DEM::particleParticleCollision(const particle *partI, const particle *partJ, const tVect& vectorDistance, IO& io, Elongation& elongation) {

    // pointers to elements
    elmt *elmtI = &elmts[partI->clusterIndex];
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry /////////////////////////////
    const double radI = partI->r;
    const double radJ = partJ->r;
    // distance norm
    const double distance = vectorDistance.norm();
    // overlap
    const double overlap = radI + radJ - distance;
    // relative velocity
    const tVect relVel = partJ->x1 - partI->x1;
    // first local unit vector (normal)
    const tVect en = vectorDistance / distance;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;
    // effective mass
    const double effMass = elmtI->m * elmtJ->m / (elmtI->m + elmtJ->m);
    // effective radius
    const double effRad = radI * radJ / (radI + radJ);

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(overlap, normRelVel, effRad, effMass);
    ASSERT(normNormalForce >= 0.0);

    //    // Overlap elastic potential energy
    //    energy.elastic+=0.5*fn*xi*xi;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radii
    const tVect vecRadI = radI*en;
    const tVect vecRadj = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistI = vecRadI + (partI->radiusVec);
    const tVect centerDistJ = vecRadj + (partJ->radiusVec);

    if (partI->particleIndex < stdParticles) {
        elmtI->FParticle = elmtI->FParticle - normalForce;
        //  moment generated in non-spherical particles
        if (elmtI->size > 1) {
            elmtI->MParticle = elmtI->MParticle - centerDistI.cross(normalForce);
        }
    }
    if (partJ->particleIndex < stdParticles) {
        elmtJ->FParticle = elmtJ->FParticle + normalForce;
        //  moment generated in non-spherical particles
        if (elmtJ->size > 1) {
            elmtJ->MParticle = elmtJ->MParticle + centerDistJ.cross(normalForce);
        }
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wI = elmtI->wpGlobal; //2.0*quat2vec( elmtI->qp1.multiply( elmtI->qp0.adjoint() ) );
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel - wI.cross(vecRadI) + wJ.cross(vecRadj);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {

        const double elong_old_norm = elongation.e.norm();


        double normElong = elongation.e.dot(en);
        const tVect normalElong = en*normElong;


        elongation.e = elongation.e - normalElong;

        if (elongation.e.norm() != 0) {
            const double scaling = elong_old_norm / elongation.e.norm();
            elongation.e = elongation.e*scaling;
        } else {
            const double scaling = 0.0;
            elongation.e = elongation.e*scaling;
        }

        //elongation.e.show();
        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, effMass, elongation, sphereMat.frictionCoefPart,
                sphereMat.linearStiff, sphereMat.viscTang);



        // torque updating
        if (partI->particleIndex < stdParticles) {
            elmtI->MParticle = elmtI->MParticle + centerDistI.cross(tangForce);
            elmtI->FSpringP = elmtI->FSpringP + tangForce;
            elmtI->FParticle = elmtI->FParticle + tangForce;

        }
        if (partJ->particleIndex < stdParticles) {
            elmtJ->MParticle = elmtJ->MParticle - centerDistJ.cross(tangForce);
            elmtJ->FSpringP = elmtJ->FSpringP - tangForce;
            elmtJ->FParticle = elmtJ->FParticle - tangForce;
        }
        //tangForce.show();
        // save particle-particle collision into statistics file
        //saveStatPPcollision(io,partI,partJ,overlap,normNormalForce,en,normTangForce,et);

    } else {
        //cout<<"normTangRelVelContact==0!!!"<<endl;
    }
    //tangRelVelContact.show();
    //cout<<normNormalForce;


    //ROLLING
    const double rolling_coeff = 0.5;

    tVect rollingMoment = rollingContact(wI, wJ, effRad, normNormalForce,
            rolling_coeff);


    if (partI->particleIndex < stdParticles) {
        elmtI->MRolling = elmtI->MRolling - rollingMoment;

    }
    if (partJ->particleIndex < stdParticles) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }

}

inline void DEM::wallParticleCollision(wall *wallI, const particle *partJ, const double& overlap, IO& io, Elongation& elongation) {

    // pointers to element
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = partJ->r;
    // first local unit vector (normal)
    const tVect en = wallI->n;
    // speed of the wall at contact point
    const tVect contactPointVelocity = wallI->getSpeed(partJ->x0); // fix this, contact point not defined
    // relative velocity
    const tVect relVel = partJ->x1 - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(2.0 * overlap, normRelVel, radJ, elmtJ->m);


    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    tVect centerDistJ = vecRadJ;
    if (elmtJ->size > 1) {
        // vectorized distance contactPoint-center of cluster
        centerDistJ = centerDistJ + (partJ->x0 - elmtJ->xp0);
    }

    // force updating
    elmtJ->FWall = elmtJ->FWall + normalForce;
    wallI->FParticle = wallI->FParticle - normalForce;
    // torque updating
    if (elmtJ->size > 1) {
        elmtJ->MWall = elmtJ->MWall + centerDistJ.cross(normalForce);
    }


    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();


    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {


        const double elong_old_norm = elongation.e.norm();

        double normElong = elongation.e.dot(en);
        const tVect normalElong = en*normElong;

        elongation.e = elongation.e - normalElong;

        if (elongation.e.norm() != 0) {
            const double scaling = elong_old_norm / elongation.e.norm();
            elongation.e = elongation.e*scaling;
        } else {
            const double scaling = 0.0;

            elongation.e = elongation.e*scaling;
        }

        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, elmtJ->m, elongation, sphereMat.frictionCoefWall,
                sphereMat.linearStiff, sphereMat.viscTang);

        // torque updating
        elmtJ->MWall = elmtJ->MWall - centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FSpringW = elmtJ->FSpringW - tangForce;
        elmtJ->FWall = elmtJ->FWall - tangForce;

        wallI->FParticle = wallI->FParticle + tangForce;
        // save wall-particle collision into statistics file
        //saveStatWPcollision(io, wallI, partJ, overlap, normNormalForce, en,normTangForce, et);

    }
    //ROLLING
    const double rolling_coeff = 2 * 0.5;

    tVect rollingMoment = rollingContact(contactPointVelocity, wJ, radJ, normNormalForce,
            rolling_coeff);



    if (partJ->particleIndex < stdParticles) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }

}

inline void DEM::cylinderParticleCollision(cylinder *cylinderI, const particle *partJ, const double& overlap, Elongation& elongation) {

    // pointers to element
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = partJ->r;
    // vectorial distance
    const tVect vecDistance = cylinderI->vecDist(partJ->x0);
    // contact point
    const tVect contactPoint = partJ->x0 - vecDistance;
    // first local unit vector (normal)
    const tVect en = vecDistance / (radJ - overlap);
    // speed of the cylinder at contact point
    const tVect contactPointVelocity = cylinderI->getSpeed(contactPoint);
    // relative velocity
    const tVect relVel = partJ->x1 - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(2.0 * overlap, normRelVel, radJ, elmtJ->m);

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ = vecRadJ + (partJ->x0 - elmtJ->xp0);

    // force updating
    elmtJ->FWall = elmtJ->FWall + normalForce;
    // wallI->FParticle=wallI->FParticle-fnv;
    // torque updating
    if (elmtJ->size > 1) {
        elmtJ->MWall = elmtJ->MWall + centerDistJ.cross(normalForce);
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ); // couldn't we just use elmtJ.w?
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {
        const double elong_old_norm = elongation.e.norm();

        double normElong = elongation.e.dot(en);
        const tVect normalElong = en*normElong;

        elongation.e = elongation.e - normalElong;

        if (elongation.e.norm() != 0) {
            const double scaling = elong_old_norm / elongation.e.norm();
            elongation.e = elongation.e*scaling;
        } else {
            const double scaling = 0.0;

            elongation.e = elongation.e*scaling;
        }



        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, elmtJ->m, elongation, sphereMat.frictionCoefWall,
                sphereMat.linearStiff, sphereMat.viscTang);

        /*// tangential force
        double normTangForce = tangentialContact(normTangRelVelContact, normNormalForce, radJ, elmtJ->m, sphereMat.frictionCoefWall);
        // second local unit vector (tangential)
        const tVect et = tangRelVelContact / tangRelVelContact.norm();
        // vectorial tangential force
        const tVect tangForce = normTangForce*et;*/

        // torque updating
        elmtJ->MWall = elmtJ->MWall - centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FWall = elmtJ->FWall - tangForce;
        //wallI->FParticle = wallI->FParticle+ftv;
    }
    //ROLLING
    const double rolling_coeff = 2 * 0.5;

    tVect rollingMoment = rollingContact(contactPointVelocity, wJ, radJ, normNormalForce,
            rolling_coeff);


    if (partJ->particleIndex < stdParticles) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }
}

inline void DEM::objectParticleCollision(object *objectI, const particle *partJ, const tVect& vectorDistance, IO& io, Elongation& elongation) {

    // pointers to element
    elmt *elmtJ = &elmts[partJ->clusterIndex];

    // NORMAL FORCE ///////////////////////////////////////////////////////////////////

    // geometry ///////////////
    // particle radius
    const double radJ = partJ->r;
    // distance from object (norm)
    const double distance = vectorDistance.norm();
    // distance before contact
    double overlap = radJ + objectI->r - distance;
    // first local unit vector (normal)
    const tVect en = (1.0 / distance) * vectorDistance;
    // speed of the wall at contact point
    const tVect contactPointVelocity = objectI->x1; // fix this, contact point not defined
    // relative velocity
    const tVect relVel = partJ->x1 - contactPointVelocity;
    // relative normal velocity (modulus)
    const double normRelVel = relVel.dot(en);
    // relative normal velocity
    const tVect normalRelVel = en*normRelVel;

    // force computation /////////////////////////////////
    const double normNormalForce = normalContact(2.0 * overlap, normRelVel, radJ, elmtJ->m);

    // Overlap elastic potential energy
    //                    energy.elastic+=fn*xj/2.5;

    // force application ///////////////////////////////////
    // vectorial normal force
    const tVect normalForce = en*normNormalForce;

    // cluster geometry
    // vectorized radius
    const tVect vecRadJ = -radJ*en;
    // vectorized distance contactPoint-center of cluster
    const tVect centerDistJ = vecRadJ + (partJ->x0 - elmtJ->xp0);

    // force updating
    elmtJ->FWall = elmtJ->FWall + normalForce;
    objectI->FParticle = objectI->FParticle - normalForce;
    // torque updating
    if (elmtJ->size > 1) {
        elmtJ->MWall = elmtJ->MWall + centerDistJ.cross(normalForce);
    }

    // TANGENTIAL FORCE ///////////////////////////////////////////////////////////////////

    // angular velocities (global reference frame)
    const tVect wJ = elmtJ->wpGlobal; //2.0*quat2vec( elmtJ->qp1.multiply( elmtJ->qp0.adjoint() ) );
    // relative velocity at contact point
    const tVect relVelContact = relVel + wJ.cross(vecRadJ);
    // tangential component of relative velocity
    const tVect tangRelVelContact = relVelContact - normalRelVel;
    // norm of tangential velocity
    const double normTangRelVelContact = tangRelVelContact.norm();
    // checking if there is any tangential motion
    if (normTangRelVelContact != 0.0) {
        const double elong_old_norm = elongation.e.norm();

        double normElong = elongation.e.dot(en);
        const tVect normalElong = en*normElong;

        elongation.e = elongation.e - normalElong;

        if (elongation.e.norm() != 0) {
            const double scaling = elong_old_norm / elongation.e.norm();
            elongation.e = elongation.e*scaling;
        } else {
            const double scaling = 0.0;

            elongation.e = elongation.e*scaling;
        }


        tVect tangForce = FRtangentialContact(tangRelVelContact, normNormalForce, elmtJ->m, elongation, sphereMat.frictionCoefWall,
                sphereMat.linearStiff, sphereMat.viscTang);



        // torque updating
        elmtJ->MWall = elmtJ->MWall - centerDistJ.cross(tangForce);
        // force updating
        elmtJ->FWall = elmtJ->FWall - tangForce;

        elmtJ->FSpringW = elmtJ->FSpringW - tangForce;
        objectI->FParticle = objectI->FParticle + tangForce;
        // save object-particle collision into statistics file
        //saveStatOPcollision(io, objectI, partJ, overlap, normNormalForce, en, normTangForce, et);
    }
    //ROLLING
    const double rolling_coeff = 2 * 0.5;

    tVect rollingMoment = rollingContact(contactPointVelocity, wJ, radJ, normNormalForce,
            rolling_coeff);


    if (partJ->particleIndex < stdParticles) {
        elmtJ->MRolling = elmtJ->MRolling + rollingMoment;
    }
}

double DEM::normalContact(const double& overlap, const double& vrelnnorm, const double& effRad, const double& effMass) const {

    // total normal force
    double fn = 0.0;

    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            // square root of effective radius
            const double sqrtEffRad = sqrt(effRad);
            // square root of overlap
            const double sqrtOverlap = sqrt(overlap);
            // normal stiffness = 2/3*(Reff^1/2)*Y/(1-nu²)
            const double kn = sphereMat.knConst * sqrtEffRad*sqrtOverlap;
            // damping
            const double gamman = 2.0 * sphereMat.dampCoeff * sqrt(kn * effMass); // 2 was missing
            const double dumpfn = -gamman*vrelnnorm;
            // elastic normal force (Hertzian contact)
            const double elasticfn = kn*overlap;
            // total normal force (limit to repulsive force)
            fn = std::max(elasticfn + dumpfn, 0.0);
            break;
        }
        case LINEAR:
        {
            // elastic normal force (linear contact)
            const double elasticfn = sphereMat.linearStiff*overlap;
            // damping
            const double gamman = 2.0 * sphereMat.dampCoeff * sqrt(sphereMat.linearStiff * effMass);
            const double dumpfn = -gamman*vrelnnorm;
            // total normal force (limit to repulsive force)
            fn = std::max(elasticfn + dumpfn, 0.0);
            break;
        }
    }

    return fn;

}

tVect DEM::FRtangentialContact(const tVect& tangRelVelContact, const double& fn,
        const double& effMass, Elongation& elongation, const double& friction, const double& tangStiff, const double& viscTang) {
    // tangential force
    tVect fs = tVect(0.0, 0.0, 0.0);
    // tangent stiffness -> not physically sound, static friction is actually missing
    const double ks = tangStiff;
    // maximum force due to static friction
    const double fsMax = friction*fn;
    const double fdMax = friction*fn;

    // viscous force
    const tVect viscousForce = 2.0 * viscTang * sqrt(effMass * ks) * tangRelVelContact;



    /*tVect et=tVect (0.0,0.0,0.0);
   if (tangRelVelContact.norm()!=0){
   et = tangRelVelContact / tangRelVelContact.norm();
   } */


    const tVect F_control = elongation.e * ks + viscousForce;

    const double F_control_norm = F_control.norm();

    tVect et = tVect(0.0, 0.0, 0.0);
    if (F_control_norm != 0) {
        et = F_control / F_control_norm;
    }

    if (!elongation.sliping) {
        if (F_control_norm < fsMax) {
            fs = F_control;
            elongation.e = elongation.e + tangRelVelContact*deltat;

        } else {
            fs = fdMax*et;
            elongation.sliping = true;
            const tVect f = fdMax * et - viscousForce;
            elongation.e = f / ks;

        }
    }
    else {

        if (F_control_norm > fdMax) {

            fs = fdMax*et;
            const tVect f = fdMax * et - viscousForce;
            elongation.e = f / ks;

        } else {
            fs = F_control;
            elongation.e = elongation.e + tangRelVelContact*deltat;
            elongation.sliping = false;

        }

    }



    return fs;
}

double DEM::tangentialContact(const double& vreltNorm, const double& fn, const double& effRad, const double& effMass, const double& friction) const {

    // tangential force
    double fs = 0.0;

    switch (sphereMat.contactModel) {
        case HERTZIAN:
        {
            static const double power = 1.0 / 3.0;
            // square root of effective radius
            const double sqrtEffRad = sqrt(effRad);
            // tangent stiffness -> not physically sound, static friction is actually missing
            const double ks = sphereMat.ksConst * sqrtEffRad * std::pow(std::abs(fn), power);
            // maximum force due to dynamic friction
            const double fsMax = friction*fn;
            // viscous force
            const double gammas = 2.0 * sphereMat.viscTang * sqrt(effMass * ks);
            const double dumpfs = gammas*vreltNorm;
            // old version
            // const double viscousForce=2.0*sphereMat.viscTang*sqrt(effMass*ks)*vreltNorm;
            // capping with dynamic friction
            fs = std::min(dumpfs, fsMax);
            //ASSERT(fsMax>=0.0);
            break;
        }
        case LINEAR:
        {
            // tangent stiffness -> not physically sound, static friction is actually missing
            const double ks = sphereMat.linearStiff;
            // maximum force due to dynamic friction
            const double fsMax = friction*fn;
            // viscous force
            const double gammas = 2.0 * sphereMat.viscTang * sqrt(effMass * ks);
            const double dumpfs = gammas*vreltNorm;
            // old version
            //const double viscousForce=2.0*sphereMat.viscTang*sqrt(effMass*ks)*vreltNorm;
            // capping with dynamic friction
            fs = std::min(dumpfs, fsMax);
            break;
        }
    }
    //if(fs<0.0){
    //    cout<<"fs="<<fs<<endl;
    //}
    ASSERT(!isnan(fs));

    return fs;

}

void DEM::deleteParticleElongation() {

    for (unsIntList::iterator iti = neighborTable.begin(); iti != neighborTable.end(); iti = iti + 2) {
        // couple of contact candidates
        unsIntList::iterator itj = iti + 1;
        // pointers to particles
        const particle *parti = &particles[*iti];
        const particle *partj = &particles[*itj];

        string s;
        if (parti->particleIndex > partj->particleIndex) {
            s = "p" + to_string(parti->particleIndex) + ":p" + to_string(partj->particleIndex);
        } else {
            s = "p" + to_string(partj->particleIndex) + ":p" + to_string(parti->particleIndex);
        }


        // checking for overlap
        const double ri = parti->r;
        const double rj = partj->r;
        const double sigij = ri + rj;
        const double sigij2 = sigij*sigij;
        // distance between centers
        const tVect vectorDistance = partj->x0 - parti->x0;
        const double distance2 = vectorDistance.norm2();
        tVect spring = parti->x0 + vectorDistance / 2;

        // check for contact
        if (distance2 > sigij2) {

            //cout<<"coppie cancellate"<<s<<endl;

            //unsigned int c=0;

            elongTable.erase(s); //before it was erased, but so the couple is eliminated in the map    
            // apply lubrication force
            //                lubrication(i, j, ri, rj, vectorDistance);
            //cout<<"numero elementi  "<<c<<endl;
        }
    }
}

void DEM::deleteWallElongation() {
    for (wallList::iterator ip = walls.begin(); ip != walls.end(); ip++) {
        wall *wallI = &*ip;
        //cout<<"muro n "<<wallI->index<<endl;
        // cycling through particle in wall neighbor list
        for (unsIntList::iterator it = nearWallTable[wallI->index].begin(); it != nearWallTable[wallI->index].end(); it++) {

            // particle
            const particle *partJ = &particles[*it];

            //cout<<"particella n "<<partJ->particleIndex<<endl;
            // radius
            const double rj = partJ->r;
            // distance from wall (norm)
            const double distance = wallI->dist(partJ->x0);


            const string s = "w" + to_string(wallI->index) + ":p" + to_string(partJ->particleIndex);


            // distance before contact
            const double overlap = rj - distance;
            if (overlap <= 0.0) {
                elongTable.erase(s);

            }
        }
    }
}

tVect DEM::rollingContact(const tVect& wI, const tVect& wJ, const double effRad, const double& fn,
        const double& rolling) {
    tVect fr = tVect(0.0, 0.0, 0.0);
    const tVect wrel = wI - wJ;
    const double wrel_norm = wrel.norm();

    tVect wIJ = tVect(0.0, 0.0, 0.0);
    if (wrel_norm != 0) {
        wIJ = wrel / wrel_norm;
    }


    fr = wIJ * rolling * fn*effRad;

    return fr;
}


// energy functions

void DEM::updateEnergy() {

    //resetting energy
    particleEnergy.reset();

    // kinetic energy
    double tKin(0.0), rKin(0.0);
    for (int n = 0; n < elmts.size(); ++n) {
        tKin += 0.5 * elmts[n].m * elmts[n].x1.norm2();
        // adjoint of orientation quaternion
        //const tQuat q0adj=elmts[n].q0.adjoint();
        // rotational velocity (body-fixed reference frame)
        //const tVect w=wLocal;//2.0*quat2vec( q0adj.multiply( elmts[n].q1 ) );
        //        Tvect w=2.0*quat2vec( elmts[i].q1.multiply( elmts[i].q0.adjoint() ) );
        const tVect wSquare = elmts[n].wLocal.compProd(elmts[n].wLocal);
        rKin += elmts[n].I.dot(wSquare);
    }

    particleEnergy.rotKin = rKin;
    particleEnergy.trKin = tKin;

    // potential energy
    // defining reference plane
    wall zeroWall;
    double g(0.0);
    const double gravityNorm = demF.norm();
    if (gravityNorm != 0.0) {
        zeroWall.n = -1.0 * demF / gravityNorm;
        zeroWall.p = tVect(0.0, 0.0, 0.0);
        if (zeroWall.n.dot(Xdirec) < 0) {
            zeroWall.p += tVect(demSize[0], 0.0, 0.0);
        }
        if (zeroWall.n.dot(Y) < 0) {
            zeroWall.p += tVect(0.0, demSize[1], 0.0);
        }
        if (zeroWall.n.dot(Z) < 0) {
            zeroWall.p += tVect(0.0, 0.0, demSize[2]);
        }
        for (int n = 0; n < elmts.size(); n++) {
            const double heigth = zeroWall.dist(elmts[n].x0);
            g += elmts[n].m * heigth*gravityNorm;
        }
    }
    particleEnergy.grav = g;

    // elastic is missing

    // total
    particleEnergy.updateTotal();

}

/////// This is for .data file

void DEM::exportDataFile(IO& io) {
    // export data particles: positions, velocities, angular positions, angular velocity. 

    static int w = io.dataFile.precision() + 5;
    io.dataFile.open(io.dataFileName.c_str(), ios::app);

    double zero = 0.0;

    // set data header in order: N (standParticles + standObjects), time, simulation domain (volume)
    io.dataFile << stdParticles + stdObjects
            << " " << setprecision(5) << demTime
            << " " << setprecision(w) << zero
            << " " << setprecision(w) << zero
            << " " << setprecision(w) << zero
            << " " << setprecision(w) << demSize[0]
            << " " << setprecision(w) << demSize[1]
            << " " << setprecision(w) << demSize[2]
            << "\n";

    tVect angPos = tVect(0.0, 0.0, 0.0); // used to match MercuryDPM outputFiles 
    tVect Obw = tVect(0.0, 0.0, 0.0); // objects angular velocity, used to match MercuryDPM outputFiles 

    for (int o = 0; o < stdObjects; ++o) {
        // export object variables
        objects[o].x0.print(io.dataFile);
        objects[o].x1.print(io.dataFile);
        io.dataFile << objects[o].r << " ";
        io.dataFile << angPos.dot(Xdirec) << " ";
        io.dataFile << angPos.dot(Y) << " ";
        io.dataFile << angPos.dot(Z) << " ";
        io.dataFile << Obw.dot(Xdirec) << " ";
        io.dataFile << Obw.dot(Y) << " ";
        io.dataFile << Obw.dot(Z) << " ";
        io.dataFile << zero << " ";
        io.dataFile << "\n";
    }

    for (int i = 0; i < stdParticles; ++i) {
        // export particle variables
        particles[i].x0.print(io.dataFile);
        particles[i].x1.print(io.dataFile);
        io.dataFile << particles[i].r << " ";
        io.dataFile << angPos.dot(Xdirec) << " ";
        io.dataFile << angPos.dot(Y) << " ";
        io.dataFile << angPos.dot(Z) << " ";
        elmts[particles[i].clusterIndex].wGlobal.print(io.dataFile);
        io.dataFile << zero << " ";
        // dataFile<<dem.particles[i].index<<" "; //keep default 0 because we don't have different species (i.e., sizes), should store information regarding the particle’s material properties (see: http://docs.mercurydpm.org/Beta/d1/d6b/analysing.html).
        io.dataFile << "\n";
    }

    io.dataFile.close();

}

/////// This is for .restart file

void DEM::exportRestartFile(IO& io) {

    double counter = io.currentTimeStep / io.saveCount;

    //  This lines are used to write one single files with restart data overtime
    // static int w = restartFile.precision() + 5;
    io.restartFile.open(io.restartFileName.c_str()); // this is file is overwritten and each time step, (ios::app is NOT used)

    io.restartFile << "restart_version " << 1 << " name" << " " << io.problemNameString
            << "\n"
            << "dataFile " << "fileType " << "ONE_FILE " << "saveCount " << io.saveCount << " " << "counter " << counter << " nextSavedTimeStep " << 0 << "\n"
            << "fStatFile " << "fileType " << "ONE_FILE " << "saveCount " << io.saveCount << " " << "counter " << counter << " nextSavedTimeStep " << 0 << "\n"
            << "eneFile " << "fileType " << "ONE_FILE " << "saveCount " << io.saveCount << " " << "counter " << counter << " nextSavedTimeStep " << 0 << "\n"
            << "restartFile " << "fileType " << "ONE_FILE " << "saveCount " << io.saveCount << " " << "counter " << counter << " nextSavedTimeStep " << 0 << "\n"
            << "statFile " << "fileType " << "ONE_FILE " << "saveCount " << io.saveCount << " " << "counter " << counter << " nextSavedTimeStep " << 0 << "\n"
            << "xMin " << 0 << " xMax " << demSize[0]
            << " yMin " << 0 << " yMax " << demSize[1]
            << " zMin " << 0 << " zMax " << demSize[2]
            << "\n"
            << "timeStep " << setprecision(10) << deltat //statisticExpTime //dem.deltaT   // is it deltaT or just the gap time between two file output??????
            << " time " << demTime // it will be the end-simulation-time as the .restart file is overwritten till the end of simulation
            << " ntimeSteps " << io.currentTimeStep
            << " timeMax " << demTime
            << "\n"
            << "systemDimensions " << 3 // Dimension simulation (Keep default 3)
            << " particleDimensions " << 3
            << " gravity " << demF.dot(Xdirec) << " " << demF.dot(Y) << " " << demF.dot(Z)
            << "\n"
            << "Species " << 1 // Keep default 1
            << "\n"
            << "LinearViscoelasticSlidingFrictionSpecies id " << 0
            << " density " << sphereMat.density
            << " stiffness " << sphereMat.linearStiff //// normal spring constant (linear stiffness) -> unused
            << " dissipation " << sphereMat.dampCoeff //// normal viscous (dampingCoeff )  -> unused
            << " slidingStiffness " << 0 // tangential spring constant (not implemented). In MercuryDPM is set as 2/7 * k(normal)
            << " slidingDissipation " << 0 // tangential viscous constant (not implemented). In MercuryDPM is set as 2/7 * disp(normal) 
            << " slidingFrictionCoefficient " << sphereMat.frictionCoefPart // Coulomb friction (constant)
            << " slidingFrictionCoefficientStatic " << 0
            << "\n"
            << "Walls " << 0 //\todo might write FiniteWall instead (merucryDPM WallHandler.cc - lines 128ith) for the 4 walls you have with normal vector found in ??    
            << "\n"
            //                << "InfiniteWall normal " << -1 << " " << -0 << " " << -0 << " position " << dem.pbcs[1].v.dot(X)     //find it in pbcs.plane 1  
            //                << "\n"
            //                << "InfiniteWall normal " << 1 << " " << 0 << " " << 0 << " position " << dem.pbcs[0].v.dot(X)
            //                << "\n"
            << "Boundaries " << 0
            << "\n"
            << "Particles " << stdParticles + stdObjects
            << "\n";

    tVect angPos = tVect(0.0, 0.0, 0.0); // used to match MercuryDPM output files 
    tVect Obw = tVect(0.0, 0.0, 0.0); // objects angular velocity, used to match MercuryDPM output files 
    double IndSpecies = 0.0;
    double zero = 0.0;

    // Coarse Graining reads the indeces based on the saving order. Objects (fixed particles) must be first. The CG code recognises those because of the inverse mass 

    for (int o = 0; o < stdObjects; ++o) {
        // export object variables
        io.restartFile << "BaseParticle id "; //BP = BaseParticle (see mercuryDPM - ParticleHandler.cc lines 223-256)
        io.restartFile << o << " ";
        io.restartFile << "indSpecies " << 0 << " ";
        io.restartFile << "position ";
        objects[o].x0.print(io.restartFile);
        io.restartFile << "orientation " << 0 << " " << 0 << " " << 0 << " " << 0 << " ";
        io.restartFile << "velocity ";
        objects[o].x1.print(io.restartFile);
        io.restartFile << "angularVelocity ";
        io.restartFile << Obw.dot(Xdirec) << " ";
        io.restartFile << Obw.dot(Y) << " ";
        io.restartFile << Obw.dot(Z) << " ";
        io.restartFile << 0 << " ";
        io.restartFile << "force " << 0 << " " << 0 << " " << 0 << " ";
        io.restartFile << "torque " << 0 << " " << 0 << " " << 0 << " ";
        io.restartFile << "radius ";
        io.restartFile << objects[o].r << " ";
        io.restartFile << "invMass ";
        io.restartFile << zero << " "; //inverse mass
        io.restartFile << "invInertia ";
        io.restartFile << zero << " "; //inverse Inertia
        io.restartFile << "\n";
    }

    for (int i = 0; i < stdParticles; ++i) {
        io.restartFile << "BaseParticle id "; //BP = BaseParticle (see mercuryDPM - ParticleHandler.cc lines 223-256)
        io.restartFile << i + stdObjects << " ";
        io.restartFile << "indSpecies " << 0 << " ";
        io.restartFile << "position ";
        particles[i].x0.print(io.restartFile);
        io.restartFile << "orientation " << 0 << " " << 0 << " " << 0 << " " << 0 << " ";
        io.restartFile << "velocity ";
        particles[i].x1.print(io.restartFile);
        io.restartFile << "angularVelocity ";
        elmts[particles[i].clusterIndex].wGlobal.print(io.restartFile);
        io.restartFile << 0 << " ";
        io.restartFile << "force " << 0 << " " << 0 << " " << 0 << " ";
        io.restartFile << "torque " << 0 << " " << 0 << " " << 0 << " ";
        io.restartFile << "radius ";
        io.restartFile << particles[i].r << " ";
        io.restartFile << "invMass ";
        io.restartFile << 1 / elmts[particles[i].clusterIndex].m << " "; //inverse mass
        io.restartFile << "invInertia ";
        io.restartFile << 1 / elmts[particles[i].clusterIndex].I.dot(Xdirec) << " "; //inverse Inertia
        io.restartFile << "\n";
    }

    io.restartFile << "Interactions " << 0;
    io.restartFile << "\n";

    io.restartFile.close();

}

void DEM::printHeaderStatFile(IO& io) {

    io.statParticleFile.open(io.statisticFileName.c_str(), ios::app);

    double zero = 0.0;

    // set statistic header in order: time, min and max simulation domain (volume), min max particle radii
    io.statParticleFile << "#"
            << " " << setprecision(5) << demTime
            << " " << zero
            << endl;
    io.statParticleFile << "#"
            << " " << zero
            << " " << zero
            << " " << zero
            << " " << demSize[0]
            << " " << demSize[1]
            << " " << demSize[2]
            << endl;
    io.statParticleFile << "#"
            << " " << particles[0].r
            << " " << particles[0].r
            << " " << zero
            << " " << zero
            << " " << zero
            << " " << zero
            << endl;
}

void DEM::saveStatPPcollision(IO& io, const particle* partI, const particle* partJ, const double& overlap, const double& normNormalForce, const tVect& en, const double normTangForce, const tVect et) const {
    // Save data for coarse graining (DEM statistics)
    if (io.saveCount != 0 && io.currentTimeStep % io.saveCount == 0
            && partI->particleIndex < stdParticles && partJ->particleIndex < stdParticles) { // \todo deal better with ghosts

        tVect centre = 0.5 * (partI->x0 + partJ->x0);
        double deltat = 0.0; // length of tangential spring, not yet implemented 

        io.statParticleFile << setprecision(3) << demTime
                << " " << partI->particleIndex + stdObjects
                << " " << partJ->particleIndex + stdObjects
                << " " << setprecision(wSt) << centre.dot(Xdirec)
                << " " << setprecision(wSt) << centre.dot(Y)
                << " " << setprecision(wSt) << centre.dot(Z)
                << " " << setprecision(wSt) << overlap
                << " " << setprecision(wSt) << deltat
                << " " << setprecision(wSt) << abs(normNormalForce)
                << " " << setprecision(wSt) << abs(normTangForce)
                << " " << setprecision(wSt) << en.dot(Xdirec)
                << " " << setprecision(wSt) << en.dot(Y)
                << " " << setprecision(wSt) << en.dot(Z)
                << " " << setprecision(wSt) << et.dot(Xdirec)
                << " " << setprecision(wSt) << et.dot(Y)
                << " " << setprecision(wSt) << et.dot(Z)
                << endl;
    }
}

void DEM::saveStatWPcollision(IO& io, const wall* wallI, const particle* partJ, const double& overlap, const double& normNormalForce, const tVect& en, const double normTangForce, const tVect et) const {

    // Save data for coarse graining (DEM statistics)
    if (io.saveCount != 0 && io.currentTimeStep % io.saveCount == 0
            && partJ->particleIndex < stdParticles) { // \todo deal better with ghosts        

        tVect centre = partJ->x0 - en * (partJ->r - overlap);
        double deltat = 0.0; // length of tangential spring, not yet implemented  
        int wallIdx = wallI->index + 1;

        io.statParticleFile << setprecision(3) << demTime
                << " " << partJ->particleIndex + stdObjects
                << " " << -wallIdx
                << " " << setprecision(wSt) << centre.dot(Xdirec)
                << " " << setprecision(wSt) << centre.dot(Y)
                << " " << setprecision(wSt) << centre.dot(Z)
                << " " << setprecision(wSt) << overlap
                << " " << setprecision(wSt) << deltat
                << " " << setprecision(wSt) << abs(normNormalForce)
                << " " << setprecision(wSt) << abs(normTangForce)
                << " " << setprecision(wSt) << en.dot(Xdirec)
                << " " << setprecision(wSt) << en.dot(Y)
                << " " << setprecision(wSt) << en.dot(Z)
                << " " << setprecision(wSt) << et.dot(Xdirec)
                << " " << setprecision(wSt) << et.dot(Y)
                << " " << setprecision(wSt) << et.dot(Z)
                << endl;
    }
}

void DEM::saveStatOPcollision(IO& io, const object* objectI, const particle* partJ, const double& overlap, const double& normNormalForce, const tVect& en, const double normTangForce, const tVect et) const {

    // Save data for coarse graining (DEM statistics)
    if (io.saveCount != 0 && io.currentTimeStep % io.saveCount == 0
            && partJ->particleIndex < stdParticles) { // \todo deal better with ghosts

        tVect centre = 0.5 * (objectI->x0 + partJ->x0);
        double deltat = 0.0; // length of tangential spring, not yet implemented  

        io.statParticleFile << setprecision(3) << demTime
                << " " << objectI->index //order here is reversed to comply with coarse graining post processing (index should be first, see StatisticsVector.hcc line 1520(ith))
                << " " << partJ->particleIndex + stdObjects
                << " " << setprecision(wSt) << centre.dot(Xdirec)
                << " " << setprecision(wSt) << centre.dot(Y)
                << " " << setprecision(wSt) << centre.dot(Z)
                << " " << setprecision(wSt) << overlap
                << " " << setprecision(wSt) << deltat
                << " " << setprecision(wSt) << abs(normNormalForce)
                << " " << setprecision(wSt) << abs(normTangForce)
                << " " << setprecision(wSt) << en.dot(Xdirec)
                << " " << setprecision(wSt) << en.dot(Y)
                << " " << setprecision(wSt) << en.dot(Z)
                << " " << setprecision(wSt) << et.dot(Xdirec)
                << " " << setprecision(wSt) << et.dot(Y)
                << " " << setprecision(wSt) << et.dot(Z)
                << endl;
    }

}
