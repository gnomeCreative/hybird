#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <csignal>

#include "getpot.h"

#include "IO.h"
#include "DEM.h"
#include "LB.h"

/*
 *      HYBIRD
 *      DEM-LBM coupling code
 *          
 *      author: Alessandro Leonardi
 *
 */

enum ExitCode { UNFINISHED = -1, SUCCESS = 0, TIME_LIMIT_REACHED, SIGNAL_CAUGHT, ERROR };

ExitCode exit_code = UNFINISHED;
ProblemName problemName=NONE;

void catchSignal(int sig) {
  std::cout << "Received signal " << sig << ": " << strsignal(sig) << " " << std::endl;
  exit_code = (SIGNAL_CAUGHT > exit_code ? SIGNAL_CAUGHT : exit_code);
}

void print_help() {
  cout << "Usage: " << " [-c configFile] [-d datadir] [-n resultFolderName (opt. time)] [options}\n";
}

void goCycle(IO& io, DEM& dem, LB& lb) {
    
    // advance one step in time
    io.outputStep(dem.demTime, lb, dem);

    ++io.currentTimeStep;
    ++lb.time;
    dem.demTime+=lb.unit.Time;

    if (io.demSolve) {
        dem.discreteElementStep(io);
    }

    if (io.lbmSolve && dem.demTime>dem.demInitialRepeat) {

        if (lb.freeSurface) {
            lb.latticeBoltzmannFreeSurfaceStep();
    //                lb.evolveBoundary();
        }
        if (io.demSolve) {
            lb.latticeBoltzmannCouplingStep(dem.newNeighborList, dem.elmts, dem.particles);
        }
        
        lb.latticeBolzmannStep(dem.elmts, dem.particles, dem.walls);
        
    }
    else if (io.lbmSolve && dem.demTime<=dem.demInitialRepeat) {
        if (io.demSolve) {
            lb.latticeBoltzmannCouplingStep(dem.newNeighborList, dem.elmts, dem.particles);
        }
    }

}

void parseCommandLine(IO& io, GetPot& command_line) {

    //io.base_directory = command_line.follow("./", "--directory");

    // print help and exit if requested
    if (command_line.search("-h") || command_line.search("--help")) {
        print_help();
        exit(0);
    }

    // create the data directory in case it doesn't exist
    io.workDirectory = "data";
    if (command_line.search("-d")) {
        io.workDirectory = command_line.next(io.workDirectory);
    }
    
    // if the name of the simulation is specified create subfolder
    std::string programName = "name";
    if (command_line.search("-n")) {
        // creating work directory
        programName=command_line.next(programName);
        // if the specified name is "time" then create a folder named with the initial time
        if (programName=="time") {
            ostringstream convertTime;
            convertTime<<std::setfill('0') << std::setw(4) <<(io.now->tm_year + 1900)
                    <<std::setfill('0') << std::setw(2) <<(io.now->tm_mon+1)
                    <<std::setfill('0') << std::setw(2) <<io.now->tm_mday
                    <<"_"
                    <<std::setfill('0') << std::setw(2) <<io.now->tm_hour
                    <<std::setfill('0') << std::setw(2) <<io.now->tm_min
                    <<std::setfill('0') << std::setw(2) <<io.now->tm_sec;
            string simulationTime=convertTime.str();
            io.workDirectory=io.workDirectory+"/"+simulationTime;
            io.workDirectory=io.workDirectory.substr(0,io.workDirectory.size());
        }
        else {
            io.workDirectory=io.workDirectory+"/"+programName;
        }
        
//        io.workDirectory = command_line.next(io.workDirectory);
    }

    // lbm configuration file
    io.lbmCfgName = "config.cfg";
    if (command_line.search("-c")) {
        io.lbmCfgName = command_line.next(io.lbmCfgName.c_str());
        cout<<"Using "<<io.lbmCfgName<<" for LBM and DEM parameters\n";
        // copy the lbm configuration file into the data folder for later reproducibility
        int a;
        a=mkdir(io.workDirectory.c_str(),0777);
        cout<<"Work directory created = "<<io.workDirectory<<". Result: "<<a<<"\n";

        std::system(("cp '" + io.lbmCfgName + "' '" + io.workDirectory + "'").c_str());
    }

    // make sure the config files can be read
    std::ifstream lbmCfgFile(io.lbmCfgName.c_str());
    if (!lbmCfgFile) {
        cout << "ERROR: Can't open config file \"" << io.lbmCfgName << "\" for reading!\n";
//        return ERROR;
    }
    lbmCfgFile.close();
}

void parseConfigFile(IO& io, DEM& dem, LB& lb, GetPot& lbmCfgFile, GetPot& command_line) {

    // PROBLEM NAME //////////////
    // necessary for hard coded sections of the code
    PARSE_CLASS_MEMBER(lbmCfgFile, io.problemNameString, "problemName","none");
    if (io.problemNameString=="DRUM") problemName=DRUM;
    else if (io.problemNameString=="SHEARCELL") problemName=SHEARCELL;
    else if (io.problemNameString=="NONE") problemName=NONE;
    else if (io.problemNameString=="AVALANCHE") problemName=AVALANCHE;
    else if (io.problemNameString=="SPLASH") problemName=SPLASH;
    else if (io.problemNameString=="BOX") problemName=BOX;
    else if (io.problemNameString=="NET") problemName=NET;
    else if (io.problemNameString=="DIFF") problemName=DIFF;
    else if (io.problemNameString=="BARRIER") problemName=BARRIER;
    else if (io.problemNameString=="demChute") problemName=demChute;
    // else if (problemNameString=="SETT") problemName=SETT;


    // GETTING SIMULATION PARAMETERS  /////////
    // DEM initial iterations
    PARSE_CLASS_MEMBER(lbmCfgFile, dem.demInitialRepeat, "demInitialRepeat",0.0);
    ASSERT(dem.demInitialRepeat >= 0);
    // LB iteration without coupling (for initial stability) - > time is frozen here
    PARSE_CLASS_MEMBER(lbmCfgFile, lb.lbmInitialRepeat, "lbmInitialRepeat",0);
    ASSERT(lb.lbmInitialRepeat >= 0);
    // maximum time variable value
    PARSE_CLASS_MEMBER(lbmCfgFile, io.maxTime, "maxTime", 0.0);
    ASSERT(io.maxTime >= 0);

    // for time integration and output
    PARSE_CLASS_MEMBER(lbmCfgFile, io.maximumTimeSteps, "maximumTimeSteps", 0);
    ASSERT(io.maximumTimeSteps>=0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.screenExpTime, "screenExpTime", 0.0);
    ASSERT(io.screenExpTime>=0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.fluidExpTime, "fluidExpTime", 0.0);
    ASSERT(io.fluidExpTime>=0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.partExpTime, "partExpTime", 0.0);
    ASSERT(io.partExpTime>=0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.objectExpTime, "objectExpTime", 0.0);
    ASSERT(io.objectExpTime>=0);
    PARSE_CLASS_MEMBER(lbmCfgFile, io.recycleExpTime, "recycleExpTime", 0.0);
    ASSERT(io.recycleExpTime>=0);
    // saveCount: how namy time steps are skipped before .restart , .data, .fstat files are saved 
    PARSE_CLASS_MEMBER(lbmCfgFile, io.saveCount, "saveCount", 0.0);
    ASSERT(io.saveCount>= 0);


    // solving particles?
    PARSE_CLASS_MEMBER(lbmCfgFile, io.demSolve, "demSolve",0);
    // solving fluid?
    PARSE_CLASS_MEMBER(lbmCfgFile, io.lbmSolve, "lbSolve",0);
    // solving free surface?
    PARSE_CLASS_MEMBER(lbmCfgFile, lb.freeSurface, "freeSurfaceSolve",0);
    // is there a force field? (shutting down improves performances)
    PARSE_CLASS_MEMBER(lbmCfgFile, lb.forceField, "forceFieldSolve",0);
    // is the fluid non-Newtonian? (shutting down improves performances)
    PARSE_CLASS_MEMBER(lbmCfgFile, lb.nonNewtonian, "nonNewtonianSolve",0);
    // are we using turbulence modeling? (shutting down improves performances)
    PARSE_CLASS_MEMBER(lbmCfgFile, lb.turbulenceOn, "turbulenceSolve",0);

    // getting LBM parameters
    lb.latticeBoltzmannGet(lbmCfgFile, command_line);
    // getting DEM parameters and particle initial status
    dem.discreteElementGet(lbmCfgFile, command_line);

    // problem specific parameters
    switch (problemName) {
        case DRUM: {
            // drum rotational speed
            PARSE_CLASS_MEMBER(lbmCfgFile, dem.drumSpeed, "drumSpeed", 0.0);
            cout<<"DRUM SPEED ="<<dem.drumSpeed<<"\n";
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.fluidMass, "fluidMass",0.0);
            cout<<"INPUT FLUID MASS="<<lb.fluidMass<<"\n";
            break;
        }
        case SHEARCELL: {
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.maxVisc, "maxVisc", 0.0);
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.maxPlasticVisc, "maxPlasticVisc", 0.0);
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.maxYieldStress, "maxYieldStress", 0.0);
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.viscSteps, "viscSteps", 0.0);
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.maxShearRate, "maxShearRate", 0.0);
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.shearRateSteps, "shearRateSteps", 0.0);
            break;
        }
        case NET:
        case BARRIER: {
            PARSE_CLASS_MEMBER(lbmCfgFile, lb.avalanchePosit, "avalanchePosit", 0.0);
            break;
        }

    }
}

void printUfo(GetPot& command_line, GetPot& lbmCfgFile){
        // warn about unused parameters
    std::vector<std::string> ufos=lbmCfgFile.unidentified_variables();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized lbm config file parameter(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
    ufos = command_line.unidentified_arguments();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized command line argument(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {

    // Checking number of processes involved
    cout<<"Program starts with threads:";
    #pragma omp parallel
    {
        #pragma omp critical
        cout<<" X ";
    }
    cout<<"\n";

    // DECLARATION OF VARIABLES - Input-Output ///////////////
    IO io;

    // DECLARATION OF VARIABLES - DEM ///////////////
    DEM dem;

    // DECLARATION OF VARIABLES - LB ///////////////
    LB lb;
    
    // print some info for restorability
    time_t t=time(0);   // get time now
    io.now=localtime(&t);
    cout<<"Binary was compiled on "<<__DATE__<<" at "<<__TIME__<<endl;
    cout<<"Invocation on "<<(io.now->tm_year + 1900)<<"-"<<(io.now->tm_mon+1)<<"-"<<io.now->tm_mday<<
            ", at "<<io.now->tm_hour<<":"<<io.now->tm_min<<":"<<io.now->tm_sec<<endl;
    cout<<"Invocation line: " << argv[0];
    for (int i=1; i<argc; ++i) {
        cout<<" "<<argv[i];
    }
    cout<<endl;

    // parsing command line and configuration file
    GetPot command_line(argc, argv);
    parseCommandLine(io, command_line);

    // parsing LBM input file
    GetPot lbmCfgFile(io.lbmCfgName);
    parseConfigFile(io, dem, lb, lbmCfgFile, command_line);
    
    printUfo(command_line, lbmCfgFile);

    // bind signal handlers
    signal(SIGHUP, SIG_IGN); // for easy running through ssh
    signal(SIGQUIT, catchSignal);
    signal(SIGTERM, catchSignal); // important for job killing on a cluster

    io.currentTimeStep = 0;

    /* /////////////////////////////////////////////////////////////////////////
    // PROGRAM CORE  ///////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////*/
    cout<<"PROBLEM ID: "<<problemName<<"\n";

    //  GLOBAL INITIALIZATION ///////////////////////////////////////////////
    // initializing input-output
    io.initialize();
    
    // initializing lattice
    lb.latticeDefinition();
    lb.LBShow();

    // initializing DEM parameters
    dem.discreteElementInit(lb.boundary, lb.lbSize, lb.unit, lb.lbF);
    if (io.lbmSolve) {
        lb.latticeBolzmannInit(dem.cylinders, dem.walls, dem.particles, dem.objects);
    }

    // setting time
    lb.time=0;
    dem.demTime=0.0;

    // CYCLE /////////////////////////////
    // integrate in time
    while (true) {
        
        if (io.realTime!=0.0 && io.realTime>io.maxTime) {
            exit_code = SUCCESS;
        }
        // exit normally if the maximum simulation time has been reached
        else if (io.maximumTimeSteps && io.currentTimeStep>=io.maximumTimeSteps) {
            exit_code = SUCCESS;
        }
        else {
            // core of the code, performs time steps
            goCycle(io, dem, lb);

//            // exit abnormally if a serious problem has occurred
//            if (io.problem) {
//                exit_code = ERROR;
//            }
        }

        if (exit_code > UNFINISHED) {
            break;
        }
    }
  return exit_code;
}
