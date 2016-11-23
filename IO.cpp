
#include "IO.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/


void IO::initialize(){

    partDirectory=workDirectory+"/particleData";
    fluidDirectory=workDirectory+"/fluidData";

    int a;
    if (demSolve) {
        a=mkdir(partDirectory.c_str(),0777);
        cout<<"Work directory created = "<<partDirectory<<". Result: "<<a<<"\n";
    }
    if (lbmSolve) {
        a=mkdir(fluidDirectory.c_str(),0777);
        cout<<"Work directory created = "<<fluidDirectory<<". Result: "<<a<<"\n";
    }

    
    //  initialising energy file
    energyFileName=workDirectory+"/energy.dat";
    static int width = energyFile.precision() + 5;
    energyFile.open(energyFileName.c_str(), ios::app);
    // Set energyFile header 
    energyFile  << "\t" << " time" 
                << "\t" << setw(width) << "ene_grav" 
                << "\t" << setw(width) << "ene_kin" 
                << "\t" << setw(width) << "ene_rot" 
                << "\t" << setw(width) << "ene_ela" 
                << "\t" << setw(width) << "X_com" 
                << "\t" << setw(width) << "Y_com" 
                << "\t" << setw(width) << "Z_com" 
                << "\n";
    energyFile.close();
	
     //  initialising .restart file
    restartFileName=workDirectory+"/"+problemNameString+".restart";
    
    //  initialising .data file 
    dataFileName=workDirectory+"/"+problemNameString+".data";
    
    //  initialising .fstat file
    statisticFileName=workDirectory+"/"+problemNameString+".fstat";
            
            
    switch (problemName) {
        case DRUM: {
            // create folder
            string drumDirectory=workDirectory+"/drumData";
            a=mkdir(drumDirectory.c_str(),0777);
            cout<<"Work directory created = "<<drumDirectory<<". Result: "<<a<<"\n";

            //  initializing drum files
            centerOfMassFileName=drumDirectory+"/centerOfMass.dat";
            positionFileName=drumDirectory+"/fluidPosition.dat";
            particleAngleFileName=drumDirectory+"/particleAngles.dat";
            particleNumberFileName=drumDirectory+"/particleNumber.dat";
            collisionNumberFileName=drumDirectory+"/particleCollisionNumber.dat";
            hydraulicNumberFileName=drumDirectory+"/particleHydraulicNumber.dat";
            collisionForceFileName=drumDirectory+"/particleCollisionForce.dat";
            hydraulicForceFileName=drumDirectory+"/particleHydraulicForce.dat";
            profileAngleFileName=drumDirectory+"/profileAngles.dat";
            fluidLevelFileName=drumDirectory+"/fluidLevelProfile.dat";
            particleLevelFileName=drumDirectory+"/particleLevelProfile.dat";
            pressureFileName=drumDirectory+"/pressureProfile.dat";
            fluidSpeedFileName=drumDirectory+"/fluidSpeedProfile.dat";
            particleSpeedFileName=drumDirectory+"/particleSpeedProfile.dat";
            particleSurfaceFileName=drumDirectory+"/particleSurfaceProfile.dat";
            fluidSurfaceFileName=drumDirectory+"/fluidSurfaceProfile.dat";


            // DRUM POSITIONS

            //  initializing position file
            positionFileName=drumDirectory+"/positions.dat";

            break;
        }
        case AVALANCHE: {
            // create folder
            string avalancheDirectory=workDirectory+"/avalancheData";
            a=mkdir(avalancheDirectory.c_str(),0777);
            cout<<"Work directory created = "<<avalancheDirectory<<". Result: "<<a<<"\n";

            //  initializing drum files
            obstacleFileName=avalancheDirectory+"/obstacleForce.dat";
            break;
        }
    }


    // build data file names 
    fluidFileFormat=fluidDirectory+"/fluid%010u.vti";
    partFileFormat=partDirectory+"/part%010u.vtu";
    recycleFileFormat=partDirectory+"/recycle%010u.vtu";
    objectFileFormat=partDirectory+"/object%010u.vtu";

    //  initializing output file
    exportFileName=workDirectory+"/export.dat";

    //  initializing max speed file
    maxSpeedFileName=workDirectory+"/maxVel.dat";
    maxSpeedFile.open(maxSpeedFileName.c_str(), ios::app);
    maxSpeedFile<<"time maxFluidSpeed\n";
    maxSpeedFile.close();

    //  initializing force file
    forceFileName=workDirectory+"/force.dat";
    forceFile.open(forceFileName.c_str(), ios::app);
    forceFile<<"time collisionForces hydraulicForces\n";
    forceFile.close();

    //  initializing plasticity file
    plasticityFileName=workDirectory+"/plasticity.dat";
    plasticityFile.open(plasticityFileName.c_str(), ios::app);
    plasticityFile<<"time percPlastic\n";
    plasticityFile.close();

    lastScreenExp=0;
    lastFluidExp=0;
    lastPartExp=0;
    lastObjectExp=0;
    lastOutputExp=0;
    lastRecycleExp=0;
}

void IO::outputStep(const double& demTime, const LB& lb, const DEM& dem){

    realTime=demTime;

    // SPECIAL: the analysis of drum particles requires this to be run every time step :(
    switch (problemName) {
        case DRUM: {
            updateDistributions(dem, dem.cylinders[0]);
            break;
        }
    }

    // FILE CREATION PHASE  ////////////////////////////////////////////////////////////////
    createParaviewFiles(lb, dem);

    // PLOTTING PHASE  ////////////////////////////////////////////////////////////////
    const unsigned int screenExpCounter = (screenExpTime > 0 ? static_cast<unsigned int>(realTime / screenExpTime)+1 : 0);
    if (screenExpCounter>lastScreenExp) {
        lastScreenExp = screenExpCounter;
        
        exportFile.open(exportFileName.c_str(), ios::app);

        // current iteration and real time
        cout<<currentTimeStep<<"; time="<<realTime<<"\t";
        exportFile<<currentTimeStep<<"; time="<<realTime<<"\t";

        if (lbmSolve) {
            exportMaxSpeedFluid(lb);
            exportTotalMass(lb);
        }

        if (lbmSolve & lb.nonNewtonian) {
            exportPlasticity(lb);
        }

        if (dem.elmts.size()) {
            exportMaxSpeedParticles(dem);
            exportForces(dem);
            //exportDiffusion(dem);
        }

        if (dem.elmts.size() && !lbmSolve) {
            exportEnergy(dem);
        }

        switch (problemName) {
            case DRUM: {
                exportDrum(lb, dem);
                break;
            }
            case SHEARCELL: {
                exportShearCell(lb, dem);
                break;
            }
            case AVALANCHE: {
                exportForceObstacle(dem);
                break;
            }
        }

        // closing file
        cout<<"\n";
        exportFile<<"\n";
        exportFile.close();
        cout.flush();
    }
}

void IO::outputFinal(){
    // drag file closure
    exportFile.close();
}

void IO::initViscFile(){
    //  initializing viscosity file
    stringstream sf;
    sf<<currentTimeStep;
    viscFileName=workDirectory+"/viscosity"+sf.str()+".dat";
    viscFile.open(viscFileName.c_str(), ios::app);
    viscFile<<"Visc = "<<currentTimeStep<<"\n";
    viscFile<<"shearRate shearStress apparentViscosity Zdiff\n";
    viscFile.close();
}

void IO::writeViscFile(const LB& lb, const wallList& walls, const elmtList& elmts) {
    // calculating apparent viscosity
    double appVisc=0.0;
    double externalShear=0.0;
    double wallStress=0.0;
    apparentViscosity(lb, walls, externalShear, wallStress, appVisc);

    tVect msdVecAv=msdVecSum/double(msdSumN);
    double msdVecAvX=msdVecAv.dot(Xdirec);
    double msdVecAvY=msdVecAv.dot(Y);
    double msdVecAvZ=msdVecAv.dot(Z);

    tVect msdVecSigma=tVect( sqrt( msdVecSum2.dot(Xdirec)/double(msdSumN)-msdVecAvX*msdVecAvX ),
                                                    sqrt( msdVecSum2.dot(Y)/double(msdSumN)-msdVecAvY*msdVecAvY ),
                                                    sqrt( msdVecSum2.dot(Z)/double(msdSumN)-msdVecAvZ*msdVecAvZ ) );

    // writing on viscosity file
    viscFile.open(viscFileName.c_str(), ios::app);
    viscFile<<externalShear<<" "<<wallStress<<" "<<appVisc<<" "<<msdVecAv.dot(Z)/realTime<<" "<<msdVecSigma.dot(Z)/realTime<<"\n";
    viscFile.close();

    // killing averages
    msdSum=msdSum2=0.0;
    msdSumN=0;
    msdVecSum.reset();
    msdVecSum2.reset();;
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// paraview files

void IO::createParaviewFiles(const LB& lb, const DEM& dem) {
    // write vtk at regular interval defined with the input file
    
    if (lbmSolve) {
        const unsigned int fluidExpCounter=(fluidExpTime>0 ? static_cast<unsigned int>(realTime/fluidExpTime)+1 : 0);
        if (fluidExpCounter>lastFluidExp) {
            lastFluidExp=fluidExpCounter;
            sprintf(filePathBuffer, fluidFileFormat.c_str(), currentTimeStep);
            //if (currentTimeStep>120) {
            exportParaviewFluidOld(lb, filePathBuffer);
            //}
        }
    }

    if (demSolve) {
        
        const unsigned int partExpCounter=(partExpTime>0 ? static_cast<unsigned int>(realTime/partExpTime)+1 : 0);
        if (partExpCounter>lastPartExp) {
            
            lastPartExp=partExpCounter;
            sprintf(filePathBuffer, partFileFormat.c_str(), currentTimeStep);
            exportParaviewParticles(dem.elmts, dem.particles, filePathBuffer);
        }
        const unsigned int recycleExpCounter=(recycleExpTime>0 ? static_cast<unsigned int>(realTime/recycleExpTime)+1 : 0);
        if (recycleExpCounter>lastRecycleExp) {
            lastRecycleExp=recycleExpCounter;
            sprintf(filePathBuffer, partFileFormat.c_str(), currentTimeStep);
            // requires the pbcShift, contained int he neighborList function
            exportRecycleParticles(dem.elmts, dem.pbcs);
        }
        const unsigned int objectExpCounter=(objectExpTime>0 ? static_cast<unsigned int>(realTime/objectExpTime)+1 : 0);
        if (objectExpCounter>lastObjectExp) {
            lastObjectExp=objectExpCounter;
            sprintf(filePathBuffer, objectFileFormat.c_str(), currentTimeStep);
            exportParaviewObjects(dem.objects, filePathBuffer);
        }
    }
}

void IO::exportParaviewParticles(const elmtList& elmts, const particleList& particles, const string& particleFile){

    const int one=1;
    const int Pnumber=particles.size();
    const tVect nx(1.0,0.0,0.0),ny(0.0,1.0,0.0),nz(0.0,0.0,1.0);
    tVect n1,n2,n3;

    // file opening
    ofstream paraviewParticleFile;
    paraviewParticleFile.open(particleFile.c_str());
    // writing on header file
    paraviewParticleFile<<"<?xml version=\"1.0\"?>\n";
    paraviewParticleFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewParticleFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewParticleFile<<"  <Piece NumberOfPoints=\""<<Pnumber<<"\" NumberOfCells=\""<<Pnumber<<"\">\n";
    paraviewParticleFile<<"   <PointData>\n";
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"radius\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        paraviewParticleFile<<particles[i].r<<"\n";
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"mass\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        paraviewParticleFile<<elmts[particles[i].clusterIndex].m<<"\n";
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"particleIndex\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        paraviewParticleFile<<particles[i].particleIndex<<"\n";
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"clusterIndex\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        paraviewParticleFile<<particles[i].clusterIndex<<"\n";
    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"tableCell\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        paraviewParticleFile<<particles[i].tableCell<<"\n";
//    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"Inertia\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        elmts[particles[i].clusterIndex].I.printLine(paraviewParticleFile);
//    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        particles[i].x1.printLine(paraviewParticleFile);
    }
        paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"w\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].wGlobal.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].FParticle.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FGrav\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].FGrav.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"MParticle\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].MParticle.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FWall\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].FWall.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"MWall\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].MWall.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FHydro\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        elmts[particles[i].clusterIndex].FHydro.printLine(paraviewParticleFile);
    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"MHydro\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        elmts[particles[i].clusterIndex].MHydro.printLine(paraviewParticleFile);
//    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FLubParticle\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        elmts[particles[i].clusterIndex].FLubParticle.printLine(paraviewParticleFile);
//    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"FLubWall\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        elmts[particles[i].clusterIndex].FLubWall.printLine(paraviewParticleFile);
//    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"n1\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        n1=project(nx,elmts[particles[i].clusterIndex].q0);
//        n1.printLine(paraviewParticleFile);
//    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"n2\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        n2=project(ny,elmts[particles[i].clusterIndex].q0);
//        n2.printLine(paraviewParticleFile);
//    }
//    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"n3\" NumberOfComponents=\"3\"/>\n";
//    for (int i=0; i<Pnumber; ++i){
//        n3=project(nz,elmts[particles[i].clusterIndex].q0);
//        n3.printLine(paraviewParticleFile);
//    }
    paraviewParticleFile<<"   </PointData>\n";
    paraviewParticleFile<<"   <CellData>\n";
    paraviewParticleFile<<"   </CellData>\n";
    paraviewParticleFile<<"   <Points>\n";
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Pnumber; ++i){
        particles[i].x0.printLine(paraviewParticleFile);
    }
    paraviewParticleFile<<"   </Points>\n";
    paraviewParticleFile<<"   <Cells>\n";
    paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i=1; i<Pnumber+1; ++i){
        paraviewParticleFile<<i-1<<"\n";
    }
    paraviewParticleFile<<"    </DataArray>\n";
    paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i=1; i<Pnumber+1; ++i){
        paraviewParticleFile<<i<<"\n";
    }
    paraviewParticleFile<<"    </DataArray>\n";
    paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i=0; i<Pnumber; ++i){
        paraviewParticleFile<<one<<"\n";
    }
    paraviewParticleFile<<"    </DataArray>\n";
    paraviewParticleFile<<"   </Cells>\n";
    paraviewParticleFile<<"  </Piece>\n";
//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewParticleFile<<" </UnstructuredGrid>\n";
    paraviewParticleFile<<"</VTKFile>";
    // header file closing
    paraviewParticleFile.close();
}

void IO::exportParaviewObjects(const objectList& objects, const string& objectFile){

    const int one=1;
    const int Onumber=objects.size();

    // file opening
    ofstream paraviewObjectFile;
    paraviewObjectFile.open(objectFile.c_str());
    // writing on header file
    paraviewObjectFile<<"<?xml version=\"1.0\"?>\n";
    paraviewObjectFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewObjectFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewObjectFile<<"  <Piece NumberOfPoints=\""<<Onumber<<"\" NumberOfCells=\""<<Onumber<<"\">\n";
    paraviewObjectFile<<"   <PointData>\n";
    paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"radius\"/>\n";
    for (int i=0; i<Onumber; ++i){
        paraviewObjectFile<<objects[i].r<<"\n";
    }
    paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"objectIndex\"/>\n";
    for (int i=0; i<Onumber; ++i){
        paraviewObjectFile<<objects[i].index<<"\n";
    }
    paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Onumber; ++i){
        objects[i].x1.printLine(paraviewObjectFile);
    }
    paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"Force\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Onumber; ++i){
        objects[i].FParticle.printLine(paraviewObjectFile);
    }
    paraviewObjectFile<<"   </PointData>\n";
    paraviewObjectFile<<"   <CellData>\n";
    paraviewObjectFile<<"   </CellData>\n";
    paraviewObjectFile<<"   <Points>\n";
    paraviewObjectFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (int i=0; i<Onumber; ++i){
        objects[i].x0.printLine(paraviewObjectFile);
    }
    paraviewObjectFile<<"   </Points>\n";
    paraviewObjectFile<<"   <Cells>\n";
    paraviewObjectFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i=1; i<Onumber+1; ++i){
        paraviewObjectFile<<i-1<<"\n";
    }
    paraviewObjectFile<<"    </DataArray>\n";
    paraviewObjectFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i=1; i<Onumber+1; ++i){
        paraviewObjectFile<<i<<"\n";
    }
    paraviewObjectFile<<"    </DataArray>\n";
    paraviewObjectFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i=0; i<Onumber; ++i){
        paraviewObjectFile<<one<<"\n";
    }
    paraviewObjectFile<<"    </DataArray>\n";
    paraviewObjectFile<<"   </Cells>\n";
    paraviewObjectFile<<"  </Piece>\n";
//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewObjectFile<<" </UnstructuredGrid>\n";
    paraviewObjectFile<<"</VTKFile>";
    // header file closing
    paraviewObjectFile.close();
}

void IO::exportRecycleParticles(const elmtList& elmts, const pbcList& pbcs){
    // exports a format readable by this code itself, allowing the re-use of computed particle data.
    // needs further work...

    int i;
    const char* charParticleFile;
    ofstream recycleParticleFile;
    string title, particleFile;
    title=partDirectory+"/recycle_";
    stringstream ss;
    ss<<currentTimeStep;
    particleFile=title+ss.str()+".dat";
    charParticleFile=particleFile.c_str();

    double zero=0.0;
    // header file opening
    recycleParticleFile.open(charParticleFile);

    // total number of element

    recycleParticleFile<<elmts.size()<<"\n";

    // ATTENTION this does not work for clusters unless also rotational velocity and position are included in the initial file :(
    for (int i=0; i<elmts.size(); ++i) {
        // import variables
        recycleParticleFile<<elmts[i].index<<" ";
        recycleParticleFile<<elmts[i].size<<" ";
        recycleParticleFile<<elmts[i].radius<<" ";
        // element center could be out of domain of there are periodic walls, fix this.
        tVect printPosition=elmts[i].x0;
        for (int b=0; b<pbcs.size(); ++b) {
            const tVect pbcVector=pbcs[b].v;
            const double leftDist=pbcs[b].pl1.dist(elmts[i].x0);
            const double rightDist=pbcs[b].pl2.dist(elmts[i].x0);
            // first plane of couple (left)
            if (leftDist<0.0) {
                printPosition=printPosition+pbcVector;
            }
            if (rightDist<0.0) {
                printPosition=printPosition-pbcVector;
            }
        }
        printPosition.print(recycleParticleFile);
        elmts[i].x1.print(recycleParticleFile);
        elmts[i].w0.print(recycleParticleFile);
        elmts[i].q0.print(recycleParticleFile);
        elmts[i].q1.printLine(recycleParticleFile);
   }
    recycleParticleFile.close();
}

void IO::exportParaviewBox(const double& time){
    int i;
    const char* charParticleFile;
    ofstream paraviewParticleFile;
    string title, particleFile;
    title=partDirectory+"/box_";
    stringstream ss;
    ss<<currentTimeStep;
    particleFile=title+ss.str()+".vtu";
    charParticleFile=particleFile.c_str();

    int one=1;
    int Pnumber=1;
    double boxVel=0.15/1.5;
    // start printing all the crap required for Paraview
    // header file opening
    paraviewParticleFile.open(charParticleFile);
    // writing on header file
    paraviewParticleFile<<"<?xml version=\"1.0\"?>\n";
    paraviewParticleFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewParticleFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewParticleFile<<"  <Piece NumberOfPoints=\""<<Pnumber<<"\" NumberOfCells=\""<<Pnumber<<"\">\n";
    paraviewParticleFile<<"   <PointData>\n";
    paraviewParticleFile<<"   </PointData>\n";
    paraviewParticleFile<<"   <CellData>\n";
    paraviewParticleFile<<"   </CellData>\n";
    paraviewParticleFile<<"   <Points>\n";
    paraviewParticleFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (i=0; i<Pnumber; ++i){
        paraviewParticleFile<<0.17<<" "<<0.25<<" "<<0.08+boxVel*time<<"\n"; //third component to be added in the future
    }
    paraviewParticleFile<<"   </Points>\n";
    paraviewParticleFile<<"   <Cells>\n";
    paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (i=1; i<Pnumber+1; ++i){
        paraviewParticleFile<<i-1<<"\n";
    }
    paraviewParticleFile<<"    </DataArray>\n";
    paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (i=1; i<Pnumber+1; ++i){
        paraviewParticleFile<<i<<"\n";
    }
    paraviewParticleFile<<"    </DataArray>\n";
    paraviewParticleFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (i=0; i<Pnumber; ++i){
        paraviewParticleFile<<one<<"\n";
    }
    paraviewParticleFile<<"    </DataArray>\n";
    paraviewParticleFile<<"   </Cells>\n";
    paraviewParticleFile<<"  </Piece>\n";
//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewParticleFile<<" </UnstructuredGrid>\n";
    paraviewParticleFile<<"</VTKFile>";
    // header file closing
    paraviewParticleFile.close();
}

void IO::exportParaviewFluidNew(const LB& lb){
    int i;
    int one=1;
    const char* charFluidFile;
    ofstream paraviewFluidFile;
    string title, fluidFile;
    title=fluidDirectory+"/fluid_";
    stringstream ss;
    ss<<currentTimeStep;
    fluidFile=title+ss.str()+".vtu";
    charFluidFile=fluidFile.c_str();
    unsigned int it=0;
    int Pnumber=lb.activeNodes.size();
    cout<<Pnumber<<"\n";
    // start printing all the crap required for Paraview
    // header file opening
    paraviewFluidFile.open(charFluidFile);
    // writing on header file
    paraviewFluidFile<<"<?xml version=\"1.0\"?>\n";
    paraviewFluidFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewFluidFile<<" <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewFluidFile<<"  <Piece NumberOfPoints=\""<<Pnumber<<"\" NumberOfCells=\""<<Pnumber<<"\">\n";
    paraviewFluidFile<<"   <PointData>\n";
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        tVect uPhysic=lb.nodes[it]->u*lb.unit.Speed;
        uPhysic.printLine(paraviewFluidFile);
    }
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        if (lb.nodes[it]->n!=0.0) {
            paraviewFluidFile<<0.3333333*(lb.nodes[it]->n-1.0)*lb.unit.Pressure<<"\n";
        }
        else {
            paraviewFluidFile<<lb.nodes[it]->n<<"\n";
        }
    }
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"viscosity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        paraviewFluidFile<<lb.nodes[it]->visc*lb.unit.KinVisc<<"\n";
    }
//    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"isInit\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
//    for (i=0; i<Pnumber; ++i){
//
//        unsigned int a=0;
//        if (lb.nodes[i]!=0) {
//            a=1;
//        }
//        paraviewFluidFile<<a<<"\n";
//    }
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"shearRate\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        paraviewFluidFile<<lb.nodes[it]->shearRate*lb.unit.AngVel<<"\n";
    }
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        paraviewFluidFile<<lb.nodes[it]->mass*lb.unit.Density<<"\n";
    }
    paraviewFluidFile<<"    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        paraviewFluidFile<<lb.types[it].getType()<<"\n";
    }
    paraviewFluidFile<<"   </PointData>\n";
    paraviewFluidFile<<"   <CellData>\n";
    paraviewFluidFile<<"   </CellData>\n";
    paraviewFluidFile<<"   <Points>\n";
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (i=0; i<Pnumber; ++i){
        it=lb.activeNodes[i];
        lb.positions[it].printLine(paraviewFluidFile);
    }
    paraviewFluidFile<<"   </Points>\n";
    paraviewFluidFile<<"   <Cells>\n";
    paraviewFluidFile<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (i=1; i<Pnumber+1; ++i){
        paraviewFluidFile<<i-1<<"\n";
    }
    paraviewFluidFile<<"    </DataArray>\n";
    paraviewFluidFile<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (i=1; i<Pnumber+1; ++i){
        paraviewFluidFile<<i<<"\n";
    }
    paraviewFluidFile<<"    </DataArray>\n";
    paraviewFluidFile<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (i=0; i<Pnumber; ++i){
        paraviewFluidFile<<one<<"\n";
    }
    paraviewFluidFile<<"    </DataArray>\n";
    paraviewFluidFile<<"   </Cells>\n";
    paraviewFluidFile<<"  </Piece>\n";
//    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewFluidFile<<" </UnstructuredGrid>\n";
    paraviewFluidFile<<"</VTKFile>";

    // data file closing
    paraviewFluidFile.close();
}

void IO::exportParaviewFluidOld(const LB& lb, const string& fluidFile){
    int i;
    double zero=0.0;

//    const char* charFluidFile;
//    
//    string title, fluidFile;
//    title=fluidDirectory+"/fluid_";
//    stringstream ss;
//    ss<<currentTimeStep;
//    fluidFile=title+ss.str()+".vti";
//    charFluidFile=fluidFile.c_str();

    int Pnumber=lb.totNodes;
    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;
    paraviewFluidFile.open(fluidFile.c_str());
    // writing on header file
    paraviewFluidFile<<"<VTKFile type=\"ImageData\" version=\"0.1\">\n";
    paraviewFluidFile<<" <ImageData WholeExtent=\"0 "<<lb.lbSize[0]-1<<" 0 "<<lb.lbSize[1]-1<<" 0 "<<lb.lbSize[2]-1<<"\" "
            <<"Origin=\"0.0 0.0 0.0\" Spacing=\""<<lb.unit.Length<<" "<<lb.unit.Length<<" "<<lb.unit.Length<<"\">\n";
    paraviewFluidFile<<"  <Piece Extent=\"0 "<<lb.lbSize[0]-1<<" 0 "<<lb.lbSize[1]-1<<" 0 "<<lb.lbSize[2]-1<<"\">\n";
    paraviewFluidFile<<"   <PointData>\n";
    paraviewFluidFile<<"    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
    for (i=0; i<Pnumber; ++i){
        if (!lb.types[i].isInsideParticle()) {
            paraviewFluidFile<<lb.types[i].getType()<<" ";
        }
        else {
            paraviewFluidFile<<1<<" ";
        }
    }
    paraviewFluidFile<<"    </DataArray>\n";
//    for (int j=1; j<lbmDirec; ++j) {
//        paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"j"<<j<<"\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//        for (i=0; i<Pnumber; ++i){
//            if (lb.curves[i]==0) {
//                Zero.print(paraviewFluidFile);
//            }
//            else {
//                tVect deltaVec=lb.unit.Length*lb.curves[i]->delta[j]*v[j];
//                deltaVec.print(paraviewFluidFile);
//            }
//        }
//        paraviewFluidFile<<"    </DataArray>\n";
//    }
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (i=0; i<Pnumber; ++i){
        if (lb.nodes[i]==0) {
            Zero.print(paraviewFluidFile);
        }
        else {
            const tVect uPhysic=lb.nodes[i]->u*lb.unit.Speed;
            uPhysic.print(paraviewFluidFile);
        }
    }
    paraviewFluidFile<<"    </DataArray>\n";
//    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"f\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//    for (i=0; i<Pnumber; ++i){
//        if (lb.nodes[i]==0) {
//            Zero.print(paraviewFluidFile);
//        }
//        else {
//            tVect fPhysic=lb.nodes[i]->hydroForce*lb.unit.Force;
//            fPhysic.print(paraviewFluidFile);
//        }
//    }
//    paraviewFluidFile<<"    </DataArray>\n";
    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
    for (i=0; i<Pnumber; ++i){
        if (lb.nodes[i]==0) {
            paraviewFluidFile<<zero<<" ";
        }
        else {
            if (lb.nodes[i]->n!=0.0) {
                paraviewFluidFile<<0.3333333*(lb.nodes[i]->n-1.0)*lb.unit.Pressure<<" ";
            }
            else {
                paraviewFluidFile<<lb.nodes[i]->n<<" ";
            }
        }
    }
    paraviewFluidFile<<"    </DataArray>\n";
    if (lb.nonNewtonian) {
        paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"dynVisc\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (i=0; i<Pnumber; ++i){
            if (lb.nodes[i]==0) {
                paraviewFluidFile<<zero<<" ";
            }
            else {
                paraviewFluidFile<<lb.nodes[i]->visc*lb.unit.DynVisc<<" ";
            }
        }
        paraviewFluidFile<<"    </DataArray>\n";
//        paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"shearRate\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
//        for (i=0; i<Pnumber; ++i){
//            if (lb.nodes[i]==0) {
//                paraviewFluidFile<<zero<<" ";
//            }
//            else {
//                paraviewFluidFile<<lb.nodes[i]->shearRate*lb.unit.AngVel<<" ";
//            }
//        }
//        paraviewFluidFile<<"    </DataArray>\n";
    }
    if (lb.freeSurface) {
        paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"AAAmass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (i=0; i<Pnumber; ++i){
            if (lb.nodes[i]==0) {
                paraviewFluidFile<<zero<<" ";
            }
            else {
                paraviewFluidFile<<lb.nodes[i]->mass*lb.unit.Density<<" ";
            }
        }
        paraviewFluidFile<<"    </DataArray>\n";
    }
    if (demSolve) {
        paraviewFluidFile<<"    <DataArray type=\"Int16\" Name=\"solidIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
        for (i=0; i<Pnumber; ++i){
            paraviewFluidFile<<lb.types[i].getSolidIndex()<<" ";
        }
        paraviewFluidFile<<"    </DataArray>\n";
    }
    paraviewFluidFile<<"   </PointData>\n";
    paraviewFluidFile<<"   <CellData>\n";
    paraviewFluidFile<<"   </CellData>\n";
    paraviewFluidFile<<"  </Piece>\n";
    paraviewFluidFile<<" </ImageData>\n";
    paraviewFluidFile<<"</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

// print export stuff

void IO::exportMaxSpeedFluid(const LB& lb) {
        // fluid max velocity
        double maxFluidSpeed=0.0;
        for (int it=0; it<lb.totNodes; ++it) {
            if (lb.types[it].isActive()) {
                maxFluidSpeed=std::max(maxFluidSpeed,lb.nodes[it]->u.norm2());
            }
        }
        maxFluidSpeed=sqrt(maxFluidSpeed);
        cout<<"MaxFSpeed= "<<maxFluidSpeed*lb.unit.Speed<<"("<<int(lbmDt*lbmDt*maxFluidSpeed*100.0/0.01)<<"%)\t";
        exportFile<<"MaxFSpeed= "<<maxFluidSpeed*lb.unit.Speed<<"("<<int(lbmDt*lbmDt*maxFluidSpeed*100.0/0.01)<<"%)\t";
        
        // printing max speed
        maxSpeedFile.open(maxSpeedFileName.c_str(), ios::app);
        maxSpeedFile<<realTime<<" "<<maxFluidSpeed*lb.unit.Speed<<"\n";
        maxSpeedFile.close();
}

void IO::exportMaxSpeedParticles(const DEM& dem) {
    // calculating speed for real time check
    // particle max velocity
    double maxPartSpeed=0.0;
    for (int n=0;n<dem.elmts.size();++n){
        maxPartSpeed=std::max(maxPartSpeed,dem.elmts[n].x1.norm2());
    }
    maxPartSpeed=sqrt(maxPartSpeed);
    cout<<"MaxPSpeed= "<<maxPartSpeed<<"\t";
    exportFile<<"MaxPSpeed= "<<maxPartSpeed<<"\t";
    
    // particle max rotational velocity
    double maxPartSpin=0.0;
    for (int n=0;n<dem.elmts.size();++n){
        maxPartSpin=std::max(maxPartSpin,dem.elmts[n].wGlobal.norm2());
    }
    maxPartSpin=sqrt(maxPartSpin);
    cout<<"MaxPSpin= "<<maxPartSpin<<"\t";
    exportFile<<"MaxPSpin= "<<maxPartSpin<<"\t";

    // cout<<"Force= ";
    // dem.elmts[0].FHydro.show();
    cout<<"\t";
}

void IO::exportTotalMass(const LB& lb) {
     // total fluid mass
    double massTot=totFluidMass(lb);
    cout<<"Volume="<<massTot*lb.unit.Volume<<"; Mass = "<<massTot*lb.unit.Mass<<"\t";
    exportFile<<"Volume="<<massTot*lb.unit.Volume<<"; Mass = "<<massTot*lb.unit.Mass<<"\t";
}

void IO::exportPlasticity(const LB& lb) {
    // fluid plasticity state
    double percPlastic=totPlastic(lb);
    cout<<"Plastic ="<<int(percPlastic)<<"%;\t";
    exportFile<<"Plastic ="<<int(percPlastic)<<"%;\t";
    // printing plasticity level
    plasticityFile.open(plasticityFileName.c_str(), ios::app);
    plasticityFile<<realTime<<" "<<percPlastic<<"\n";
    plasticityFile.close();

}

void IO::exportShearCell(const LB& lb, const DEM& dem) {
    // apparent viscosity from shear cell
    double appVisc=0.0;
    double externalShear=0.0;
    double wallStress=0.0;
    apparentViscosity(lb, dem.walls, externalShear, wallStress, appVisc);
    cout<<"App Visc= "<<appVisc<<"\t";
    exportFile<<"App Visc= "<<appVisc<<"\t";

//            double volumeParticles=totVolumeParticles(nodes);
//            cout<<"Volume particles="<<volumeParticles*unit.Length*unit.Length*unit.Length<<"\t";
//        cout<<" Vel= ";
//        dem.elmts[0].x1.show();
//        cout<<" Hydro= ";
//        dem.elmts[0].FHydro.show();
//        cout<<" Wall= ";
//        dem.elmts[0].FWall.show();
//        dem.walls[0].force.show();
    tVect xDirec=tVect(1.0,0.0,0.0);
    cout<<"wallDown = "<<dem.walls[0].FHydro.dot(xDirec)<<" wallUp = "<<dem.walls[1].FHydro.dot(xDirec)<<"\t";
    cout<<"wallDown = "<<dem.walls[0].FParticle.dot(xDirec)<<" wallUp = "<<dem.walls[1].FParticle.dot(xDirec)<<"\t";
}

void IO::exportEnergy(const DEM& dem) {

    cout<<"eKin = "<<dem.particleEnergy.kin<<"\t";
    exportFile<<"eKin = "<<dem.particleEnergy.kin<<"\t";
    cout<<"eGrav = "<<dem.particleEnergy.grav<<"\t";
    exportFile<<"eGrav = "<<dem.particleEnergy.grav<<"\t";
    cout<<"eTot = "<<dem.particleEnergy.total<<"\t";
    exportFile<<"eTot = "<<dem.particleEnergy.total<<"\t";
}

void IO::exportForces(const DEM& dem) {
    // printing force
    // hydraulic and collision forces
    forceFile.open(forceFileName.c_str(), ios::app);
    double collisionForce=collisionForceTot(dem.elmts);
    double hydraulicForce=hydraulicForceTot(dem.elmts);
    forceFile<<realTime<<" "<<collisionForce<<" "<<hydraulicForce<<"\n";
    forceFile.close();
}

void IO::exportDrum(const LB& lb, const DEM& dem) {
    // printing centers of mass (fluid, particle and total)
    drumMassCenter(lb, dem.elmts, dem.sphereMat.density, dem.cylinders[0]);

    // printing phase positions
    drumPositions(lb, dem.elmts, dem.cylinders[0]);

    if (lbmSolve) {
        drumProfileProperties(lb, dem.elmts, dem.cylinders[0]);
    }

    if (demSolve) {
        drumParticleProperties(dem);
    }


}

void IO::exportForceObstacle(const DEM& dem) {
    // printing total force on obstacle
    obstacleFile.open(obstacleFileName.c_str(), ios::app);
    tVect totObstacleForce=totForceObject(dem.objects);
    obstacleFile<<realTime<<" ";
    totObstacleForce.printLine(obstacleFile);
    obstacleFile.close();
}

// data elaboration

double IO::totPlastic(const LB& lb) const {
    // prints the total mass in the free fluid domain
    static double maxVisc=(maxTau-0.5)/3/lbmDt;
    unsigned int totPlastic=0;
    unsigned int totActive=0;

    for (int it=0; it<lb.nodes.size(); ++it){
        if (lb.types[it].isActive()) {
            ++totActive;
            if (lb.nodes[it]->visc>0.95*maxVisc) {
                ++totPlastic;
            }
        }
    }
    double perc=100.0*double(totPlastic)/double(totActive);
    return perc;
}

double IO::totFluidMass(const LB& lb) const {
    // prints the total mass in the free fluid domain
    double mass=0.0;

    if (lbmSolve) {
        for (int it=0; it<lb.nodes.size(); ++it){
            if (lb.types[it].isActive() && !lb.types[it].isInsideParticle()) {
                mass+=lb.nodes[it]->mass;
            }
        }
    }
    return mass;
}

double IO::totParticleMass(const elmtList& elmts) const {
    // prints the total mass in the free fluid domain
    double mass=0.0;

    for (int n=0; n<elmts.size(); ++n){
        // calculate mass
        mass+=elmts[n].m;
    }

    return mass;
}

void  IO::apparentViscosity(const LB& lb, const wallList& walls, double& externalShear, double& wallStress, double& appVisc) const {
    // computes the apparent viscosity of a shear cell, as the ratio between shear at the wall and imposed shear rate
    // the cell is sheared in direction x and has moving wall in z
    appVisc=0.0;
    tVect xDirec=tVect(1.0,0.0,0.0);

    double cellSize=double(lb.lbSize[2]-2)*lb.unit.Length;
    
    externalShear=0.5*(walls[1].vel.dot(xDirec)-walls[0].vel.dot(xDirec))/cellSize;
    
    double plateSize=double(lb.lbSize[0]-2)*double(lb.lbSize[1]-2)*lb.unit.Length*lb.unit.Length;

    wallStress=-0.5*(walls[1].FHydro.dot(xDirec)+walls[1].FParticle.dot(xDirec)-walls[0].FHydro.dot(xDirec)-walls[0].FParticle.dot(xDirec))/plateSize;

    appVisc=0.5*wallStress/externalShear;
}

tVect IO::centerOfMassFluid(const LB& lb) const {
    // prints the total mass in the free fluid domain
    tVect center(0.0,0.0,0.0);
    double totMass(0.0);

    if (lbmSolve) {
        for (int it=0; it<lb.nodes.size(); ++it) {
            if (lb.types[it].isActive() && !lb.types[it].isInsideParticle()) {
                center+=lb.nodes[it]->mass*lb.positions[it];
                totMass+=lb.nodes[it]->mass;
            }
        }
        return center/totMass;
    }
    else {
        return Zero;
    }

    
}

tVect IO::centerOfMassParticles(const elmtList& elmts) const {
    // prints the total mass in the free fluid domain
    tVect center(0.0,0.0,0.0);
    double totMass(0.0);

    for (int p=0; p<elmts.size(); ++p){
        center+=elmts[p].m*elmts[p].x0;
        totMass+=elmts[p].m;
    }

    return center/totMass;
}

double IO::hydraulicForceTot(const elmtList& elmts) const {
    double totForce=0.0;

    for (int n=0; n<elmts.size(); ++n) {
        totForce+=elmts[n].FHydro.norm();
    }

    return totForce;
}

double IO::collisionForceTot(const elmtList& elmts) const {
    double totForce=0.0;

    for (int n=0; n<elmts.size(); ++n){
        totForce+=elmts[n].FParticle.norm();
        totForce+=elmts[n].FWall.norm();
    }

    return totForce;
}

double IO::totVolumeParticles(const LB& lb) const {
    // prints the total mass of the particles
    double volume=0.0;

        for (int it=0; it<lb.nodes.size(); ++it){
            if (lb.types[it].isInsideParticle()) {
                volume+=1;
            }
        }
    return volume;
}

tVect IO::totForceObject(const objectList& objects) const {
    // prints the total mass in the free fluid domain
    tVect totForce(0.0,0.0,0.0);

    for (int o=0; o<objects.size(); ++o){
        totForce+=objects[o].FParticle;
    }

    return totForce;
}

// DRUM

void IO::drumMassCenter(const LB& lb, const elmtList& elmts, const double& particleDensity, const cylinder& drumCylinder) {

    // FLUID
    // calculate fluid mass
    const double fluidMass=totFluidMass(lb)*lb.unit.Mass;
    // calculate fluid center of mass
    const tVect fluidCenter=centerOfMassFluid(lb)*lb.unit.Length;
    
    // calculate fluid volume
    const double fluidVolume=fluidMass/lb.unit.Mass*lb.unit.Volume;
    // calculate fluid center of mass angle
    const double fluidAngle=drumCylinder.getTheta(fluidCenter, Zm, Xdirec);

    // PARTICLES
    // calculate particle center of mass
    const tVect particleCenter=centerOfMassParticles(elmts);
    // calculate particle mass
    const double particleMass=totParticleMass(elmts);
    // calculate particle volume
    const double particleVolume=particleMass/particleDensity;
    // calculate particle center of mass angle
    const double particleAngle=drumCylinder.getTheta(particleCenter, Zm, Xdirec);

    // MIXTURE
    // calculate mixture mass
    const double mixtureMass=fluidMass+particleMass;
    // calculate mixture center of mass
    const tVect mixtureMassCenter=(fluidMass*fluidCenter+particleMass*particleCenter)/mixtureMass;
    // calculate mixture center of mass angle
    const double mixtureMassAngle=drumCylinder.getTheta(mixtureMassCenter, Zm, Xdirec);
    // calculate mixture center of volume
    const double mixtureVolume=fluidVolume+particleVolume;
    // calculate mixture center of volume
    const tVect mixtureVolumeCenter=(fluidVolume*fluidCenter+particleVolume*particleCenter)/mixtureVolume;
    // calculate mixture center of volume angle
    const double mixtureVolumeAngle=drumCylinder.getTheta(mixtureVolumeCenter, Zm, Xdirec);

    // print all the stuff
    centerOfMassFile.open(centerOfMassFileName.c_str(), ios::app);

    centerOfMassFile<<realTime<<" ";
    // fluid
    centerOfMassFile<<fluidMass<<" ";
    fluidCenter.print(centerOfMassFile);
    centerOfMassFile<<fluidAngle<<" ";
    // particles
    centerOfMassFile<<particleMass<<" ";
    particleCenter.print(centerOfMassFile);
    centerOfMassFile<<particleAngle<<" ";
    // mixture mass
    centerOfMassFile<<mixtureMass<<" ";
    mixtureMassCenter.print(centerOfMassFile);
    centerOfMassFile<<mixtureMassAngle<<" ";
    // mixture volume
    centerOfMassFile<<mixtureVolume<<" ";
    mixtureVolumeCenter.print(centerOfMassFile);
    centerOfMassFile<<mixtureVolumeAngle<<"\n";

    centerOfMassFile.close();
}

void IO::drumPositions(const LB& lb, const elmtList& elmts, const cylinder& drumCylinder) {
    // finds position of particle and fluid front and rear in the drum

    // PARTICLES
    static const double minAngle=-1.0;
    static const double maxAngle=1.0;
    static const double deltaAngle=0.01;
    static const unsigned int totIndices=floor((maxAngle-minAngle)/deltaAngle);

/*/    // PARTICLE DENSITY
//    if (elmts.size()) {
//
//        unsigned int particleDensity[totIndices];
//        for (int i=0; i<totIndices; ++i) {
//            particleDensity[i]=0;
//        }
//
//        for (int n=0; n<elmts.size(); ++n) {
//            const double particleAngle=drumCylinder.getTheta(elmts[n].x0, Zm, X);
//            const unsigned int particleAngleIndex=floor((particleAngle-minAngle)/deltaAngle);
//            ++particleDensity[particleAngleIndex];
//        }

//        // print all the stuff
//        particleDensityFile.open(particleDensityFileName.c_str(), ios::app);
//
//        particleDensityFile<<realTime<<" "; //particlePositionFile<<"time minAngle maxAngle\n";
//        // particles
//        for (int i=0; i<totIndices; ++i) {
//            particleDensityFile<<particleDensity[i]<<" ";
//        }
//        particleDensityFile<<"\n";
//        particleDensityFile.close();

    }

/    // PARTICLES (old)
//    if (elmts.size()) {
//        double particleFront(1.0e20), particleRear(-1.0e20);
//        unsigned int particleFrontIndex(0), particleRearIndex(0);
//        for (int n=0; n<elmts.size(); ++n) {
//            const double x=elmts[n].x0.dot(X);
//            if (x<particleFront) {
//                particleFront=x;
//                particleFrontIndex=n;
//            }
//            else if (x>particleRear) {
//                particleRear=x;
//                particleRearIndex=n;
//            }
//        }
//
//        const double particleFrontAngle=drumCylinder.getTheta(elmts[particleFrontIndex].x0, Zm, X);
//        const double particleRearAngle=drumCylinder.getTheta(elmts[particleRearIndex].x0, Zm, X);
//
//        // print all the stuff
//        particlePositionFile.open(particlePositionFileName.c_str(), ios::app);
//
//        particlePositionFile<<realTime<<" "; //particlePositionFile<<"time minAngle maxAngle\n";
//        // particles
//        particlePositionFile<<particleFrontAngle<<" ";
//        particlePositionFile<<particleRearAngle<<"\n";
//
//        particlePositionFile.close();
//
//    } */

    // FRONTS
    if (lbmSolve) {
        static const unsigned int sections=10;
        // size of the drum
        static const unsigned int drumWidth=lb.lbSize[1];
        // spacing between sections
        static const unsigned int spacing=int(double(drumWidth)/3.0/double(sections));
        // definition of sections
        static unsigned int rowNumber[sections];
        // first is at 1/3 of the drum
        rowNumber[0]=int(double(drumWidth)/3.0);
        for (int i=1; i<sections-1; ++i) {
            rowNumber[i]=rowNumber[0]+spacing*i;
        }

        unsigned int fluidRear[sections], fluidFront[sections], fluidRearIndex[sections], fluidFrontIndex[sections]; // <- Must be initialized!!!
        unsigned int particleRear[sections], particleFront[sections], particleRearIndex[sections], particleFrontIndex[sections]; // <- Must be initialized!!!

        // initialize
        for (int i=0; i<sections; ++i) {
            fluidRear[i]=0;
            fluidFront[i]=1000000;
            fluidRearIndex[i]=0;
            fluidFrontIndex[i]=0;
            //
            particleRear[i]=0;
            particleFront[i]=1000000;
            particleRearIndex[i]=0;
            particleFrontIndex[i]=0;
        }

        // we go through all active nodes, and check nodes belonging to the sections
        for (int it=0; it<lb.totNodes; ++it) {
            if (lb.types[it].isActive()) {
                const unsigned int y=lb.getY(it);
                for (int i=0; i<sections; ++i) {
                    if (y==rowNumber[i]) {
                        const unsigned int x=lb.getX(it);
                        if (x<fluidFront[i]) {
                            fluidFront[i]=x;
                            fluidFrontIndex[i]=it;
                        }
                        else if (x>fluidRear[i]) {
                            fluidRear[i]=x;
                            fluidRearIndex[i]=it;
                        }
                    }
                }
            }
            if (lb.types[it].isInsideParticle()) {
                const unsigned int y=lb.getY(it);
                for (int i=0; i<sections; ++i) {
                    if (y==rowNumber[i]) {
                        const unsigned int x=lb.getX(it);
                        if (x<particleFront[i]) {
                            particleFront[i]=x;
                            particleFrontIndex[i]=it;
                        }
                        else if (x>particleRear[i]) {
                            particleRear[i]=x;
                            particleRearIndex[i]=it;
                        }
                    }
                }
            }
        }

        // cylindrical coordinates
        double fluidFrontAngle[sections], fluidRearAngle[sections];
        double particleFrontAngle[sections], particleRearAngle[sections];
        for (int i=0; i<sections; ++i) {
            const tVect fluidFrontPosition=lb.positions[fluidFrontIndex[i]]*lb.unit.Length;
            const tVect fluidRearPosition=lb.positions[fluidRearIndex[i]]*lb.unit.Length;
            fluidFrontAngle[i]=drumCylinder.getTheta(fluidFrontPosition, Zm, Xdirec);
            fluidRearAngle[i]=drumCylinder.getTheta(fluidRearPosition, Zm, Xdirec);
            //
            const tVect particleFrontPosition=lb.positions[particleFrontIndex[i]]*lb.unit.Length;
            const tVect particleRearPosition=lb.positions[particleRearIndex[i]]*lb.unit.Length;
            particleFrontAngle[i]=drumCylinder.getTheta(particleFrontPosition, Zm, Xdirec);
            particleRearAngle[i]=drumCylinder.getTheta(particleRearPosition, Zm, Xdirec);
        }

        // computing averages
        double averageFluidFront(0.0), averageParticleFront(0.0), averageFluidRear(0.0), averageParticleRear(0.0);
        for (int i=0; i<sections; ++i) {
            averageFluidFront+=fluidFrontAngle[i]/double(sections);
            averageParticleFront+=particleFrontAngle[i]/double(sections);
            averageFluidRear+=fluidRearAngle[i]/double(sections);
            averageParticleRear+=particleRearAngle[i]/double(sections);
        }

        //  writing in position file
        positionFile.open(positionFileName.c_str(), ios::app);
        positionFile<<realTime<<" "<<averageFluidFront<<" "<<averageFluidRear<<" "<<averageParticleFront<<" "<<averageParticleRear<<"\n";
        positionFile.close();

    }

}

//void IO::drumSurfaceProperties(const LB& lb, const elmtList& elmts, const cylinder& drum) {
//
//    string title;
//    const char* charSurfaceFile;
//    ofstream outputSurfaceFile;
//    string surfaceFile;
//    title=workDirectory+"/drumSurface/";
//    stringstream sf;
//    sf<<currentTimeStep;
//    surfaceFile=title+sf.str()+".dat";
//    charSurfaceFile=surfaceFile.c_str();
//    outputSurfaceFile.open(charSurfaceFile);
//
//    // list for raw data
//    doubleList fluidHeigth, solidHeigth, particleHeigth;
//    vecList fluidVelocity, particleVelocity;
//    // list for post-processed data
//    doubleList angle, fluidProjectedHeigth, particleProjectedHeigth, totalProjectedHeigth, fluidProjectedVelocity, particleProjectedVelocity;
//
//    static unsigned int yCenter=int(double(lb.lbSize[1])/2.0);
//    static unsigned int xSize=lb.lbSize[0];
//
//    // initialize containers
//    fluidHeigth.resize(xSize);
//    solidHeigth.resize(xSize);
//    particleHeigth.resize(xSize);
//    angle.resize(xSize);
//    fluidProjectedHeigth.resize(xSize);
//    particleProjectedHeigth.resize(xSize);
//    totalProjectedHeigth.resize(xSize);
//    fluidVelocity.resize(xSize);
//    particleVelocity.resize(xSize);
//    fluidProjectedVelocity.resize(xSize);
//    particleProjectedVelocity.resize(xSize);
//
//    // getting raw data (Cartesian coordinate)
//    for (int it=0; it<lb.totNodes; ++it){
//        if (lb.getY(it)==yCenter) {
//            // height of the particles
//            if (lb.types[it].isInsideParticle()) {
//                particleHeigth[lb.getX(it)]+=1.0*lb.unit.Length;
//            }
//            // height of the fluid
//            else if (lb.types[it].isActive()) {
//                fluidHeigth[lb.getX(it)]+=lb.nodes[it]->liquidFraction()*lb.unit.Length;
//            }
//            // height of the drum
//            else if (lb.types[it].isCurvedWall()) {
//                solidHeigth[lb.getX(it)]+=1.0*lb.unit.Length;
//            }
//
//            // fluid velocity at the surface
//            if (lb.types[it].isActive()) {
//                // checking if we are indeed at the surface
//                const unsigned int upLink=lb.neighbors[it].d[5];
//                if (lb.types[upLink].isGas()) {
//                    // getting velocity
//                    fluidVelocity[lb.getX(it)]=lb.nodes[it]->u*lb.unit.Speed;
//                }
//            }
//            // particle velocity at the surface (if it exists)
//            if (lb.types[it].isInsideParticle()) {
//                // checking if we are indeed at the surface
//                const unsigned int upLink=lb.neighbors[it].d[5];
//                if (lb.types[upLink].isGas()) {
//                    // particle index
//                    const unsigned int solidIndex=lb.types[it].getSolidIndex();
//                    // getting velocity
//                    particleVelocity[lb.getX(it)]=elmts[solidIndex].x1;
//                }
//            }
//        }
//    }
//
//    // elaborating data (passing to cylindrical coordinate)
//    for (int is=1; is<xSize; ++is) {
//        // coordinate Z of the fluid surface in Cartesian coordinates
//        const double totalHeight=solidHeigth[is]+fluidHeigth[is]+particleHeigth[is];
//        // vectorial position of the fluid surface in Cartesian coordinates
//        const tVect heightPosition=tVect(double(is)*lb.unit.Length, double(yCenter)*lb.unit.Length, totalHeight);
//        // height of the fluid surface in cylindrical coordinates
//        const double angleHeight=drum.dist(heightPosition);
//        // coordinate theta of the fluid surface in cylindrical coordinates
//        angle[is]=drum.getTheta(heightPosition, Zm, X);
//        // we need a trick to transform the particle and fluid height in the new system
//        // we suppose a change in the height proportional to the difference between
//        // the surface level in cylindrical and Cartesian coordinates
//        const double prop=angleHeight/totalHeight;
//        fluidProjectedHeigth[is]=fluidHeigth[is]*prop;
//        particleProjectedHeigth[is]=particleHeigth[is]*prop;
//        totalProjectedHeigth[is]=fluidProjectedHeigth[is]+particleProjectedHeigth[is];
//
//        // getting versor pointing towards cylinder axes
//        const tVect radiusVersor=drum.centerVersor(heightPosition);
//        // getting versor on theta direction
//        const tVect thetaVersor=radiusVersor.cross(Ym);
//        // get component of velocity on theta direction
//        fluidProjectedVelocity[is]=fluidVelocity[is].dot(thetaVersor);
//        particleProjectedVelocity[is]=particleVelocity[is].dot(thetaVersor);
//
//
//
//        outputSurfaceFile<<angle[is]<<" "<<fluidProjectedHeigth[is]<<" "<<particleProjectedHeigth[is]<<" "<<totalProjectedHeigth[is]<<" "<<fluidProjectedVelocity[is]<<" "<<particleProjectedVelocity[is]<<"\n";
//
//    }
//    outputSurfaceFile.close();
////    outputSolidShapeFile.close();
//}
//
//void IO::drumBaseProperties(const LB& lb, const elmtList& elmts, const cylinder& drum) {
//    // prints the total mass in the free fluid domain
//
//    string title;
//    const char* charBaseFile;
//    ofstream outputBaseFile;
//    string baseFile;
//    title=workDirectory+"/drumBase/";
//    stringstream sf;
//    sf<<currentTimeStep;
//    baseFile=title+sf.str()+".dat";
//    charBaseFile=baseFile.c_str();
//    outputBaseFile.open(charBaseFile);
//
//    // list for raw data
//    doubleList porePressure, impactPressure, totalPressure, solidHeight;
//    // list for post-processed data
//    doubleList angle;
//
//    int range=6;
//    static unsigned int yc=int(double(lb.lbSize[1])/2.0);
//    static unsigned int pressureSize=lb.lbSize[0];
//
//    porePressure.resize(pressureSize);
//    impactPressure.resize(pressureSize);
//    totalPressure.resize(pressureSize);
//    solidHeight.resize(pressureSize);
//    //
//    angle.resize(pressureSize);
//
//    for (int is=1; is<pressureSize; ++is) {
//        porePressure[is]=0.0;
//        impactPressure[is]=0.0;
//        totalPressure[is]=0.0;
//    }
//
//    unsIntList particleCount;
//    particleCount.resize(elmts.size());
//
//    unsIntList solidCount, fluidCount, solidIndices;
//    solidCount.resize(pressureSize);
//    solidIndices.resize(pressureSize);
//    fluidCount.resize(pressureSize);
//
//    // getting raw data (Cartesian coordinate)
//    for (int it=0; it<lb.totNodes; ++it) {
////        if ( (nodes[it].y>=(yc-range)) && (nodes[it].y<=(yc+range)) ) {
////        if (nodes[it].y==yc) {
//            // link to lower node
//            unsigned int low=lb.neighbors[it].d[6];
//            // if the lower node is a wall we assume we have the drum
//            if (lb.types[low].isWall()) {
//                solidHeight[lb.getX(it)]=double(lb.getZ(low))*lb.unit.Length;
//                // if we have fluid we can record the pore pressure
//                if (lb.types[it].isActive()) {
//                    porePressure[lb.getX(it)]+=0.3333333*(lb.nodes[it]->n-1.0)*lb.unit.Pressure;
//                    fluidCount[lb.getX(it)]++;
//                }
//                // if we have particle we can record the collision force
//                else if (lb.types[it].isInsideParticle()) {
//                    // getting solid index (we need to average the force on all solid nodes sharing the element)
//                    unsigned int solidIndex=lb.types[it].getSolidIndex();
//                    solidIndices[lb.getX(it)]=solidIndex;
//                    particleCount[solidIndex]++;
//                    solidCount[lb.getX(it)]++;
//                    // getting radius (we need to project on the drum surface
//                    tVect radius=drum.centerDist(lb.unit.Length*lb.positions[it]);
//                    tVect normRadius=-1.0*radius/radius.norm();
//                    // finally getting the pressure (still needs to be averaged)
//                    impactPressure[lb.getX(it)]+=(elmts[solidIndex].FWall.dot(normRadius))/(lb.unit.Length*lb.unit.Length);
//                }
//            }
////        }
//    }
//
//    // elaborating data (passing to cylindrical coordinate)
//    for (int is=1; is<pressureSize; ++is) {
//        // vectorial position of the drum surface in Cartesian coordinates
//        tVect heightPosition=tVect(double(is)*lb.unit.Length, double(yc)*lb.unit.Length, solidHeight[is]);
//        // coordinate theta of the drum surface in cylindrical coordinates
//        angle[is]=drum.getTheta(heightPosition, Zm, X);
//        // we have to average over
//        // 1. the number of cell in the y and x direction
//        // 2. the number of solid nodes (only for collisional)
//        if (impactPressure[is]!=0.0) {
//            impactPressure[is]/=double(particleCount[solidIndices[is]]);
//            impactPressure[is]/=double(solidCount[is]);
//        }
//        if (porePressure[is]!=0.0) {
//            porePressure[is]/=double(fluidCount[is]);
//        }
//        totalPressure[is]=porePressure[is]+impactPressure[is];
//        outputBaseFile<<angle[is]<<" "<<porePressure[is]<<" "<<impactPressure[is]<<" "<<totalPressure[is]<<"\n";
//
//    }
//    outputBaseFile.close();
////    outputSolidShapeFile.close();
//}

void IO::drumProfileProperties(const LB& lb, const elmtList& elmts, const cylinder& drum) {
    // list for raw data
    doubleList pressure, fluidHeight, solidHeight, particleHeight, particleSurface;
    // raw data with occurrences
    std::vector<vecList> fluidVelocity, particleVelocity;
    std::vector<doubleList> fluidSurface;
    // list for post-processed data
    doubleList angle, fluidProjectedVelocity, particleProjectedVelocity;
    // stuff to average velocity
//    std::vector<intList> particleVelocityNumber, fluidVelocityNumber, particleSurfaceNumber;


    int range=6;
    static unsigned int yCenter=int(double(lb.lbSize[1])/2.0);
    static unsigned int xSize=lb.lbSize[0];

    // number of sections
    static const unsigned int sections=10;
    // size of the drum
    static const unsigned int drumWidth=lb.lbSize[1];
    // spacing between sections
    static const unsigned int spacing=int(double(drumWidth)/3.0/double(sections));
    // definition of sections
    static unsigned int rowNumber[sections];
    // first is at 1/3 of the drum
    rowNumber[0]=int(double(drumWidth)/3.0);
    for (int i=1; i<sections; ++i) {
        rowNumber[i]=rowNumber[0]+spacing*i;
    }

    // initialize containers (raw)
    pressure.resize(xSize);
    fluidHeight.resize(xSize);
    solidHeight.resize(xSize);
    particleHeight.resize(xSize);
    // initialize containers (with occurrences)
    fluidVelocity.resize(xSize);
    particleVelocity.resize(xSize);
    fluidSurface.resize(xSize);
    particleSurface.resize(xSize);
    for (int i=0; i<xSize; ++i) {
        fluidVelocity[i].resize(sections);
        particleVelocity[i].resize(sections);
        fluidSurface[i].resize(sections);
    }
    // initialize containers (processed)
    angle.resize(xSize);
    fluidProjectedVelocity.resize(xSize);
    particleProjectedVelocity.resize(xSize);

    for (int i=0; i<xSize; ++i) {
        // measures
        fluidHeight[i]=0.0;
        solidHeight[i]=0.0;
        particleHeight[i]=0.0;
        angle[i]=0.0;
        fluidProjectedVelocity[i]=0.0;
        particleProjectedVelocity[i]=0.0;
        particleSurface[i]=0.0;
        // with occurrences
        for (int j=0; j<sections; ++j) {
            fluidVelocity[i][j]=Zero;
            particleVelocity[i][j]=Zero;
            fluidSurface[i][j]=0;  
        }
    }

    // getting raw data (Cartesian coordinate)
    for (int it=0; it<lb.totNodes; ++it) {
        const unsigned int y=lb.getY(it);
        const unsigned int x=lb.getX(it);
        for (int i=0; i<sections; ++i) {
            if (y==rowNumber[i]) {
                // link to lower and upper node
                const unsigned int lowLink=lb.neighbors[it].d[6];
                const unsigned int upLink=lb.neighbors[it].d[5];
                if (lb.types[lowLink].isCurvedWall()) {
                    solidHeight[x]=(double(lb.getZ(lowLink))+0.5)*lb.unit.Length;
                }
                if (lb.types[it].isActive()) {
                    // height of the fluid
                    if (!lb.types[it].isInsideParticle()) {
                        fluidHeight[x]+=lb.nodes[it]->liquidFraction()*lb.unit.Length;
                    }
                    // if the lower node is a wall we assume we have the drum
                    if (lb.types[lowLink].isCurvedWall()) {
                        pressure[x]+=0.3333333*(lb.nodes[it]->n-1.0)*lb.unit.Pressure;
                    }
                    // if the upper node is gas we assume we are at the fluid surface
                    if (lb.types[upLink].isGas()) {
                        // getting velocity
                        fluidVelocity[x][i]=lb.nodes[it]->u*lb.unit.Speed;
                        // getting surface position (I assume the drum base level is the solidHeight)
//                        cout<<"A "<<fluidSurface[x]<<" "<<(double(lb.getZ(it))+lb.nodes[it]->liquidFraction())*lb.unit.Length<<endl;
                        fluidSurface[x][i]=(double(lb.getZ(it))+lb.nodes[it]->liquidFraction())*lb.unit.Length;
                    }
                }
                if (lb.types[it].isInsideParticle()) {
                    // height of the particles
                    particleHeight[x]+=1.0*lb.unit.Length;
                    // if the upper node is not inside the particle we check if we have a surface
                    if (!lb.types[upLink].isInsideParticle()) {
                        // getting surface position
                        const double surfacePosition=(double(lb.getZ(it))+0.5)*lb.unit.Length;
                        if (surfacePosition>particleSurface[x]) {
                            // particle index
                            const unsigned int solidIndex=lb.types[it].getSolidIndex();
                            // getting velocity
                            particleVelocity[x][i]=elmts[solidIndex].x1;
                            particleSurface[x]=surfacePosition;
                        }
                    }
                }
            }
        }
    }

    // declare and intialize average for measures with occurrences
    vecList avFluidVelocity, avParticleVelocity;
    doubleList avFluidSurface;

    avFluidVelocity.resize(xSize);
    avParticleVelocity.resize(xSize);
    avFluidSurface.resize(xSize);

    // averaging through the number of lines (and occurrences if necessary)
    for (int i=1; i<xSize; ++i) {
        // simple averages
        fluidHeight[i]/=double(sections);
        particleHeight[i]/=double(sections);
        pressure[i]/=double(sections);

        // average for velocities
        avFluidVelocity[i]=Zero;
        avParticleVelocity[i]=Zero;
        unsigned int fluidVelocityOccurrences(0), partVelocityOccurrences(0);
        for (int j=0; j<sections; ++j) {
            if (fluidVelocity[i][j].norm2()!=0) {
                ++fluidVelocityOccurrences;
                avFluidVelocity[i]+=fluidVelocity[i][j];
            }
            if (particleVelocity[i][j].norm2()!=0) {
                ++partVelocityOccurrences;
                avParticleVelocity[i]+=particleVelocity[i][j];
            }
        }
        avFluidVelocity[i]/=double(fluidVelocityOccurrences);
        avParticleVelocity[i]/=double(partVelocityOccurrences);
        
        // average for surface
        avFluidSurface[i]=0.0;
        for (int j=0; j<sections; ++j) {
            if (fluidSurface[i][j]!=0.0) {
                avFluidSurface[i]+=fluidSurface[i][j]-solidHeight[i];
            }
        }
        avFluidSurface[i]/=double(sections);
    }

    // elaborating data (passing to cylindrical coordinate)
    for (int i=1; i<xSize; ++i) {
        // vectorial position of the drum surface in Cartesian coordinates
        tVect heightPosition=tVect(double(i)*lb.unit.Length, double(yCenter)*lb.unit.Length, solidHeight[i]);
        // coordinate theta of the drum surface in cylindrical coordinates
        angle[i]=drum.getTheta(heightPosition, Zm, Xdirec);

        // getting versor pointing towards cylinder axes
        const tVect radiusVersor=drum.centerVersor(heightPosition);
        // getting versor on theta direction
        const tVect thetaVersor=radiusVersor.cross(Ym);
        // get component of velocity on theta direction
        fluidProjectedVelocity[i]=avFluidVelocity[i].dot(thetaVersor);
        particleProjectedVelocity[i]=avParticleVelocity[i].dot(thetaVersor);
    }

    // PRINTING DATA ON FILES

    // if first time, then just printout the angle file
    if (currentTimeStep==0) {
        // angle profile
        profileAngleFile.open(profileAngleFileName.c_str(), ios::app);
        for (int i=1; i<xSize; ++i) {
            profileAngleFile<<angle[i]<<" ";
//            cout<<i<<" "<<solidHeight[i]<<" "<<angle[i]<<endl;
        }
        profileAngleFile<<"\n";
        profileAngleFile.close();
    }
    else {
        // fluid level profile
        fluidLevelFile.open(fluidLevelFileName.c_str(), ios::app);
        fluidLevelFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            fluidLevelFile<<fluidHeight[i]<<" ";
        }
        fluidLevelFile<<"\n";
        fluidLevelFile.close();

        // particle level profile
        particleLevelFile.open(particleLevelFileName.c_str(), ios::app);
        particleLevelFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            particleLevelFile<<particleHeight[i]<<" ";
        }
        particleLevelFile<<"\n";
        particleLevelFile.close();

        // base pressure profile
        pressureFile.open(pressureFileName.c_str(), ios::app);
        pressureFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            pressureFile<<pressure[i]<<" ";
        }
        pressureFile<<"\n";
        pressureFile.close();

        // fluid surface speed profile
        fluidSpeedFile.open(fluidSpeedFileName.c_str(), ios::app);
        fluidSpeedFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            fluidSpeedFile<<fluidProjectedVelocity[i]<<" ";
        }
        fluidSpeedFile<<"\n";
        fluidSpeedFile.close();

        // fluid surface speed profile
        particleSpeedFile.open(particleSpeedFileName.c_str(), ios::app);
        particleSpeedFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            particleSpeedFile<<particleProjectedVelocity[i]<<" ";
        }
        particleSpeedFile<<"\n";
        particleSpeedFile.close();

        // particle surface profile
        particleSurfaceFile.open(particleSurfaceFileName.c_str(), ios::app);
        particleSurfaceFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            if (particleSurface[i]!=0) {
                particleSurfaceFile<<particleSurface[i]-solidHeight[i]<<" ";
            }
            else {
                particleSurfaceFile<<0.0<<" ";
            }
        }
        particleSurfaceFile<<"\n";
        particleSurfaceFile.close();

        // particle surface profile
        fluidSurfaceFile.open(fluidSurfaceFileName.c_str(), ios::app);
        fluidSurfaceFile<<realTime<<" ";
        for (int i=1; i<xSize; ++i) {
            fluidSurfaceFile<<avFluidSurface[i]<<" ";
        }
        fluidSurfaceFile<<"\n";
        fluidSurfaceFile.close();
    }

}

void IO::drumParticleProperties(const DEM& dem) {

    static const double minAngle=-1.0;
    static const double maxAngle=1.0;
    static const double deltaAngle=0.01;
    static const unsigned int totIndices=int((maxAngle-minAngle)/deltaAngle);

    // PRINTING DATA ON FILES

    // if first time, then just printout the angle file
    if (currentTimeStep==0) {
        // angle profile
        particleAngleFile.open(particleAngleFileName.c_str(), ios::app);
        for (int i=0; i<totIndices; ++i) {
            particleAngleFile<<minAngle+deltaAngle*double(i)<<" ";
        }
        particleAngleFile<<"\n";
        particleAngleFile.close();
    }
    else {

        // particle number
        particleNumberFile.open(particleNumberFileName.c_str(), ios::app);
        particleNumberFile<<realTime<<" ";
        for (int i=0; i<totIndices; ++i) {
            particleNumberFile<<double(particleNumberAngle[i])/double(timeStepCount)<<" ";
        }
        particleNumberFile<<"\n";
        particleNumberFile.close();

        // particle collision number
        collisionNumberFile.open(collisionNumberFileName.c_str(), ios::app);
        collisionNumberFile<<realTime<<" ";
        for (int i=0; i<totIndices; ++i) {
            collisionNumberFile<<double(collisionNumberAngle[i])/double(timeStepCount)<<" ";
        }
        collisionNumberFile<<"\n";
        collisionNumberFile.close();

        // particle hydraulic number
        hydraulicNumberFile.open(hydraulicNumberFileName.c_str(), ios::app);
        hydraulicNumberFile<<realTime<<" ";
        for (int i=0; i<totIndices; ++i) {
            hydraulicNumberFile<<double(hydraulicNumberAngle[i])/double(timeStepCount)<<" ";
        }
        hydraulicNumberFile<<"\n";
        hydraulicNumberFile.close();

        // particle collision force
        collisionForceFile.open(collisionForceFileName.c_str(), ios::app);
        collisionForceFile<<realTime<<" ";
        for (int i=0; i<totIndices; ++i) {
            collisionForceFile<<double(collisionForceAngle[i])/double(timeStepCount)<<" ";
        }
        collisionForceFile<<"\n";
        collisionForceFile.close();

        // particle hydraulic force
        hydraulicForceFile.open(hydraulicForceFileName.c_str(), ios::app);
        hydraulicForceFile<<realTime<<" ";
        for (int i=0; i<totIndices; ++i) {
            hydraulicForceFile<<double(hydraulicForceAngle[i])/double(timeStepCount)<<" ";
        }
        hydraulicForceFile<<"\n";
        hydraulicForceFile.close();

        // resetting the vectors
        thetaList.resize(totIndices);
        for (int i=0; i<totIndices; i++) {
            thetaList[i]=minAngle+double(i)*deltaAngle;
        }

        timeStepCount=0;
    }

}

void IO::updateDistributions(const DEM& dem, const cylinder& drum){
    // be sure that the first time this is executed, an instance of  <- FIX THIS!!!!!!!!


    static const double minAngle=-1.0;
    static const double maxAngle=1.0;
    static const double deltaAngle=0.01;
    static const unsigned int totIndices=floor((maxAngle-minAngle)/deltaAngle);

    if (timeStepCount==0 || currentTimeStep==0) {

        collisionForceAngle.clear();
        hydraulicForceAngle.clear();
        collisionNumberAngle.clear();
        hydraulicNumberAngle.clear();
        particleNumberAngle.clear();

        collisionForceAngle.resize(totIndices);
        hydraulicForceAngle.resize(totIndices);
        collisionNumberAngle.resize(totIndices);
        hydraulicNumberAngle.resize(totIndices);
        particleNumberAngle.resize(totIndices);

        for (int i=0; i<totIndices; ++i) {
            collisionForceAngle[i]=0.0;
            hydraulicForceAngle[i]=0.0;
            collisionNumberAngle[i]=0;
            hydraulicNumberAngle[i]=0;
            particleNumberAngle[i]=0;
        }

    }

    for (int n=0; n<dem.elmts.size(); ++n) {
        const double thisAngle=drum.getTheta(dem.elmts[n].x0, Zm, Xdirec);
        const int thetaIndex=int((thisAngle-minAngle)/deltaAngle);
        const double hydroForce=dem.elmts[n].FHydro.norm2();
        const double collisionForce=dem.elmts[n].FParticle.norm2();

        ++particleNumberAngle[thetaIndex];
        if (collisionForce!=0.0) {
            collisionForceAngle[thetaIndex]+=sqrt(collisionForce);
            ++collisionNumberAngle[thetaIndex];
        }
        if (hydroForce!=0.0) {
            hydraulicForceAngle[thetaIndex]+=sqrt(hydroForce);
            ++hydraulicNumberAngle[thetaIndex];
        }
    }

    ++timeStepCount;

}

// AVALANCHE


