
#include "node.h"

// basic functions

//tVect node::position(){
//    return tVect(double(x), double(y), double(z));
//}

// liquid fraction

double node::liquidFraction() const {
    return mass/n;
}

// functions working on distribution functions

void node::showPopulations() const {
    cout<<"(";
    for (int j=0; j<lbmDirec; ++j) {
        cout<<j<<":"<<f[j]<<", ";
    }
    cout<<")";
}

void node::initialize(const double& density,const tVect& velocity, const double& massFunction, const double& viscosity, const tVect& lbF){

    // just for initializing
    shearRate=0.0;
    //macroscopic properties of the initial condition
    n=density;
    mass=massFunction;
    u=lbmDt*lbF/2.0/n+velocity; //FORCE PART
    visc=viscosity;
//    force=tVect(0.0,0.0,0.0);

    setEquilibrium(lbF);
}

void node::setEquilibrium(const tVect& F){
    static double C1=3.0*lbmDt*lbmDt;
    static double C2=4.5*lbmDt*lbmDt*lbmDt*lbmDt;
    static double C3=1.5*lbmDt*lbmDt;
//    static double F1=3.0*dt*dt;
//    static double F2=9.0*dt*dt*dt*dt;

//    double tauF=3.0*visc*dt; // = tau-0.5

    const double usq=u.norm2();
//    tVect vmu, forcePart;
    for (int j=0; j<lbmDirec; ++j){
      // the following lines initialize f to be the local equilibrium values
        const double vu=u.dot(v[j]);
//        vmu=v[j]; // FORCE PART
//        vmu-=u; // FORCE PART
//        forcePart=F1*vmu+F2*vu*v[j]; //FORCE PART
//        dummyNode.f[j]=dummyNode.fs[j]=coeff[j]*dummyNode.n*(1.0+3.0*vu+4.5*vu*vu-1.5*usq)-2.0*(tau-0.5)*coeff[j]*forcePart.dot(F); //FORCE PART ((tau-0.5)=omegaf/omega)
//        f[j]=fs[j]=coeff[j]*n*(1.0+C1*vu+C2*vu*vu-C3*usq)-coeff[j]*forcePart.dot(F); //FORCE PART ((tau-0.5)=omegaf/omega)
        f[j]=fs[j]=coeff[j]*n*(1.0+C1*vu+C2*vu*vu-C3*usq);//+tauF*coeff[j]*forcePart.dot(F); //FORCE PART ((tau-0.5)=omegaf/omega)
    }
}

void node::reconstruct(){
    // reconstruction of macroscopical physical variables

    // density
    n=0.0;
    for (int j=0; j<lbmDirec; ++j) {
        n+=f[j];
    }

    // momentum
    u.reset();
    for (int j=0; j<lbmDirec; ++j) {
        u+=f[j]*v[j];
    }
//    if (u.norm2()>0.004) {
//        u/=1.25;
//    }

    // velocity
    u/=n;
}

void node::shiftVelocity(const tVect& F){
    const tVect totalForce=F+hydroForce;
    u+=0.5*lbmDt*totalForce/n;
}

void node::computeEquilibrium(double feq[]){
    static double C1=3.0*lbmDt*lbmDt;
    static double C2=4.5*lbmDt*lbmDt*lbmDt*lbmDt;
    static double C3=1.5*lbmDt*lbmDt;

    const double usq=u.norm2();

    for (int j=0; j<lbmDirec; ++j) {
        const double vu=u.dot(v[j]);
        feq[j]=coeff[j]*n*(1.0+C1*vu+C2*vu*vu-C3*usq);
    }
}

void node::computeShearRate(const double feq[], const double& turbConst, const double& plasticVisc, const double& yieldStress, const bool& nonNewtonian, const bool& turbulenceOn){

    // minimum and maximum viscosity
    static double minVisc=(minTau-0.5)/3/lbmDt;
    static double maxVisc=(maxTau-0.5)/3/lbmDt;
    const double tau=0.5+3.0*visc*lbmDt;

    tMat gamma(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

    for (int j=0; j<lbmDirec; ++j) {
        gamma+=(f[j]-feq[j])*vv[j];
    }
    gamma*=1.5*lbmDt/(tau*n);

    // shear rate (second invariant)
    shearRate=gamma.magnitude();

    // Smagorinsky turbulence
    double nuTurb=0.0;
    if (turbulenceOn) {
        nuTurb=turbConst*shearRate;
    }

    // Bingham model
    double nuApp=0.0;
    if (nonNewtonian) {
        nuApp=plasticVisc+yieldStress/(n*2.0*shearRate);
    }
    else {
        nuApp=visc;
    }
    // fixing plastic behavior at surface
//    if (isInterface()) {
//        visc=maxVisc;
//    }
//    else {
        // limiting for stability
        visc=std::max(minVisc,std::min(maxVisc,nuApp+nuTurb));
//        if (flag) {
//            visc=plasticVisc;
//        }
//    }
}

void node::resetViscosity(const double& initViscosity){

//    else {
        // limiting for stability
        visc=initViscosity;
//        if (flag) {
//            visc=plasticVisc;
//        }
//    }
}

void node::solveCollision(const double feq[]){
    // relaxation frequency
    const double omega=1.0/(0.5+3.0*visc*lbmDt);

    for (int j=0; j<lbmDirec; j++){
        f[j]+=omega*(feq[j]-f[j]);
    }
}

void node::addForce(const tVect& F){
    static double F1=3.0*lbmDt*lbmDt;
    static double F2=9.0*lbmDt*lbmDt*lbmDt*lbmDt;

    tVect vmu, forcePart;

    const double omegaf=1.0-1.0/(1.0+6.0*visc*lbmDt);
    const tVect totalForce=F+hydroForce;

    for (int j=0; j<lbmDirec; j++){
        const double vu=u.dot(v[j]);
        const tVect vmu=v[j]-u;
//        vmu-=u;
        const tVect forcePart=F2*vu*v[j]+F1*vmu;
//        forcePart+=F1*vmu;
        f[j]+=lbmDt*omegaf*coeff[j]*forcePart.dot(totalForce);
    }
}

/*
void node::collideFNN(const double& turbConst, const tVect& F, const double& plasticVisc, const double& yieldStress) {

    // equilibrium distributions
    double feq[lbmDirec];

    // shift velocity field to F/2
    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
    computeShearRate(feq, turbConst, plasticVisc, yieldStress);
    // compute new distributions
    solveCollision(feq);
    // add force term to new distributions
    addForce(F);
}

void node::collideNN(const double& turbConst, const double& plasticVisc, const double& yieldStress) {

    // equilibrium distributions
    double feq[lbmDirec];

//    // shift velocity field to F/2
//    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
    computeShearRate(feq, turbConst, plasticVisc, yieldStress);
    // compute new distributions
    solveCollision(feq);
//    // add force term to new distributions
//    addForce(F);

}

void node::collideF(const tVect& F, const double& initVisc) {

    // equilibrium distributions
    double feq[lbmDirec];

    // shift velocity field to F/2
    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
//    computeShearRate(feq, plasticVisc, yieldStress, minVisc, maxVisc);
    visc=initVisc; // for NEWTONIAN - viscosity part
    // compute new distributions
    solveCollision(feq);
    // add force term to new distributions
    addForce(F);

}

void node::collide(const double& initVisc) {

    // equilibrium distributions
    double feq[lbmDirec];

    // shift velocity field to F/2
//    shiftVelocity(F);
    // compute equilibrium distributions
    computeEquilibrium(feq);
    // compute shear rate tensor, find invariant calculate viscosity (Bingham)
//    computeShearRate(feq, plasticVisc, yieldStress, minVisc, maxVisc);

//    double tau=0.5+3.0*visc*dt;
//    tMat gamma(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
//    for (int j=0; j<direc; ++j) {
//        gamma+=(f[j]-feq[j])*vv[j];
//    }
//    gamma*=1.5*dt/(tau*n);
//    // shear rate (second invariant)
//    shearRate=gamma.magnitude();

    visc=initVisc; // for NEWTONIAN - viscosity part
    // compute new distributions
    solveCollision(feq);
    // add force term to new distributions
//    addForce(F);

}
*/

void node::store() {
    for (int j=0; j<lbmDirec; ++j) {
            fs[j]=f[j];
    }
}

//void node::stream(unsigned int& fdir, node& linknode, unsigned int& linkdir, double& prod, double& sum){
//    // this is a general expression for the streaming, where
//    // a population can change direction and size (Bounce-Back & Free Surface)
//    f[fdir]=prod*linknode.fs[linkdir]+sum;
//}

tVect node::bounceBackForce(const unsigned int& j, const double staticPres[], const double& BBi) const {
//    tVect force=(2.0*fs[j]-2.0*coeff[j])*dt*v[j];
//    cout<<"AAAA "<<force.norm()<<" "<<n<<" "<<n*force.norm()<<"\n";
//    return (2.0*fs[j]-BBi)*dt*v[j];//*(n-1.0)/n;
    return (2.0*(fs[j]-staticPres[j])-BBi)*lbmDt*v[j];
//    return (2.0*fs[j]-BBi)*mass*mass/n*dt*v[j];
//    return force;
//    return -staticPres[j]*v[j];
}

double node::massStream(const unsigned int& sdir) const {
    return (f[opp[sdir]]-fs[sdir]);
}

// functions for identification of type

bool nodeType::isActive() const {
    if (t==0) {
        return true;
    }
    else if (t==3) {
        return true;
    }
    else
    return false;
}

bool nodeType::isSlipStatWall() const {
    return t==5;
}

bool nodeType::isSlipDynWall() const {
    return t==6;
}

bool nodeType::isWall() const {
    return ((t==5)||(t==6)||(t==7)||(t==8)||(t==9));
}

bool nodeType::isStatWall() const {
    return t==7;
}

bool nodeType::isDynWall() const {
    return t==8;
}

bool nodeType::isCurvedWall() const {
    return t==9;
}

bool nodeType::isInterface() const {
    return t==3;
}

bool nodeType::isFluid() const {
    return t==0;
}

bool nodeType::isGas() const {
    return t==2;
}

bool nodeType::isPeriodic() const {
    return t==4;
}

// functions for change of type

void nodeType::setSlipStatWall(){
    t=5;
}

void nodeType::setSlipDynWall(){
    t=6;
}

void nodeType::setStatWall(){
    t=7;
}

void nodeType::setDynWall(){
    t=8;
}

void nodeType::setCurvedWall(){
    t=9;
}

void nodeType::setInterface(){
    t=3;
}

void nodeType::setFluid(){
    t=0;
}

void nodeType::setGas(){
    t=2;
}

void nodeType::setPeriodic(){
    t=4;
}

void nodeType::setType(unsigned int& typ){
    // setting type according to identification number
    if (typ==0){
        setFluid();
    }
    else if (typ==1){
        cout<<"some node has index 1"<<endl;
        exit(1);
    }
    else if (typ==2){
        setGas();
    }
    else if (typ==3){
        setInterface();
    }
    else if (typ==4){
        setPeriodic();
    }
    else if (typ==5){
        setSlipStatWall();
    }
    else if (typ==6){
        setSlipDynWall();
    }
    else if (typ==7){
        setStatWall();
    }
    else if (typ==8){
        setDynWall();
    }
    else if (typ==9){
        setCurvedWall();
    }
    else {
        cout<<"type not defined\n";
        exit(0);
    }
}

// particle flag

bool nodeType::isInsideParticle() const {
    return p;
}

void nodeType::setInsideParticle() {
    p=true;
}

void nodeType::setOutsideParticle() {
    p=false;
}

// get type

unsigned int nodeType::getType() const {
    return t;
}

// solid index functions
unsigned int nodeType::getSolidIndex() const {
    return solidIndex;
}

void nodeType::setSolidIndex(unsigned int& ind){
    solidIndex=ind;
}

// CURVED ////////////////////////////////

void curve::computeCoefficients() {
    for (int j=1; j<lbmDirec; ++j) {
        m1[j]=(delta[j]-1)/delta[j];
        m2[j]=1/delta[j];
    }
}

double curve::getChi(const unsigned int& j, const double& tau) const {
    if (delta[j]>=0.5) {
        return (2.0*delta[j]-1.0)/tau;
    }
    else {
        return (2.0*delta[j]-1.0)/(tau-2.0);
    }
}

// UNITS ///////////////////////////////////

void measureUnits::setComposite(){
    Volume=Length*Length*Length;
    Speed=Length/Time; //+2
    Accel=Length/Time/Time; // +6
    AngVel=1.0/Time;
    KinVisc=Length*Length/Time; // 0
    DynVisc=Density*Length*Length/Time; // 0
    Force=Density*Length*Length*Length*Length/Time/Time; // +3
    Torque=Density*Length*Length*Length*Length*Length/Time/Time; // +1
    Mass=Density*Length*Length*Length; // -3
    Stress=Density*Length*Length/Time/Time;
    Pressure=Density*Length*Length/Time/Time;
}