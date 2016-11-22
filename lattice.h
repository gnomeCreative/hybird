/* 
 * File:   lattice.h
 * Author: aleonard
 *
 * Created on November 6, 2013, 11:43 AM
 */

#ifndef LATTICE_H
#define	LATTICE_H

#include "vector.h"

// dimensions
#define lbmDim 3
// directions
#define lbmDirec 19
#define lbmMainDirec 7

// LATTICE PARAMETERS  ////////////////
// definition of parameter for the lattice
// the implemented lattice is the D3Q19
// changing this would, theoretically, change the lattice structure
// with no consequences on the rest of the program

// time step size
# define lbmDt 1.0

// minimum and maximum relaxation time
#define minTau 0.501
#define maxTau 1.8


// direction vectors of the lattice D3Q19
// order is O,x,y,z,xy,yz,zx.

const tVect v[lbmDirec]={     (1/lbmDt)*tVect(0.0,0.0,0.0), //0
                                  //
                                  (1/lbmDt)*tVect(1.0,0.0,0.0), // 1
                                  (1/lbmDt)*tVect(-1.0,0.0,0.0), // 2
                                  //
                                  (1/lbmDt)*tVect(0.0,1.0,0.0), // 3
                                  (1/lbmDt)*tVect(0.0,-1.0,0.0), // 4
                                  //
                                  (1/lbmDt)*tVect(0.0,0.0,1.0), // 5
                                  (1/lbmDt)*tVect(0.0,0.0,-1.0), // 6
                                  //
                                  (1/lbmDt)*tVect(1.0,1.0,0.0), // 7
                                  (1/lbmDt)*tVect(-1.0,-1.0,0.0), // 8
                                  (1/lbmDt)*tVect(-1.0,1.0,0.0), // 9
                                  (1/lbmDt)*tVect(1.0,-1.0,0.0), // 10
                                  //
                                  (1/lbmDt)*tVect(0.0,1.0,1.0), // 11
                                  (1/lbmDt)*tVect(0.0,-1.0,-1.0), // 12
                                  (1/lbmDt)*tVect(0.0,-1.0,1.0), // 13
                                  (1/lbmDt)*tVect(0.0,1.0,-1.0), // 14
                                  //
                                  (1/lbmDt)*tVect(1.0,0.0,1.0), // 15
                                  (1/lbmDt)*tVect(-1.0,0.0,-1.0), // 16
                                  (1/lbmDt)*tVect(1.0,0.0,-1.0), // 17
                                  (1/lbmDt)* tVect(-1.0,0.0,1.0) }; // 18

// direction vectors of the lattice D2Q9
// tensor v_i x v_i
const tMat vv[lbmDirec]={    tMat(v[0],v[0]),
                                    tMat(v[1],v[1]),
                                    tMat(v[2],v[2]),
                                    tMat(v[3],v[3]),
                                    tMat(v[4],v[4]),
                                    tMat(v[5],v[5]),
                                    tMat(v[6],v[6]),
                                    tMat(v[7],v[7]),
                                    tMat(v[8],v[8]),
                                    tMat(v[9],v[9]),
                                    tMat(v[10],v[10]),
                                    tMat(v[11],v[11]),
                                    tMat(v[12],v[12]),
                                    tMat(v[13],v[13]),
                                    tMat(v[14],v[14]),
                                    tMat(v[15],v[15]),
                                    tMat(v[16],v[16]),
                                    tMat(v[17],v[17]),
                                    tMat(v[18],v[18]) };

// opposed directions (used for bounce back)
const unsigned int opp[lbmDirec]={0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};

// slip directions (used for slip walls)
const unsigned int slip1CheckSup[lbmDirec]={0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 3, 4, 5, 6, 1, 2, 6, 5};
const unsigned int slip1Sup[lbmDirec]={0, 0, 0, 0, 0, 0, 0, 9, 10, 8, 7, 13, 14, 12, 11, 18, 17, 15, 16};
const unsigned int slip2CheckSup[lbmDirec]={0, 0, 0, 0, 0, 0, 0, 3, 4, 2, 1, 5, 6, 4, 3, 5, 6, 1, 2};
const unsigned int slip2Sup[lbmDirec]={0, 0, 0, 0, 0, 0, 0, 10, 9, 7, 8, 14, 13, 11, 12, 17, 18, 16, 15};
const unsIntList slip1Check(slip1CheckSup,slip1CheckSup+19);
const unsIntList slip1(slip1Sup,slip1Sup+19);
const unsIntList slip2Check(slip2CheckSup,slip2CheckSup+19);
const unsIntList slip2(slip2Sup,slip2Sup+19);

// weight coefficients for the D2Q9 lattice (universal principle of laziness)
const double coeff[lbmDirec]={       12.0/36.0,
                                                              2.0/36.0, 2.0/36.0, 2.0/36.0, 2.0/36.0, 2.0/36.0, 2.0/36.0,
                                                              1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
                                                              1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };


#endif	/* LATTICE_H */

