//
// Created by xos81802 on 13/07/2018.
//

#ifndef SUGATAMA_FUNCTIONS_H
#define SUGATAMA_FUNCTIONS_H


#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <numeric>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <sstream>

typedef std::vector<float> vector1D;
typedef std::vector<vector1D> vector2D;
typedef std::vector<vector2D> vector3D;

typedef std::vector< std::complex<float> > vector1DC;
typedef std::vector< vector1DC> vector2DC;
typedef std::vector< vector2DC> vector3DC;



float asf ( unsigned int atomicNumber, float q);

float assignOccupancy (  std::string  * neighboringAtom, std::string * neighboringResi);
//float dmaxFromPDB(float x[], float y[], float z[], int numberOfElements, float *pcenterX, float *pcenterY, float *pcenterZ);

//std::vector < std::vector < float > > calcFQ ( std::vector < float >  qvalues, std::set<int> atomList );
void calcFQ (const std::vector<float> &qvalues, std::set<int> atomList, float * f_q_array, int atomSize, float fprime[] );

unsigned int convertToAtomicNumber(const std::string atomType,  const std::string resiName);


void dmaxFromPDB(const std::vector <float>& x,
                 const std::vector <float>& y,
                 const std::vector <float>& z,
                 float * dmax,
                 float * centeredX,
                 float * centeredY,
                 float * centeredZ,
                 size_t const numAtoms
);

std::string newFilename(std::string currentFile, int loop, std::string extension);

double residueToVolume(std::string atomType, std::string residue);

void sub_select(std::vector <float> & qvalue, std::vector<std::vector <float> > & iobs, float subset_fraction, int * qvaluesSize, int lmax,  float qmin, float qmax);

void dmaxFromCenteredCoordinates (float * centeredX,
                                  float * centeredY,
                                  float * centeredZ,
                                  size_t const numAtoms,
                                  float * dmax);
void waterLattice(
        float const * const x,
        float const * const y,
        float const * const z,
        float lowerBound,
        float upperBound,
        float waters[][5],
        //                  vector2D &waters,
        size_t const xSize,
        int * numWaters
);

float * xyz_to_rtp(float const &tempx, float const &tempy, float const &tempz);
float kurtosis(float * residuals, int residualsSize);

float median(std::vector<float> * scores);

void printInfo(std::string text);
void printError(std::string text);

#endif //SUGATAMA_FUNCTIONS_H
