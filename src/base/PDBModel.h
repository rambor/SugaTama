//
// Created by xos81802 on 13/07/2018.
//

#ifndef SUGATAMA_PDBMODEL_H
#define SUGATAMA_PDBMODEL_H


#include <boost/regex.hpp>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <numeric>
#include <functional>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "functions.h"
#include <math.h> // debug
#include <fstream>
#include <sstream>
//#include "Waters.h"
#include "Coords.h"
#include "vector3.h"
#include <dirent.h>
#include <boost/lexical_cast.hpp>

const float invFourPI = (float)(1.0f/(4.0f*M_PI));

class PDBModel {

private:
    std::string filename;
    float edge_radius=0.0;

    bool ifRNA, found_edge_radius=false;
    bool ifWaters;
    unsigned int totalAtoms, totalResidues, waterCount, watersPerResidue;  // maximum number of waters per residue, should match declaration in waters

    // Following vectors will all be same length
    std::vector < std::string > atomType, resi, trimmedAtomType, trimmedResi, chainID, waterLines;
    // waters

    std::vector < unsigned int > resID;
    std::vector < float > x,y,z, atomVolume, atomicRadii;

    Coords * keptWaters = nullptr;

    std::set <unsigned int> atomList;

    float volume, dmax, cutoff;

    // setup for dynamically allocated arrays
    float * rThetaPhiAtomType = nullptr;
    float * atomicExpVolume = nullptr;
    unsigned int * atomicNumbers = nullptr;

    float * rWater = nullptr;
    float * phiWater = nullptr;
    float * waterCosTheta = nullptr;

    float * centeredX = nullptr; // new Array declared on heap
    float * centeredY = nullptr;
    float * centeredZ = nullptr;

    float centeringX=0.0f, centeringY=0.0f, centeringZ=0.0f;
    vector3 centeringVector;

    void setWaterRPhi(int index, float rValue, float phiValue);

public:
    PDBModel(std::string filename, bool rna, bool discardWaters, float lower);

    PDBModel(const PDBModel &e2);

    ~PDBModel(){

        delete[] rThetaPhiAtomType;
        rThetaPhiAtomType = nullptr;
        delete[] atomicExpVolume;
        atomicExpVolume = nullptr;
        delete[] atomicNumbers;
        atomicNumbers= nullptr;

        delete[] centeredX;
        centeredX = nullptr;
        delete[] centeredY;
        centeredY = nullptr;
        delete[] centeredZ;
        centeredZ = nullptr;

        delete[] keptWaters;
        keptWaters = nullptr;

        if (ifWaters){
            delete[] rWater;
            rWater=nullptr;
            delete[] phiWater;
            phiWater=nullptr;
            delete[] waterCosTheta;
            waterCosTheta=nullptr;
        }
    } // destructor defined in header file


    void writeCenteredCoordinates() const;
    void writeKeptWaters() const;
    //void addWaters(Waters & waterModel);

    unsigned int getWaterCount() const { return waterCount;}

    unsigned int getTotalAtoms() const { return totalAtoms;}
    float getDmax() const { return dmax; }
    float getVolume() const { return volume; }
    const std::set <unsigned int> &getAtomList() const { return atomList;}

    const float * getrThetaPhiAtomType() const {return rThetaPhiAtomType;}
    const float * getAtomicExpVolume() const {return atomicExpVolume;}

    const float * getAtomVolume() const {return &atomVolume[0];}

    const unsigned int * getAtomicNumbers() const {return atomicNumbers;}
    const float * getCenteredX() const {return centeredX;}
    const float * getCenteredY() const {return centeredY;}
    const float * getCenteredZ() const {return centeredZ;}
    const float * getWaterR() const {return &rWater[0];}
    const float * getWaterPhi() const {return &phiWater[0];}
    const float * getWaterCosTheta() const {return &waterCosTheta[0];}

    float getXCoord(int i){ return x[i];}
    float getYCoord(int i){ return y[i];}
    float getZCoord(int i){ return z[i];}

    float getAtomicRadius(unsigned int i){return atomicRadii[i];}

    // resID;
    const std::vector<unsigned int>::const_iterator getResIDIterator() const { return resID.cbegin(); }
    const std::vector<std::string>::const_iterator getAtomTypeIterator() const { return trimmedAtomType.cbegin(); }
    const std::vector<std::string>::const_iterator getChainIDIterator() const { return chainID.cbegin(); }

    //void setCutoff(float number){ cutoff = number*number;} // in angstroms
    const Coords * getKeptWaters() const {return &keptWaters[0];}

    void writeCenteredCoordinatesToFile(std::string name) const;

    std::string getAtomTypeByIndex(int index){ return atomType[index];}

    bool getEdgeRadiusStatus(){return found_edge_radius;}

    float getEdgeRadius(){return edge_radius;}

    const vector3 & getCenteringVec(){ return centeringVector;}

};
#endif //SUGATAMA_PDBMODEL_H
