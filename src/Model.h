//
// Created by xos81802 on 10/07/2018.
//

#ifndef SUGATAMA_MODEL_H
#define SUGATAMA_MODEL_H

#include <string>
#include <vector>
#include <math.h>
#include <regex>
#include <iostream>
#include <Bead.h>
#include <set>
#include <Eigen/Core>
#include <SubUnit.h>

const float sqrt6 = std::sqrt(6.0d);//2.449489742783178;
const float invsqrt6 = (float)(1.0d/std::sqrt(6.0d));//0.4082482904638631;
const float sqrt3 = (float)std::sqrt(3.0d);//1.7320508075688772;
const float inv3 = (float)(1.0d/3.0d);//0.3333333333333333;

class Data;
class Anneal;

class Model {

    float radius_of_universe, bead_radius, cutOffNeighbor;
    float xaxis, yaxis, zaxis;
    float inv_bead_radius, helicalVolume;
    float beadAverage, beadStDev, cvx_volume;
    float averageNumberOfContactsInModel;

    unsigned int sizeOfNeighborhood = 12;
    unsigned int number_of_beads, total_in_reduced_seed, total_in_working_universe, neighborLimit, total_in_seed; // there will be no bead with this index
    unsigned long int totalDistances, startingWorkingLimit, workingLimit, deadLimit;
    std::vector<unsigned short int> bins; // linear array (n*(n-1)/2) should this be unsigned short int, bins will never be greater than 100
    float sasa_volume_start=0, sasa_volume=0;
    float sasa_surface_area=0;
    float surface_to_volume=0.0f;

    float radial_limit, bead_volume;

    bool useDirectMethod = false, useCylindricalSearchSpace=false;

    std::string symmetry="C1", helical_point_group_symmetry = "C1";
    float risePerSubUnit, rotationPerSubunit;

    unsigned int maxHelicalIndex;
    unsigned int numberOfSubUnits=1, numberOfHelicalSubUnits, numberOfHelicalSubUnitsInSingleFilament;
    unsigned int symmetryIndex=1;
    std::string symmetryGroup="";
    std::vector<unsigned int> neighboringSubunits; // subunits that correspond to neighbor of the asymmetry unit (zero order)
    std::vector<Bead> beads;
    std::vector<unsigned int> bead_indices;
    std::vector<unsigned int> neighbors; // 12 * bead_indices size
    std::vector<unsigned int> startingSetVector;
    std::vector<unsigned int> seed_indices;  // contains true model of PDB model converted to lattice model
    std::set<unsigned int> reduced_seed;

    std::vector<float> subUnitAnglesCos;
    std::vector<float> subUnitAnglesSin;

    std::vector<SubUnit> subUnits;

    void createUniverse(bool useSphericalModel);
    void createRectangularUniverse(bool useSphericalModel);

    void createUniverseFromMask(std::string filename);
    void createUniverseFromMaskSym(std::string filename);
    void createDistancesAndConvertToSphericalCoordinates();
    void setSymmetryParameters(std::string sym);
    void createSubUnitAngles();
    void pruneUniverseToAsymmetricUnit();

    void logger(std::string description, std::string value) {

        unsigned int len = 40 - description.size();
        std::string firsthalf = std::string(len, ' ');
        printf("%s%s : %s\n", firsthalf.c_str(), description.c_str(), value.c_str());

    }

    std::string formatNumber(float number, int decimals = 2) {
        char buffer [50];
        switch(decimals){
            case 1 :
                sprintf (buffer, "%.1f", number);
                break;
            case 2 :
                sprintf (buffer, "%.2f", number);
                break;
            case 3 :
                sprintf (buffer, "%.3f", number);
                break;
            default:
                sprintf (buffer, "%.4f", number);
        }

        return std::string(buffer);
    }


public:

    Model();
//    Model(float size, float bead_r, bool fastmode, bool shape, int xaxis, int yaxis, int zaxis);
//    Model(float size, float bead_r, bool fastmode, bool shape, int xaxis, int yaxis, int zaxis, std::string sym);

    Model(float beadradius, float xaxis, float yaxis, float zaxis, std::string sym);
    Model(float beadradius, float height, float radius, std::string sym);
    Model(float searchSpace, float beadradius, std::string sym);
    Model(std::string maskfile, float searchSpace, float beadradius, std::string sym, bool maskIt);

    Model(float beadradius, float xaxis, float yaxis, float zaxis);
    Model(float beadradius, float xmin, float ymin, float zmin, float xmax, float ymax, float zmax, std::string sym);
    Model(float beadradius, float height, float radius);
    Model(float searchSpace, float beadradius);
    Model(std::string maskfile, float searchSpace, float beadradius);

    Model(std::string helicalfile, float beadradius);

    bool isUseDirectMethod() const;

    unsigned short int populateBins(Data * pData);
    unsigned short int * getPointerToBins(){ return &bins[0];}
    unsigned int getSizeOfNeighborhood(){ return sizeOfNeighborhood;}

    unsigned int getSizeOfSubUnitNeighbors(){ return neighboringSubunits.size();}
    std::vector<unsigned int> * const getPointerToSubUnitNeighbors(){ return neighboringSubunits.size() > 0 ? &neighboringSubunits : nullptr;}

    unsigned int getNeighborLimit() { return neighborLimit;}
    unsigned int getTotalInSeed(){ return this->total_in_seed; }
    void createSeedFromPDB(std::string filename, Data * pData, unsigned int totalBins, std::vector<double> * pdbPr);

    float getBeadRadius(){ return bead_radius;}
    float getBeadVolume() const {return bead_volume;}
    Bead * getBead(int i) { return &beads[i];}

    //const float * getrThetaPhiAtomType() const {return rThetaPhiAtomType;}
    unsigned int getTotalNumberOfBeadsInUniverse() const {return number_of_beads;}
    unsigned int getNumberOfSubUnits(){return numberOfSubUnits;}
    unsigned int getTotalNumberOfHelicalSubUnits(){return numberOfHelicalSubUnits;}
    unsigned int getTotalNumberOfSubUnitsInSingleFiber(){return numberOfHelicalSubUnitsInSingleFilament;}

    void setNumberOfSubUnitsHelical(unsigned int value){numberOfHelicalSubUnits = value*numberOfSubUnits; numberOfHelicalSubUnitsInSingleFilament = value;}

    void updateUniverse(std::string filename);

    std::vector<unsigned int>::iterator getPointerToNeighborhood(unsigned int index);
    const unsigned int * getDirectPointerToNeighborhood();
    std::string createHeader(float dkl, Anneal *annealedObject, Data *pData, unsigned int steps, unsigned int workingNumber, float volume, float averageContacts);
    void writeModelToFile(const unsigned int workingLimit, std::vector<unsigned int> &beads, std::string & nameOf, unsigned int steps);
    void writeModelToFileFlipped(const unsigned int workingLimit, std::vector<unsigned int> &beads, std::string & nameOf, unsigned int steps);
    void writeSubModelToFile(unsigned int startIndex, unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string nameOf);
    std::string writeSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::vector<unsigned int> &pofrModel, std::string name, Anneal *annealedObject, Data *pData, unsigned int totalSteps, float volume, float averageContacts);
    std::string writeBasicSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::string name);
    std::string writeHelicalSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::vector<unsigned int> &pofrModel, std::string name, Anneal *annealedObject, Data *pData, unsigned int totalSteps, float volume, float averageContacts);
    void writeSetToFile(std::set<unsigned int> &selectedBeads, std::string & nameOf);

    void centerLatticeModel(unsigned int * workingLimit, std::vector<unsigned int> & indices, std::set<unsigned int> & hullPts);

    std::string writeModelToFileBare(float dkl, unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::vector<unsigned int> &pofrModel, std::string nameOf, Anneal *annealedObject, unsigned int steps, float volume, float averageContacts);
    std::string writeModelToFile2(float dkl, unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::vector<unsigned int> &pofrModel, std::string nameOf, Anneal *annealedObject, Data *pData, unsigned int steps, float volume, float averageContacts);
    std::string writeSymXModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::vector<SubUnit> &subunits, std::vector<unsigned int> &pofrModel, std::string name, Anneal *annealedObject, Data *pData, int totalSteps, float volume, float averageContacts);
    std::string getSymmetryString(){ return symmetry;}
    /**
     *
     */
    void setStartingSet(std::vector<unsigned int> & beads){
        total_in_working_universe = beads.size();
        startingSetVector.resize(total_in_working_universe);
        std::copy(beads.begin(), beads.end(), startingSetVector.begin());
    }

    void setStartingWorkingLimit(unsigned int limit){ startingWorkingLimit = limit;}
    unsigned int getStartingWorkingLimit(){return startingWorkingLimit;}

    void setBeadAverageAndStdev(float volume, float stdev);
    float getVolumeAverage(){return beadAverage;}
    float getVolumeStdev(){return beadStDev;}
    void setCVXHullVolume(float volume){ this->cvx_volume = volume; }
    float getCVXHullVolume () const { return this->cvx_volume; }
    void setAverageNumberOfContactsInModel(float number){ this->averageNumberOfContactsInModel = number; }

    inline void printAtomLine(FILE * pFile, unsigned int index, std::string chain, unsigned int residue_index, float x, float y, float z){
        std::string resid = std::to_string(residue_index);
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", index," CA ", "ALA", chain.c_str(), resid.c_str(), x, y, z );
    }

    /*
     * this is sorted!
     */
    void copyStartingModelIntoVector(std::vector<unsigned int> & vec){std::copy(startingSetVector.begin(), startingSetVector.end(), vec.begin());}

    bool isSubUnitNeighbor(unsigned int value);
    void printBeadsFromSet(std::set<unsigned int> &beadIDs);

    const std::vector<unsigned int>::const_iterator getSeedBegin() const { return seed_indices.cbegin(); }
    const std::vector<unsigned int>::const_iterator getSeedEnd() const { return seed_indices.cend(); }
    const std::set<unsigned int>::const_iterator getReducedSeedBegin() const { return reduced_seed.cbegin(); }
    const std::set<unsigned int>::const_iterator getReducedSeedEnd() const { return reduced_seed.cend(); }

    void transformCoordinatesBySymmetryPreCalc(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates);
    void transformCoordinatesByHelicalSymmetry(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates);
    void calculateSymmetryMatesForSingleCoordinatePreCalc(const unsigned int beadIndex, unsigned int &startCount, std::vector<Eigen::Vector3f> &coordinates);

    vector3 transformVectorBySymmetry(unsigned int subunitIndex, const vector3 & vec);
    void transformCoordinatesBySymmetry(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates);

    unsigned short int getMaxBin(Data * pData);

    void createSeedFromPDBSym(std::string filename, Data * pData, unsigned int totalBins, std::vector<double> * pdbPr);

    void setReducedSeed(unsigned int limit, std::vector<unsigned int> &reduced){
        //reduced_seed.resize(limit);
        total_in_reduced_seed = limit;

        reduced_seed.clear();
        for(unsigned int i=0; i<limit; i++){
            reduced_seed.insert(reduced[i]);
        }
    }

    float getRadialLimit(){return radial_limit;}
    unsigned int getTotalInReducedSeed(){ return this->total_in_reduced_seed; }
    void copyReducedModelIntoVector(std::vector<unsigned int> & vec){std::copy(reduced_seed.begin(), reduced_seed.end(), vec.begin());}


    float getDiameterOfUniverse(){ return 2.0*radial_limit;}

    void copySubUnits(std::vector<SubUnit> & copyTo) {
        unsigned int total = subUnits.size();
        copyTo.resize(total);
        std::copy(subUnits.begin(), subUnits.end(), copyTo.begin());
    }

    void updateSubUnits(std::vector<SubUnit> & subunits);
    void updateBeadIndices(unsigned int workingLimit, std::vector<unsigned int> &indices);

    float getHelicalVolume(){return helicalVolume;}
    float getRisePerSubUnit(){ return risePerSubUnit;}
    float getRotationPerSubUnit() { return rotationPerSubunit;}
    void setRiseAndRotationPerSubUnit(float rise, float theta){ this->risePerSubUnit = rise; this->rotationPerSubunit = theta;}
    float getZaxis(){ return zaxis;}
    float getXaxis(){ return xaxis;}

    float getNieghborCutOffLimit(){ return cutOffNeighbor;}
};


#endif //SUGATAMA_MODEL_H
