//
// Created by xos81802 on 03/02/2019.
//

#ifndef SUGATAMA_KDE_H
#define SUGATAMA_KDE_H

#include "../thirdparty/Eigen/Core"
#include "../thirdparty/Eigen/Geometry"
#include "../thirdparty/power_sasa.h"

#include <string>
#include <vector>
#include <algorithm>
#include <vector3.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iostream>
#include <cstdio>
#include <utility>
#include <limits>
#include <fstream>
#include "../Model.h"
#include "../base/Data.h"

class KDE {

    struct Trial{ double score; std::vector<unsigned int> indices;};

    bool priorSet = false;
    std::string prior_name;
    unsigned int minN=10000, maxN=0, kdeSUM;
    std::string listofPDBFiles;
    std::vector<vector3> coordinates;
    std::vector<vector3> prior;
    std::vector<vector3> centered_prior;
    std::vector<vector3> centered_coordinates;
    std::vector<std::string> pdbfiles;
    std::map<unsigned int, std::vector<unsigned int> > bead_to_vec_map;
    std::map<std::string, float > kernel_mask_mapping;
    std::map<std::string, float > kernel_mapping;
    vector3 centering_vec;
    unsigned int totalCoordinates, totalFiles, trialSize, topNSize;
    float bead_radius, topNpercent;
    float minx =  std::numeric_limits<float>::max();
    float miny =  std::numeric_limits<float>::max();
    float minz =  std::numeric_limits<float>::max();
    float maxx =  -std::numeric_limits<float>::max();
    float maxy =  -std::numeric_limits<float>::max();
    float maxz =  -std::numeric_limits<float>::max();
    float cminx, cminy, cminz, cmaxx, cmaxy, cmaxz;
    float kdeCount, grid_spacing;
    float bandwidth, delta_r = 1.8;
    float average_dmin, dmin_supremum, dmin_infimum, stdev_dmin;
    float solventContrast = 0.334;
    float proteinContrast = 0.414;
    float nucleicacidContrast = 0.621;

    double mapAverage, stdev;

    float remapToNewHCPLattice(Model * pModel, float cutoff);
    float remapToNewHCPLatticeSym(Model * pModel, float cutoff);
    void createMaskOnCenteredObject(std::vector<unsigned int> & mask, float limit);
    std::string headerParametersForXPLOR(int & na, int & nb, int & nc, float & adjustedCMINX, float & adjustedCMINY, float & adjustedCMINZ);
    std::string headerParametersForXPLORSYM(int & na, int & nb, int & nc, float & adjustedCMINX, float & adjustedCMINY, float & adjustedCMINZ, std::string sym);
    std::string mirrorTempHeader;
    void logger(std::string description, std::string value);
    std::string formatNumber(float number, int decimals);

    unsigned int estimateRounds(unsigned int lowerBound, unsigned int totalInSet);
public:
    struct ProbabilityBead{ double prob; unsigned int index; double prior;
        bool operator<(const ProbabilityBead & a) const{
            return prob < a.prob;
        }
    };

    KDE(std::string filename, float bead_radius, float topNpercent);

    std::string getFileExt(const std::string &s);

    bool checkKDEFile(std::string basic_string);
    bool createKDE(Model * pModel);

    void extractCoordinates();

    float getCenteredMaxZ() {return std::abs(cminz) > std::abs(cmaxz) ? std::abs(cminz) : std::abs(cmaxz);}
    float getCenteredMaxX() {return std::abs(cminx) > std::abs(cmaxx) ? std::abs(cminx) : std::abs(cmaxx);}
    float getCenteredMaxY() {return std::abs(cminy) > std::abs(cmaxy) ? std::abs(cminy) : std::abs(cmaxy);}

    float kernel(float value);

    float map_refine(Model *pModel, Data *pData, std::string outname);
    float map_refineSym(Model *pModel, Data *pData, std::string outname);

    void generateSeedFromInputFile(Model *pModel, Data *pData, std::string outname);

    void writeCoordinatesToFile(std::string name, std::vector<vector3> & coords) const;

    float getMinx(){return minx;}
    float getMiny(){return miny;}
    float getMinz(){return minz;}
    float getMaxx(){return maxx;}
    float getMaxy(){return maxy;}
    float getMaxz(){return maxz;}

    void add_prior(std::string);
    float getAverageDmin(){ return average_dmin;}
    float getSupremumDmin(){ return dmin_supremum;}
    float surfaceVolumeCalculation(unsigned int workingLimit, std::vector<unsigned int> &mask, Model *pModel);
};


#endif //SUGATAMA_KDE_H
