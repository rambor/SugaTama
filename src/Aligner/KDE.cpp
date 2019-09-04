//
// Created by xos81802 on 03/02/2019.
//


#include <boost/filesystem.hpp>
#include <PDBModel.h>
#include "KDE.h"
#include "../Anneal.h"

/**
 * Assume files are aligned as appropriate
 * We will construct a bounding box based on the input coordinates
 * HCP lattice will be used to establish the kernel function for writing out a map
 *
 * check center, but will
 *
 *
 * @param filename
 * @param bead_radius
 */
KDE::KDE(std::string filename, float bead_radius, float topN) : listofPDBFiles(filename), bead_radius(bead_radius), topNpercent(topN) {
    // determine if dat file, it is IofQ or PofR
    std::string ext = this->getFileExt(filename);
    std::cout << "** CHECKING FILE " << ext << std::endl;
    // compile list of coordinates into single array
    if (ext == "inp"){
        if (this->checkKDEFile(filename)){
            // extract coordinates
            this->extractCoordinates();
        } else {
            throw std::invalid_argument("** ERROR => CORRUPTED INPUT FILE : " + filename);
        }
    } else {
        throw std::invalid_argument("** ERROR => INCORRECT INPUT FILE MUST be .inp : " + filename);
    }
}


std::string KDE::getFileExt(const std::string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
        return(s.substr(i+1, s.length() - i));
    }

    return("");
}

/*
 * Assume first file in list is the reference
 */
void KDE::extractCoordinates(){
    float sumx=0, sumy=0, sumz=0;
    for(auto file : pdbfiles){

        PDBModel tempmodel (file, false, true, bead_radius);
        unsigned int totalCoords = tempmodel.getTotalAtoms();

        minN = (totalCoords < minN) ? totalCoords : minN;
        maxN = (totalCoords > maxN) ? totalCoords : maxN;

        for(unsigned int i = 0; i<totalCoords; i++){
            float pX = tempmodel.getXCoord(i);
            float pY = tempmodel.getYCoord(i);
            float pZ = tempmodel.getZCoord(i);
            sumx += pX;
            sumy += pY;
            sumz += pZ;

            coordinates.emplace_back(vector3(pX, pY, pZ));
            if (pX > maxx){
                maxx = pX;
            }
            if (pY > maxy){
                maxy = pY;
            }
            if (pZ > maxz){
                maxz = pZ;
            }
            if (pX < minx){
                minx = pX;
            }
            if (pY < miny){
                miny = pY;
            }
            if (pZ < minz){
                minz = pZ;
            }
        }
    }

    float invT = 1.0f/(float)coordinates.size();
    sumx *= invT;
    sumy *= invT;
    sumz *= invT;
    centering_vec = vector3(sumx, sumy, sumz);

    const vector3 * pVec = coordinates.data();

    for (unsigned int n=0; n < coordinates.size(); n++) {
        const vector3 &pVec1 = *(pVec + n);
        centered_coordinates.emplace_back(pVec1 - centering_vec);
    }


    std::vector<float> dmin;
    float dminSum =0.0f;
    float dminSumSquare =0.0f;
    for (unsigned int n=0; n < coordinates.size(); n++) {
        const vector3 &pVec1 = *(pVec + n);
        float tempDmin = 100;
        for (unsigned int m=0; m < n; m++) {
            const vector3 &pVec2 = *(pVec + m);
            float dis = (pVec1 - pVec2).length();
            if (dis < tempDmin) {
                tempDmin = dis;
            }
        }

        for (unsigned int m=n+1; m < coordinates.size(); m++) {
            const vector3 &pVec2 = *(pVec + m);
            float dis = (pVec1 - pVec2).length();
            if (dis < tempDmin) {
                tempDmin = dis;
            }
        }

        dmin.emplace_back(tempDmin);
        dminSum += tempDmin;
        dminSumSquare += tempDmin*tempDmin;
    }

    std::sort(dmin.begin(), dmin.end());
    dmin_supremum = dmin[dmin.size()-1];
    dmin_infimum = dmin[0];
    average_dmin = dminSum/(float)dmin.size();
    stdev_dmin = std::sqrt(dminSumSquare/(float)dmin.size() - average_dmin*average_dmin);

    logger("min dmin", std::to_string(dmin_infimum));
    logger("max dmin", std::to_string(dmin_supremum));
    logger("<dmin>", std::to_string(average_dmin));
    logger("stdev", std::to_string(stdev_dmin));
    logger("width", std::to_string(2*stdev_dmin));

    /*
     * set bounding box
     */
    minx -= bead_radius;
    miny -= bead_radius;
    minz -= bead_radius;
    maxx += bead_radius;
    maxy += bead_radius;
    maxz += bead_radius;


    cminx = minx - centering_vec.x;
    cminy = miny - centering_vec.y;
    cminz = minz - centering_vec.z;
    cmaxx = maxx - centering_vec.x;
    cmaxy = maxy - centering_vec.y;
    cmaxz = maxz - centering_vec.z;

    this->writeCoordinatesToFile("centered", centered_coordinates);
    logger("Total Coordinates", std::to_string(coordinates.size()));
    logger("Total Centered Coordinates", std::to_string(centered_coordinates.size()));
    logger("Total PDB Files", std::to_string( pdbfiles.size()));
    totalCoordinates = centered_coordinates.size();
}

/*
 * should check of list of files that exists
 *
 */
bool KDE::checkKDEFile(std::string file) {

    bool returnMe = false;
    // read in file
    std::ifstream data (file, std::ifstream::in);

    if (data.is_open()) {

        boost::regex format("pdb");
        std::string line;

        while(!data.eof()) //
        {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if(boost::regex_search(tempLine[0], format)){
                std::cout << "  =>   CHECKING PDB COORDINATE FILE : " << tempLine[0] <<std::endl;
                if (!boost::filesystem::exists(tempLine[0])){
                    throw std::invalid_argument("** ERROR => FILE NOT FOUND : " + tempLine[0]);
                }
                pdbfiles.push_back(tempLine[0]);
            } else {
                std::cout << "NO FILE FOUND FOR LINE  =>  " << line <<std::endl;
            }
        }
    }

    totalFiles = pdbfiles.size();
    data.close();
    returnMe = true;
    return returnMe;
}


bool KDE::createKDE(Model *pModel) {

    float cutoff = std::sqrt(3.0d/2.0d)*average_dmin*1.02; // midpoint of a tetrahedron
    //grid_spacing = 2.0/3.0*bead_radius; // spacing in rectangular grid for calculating map
    grid_spacing = 2*average_dmin; // spacing in rectangular grid for calculating map
    bandwidth = grid_spacing*std::sqrt(3.0d); // diagonal of a square 1.73*r

    const vector3 * pVec = centered_coordinates.data();
    float limit = remapToNewHCPLattice(pModel, cutoff);

    std::vector<unsigned int> mask;
    createMaskOnCenteredObject(mask, limit);

    //write to file
    unsigned int workingLimit = mask.size();
    std::string name = "ave_mask";
    float maskvolume = surfaceVolumeCalculation(workingLimit, mask, pModel);
    pModel->writeModelToFile(workingLimit, mask, name, 0);

    int na, nb, nc;
    float adjustedCMINX, adjustedCMINY, adjustedCMINZ;
    std::string tempHeader = headerParametersForXPLOR(na, nb, nc, adjustedCMINX, adjustedCMINY, adjustedCMINZ);


    float inv2Pi = 1.0f/(float)std::sqrt(2*M_PI);
    float invhd = 1.0f;///(bandwidth*bandwidth*bandwidth);
    float mapSum = 0.0f, mapSumSquared = 0.0f;
    float mapCount = 0.0f;
    float averageInMask = 0.0f;
    float maskCount = 0.0f;

    for(int z=0; z<nc; z++){
        float zsection = adjustedCMINZ+z*grid_spacing;
        for(int y=0; y<nb; y++){
            float ysection = adjustedCMINY+y*grid_spacing;
            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                vector3 locale(adjustedCMINX+x*grid_spacing, ysection, zsection);
                float kernelSum = 0.0f;
                float incount = 0.0f;
                for(unsigned int c=0; c<totalCoordinates; c++){ // calculate distances to vectors in compiled model
                    vector3 vec1 = *(pVec+c) - locale;
                    float length = vec1.length();
                    kernelSum += kernel(length/bandwidth);
                    incount += 1.0f;
                }

                kernelSum *= (incount > 0) ? inv2Pi/incount*invhd : 0;
                kernel_mapping.emplace(key, kernelSum);
                /*
                 * need to determine if grid point is within the mask
                 * if so, determine average within the mask
                 *
                 */
                for(auto mit : mask){
                    Bead * pBead = pModel->getBead(mit);
                    vector3 vec1 = pBead->getVec() - locale;
                    if (vec1.length() < bandwidth){
                        kernel_mask_mapping.emplace(key, kernelSum);
                        averageInMask += kernelSum;
                        maskCount += 1.0;
                    }
                }

                mapSum += kernelSum;
                mapSumSquared += kernelSum*kernelSum;
                mapCount += 1.0;
            }
        }
    }

    /*
     * ADJUST THE MAP SUCH THAT THE KERNEL THAT INTERSECTS THE MASK IS SCALE TO SOMETHING SENSIBLE ON AVERAGE
     * ESSENTIALLY EVERTHING OUTSIDE MASK STAYS THE SAME (LOW VALUE NEAR ZERO) WHEREAS EVERYTHING ELSE IS SCALED TO 0.414
     * SHOULD scale max value and then set average?
     */
    float scale = 0.414/(averageInMask/maskCount);
    char buffer[80];

    for(int z=0; z<nc; z++){

        std::snprintf(buffer, 80, "%8i\n", z);
        tempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){
            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                auto allIt = kernel_mapping.find(key);
                auto maskIt = kernel_mask_mapping.find(key);
                /*
                 * write to file
                 */
                float kernelSum = (*allIt).second;
                if (maskIt != kernel_mask_mapping.end()){
                    kernelSum = (*maskIt).second*scale;
                }

                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                tempHeader.append(buffer);
                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    tempHeader += "\n";
                }
            }
        }

        if (stringIndex % 6 != 0){
            tempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    tempHeader.append(buffer);
    float ave = mapSum/mapCount;
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, std::sqrt(mapSumSquared/mapCount - ave*ave));
    tempHeader.append(buffer);

    // write to file
    std::string map = "map.xplor";
    const char * outputFileName = map.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    return false;
}

/*
 * Gaussian with unit variance
 */
float KDE::kernel(float value){
    return std::exp(-value*value*0.5);
}


void KDE::writeCoordinatesToFile(std::string name, std::vector<vector3> & coords) const {

    char alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ;
    std::vector<char> alphabet( alpha, alpha+sizeof(alpha)-1 ) ;

    std::string residue_index;
    std::string remarkinfo = "KDE " + name;
    name = name + ".pdb";
    std::string atom = "CA";
    std::string resi = "ALA";
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK COORDINATES FOR KDE : %s\n", remarkinfo.c_str());
    std::string chain;
    unsigned int chainIndex = 0;
    unsigned int residueIndex = 1;
    const vector3 * pVec = coords.data();

    for (unsigned int n=0; n < coords.size(); n++) {
        if (n>0 && n % 999 == 0){
            chainIndex++;
            residueIndex = 1;
        }
        chain = alphabet[chainIndex];
        vector3 vec = *(pVec + n);
        residue_index = boost::lexical_cast<std::string>(residueIndex);
        fprintf(pFile, "ATOM  %5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", n+1, atom.c_str(), resi.c_str(), chain.c_str(), residue_index.c_str(), vec.x, vec.y, vec.z );
        residueIndex++;
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}

float KDE::surfaceVolumeCalculation(unsigned int workingLimit, std::vector<unsigned int> & mask, Model * pModel){

    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
    std::vector<float> weights(workingLimit);
    std::fill(weights.begin(), weights.end(), bead_radius + delta_r);
    std::vector<Eigen::Vector3f> coordinates(workingLimit);

    for(unsigned int i=0; i<workingLimit; i++){
        vector3 pBeadVec = pModel->getBead(mask[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
        unsorted_bead_indices[i] = mask[i];
    }

    POWERSASA::PowerSasa<float, Eigen::Vector3f> *ps =
            new POWERSASA::PowerSasa<float, Eigen::Vector3f>(coordinates, weights, 1, 1, 1, 1);

    ps->calc_sasa_all();
    float volume = 0.0;
    float sasa = 0.0;
    for (unsigned int i = 0; i < workingLimit; ++i) {
        sasa += ps->getSasa()[i];
        volume += ps->getVol()[i];
    }

    delete ps; //maximizes surface area while minimizing number of regions
    return volume;
}


/**
 * maps are not centered for models with symmetry
 * First model in the kde.list is the reference to which symmtry operators will be applied
 * @param pModel
 * @param cutoff
 * @return value to use to determine output mask, everything >= limit is kept
 */
float KDE::remapToNewHCPLattice(Model * pModel, float cutoff){

    const vector3 * pVec = centered_coordinates.data();
    /*
     * go through each bead in pModel and assign it a bead from vector3
     * assignment is based on cutoff
     */
    unsigned int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();
    for(unsigned int i = 0; i<totalBeads; i++){
        Bead * pBead = pModel->getBead(i);
        auto pPosition = bead_to_vec_map.find(i);

        for(unsigned int c=0; c<totalCoordinates; c++){
            vector3 vec1 = *(pVec+c) - pBead->getVec();

            if (vec1.length() < cutoff){
                if (pPosition != bead_to_vec_map.end()){
                    (*pPosition).second.push_back(c);
                } else {
                    auto temp = bead_to_vec_map.emplace(i, std::vector<unsigned int>());
                    (*temp.first).second.push_back(c);
                    pPosition = bead_to_vec_map.find(i);
                }
            }
        }
    }
    logger("SIZE OF UNIVERSE", std::to_string(totalBeads));
    logger("TOTAL IN LATTICE", std::to_string(bead_to_vec_map.size()));

    // calculate average per bead and variance of remapped lattice
    kdeSUM = 0.0f;
    float squaredSum = 0.0f;
    kdeCount = 0.0f;
    for(auto mit : bead_to_vec_map){
        float temp = mit.second.size();
        squaredSum += temp*temp;
        kdeSUM += temp;
        kdeCount += 1.0f;
    }


    mapAverage = kdeSUM/kdeCount;
    float variance = squaredSum/kdeCount - mapAverage*mapAverage;
    stdev = std::sqrt(variance);
    float limit = mapAverage;// - 0.25*stdev;

    logger("TOTAL IN REMAPPED HCP LATTICE", std::to_string(bead_to_vec_map.size()));
    logger("AVERAGE PER BEAD", formatNumber(mapAverage, 2));
    logger("VARIANCE", formatNumber(variance, 2));
    logger("STDEV", formatNumber(stdev, 2));

    return limit;
}


/**
 * maps are not centered for models with symmetry
 * First model in the kde.list is the reference to which symmetry operators will be applied
 * dmax of the search space must support the Model
 *
 * Mapping consists of going through the new Universe (pModel)
 * Determining which beads in aligned set of coordinates overlap with pModel
 *
 * @param pModel
 * @param cutoff
 * @return
 */
float KDE::remapToNewHCPLatticeSym(Model * pModel, float cutoff){

    const vector3 * pVec = coordinates.data();
    /*
     * go through each bead in pModel
     * if point in coordinates overlaps with pModel
     * add index of pModel to bead_to_vec_map
     * assignment is based on cutoff
     */
    unsigned int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();
    for(unsigned int i = 0; i<totalBeads; i++){
        Bead * pBead = pModel->getBead(i);
        auto pPosition = bead_to_vec_map.find(i);

        for(unsigned int c=0; c<totalCoordinates; c++){ // iterate through the points that make up the averaging cloud
            vector3 vec1 = *(pVec+c) - pBead->getVec();

            if (vec1.length() < cutoff){
                if (pPosition != bead_to_vec_map.end()){
                    (*pPosition).second.push_back(c);
                } else {
                    auto temp = bead_to_vec_map.emplace(i, std::vector<unsigned int>());
                    (*temp.first).second.push_back(c);
                    pPosition = bead_to_vec_map.find(i);
                }
            }
        }
    }
    // calculate average per bead and variance of remapped lattice
    float sumIt = 0.0f;
    float squaredSum = 0.0f;
    float countIt = 0.0f;
    for(auto mit : bead_to_vec_map){
        float temp = mit.second.size();
        squaredSum += temp*temp;
        sumIt += temp;
        countIt += 1.0f;
    }

    mapAverage = sumIt/countIt;
    float variance = squaredSum/countIt - mapAverage*mapAverage;
    stdev = std::sqrt(variance);
    float limit = mapAverage - 0.25*stdev;


    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************            REMAPPED LATTICE            *******************" << std::endl;

    printf("       TOTAL IN REMAPPED HCP LATTICE : %-i  \n", bead_to_vec_map.size());
    printf("                    AVERAGE PER BEAD : %0.1f  \n", mapAverage);
    printf("                            VARIANCE : %.2f  \n", variance);
    printf("                               STDEV : %.2f  \n", stdev);

    return limit;
}

/**
 * take map and refined model and resample map
 * @return
 */
float KDE::map_refineSym(Model *pModel, Data *pData, std::string outname) {

    std::cout << " SUBUNIT BOUNDING BOX "<< std::endl;
    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    printf("  => X %8.4f %8.4f %8.4f \n", minx, maxx, (maxx-minx));
    printf("  => Y %8.4f %8.4f %8.4f \n", miny, maxy, (maxy-miny));
    printf("  => Z %8.4f %8.4f %8.4f \n", minz, maxz, (maxz-minz));
    /*
     * use coordinates to build symmetric model and determine size of box
     * update size of the bounding box
     */
    unsigned int totalSubunits = pModel->getNumberOfSubUnits();
    for(unsigned int s = 0; s<totalSubunits; s++){
        for(auto coords : coordinates){
            vector3 tempVec = pModel->transformVectorBySymmetry(s, coords);
            float pX = tempVec.x;
            float pY = tempVec.y;
            float pZ = tempVec.z;
            if (pX > maxx){
                maxx = pX;
            }
            if (pY > maxy){
                maxy = pY;
            }
            if (pZ > maxz){
                maxz = pZ;
            }
            if (pX < minx){
                minx = pX;
            }
            if (pY < miny){
                miny = pY;
            }
            if (pZ < minz){
                minz = pZ;
            }
        }
    }
    /*
     * set bounding box
     */
    minx -= bead_radius;
    miny -= bead_radius;
    minz -= bead_radius;
    maxx += bead_radius;
    maxy += bead_radius;
    maxz += bead_radius;

    std::cout << " ENTIRE STRUCTURE BOUNDING BOX "<< std::endl;
    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    printf("  => X %8.4f %8.4f %8.4f \n", minx, maxx, (maxx-minx));
    printf("  => Y %8.4f %8.4f %8.4f \n", miny, maxy, (maxy-miny));
    printf("  => Z %8.4f %8.4f %8.4f \n", minz, maxz, (maxz-minz));

    float cutoff = bead_radius + 0.5*average_dmin; // sphere of influence, looking for all vectors in centered coordnates that are within...
    cutoff = std::sqrt(3.0d/2.0d)*bead_radius*1.02;

    grid_spacing = average_dmin; // spacing in rectangular grid for calculating map
    bandwidth = grid_spacing*std::sqrt(3.0) > (std::sqrt(3.0d/2.0d)*bead_radius) ? std::sqrt(3.0d/2.0d)*bead_radius : grid_spacing*std::sqrt(3.0);


    unsigned int highTempRounds = 23;
    std::string prefix = "CE";

    float limit = remapToNewHCPLatticeSym(pModel, cutoff);

    std::vector<ProbabilityBead> lattice;
    //double prob = 1.0d/bead_to_vec_map.size();
    std::vector<vector3> remapped;
    // bias the initial map based on point density
    for(auto & index : bead_to_vec_map){
        if (index.second.size() >= limit){
            lattice.emplace_back(ProbabilityBead{0.5, index.first, 0.5});
        } else {
            lattice.emplace_back(ProbabilityBead{0.1, index.first, 0.1});
        }
        remapped.emplace_back(pModel->getBead(index.first)->getVec());
    }


    // create vector of indices and a vector to be cdf
    unsigned int pot = (double)estimateRounds(minN,lattice.size());
    trialSize = (pot < 23571) ? 23571 : pot;
    logger("estimated pool", std::to_string(pot));
    logger("estimated eqn", std::to_string((unsigned int)(1.0/(double)minN*lattice.size()*std::log(lattice.size()))));
    unsigned int totalpairs = (lattice.size()*(lattice.size()-1))/2;
    unsigned int pairs = (minN*(minN-1))/2;
    logger("estimated eqn pairs", std::to_string((unsigned int)(1.0/(double)pairs*totalpairs*std::log(totalpairs))));

    trialSize = 10000;//minN*4*40;
    topNSize = (unsigned int)(topNpercent*trialSize);

    Anneal mainAnneal(0,
                      0,
                      highTempRounds,
                      prefix,
                      0.0731, //alpha
                      0.01,  //beta
                      0,    //eta
                      0.01, //lambda
                      0,    //mu
                      trialSize,
                      0,    //accept rate
                      2*average_dmin // interconnectivity
    );

    writeCoordinatesToFile("uncentered", coordinates);


    // if input refined model is added update probabilites of lattice
    if (priorSet){
        std::cout << "       UPDATING MAP WITH PRIOR : " << prior_name << std::endl;
        double sum=0.0d;
        double alpha = 0.7;
        std::vector<double> temp;
        for(auto & lat : lattice){

            double kernelSum = 0.0d;
            for(unsigned int i=0; i<prior.size(); i++){
                vector3 lVec = (prior[i] - pModel->getBead(lat.index)->getVec());
                kernelSum += kernel(lVec.length()/cutoff);
            }
            sum+=kernelSum;
            temp.emplace_back(kernelSum);
        }
        double inv = 1.0d/sum;
        for(unsigned int i=0; i<lattice.size(); i++){ // update lattice probabilities
            ProbabilityBead * pPoint = &lattice[i];
            double old = pPoint->prob;
            pPoint->prob = old*(1.0-alpha) + alpha*temp[i]*inv;
        }
    }

    writeCoordinatesToFile("remapped", remapped);
    unsigned int lowerBound = (unsigned int)std::floor(0.8*minN);

    mainAnneal.ceMapOptimizationSym(pModel, pData, topNSize, lowerBound, maxN, lattice);

    // create map
    unsigned int totalInLattice = lattice.size();
    double sumIt=0.0d;
    double var=0.0d;
    for(auto & point : lattice){
        sumIt+=point.prob;
        var += point.prob*point.prob;
    }

    double average = sumIt/(double)totalInLattice;
    double variance = var/(double)totalInLattice - average*average;

    logger("MAP AVERAGE", formatNumber(average, 4));
    logger("MAP SIGMA", formatNumber(std::sqrt(variance), 4));

    std::vector<unsigned int> mask;
    // make a kept list that is essentially contoured to >= average
    createMaskOnCenteredObject(mask, limit);
    std::string name = outname+"_mask_sym";
    std::string subname = outname+"_mask_sub";
    pModel->writeModelToFile(mask.size(), mask, subname, 0);
    pModel->writeBasicSymModelToFile(0.0f, mask.size(), mask, name);

    int na, nb, nc;
    float adjustedCMINX, adjustedCMINY, adjustedCMINZ;
    std::string tempHeader = headerParametersForXPLORSYM(na, nb, nc, adjustedCMINX, adjustedCMINY, adjustedCMINZ, pModel->getSymmetryString());

    // draw map over entire symmetry space
    float inv2Pi = 1.0f/(float)std::sqrt(2*M_PI);
    float invhd = 1.0f; //dimension of the kernel
    float mapSum = 0.0f, mapSumSquared = 0.0f;
    float mapCount = 0.0f;
    float averageInMask = 0.0f;
    float maskCount = 0.0f;

    logger("","WRITING MAP WITH SYMMETRY");
    for(int z=0; z<nc; z++){ // rectangular grid covers entire molecule i.e., all subunits
        float zsection = adjustedCMINZ+z*grid_spacing;
        for(int y=0; y<nb; y++){
            float ysection = adjustedCMINY+y*grid_spacing;
            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                const vector3 locale(adjustedCMINX+x*grid_spacing, ysection, zsection);
                float kernelSum = 0.0f;
                float incount = 0.0f;

                for(auto & point : lattice){
                    const vector3 * pVecPrimary = &(pModel->getBead(point.index)->getVec());
                    vector3 vec1 = *pVecPrimary - locale;
                    float length = vec1.length();
                    kernelSum += kernel(length/bandwidth)*point.prob;
                    for(unsigned int s = 1; s<totalSubunits; s++){ // add symmtery mates
                        vector3 vec2 = pModel->transformVectorBySymmetry(s, *pVecPrimary);
                        vec1 = vec2 - locale;
                        length = vec1.length();
                        kernelSum += kernel(length/bandwidth)*point.prob;
                        incount += 1.0f;
                    }
                    incount += 1.0f;
                }

                kernelSum *= (incount > 0) ? inv2Pi/incount*invhd : 0;
                kernel_mapping.emplace(key, kernelSum);
                /*
                 * need to determine if grid point is within the mask
                 * if so, determine average within the mask
                 */
//                for(auto mit : mask){
//                    Bead * pBead = pModel->getBead(mit);
//                    vector3 vec1 = pBead->getVec() - locale;
//                    if (vec1.length() < bandwidth){
//                        kernel_mask_mapping.emplace(key, kernelSum);
//                        averageInMask += kernelSum;
//                        maskCount += 1.0;
//                    }
//                }

                mapSum += kernelSum;
                mapSumSquared += kernelSum*kernelSum;
                mapCount += 1.0;
            }
        }
    }

    /*
     * ADJUST THE MAP SUCH THAT THE KERNEL THAT INTERSECTS THE MASK IS SCALE TO SOMETHING SENSIBLE ON AVERAGE
     * ESSENTIALLY EVERTHING OUTSIDE MASK STAYS THE SAME (LOW VALUE NEAR ZERO) WHEREAS EVERYTHING ELSE IS SCALED TO 0.414
     * SHOULD scale max value and then set average?
     */
    logger("","RESCALING MAP");
    //float scale = 0.414/(averageInMask/maskCount);
    float averageOriginalMap = mapSum/mapCount;
    float stdev = std::sqrt(mapSumSquared/mapCount - averageOriginalMap*averageOriginalMap);
    float scale = 1.0f/stdev;
    char buffer[80];

    for(int z=0; z<nc; z++){

        std::snprintf(buffer, 80, "%8i\n", z);
        tempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){
            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                auto allIt = kernel_mapping.find(key);
          //      auto maskIt = kernel_mask_mapping.find(key);
                /*
                 * write to file
                 */
                float kernelSum = (*allIt).second*scale;
                //    if (maskIt != kernel_mask_mapping.end()){
                //        kernelSum = (*maskIt).second*scale;
                //    }

                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                tempHeader.append(buffer);
                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    tempHeader += "\n";
                }
            }
        }

        if (stringIndex % 6 != 0){
            tempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    tempHeader.append(buffer);
    float ave = averageOriginalMap*scale;
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, 1.0);
    tempHeader.append(buffer);

    // write to file
    std::string map = outname + "_cemap_sym.xplor";
    const char * outputFileName = map.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);


    /*
     * flip the map x -> -x
     * reverse the order for each section
     */
    name = outname + "_cemask_flipped";
    pModel->writeModelToFileFlipped(mask.size(), mask, name, 0);

    for(int z=0; z<nc; z++){

        std::snprintf(buffer, 80, "%8i\n", z);
        mirrorTempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){
            for(int x=(na-1); x>=0; x--){ // reverse order

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                auto allIt = kernel_mapping.find(key);
                /*
                 * write to file
                 */
                float kernelSum = (*allIt).second*scale;

                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                mirrorTempHeader.append(buffer);
                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    mirrorTempHeader += "\n";
                }
            }
        }

        if (stringIndex % 6 != 0){
            mirrorTempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    mirrorTempHeader.append(buffer);
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, 1.0);
    mirrorTempHeader.append(buffer);

    // write to file
    map = outname+"_cemap_flipped.xplor";
    outputFileName = map.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, mirrorTempHeader.c_str());
    fclose(pFile);


    return 0;
}


void KDE::generateSeedFromInputFile(Model *pModel, Data *pData, std::string outname){
    float cutoff = bead_radius + 0.5*dmin_supremum; // sphere of influence, looking for all vectors in centered coordnates that are within...
    //cutoff = std::sqrt(3.0d/2.0d)*average_dmin*1.02; // midpoint of a tetrahedron
    //grid_spacing = 2.0/3.0*bead_radius; // spacing in rectangular grid for calculating map
    grid_spacing = average_dmin; // spacing in rectangular grid for calculating map
    bandwidth = grid_spacing*std::sqrt(3.0d); // diagonal of a square 1.73*r

    std::string prefix = outname+"_seed";
    logger("REMAPPING TO NEW LATTICE",outname);
    remapToNewHCPLattice(pModel, cutoff);
    std::vector<unsigned int> mask; // contains indices to write out with respect to Model (Universe)
    // make a kept list that is essentially contoured to >= average
    createMaskOnCenteredObject(mask, 1);

    EulerTour eulerTour (mask, mask.size(), pModel);
    if (eulerTour.getNumberOfComponents() > 1){ // prune the mask
        std::cout << "               EULER TOUR SIZE : TOO MANY TOURS -- ADJUSTING "  <<std::endl;
        const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
        const unsigned int neighborLimit = pModel->getNeighborLimit();
        const unsigned int * ptr = pModel->getDirectPointerToNeighborhood();

        auto pLargest = eulerTour.getPointerToLargestTour(); // assuming this is by far the majority of the points
        auto pETours = eulerTour.getTours();


        std::set<unsigned int> largeTour;

        for(auto it = pLargest->begin(); it != pLargest->end(); ++it){
            largeTour.insert((*it)->getKey());
        }
        std::set<unsigned int> beads_in_use(mask.begin(), mask.end());


        for(auto & tour : *pETours){
            if (&tour.second != pLargest){
                // try finding a common neighbor betweeen the largest and current tours
                std::set<unsigned int> smallTour;
                for(auto it = tour.second.begin(); it != tour.second.end(); ++it){
                    smallTour.insert((*it)->getKey());
                }

                for(auto index : smallTour){

                    for (unsigned int i=0; i< totalNeighbors; i++){
                        //unsigned int neighbor = *(it+i);
                        unsigned int neighbor = ptr[totalNeighbors*index + i]; // is this neighbor

                        if (neighbor < neighborLimit && smallTour.find(neighbor) == smallTour.end()){ // not in small tour
                            // check if a neighbor of largeTour
                            for(auto & lindex : largeTour){
                                for (unsigned int j=0; j< totalNeighbors; j++) {
                                    unsigned int largeNeighbor = ptr[totalNeighbors * lindex + j]; // is this neighbor
                                    if (largeNeighbor == neighbor){
                                        beads_in_use.insert(neighbor); // add neighbor that is in common
                                        eulerTour.addNode(neighbor, pModel);
                                        if (eulerTour.getNumberOfComponents() == 1){
                                            goto alldone;
                                        }
                                        goto next_neighbor;
                                    } else if (largeNeighbor == neighborLimit){
                                        break;
                                    }
                                }
                            }
                        } else if (neighbor == neighborLimit){
                            next_neighbor:
                            break;
                        }
                    }
                }
            }
        }

        alldone:
        mask.clear();
        for(auto ind : beads_in_use){
            mask.emplace_back(ind);
        }
    }

    logger("EULER TOUR SIZE", std::to_string(eulerTour.getNumberOfComponents()));
    logger("MASK SIZE", std::to_string(mask.size()));
    logger("OUTFILE NAME", prefix);

    pModel->writeModelToFile(mask.size(), mask, prefix, 0);

}

/**
 * take map and refined model and resample map
 * if using refined SEED to bias the map, the seed must be aligned to the first file in the kde.inp
 *
 * @return
 */
float KDE::map_refine(Model *pModel, Data *pData, std::string outname) {

    float cutoff = bead_radius + 0.5*average_dmin; // sphere of influence, looking for all vectors in centered coordnates that are within...
    cutoff = std::sqrt(3.0d/2.0d)*bead_radius*1.02; // midpoint of a tetrahedron
    //grid_spacing = 2.0/3.0*bead_radius; // spacing in rectangular grid for calculating map
    grid_spacing = average_dmin; // spacing in rectangular grid for calculating map
    bandwidth = grid_spacing*std::sqrt(3.0) > (std::sqrt(3.0d/2.0d)*bead_radius) ? std::sqrt(3.0d/2.0d)*bead_radius : grid_spacing*std::sqrt(3.0);


    logger("DATA LATTICE SPACING", formatNumber(bead_radius,2));
    logger("CUTOFF", formatNumber(cutoff,3));
    logger("GRID SPACING", formatNumber(grid_spacing,2));
    logger("BANDWIDTH", formatNumber(bandwidth,2));

    unsigned int highTempRounds = 31;
    std::string prefix = "CE";

    float limit = remapToNewHCPLattice(pModel, cutoff); // only map the pModel beads that intersect input models
    // create vector of indices and a vector to be cdf

    std::vector<ProbabilityBead> lattice; // related to pModel
    std::vector<vector3> remapped;
    /*
     * probabilities are calculated as the number of beads overlapping with new lattice
     * We can't treat each bead equally, as some may barely intersect.
     */
    for(auto & index : bead_to_vec_map){
        if (index.second.size() >= limit){
            lattice.emplace_back(ProbabilityBead{0.67, index.first, 0.67});
            remapped.emplace_back(pModel->getBead(index.first)->getVec());
        } else if (index.second.size() < limit && index.second.size() >= (limit - 0.5*stdev)) {
            lattice.emplace_back(ProbabilityBead{0.1, index.first, 0.1});
        } else if (index.second.size() < (limit - 0.5*stdev)){
            lattice.emplace_back(ProbabilityBead{0.02, index.first, 0.02});
        }
    }

    writeCoordinatesToFile("filtered", remapped);

    unsigned int pot = (double)estimateRounds(lattice.size()/2,lattice.size())*87;

    trialSize = (pot < 37571) ? 37571 : pot;
//    trialSize = pot;
//    trialSize = (unsigned int)((double)estimateRounds(minN,lattice.size())/topNpercent);
    topNSize = (unsigned int)(topNpercent*trialSize);

    Anneal mainAnneal(0,
                      0,
                      highTempRounds,
                      prefix,
                      0.03, //alpha
                      0.001,//beta
                      0,    //eta
                      0.01, //lambda
                      0,    //mu
                      trialSize,
                      0,
                      2*bead_radius
    );

//    std::vector<unsigned int> mask; // contains indices to write out with respect to Model (Universe)
//    // make a kept list that is essentially contoured to >= average
//    createMaskOnCenteredObject(mask, 1);
//    std::string masknem = "entiremaskfromlattice_seed";
//    pModel->writeModelToFile(mask.size(), mask, masknem, 0);

    if (priorSet){
        std::cout << "       UPDATING MAP WITH PRIOR : " << prior_name << std::endl;
        double sum=0.0d;
        double alpha = 0.7;
        std::vector<double> temp;
        std::set<unsigned int> indicesToUpdate;
        for(auto & lat : lattice){ //determine which lattice points overlap with prior

            double kernelSum = 0.0d;
            double minL = 1000.0d;
            unsigned int indexToUpdate=0;
            for(unsigned int i=0; i<centered_prior.size(); i++){
                vector3 lVec = (centered_prior[i] - pModel->getBead(lat.index)->getVec());
                kernelSum += kernel(lVec.length()/cutoff);
                if (lVec.length() < minL){
                    minL = lVec.length();
                    indexToUpdate = lat.index;
                }
            }

            if (minL < cutoff){
                indicesToUpdate.insert(indexToUpdate);
            }
            sum+=kernelSum;
            temp.emplace_back(kernelSum);
        }


        double inv = 1.0d/sum;
//        double invC = 1.0d/(double)indicesToUpdate.size();
        for(unsigned int i=0; i<lattice.size(); i++){ // update lattice probabilities
            ProbabilityBead * pPoint = &lattice[i];
            double old = pPoint->prob;
            pPoint->prob = old*(1.0-alpha) + alpha*temp[i]*inv;
//            if (indicesToUpdate.find(pPoint->index) != indicesToUpdate.end()){
//                pPoint->prob = old*(1.0-alpha) + alpha*1.0d*invC;
//            } else {
//                pPoint->prob = old*(1.0-alpha);
//            } // probabilities should sum to one
        }
    }

//    writeCoordinatesToFile("remapped", remapped);
    unsigned int lowerBound = minN*mapAverage;
    auto rescaledMaxN = (unsigned int)(maxN*mapAverage);
    if (rescaledMaxN > bead_to_vec_map.size()){
        rescaledMaxN = bead_to_vec_map.size();
    }
    mainAnneal.ceMapOptimization(pModel, pData, topNSize, lowerBound, rescaledMaxN, lattice);

    /*
     * determine average and variance of the refined lattice
     */
    unsigned int totalInLattice = lattice.size();
    double sumIt=0.0d;
    double var=0.0d;
    for(auto & point : lattice){
        sumIt+=point.prob;
        var += point.prob*point.prob;
    }
    double average = sumIt/(double)totalInLattice;
    double variance = var/(double)totalInLattice - average*average;
    logger("MAP AVERAGE", formatNumber(average, 2));
    logger("MAP SIGMA", formatNumber(std::sqrt(variance), 2));

    std::vector<unsigned int> mask;
    for(auto & point : lattice){
        if (point.prob >= (average)){
            mask.emplace_back(point.index);
        }
    }
    // make a kept list that is essentially contoured to >= average
    //createMaskOnCenteredObject(mask, limit);
    std::string name = outname + "_cemask";
    pModel->writeModelToFile(mask.size(), mask, name, 0);

    int na, nb, nc;
    float adjustedCMINX, adjustedCMINY, adjustedCMINZ;
    std::string tempHeader = headerParametersForXPLOR(na, nb, nc, adjustedCMINX, adjustedCMINY, adjustedCMINZ);


    float inv2Pi = 1.0f/(float)std::sqrt(2*M_PI);
    float invhd = 1.0f; //dimension of the kernel
    float mapSum = 0.0f, mapSumSquared = 0.0f;
    float mapCount = 0.0f;
    float averageInMask = 0.0f;
    float maskCount = 0.0f;

    for(int z=0; z<nc; z++){ // rectangular grid so everything is incremented by bead_radius
        float zsection = adjustedCMINZ+z*grid_spacing;
        for(int y=0; y<nb; y++){
            float ysection = adjustedCMINY+y*grid_spacing;
            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                const vector3 locale(adjustedCMINX+x*grid_spacing, ysection, zsection);
                float kernelSum = 0.0f;
                float incount = 0.0f;

                for(auto & point : lattice){
                    vector3 vec1 = pModel->getBead(point.index)->getVec() - locale;
                    float length = vec1.length();
                    kernelSum += kernel(length/bandwidth)*point.prob;
                    incount += 1.0f;
                }

                kernelSum *= (incount > 0) ? inv2Pi/incount*invhd : 0;
                kernel_mapping.emplace(key, kernelSum);
                /*
                 * need to determine if grid point is within the mask
                 * if so, determine average within the mask
                 *
                 */
                for(auto mit : mask){
                    Bead * pBead = pModel->getBead(mit);
                    vector3 vec1 = pBead->getVec() - locale;
                    if (vec1.length() < bandwidth){
                        kernel_mask_mapping.emplace(key, kernelSum);
                        averageInMask += kernelSum;
                        maskCount += 1.0;
                    }
                }

                mapSum += kernelSum;
                mapSumSquared += kernelSum*kernelSum;
                mapCount += 1.0;
            }
        }
    }
    /*
     * ADJUST THE MAP SUCH THAT THE KERNEL THAT INTERSECTS THE MASK IS SCALE TO SOMETHING SENSIBLE ON AVERAGE
     * ESSENTIALLY EVERTHING OUTSIDE MASK STAYS THE SAME (LOW VALUE NEAR ZERO) WHEREAS EVERYTHING ELSE IS SCALED TO 0.414
     * SHOULD scale max value and then set average?
     */
    //float scale = 0.414/(averageInMask/maskCount);
    float averageOriginalMap = mapSum/mapCount;
    float stdev = std::sqrt(mapSumSquared/mapCount - averageOriginalMap*averageOriginalMap);
    float scale = 1.0f/stdev;
    char buffer[80];

    for(int z=0; z<nc; z++){

        std::snprintf(buffer, 80, "%8i\n", z);
        tempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){
            for(int x=0; x<na; x++){

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                auto allIt = kernel_mapping.find(key);
           //     auto maskIt = kernel_mask_mapping.find(key);
                /*
                 * write to file
                 */
                float kernelSum = (*allIt).second*scale;
            //    if (maskIt != kernel_mask_mapping.end()){
            //        kernelSum = (*maskIt).second*scale;
            //    }
                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                tempHeader.append(buffer);
                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    tempHeader += "\n";
                }
            }
        }

        if (stringIndex % 6 != 0){
            tempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    tempHeader.append(buffer);
    float ave = averageOriginalMap*scale;
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, 1.0);
    tempHeader.append(buffer);

    // write to file
    std::string map = outname+"_cemap.xplor";
    const char * outputFileName = map.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);

    /*
     * flip the map x -> -x
     * reverse the order for each section
     */
    name = outname + "_cemask_flipped";
    pModel->writeModelToFileFlipped(mask.size(), mask, name, 0);

    for(int z=0; z<nc; z++){

        std::snprintf(buffer, 80, "%8i\n", z);
        mirrorTempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){
            for(int x=(na-1); x>=0; x--){ // reverse order

                std::string key = std::to_string(x)+"-"+std::to_string(y)+"-"+std::to_string(z);

                auto allIt = kernel_mapping.find(key);
                /*
                 * write to file
                 */
                float kernelSum = (*allIt).second*scale;

                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                mirrorTempHeader.append(buffer);
                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    mirrorTempHeader += "\n";
                }
            }
        }

        if (stringIndex % 6 != 0){
            mirrorTempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    mirrorTempHeader.append(buffer);
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, 1.0);
    mirrorTempHeader.append(buffer);

    // write to file
    map = outname+"_cemap_flipped.xplor";
    outputFileName = map.c_str();

    pFile = fopen(outputFileName, "w");
    fprintf(pFile, mirrorTempHeader.c_str());
    fclose(pFile);


    return 0;
}

void KDE::createMaskOnCenteredObject(std::vector<unsigned int> & mask, float limit){
    // make a kept list that is essentially contoured to >= average
    for(auto mit : bead_to_vec_map){
        float temp = mit.second.size();
        if (temp >= limit){
            mask.push_back(mit.first);
        }
    }
}

std::string KDE::headerParametersForXPLOR(int & na, int & nb, int & nc, float & adjustedCMINX, float & adjustedCMINY, float & adjustedCMINZ){
    std::cout << " BOUNDING BOX "<< std::endl;
    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    printf("  => X %8.4f %8.4f %8.4f \n", cminx, cmaxx, (cmaxx-cminx));
    printf("  => Y %8.4f %8.4f %8.4f \n", cminy, cmaxy, (cmaxy-cminy));
    printf("  => Z %8.4f %8.4f %8.4f \n", cminz, cmaxz, (cmaxz-cminz));

    auto startingNA = (int)std::round(cminx/grid_spacing);
    auto startingNB = (int)std::round(cminy/grid_spacing);
    auto startingNC = (int)std::round(cminz/grid_spacing);

    auto stoppingNA = (int)std::round(cmaxx/grid_spacing);
    auto stoppingNB = (int)std::round(cmaxy/grid_spacing);
    auto stoppingNC = (int)std::round(cmaxz/grid_spacing);

    adjustedCMINX = startingNA*grid_spacing;
    adjustedCMINY = startingNB*grid_spacing;
    adjustedCMINZ = startingNC*grid_spacing;

    float adjustedCMAXX = stoppingNA*grid_spacing;
    float adjustedCMAXY = stoppingNB*grid_spacing;
    float adjustedCMAXZ = stoppingNC*grid_spacing;

    float a_side = adjustedCMAXX - adjustedCMINX;
    float b_side = adjustedCMAXY - adjustedCMINY;
    float c_side = adjustedCMAXZ - adjustedCMINZ;

    std::cout << " ADJUSTED BOUNDING BOX "<< std::endl;
    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    printf("  => X %8.4f %8.4f %8.4f \n", adjustedCMINX, adjustedCMAXX, a_side);
    printf("  => Y %8.4f %8.4f %8.4f \n", adjustedCMINY, adjustedCMAXY, b_side);
    printf("  => Z %8.4f %8.4f %8.4f \n", adjustedCMINZ, adjustedCMAXZ, c_side);

    na = std::abs(startingNA) + std::abs(stoppingNA) + 1;
    nb = std::abs(startingNB) + std::abs(stoppingNB) + 1;
    nc = std::abs(startingNC) + std::abs(stoppingNC) + 1;

    mirrorTempHeader = "\n";
    std::string tempHeader = "\n";
    tempHeader += "        4 !NTITLE\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";

    mirrorTempHeader += tempHeader;
    char buffer[80];
    std::snprintf(buffer, 80, "%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", na, startingNA, stoppingNA, nb, startingNB, stoppingNB, nc, startingNC, stoppingNC);
    tempHeader.append(buffer);

    std::snprintf(buffer, 80, "%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", na, -stoppingNA, -startingNA, nb, startingNB, stoppingNB, nc, startingNC, stoppingNC);
    mirrorTempHeader.append(buffer);

    std::snprintf(buffer, 80, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E \n",a_side, b_side, c_side, 90.0, 90.0, 90.0);
    tempHeader.append(buffer);
    mirrorTempHeader.append(buffer);

    // write the matrix
    tempHeader += "ZYX\n";
    mirrorTempHeader += "ZYX\n";
    return tempHeader;
}


std::string KDE::headerParametersForXPLORSYM(int & na, int & nb, int & nc, float & adjustedCMINX, float & adjustedCMINY, float & adjustedCMINZ, std::string sym){

    auto startingNA = (int)std::round(minx/grid_spacing);
    auto startingNB = (int)std::round(miny/grid_spacing);
    auto startingNC = (int)std::round(minz/grid_spacing);

    auto stoppingNA = (int)std::round(maxx/grid_spacing);
    auto stoppingNB = (int)std::round(maxy/grid_spacing);
    auto stoppingNC = (int)std::round(maxz/grid_spacing);

    adjustedCMINX = startingNA*grid_spacing;
    adjustedCMINY = startingNB*grid_spacing;
    adjustedCMINZ = startingNC*grid_spacing;

    float adjustedCMAXX = stoppingNA*grid_spacing;
    float adjustedCMAXY = stoppingNB*grid_spacing;
    float adjustedCMAXZ = stoppingNC*grid_spacing;

    float a_side = adjustedCMAXX - adjustedCMINX;
    float b_side = adjustedCMAXY - adjustedCMINY;
    float c_side = adjustedCMAXZ - adjustedCMINZ;

    logger("ADJUSTED BOUNDING BOX", "");
    //std::cout << " ADJUSTED BOUNDING BOX "<< std::endl;
    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    printf("  => X %8.4f %8.4f %8.4f \n", adjustedCMINX, adjustedCMAXX, a_side);
    printf("  => Y %8.4f %8.4f %8.4f \n", adjustedCMINY, adjustedCMAXY, b_side);
    printf("  => Z %8.4f %8.4f %8.4f \n", adjustedCMINZ, adjustedCMAXZ, c_side);

    na = std::abs(startingNA) + std::abs(stoppingNA) + 1;
    nb = std::abs(startingNB) + std::abs(stoppingNB) + 1;
    nc = std::abs(startingNC) + std::abs(stoppingNC) + 1;

    std::string tempHeader = "\n";
    tempHeader += "        4 !NTITLE\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265        SYMMETRY : "+sym+"\n";
    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";
    char buffer[80];
    std::snprintf(buffer, 80, "%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", na, startingNA, stoppingNA, nb, startingNB, stoppingNB, nc, startingNC, stoppingNC);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E \n",a_side, b_side, c_side, 90.0, 90.0, 90.0);
    tempHeader.append(buffer);
    // write the matrix
    tempHeader += "ZYX\n";

    return tempHeader;
}

/*
 * Prior should be aligned to the base model used to make the averaged set
 * In future, will need to auto-align to first file in the kde.inp list
 *
 */
void KDE::add_prior(std::string file) {

    logger("ADDING PRIOR MODEL", file);
    PDBModel tempmodel (file, false, true, bead_radius);
    unsigned int totalCoords = tempmodel.getTotalAtoms();

    for(unsigned int i = 0; i<totalCoords; i++){
        float pX = tempmodel.getXCoord(i);
        float pY = tempmodel.getYCoord(i);
        float pZ = tempmodel.getZCoord(i);
        prior.emplace_back(vector3(pX, pY, pZ));
        centered_prior.emplace_back(vector3(pX-centering_vec.x, pY-centering_vec.y, pZ-centering_vec.z));
    }
    priorSet = true;
    prior_name = file;
    writeCoordinatesToFile("centered_seed", centered_prior);

    /*
     * flip the model
     */
    for(auto & coords : centered_prior){
        coords.x = coords.x*(-1.0f);
    }
    writeCoordinatesToFile("centered_seed_flipped", centered_prior);

    /*
     * return to original
     */
    for(auto & coords : centered_prior){
        coords.x = coords.x*(-1.0f);
    }

}


void KDE::logger(std::string description, std::string value) {

    unsigned int len = 40 - description.size();
    std::string firsthalf = std::string(len, ' ');
    printf("%s%s : %s\n", firsthalf.c_str(), description.c_str(), value.c_str());

}


std::string KDE::formatNumber(float number, int decimals = 2) {
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


unsigned int KDE::estimateRounds(unsigned int lowerBound, unsigned int totalInSet){

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<unsigned int> all(totalInSet);
    for(unsigned int i=0; i<totalInSet; i++){
        all[i] = i;
    }

    std::set<unsigned int> selected;

    unsigned int round=0;
    while (selected.size() < totalInSet){
        std::shuffle(all.begin(), all.end(), gen);
        for(unsigned int i=0; i<lowerBound; i++){
            selected.insert(all[i]);
        }
        round += 1;
    }

    return round;
}