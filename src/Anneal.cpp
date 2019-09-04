//
// Created by xos81802 on 10/07/2018.
//

#include "Anneal.h"

Anneal::Anneal(float highT,
               float percent,
               unsigned int highTempRounds,
               std::string fileprefix,
               float alpha,
               float beta,
               float eta,
               float lambda,
               float mu,
               unsigned int multiple,
               float accRate,
               float intercon) : highT(highT),
                                percentAddRemove(percent),
                                highTempRounds(highTempRounds),
                                filenameprefix(fileprefix),
                                alpha(alpha),
                                beta(beta),
                                eta(eta),
                                lambda(lambda),
                                mu(mu),
                                ccmultiple(multiple),
                                asaAcceptanceRate (accRate),
                                interconnectivityCutOff(intercon)
{


    //this->highTempStartForCooling = 0.00001; //0.00001
    this->highTempStartForCooling = 0.000016; //0.00001
    this->asaAcceptanceRate = accRate;
    complementASAAcceptanceRate = 1.0f - accRate;
    intASAAcceptanceRate = (int)1000*accRate;
    intComplementASAAcceptanceRate = (int)1000*complementASAAcceptanceRate;



    contactsDistribution.resize(13);
    /*
     * this distribution is set to q-max 0.4 and bead_radius = 0.5*bin_width
     */
    distributionlimit=13;
    contactsDistribution[0] = 0.0d;
    contactsDistribution[1] = 0.25419d;
    contactsDistribution[2] = 0.317564d;
    contactsDistribution[3] = 0.241696d;
    contactsDistribution[4] = 0.117231d;
    contactsDistribution[5] = 0.0447598d;
    contactsDistribution[6] = 0.0190639d;
    contactsDistribution[7] = 0.00549451d;
    contactsDistribution[8] = 0.0d;
    contactsDistribution[9] = 0.0d;
    contactsDistribution[10] = 0.0d;
    contactsDistribution[11] = 0.0d;
    contactsDistribution[12] = 0.0d;

    double totaltemp =0;
    for(unsigned int i=0; i<distributionlimit; i++){
        totaltemp += contactsDistribution[i];
    }

    for(unsigned int i=0; i<distributionlimit; i++){ // normalize
        contactsDistribution[i] *= 1.0/totaltemp;
    }
}

/**
 *
 * @param flags
 * @param bead_indices
 * @param upTo is the workingLimit of beads_indices
 * @param pModel
 * @return
 */
float Anneal::calculateCVXHULLVolume(char *flags, std::vector<unsigned int> *bead_indices, const unsigned int upTo, Model *pModel) {

    unsigned int numpoints = 3*upTo;
    coordT points[numpoints];

    const unsigned int * ptr = (*bead_indices).data();

    for (unsigned int i=0; i<upTo; i++){
        beadToPoint(&points[i*3], pModel->getBead(ptr[i]));
        //  beadToPoint(&(pointSet[i*3]), pModel->getBead((*bead_indices)[i]));
    }


    qh_new_qhull (3, upTo, points, 0, flags, nullptr, nullptr);
    float volume_test = (float)(qh totvol);

    //qh totarea;
    qh_freeqhull(true);
    //float calcVol = pModel->getBeadVolume()*upTo;
    return volume_test;///calcVol;
}

/**
 *
 * @param flags
 * @param bead_indices
 * @param upTo is the workingLimit of beads_indices
 * @param pModel
 * @return
 */
float Anneal::calculateCVXHULLVolumeSet(char *flags, std::set<unsigned int> *bead_indices_tree, const unsigned int upTo, Model *pModel) {

    unsigned int numpoints = 3*upTo;
    coordT points[numpoints];

    unsigned int i=0;
    for (const auto & bead : *bead_indices_tree) {
        beadToPoint(&points[i*3], pModel->getBead(bead));
        i++;
    }

    qh_new_qhull (3, upTo, points, 0, flags, nullptr, nullptr);
    float volume_test = (float)(qh totvol);

    //qh totarea;
    qh_freeqhull(true);
    //float calcVol = pModel->getBeadVolume()*upTo;
    return volume_test;///calcVol;
}


float Anneal::changeInSurfaceToVolumeAdd(unsigned int beadToAdd, std::set<unsigned int> *beads_in_use, Model * pModel  ){

    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(beadToAdd);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    auto endOfSet = beads_in_use->end();
    unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    std::vector<Eigen::Vector3f> coordinates;
    coordinates.reserve(pModel->getSizeOfNeighborhood());

    unsigned int count=0;
    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != endOfSet){ // if end of set, means not in use
            vector3 pBeadVec = pModel->getBead(neighbor)->getVec();
            coordinates.emplace_back(Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z));
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    std::vector<float> weights(count);
    float probe = pModel->getBeadRadius() + delta_r;
    std::fill(weights.begin(), weights.end(), probe);

    POWERSASA::PowerSasa<float, Eigen::Vector3f> *ps =
            new POWERSASA::PowerSasa<float, Eigen::Vector3f>(coordinates, weights, 1, 1, 1, 1);

    ps->calc_sasa_all();
    float volume = 0.0;
    float sasa = 0.0;
    for (unsigned int i = 0; i < count; ++i) {
        sasa += ps->getSasa()[i];
        volume += ps->getVol()[i];
    }

    // collect beads that are in contact
    //delete ps; //maximizes surface area while minimizing number of regions

    auto currentStoV = sasa;//(float)(36.0f*M_PI*(volume*volume)/(sasa*sasa*sasa));
    /*
     * Add new bead and recalculate
     */
    vector3 pBeadVec = pModel->getBead(beadToAdd)->getVec();
    coordinates.emplace_back(Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z));
    weights.push_back(probe);

    ps->update_coords(coordinates, weights);
    ps->calc_sasa_all();
    volume = 0.0;
    sasa = 0.0;

    count++;
    for (unsigned int i = 0; i < count; ++i) {
        sasa += ps->getSasa()[i];
        volume += ps->getVol()[i];
    }

    // collect beads that are in contact
    delete ps; //maximizes surface area while minimizing number of regions
    return (currentStoV - sasa); //sphericity : minimizing moves away from a sphere
    // lots of surface area supports an extended model
    // minimal surface area supports a compact object
    //return (sasa/(float)total)/surfaceArea;
};


float Anneal::surfaceToVolume(const unsigned int total, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates){

    POWERSASA::PowerSasa<float, Eigen::Vector3f> *ps =
            new POWERSASA::PowerSasa<float, Eigen::Vector3f>(coordinates, weights, 1, 1, 1, 1);

    ps->calc_sasa_all();
    float volume = 0.0;
    float sasa = 0.0;
    for (unsigned int i = 0; i < total; ++i) {
        sasa += ps->getSasa()[i];
        volume += ps->getVol()[i];
    }

    delete ps; //maximizes surface area while minimizing number of regions
    return (36*M_PI*(volume*volume)/(sasa*sasa*sasa)); //sphericity : minimizing moves away from a sphere
    // lots of surface area supports an extended model
    // minimal surface area supports a compact object
    //return sasa;
};



/**
 * Creates a distribution of contacts from selected lattice points specified by beads_in_use
 *
 * @param distribution
 * @param beads_in_use
 * @param pModel
 */
void Anneal::populateContactsDistribution(std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel){

    for(int i=0; i<13; i++){
        distribution[i]=0.0d;
    }

    for (auto const & it : *beads_in_use){
        ++distribution[ numberOfContactsFromSet(beads_in_use, pModel, it) ];
    }
}


double Anneal::removeFromContactsPotential(unsigned int removeMe, double currentContactsSum, std::set<unsigned int> *beads_in_use, Model *pModel){

    double currentSum = currentContactsSum;
    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(removeMe);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    double oldPot = 0.0d;
    double newPot = 0.0d;

    for(unsigned int i=0; i<count; i++){
        unsigned int tempNumber=numberOfContactsFromSet(beads_in_use, pModel, neighborsInuse[i]);
        oldPot += contactPotentialFunction(tempNumber);
        newPot += contactPotentialFunction(tempNumber-1);
    }

    currentSum += newPot - oldPot - contactPotentialFunction(numberOfContactsFromSet(beads_in_use, pModel, removeMe));

    return currentSum;
}

double Anneal::addToContactsPotential(unsigned int addMe, double currentContactsSum, std::set<unsigned int> *beads_in_use, Model *pModel){

    double currentSum = currentContactsSum;
    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(addMe);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    double delta = 0;
    for(unsigned int i=0; i<count; i++){
        unsigned int tempNumber=numberOfContactsFromSet(beads_in_use, pModel, neighborsInuse[i]);
        double oldPot = contactPotentialFunction(tempNumber-1);
        double newPot = contactPotentialFunction(tempNumber);
        delta += (newPot - oldPot);
    }

    currentSum += contactPotentialFunction(numberOfContactsFromSet(beads_in_use, pModel, addMe)) + delta;

    return currentSum;
}

inline double Anneal::contactPotentialFunction(unsigned int nc){
    double diff, value = 0.0d;
    double target = 3.0d;
   // double slope = -0.0011;//(0.0001d - 0.01d)/(12.0d-target);
   // double intercept = 0.0133;//0.01 - slope*target;
    double slope = 0.011;   // minimum at 3, rises to 0.01 at nc of 12
    double intercept = -0.032;// minimum at 3, rises to 0.01 at nc of 12

    // double slope = -0.0011;//(0.0001d - 0.01d)/(12.0d-target);
    // double intercept = 0.0133;//0.01 - slope*target;

    if (nc == 1){
        diff = (double)nc-target;
        value += 13.0d*diff*diff + 0.01d;
    } else if (nc < target){
        diff = (double)nc-target;
        value += 3.0d*diff*diff + 0.01d;
    } else {
        value += -0.00110d*(double)nc + 0.01330d;
    }

//
//    if (nc == 1){
//        diff =(double)nc-target;
//        value += 37.0d*diff*diff + 0.001d;
//    } else if (nc < target) {
//        diff =(double)nc-target;
//        value += 7*diff*diff + 0.001d;
//    } else {
//        value += 0.001;
//    }
    return value;
}

double Anneal::contactPotential(std::set<unsigned int> *beads_in_use, Model *pModel){
    double value = 0.0d;//, inv = 1.0d/(double)beads_in_use->size();

    for (auto it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {
        unsigned int nc = numberOfContactsFromSet(beads_in_use, pModel, *it);
        value += contactPotentialFunction(nc);
    }

    return value;//*inv;
}


/**
 * beads_in_use must already contain addMe
 * @param addMe
 * @param distribution
 * @param beads_in_use
 * @param pModel
 */
void Anneal::addToContactsDistribution(unsigned int addMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel){
    // addMe will have n-number of contacts and for each
    // get all neighbors for addMe
    double * const pDistribution = distribution.data(); // initialized as empty in Model class
    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(addMe);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // now adjust counts
    for(unsigned int i=0; i<count; i++){
        unsigned int tempNumber=numberOfContactsFromSet(beads_in_use, pModel, neighborsInuse[i]);
        --pDistribution[ tempNumber-1 ];
        ++pDistribution[ tempNumber ];
//        --distribution[ tempNumber-1 ]; // remove prior contribution
//        ++distribution[ tempNumber ];    // add new contribution
    }
    // now do addMe
    ++pDistribution[numberOfContactsFromSet(beads_in_use, pModel, addMe)];
}

/**
 * beads_in_use must already contain removeMe
 * @param addMe
 * @param distribution
 * @param beads_in_use
 * @param pModel
 */
void Anneal::removeFromContactsDistribution(unsigned int removeMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel){
    // get all neighbors for removeMe
    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(removeMe);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // now adjust counts

    for(unsigned int i=0; i<count; i++){
        unsigned int currentNumber=numberOfContactsFromSet(beads_in_use, pModel, neighborsInuse[i]);
        --distribution[ currentNumber ];      // remove prior contribution
        ++distribution[ currentNumber-1 ];    // add new contribution
    }

    // now remove contributions from removeMe bead
    --distribution[ numberOfContactsFromSet(beads_in_use, pModel, removeMe)];
}

/*!
 *
 */
double Anneal::connectivityPotential(unsigned int numberOfComponents){

    switch(numberOfComponents) {
        case 1:
            return 0.0d;
        case 2:
            return 0.1d;
        case 3:
            return 0.30d;
        case 4:
            return 0.70d;
        case 5:
            return 100000.0d;
        case 6:
            return 1000000.0d;
        default:
            return 1000000.0d*numberOfComponents;
    }
}

/*
* add new position as last element in each container
* reverse is a simple pop-off
*/
void Anneal::addForSurfaceToVolume(unsigned int toAdd, std::vector<unsigned int> & unsortedIndices, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates, Model * pModel){
    unsortedIndices.push_back(toAdd);
    weights.push_back(pModel->getBeadRadius() + delta_r);
    const vector3 & pVec = pModel->getBead(toAdd)->getVec();
    coordinates.push_back(Eigen::Vector3f(pVec.x, pVec.y, pVec.z));
}


void Anneal::removeForSurfaceToVolume(unsigned int toRemove, unsigned int workingLimit, std::vector<unsigned int> & unsortedIndices, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates){
    unsigned int last = workingLimit - 1;
    auto it = std::find(unsortedIndices.begin(), unsortedIndices.end(), toRemove);
    unsigned int locale = std::distance(unsortedIndices.begin(), it);

    std::iter_swap(unsortedIndices.begin() + locale, unsortedIndices.begin() + last);
    std::iter_swap(weights.begin() + locale, weights.begin() + last);
    std::iter_swap(coordinates.begin() + locale, coordinates.begin() + last);
    unsortedIndices.pop_back();
    weights.pop_back();
    coordinates.pop_back();
}


void Anneal::swapForSurfaceToVolume(unsigned int indexToFind, unsigned int indexToAdd, std::vector<unsigned int> & unsortedIndices, std::vector<Eigen::Vector3f> & coordinates, Model * pModel){
    auto it = std::find(unsortedIndices.begin(), unsortedIndices.end(), indexToFind);
    unsigned int locale = std::distance(unsortedIndices.begin(), it);
    *it = indexToAdd;

    const vector3 & pVec = pModel->getBead(indexToAdd)->getVec();
    coordinates[locale] = Eigen::Vector3f(pVec.x, pVec.y, pVec.z);
}


/**
 * from Vincent A. Cicirello
 * On the Design of an Adaptive Simulated Annealing Algorithm
 *
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void Anneal::updateASATemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

//    bool changed = false;
    double stepEval = (double)index/evalMax;
    double lamRate=asaAcceptanceRate;

    if (stepEval < 0.15) {
        //lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
        lamRate = asaAcceptanceRate+complementASAAcceptanceRate*pow(intComplementASAAcceptanceRate, -stepEval*6.666667);
    } else if (stepEval >= 0.15 && stepEval < 0.65){
        lamRate = asaAcceptanceRate;
    } else if (stepEval >= 0.65){
        lamRate = asaAcceptanceRate*pow(intASAAcceptanceRate, -(stepEval - 0.65)*2.857142857);
    }

    if (acceptRate > lamRate){
        temp = 0.999*temp;
//        changed=true;
    } else {
        temp = temp*1.001001001001;
//        changed=true;
    }

//    if (changed){
        inv_temp = 1.0/temp;
//    }
}


/**
 * from Vincent A. Cicirello
 * On the Design of an Adaptive Simulated Annealing Algorithm
 *
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void Anneal::updateASAConstantTemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp){


    double stepEval = (double)index/evalMax;
    double lamRate=asaAcceptanceRate;

    if (stepEval < 0.15) {
        //lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
        lamRate = asaAcceptanceRate+complementASAAcceptanceRate*pow(intComplementASAAcceptanceRate, -stepEval*6.666667);
    } else if (stepEval >= 0.15){
        lamRate = asaAcceptanceRate;
    }

    if (acceptRate > lamRate){
        temp = 0.999*temp;
    } else {
        temp = temp*1.001001001001;
    }

    inv_temp = 1.0/temp;
}


void Anneal::updateASAModTemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

//    bool changed = false;
    double stepEval = (double)index/evalMax;
    double lamRate=0.27;
    double finalEval = 0.65;

    if (stepEval < 0.10) { // 80% acceptance rate
        lamRate = 0.80+0.20*pow(200, -stepEval*6.666667); // exponential increase to get to
    } else if (stepEval >= 0.10 && stepEval < 0.20){
        lamRate = 0.70; // 70% acceptance rate
    } else if (stepEval >= 0.20 && stepEval < 0.30){
        lamRate = 0.60; // 60% acceptance rate
    } else if (stepEval >= 0.30 && stepEval < 0.40){
        lamRate = 0.50; // 50% acceptance rate
    } else if (stepEval >= 0.40 && stepEval < 0.55){
        lamRate = 0.44; // 44% acceptance rate
    } else if (stepEval >= 0.55 && stepEval < finalEval){
        lamRate = 0.27; // 27% acceptance rate
    } else if (stepEval >= finalEval){ // exponential decay from 0.27
        lamRate = 0.27*pow(270, -(stepEval - finalEval)*2.857142857);
    }


    if (acceptRate > lamRate){
        temp = 0.999*temp;
        inv_temp = 1.0/temp;
//        changed=true;
    } else {
        temp = temp*1.001001001001;
        inv_temp = 1.0/temp;
//        changed=true;
    }

//    if (changed){
//        inv_temp = 1.0/temp;
//    }
}



void Anneal::printParameters(std::vector<float> * accept, std::vector<double> * temp, std::vector<float> * divergence, std::vector<unsigned int> * wl){

    int total = accept->size();

    FILE * pFile;
    const char *outputFileName;
    std::string nameOf = "run_parameters.txt";
    outputFileName = nameOf.c_str() ;
    pFile = std::fopen(outputFileName, "w");

    std::string index;

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );

    for (int i=0; i < total; i++){
        //residue_index = std::to_string(selectedBeads[i]);
        //index = std::to_string(i + 1);
        fprintf(pFile, "%i %.5E %0.9f %0.9f %i\n", (i+1), (*accept)[i], (*temp)[i], (*divergence)[i], (*wl)[i] );
    }

    fclose(pFile);

}

/**
 * Calcualte P(r) distribution using selected set of bead_indices (bead indices must be sorted!)
 * binCount is a small vector
 * Should not use if max index of pModelBin is greater than unsigned_int max
 * @param bead_indices
 * @param binCount
 * @param workingLimit
 * @param totalBeadsInSphere
 * @param pModel
 * @param pData
 * @return
 */
void Anneal::calculateModelPrDistribution(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, unsigned int workingLimit, unsigned int totalBeadsInSphere, Model *pModel, Data *pData) {

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);

    unsigned int row;
    unsigned int row2;
    const unsigned short int * pModelBin = pModel->getPointerToBins();

    //const int * pBeadIndices = &(*bead_indices)[0];
    const unsigned int * ptr = (*bead_indices).data();

    // calculate P(r) for beads
    for(unsigned int m=0; m < workingLimit; m++){ // set row

        //row = (unsigned int)*(pBeadIndices + m);
        row = ptr[m];
        /*
         * if row is zero, row2 can not be less than zero.
         * move -1 into the loop, insures number is never assigned non-negative
         * Underflow error possible
         */
        row2 = (row*totalBeadsInSphere) - (row*(row+1))/2 - row;// - 1;

        for(unsigned int n=(m+1); n < workingLimit; n++){ // iterate over columns
            (*binCount)[ *( pModelBin + ((row2 + ptr[n]) - 1) ) ]++;
        }
    }

}



/**
 * Calcualte P(r) distribution using selected set of bead_indices (bead indices must be sorted!)
 * binCount is a small vector
 * Should not use if max index of pModelBin is greater than unsigned_int max
 * @param bead_indices
 * @param binCount
 * @param workingLimit
 * @param totalBeadsInSphere
 * @param pModel
 * @param pData
 * @return
 */
void Anneal::calculateModelPrDistributionDirect(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, const unsigned int workingLimit, Model *pModel, Data *pData){

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);
    // calculate P(r) for beads
    // calculate distances and determine bin
    const unsigned int * const pBead = (*bead_indices).data();

    for(unsigned int m=0; m < workingLimit; m++){ // set row

        const vector3 * firstVec = &(pModel->getBead(pBead[m])->getVec());

        for(unsigned int n=(m+1); n < workingLimit; n++){ // iterate over columns
            (*binCount)[ pData->convertToBin((*firstVec - pModel->getBead(pBead[n])->getVec()).length()) ]++;
        }
    }
}


/**
 * Given the selectedIndex, compile list of indices that are direct neighbors in use
 * Pick one at random
 *
 */
unsigned int Anneal::getConnectedNeighborFromSet(std::set<unsigned int> *beads_in_use,
                                               Model *pModel,
                                               unsigned int & selectedIndex){

    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return neighborLimit;
    } else if (count == 1){
        return neighborsInuse[0];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);

    return neighborsInuse[randomIndex(gen)];
}





/**
 * pick random point that is in use
 * get set of neighbors that are available
 * randomly select one
 */
unsigned int Anneal::getUseableNeighborFromSet(std::set<unsigned int> *beads_in_use,
                                              Model *pModel,
                                              unsigned int & selectedIndex){


    const unsigned int * ptr = pModel->getDirectPointerToNeighborhood();

    //const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        //unsigned int neighbor = *(it+i);
        unsigned int neighbor = ptr[totalNeighbors*selectedIndex + i];
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) == beads_in_use->end()){ // if end of set, means not in use
            possibleNeighbors[count] = neighbor;
            //std::cout << count << " <= " << i << " - " << neighbor << std::endl;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return neighborLimit;
    } else if (count == 1){
        return possibleNeighbors[0];
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);
    //if (count == totalNeighbors){
    //    std::cout << count <<  " limit neighbors for : " << randomIndex(gen) << " " << possibleNeighbors[randomIndex(gen)] << std::endl;
    //}

    return possibleNeighbors[randomIndex(gen)];
}

unsigned int Anneal::numberOfContactsFromSet(std::set<unsigned int> *beads_in_use,
                                   Model *pModel,
                                   unsigned int const selectedIndex){

    auto it = pModel->getPointerToNeighborhood(selectedIndex);
    unsigned int neighborContacts = 0;

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighhborLimit = pModel->getNeighborLimit();

    for (unsigned int i=0; i< totalNeighbors; i++){

        unsigned int neighbor = *(it+i);

        if ((neighbor < neighhborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){
            neighborContacts += 1;
        } else if (neighbor == neighhborLimit) {
            break;
        }
    }

    return neighborContacts;
}


/**
 * pick random point that is in use
 * get set of neighbors that are available
 * randomly select one
 */
unsigned int Anneal::getUseableNeighborFromRestrictedSet(std::set<unsigned int> *available_beads,
                                               Model *pModel,
                                               unsigned int & selectedIndex){

    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    auto endOfSet = available_beads->end();
    unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    unsigned int count=0;
    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && available_beads->find(neighbor) != endOfSet){ // if end of set, means not in use
            possibleNeighbors[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return neighborLimit;
    } else if (count == 1){
        return possibleNeighbors[0];
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);

    return possibleNeighbors[randomIndex(gen)];
}



/**
 * pre-calculation, populate vector of bins for all pairwise distances
 *
 * @param pModel
 * @param pData
 */
void Anneal::fillPrBinsAndAssignTotalBin(Model * pModel, Data * pData){

    maxbin = pModel->populateBins(pData); // populate the entire

    maxbin += 1;
    totalBins = pData->getShannonBins(); // set by experimental input Data file
    // maxbin and totalBins will be differenct
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

}


///**
// * The new position to Add MUST BE WITHIN THE WORKING LIMIT
// * beadsInUse must be sorted up to workingLimit
// *
// * @param addMeSubUnitIndex
// * @param beadsInUse
// * @param workingLimit
// * @param prBins
// * @param pModel
// * @param pData
// * @return
// */
//inline unsigned int Anneal::addToPrSym(const unsigned int addMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){
//    // each bead has a sym related partner in pModel
//    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
//    unsigned int findIt=0;
//    const unsigned int * ptr = &beadsInUse.front();
//    for(unsigned int i=0; i<workingLimit; i++){
//        if (ptr[i] == addMeSubUnitIndex){
//            findIt = i;
//            break;
//        }
//    }
//
//    const unsigned int stopAt = findIt;
//    //
//    // add interdomain distances
//    //
//    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
//    unsigned int violations=0;
//    float distance_to;
//
//    const unsigned int subUnits = pModel->getNumberOfSubUnits();
//    const unsigned int totalCoordinates = subUnits*workingLimit;
//    std::vector<vector3> coordinates(totalCoordinates);
//    vector3 * tempVec1;
//
//    // create first subunit from selected indices and populate coordinates
//    for (unsigned int i=0; i < workingLimit; i++){
//        Bead * tempBead = pModel->getBead(beadsInUse[i]);
//        tempVec1 = &coordinates[i];
//        (*tempVec1).x = tempBead->getX();
//        (*tempVec1).y = tempBead->getY();
//        (*tempVec1).z = tempBead->getZ();
//    }
//
//    unsigned int count = workingLimit;
//    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
//        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
//    }
//
//    // adjust Pr
//    for (unsigned int s=0; s < subUnits; s++){  //
//        unsigned int basis = stopAt+s*workingLimit;
//        const vector3 * prefVec = &coordinates[basis];
//
//        for (unsigned int next_i = 0; next_i < basis; next_i++) {
//            // calculate distance and convert to bin
//            distance_to = ((*prefVec) - coordinates[next_i]).length();
//            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
//
//            if (distance_to < violation_limit){
//                ++violations;
//            }
//        }
//
//        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
//            // calculate distance and convert to bin
//            distance_to = ((*prefVec) - coordinates[next_i]).length();
//            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
//
//            if (distance_to < violation_limit){
//                ++violations;
//            }
//        }
//    }
//
//
//    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
//        unsigned int basis = stopAt+s*workingLimit;
//        const vector3 * prefVec = &coordinates[basis];
//
//        for (unsigned int ss=s+1; ss < subUnits; ss++){
//            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
//            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
//            if (distance_to < violation_limit){
//                --violations;
//            }
//        }
//    }
//
//    return violations;
//}


/**
 * modelPR and targetPR are the same size
 * targetPR is derived from PDB
 */
float Anneal::calculateKLDivergenceAgainstPDBPR(std::vector<unsigned int> &modelPR, std::vector<double> &targetPR){

    double totalCounts = 0.0;
    double kl=0.0;
    double prob, * value;

    unsigned int totalm = modelPR.size();
    std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    //last nonzero bin
    unsigned int last=0;
    for (unsigned int i=0; i < totalm; i++){
        if (targetPR[last] <= 0){
            break;
        }
        last++;
    }

    for (unsigned int i=0; i < totalm; i++){
        totalCounts += modelPR_float[i];
    }

    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    for (unsigned int i=0; i < last; i++){
        prob = targetPR[i];  // get experimental bin value
        value = &modelPR_float[i]; // get model value
        if (prob > 0 && *value > 0){
            kl += prob * log(prob/(*value) * totalCounts);
        } else if (prob > 0 && *value <= 0){ // severely penalize any model bin that is zero
            kl += 10.0;
        }
    }

    return (float)(kl/(float)last);  // returns value per bin
}



bool Anneal::checkForRepeats(std::vector<unsigned int> beads) {
    bool state = false;
    unsigned int beadSize = beads.size();

    std::set<int> testSet(beads.begin(), beads.end());
    std::cout << "______________________________________________________________________________" << std::endl;
    std::cout << "*******************                 TEST                   *******************" << std::endl;
    std::cout << "*******************              -----------               *******************" << std::endl;
    std::cout << " TEST SET " << testSet.size() << " vector set " << beadSize << std::endl;
    if (testSet.size() != beadSize){
        for(int i=0; i<10; i++){
            std::cout << "                 !!!!!!!!!!!!!!!!DUPLICATE ENTRIES FOUND! " << testSet.size() << " " << beadSize <<  std::endl;
        }
        state = true;
    }
    return state;
}

/**
 * use to update distribution model
 * Input distribution must contain zereos for values not specified
 *
 * @param contactsDistributionSeed
 */
void Anneal::setContactsDistribution(std::vector<unsigned int> & contactsDistributionSeed){

    distributionlimit=contactsDistributionSeed.size();
    contactsDistribution.resize(distributionlimit);

    std::copy(contactsDistributionSeed.begin(), contactsDistributionSeed.end(), contactsDistribution.begin());

    double totaltemp =0;
    for(auto & item : contactsDistribution){
        totaltemp += item;
    }

    for(unsigned int i=0; i<distributionlimit; i++){ // normalize
        contactsDistribution[i] *= 1.0/totaltemp;
    }
}

bool Anneal::checkSetAndVector(unsigned int workingLimit, std::vector<unsigned int> * indices, std::set<unsigned int> * beads_in_use){

    // check all beadss in set
    for(auto it = beads_in_use->begin(); it != beads_in_use->end(); ++it){

        auto vit = std::find(indices->begin(), indices->end(), *it);

        if (vit == indices->end()){
            std::cout << "************** ==> SET INDEX exceeds VECTOR " << *it << " LIMIT:   "<< indices->size() << std::endl;
            std::cout << "************** ==> SET SIZE " << beads_in_use->size() << " VECTOR SIZE " << indices->size() << std::endl;
            exit(0);
            return true;
        }else{
            unsigned int dis = std::distance(indices->begin(), vit);

            if (dis > workingLimit){
                std::cout << "************** ==> SET INDEX exceeds workinglimit " << *it << std::endl;
                exit(0);
                return true;
            }
        }
    }


    for(unsigned int i=0;i<workingLimit; i++){
        unsigned int index = (*indices)[i];
        auto vit = beads_in_use->find(index);
        if (vit == beads_in_use->end()){
            std::cout << "************** ==> VECTOR INDEX NOT IN SET " << index << std::endl;
            exit(0);
            return true;
        }
    }

    return false;
}


/**
 * bead_indices should be sorted to make insertion into hull as efficient as possible
 *
 * @param workingLimit
 * @param bead_indices
 * @param hull
 * @param pModel
 * @param totalBeadsInSphere
 */
void Anneal::recalculateDeadLimit(unsigned int workingLimit, std::vector<unsigned int> &bead_indices, std::set<unsigned int> &hull, Model * pModel, unsigned int totalBeadsInSphere ){

    hull.clear();

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[] = "qhull FA";

    coordT hullPoints2[3*workingLimit];
    // calculate CVX Hull from selected indices
    const auto num = (unsigned int)bead_indices.size();
    const unsigned int * ptr = (num > 0) ? bead_indices.data() : nullptr;
    std::vector<unsigned int> active_indices(workingLimit);

    for (unsigned int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints2[i*3], pModel->getBead(ptr[i]));
        active_indices[i] =ptr[i];
    }

    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, nullptr, nullptr);

    for(unsigned int i = workingLimit; i < num; i++) {
        unsigned int value = ptr[i];
        beadToPoint(testPoint, pModel->getBead(value));
        // exclude HULL points, for each bead, determine if outside HULL
        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);

        if (!isoutside){
            hull.insert(value);
        }
    }

    std::set<unsigned int> beads_in_use(bead_indices.begin(), bead_indices.begin() + workingLimit);

    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    unsigned int neighbor;
    vertexT * vertices = qh vertex_list;
    auto totalV = (unsigned int)qh num_vertices;

    // only move CVX hull points
    for (unsigned int v = 0; v < totalV; v++) { //

        auto pptr = pModel->getDirectPointerToNeighborhood();
        unsigned int location = totalNeighbors*active_indices[qh_pointid( vertices->point)];

        for (unsigned int j=0; j < totalNeighbors; j++){
            neighbor = pptr[location+j];
            if (neighbor < neighborLimit && beads_in_use.find(neighbor) == beads_in_use.end()){
                hull.insert(neighbor);
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
        vertices = vertices->next;
    }


    qh_freeqhull(true);
}


void Anneal::getHullPoints(std::set<unsigned int> &hullpts, std::set<unsigned int> &beads_in_use, Model * pModel){
    unsigned int workingLimit = beads_in_use.size();
    char flags[] = "qhull FA";
    hullpts.clear();

    coordT hullPoints2[3*workingLimit];
    std::vector<unsigned int> active_indices(workingLimit);
    // calculate CVX Hull from selected indices
    unsigned int count=0;
    for(const auto & ind : beads_in_use){
        beadToPoint(&hullPoints2[count*3], pModel->getBead(ind));
        active_indices[count] =ind;
        count++;
    }
    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, nullptr, nullptr);

    vertexT * vertices = qh vertex_list;
    auto totalV = (unsigned int)qh num_vertices;

    // only move CVX hull points
    for (unsigned int v = 0; v < totalV; v++) { //
        hullpts.insert(active_indices[qh_pointid( vertices->point)]);
        vertices = vertices->next;
    }

    qh_freeqhull(true);
}



/**
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within hull sett
 */
void Anneal::populateLayeredDeadlimitUsingSet(std::set<unsigned int> & beads_in_use, std::set<unsigned int> & hull, Model * pModel) {

    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    unsigned int neighbor;
    hull.clear();

    for(const auto & it : beads_in_use){
        //hull.insert(it);
//        auto neighborhood = pModel->getPointerToNeighborhood(it);
        const auto ptr = pModel->getDirectPointerToNeighborhood();
        unsigned int location = totalNeighbors*it;

        for (unsigned int j=0; j < totalNeighbors; j++){
            neighbor = ptr[location+j];

            if (neighbor < neighborLimit && beads_in_use.find(neighbor) == beads_in_use.end()){
                hull.insert(neighbor);
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
    }
}

const std::vector<double> & Anneal::getContactsDistribution() const {
    return contactsDistribution;
}

void Anneal::updateContactsDistributionToFile(std::string filename){

    contactsDistribution.clear();
    contactsDistribution.resize(13);

    if (!boost::filesystem::exists(filename)){
        throw std::invalid_argument("** ERROR FILE => HISTOGRAM FILE NOT FOUND : " + filename);
    }

    std::ifstream data (filename, std::ifstream::in);

    if (data.is_open()) {

        boost::regex dataFormat("[0-9]+[[:blank:]]+[0-9]+.[0-9]+");
        boost::regex remarkFormat("REMARK", boost::regex::icase);

        std::string line;

        while (!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */

            if (isspace(line[0])) {
                line.erase(line.begin(),
                           std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }


            if (boost::regex_search(line, dataFormat)) {
                std::vector<std::string> tempLine;
                boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

                int contacts = std::stoi(tempLine[0]);
                contactsDistribution[contacts] = std::stod(tempLine[1]);


            } else if (boost::regex_search(line, remarkFormat)) {
//                std::vector<std::string> rgLine;
//                boost::split(rgLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
//                this->rg = std::strtof(rgLine[5].c_str(), nullptr);

            }
        }

        data.close();
    }

    logger("UPDATED CONTACTS DISTRIBUTION","");
    for(unsigned int i=0; i<contactsDistribution.size(); i++){
        logger(std::to_string(i),std::to_string(contactsDistribution[i]));
    }
}

/**
 * write subset of selected bead to file
 */
void Anneal::writeContactsDistributionToFile(float bin_width, float bead_radius, float qmax){

    std::string nameOf = "histogram.histo";
    //const char * outputFileName = nameOf.c_str();
    FILE * pFile;
    pFile = fopen(nameOf.c_str(), "w");

    char buffer[80];

    std::string tempHeader = "REMARK 265\n";
    tempHeader += "REMARK 265 HISTOGRAM DETAILS\n";
    std::snprintf(buffer, 80, "REMARK 265         QMAX : %.4f INV ANGSTROMS\n", qmax);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265  EDGE LENGTH : %.3f ANGSTROMS\n", 2*bead_radius);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265     BINWIDTH : %.3f ANGSTROMS\n", bin_width);
    tempHeader.append(buffer);
    tempHeader.append("REMARK 265 \n");
    tempHeader.append("REMARK 265 NORMALIZED (AREA INTEGRATES TO UNITY) \n");
    tempHeader.append("REMARK 265 \n");
    std::snprintf(buffer, 80, "REMARK 265     COL %i => NUMBER OF CONTACTS \n", 1);
    tempHeader.append(buffer);
    std::snprintf(buffer, 80, "REMARK 265     COL %i => DISTRIBUTION \n", 2);
    tempHeader.append(buffer);
    tempHeader.append("REMARK 265 \n");

    unsigned int index=0;
    for(auto cc : contactsDistribution){
        std::snprintf(buffer, 80, " %-3d %.5f\n", index, cc);
        tempHeader.append(buffer);
        index++;
    }

    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}

