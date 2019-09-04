//
// Created by xos81802 on 14/07/2018.
//

#include "Component.h"
#include "../Model.h"

Component::Component(std::string id, float volume, Model *pModel) : id(id), targetVolume(volume), pModel(pModel) {
    std::cout <<" CREATING NEW COMPONENT " << id << " VOL: " << volume << std::endl;
    this->invTargetVolume = 1.0f/targetVolume;
    contactsDistributionOfModel.resize(13);
}


Component::Component(std::string id, float volume, Model *pModel, bool contiguous) : Component(id, volume, pModel){
    this->invTargetVolume = 1.0f/targetVolume;
    contactsDistributionOfModel.resize(13);
}


void Component::addResid(unsigned int index, std::string chain) {
    std::cout << " ADDING RESID : " << index << " CHAIN => " << chain << std::endl;
    resids.push_back(index);
    chains.push_back(chain);
}

void Component::addCenteredAnchors(unsigned int index){
    centeredAnchors.insert(index);
}

/**
 * value should be in units of beads per volume
 */
void Component::setTargetNumberOfLatticePoints(float value){
    this->targetNumberOfBeads = (unsigned int)ceil(targetVolume/value);
    if (anchors.size() > 0){
        this->targetNumberOfBeads += anchors.size();
    }
    // populate potential
    this->populatePotential(0.18);
    std::cout << " SET COMPONENT " << id << " TARGET NUMBER => " << targetNumberOfBeads << " " << (targetVolume/value) << std::endl;
}

unsigned int Component::getTargetNumberOfLatticePoints(){
    std::cout << " GET COMPONENT " << id << " TARGET NUMBER => " << this->targetNumberOfBeads << " vol " << targetVolume << std::endl;
    return this->targetNumberOfBeads;
}


// Euler tour should include the anchor points
unsigned int Component::addLatticePoint(unsigned int index){

//    if (beads_in_use.find(index) != beads_in_use.end()){
//        std::cout << " ERROR ADDING BEAD ALREADY PRESENT " << std::endl;
//        exit(0);
//    }

    beads_in_use.insert(index);
    totalComponents = tour.addNode(index, pModel);
    if (isAnchor(index)){
        anchorCount++;
    }
    return totalComponents;
}

unsigned int Component::removeLatticePoint(unsigned int index){

    //unsigned int before = this->getTotalNumberOfBeads();
    beads_in_use.erase(index);
    totalComponents = tour.removeNode(index);

    if (isAnchor(index)){
        anchorCount--;
    }

    return totalComponents;
}

// ADD ANCHOR OCCURS DURING INITIALIZATION FROM THE FILE
// IF NOT ANCHORS, THEN NOTHIHNG ADDED TO TOUR
void Component::addAnchor(unsigned int index){

    if (anchors.find(index) == anchors.end()){
        std::cout << "        ADDING ANCHOR " << index << " " << anchors.size() << std::endl;
        anchors.insert(index);
        anchorCount++;

//        beads_in_use.insert(index);
//        totalComponents = tour.addNode(index, pModel);

        empty=false;
    }
}

// if anchor point, returns true
bool Component::isAnchor(unsigned int index){
    return !(anchors.find(index) == anchors.end());
}

/**
 * Test if Index is an anchor
 * returns true if centered anchor point
 * @param index
 * @return
 */
bool Component::isCenteredAnchor(unsigned int index){
    return (centeredAnchors.find(index) != centeredAnchors.end());
}


void Component::printAnchors(){
    std::cout << anchors.size() << std::endl;
}


void Component::printCenteredAnchors(){
    for(auto anc : centeredAnchors){
        std::cout << "ID: " << this->id << " Centered Anchor " << anc << std::endl;
    }
    std::cout << " TOTAL ANCHORS " << anchors.size() << std::endl;
}


void Component::populatePotential(float percentage){

    percentageStep = 1.0f/((float)targetNumberOfBeads*percentage);

    potentialFunction.resize(6);
    potentialFunction[0] = 0;
    potentialFunction[1] = 0.001;
    potentialFunction[2] = 0.005;
    potentialFunction[3] = 0.01;
    potentialFunction[4] = 0.1;
    potentialFunction[5] = 1.1;
}

/**
 * Calculate the difference between expected and actual number of beads as a percentage and then take as a ratio to percentage window
 * so if you are within 10% and have a window of 10%, then index is 0
 * @return
 */
float Component::potential() {

    auto totalInUse = (unsigned int)beads_in_use.size();
    float increment = std::abs((float)totalInUse - targetNumberOfBeads) * percentageStep;
    unsigned int locale = (increment < 1) ? 0 : (unsigned int)std::ceil(increment - 1);
    locale = locale > 5 ? 5 : locale;

    if (totalInUse < 6){
        return 100.5;
    }

    return potentialFunction[locale];
}

float Component::calculateCVXVolume() {
    auto totalInUse = (unsigned int)beads_in_use.size();
    char flags[] = "qhull FA";
    unsigned int numpoints = 3*totalInUse;
    coordT points[numpoints];

    std::vector<unsigned int> active_indices(totalInUse);

    unsigned int next, countIt=0;
    for (auto it = beads_in_use.begin(); it != beads_in_use.end(); ++it) {
        next = 3*countIt;
        Bead * pBead = pModel->getBead(*it);
        points[next] = pBead->getX();
        points[next+1] = pBead->getY();
        points[next+2] = pBead->getZ();
        active_indices[countIt] = *it;
        countIt++;
    }

    // needs to be optimized
    qh_new_qhull(3, countIt, points, 0, flags, nullptr, nullptr);
//    vertexT * vertices = qh vertex_list;
//    int totalV = qh num_vertices;
//    cvxPoints.resize(totalV);
//
//    for (int v = 0; v < totalV; v++) { //
//        cvxPoints[v] = active_indices[qh_pointid(vertices->point)];
//        vertices = vertices->next;
//    }
//
    float tempVolume = qh totvol;
    qh_freeqhull(true);
    return tempVolume;
}

void Component::printSet(){
    pModel->printBeadsFromSet(beads_in_use);
}

void Component::printBest(){
    pModel->printBeadsFromSet(best);
}


void Component::setBest(){

    //best = std::set<unsigned int>(beads_in_use.begin(), beads_in_use.end());
    best.clear();
    for(auto bead : beads_in_use){
        best.insert(bead);
    }
}


void Component::copyBestToInUse(){

    beads_in_use.clear();
    for(auto bead : best){
        beads_in_use.insert(bead);
    }

    std::vector<unsigned int> temp;

    anchorCount=0;
    for (auto it = best.begin(); it != best.end(); ++it) {
        //beads_in_use.insert(*it);
        temp.push_back(*it);
        if (isAnchor(*it)){
            anchorCount++;
        }
    }

    std::sort(temp.begin(), temp.end());
    totalComponents = tour.newTour(temp.begin(), temp.size(), pModel);
}


/**
 * write subset of selected bead to file
 */
void Component::writeToFile(std::string nameOf){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    std::string residue_index;

    std::vector<unsigned int> indices(beads_in_use.begin(), beads_in_use.end());
//    for(auto it = beads_in_use.begin(); it != beads_in_use.end(); ++it){
//        indices.push_back(*it);
//    }

    std::sort(indices.begin(), indices.end());
    fprintf(pFile, "REMARK 265         NUMBER LATTICE POINTS : %s\n", nameOf.c_str());
    fprintf(pFile, "REMARK 265                        TARGET : %i\n", targetNumberOfBeads );
    fprintf(pFile, "REMARK 265                        ACTUAL : %i\n", (int)beads_in_use.size() );
    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    unsigned int residueCnt = 1;
    for (unsigned int i=0; i<indices.size(); i++){
        currentBead = pModel->getBead(indices[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(residueCnt);
        fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}


/**
 * write subset of selected bead to file
 */
void Component::writeAnchorsToFile(std::string nameOf){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    std::string residue_index;

    std::vector<unsigned int> indices;
    for(auto it = anchors.begin(); it != anchors.end(); ++it){
        indices.push_back(*it);
    }

    for(auto it = centeredAnchors.begin(); it != centeredAnchors.end(); ++it){
        indices.push_back(*it);
    }

    std::sort(indices.begin(), indices.end());

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    unsigned int residueCnt = 1;
    for (unsigned int i=0; i<indices.size(); i++){
        currentBead = pModel->getBead(indices[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(residueCnt);
        fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}


void Component::printConstraints(){
    std::cout<< " TARGET NUMBER BEADS => " << targetNumberOfBeads << "( " << beads_in_use.size() << " )" << std::endl;
}


void Component::setCurrentCVXVolume(){
    this->currentCVXVolume = this->calculateCVXVolume();
}

void Component::printContactsPerCenteredAnchor(){

    auto endOfSet = beads_in_use.end();
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighhborLimit = pModel->getNeighborLimit();

    for(auto canchor : centeredAnchors){

        auto it = pModel->getPointerToNeighborhood(canchor);
        unsigned int neighborContacts = 0;

        for (unsigned int i=0; i< totalNeighbors; i++){

            unsigned int neighbor = *(it+i);

            if ((neighbor < neighhborLimit ) && beads_in_use.find(neighbor) != endOfSet){
                neighborContacts += 1;
            } else if (neighbor == neighhborLimit) {
                break;
            }
        }

        std::cout << "COMP ID: " << this->getID() << " CENTERED ANCHOR " << canchor << " TOTAL CONTACTS " << neighborContacts << std::endl;
    }

}

/**
 * write subset of selected bead to file
 */
void Component::writeCenteredAnchorsToFile(){
    FILE * pFile;

    const char *outputFileName;
    std::string nameOf = this->getID()+"_centeredAnchors.pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    std::string residue_index;

    std::vector<unsigned int> indices;


    for(auto it : centeredAnchors){
        indices.push_back(it);
    }

    std::sort(indices.begin(), indices.end());

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    unsigned int residueCnt = 1;
    for (unsigned int i=0; i<indices.size(); i++){
        currentBead = pModel->getBead(indices[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(residueCnt);
        fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}
