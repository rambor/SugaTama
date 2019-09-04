//
// Created by xos81802 on 19/09/2018.
//


#include "SubUnit.h"
#include "../Model.h"

SubUnit::SubUnit(){}; //


SubUnit::SubUnit(unsigned int pivotIndex, unsigned int subUnitIndex, unsigned int translateTo) : pivotIndex(pivotIndex) , subUnitIndex(subUnitIndex), translateTo(translateTo) {

}

/**
 *
 * @param subUnitIndex
 * @param pivotIndex
 * @param indices (must be sorted)
 * @param pModel
 */
SubUnit::SubUnit(unsigned int subUnitIndex, unsigned int pivotIndex, unsigned int translateTo, vector3 * translationVec, unsigned int workingLimit, std::vector<unsigned int> &indices, Model *pModel) : SubUnit(pivotIndex, subUnitIndex, translateTo){

    // subtract cofmVector and add TranslateTo

    this->workingLimit = workingLimit;
    capacity = 2*workingLimit;

    indicesInUse.resize(capacity);
    backedUpindicesInUse.resize(capacity);

    coordinates.resize(capacity);
    backedUpCoordinates.resize(capacity);

    transformedCoordinates.resize(capacity);
    backedUpTransformedCoordinates.resize(capacity);

    std::copy(indices.begin(), indices.begin()+workingLimit, indicesInUse.begin());

    this->pPivotPoint = pModel->getBead(this->pivotIndex);
//    this->pTranslateToPoint = pModel->getBead(this->translateTo); // original position in universe

    this->pivotVector = vector3(this->pPivotPoint->getVec().x, this->pPivotPoint->getVec().y, this->pPivotPoint->getVec().z);
    this->translateToVector = vector3(translationVec->x, translationVec->y, translationVec->z);

    backedUpCenteringVector = vector3(0,0,0);
    backedUpTranslationVector = vector3(0,0,0);

    this->backedUpPivotIndex = this->pivotIndex;

    // convert indices to coordinates and translate to cofm
    for(unsigned int i=0; i < this->workingLimit; i++){
        Bead * tempBead = pModel->getBead(this->indicesInUse[i]);
        coordinates[i] = vector3(tempBead->getVec().x, tempBead->getVec().y, tempBead->getVec().z);
    }

    // random rotation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randomIndex(0,360-1); // guaranteed unbiased

    alpha_x = randomIndex(gen)*convert;
    gamma_z = randomIndex(gen)*convert;
    std::uniform_int_distribution<int> randomBeta(0,180-1); // guaranteed unbiased
    beta_y = randomBeta(gen)*convert;

    this->rotate(alpha_x, beta_y, gamma_z);
}

/**
 *
 * @param subUnitIndex
 * @param pivotIndex
 * @param indices (must be sorted)
 * @param pModel
 */
SubUnit::SubUnit(unsigned int subUnitIndex, unsigned int pivotIndex, unsigned int translateTo, unsigned int workingLimit, std::vector<unsigned int> &indices, Model *pModel) : SubUnit(pivotIndex, subUnitIndex, translateTo){

    // subtract cofmVector and add TranslateTo

    this->workingLimit = workingLimit;
    capacity = 2*workingLimit;

    indicesInUse.resize(capacity);
    backedUpindicesInUse.resize(capacity);

    coordinates.resize(capacity);
    backedUpCoordinates.resize(capacity);

    transformedCoordinates.resize(capacity);
    backedUpTransformedCoordinates.resize(capacity);

    std::copy(indices.begin(), indices.begin()+workingLimit, indicesInUse.begin());

    pPivotPoint = pModel->getBead(this->pivotIndex);
//    pTranslateToPoint = pModel->getBead(this->translateTo);
    vector3 pTempTranslateToPoint = pModel->getBead(this->translateTo)->getVec();

    pivotVector = vector3(this->pPivotPoint->getVec().x, this->pPivotPoint->getVec().y, this->pPivotPoint->getVec().z);
//    translateToVector = vector3(this->pTranslateToPoint->getVec().x, this->pTranslateToPoint->getVec().y, this->pTranslateToPoint->getVec().z);
    translateToVector = vector3(pTempTranslateToPoint.x, pTempTranslateToPoint.y, pTempTranslateToPoint.z);

    backedUpCenteringVector = vector3(0,0,0);
    backedUpTranslationVector = vector3(0,0,0);

    this->backedUpPivotIndex = this->pivotIndex;

    // convert indices to coordinates
    for(unsigned int i=0; i < this->workingLimit; i++){
        Bead * tempBead = pModel->getBead(this->indicesInUse[i]);
        coordinates[i] = vector3(tempBead->getVec().x, tempBead->getVec().y, tempBead->getVec().z );
    }

    // random rotation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randomIndex(0,360); // guaranteed unbiased

    alpha_x = randomIndex(gen)*convert;
    gamma_z = randomIndex(gen)*convert;
    std::uniform_int_distribution<int> randomBeta(0,180); // guaranteed unbiased
    beta_y = randomBeta(gen)*convert;

    this->rotate(alpha_x, beta_y, gamma_z);
}

void SubUnit::makeBackUp() {

    this->backedUpPivotIndex = pivotIndex;
    this->backedUpTranslateToIndex = translateTo;

    pBackUpPivotPoint = &(*pPivotPoint);
//    pBackUpTranslateToPoint = &(*pTranslateToPoint);

    this->alpha_x_old = this->alpha_x;
    this->beta_y_old = this->beta_y;
    this->gamma_z_old = this->gamma_z;

    backedUpCenteringVector.x = pivotVector.x;
    backedUpCenteringVector.y = pivotVector.y;
    backedUpCenteringVector.z = pivotVector.z;

    backedUpTranslationVector.x = translateToVector.x;
    backedUpTranslationVector.y = translateToVector.y;
    backedUpTranslationVector.z = translateToVector.z;

    this->x1_old = x1;
    this->x2_old = x2;
    this->x3_old = x3;
    this->y1_old = y1;
    this->y2_old = y2;
    this->y3_old = y3;
    this->z1_old = z1;
    this->z2_old = z2;
    this->z3_old = z3;

    //std::copy(coordinates.begin(), coordinates.end(), backedUpCoordinates.begin());
    std::copy(transformedCoordinates.begin(), transformedCoordinates.end(), backedUpTransformedCoordinates.begin());
    //std::copy(indicesInUse.begin(), indicesInUse.end(), backedUpindicesInUse.begin());
}


/**
 * Undoing Add means decrementing workingLimit
 * @param indexToAdd
 * @param pModel
 */
void SubUnit::addIndex(unsigned int indexToAdd, Model *pModel ){

    if (workingLimit >= capacity){
        capacity += workingLimit;
        indicesInUse.resize(capacity);
        coordinates.resize(capacity);
        transformedCoordinates.resize(capacity);
        backedUpTransformedCoordinates.resize(capacity);
        backedUpindicesInUse.resize(capacity);
        backedUpCoordinates.resize(capacity);
    }

    this->indicesInUse[workingLimit] = indexToAdd; // will have repeats
    // constant point to a const vector
    const vector3 * const tempBeadVec = &(pModel->getBead(indexToAdd)->getVec());

    this->coordinates[workingLimit] = vector3(tempBeadVec->x, tempBeadVec->y, tempBeadVec->z);
    // needs to be rotated and translated
    vector3 * pVec = &coordinates[workingLimit];

    // rotate and translate
    this->transformPoint(pVec, &transformedCoordinates[workingLimit]);
    //this->rotate(alpha_x, beta_y, gamma_z);
    ++workingLimit;
}


/**
 * newly added points are always last, just have to track workinglimit
 */
void SubUnit::undoAdd(){
    workingLimit--;
}


void SubUnit::undoRemove(){
    workingLimit++;
}


void SubUnit::removeIndex(unsigned int indexToRemove){

    auto locale = std::find(indicesInUse.begin(), indicesInUse.begin() + workingLimit, indexToRemove);
    int distance = std::distance(indicesInUse.begin(), locale);

    workingLimit--;
    // swap the end and decrement workingLimit?
    std::iter_swap(locale, indicesInUse.begin() + workingLimit);
    std::iter_swap(coordinates.begin()+distance, coordinates.begin() + workingLimit);
    std::iter_swap(transformedCoordinates.begin()+distance, transformedCoordinates.begin() + workingLimit);
}

/**
 * Angles are in radians
 * Rotation is RxRyRz
 * @param alpha
 * @param beta
 * @param gamma
 */
void SubUnit::rotate(float alpha, float beta, float gamma) {

    float cosAlpha = cos(alpha);
    float sinAlpha = sin(alpha);
    float cosBeta = cos(beta);
    float sinBeta = sin(beta);
    float cosGamma = cos(gamma);
    float sinGamma = sin(gamma);

    x1 = cosGamma*cosBeta;
    x2 = sinGamma*cosBeta;
    x3 = -sinBeta;
    y1 = -sinGamma*cosAlpha + cosGamma*sinBeta*sinAlpha;
    y2 = cosGamma*cosAlpha + sinGamma*sinBeta*sinAlpha;
    y3 = cosBeta*sinAlpha;
    z1 = sinGamma*sinAlpha + cosGamma*sinBeta*cosAlpha;
    z2 = -cosGamma*sinAlpha + sinGamma*sinBeta*cosAlpha;
    z3 = cosBeta*cosAlpha;

    vector3 * tempVec1;

    for(unsigned int i=0; i<workingLimit; i++){
        tempVec1 = &transformedCoordinates[i];
        const vector3 baseVec = coordinates[i] - pivotVector;

        (*tempVec1).x = x1*(baseVec).x + y1*(baseVec).y + z1*(baseVec).z + translateToVector.x;
        (*tempVec1).y = x2*(baseVec).x + y2*(baseVec).y + z2*(baseVec).z + translateToVector.y;
        (*tempVec1).z = x3*(baseVec).x + y3*(baseVec).y + z3*(baseVec).z + translateToVector.z;
    }
}

/**
 * Transform lattice point using current alpha, beta and gamma to center-of-mass
 * @param baseVec
 * @param newVewc
 */
void SubUnit::transformPoint(const vector3 *base, vector3 *newVewc) {

    const vector3 baseVec = *base - pivotVector; // when to update centering vector?
    // translation Vector is with respect to the designated neighboring subUnit

    (*newVewc).x = x1*baseVec.x + y1*baseVec.y + z1*baseVec.z + translateToVector.x;
    (*newVewc).y = x2*baseVec.x + y2*baseVec.y + z2*baseVec.z + translateToVector.y;
    (*newVewc).z = x3*baseVec.x + y3*baseVec.y + z3*baseVec.z + translateToVector.z;
}

vector3 SubUnit::translateAndRotatePoint(const vector3 * const originalVec) {
    const vector3 baseVec = *originalVec - pivotVector; // when to update centering vector?
    // translation Vector is with respect to the designated neighboring subUnit
    return vector3(x1*baseVec.x + y1*baseVec.y + z1*baseVec.z + translateToVector.x,
                   x2*baseVec.x + y2*baseVec.y + z2*baseVec.z + translateToVector.y,
                   x3*baseVec.x + y3*baseVec.y + z3*baseVec.z + translateToVector.z);
}

unsigned int SubUnit::getIndexOfBeadtoRemove(unsigned int indexToRemove) {
    auto locale = std::find(indicesInUse.begin(), indicesInUse.begin() + workingLimit, indexToRemove);
    return std::distance(indicesInUse.begin(), locale);
}


/**
 * translate to is a Index within the subUnit
 *
 * @param newPivot
 * @param translatedTo
 * @param pModel
 */
void SubUnit::translateToTransformedPosition(unsigned int translatedTo, Model * pModel){
    this->makeBackUp();

//    const vector3 * newPivotVector = &pModel->getBead(newPivot)->getVec();
    /*
     * UNDO TRANFORMATION
     * subtract translateToVector -> should be centered on the pivot
     * add back the pvot and now have an object within Universe but rotated
     * Subtract the newPivot to center
     * Add the new translation vector
     */
    // reassign CenterOfRotation/MASS
//    pivotVector = vector3(*newPivotVector);
//    pivotIndex = newPivot;

//    this->pPivotPoint = pModel->getBead(pivotIndex);

    this->translateTo = translatedTo;
    this->translateToVector = vector3(pModel->getBead(translatedTo)->getVec()); // translate to transformed position within molecule

    vector3 * tempVec1;

    /*
     * re-transform coordinates to new translation point
     */
    for(unsigned int i=0; i<workingLimit; i++){
        tempVec1 = &transformedCoordinates[i];
        const vector3 baseVec = coordinates[i] - pivotVector;

        (*tempVec1).x = x1*(baseVec).x + y1*(baseVec).y + z1*(baseVec).z + translateToVector.x;
        (*tempVec1).y = x2*(baseVec).x + y2*(baseVec).y + z2*(baseVec).z + translateToVector.y;
        (*tempVec1).z = x3*(baseVec).x + y3*(baseVec).y + z3*(baseVec).z + translateToVector.z;
    }
}

void SubUnit::undoTranslateTo(Model *pModel) {

    this->pivotIndex = backedUpPivotIndex;
    this->translateTo = backedUpTranslateToIndex;

    pivotVector.x = backedUpCenteringVector.x;
    pivotVector.y = backedUpCenteringVector.y;
    pivotVector.z = backedUpCenteringVector.z;

    translateToVector.x = backedUpTranslationVector.x;
    translateToVector.y = backedUpTranslationVector.y;
    translateToVector.z = backedUpTranslationVector.z;

    this->pPivotPoint = pModel->getBead(backedUpPivotIndex);

    std::copy(backedUpTransformedCoordinates.begin(), backedUpTransformedCoordinates.end(), transformedCoordinates.begin());
}


/*
 * new Pivot Point must exist within selected indices, the pivot is a point within the selected Indices
 * tests should be withihn the funciton calling rotateAt()
 *
 */
void SubUnit::rotateAt(unsigned int newPivotPoint, float angleInDegreesX, float angleInDegreesY, float angleInDegreesZ, Model * pModel){

    this->makeBackUp();

    pivotIndex = newPivotPoint;

    pPivotPoint = pModel->getBead(newPivotPoint); // reference to a bead in universe

    pivotVector = vector3(pPivotPoint->getVec().x, pPivotPoint->getVec().y, pPivotPoint->getVec().z);

    this->alpha_x += angleInDegreesX*convert;
    this->beta_y += angleInDegreesY*convert;
    this->gamma_z += angleInDegreesZ*convert;

    this->rotate(alpha_x, beta_y, gamma_z);
}

void SubUnit::undoRotateTo() {

    this->alpha_x = alpha_x_old;
    this->beta_y = beta_y_old;
    this->gamma_z = gamma_z_old;

    this->x1 = x1_old;
    this->x2 = x2_old;
    this->x3 = x3_old;
    this->y1 = y1_old;
    this->y2 = y2_old;
    this->y3 = y3_old;
    this->z1 = z1_old;
    this->z2 = z2_old;
    this->z3 = z3_old;

    pivotVector.x = backedUpCenteringVector.x;
    pivotVector.y = backedUpCenteringVector.y;
    pivotVector.z = backedUpCenteringVector.z;

    pPivotPoint = &(*pBackUpPivotPoint);
    pivotIndex = backedUpPivotIndex;

    std::copy(backedUpTransformedCoordinates.begin(), backedUpTransformedCoordinates.end(), transformedCoordinates.begin());
}

void SubUnit::swapIndices(unsigned int indexOfCurrent, unsigned int indexOfNew, Model * pModel) {
// find position to swap
    auto position1 = std::find(indicesInUse.begin(), indicesInUse.begin() + workingLimit, indexOfCurrent);
    swapLocale = std::distance(indicesInUse.begin(), position1);

    vector3 const * pNewVec = &(pModel->getBead(indexOfNew)->getVec()); // new value
    vector3 * pOldVec = &coordinates[swapLocale]; // pointing to old value

    //change value to vew value
    //*position1 = indexOfNew;
    indicesInUse[swapLocale] = indexOfNew;
    pOldVec->x = pNewVec->x;
    pOldVec->y = pNewVec->y;
    pOldVec->z = pNewVec->z;

    // update transformed coordinate position
    // rotate and translate
    this->transformPoint(pOldVec, &transformedCoordinates[swapLocale]);
    unsigned int finalLocale = workingLimit - 1; // move new vector to end of selected set

    /*
     * move newly transformed position to => workingLimit - 1
     * AddToPrSymX assumes newly added position is at end
     */
    std::iter_swap(transformedCoordinates.begin() + swapLocale, transformedCoordinates.begin() + finalLocale);
    std::iter_swap(coordinates.begin() + swapLocale, coordinates.begin() + finalLocale);
    std::iter_swap(indicesInUse.begin() + swapLocale, indicesInUse.begin() + finalLocale);
}