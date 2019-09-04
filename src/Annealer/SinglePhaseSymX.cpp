//
// Created by xos81802 on 23/07/2018.
//

#include <SubUnit.h>
#include <Logger.h>
#include "../Anneal.h"
#include "Data.h"
#include "../Model.h"
#include "PDBModel.h"
#include "../EulerTour/EulerTour.h"
#include "CEPotential.h"
#include "CEConfiguration.h"

using namespace std;

/**
 * create initial model of complete object
 */
bool Anneal::createInitialModelSymmetryX(Model *pModel, Data *pData) {

    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::clock_t startTime;
    unsigned int counter = 0;

    char addRemoveText[50];
    char titleText[50];

    contactCutOff = interconnectivityCutOff;
    violation_limit = 1.5d*pModel->getBeadRadius();

    //float inv_kb_temp = 1.0/highT;
    double lowTempStop = 0.01;
    double inv_kb_temp = 1.0/lowTempStop;
    float acceptRate = 0.5, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    unsigned int volumeCount=0;
    float volumeSum=0, workingLimitSum=0, sum_x_squared=0;

    // convert distances to ShannonBin membership
    double violations, temp_violations;
    maxbin = pModel->populateBins(pData);

    // adjust this part
    maxbin = pData->convertToBin( pModel->getDiameterOfUniverse()*3 );
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    std::cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    std::cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    std::cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // as bead indices are discarded, set upper limit of vector
    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> lowest_subUnit_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    std::vector<unsigned int> cvxIndices;
    std::vector<unsigned int> tempCVXIndices;

    /*
     * set expected volume constraints
     */
    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.21);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

//    const float targetSubUnitVolume = (lowerV+upperV)*0.5/(float)pModel->getNumberOfSubUnits();
    const float targetSubUnitVolume = (float)lowerV/(float)pModel->getNumberOfSubUnits();

    const float radius_larger = (float)(pModel->getBeadRadius()*std::sqrt(7.0)/2.0);
    auto lowerN = (unsigned int)std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits());
    auto upperN = (unsigned int)std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits());
    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN);   // number of beads in ASU

    auto subUnitWorkingLimit = (unsigned int)number_of_beads_to_use(gen);
    unsigned int lowestWorkingLimit = subUnitWorkingLimit;

    std::cout << "Bead Search Limited to: " << lowerN << " <= N <= " << upperN << std::endl;
    std::cout << "      INITIAL MODEL WL: " << subUnitWorkingLimit << std::endl;
    std::cout << "              SYMMETRY: " << pModel->getSymmetryString() << std::endl;
    std::cout << "        TOTAL SUBUNITS: " << pModel->getNumberOfSubUnits() << std::endl;

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    unsigned int * ptr = (totalBeadsInSphere != 0) ? &subUnit_indices.front() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::cout << "        CREATING INITIAL RANDOM MODEL " << std::endl;
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit);
    std::set<unsigned int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomInUniverse (0, totalBeadsInSphere-1);


    // determine cvxindices
    updateCVXIndices(subUnitWorkingLimit, subUnit_indices, cvxIndices, pModel);
    std::vector<SubUnit> subunits; // this is 1 less than number of subUnits in Model, parent subUnit is not included
    unsigned int translateTo;// = getUseableNeighborFromSet(&beads_in_use_tree, pModel, cvxIndices[0] );
    for(unsigned int i=0; i<pModel->getNumberOfSubUnits(); i++){  // pick random beads, do rotations and translations with respect to cm
        translateTo = subUnit_indices[randomInUniverse(gen)];
        subunits.emplace_back(SubUnit(i, subUnit_indices[randomIndex(gen)], translateTo, subUnitWorkingLimit, subUnit_indices, pModel));
    }


    std::uniform_int_distribution<unsigned int> randomSubUnitIndex (0, subunits.size()-1);
    std::uniform_int_distribution<unsigned int> randomDegree (0, 10);
    const unsigned int maxAngle = randomDegree.max();
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::uniform_real_distribution<float> randomDegreeIncrement(0,2);

/*
 * calculate energy
 */
    // fill binCount for first time
    calculateModelPrDistributionSymX(&subUnit_indices, &subunits, &binCount, subUnitWorkingLimit, violations, pModel, pData );
    float hlambda = 0.0071;
    float testKL, currentKL = pData->getScore(binCount);

    // intra-subunit connectivity
    double tempIntraContacts, intraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);
    float gammaConstant = 0.01;//0.1*currentKL/((1.0-0.1)*intraContacts);

    // connectivity
    EulerTour eulerTour(subUnit_indices, subUnitWorkingLimit, pModel);
    unsigned int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    // volume
    unsigned int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    char flags[25];
    std::sprintf(flags, "qhull s FA");
    float invTarVolSubUnitmuConstant = 0.0001f/targetSubUnitVolume;
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
    float temp_subunit_volume_energy, current_subunit_volume_energy = std::abs(current_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
    float prior_volume = current_volume;
    //float muConstant = mu*currentKL/((1.0f-mu)*current_volume);


    double current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + current_subunit_volume_energy + gammaConstant*intraContacts;

    double totalViolations = violations;
    beta = 0.001f/(float)pModel->getNumberOfSubUnits(); // has effects on lowerLimit
    current_energy += beta * totalViolations;


    std::cout << " INITIAL => ENERGY :  " << current_energy << endl;
    std::cout << " INITIAL =>   D_KL : " <<  currentKL << endl;
    std::cout << " INITIAL => VOLUME : " << current_volume << endl;
    std::cout << "        VIOLATIONS : " << violations << endl;
    std::cout << " INITIAL        WL : " << subUnitWorkingLimit << endl;

    double test_energy, lowest_energy = current_energy;


    float addRemovePosProb = 0.6537;
    bool updated = false, updateCVX = false, lowestComp = false;
    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int high;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup cop
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());

    for (high=0; high < highTempRounds; high++){

        if (distribution(gen) < addRemovePosProb ) { // add/remove/positional

            startTime = std::clock();

            if (distribution(gen) < 0.376551) { // positional
                std::sprintf(titleText, "POSITIONAL");

                auto itIndex = subUnit_indices.begin() + randomIndex(gen);
                unsigned int swap1 = *itIndex;

                if ( (distribution(gen) < 0.353) || currentNumberOfComponents > 1) { // only move CVX points if true
                    for (unsigned int i = 0; i < subUnitWorkingLimit; i++) {
                        beadToPoint(&points[i*3], pModel->getBead(subUnit_indices[i]));
                        active_indices[i] = subUnit_indices[i];
                    }
                    // needs to be optimized
                    qh_new_qhull(3, subUnitWorkingLimit, points, 0, flags, nullptr, nullptr);
                    vertexT * vertices = qh vertex_list;
                    auto totalV = (unsigned int)(qh num_vertices);

                    // only move CVX hull points
                    std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
                    for (unsigned int v = 0; v < totalV; v++) { //
                        indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                        vertices = vertices->next;
                    }
                    qh_freeqhull(1);

                    std::uniform_int_distribution<unsigned int> randomCVXIndex (0, totalV-1);
                    swap1 = indices_to_check[randomCVXIndex(gen)];
                    itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit, swap1);
                }

                if (std::distance(subUnit_indices.begin(), itIndex) >= subUnitWorkingLimit){
                    exit(0);
                }

                /**
                 * undoing remove may not work correctly
                 * removing from Pr means moving the lattice position to the end in the non-base subunits
                 * undoing is achieved by incrementing workinglimit
                 */
                double tempValue = totalViolations - removeFromPrSymX(swap1, subUnit_indices, subunits, subUnitWorkingLimit, binCount, pModel, pData);

                beads_in_use_tree.erase(swap1);
                eulerTour.removeNode(swap1);

                unsigned int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                while (neighbor == noNeigborIndex || neighbor == swap1){
                    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                }

                // make the swap, sort and update P(r)
                auto pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.end(), neighbor);
                std::iter_swap(itIndex, pSwap2);
                std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit); // bead_indices needs to be sorted

                swapInSubUnits(swap1, neighbor, subunits, pModel);

                temp_violations = tempValue + addToPrSymX( subUnit_indices, &subunits, subUnitWorkingLimit, binCount, pModel, pData);
                testKL = pData->getScore(binCount);

                tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

                test_volume = updateCVXIndices(subUnitWorkingLimit, subUnit_indices, tempCVXIndices, pModel);
                temp_subunit_volume_energy = std::abs(test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;


                tempIntraContacts = intraSubUnitContact(&subunits, &tempCVXIndices, pModel);

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_subunit_volume_energy + beta*temp_violations + gammaConstant*tempIntraContacts;

                if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                    beads_in_use_tree.insert(neighbor);
                    updateCVX = true;
                    updated = true;
                    std::sprintf(addRemoveText, "SWAPPED : %i <=> %i", swap1, neighbor);
                } else {
                    // reverse changes; find index in sorted list and replace
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());

                    eulerTour.removeNode(neighbor);
                    eulerTour.addNode(swap1, pModel);
                    beads_in_use_tree.insert(swap1);

                    // undo swap
                    swapInSubUnits(neighbor, swap1, subunits, pModel);
                    std::sprintf(addRemoveText, "FAILED");
                }

            } else {

                if ((distribution(gen) < 0.5) && subUnitWorkingLimit > lowerN) { // REMOVE beads from sorted list into useable range < deadLimit
                    // randomly swap positions with end of workingLength, could remove CVX Hull Point
                    std::sprintf(titleText, "  REMOVE  ");
                    startTime = std::clock();

                    unsigned int testIndex = subUnit_indices[randomIndex(gen)]; // could include CVX point

                    temp_violations = totalViolations - removeLatticePositionToModelSymX(subUnit_indices, subunits,  binCount, &subUnitWorkingLimit, &testIndex, pModel, pData);
                    testKL = pData->getScore(binCount);

                    // subUnits still has the index to remove, must remove it
                    removeFromSubUnits(testIndex, subunits);

                    test_volume = updateCVXIndices(subUnitWorkingLimit, subUnit_indices, tempCVXIndices, pModel);
                    temp_subunit_volume_energy = std::abs(test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);

                    tempIntraContacts = intraSubUnitContact(&subunits, &tempCVXIndices, pModel);

                    test_energy = testKL + hlambda * (tempNumberOfComponents - 1) * (tempNumberOfComponents - 1) +
                                  temp_subunit_volume_energy + beta*temp_violations + gammaConstant*tempIntraContacts;

                    if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                        beads_in_use_tree.erase(testIndex);
                        updateCVX=true;
                        updated = true;
                        std::sprintf(addRemoveText, "REMOVED ");
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        auto beginIt = subUnit_indices.begin();
                        auto beginBinCount = binCount.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt, &subUnitWorkingLimit, &binCountBackUp,
                                                              &beginBinCount);
                        undoRemoveFromSubUnits(subunits);
                        std::sprintf(addRemoveText, "FAILED");
                    }

                } else { // ADD beads
                    // randomly add from deadLimit
                    std::sprintf(titleText, "    ADD   ");
                    startTime = std::clock();

                    unsigned int testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                    while ( testIndex == noNeigborIndex) { // find a new position
                        testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                    }

                    auto itIndex = std::find(subUnit_indices.begin(), subUnit_indices.end(), testIndex);

                    addLatticPositionToModel(&subUnit_indices, &subUnitWorkingLimit, &itIndex);
                    // add new position to subUnits
                    addToSubUnits(testIndex, subunits, pModel);

                    temp_violations = totalViolations + addToPrSymX(subUnit_indices, &subunits, subUnitWorkingLimit, binCount, pModel, pData);
                    testKL = pData->getScore(binCount);

                    test_volume = updateCVXIndices(subUnitWorkingLimit, subUnit_indices, tempCVXIndices, pModel);
                    temp_subunit_volume_energy = std::abs(test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                    tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                    tempIntraContacts = intraSubUnitContact(&subunits, &tempCVXIndices, pModel);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1)
                                  + temp_subunit_volume_energy + beta*temp_violations + gammaConstant*tempIntraContacts;

                    if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                        beads_in_use_tree.insert(testIndex);
                        updateCVX=true;
                        updated = true;
                        std::sprintf(addRemoveText, "ADDED   ");
                    } else { // undo changes (rejecting)
                        auto beginBinCount = binCount.begin();
                        restoreAddingFromBackUp(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                        currentNumberOfComponents = eulerTour.removeNode(testIndex);

                        undoAddToSubUnits(subunits);
                        std::sprintf(addRemoveText, "FAILED");
                    }
                }
            }



        } else { // rotate/translate
            // randomly pick subunit to alter
            unsigned int subUnitToAlter = randomSubUnitIndex(gen);
            SubUnit * pSubUnit = &subunits[subUnitToAlter];
            startTime = std::clock();

            if (distribution(gen) < 0.367){ //random rotate, pick a pivot point inside molecule

                std::sprintf(titleText, "  ROTATE  ");

                if (distribution(gen) > acceptRate){
                    float anglex = 3.0f*(( (randomDegree(gen) - (randomDegreeIncrement(gen))*maxAngle) ));
                    float angley = 3.0f*(( (randomDegree(gen) - (randomDegreeIncrement(gen))*maxAngle) ));
                    float anglez = 3.0f*(( (randomDegree(gen) - (randomDegreeIncrement(gen))*maxAngle) ));
                    pSubUnit->rotateAt(subUnit_indices[randomIndex(gen)], anglex, angley, anglez, pModel);
                } else {
                    float anglex = 17.0f*(( (randomDegree(gen) - (randomDegreeIncrement(gen))*maxAngle) ));
                    float angley = 17.0f*(( (randomDegree(gen) - (randomDegreeIncrement(gen))*maxAngle) )); // max is 170 in steps of 17
                    float anglez =  8.9f*(( (randomDegree(gen) - (randomDegreeIncrement(gen))*maxAngle) )); // max 89 degrees in steps 8.9
                    pSubUnit->rotateAt(subUnit_indices[randomIndex(gen)], anglex, angley, anglez, pModel);
                }
                // translate all other subUnits
                // calculate new energy
                calculateModelPrDistributionSymX(&subUnit_indices, &subunits, &binCount, subUnitWorkingLimit, temp_violations, pModel, pData );
                testKL = pData->getScore(binCount);

                // randomly add from deadLimit
                tempIntraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);
                tempNumberOfComponents = currentNumberOfComponents;
                test_volume = current_volume;
                temp_subunit_volume_energy = current_subunit_volume_energy;

                test_energy = testKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1)
                              + current_subunit_volume_energy + beta*temp_violations + gammaConstant*tempIntraContacts;

                if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                    updated = true;
                    std::sprintf(addRemoveText, "ROTATED ");
                } else { // undo changes (rejecting)
                    pSubUnit->undoRotateTo();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    std::sprintf(addRemoveText, "FAILED");
                }

            } else { // translate center-of-mass to neighbor that is not in use
                std::sprintf(titleText, " TRANSLATE");

                translateTo = (distribution(gen) > acceptRate) ? subUnit_indices[randomIndex(gen)] : subUnit_indices[randomInUniverse(gen)];

                /*
                 * Any random place not selected
                 */
                // translate the subunit
                pSubUnit->translateToTransformedPosition(translateTo, pModel);
                // calculate new energy
                calculateModelPrDistributionSymX(&subUnit_indices, &subunits, &binCount, subUnitWorkingLimit, temp_violations, pModel, pData );
                testKL = pData->getScore(binCount);
                // randomly add from deadLimit
                tempIntraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);
                tempNumberOfComponents = currentNumberOfComponents;
                test_volume = current_volume;
                temp_subunit_volume_energy = current_subunit_volume_energy;

                test_energy = testKL +
                              hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + current_subunit_volume_energy +
                              beta*temp_violations + gammaConstant*tempIntraContacts;

                if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                    updated = true;
                    std::sprintf(addRemoveText, "TRANSLATED");
                } else { // undo changes (rejecting)
                    pSubUnit->undoTranslateTo(pModel);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    std::sprintf(addRemoveText, "FAILED");
                }
            }
        }



        if (updated){
            currentNumberOfComponents = tempNumberOfComponents;
            currentKL = testKL;
            current_volume = test_volume;
            totalViolations = temp_violations;
            current_energy = test_energy;
            current_subunit_volume_energy = temp_subunit_volume_energy;
            intraContacts = tempIntraContacts;

            acceptRate = inv500slash499*acceptRate+inv500;

            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy

            randomIndex = std::uniform_int_distribution<unsigned int>(0,subUnitWorkingLimit-1); // guaranteed unbiased

            if (updateCVX){ // update CVX indices
                //updateCVXIndices(subUnitWorkingLimit, subUnit_indices, cvxIndices, pModel);
                cvxIndices.resize(tempCVXIndices.size());
                std::copy(tempCVXIndices.begin(), tempCVXIndices.end(), cvxIndices.begin());
            }

            updateCVX=false;
            updated = false;
//            pModel->writeSymXModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, subunits, binCount, "file_" + to_string(counter), this, pData, high, current_volume, this->contactsPerBead);
            counter++;

            if (currentNumberOfComponents == 1){
                addRemovePosProb *= 0.9995;
            }
        } else {
            acceptRate = inv500slash499*acceptRate;
        }

        updateASAModTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);
        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts
        std::printf("*******************             %s                 ******************* \n", titleText);
        std::printf("      TEMP : %-.4E ACCEPT : %.5f ( %.5f ) \n", lowTempStop, acceptRate, addRemovePosProb);
        std::printf("         MAXSTEPS => %i (%i) TIME : %.5f (SECONDS) \n", highTempRounds, high, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
        std::printf("            GRAPH => %i %s \n", currentNumberOfComponents, addRemoveText);
        std::printf("            UPPER => %i LOWER => %i LIMIT : %5i \n", upperN, lowerN, subUnitWorkingLimit);
        std::printf("           VOLUME => %.0f MU*VOL => %.3E\n", current_volume, current_subunit_volume_energy);
        std::printf("       VIOLATIONS => %.3E (%.3E)  BETA => %.4E  \n", totalViolations, intraContacts, beta);
        std::printf("              D_KL : %.4E ENRGY: %.4E \n", currentKL, current_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            /*
             * make copy of lowest configuraiton, must include subUnits
             */
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            lowestWorkingLimit = subUnitWorkingLimit;
            lowestComp = true;

            //gammaConstant = 0.1*currentKL/((1.0-0.1)*intraContacts);

            current_energy = currentKL + current_subunit_volume_energy + beta*totalViolations + gammaConstant*intraContacts;

            lowest_energy = current_energy;
            workingLimitSum += subUnitWorkingLimit;
            sum_x_squared += subUnitWorkingLimit*subUnitWorkingLimit;
            volumeSum += current_volume;
            volumeCount++;
        }
    }


    pModel->updateBeadIndices(lowestWorkingLimit, lowest_subUnit_indices);

    // calculate average volume and standard deviation
    float volumeAverage = workingLimitSum/(float)volumeCount;
    float volumeStdev = sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage;
    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);

    // remove points close to hull
    std::string unname = "lowest_CVX_subunit_" + filenameprefix;
    pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, unname, high);
    unname = "initial_CVX_subunit_" + filenameprefix;
    pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, unname, high);
    pModel->writeSymXModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, subunits, binCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, current_volume, 0);

    pModel->setStartingSet(subUnit_indices);
    pModel->setStartingWorkingLimit(subUnitWorkingLimit);
    pModel->updateSubUnits(subunits);

    std::cout << "AVERAGE # BEADS EST HIGH TEMP SELECTION: " << (int)volumeAverage << " SIGMA: " << (int)volumeStdev << std::endl;
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    std::printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    std::printf("    LOWEST => %d \n", lowestWorkingLimit);
    std::cout << "*******************                                        *******************" << std::endl;


    if (lowestComp) {
        return true;
    } else {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return false;
    }

}


void Anneal::estimateInterSubUnitDistancesX(Model *pModel, Data *pData){

    srand(time(0));
    std::random_device rd;
    char outBuffer[100];

    std::mt19937 gen(rd());

    unsigned int totalTrialsPerRound = 51; // should be proportional to 5^n where n is number of subunits

//    auto topN = (unsigned int)std::ceil(totalTrialsPerRound*0.01);
    unsigned int topN = 6;

    unsigned int totalSubUnits = pModel->getNumberOfSubUnits();

    unsigned int totalDistanceBins = 5;
    float myBinWidth = pData->getDmax()/(float)totalDistanceBins;

    //std::vector<SubUnit> & subunits = pModel->getSubUnits();
    // build distribution potentials for each intersubunit distance
//    for (unsigned int i=0; i < subUnitsInPotential.size(); i++){
//
//        SubUnit & pSubUnit1 = subunits[subUnitsInPotential[i]];
//        unsigned int next = subUnitsInPotential[i] + 1;
//
//        for (unsigned int j=next; j < subunits.size(); j++){
//            cePotentials.push_back(CEPotential(pSubUnit1, subunits[j], myBinWidth, totalDistanceBins));
//        }
//    }
    std::vector<SubUnit> subunits;
    std::vector<SubUnit> bestSubunits(totalSubUnits);
    pModel->copySubUnits(subunits);

    const auto lastIndex = (unsigned int)(subunits.size()-1); // we don't know how first and last connect, but we assume sequential, 1->2, 2->3, etc

    /*
     * cePotentials uses references to the subUnits object.
     * If subUnits is recreated/copied, reference is dead
     * do all neighboring pairwise distances, should be n-1
     */
    cePotentials.clear();
    cePotentials.push_back(CEPotential(subunits[0], subunits[lastIndex], myBinWidth, totalDistanceBins)); // only applying potential between 1st and last subunit
    for(unsigned int i=0; i < lastIndex; i++){
        cePotentials.push_back(CEPotential(subunits[i], subunits[i+1], myBinWidth, totalDistanceBins)); // only applying potential between 1st and last subunit
    }
    // add constraint for first and last
    cePotentials.push_back(CEPotential(subunits[0], subunits[lastIndex], myBinWidth, totalDistanceBins)); // only applying potential between 1st and last subunit
    auto totalInterSubUnitDistancesInPotential = (unsigned int) cePotentials.size();


    unsigned int totalRounds = 3;

    // 10^-3
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50

    // build potentials for each intersubunit distance
    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(pModel->getTotalNumberOfBeadsInUniverse());   // large vector ~1000's

    pModel->copyStartingModelIntoVector(bead_indices);
    std::vector<unsigned int> universe(bead_indices.begin(), bead_indices.end());   // large vector ~1000's

    std::fill(binCount.begin(), binCount.end(), 0);

    double violations;

    unsigned int workingLimit = pModel->getStartingWorkingLimit();
    /*
     * Calculate Center-of-Mass
     */
    vector3 centerOfMass = vector3(0,0,0);
    float invWL = 1.0f/(float)workingLimit;
    for(unsigned int i=0; i<workingLimit; i++){
        centerOfMass += pModel->getBead(bead_indices[i])->getVec();
    }
    centerOfMass *= invWL;

    std::cout << " center of mass " << centerOfMass.x << " " << centerOfMass.y << " " << centerOfMass.z << std::endl;
    for(unsigned int i=0; i < (totalInterSubUnitDistancesInPotential-1); i++){
        CEPotential * pCEPotential = &cePotentials[i];
        // estimate pairwise distance and initialize probabilities
        float distest = (subunits[pCEPotential->getSubUnit1()].translateAndRotatePoint(&centerOfMass) - subunits[pCEPotential->getSubUnit2()].translateAndRotatePoint(&centerOfMass)).length();
        //pCEPotential->seedProbabilities(distest);
        pCEPotential->setTest(distest);
    }

    CEPotential * pCEPotential = &cePotentials[totalInterSubUnitDistancesInPotential-1];
    float distest = (subunits[pCEPotential->getSubUnit1()].translateAndRotatePoint(&centerOfMass) - subunits[pCEPotential->getSubUnit2()].translateAndRotatePoint(&centerOfMass)).length();
    pCEPotential->seedProbabilities(distest);



    // fill binCount for first time
    calculateModelPrDistributionSymX(&bead_indices, &subunits, &binCount, workingLimit, violations, pModel, pData );
    float currentKL = pData->getScore(binCount);

    double coupons = pModel->getNumberOfSubUnits()*(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    float step_limit = (float)(ccmultiple*coupons);

    std::uniform_int_distribution<unsigned int> randomInUniverse (0, pModel->getTotalNumberOfBeadsInUniverse()-1);
    std::uniform_int_distribution<unsigned int> randomInSubUnit (0, workingLimit-1);
    std::uniform_int_distribution<unsigned int> randomSubUnit (0, totalSubUnits-1);

    const float base_energy = estimateInterSubUnitEnergyX(subunits, pModel, pData);
     //(float)(currentKL + beta*totalViolations + zetaConstant*rgPotential + gammaConstant*intraContacts);
    float absolute_best = base_energy;

    std::cout << " STARTING ENERGY " << absolute_best << std::endl;
    std::cout << " ESTIMATING DISTANCES " << std::endl;

    for(unsigned int round=0; round < totalRounds; round++){

        std::vector<CEConfiguration> acceptedConfigurations;
        acceptedConfigurations.resize(topN);

        unsigned int topCount=0;
        unsigned int roundTopCount=0;

        float best_energy = base_energy;

        while(topCount < topN){ // need configurations that are at least as good as start

            subunits.clear();
            subunits.resize(pModel->getNumberOfSubUnits());
            pModel->copySubUnits(subunits);

//            for(unsigned int i=0; i<subunits.size(); i++){ // lack of memory problem
//                SubUnit * pSubUnit = &(subunits[i]);
//                int translateTo = pSubUnit->getTranslateToIndex();
//                translateTo = getTranslateToPointNearCOM(pModel, translateTo);
//                pSubUnit->translateToTransformedPosition(translateTo, pModel);
//                pSubUnit->translateToTransformedPosition(bead_indices[randomInSubUnit(gen)], pModel);
//                pSubUnit->translateToTransformedPosition(bead_indices[randomInUniverse(gen)], pModel);
//                 based on distance, find a bead within the distance?
//            }

            /*
             * create configuration and set which bin is in use
             */
            acceptedConfigurations[topCount] = CEConfiguration(myBinWidth, totalInterSubUnitDistancesInPotential);

            //calculate total distances sum
            float startingMax = totalInterSubUnitDistances(&centerOfMass, subunits, pModel);

            for(int max=0; max<1000; max++){

                //pick random subunit
                SubUnit * pSubUnit = &(subunits[randomSubUnit(gen)]);
                unsigned int oldposition = pSubUnit->getTranslateToIndex();
                pSubUnit->translateToTransformedPosition(bead_indices[randomInUniverse(gen)], pModel);

                float tempMax = totalInterSubUnitDistances(&centerOfMass, subunits, pModel);

                if (tempMax > startingMax){
                    startingMax = tempMax;
                } else {
                    pSubUnit->translateToTransformedPosition(oldposition, pModel);
                }
            }


            pModel->writeSymXModelToFile(currentKL,
                                         workingLimit,
                                         bead_indices,
                                         subunits,
                                         binCount,
                                         "max_" + std::to_string(round),
                                         this,
                                         pData,
                                         1000,
                                         10000,
                                         0);



//            CEConfiguration * pConf = &acceptedConfigurations[topCount];
            // setting bin in use by setConfiguration function in cePotentials
//            SubUnit * pSubUnit1;
//            CEPotential * pCEPotential;
//            bool notset = true;
//            while(notset){
//
//                std::shuffle(universe.begin(), universe.end(), gen);
//                unsigned int stopAt = pModel->getTotalNumberOfBeadsInUniverse()-1;
//                unsigned int randomBase = universe[stopAt];
//                subunits[0].translateToTransformedPosition(randomBase, pModel);
//                cePotentials[totalInterSubUnitDistancesInPotential-1].setConfiguration(); // set configuration for first and last
//
//                std::shuffle(universe.begin(), universe.begin() + stopAt, gen);
//                unsigned int swapTo=stopAt-1;
//
//                for (unsigned int p=0; p < (totalSubUnits-1); p++){
//
//                    pCEPotential = &cePotentials[p];
//                    pConf->populateAssignment(p, pCEPotential->setConfiguration());
//
//                    pSubUnit1 = &subunits[pCEPotential->getSubUnit1()];
//                    vector3 pVecBase = pSubUnit1->translateAndRotatePoint(&centerOfMass);
//
//                    unsigned int randomlocation = universe[0];
//                    bool foundIt = false;
//                    for(unsigned int u=0; u< (swapTo+1); u++){
//                        randomlocation = universe[u];
//                        float dis = (pVecBase - pModel->getBead(randomlocation)->getVec()).length();
//                        if (pCEPotential->keepIt(dis)){
//                            std::iter_swap(universe.begin()+u, universe.begin()+swapTo);
//                            swapTo--;
//                            foundIt = true;
//                            break;
//                        }
//                    } // translate the second subUnit
//
//                    if (!foundIt){ // try again
//                        break;
//                    }
//
//                    subunits[pCEPotential->getSubUnit2()].translateToTransformedPosition(randomlocation, pModel);
//                }
//
//                float distest = (subunits[0].translateAndRotatePoint(&centerOfMass) - subunits[lastIndex].translateAndRotatePoint(&centerOfMass)).length();
//
//                if (cePotentials[totalInterSubUnitDistancesInPotential-1].keepIt(distest )) {
//                    notset = false;
//                }
//            }


            float current_energy = symXRotateTranslate(round, topCount, step_limit, subunits, cePotentials, pData, pModel, true, "rt_search");

            //float current_energy = estimateInterSubUnitEnergyX(subunits, pModel, pData);
            //std::cout << "FOUND suitable configuration " << current_energy << " " << base_energy << std::endl;


            //double floored = floor(current_energy);
//            if (current_energy < base_energy && round < 2) {
//                // update acceptedConfigurations
//                acceptedConfigurations[topCount].setScore(current_energy);
                std::sprintf(outBuffer, "%d INITIAL TOP N => %d ( %d) SCORE %.3E  BASE : %.2E ( %.2E )\n", round, topCount, topCount, current_energy, base_energy, currentKL);
                logToCeFile(outBuffer);
//                topCount++;
//                roundTopCount = 0;
//            } else if (round > 1) {
                acceptedConfigurations[topCount].setScore(current_energy);
                topCount++;
//            }

            if (current_energy < best_energy){
                best_energy = current_energy;
                std::copy(subunits.begin(), subunits.end(), bestSubunits.begin());
            }

            roundTopCount++;
        }


        std::sort(acceptedConfigurations.begin(), acceptedConfigurations.end(), [](const CEConfiguration & lhs, const CEConfiguration & rhs) {
            return lhs.getScore() < rhs.getScore();
        });

        float last = acceptedConfigurations[topN-1].getScore();


        for(unsigned int trial=topN; trial< totalTrialsPerRound; trial++){
            // generate configuration
            subunits.clear();
            subunits.resize(pModel->getNumberOfSubUnits());
            pModel->copySubUnits(subunits);

//            for(unsigned int i=0; i<subunits.size(); i++){
//                SubUnit * pSubUnit = &(subunits[i]);
//                int translateTo = pSubUnit->getTranslateToIndex();
//                translateTo = getTranslateToPointNearCOM(pModel, translateTo);
//                pSubUnit->translateToTransformedPosition(translateTo, pModel);
//                pSubUnit->translateToTransformedPosition(bead_indices[randomInSubUnit(gen)], pModel);
//                pSubUnit->translateToTransformedPosition(bead_indices[randomInUniverse(gen)], pModel);
//            }


            //calculate total distances sum
            float startingMax = totalInterSubUnitDistances(&centerOfMass, subunits, pModel);

            for(int max=0; max<1000; max++){

                //pick random subunit
                SubUnit * pSubUnit = &(subunits[randomSubUnit(gen)]);
                unsigned int oldposition = pSubUnit->getTranslateToIndex();
                pSubUnit->translateToTransformedPosition(bead_indices[randomInUniverse(gen)], pModel);

                float tempMax = totalInterSubUnitDistances(&centerOfMass, subunits, pModel);

                if (tempMax > startingMax){
                    startingMax = tempMax;
                } else {
                    pSubUnit->translateToTransformedPosition(oldposition, pModel);
                }
            }



            // copy Starting_Set from initial model
            // end of steps
            // update acceptedConfigruations
            CEConfiguration ceConf = CEConfiguration(myBinWidth, totalInterSubUnitDistancesInPotential);
            for (unsigned int p=0; p<totalInterSubUnitDistancesInPotential; p++){
                ceConf.populateAssignment(p, cePotentials[p].setConfiguration());
            }

            //float current_energy = estimateInterSubUnitEnergyX(subunits, pModel, pData);
            float current_energy = symXRotateTranslate(round, trial, step_limit, subunits, cePotentials, pData, pModel, false, "rt_search");

            if (current_energy < last){
                acceptedConfigurations[topN-1] = ceConf;
                acceptedConfigurations[topN-1].setScore(current_energy);

                std::sort(acceptedConfigurations.begin(), acceptedConfigurations.end(), [](const CEConfiguration & lhs, const CEConfiguration & rhs) {
                    return lhs.getScore() < rhs.getScore();
                });

                std::sprintf(outBuffer, "%d     NEW TOP N => %d       SCORE %.3E  <  %.3E [ %.3E ] \n", round,  trial, current_energy, last, best_energy);
                logToCeFile(outBuffer);

                last = acceptedConfigurations[topN-1].getScore();
            }

            if (current_energy < best_energy){
                best_energy = current_energy;
                absolute_best = best_energy;
                std::copy(subunits.begin(), subunits.end(), bestSubunits.begin());
                pModel->writeSymXModelToFile(current_energy,
                                             workingLimit,
                                             bead_indices,
                                             subunits,
                                             binCount,
                                             "best_CVX_sym_" + std::to_string(round),
                                             this,
                                             pData,
                                             1000,
                                             100000,
                                             0);
            }
        }

        /*
         * update probabilities based on topN
         * totalInPotential
         */
        for(unsigned int i=0; i < totalInterSubUnitDistancesInPotential; i++){ // go through each SubUnit-to-SubUnit potential
            //for each bin determine distribution count
            float totalCounts = 0;
            std::vector<unsigned int> counts(totalBins);
            std::fill(counts.begin(), counts.end(), 0);

            for(unsigned int top=0; top<topN; top++){
                CEConfiguration * pConf = &acceptedConfigurations[top]; // bin assigments for each subunit-to-subunit distance
                std::sprintf(outBuffer, " TOP N => %d  SCORE %.3E  BIN %d (BASE : %.2E)\n", top, pConf->getScore(), pConf->getBinAssignment(i), best_energy);
                logToCeFile(outBuffer);
                ++counts[pConf->getBinAssignment(i)];
                ++totalCounts;
            }

            // update each bin in CEPotential
            CEPotential * pCEP = &cePotentials[i];
            pCEP->updateProbabilities(counts, totalCounts);
            pCEP->printDistribution(round);
        }
    }

    // which bin is largest
//    for(unsigned int i=0; i<totalInterSubUnitDistancesInPotential; i++){
//        cePotentials[i].setMostProbableDistance();
//        //cePotentials[i].setTest(81);
//    }

    cePotentials[totalInterSubUnitDistancesInPotential-1].setMostProbableDistance();
    float distance = cePotentials[totalInterSubUnitDistancesInPotential-1].getMostProbableDistance();
    cePotentials.clear();

    cePotentials.push_back(CEPotential(subunits[0], subunits[lastIndex], myBinWidth, totalDistanceBins)); // only applying potential between 1st and last subunit
    cePotentials[0].setTest(distance);
    // score must be better that best
//    while (symXRotateTranslate(1234, 1234, bestSubunits, cePotentials, pData, pModel, true) > absoluteBest){
//        for(int i=0; i<subunits.size(); i++){
//            SubUnit * pSubUnit = &(bestSubunits[i]);
//            int translateTo = pSubUnit->getTranslateToIndex();
//            translateTo = getTranslateToPointNearCOM(pModel, translateTo);
//            pSubUnit->translateToTransformedPosition(translateTo, pModel);
//            //pSubUnit->translateToTransformedPosition(bead_indices[randomInSubUnit(gen)], pModel);
//        }
//    }
//

//    pModel->updateSubUnits(subunits);
    std::cout << " UPDATING SUBUNITS " << std::endl;

}


/**
 * The refinement alternates between :
 *     1. rotating and translating the subunits
 *     2. add/remove/positional refinement of a subunit with no rotation or translation
 *
 *     The acceptance rate starts at 0.44 for each round is gradually reduced to 0.1 in the final round
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineHomogenousBodyASAHybridSymmetryX(Model *pModel, Data *pData, std::string nameTo) {

    std::cout << " refineHomogenousBodyASAHybridSymmetryX " << std::endl;

    char outBuffer[80];
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY IDENTICAL SUBUNITS" << std::endl;

    std::string tempName = nameTo + "_refinementSymX";
    Logger logger(tempName);
    logger.setDataFilename(pData->getDataFilename());
    char stringbuffer[90];


    // set deadLimit of the selected set
    // convert distances in Search Sphere to ShannonBin membership
    maxbin = pModel->populateBins(pData);maxbin = pData->convertToBin(pModel->getDiameterOfUniverse()*3);
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    pData->setScoringFunction(maxbin);

    // std::vector<SubUnit> & subunits = pModel->getSubUnits();
    // cycle between rounds of refinment of the subUnit and rotation/translations
    std::vector<SubUnit> subunits;
    std::vector<SubUnit> backUpsubunits;
    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    std::uniform_int_distribution<unsigned int> randomInUniverse (0, totalBeadsInSphere-1);

    const unsigned int totalRounds = 41;
    float startAcc = 0.44, finalAcc = 0.1, factor = std::pow(finalAcc/startAcc, 1.0f/(float)totalRounds);
    float bestARPEnergy=0, bestRotaTrans=0;

    for(unsigned int alternatingRounds=0; alternatingRounds < totalRounds; alternatingRounds++) {
        /*
         * lower acceptance rate per round
         */
        auto accRate = (float)(startAcc*std::pow(factor, alternatingRounds));
        this->asaAcceptanceRate = accRate;
        complementASAAcceptanceRate = 1.0f - accRate;
        intASAAcceptanceRate = 1000*accRate;
        intComplementASAAcceptanceRate = 1000*complementASAAcceptanceRate;

        pModel->copySubUnits(subunits);
        pModel->copySubUnits(backUpsubunits);

        if (true){ // do rotate/translate
//        if (alternatingRounds % 2 == 0){ // do rotate/translate
            std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
            // copy Starting_Set from initial model
            pModel->copyStartingModelIntoVector(bead_indices);
            unsigned int workingLimit = pModel->getStartingWorkingLimit();

            // make small random displacements of starting model?
//            for(unsigned int i=0; i<subunits.size(); i++){
//                SubUnit * pSubUnit = &(subunits[i]);
//                unsigned int translateTo = pSubUnit->getTranslateToIndex();
//                translateTo = (accRate > 0.17) ?  getTranslateToPointNearCOM(pModel, translateTo) : getUseableNeighborForCOM(pModel, translateTo);
//                //translateTo = (accRate > 0.20) ?  bead_indices[randomInUniverse(gen)] : getUseableNeighborForCOM(pModel, translateTo);
//                pSubUnit->translateToTransformedPosition(translateTo, pModel);
//            }

            auto coupons = (float)(pModel->getNumberOfSubUnits()*(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5));
            float step_limit = ccmultiple*coupons;

            float rotaTransE = symXRotateTranslate(alternatingRounds, alternatingRounds, step_limit, subunits, cePotentials, pData, pModel, true, "rottrans_round");

            if (alternatingRounds < 2){
                bestRotaTrans = rotaTransE;
                pModel->updateSubUnits(subunits);
            } else if (rotaTransE < bestRotaTrans){
                pModel->updateSubUnits(subunits);

                std::string tempString = "";
                std::snprintf(stringbuffer, sizeof(stringbuffer), "ACCEPT RT %4d : %.7f < %.7f [ %.3f ]\n", alternatingRounds, rotaTransE, bestRotaTrans, accRate);
                logger.logIt(std::string(stringbuffer));

                bestRotaTrans = rotaTransE;
                logToCeFile(outBuffer);
            } else {
                std::string tempString = "";
                std::snprintf(stringbuffer, sizeof(stringbuffer), "FAILED RT %4d : %.7f > %.7f [ %.3f ]\n", alternatingRounds, rotaTransE, bestRotaTrans, accRate);
                logger.logIt(std::string(stringbuffer));
            }

            std::cout << "FNISHED ROTATE/TRANSLATE " << alternatingRounds << std::endl;

        } else { // even
            /*
             * make adjustments to subUnit without rotation/translation
             */
            // number of shannon bins for the model is calculated over the Universe (not the data)
            std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
            std::vector<unsigned int> workingBinCount(maxbin);   // smallish vector, typically < 50
            std::vector<unsigned int> testBinCount(maxbin);      // smallish vector, typically < 50
            std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

            std::cout << "   TOTAL EXP N_S BINS : " << totalBins << std::endl;
            std::cout << "   MAX MODEL N_S BINS : " << maxbin << std::endl;
            std::cout << "             BINWIDTH : " << pData->getBinWidth() << std::endl; // lattice should be limited by binwidth

            auto beginBinCount = binCount.begin();
            auto endBinCount = binCount.end();

            std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
            std::vector<unsigned int> backUpState(totalBeadsInSphere);
            std::vector<unsigned int> cvxIndices;
            std::vector<unsigned int> tempCVXIndices;

            // copy Starting_Set from initial model
            pModel->copyStartingModelIntoVector(bead_indices);;

            unsigned int workingLimit = pModel->getStartingWorkingLimit();
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
            std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

            // coupon collector's problem
            auto coupons = (unsigned int)std::lround(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
            unsigned int updateCount = 3*coupons;
            auto step_limit = (float) updateCount;

            double violations, temp_violations;
            calculateModelPrDistributionSymX(&bead_indices, &subunits, &binCount, workingLimit, violations, pModel, pData );
            float testKL, currentKL = pData->getScore(binCount);

//            char flags[] = "qhull FA"; // CVX HULL STUFF
            float temp_volume, current_volume = updateCVXIndices(workingLimit, bead_indices, cvxIndices, pModel);

            std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;
            std::uniform_real_distribution<float> distribution(0.0,1.0);


            EulerTour eulerTour(bead_indices, workingLimit, pModel);
            unsigned int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

            bool isUpdated = false, updateCVX = false;
            float acceptRate = 0.5;
            const float inv500 = 1.0f/500.0f, inv500slash499 = 499.0f/500.0f;
            unsigned int original, counter=1, failures=0, updated=0;
            double lowTempStop = (double)highTempStartForCooling, inv_kb_temp = 1.0/lowTempStop ;

            float this_energy, lowestKL = currentKL;
            char addRemoveText[50];
            /*
             * estimates constants, eta and mu
             */
            float muConstant = 0.000001f/current_volume;//mu*currentKL/((1.0-mu)*current_volume); // volume
            const float targetSubUnitVolume = current_volume;

            float temp_subunit_volume_energy, current_subunit_volume_energy = 0;

            double tempIntraContacts, intraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel); // subunit to subunit contacts
            float gammaConstant = 0.001; // between subUnit contacts
            // Add components to estimate total energy ( not including totalContact energy

            float current_energy = currentKL +
                                   lambda*connectivityPotential(currentNumberOfComponents) +
                                   gammaConstant*intraContacts + current_subunit_volume_energy;

            double totalViolations = violations;
            current_energy += beta*totalViolations;

            float lowest_total_energy = current_energy;

            if (alternatingRounds < 2){
                bestARPEnergy = current_energy;
            }

            //float initial_energy = current_energy;
//            std::string tempString = "";
//            std::snprintf(stringbuffer, sizeof(stringbuffer), "    VIOLATIONS POTENTIAL : %.3E", violations);
//            logger.logIt(std::string(stringbuffer));
//            std::snprintf(stringbuffer, sizeof(stringbuffer),"           CURRENT ENRGY : %.3E", current_energy);
//            logger.logIt(std::string(stringbuffer));

            std::clock_t startTime;
            std::uniform_int_distribution<unsigned int> randomIndex (0,workingLimit-1); // guaranteed unbiased
            std::uniform_int_distribution<unsigned int> randomSubUnitIndex (0, subunits.size()-1);

            std::vector<unsigned int> selections(workingLimit);
            for(unsigned int i=0; i < workingLimit; i++){
                selections[i] = i; // some distances will exceed dmax
            }

            unsigned int numberOfCoolingTempSteps=0, noNeighborLimit = pModel->getNeighborLimit();

            for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

                if (distribution(gen) < percentAddRemove ){ //add or remove bead within working Set (exclude deadzone)

                    std::cout << "______________________________________________________________________________" << std::endl;
                    std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
                    // additional points to expand deadlimit will occur via enlarging CVX Hull
                    startTime = std::clock();
                    if (distribution(gen) < 0.5 || workingLimit < 10){ // ADD BEAD?
                        std::cout << "*******************                  ADD                   *******************" << std::endl;

                        original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                        while(original == noNeighborLimit){
                            original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                        }

                        // rather than try to find a bead to add, see if the one random one is acceptable
                        auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);
                        // check if available neighbor can be added
                        addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);
                        // add new position to subUnits
                        addToSubUnits(original, subunits, pModel);

                        beads_in_use_tree.insert(original);

                        temp_violations = totalViolations + addToPrSymX(bead_indices, &subunits, workingLimit, binCount, pModel, pData);
                        testKL = pData->getScore(binCount);

                        tempNumberOfComponents = eulerTour.addNode(original, pModel);

                        temp_volume = updateCVXIndices(workingLimit, bead_indices, tempCVXIndices, pModel);
                        temp_subunit_volume_energy = muConstant*std::abs(temp_volume-targetSubUnitVolume);

                        tempIntraContacts = intraSubUnitContact(&subunits, &tempCVXIndices, pModel);

                        this_energy = testKL +
                                      lambda*connectivityPotential(tempNumberOfComponents) +
                                      beta*temp_violations +
                                      gammaConstant*tempIntraContacts + temp_subunit_volume_energy;


                        if ((this_energy < current_energy) || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                            isUpdated = true;
                            updateCVX = true;
                            std::sprintf(addRemoveText, "     ADD => %i", 1);
                        } else { // undo changes (rejecting)
                            beads_in_use_tree.erase(original);
                            eulerTour.removeNode(original);
                            sprintf(addRemoveText, "     ADD => %i", 0);
                            beginBinCount = binCount.begin();
                            restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                            undoAddToSubUnits(subunits);
                        }

                    } else { // REMOVE BEADS?

                        std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                        // test for deletion
                        std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                        original = bead_indices[selections[0]];
                        // test for deletion
                        for(unsigned int i=1; i<workingLimit; i++){
                            tempNumberOfComponents = eulerTour.removeNode(original);
                            if (tempNumberOfComponents <= currentNumberOfComponents){
                                break;
                            } else {
                                eulerTour.addNode(original, pModel);
                                original = bead_indices[selections[i]]; // potential to select the same thing twice
                            }
                        }

                        temp_violations = totalViolations - removeLatticePositionToModelSymX(bead_indices, subunits,  binCount, &workingLimit, &original, pModel, pData);

                        beads_in_use_tree.erase(original);
                        testKL = pData->getScore(binCount);
                        tempNumberOfComponents = 1;//eulerTour.removeNode(original);

                        temp_volume = updateCVXIndices(workingLimit, bead_indices, tempCVXIndices, pModel);
                        temp_subunit_volume_energy = muConstant*std::abs(temp_volume-targetSubUnitVolume);

                        tempIntraContacts = intraSubUnitContact(&subunits, &tempCVXIndices, pModel);

                        this_energy = testKL +
                                      lambda*connectivityPotential(tempNumberOfComponents) +
                                      beta*temp_violations +
                                      gammaConstant*tempIntraContacts + temp_subunit_volume_energy;

                        if ((this_energy < current_energy) || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                            isUpdated = true;
                            updateCVX = true;
                            removeFromSubUnits(original, subunits);
                            sprintf(addRemoveText, " REMOVED => %i", original);
                        } else { // undo changes and move to next bead (rejecting)
                            beads_in_use_tree.insert(original);
                            auto beginIt = bead_indices.begin();
                            restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                            eulerTour.addNode(original, pModel);
                            sprintf(addRemoveText, "  FAILED => %i", original);
                        }
                    }

                } else { // positional refinement
                    std::cout << "______________________________________________________________________________" << std::endl;
                    std::cout << "*******************               POSITIONAL               *******************" << std::endl;
                    std::cout << "*******************                                        *******************" << std::endl;
                    // float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
                    // find bead to swap in active set
                    // only refine a single position with 2 or less contacts
                    startTime = std::clock();
                    // randomly select an index to move
                    // select only node I can move?
                    std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                    unsigned int position = selections[0];
                    unsigned int swap1 = bead_indices[ position ];

                    for(unsigned int i=1; i<workingLimit; i++){
                        if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            break;
                        } else { // if last node, will break out
                            eulerTour.addNode(swap1, pModel);
                            position = selections[i];
                            swap1 = bead_indices[ position ]; // what happens if I fail all the way?
                        }
                    }

                    // swap
                    auto itIndex = bead_indices.begin() + position;
                    // remove selected index from P(r)
                    double tempValue = totalViolations - removeFromPrSymX(swap1, bead_indices, subunits, workingLimit, binCount, pModel, pData);

                    beads_in_use_tree.erase(swap1);

                    // find better position
                    unsigned int originalSwap2Value;
                    for(unsigned int i=0;i<workingLimit; i++){
                        unsigned int newPosition = bead_indices[selections[i]];
                        if (newPosition != swap1){
                            originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, newPosition);
                            if (originalSwap2Value != noNeighborLimit && originalSwap2Value != swap1){
                                break;
                            }
                        }
                    }

                    // get available neighbor position
                    auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
                    std::iter_swap(itIndex, pSwap2);
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                    swapInSubUnits(swap1, originalSwap2Value, subunits, pModel);
                    //addToSubUnits(originalSwap2Value, subunits, pModel);
                    temp_violations = tempValue + addToPrSymX(bead_indices, &subunits,  workingLimit, binCount, pModel, pData);

                    // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                    testKL = pData->getScore(binCount);

                    beads_in_use_tree.insert(originalSwap2Value);

                    tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);

                    temp_volume = updateCVXIndices(workingLimit, bead_indices, tempCVXIndices, pModel);
                    temp_subunit_volume_energy = muConstant*std::abs(temp_volume-targetSubUnitVolume);

                    tempIntraContacts = intraSubUnitContact(&subunits, &tempCVXIndices, pModel);

                    this_energy = testKL +
                                  lambda*connectivityPotential(tempNumberOfComponents) +
                                  beta*temp_violations +
                                  gammaConstant*tempIntraContacts +
                                  temp_subunit_volume_energy;

                    if ((this_energy < current_energy) || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                        isUpdated = true;
                        updateCVX = true;
                        std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                    } else {
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                        beads_in_use_tree.erase(originalSwap2Value);
                        beads_in_use_tree.insert(swap1); // add back swap one to tree

                        std::sprintf(addRemoveText, "      FAILED => %i", swap1);

                        eulerTour.removeNode(originalSwap2Value);
                        currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                        swapInSubUnits(originalSwap2Value, swap1, subunits, pModel);
                    }
                } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

                printf("       TEMP : %-.4E MAXSTEPS => %.0f (%4i) ROUND => %d \n", lowTempStop, step_limit, numberOfCoolingTempSteps, alternatingRounds);
                printf("     ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
                printf("       TIME : %.5f (SECONDS)  %s\n", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC), addRemoveText);
                printf("      LIMIT : %5i\n", workingLimit);
                printf(" VIOLATIONS => %.3E  INTRA => %.3E COMP => %d \n", totalViolations, intraContacts, currentNumberOfComponents);
                printf("       D_KL => %-5.4E ( %.4E ) ENRGY : %.4E\n", currentKL, bestARPEnergy, current_energy);


//        if (!symXBinTest(&bead_indices, &subunits, &binCount, workingLimit, pModel, pData )){
//            std::cout << " BIN MISMATCH " << std::endl;
//            exit(0);
//        }

                //cout << "______________________________________________________________________________" << endl;
                //cout << "*******************                 TEST                   *******************" << endl;
                //cout << "*******************              -----------               *******************" << endl;
//        float testKL1 = calculateKLEnergySymmetry(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//        if (currentKL != testKL1){// || checkForRepeats(bead_indices)){
//            testKL = pData->calculateKLDivergence(binCount);
//            std::cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << " -- " << testKL << std::endl;
//            return "stopped";
//        }

//        checkSubUnits(workingLimit, subunits, "After Change");

                // update run time parameters
//                logger.update(lowTempStop, currentKL, workingLimit);

                // Adaptive simulated annealing part
                if (isUpdated){
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    intraContacts = tempIntraContacts;
                    current_volume = temp_volume;
                    totalViolations = temp_violations;
                    current_subunit_volume_energy = temp_subunit_volume_energy;

                    acceptRate = inv500slash499*acceptRate+inv500;
                    isUpdated = false;
                    failures=0;
                    updated++;

                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased
                    //randomDeadIndex = std::uniform_int_distribution<int>(workingLimit, deadLimit-1);
                    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
                    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

                    selections.resize(workingLimit);
                    unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
                    for(unsigned int i=0; i < workingLimit; i++){
                        *(pSelections+i) = i; // some distances will exceed dmax
                    }

                    if (updateCVX){ // update CVX indices
                        //updateCVXIndices(workingLimit, bead_indices, cvxIndices, pModel);
                        cvxIndices.resize(tempCVXIndices.size());
                        std::copy(tempCVXIndices.begin(), tempCVXIndices.end(), cvxIndices.begin());
                        updateCVX = false;
                    }

                } else {
                    acceptRate = inv500slash499*acceptRate;
                    failures++;
                }

                updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);

                //update for running average
                counter++;
            } // end of steps

            float tempAverageContacts=0.0;

            // calculate average_number_of_contacts
//            for (unsigned int i=0; i< workingLimit; i++){
//                unsigned int temp = bead_indices[i];
//                tempAverageContacts += numberOfContacts(temp, &bead_indices, workingLimit, pModel);
//            }
//            tempAverageContacts *= 1.0/(double)workingLimit;

            if (current_energy < bestARPEnergy){
                std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "best_subunit_round_" + std::to_string(alternatingRounds), this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
                pModel->writeSymXModelToFile(currentKL, workingLimit, bead_indices, subunits, binCount, "best_round_" + std::to_string(alternatingRounds), this, pData, numberOfCoolingTempSteps, current_volume, 0);

                pModel->updateBeadIndices(workingLimit, bead_indices);
                pModel->setStartingSet(bead_indices);
                pModel->setStartingWorkingLimit(workingLimit);

                pModel->updateSubUnits(subunits);

                std::string tempString = "";
                std::snprintf(stringbuffer, sizeof(stringbuffer), "ACCEPT AP %4d : %.7f < %.7f [ %.3f ]\n", alternatingRounds, current_energy, bestARPEnergy, accRate);
                logger.logIt(std::string(stringbuffer));


                logToCeFile(outBuffer);
                bestARPEnergy = current_energy;
            } else {
                std::copy(backUpsubunits.begin(), backUpsubunits.end(), subunits.begin());

                std::string tempString = "";
                std::snprintf(stringbuffer, sizeof(stringbuffer), "REJECT AP %4d : %.5f > %.5f [ %.3f ]\n", alternatingRounds, current_energy, bestARPEnergy, accRate);
                logger.logIt(std::string(stringbuffer));

                logToCeFile(outBuffer);
            }

        } // end of add remove positional

        // write out best model?


    }

    logger.writeToFile();

    return "stop";
}


float Anneal::symXRotateTranslate(unsigned int round, unsigned int index, float step_limit, std::vector<SubUnit> & subunits, std::vector<CEPotential> & cePotentials, Data * pData,  Model * pModel, bool printIt, std::string name) {

    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    char addRemoveText[50];

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    unsigned int workingLimit = pModel->getStartingWorkingLimit();
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    std::uniform_int_distribution<unsigned int> randomInUniverse (0, totalBeadsInSphere-1);
    std::uniform_int_distribution<unsigned int> randomDegree (0, 10);
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> sign(0,1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> stepDegree(0,2); // guaranteed unbiased
    float maxAngle = randomDegree.max();

    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    const unsigned int totalSubUnits = pModel->getNumberOfSubUnits();
    std::uniform_int_distribution<unsigned int> randomSubUnitIndex (0, totalSubUnits-1);
    double violations, temp_violations;

    // build potentials for each intersubunit distance
    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    std::vector<unsigned int> cvxIndices;
    std::vector<unsigned int> tempCVXIndices;

    pModel->copyStartingModelIntoVector(bead_indices);
    std::fill(binCount.begin(), binCount.end(), 0);
    // reset iterators to internal bead_indices
    /*
     * should I randomize?
     */
//    for(int i=0; i<subunits.size(); i++){ // randomize subUnits
//        SubUnit * pSubUnit = &(subunits[i]);
//        pSubUnit->translateToTransformedPosition(bead_indices[randomInUniverse(gen)], pModel);
//        //pSubUnit->translateToTransformedPosition(bead_indices[randomIndex(gen)], pModel);
//    }

    float current_volume = updateCVXIndices(workingLimit, bead_indices, cvxIndices, pModel);
    double tempIntraContacts, intraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);

    bool isUpdated = false;
    float acceptRate = 0.5, inv500 = 1.0f/500.0f, inv500slash499 = 499.0f/500.0f;
    unsigned int failures=0, updated=0;
    double lowTempStop = 0.001;
    double inv_kb_temp = 1.0/lowTempStop;

    float gammaConstant = 0.001; // contacts between subUnit contacts

    float this_energy;

    // initialize for new round
    // fill binCount for first time
    calculateModelPrDistributionSymX(&bead_indices, &subunits, &binCount, workingLimit, violations, pModel, pData );
    float testKL, currentKL = pData->getScore(binCount);


    const float expRg = pData->getRg();
    double zeta = 0.53721;
    //float diff = expRg - pData->calculateRgFromExperimentalDistribution(binCount);
    double tempRgPotential, rgPotential = moments(&subunits, &bead_indices, workingLimit, pData->getRg(), pModel);
    float zetaConstant = 0.0001;//(rgPotential == 0) ? 0.0001f : (float)(zeta*currentKL/((1.0f-zeta)*rgPotential)); // rg potential

    float current_energy = currentKL +
                           gammaConstant*intraContacts +
                           zetaConstant*rgPotential;

    double totalViolations = violations;
    beta=0.001;
    current_energy += beta*totalViolations;

    /*
     * Calculate Center-of-Mass
     */
    vector3 cofm = vector3(0,0,0);
    float invWL = 1.0f/(float)workingLimit;
    for(unsigned int i=0; i<workingLimit; i++){
        cofm += pModel->getBead(bead_indices[i])->getVec();
    }
    cofm *= invWL;

    float cofmdist = 1000;
    unsigned int cofmIndex=0;
    for(unsigned int i=0; i<totalBeadsInSphere; i++){
        float dis = (cofm - pModel->getBead(bead_indices[i])->getVec()).length();
        if (dis < cofmdist){
            cofmdist = dis;
            cofmIndex = i;
        }
    }


    float tempCOMEnergy, comEnergy = calculateCOMPotential(&cofm, subunits, cePotentials);
    float epsilonConstant = 0.01f;// (comEnergy == 0) ? 0.001 : epsilon*currentKL/((1.0-epsilon)*comEnergy);

    current_energy += epsilonConstant*comEnergy;
    float lowest_energy = current_energy;

    int wholeCount=0;
    auto sixtyfivepercent = (unsigned int)(step_limit*0.65);

    pModel->writeSymXModelToFile(currentKL,
                                 workingLimit,
                                 bead_indices,
                                 subunits,
                                 binCount,
                                 "before_" + std::to_string(round),
                                 this,
                                 pData,
                                 1000,
                                 current_volume,
                                 0);


    //perform minimization
    for(unsigned int numberOfCoolingTempSteps=0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {
        // randomly pick subunit to alter
        unsigned int subUnitToAlter = randomSubUnitIndex(gen);
        SubUnit * pSubUnit = &(subunits[subUnitToAlter]);
        std::cout << "______________________________________________________________________________" << std::endl;

        if (distribution(gen) < 0.312167){ //random rotate
            std::cout << "*******************               ROTATIONAL               *******************" << std::endl;

            if (numberOfCoolingTempSteps > sixtyfivepercent){
//            float anglex = 1.75f*((float)( (randomDegree(gen) - (rand()%2)*maxAngle) )); // 30 degrees
//            float angley = 1.75f*((float)( (randomDegree(gen) - (rand()%2)*maxAngle) ));
//            float anglez = 1.75f*((float)( (randomDegree(gen) - (rand()%2)*maxAngle) ));
                float plusminus = (sign(gen) == 1) ? 1.0f : -1.0f;
                float fraction = (sign(gen) == 1) ? 1.0f : 0.5f;

                float anglex = 10.89f*(plusminus)*(stepDegree(gen))*fraction; // 30 degrees
                float angley = 10.89f*(plusminus)*(stepDegree(gen))*fraction;
                float anglez = 10.89f*(plusminus)*(stepDegree(gen))*fraction;
                //pSubUnit->rotateAt(bead_indices[randomIndex(gen)], anglex, angley, anglez, pModel);
                pSubUnit->rotateAt(bead_indices[cofmIndex], anglex, angley, anglez, pModel);
            } else {
                float anglex = 17.0f*(( (randomDegree(gen) - (rand()%2)*maxAngle) ));
                float angley =  8.9f*(( (randomDegree(gen) - (rand()%2)*maxAngle) ));  // max is 170 in steps of 17
                float anglez = 17.0f*(( (randomDegree(gen) - (rand()%2)*maxAngle) )); // max 89 degrees in steps 8.9
                pSubUnit->rotateAt(bead_indices[cofmIndex], anglex, angley, anglez, pModel);
            }

            // calculate new energy
            calculateModelPrDistributionSymX(&bead_indices, &subunits, &binCount, workingLimit, temp_violations, pModel, pData );
            testKL = pData->getScore(binCount);

            // contact energy within subUnit is the same so does not get included in the sum
            tempIntraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);
            tempRgPotential = moments(&subunits, &bead_indices, workingLimit, pData->getRg(), pModel);
            tempCOMEnergy = calculateCOMPotential(&cofm, subunits, cePotentials);

            this_energy = testKL +
                          beta*temp_violations +
                          gammaConstant*tempIntraContacts +
                          zetaConstant*tempRgPotential +
                          epsilonConstant*tempCOMEnergy;

            if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                isUpdated=true;
                std::sprintf(addRemoveText, "  ROTATED ");
            } else { // undo changes (rejecting)
                pSubUnit->undoRotateTo();
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::sprintf(addRemoveText, " NO ROTATED");
            }

        } else { // translate center-of-mass to neighbor that is not in use
            std::cout << "*******************             TRANSLATIONAL              *******************" << std::endl;
            // translate the subunit, either pick a place within the subUnit (small translation) or anywhere
            // based on
            bool single = false;
            bool whole = false;
            if (distribution(gen) > acceptRate){ // grab neighbor to current translation point
                unsigned int translateTo = pSubUnit->getTranslateToIndex();
                translateTo = getUseableNeighborForCOM(pModel, translateTo);
                pSubUnit->translateToTransformedPosition(translateTo, pModel);
                single = true;
            } else {
                // pSubUnit->translateToTransformedPosition(bead_indices[randomInUniverse(gen)], pModel);
//
                for(unsigned int i=0; i<subunits.size(); i++){
                    pSubUnit = &(subunits[i]);
                    unsigned int translateTo = pSubUnit->getTranslateToIndex();
//                    translateTo = getUseableNeighborForCOM(pModel, translateTo);
                    translateTo = getTranslateToPointNearCOM(pModel, translateTo);
                    pSubUnit->translateToTransformedPosition(translateTo, pModel);
                }
                whole = true;
//                single=true;
            }

            // calculate new energy
            calculateModelPrDistributionSymX(&bead_indices, &subunits, &binCount, workingLimit, temp_violations, pModel, pData );
            testKL = pData->getScore(binCount);

            tempRgPotential = moments(&subunits, &bead_indices, workingLimit, pData->getRg(), pModel);
            tempCOMEnergy = calculateCOMPotential(&cofm, subunits, cePotentials);
            tempIntraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);

            this_energy = testKL +
                          beta*temp_violations +
                          gammaConstant*tempIntraContacts +
                          zetaConstant*tempRgPotential +
                          epsilonConstant*tempCOMEnergy;

            if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                isUpdated=true;
                std::sprintf(addRemoveText, "TRANSLATED");
                if (whole){
                    ++wholeCount;
                }
            } else { // undo changes (rejecting)
                if (single){
                    pSubUnit->undoTranslateTo(pModel);
                } else {
                    for(unsigned int i=0; i<subunits.size(); i++){
                        pSubUnit = &(subunits[i]);
                        pSubUnit->undoTranslateTo(pModel);
                    }
                }

                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::sprintf(addRemoveText, " NOT TRANSD ");
            }
        }

        std::cout << "*******************                                        *******************" << std::endl;
        std::printf("      ROUND : %5d         INDEX => %d ( %d )\n", round, index, wholeCount);
        std::printf("       TEMP : %-.4E MAXSTEPS => %.0f ( %i ) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        std::printf("     ACCEPT : %.5f  FAILURES => %i %s \n", acceptRate, failures, addRemoveText);
        std::printf(" VIOLATIONS : %.3E COM : %-4.3E RG : %.3E \n", totalViolations, comEnergy, rgPotential);
        std::printf("   INTERCON : %.3E  \n", intraContacts);
        std::printf("       D_KL => %-5.4E ENRGY : %.4E ( %.4E ) \n", currentKL, current_energy, lowest_energy);

        // Adaptive simulated annealing part
        if (isUpdated){

            current_energy = this_energy;
            currentKL = testKL;
            totalViolations = temp_violations;
            rgPotential = tempRgPotential;
            comEnergy = tempCOMEnergy;
            intraContacts = tempIntraContacts;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);

        if ( current_energy < lowest_energy){ // if counter too small, add/remove may not sample sufficiently
            // if outside window, renormalize
            // if less than, reset, the larger eta is the more it dominates the minimization
//            if ( (zetaConstant*rgPotential)/currentKL < zeta) {
//                current_energy -= zetaConstant*rgPotential;
//                //zetaConstant = currentKL*zeta/((1.0 - zeta)*rgPotential);
//                zetaConstant = (rgPotential == 0) ? 0.0001 : zeta*currentKL/((1.0-zeta)*rgPotential); // rg potential
//                current_energy += zetaConstant*rgPotential;
//            }
            lowest_energy = current_energy;
        }
    }

    if (printIt){
        pModel->writeSymXModelToFile(currentKL,
                                     workingLimit,
                                     bead_indices,
                                     subunits,
                                     binCount,
                                     name + "_" + std::to_string(round),
                                     this,
                                     pData,
                                     1000,
                                     current_volume,
                                     0);
    }

//    return (currentKL + beta*totalViolations + gammaConstant*intraContacts);
    return (float)(currentKL + beta*totalViolations + zetaConstant*rgPotential + gammaConstant*intraContacts);
}


/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSymX(std::vector<unsigned int> *subUnit_indices, std::vector<SubUnit> *subUnits, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, double &potential, Model *pModel, Data *pData) {

    std::fill(binCount->begin(), binCount->end(), 0);
    potential=0.0d;
    const unsigned int totalSubUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalNonParentSubUnits = pModel->getNumberOfSubUnits()-1;
    float distance_to;
    float diff;//, potential=0;

    // float max=0;
    // all subUnits are identical
    for (unsigned int i=0; i < indicesWorkingLimit; i++){

        const vector3 * const tempVec2 = &(pModel->getBead((*subUnit_indices)[i])->getVec());

        for(unsigned int  j=(i+1); j < indicesWorkingLimit; j++){
            distance_to =((*tempVec2) - pModel->getBead((*subUnit_indices)[j])->getVec()).length();
            (*binCount)[pData->convertToBinUnsignedInt(distance_to)] += totalSubUnits; // some distances will exceed dmax
        }
    }



    for (unsigned int s=0; s < totalNonParentSubUnits; s++){  //

        const vector3 * pSubUnit = &(*subUnits)[s].getStartToTransformedCoordinates(); //

        for(unsigned int  j=0; j < indicesWorkingLimit; j++){ // for each subUnit, go through its coordinates
            const vector3 * const pVec1 = pSubUnit + j; // constant pointer to a constant vector

            for (unsigned int ss=s+1; ss < totalSubUnits; ss++){  // create sym related subunits and add to coordinates vector

                const vector3 * pSubUnit2 = &(*subUnits)[ss].getStartToTransformedCoordinates(); // initialized as emptyin Model class

                for(unsigned int  jj=0; jj < indicesWorkingLimit; jj++){
                    distance_to =(*pVec1 - *(pSubUnit2+jj)).length();

                    ++(*binCount)[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax
                    // binCount is a vector of unsigned ints and access
                    // convert to Bin is return an unsigned integer or short
                    // seems like some incompatibility with short and regular int for finding bin locations
                    if (distance_to < violation_limit){
                        diff = (distance_to - contactCutOff);
                        potential += diff*diff*diff*diff*diff*diff;
                        // violation++;
                    }

                }
            }
        }
    }

}



double Anneal::intraSubUnitContact(std::vector<SubUnit> * subUnits, std::vector<unsigned int> * cvxIndices, Model *pModel){
    //calculate contacts per bead for selectedIndex
    const unsigned int totalSubUnits = pModel->getNumberOfSubUnits();
    float contactCutOffLimit = 1.721f*contactCutOff; // contactCutOff is 2x bead_radius
    // float midpoint = (contactCutOffLimit - violation_limit)*0.5 + violation_limit;
    // only look at CVX hull points
    // create first subunit from selected indices and populate coordinates
    // determine contacts
    // For each subUnit, check that it is in contact with at least another subunit
    unsigned int totalCVXIndices = cvxIndices->size();

    float distance_to, closestPoint;// = pModel->getDiameterOfUniverse();
    // check first subunit, must be in contact with second
    float potential = 0;

    for (unsigned int s=0; s< (totalSubUnits-1); s++){

        SubUnit & pSubUnit1 = (*subUnits)[s];
        SubUnit & pSubUnit2 = (*subUnits)[s+1];

        closestPoint = pModel->getDiameterOfUniverse();

        for(unsigned int subIndex1=0; subIndex1 < totalCVXIndices; subIndex1++){

            const vector3 * const pVec1 = &pSubUnit1.getPointerToTransformedCoordinate((*cvxIndices)[subIndex1]);

            for(unsigned int subIndex2=0; subIndex2 < totalCVXIndices; subIndex2++){

                distance_to = (*pVec1 - (pSubUnit2.getPointerToTransformedCoordinate((*cvxIndices)[subIndex2]))).length();

                if (distance_to < contactCutOffLimit && distance_to > violation_limit){ /*!< distance between lattice points to be counted as a neighbor */
                    // if no distance can be considered a contact between neighbors must consider not touching
                    //closestPoint = contactCutOff;
                    closestPoint = contactCutOffLimit;
                    goto loop2;
                } else if (distance_to<closestPoint){
                    closestPoint = distance_to;
                }
            }
        }

        loop2:
        float diff = closestPoint - contactCutOffLimit;
        potential += diff*diff;
    }

    auto totalContraints = (unsigned int)intraSubUnitContacts.size();
    if (intraSubUnitContacts.size() > 0){

        for(unsigned int m=0; m<totalContraints; m++){
            SubUnit & pSubUnit1 = (*subUnits)[intraSubUnitContacts[m].getSubUnit1()];
            SubUnit & pSubUnit2 = (*subUnits)[intraSubUnitContacts[m].getSubUnit2()];

            closestPoint = pModel->getDiameterOfUniverse();

            for(unsigned int subIndex1=0; subIndex1 < totalCVXIndices; subIndex1++){

                const vector3 * const pVec1 = &pSubUnit1.getPointerToTransformedCoordinate((*cvxIndices)[subIndex1]);

                for(unsigned int subIndex2=0; subIndex2 < totalCVXIndices; subIndex2++){

                    distance_to = (*pVec1 - (pSubUnit2.getPointerToTransformedCoordinate((*cvxIndices)[subIndex2]))).length();

                    if (distance_to < contactCutOffLimit && distance_to > violation_limit){ /*!< distance between lattice points to be counted as a neighbor */
                        // if no distance can be considered a contact between neighbors must consider not touching
                        closestPoint = contactCutOffLimit;
                        goto loop3;
                    } else if (distance_to<closestPoint){
                        closestPoint = distance_to;
                    }
                }
            }

            loop3:
            float diff = closestPoint - contactCutOffLimit;
            potential += diff*diff;
        }
    }

    return potential/(double)(totalSubUnits + totalContraints - 1) + 0.00001;
}


/**
 * beadsInUse does not have to be sorted
 *
 * prBins is the model P(r) distribution
 * Index to remove should still be in beadsInuse for all subUnits
 */
inline double Anneal::removeFromPrSymX(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, std::vector<SubUnit> & subUnits, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    //
    // remove interdomain distances
    //
    unsigned int violations=0;
    double diff, potential=0.0d;
    float distance_to;

    const unsigned int totalSubUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = totalSubUnits*workingLimit;
    std::vector<const vector3 *> pCoordinates(totalCoordinates);

    // vector3 * ptr = &pCoordinates.front();
    // compile coordinates
    unsigned int count = 0;
    for (unsigned int s=0; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector

        const vector3 * pSubUnit = &subUnits[s].getStartToTransformedCoordinates(); // return address to start

        for(unsigned int  j=0; j < workingLimit; j++){ //new position is not in same index location as parent subUnit
            pCoordinates[count] = (pSubUnit+j);
            ++count;
        }
    }

    // get location of bead to remove in subUnits
    const unsigned  int indexInSubUnit = subUnits[0].getIndexOfBeadtoRemove(removeMeSubUnitIndex);
    // now go through each subunit, adjust Pr
    for (unsigned int s=0; s < totalSubUnits; s++){  //

        unsigned int basis = s*workingLimit+indexInSubUnit; // location of position to remove

        const vector3 * prefVec = pCoordinates[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - (*pCoordinates[next_i])).length();
            --prBins[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                diff = (distance_to - contactCutOff);
                potential += diff*diff*diff*diff*diff*diff;
                ++violations;
            }
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - (*pCoordinates[next_i])).length();
            --prBins[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                diff = (distance_to - contactCutOff);
                potential += diff*diff*diff*diff*diff*diff;
                ++violations;
            }
        }
    }


    /**
     *
     * we double remove the symmetry related vectors, need to add back
     */
    // remove double counts
    // start with parent subUnit
    for (unsigned int s=0; s < totalSubUnits; s++){

        const vector3 * prefVec1 = pCoordinates[s*workingLimit+indexInSubUnit];

        for (unsigned int ss=s+1; ss < totalSubUnits; ss++){
            distance_to = ((*prefVec1) - (*pCoordinates[ss*workingLimit+indexInSubUnit])).length();
            ++prBins[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                diff = (distance_to - contactCutOff);
                potential -= diff*diff*diff*diff*diff*diff;
                --violations;
            }
        }
    }

    return potential;
}


void Anneal::swapInSubUnits(unsigned int indexOfCurrent, unsigned int indexOfNew, std::vector<SubUnit> & subUnits, Model * pModel){

    unsigned int total = subUnits.size();
    for(unsigned int i=0; i<total; i++){
        subUnits[i].swapIndices(indexOfCurrent, indexOfNew, pModel);
    }
}

void Anneal::removeFromSubUnits(unsigned int indexOfPositionToRemove, std::vector<SubUnit> & subUnits){
    unsigned int total = subUnits.size();
    for(unsigned int i=0; i<total; i++){
        subUnits[i].removeIndex(indexOfPositionToRemove);
    }
}

/**
 * The new position to Add MUST BE WITHIN THE WORKING LIMIT
 * beadsInUse must be sorted up to workingLimit
 *
 * for subUnits, newly added bead is in the last position
 *
 * @param addMeSubUnitIndex
 * @param beadsInUse
 * @param workingLimit
 * @param prBins
 * @param pModel
 * @param pData
 * @return
 */
double Anneal::addToPrSymX(std::vector<unsigned int> & beadsInUse, std::vector<SubUnit> *subUnits, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){
    // each bead has a sym related partner in pModel
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin

    unsigned int violations=0;
    double potential=0, diff;
    float distance_to;

    const unsigned int totalSubUnits = pModel->getNumberOfSubUnits(); // specific to subUnits vector
    const unsigned int totalCoordinates = totalSubUnits*workingLimit;
    std::vector<const vector3 *> pCoordinates(totalCoordinates);
    /*
     * compile coordinates of the complex
     */
    unsigned int count = 0;
    for (unsigned int s=0; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector

        const vector3 * const pSubUnit = &(*subUnits)[s].getStartToTransformedCoordinates(); // initialized as emptyin Model class

        for(unsigned int j=0; j < workingLimit; j++){ //new position is not in same index location as parent subUnit
            pCoordinates[count] = (pSubUnit+j);
            ++count;
        }
    }


    // now go through each subunit, adjust Pr
    for (unsigned int s=0; s < totalSubUnits; s++){  //

        unsigned int basis = (s+1)*workingLimit-1; // location of newly added position

        const vector3 * const prefVec = pCoordinates[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - (*pCoordinates[next_i])).length();
            ++prBins[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                diff = (distance_to - contactCutOff);
                potential += diff*diff*diff*diff*diff*diff;
                ++violations;
            }
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - (*pCoordinates[next_i])).length();
            ++prBins[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                diff = (distance_to - contactCutOff);
                potential += diff*diff*diff*diff*diff*diff;
                ++violations;
            }
        }
    }


    // remove double counts
    // start with parent subUnit
    for (unsigned int s=0; s < totalSubUnits; s++){

        const vector3 * const prefVec1 = pCoordinates[(s+1)*workingLimit-1]; // location of the newly added position

        for (unsigned int ss=s+1; ss < totalSubUnits; ss++){
            distance_to = ((*prefVec1) - (*pCoordinates[(ss+1)*workingLimit-1])).length();
            --prBins[pData->convertToBinUnsignedInt(distance_to)]; // some distances will exceed dmax
            if (distance_to < violation_limit){
                diff = (distance_to - contactCutOff);
                potential -= diff*diff*diff*diff*diff*diff;
                --violations;
            }
        }
    }

    return potential;
}

/**
 *
 *
 */
double Anneal::removeLatticePositionToModelSymX(
        std::vector<unsigned int> & bead_indices,
        std::vector<SubUnit> & subUnits,
        std::vector<unsigned int> & modelPrBins,
        unsigned int * pWorkingLimit,
        const unsigned int * pLatticePointToRemove, Model * pModel, Data *pData){

    auto pBeginIt = bead_indices.begin();
    auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    double violations = removeFromPrSymX(*pLatticePointToRemove, bead_indices, subUnits, *pWorkingLimit, modelPrBins, pModel, pData);

    // reduce the workingLimit
    // if wl = 10
    // 0 1 2 3 4 5 6 7 8 9 10
    // remove 4
    // 0 1 2 3 9 5 6 7 8 4 10
    //
    // 0 1 2 3 5 6 7 8 9 4 10
    //
    *pWorkingLimit -= 1;
    // swap selected point to the workingLimit
    std::iter_swap(itIndex, pBeginIt + *pWorkingLimit);
    // still need to sort, swap changes the order
    // to save some time, sort should only be from swapped point to workingLimit
    std::sort(pBeginIt, pBeginIt + *pWorkingLimit);
    return violations;
}


void Anneal::undoRemoveFromSubUnits(std::vector<SubUnit> &subUnits) {
    unsigned int total = subUnits.size();
    for(unsigned int i=0; i<total; i++){
        subUnits[i].undoRemove();
    }
}

void Anneal::addToSubUnits(unsigned int indexOfPositionToAdd, std::vector<SubUnit> & subUnits, Model * pModel){
    unsigned int total = subUnits.size();
    for(unsigned int i=0; i<total; i++){
        subUnits[i].addIndex(indexOfPositionToAdd, pModel);
    }
}

void Anneal::undoAddToSubUnits(std::vector<SubUnit> &subUnits) {
    unsigned int total = subUnits.size();
    for(unsigned int i=0; i<total; i++){
        subUnits[i].undoAdd();
    }
}



/**
 * apply Parallel Axis theorem to calculate Rg based on CVX points
 * @param subUnits
 * @param cvxIndices
 * @param pModel
 * @return
 */
double Anneal::moments(std::vector<SubUnit> *subUnits, std::vector<unsigned int> *indices, unsigned int totalIndices, float targetRg, Model *pModel) {

    const int totalSubUnits = pModel->getNumberOfSubUnits();
    // only look at CVX hull points
    // create first subunit from selected indices and populate coordinates
    // determine contacts
    // For each subUnit, check that it is in contact with at least another subunit
    const float inv = 1.0f/(float)totalIndices;
    float sum_x=0, sum_y=0, sum_z=0;
    for (unsigned int i=0; i<totalIndices; i++){
        const vector3 * const pVec = &pModel->getBead((*indices)[i])->getVec();
        sum_x += pVec->x;
        sum_y += pVec->y;
        sum_z += pVec->z;
    }
    // center coordinates
    const vector3 centerOfMass = vector3(sum_x*inv, sum_y*inv, sum_z*inv);

    // calculate Rg of SubUnit
    float squared = 0;
    for (unsigned int i=0; i<totalIndices; i++){
        const vector3 pVec1 = (pModel->getBead((*indices)[i])->getVec() - centerOfMass);
        squared += pVec1.dot(pVec1); // squared distance
    }

    const float massFraction = 1.0f/(float)pModel->getNumberOfSubUnits();
    float squaredRg = totalSubUnits*(squared*inv);
    const float inv2 = massFraction*massFraction;

    // distance between subUnits
    for (unsigned int s=0; s< (totalSubUnits-1); s++){

        SubUnit & pSubUnit1 = (*subUnits)[s];

        const vector3 pVec1 = pSubUnit1.translateAndRotatePoint(&centerOfMass);

        for (unsigned int ss=s+1; ss< totalSubUnits; ss++){

            SubUnit & pSubUnit2 = (*subUnits)[ss];
            const vector3 pVec2 = pSubUnit2.translateAndRotatePoint(&centerOfMass);
            float dis = (pVec1-pVec2).length();
            squaredRg += inv2*dis*dis;
        }
    }

    //direct Rg calculation
    unsigned int totalCoordinates = totalSubUnits*totalIndices;
    std::vector<const vector3 *> pCoordinates(totalCoordinates);
    vector3 com = vector3(0,0,0);
    unsigned int count = 0;
    for (unsigned int s=0; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector
        SubUnit & pSubUnit1 = (*subUnits)[s];

        const vector3 * pSubUnit = &(pSubUnit1.getStartToTransformedCoordinates()); // return address to start

        for(unsigned int  j=0; j < totalIndices; j++){ //new position is not in same index location as parent subUnit
            pCoordinates[count] = (pSubUnit+j);
            com += *(pSubUnit+j);
            ++count;
        }
    }

    const float invSum = 1.0f/(float)count;
    com.x = com.x*invSum;
    com.y = com.y*invSum;
    com.z = com.z*invSum;

    float actualSecondMoment=0;
    for(unsigned int i=0; i<totalCoordinates; i++){
        float diff = (*pCoordinates[i] - com).length();
        actualSecondMoment += diff*diff;
    }

    actualSecondMoment *= 1.0/(float)totalCoordinates;
    float expRg = sqrt(actualSecondMoment);
    float diff = (targetRg-expRg)/targetRg;

    //std::cout << " Rg " << expRg << " " << targetRg << " DIFF " << diff << std::endl;
    return (diff*diff);
}


void Anneal::logToCeFile(char * text){
    FILE *pFile;
    pFile = fopen("ce_probabilities.txt", "a");
    fprintf(pFile, text);
    fclose(pFile);
}


float Anneal::calculateCOMPotential(const vector3 * const com, std::vector<SubUnit> & subunits, std::vector<CEPotential> & cePotentials){
    float energy = 0;
    float distance;

    unsigned int totalInP = cePotentials.size();
    SubUnit * pSubUnit1, * pSubUnit2;
    CEPotential * pCEPotential;

    for(unsigned int i=0; i<totalInP; i++){
        pCEPotential = &cePotentials[i];

        pSubUnit1 = &subunits[pCEPotential->getSubUnit1()];
        pSubUnit2 = &subunits[pCEPotential->getSubUnit2()];

        distance = (pSubUnit1->translateAndRotatePoint(com) - pSubUnit2->translateAndRotatePoint(com)).length();
        energy += pCEPotential->calculatePotentialEnergy(distance);
    }

    return energy;
}



/**
 * based on COM, pick a random point that is a neighbor of COM
 * returned point may or may not be in part of the reconstructed model (that is, already in use)
 *
 */
inline unsigned int Anneal::getUseableNeighborForCOM(Model *pModel,
                                            unsigned int & com){

    auto it = pModel->getPointerToNeighborhood(com);
    unsigned int noNeighborLimit = pModel->getNeighborLimit();
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    unsigned int count=0;
    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ( neighbor < noNeighborLimit ){ // if end of set, means not in use
            possibleNeighbors[count] = neighbor;
            count++;
        } else if (neighbor == noNeighborLimit) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return noNeighborLimit;
    } else if (count == 1){
        return possibleNeighbors[0];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);

    return possibleNeighbors[randomIndex(gen)];
}


/**
 * perform random walk starting from COM
 *
 * randomly select one
 */
inline unsigned int Anneal::getTranslateToPointNearCOM(Model *pModel, unsigned int & com){
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    unsigned int neighborToMoveTo = getUseableNeighborForCOM(pModel, com);
    for(int i=0; i<5; i++){
        neighborToMoveTo = getUseableNeighborForCOM(pModel, neighborToMoveTo);
        if (distribution(gen) < 0.3335){
            break;
        }
    }

    return neighborToMoveTo;
}



/*
 * does not include COMPotential
 */
float Anneal::estimateInterSubUnitEnergyX(std::vector<SubUnit> & subunits, Model *pModel, Data *pData) {

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    const unsigned int totalSubUnits = pModel->getNumberOfSubUnits();
    std::uniform_int_distribution<unsigned int> randomSubUnitIndex (0, totalSubUnits-1);
    double violations;

    // build potentials for each intersubunit distance
    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    std::vector<unsigned int> cvxIndices;
    std::vector<unsigned int> tempCVXIndices;

    pModel->copyStartingModelIntoVector(bead_indices);
    std::fill(binCount.begin(), binCount.end(), 0);

    updateCVXIndices(workingLimit, bead_indices, cvxIndices, pModel);
    double intraContacts = intraSubUnitContact(&subunits, &cvxIndices, pModel);

    const float gammaConstant = 0.001; // contacts between subUnit contacts

    // initialize for new round
    // fill binCount for first time
    calculateModelPrDistributionSymX(&bead_indices, &subunits, &binCount, workingLimit, violations, pModel, pData );
    float currentKL = pData->getScore(binCount);

    double rgPotential = moments(&subunits, &bead_indices, workingLimit, pData->getRg(), pModel);
    const float zetaConstant = 0.0001;

    float current_energy = currentKL +
                           gammaConstant*intraContacts +
                           zetaConstant*rgPotential;

    beta=0.001;
    current_energy += beta*violations;

    //(float)(currentKL + beta*totalViolations + zetaConstant*rgPotential + gammaConstant*intraContacts);

    return current_energy;
}



float Anneal::totalInterSubUnitDistances(const vector3 * const comvector, std::vector<SubUnit> & subunits, Model *pModel){

    unsigned int totalInterSubUnitDistancesInPotential = subunits.size();
    float total = 0;
    for(unsigned int i=0; i < totalInterSubUnitDistancesInPotential; i++){

        vector3 pVec1 = subunits[i].translateAndRotatePoint(comvector);

        for(unsigned int j=i+1; j < totalInterSubUnitDistancesInPotential; j++){

            total += (pVec1 - subunits[j].translateAndRotatePoint(comvector)).length();
        }
    }

    return total;
}