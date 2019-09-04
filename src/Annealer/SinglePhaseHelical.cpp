//
// Created by xos81802 on 09/10/2018.
//

#include <random>
#include "../Anneal.h"


bool Anneal::createInitialModelHelicalSymmetry(Model *pModel, Data *pData) {

    //highTempRounds*=1.62;

    std::cout << " --------------------------------------------------- " << std::endl;
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    std::bernoulli_distribution binaryDistribution(0.50);

    contactCutOff = interconnectivityCutOff;
    violation_limit = 1.7*pModel->getBeadRadius();

    float priorTheta, newTheta;

    float delta_theta = 0.05f*pModel->getRotationPerSubUnit();
    float priorRise, newRise;
    float delta_rise = pModel->getRisePerSubUnit()*0.025f;
    const float lowerRiseThreshold = pModel->getRisePerSubUnit()*0.80f;
    const float upperRiseThreshold = pModel->getRisePerSubUnit()*1.2f;
    const float lowerThetaThreshold = pModel->getRotationPerSubUnit()*0.98f;
    const float upperThetaThreshold = pModel->getRotationPerSubUnit()*1.02f;
    // convert distances to ShannonBin membership
    maxbin = pModel->getMaxBin(pData);

    float hyp = 0.5f*(pData->getDmax() + 2*pModel->getZaxis());
    float longest = 2*std::sqrt(hyp*hyp + pModel->getXaxis()*pModel->getXaxis());
    auto maxbintemp = pData->convertToBin(longest);

    maxbin = maxbintemp;  // maximum bin index, so vector size must be +1
    std::cout << "         INITIAL MAX BIN : " << maxbin << std::endl;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "          MAX BIN SEARCH : " << maxbintemp << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's
    //std::clock_t start;
    // c-style is slightly faster for large vector sizes
    // start = std::clock();
    unsigned int *ptr = (totalBeadsInSphere != 0) ? &bead_indices.front() : nullptr;
    for (unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    // cout << "C-STYLE " << (std::clock() - start)/(double) CLOCKS_PER_SEC << endl;
    std::random_device rd;
    std::mt19937 gen(rd());

    float invBeadVolume = 1.0f / pModel->getBeadVolume();
    lowerV = (unsigned int) (pModel->getHelicalVolume() - pModel->getHelicalVolume() * 0.29);
    upperV = (unsigned int) (pModel->getHelicalVolume() + pModel->getHelicalVolume() * 0.05);
    const float targetSubUnitVolume = (float) upperV / (float) pModel->getNumberOfSubUnits();

    auto lowerN = (unsigned int) std::round(lowerV * invBeadVolume);
    auto upperN = (unsigned int) std::round(upperV * invBeadVolume);

    // pick random number between lowerN and upperN
    std::uniform_int_distribution<unsigned int> number_of_beads_to_use(lowerN, upperN);
    unsigned int workingLimit = number_of_beads_to_use(gen);
    unsigned int lowestWorkingLimit = workingLimit;
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    // create initial model
    // randomize bead indices
    // sort to workingLimit
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    // make copy as initial best model
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

    std::cout << "        WORKINGLIMIT SET : " << workingLimit << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";
    coordT points[3 * (upperN + (unsigned int) (0.50 * upperN))];

    float subunit_test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float lowest_volume = current_volume;
    float initialVolume = current_volume;

    std::cout << "              CVX VOLUME : " << current_volume << std::endl;
    EulerTour eulerTour(bead_indices.cbegin(), workingLimit, pModel);

    unsigned int lowestNumberOfComponents = 0, tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();
    std::cout << "      CREATED EULER TOUR : " << currentNumberOfComponents << std::endl;

    pModel->setNumberOfSubUnitsHelical((unsigned int) std::ceil(pData->getDmax() / pModel->getRisePerSubUnit()));

    unsigned int violations, temp_violations;
    calculateModelPrDistributionSymHelical(&bead_indices, &binCount, workingLimit, violations, pModel, pData);

    // fill binCount for first time
    float testKL, lowestKL, currentKL = pData->getScore(binCount);

    float hlambda = pData->getIsIntensity() ? 10.0f : 0.001f;
    float muConstant = pData->getIsPr() ? 0.00001f : 0.1f; // chi tends to start in 100's whereas DKL is 0.1;

    const float invTarVolSubUnitmuConstant = muConstant / targetSubUnitVolume;

    unsigned int failures = 0;
    bool updateCVX = false;
    //float muConstant = 0.0001;//mu*0.00001f/((1.0f-mu)*upperV/(float)pModel->getNumberOfSubUnits());
    // high temp search
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");
    float test_energy, current_energy = currentKL + hlambda * (currentNumberOfComponents - 1) * (currentNumberOfComponents - 1);

    std::cout << "       INITIAL => ENERGY : " << current_energy << std::endl;
    std::cout << "       INITIAL =>   D_KL : " << currentKL << std::endl;
    std::cout << "       INITIAL => VOLUME : " << current_volume << std::endl;
    std::cout << "              VIOLATIONS : " << violations << std::endl;
    std::cout << "       INITIAL        WL : " << workingLimit << std::endl;


    char addRemoveText[50];
    float sum_x_squared = (workingLimit*workingLimit);
    unsigned int volumeCount = 1, testIndex;
    float volumeSum = current_volume, workingLimitSum = workingLimit;

    unsigned int totalViolations = violations; //(double)subUnitWorkingLimit;
    // what is tolerable number of violations per subsunit?
    current_energy += beta * totalViolations;

    float lowest_energy = current_energy;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup cop

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int high;
    auto lowTempStop =  (double)highT;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    float temp_subunit_volume_energy, current_subunit_volume_energy =
            std::abs(current_volume - targetSubUnitVolume) * invTarVolSubUnitmuConstant;

    current_energy += current_subunit_volume_energy;

//    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "initial_subunit_annealed", this, pData, 0, current_volume, 0);
//    pModel->writeHelicalSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "initial_hsym_annealed", this, pData, 0, current_volume, 0);


    for (high = 0; high < highTempRounds; high++) { // iterations during the high temp search

        if (distribution(gen) > 0.5){

            if (distribution(gen) > 0.17) {
                // create points from workSet to determine HULL
                std::sprintf(addRemoveText, "POSITIONAL");
                unsigned int swap1 = bead_indices[randomIndex(gen)];

                // setup parameters for hull
                if ( (distribution(gen) < 0.353) || currentNumberOfComponents > 1){ // only move CVX points if true

                    for (unsigned int i=0; i < workingLimit; i++) {
                        beadToPoint(&points[i*3], pModel->getBead(bead_indices[i]));
                        active_indices[i] = bead_indices[i];
                    }

                    // needs to be optimized
                    qh_new_qhull(3, workingLimit, points, 0, flags, nullptr, nullptr);
                    vertexT * vertices = qh vertex_list;
                    auto totalV = (unsigned int)(qh num_vertices);

                    // only move CVX hull points
                    std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
                    for (unsigned int v=0; v < totalV; v++) { //
                        indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                        vertices = vertices->next;
                    }

                    qh_freeqhull(true);

                    std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
                    unsigned int randomInt = randomVertices(gen);
                    swap1 = indices_to_check[randomInt];
                } // recalculate

                auto itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                // remove selected index from P(r)
                removeFromPrSymHelical(swap1, bead_indices, workingLimit, binCount, pModel, pData);

                beads_in_use_tree.erase(swap1);
                eulerTour.removeNode(swap1);

                unsigned int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while (neighbor == noNeigborIndex || neighbor == swap1){
                    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }
                // make the swap, sort and update P(r)
                auto pSwap2 = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), neighbor);
                std::iter_swap(itIndex, pSwap2);
                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                temp_violations = addToPrSymHelical(neighbor, bead_indices, workingLimit, binCount, pModel, pData);

                testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry

                tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel); // calculate volume of the subUnit
                temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

                if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.insert(neighbor);
                    updateCVX = true;
                } else {
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    beads_in_use_tree.insert(swap1);
                    eulerTour.removeNode(neighbor);
                    eulerTour.addNode(swap1, pModel);
                }

            } else {

                if ( ((distribution(gen) < 0.5)) || workingLimit > upperN) { // REMOVE beads from sorted list into useable range < deadLimit
                    // randomly swap positions with end of workingLength, could remove CVX Hull Point
                    std::sprintf(addRemoveText, "  REMOVE  ");

                    testIndex = bead_indices[randomIndex(gen)];

                    removeLatticePositionToModelSymHelical(bead_indices, binCount, &workingLimit, &testIndex, pModel, pData);
                    temp_violations = getViolationsHelical(&bead_indices, workingLimit, pModel, pData);

                    testKL = pData->getScore(binCount);

                    subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel); // calculate volume of the subUnit
                    temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

                    if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.erase(testIndex);
                        updateCVX=true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        auto beginIt = bead_indices.begin();
                        auto beginBinCount = binCount.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt,
                                                              &workingLimit,
                                                              &binCountBackUp,
                                                              &beginBinCount);
                    }

                } else { // ADD beads
                    std::sprintf(addRemoveText, "    ADD   ");

                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    while (testIndex == noNeigborIndex) { // find a new position
                        testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    }

                    auto itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);
                    addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);

                    temp_violations = addToPrSymHelical(testIndex, bead_indices, workingLimit, binCount, pModel, pData);

                    testKL = pData->getScore(binCount);

                    subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel); // calculate volume of the subUnit
                    temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                    tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

                    if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.insert(testIndex);
                        updateCVX=true;
                    } else { // undo changes (rejecting)
                        auto beginBinCount = binCount.begin();
                        restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                        eulerTour.removeNode(testIndex);
                    }
                }
            }

        } else {

            std::sprintf(addRemoveText, "  ADJUST  ");
            priorRise = pModel->getRisePerSubUnit();
            priorTheta = pModel->getRotationPerSubUnit();

            //adjust delta and theta
            float plusMinus =  (binaryDistribution(gen) > 0) ? -1.0f : 1.0f;

            newRise = (binaryDistribution(gen) > 0) ? priorRise + plusMinus*delta_rise : priorRise;
            plusMinus =  (binaryDistribution(gen) > 0) ? -1.0f : 1.0f;
            newTheta = (binaryDistribution(gen) > 0) ? priorTheta + plusMinus*delta_theta : priorTheta;
            // determine
            if (newRise < lowerRiseThreshold){
                newRise = lowerRiseThreshold;
            } else if (newRise > upperRiseThreshold){
                newRise = upperRiseThreshold;
            }

            if (newTheta < lowerThetaThreshold){
                newTheta = lowerThetaThreshold;
            } else if (newTheta > upperThetaThreshold){
                newTheta = upperThetaThreshold;
            }

            pModel->setRiseAndRotationPerSubUnit(newRise, newTheta);
            pModel->setNumberOfSubUnitsHelical((unsigned int) std::ceil(pData->getDmax() / pModel->getRisePerSubUnit()));

            subunit_test_volume = current_volume;
            temp_subunit_volume_energy = current_subunit_volume_energy;
            tempNumberOfComponents = currentNumberOfComponents;

            calculateModelPrDistributionSymHelical(&bead_indices, &binCount, workingLimit, temp_violations, pModel, pData);
            testKL = pData->getScore(binCount);

            test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

            if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                pModel->setRiseAndRotationPerSubUnit(priorRise, priorTheta);
                pModel->setNumberOfSubUnitsHelical((unsigned int) std::ceil(pData->getDmax() / pModel->getRisePerSubUnit()));
            }

        }




        if (updateCVX){

            current_energy = test_energy;
            currentKL = testKL;
            current_volume = subunit_test_volume;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            current_subunit_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASAConstantTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts
        std::printf("*******************             %s                 ******************* \n", addRemoveText);
        std::printf("         MAXSTEPS => %i (%4i) ACCEPTRATE : %.3f TEMP : %.3E\n", highTempRounds, high, acceptRate, lowTempStop);
        std::printf("            GRAPH => %3i\n", currentNumberOfComponents);
        std::printf("            UPPER => %i LOWER => %i  LIMIT : %5i \n", upperN, lowerN, workingLimit);
        std::printf("           VOLUME => %.0f (%.0f)  MU => %.4E  MU*VOL => %.3E\n", current_volume, targetSubUnitVolume, muConstant, current_subunit_volume_energy);
        std::printf("            HELIX : RISE => %.1f THETA %.1f TOTALSUBS %i \n", priorRise, priorTheta, pModel->getTotalNumberOfSubUnitsInSingleFiber());
        std::printf("       VIOLATIONS => %d  BETA => %.4E (%i)  \n", totalViolations, beta, volumeCount);
        std::printf("              D_KL : %.4E ENRGY: %.4E (%.4E) \n", currentKL, current_energy, lowest_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy ){
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            lowestWorkingLimit = workingLimit;
            lowest_energy = current_energy;
            lowestKL = currentKL;
            lowestNumberOfComponents = currentNumberOfComponents;
            workingLimitSum += workingLimit;
            sum_x_squared += workingLimit*workingLimit;
            volumeSum += current_volume;
            volumeCount++;

        }

//        calculateModelPrDistributionSymHelical(&bead_indices, &testBinCount, workingLimit, violations, pModel, pData);
//        // fill binCount for first time
//        float testKL = pData->getScore(testBinCount);
//
//        if (testKL != currentKL){
//            exit(0);
//        }

    }


    highTempStartForCooling = (float)lowTempStop;
    // pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage = workingLimitSum/(float)volumeCount;
    float volumeStdev = std::sqrt(sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage);

    // remove points close to hull

    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);


    std::string nameOfModel = pModel->writeModelToFile2(lowestKL, workingLimit, bead_indices, binCount, "subunit_hiT_", this, pData, 0, current_volume, 0);
    pModel->writeHelicalSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "hsym_hiT_", this, pData, 0, current_volume, 0);

    if (lowestNumberOfComponents != 1){
        return false;
    }
    return true;
}


/**
 *
 */
std::string Anneal::refineSymModelHelical(Model *pModel, Data *pData, std::string nameTo){

    std::cout << " --------------------------------------------------- " << std::endl;
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    std::bernoulli_distribution binaryDistribution(0.50);
    contactCutOff = interconnectivityCutOff;
    violation_limit = 1.71f*pModel->getBeadRadius();

    float hyp = 0.5f*(pData->getDmax() + 2*pModel->getZaxis());
    float longest = 2*std::sqrt(hyp*hyp + pModel->getXaxis()*pModel->getXaxis());
    auto maxbintemp = pData->convertToBin(longest);

    float priorTheta, newTheta;
    float delta_theta = 0.05f*pModel->getRotationPerSubUnit();
    float priorRise, newRise;
    float delta_rise = pModel->getRisePerSubUnit()*0.025f;
    const float lowerRiseThreshold = pModel->getRisePerSubUnit()*0.80f;
    const float upperRiseThreshold = pModel->getRisePerSubUnit()*1.2f;
    const float lowerThetaThreshold = pModel->getRotationPerSubUnit()*0.9f;
    const float upperThetaThreshold = pModel->getRotationPerSubUnit()*1.1f;

    maxbin = maxbintemp;  // maximum bin index, so vector size must be +1
    std::cout << "         INITIAL MAX BIN : " << maxbin << std::endl;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's

    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // cout << "C-STYLE " << (std::clock() - start)/(double) CLOCKS_PER_SEC << endl;
    std::random_device rd;
    std::mt19937 gen(rd());

    float invBeadVolume = 1.0f / pModel->getBeadVolume();
    lowerV = (unsigned int) (pModel->getHelicalVolume() - pModel->getHelicalVolume() * 0.29);
    upperV = (unsigned int) (pModel->getHelicalVolume() + pModel->getHelicalVolume() * 0.05);
    const float targetSubUnitVolume = (float) upperV / (float) pModel->getNumberOfSubUnits();

    auto lowerN = (unsigned int) std::round(lowerV * invBeadVolume);
    auto upperN = (unsigned int) std::round(upperV * invBeadVolume);

    // pick random number between lowerN and upperN
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    // create initial model
    // randomize bead indices
    // sort to workingLimit
    std::cout << "        WORKINGLIMIT SET : " << workingLimit << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";
    float subunit_test_volume, subunit_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float current_volume = subunit_volume;

    std::cout << "              CVX VOLUME : " << current_volume << std::endl;
    EulerTour eulerTour(bead_indices.cbegin(), workingLimit, pModel);

    unsigned int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();
    std::cout << "              EULER TOUR : " << currentNumberOfComponents << std::endl;

    pModel->setNumberOfSubUnitsHelical((unsigned int) std::ceil(pData->getDmax() / pModel->getRisePerSubUnit()));

    unsigned int violations, temp_violations;
    calculateModelPrDistributionSymHelical(&bead_indices, &binCount, workingLimit, violations, pModel, pData);

    // fill binCount for first time
    float lowestKL, currentKL = pData->getScore(binCount);

    float hlambda = pData->getIsIntensity() ? 10.0f : 0.001f;
    float muConstant = pData->getIsPr() ? 0.000001f : 0.1f; // chi tends to start in 100's whereas DKL is 0.1;

    const float invTarVolSubUnitmuConstant = muConstant / targetSubUnitVolume;

    unsigned int failures = 0;
    bool updateCVX = false;

    float testKL, test_energy; // sets alpha as a constant during high temp e
    float current_energy = currentKL + hlambda * (currentNumberOfComponents - 1) * (currentNumberOfComponents - 1);

    std::cout << "               => ENERGY : " << current_energy << std::endl;
    std::cout << "                    D_KL : " << currentKL << std::endl;
    std::cout << "                  VOLUME : " << current_volume << std::endl;
    std::cout << "              VIOLATIONS : " << violations << std::endl;
    std::cout << "                      WL : " << workingLimit << std::endl;


    char addRemoveText[50];
    unsigned int testIndex, totalViolations = violations; //(double)subUnitWorkingLimit;
    // what is tolerable number of violations per subsunit?
    current_energy += beta * totalViolations;

    float lowest_energy = current_energy;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup cop

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int high;
    auto lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000000001 : (double)highTempStartForCooling;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    float temp_subunit_volume_energy, current_subunit_volume_energy =
            std::abs(subunit_volume - targetSubUnitVolume) * invTarVolSubUnitmuConstant;

    current_energy += current_subunit_volume_energy;

//    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "initial_subunit_annealed", this, pData, 0, current_volume, 0);
//    pModel->writeHelicalSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "initial_hsym_annealed", this, pData, 0, current_volume, 0);
    std::vector<unsigned int> selections(workingLimit);
    for(unsigned int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    auto coupons = (unsigned int)std::lround((workingLimit*std::log((double)workingLimit) + 0.5772156649f*workingLimit + 0.5f));
    const unsigned int updateCount = 5*ccmultiple*coupons;
    float step_limit = (updateCount < 51357) ? 51357 : (float)updateCount;

    for (high = 0; high < step_limit; high++) { // iterations during the high temp search

        std::cout << "______________________________________________________________________________" << std::endl;

        if (distribution(gen) > 0.5) {
            std::sprintf(addRemoveText, "  ADJUST  ");
            priorRise = pModel->getRisePerSubUnit();
            priorTheta = pModel->getRotationPerSubUnit();

            //adjust delta and theta
            float plusMinus =  (binaryDistribution(gen) > 0) ? -1.0f : 1.0f;
            newRise = (binaryDistribution(gen) > 0) ? priorRise + plusMinus*delta_rise : priorRise;
            plusMinus =  (binaryDistribution(gen) > 0) ? -1.0f : 1.0f;
            newTheta = (binaryDistribution(gen) > 0) ? priorTheta + plusMinus*delta_theta : priorTheta;

            if (newRise < lowerRiseThreshold){
                newRise = lowerRiseThreshold;
            } else if (newRise > upperRiseThreshold){
                newRise = upperRiseThreshold;
            }


            if (newTheta < lowerThetaThreshold){
                newTheta = lowerThetaThreshold;
            } else if (newTheta > upperThetaThreshold){
                newTheta = upperThetaThreshold;
            }


            // determine
            pModel->setRiseAndRotationPerSubUnit(newRise, newTheta);
            pModel->setNumberOfSubUnitsHelical((unsigned int) std::ceil(pData->getDmax() / pModel->getRisePerSubUnit()));

            subunit_test_volume = current_volume;
            temp_subunit_volume_energy = current_subunit_volume_energy;
            tempNumberOfComponents = currentNumberOfComponents;

            calculateModelPrDistributionSymHelical(&bead_indices, &binCount, workingLimit, temp_violations, pModel, pData);
            testKL = pData->getScore(binCount);

            test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

            if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                pModel->setRiseAndRotationPerSubUnit(priorRise, priorTheta);
                pModel->setNumberOfSubUnitsHelical((unsigned int) std::ceil(pData->getDmax() / pModel->getRisePerSubUnit()));
            }

        } else {
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            if ( distribution(gen) < percentAddRemove ){
                // additional points to expand deadlimit will occur via enlarging CVX Hull
                if (distribution(gen) < 0.5637) { // ADD BEAD?
                    std::sprintf(addRemoveText, "    ADD   ");
                    unsigned int selectedIndex = 0;
                    unsigned int randomIndexToUse = bead_indices[selections[selectedIndex]];
                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomIndexToUse);
//
                    for(;selectedIndex < workingLimit; selectedIndex++){ // only add to bead that has fewer contacts?
                        if ((testIndex != noNeigborIndex)){ // && numberOfContactsFromSetSym(&beads_in_use_tree, pModel, randomIndexToUse) < 3){
                            break;
                        } else {
                            randomIndexToUse = bead_indices[selections[selectedIndex]];
                            testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomIndexToUse);
                        }
                    }

                    auto itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);
                    addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);

                    temp_violations = addToPrSymHelical(testIndex, bead_indices, workingLimit, binCount, pModel, pData);

                    testKL = pData->getScore(binCount);

                    subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel); // calculate volume of the subUnit
                    temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                    tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

                    if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.insert(testIndex);
                        updateCVX=true;
                    } else { // undo changes (rejecting)
                        auto beginBinCount = binCount.begin();
                        restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                        eulerTour.removeNode(testIndex);
                    }

                } else {
                    std::sprintf(addRemoveText, "  REMOVE  ");

                    testIndex = bead_indices[selections[0]];

                    // test for deletion
                    for(unsigned int i=1; i<workingLimit; i++){
                        tempNumberOfComponents = eulerTour.removeNode(testIndex);
                        if (tempNumberOfComponents <= currentNumberOfComponents){
                            break;
                        } else {
                            eulerTour.addNode(testIndex, pModel);
                            testIndex = bead_indices[selections[i]]; // potential to select the same thing twice
                        }
                    }

                    removeLatticePositionToModelSymHelical(bead_indices, binCount, &workingLimit, &testIndex, pModel, pData);
                    temp_violations = getViolationsHelical(&bead_indices, workingLimit, pModel, pData);

                    testKL = pData->getScore(binCount);

                    subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel); // calculate volume of the subUnit

                    temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;
                    if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.erase(testIndex);
                        updateCVX=true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        auto beginIt = bead_indices.begin();
                        auto beginBinCount = binCount.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt,
                                                              &workingLimit,
                                                              &binCountBackUp,
                                                              &beginBinCount);
                    }

                }

            } else {
                std::sprintf(addRemoveText, "POSITIONAL");

                unsigned int position = selections[0];
                unsigned int swap1 = bead_indices[position];

                if (distribution(gen) < (1.0f-acceptRate)){ // try to only move beads that are singly connected
                    bool found=true;
                    for(unsigned int i=1;i<workingLimit; i++){
                        if ( numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1) < 3 ){
                            if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                                found = false;
                                break;

                            }
                            //currentNumberOfComponents = eulerTour.removeNode(swap1);
                            eulerTour.addNode(swap1, pModel);
                        }

                        position = selections[i];
                        swap1 = bead_indices[position];
                    }

                    if (found){
                        position = selections[0];
                        swap1 = bead_indices[position];
                        for(unsigned int i=1;i<workingLimit; i++){
                            if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                                break;
                            } else { // if last node, will break out
                                eulerTour.addNode(swap1, pModel);
                                position = selections[i];
                                swap1 = bead_indices[position]; // what happens if I fail all the way?
                            }
                        }
                    }

                } else { // any position
                    swap1 = bead_indices[position];
                    for(unsigned int i=1;i<workingLimit; i++){
                        if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                        position = selections[i];
                        swap1 = bead_indices[position]; // what happens if I fail all the way?
                    }
                }

                auto itIndex = bead_indices.begin() + position;
                removeFromPrSymHelical(swap1, bead_indices, workingLimit, binCount, pModel, pData);

                beads_in_use_tree.erase(swap1);

                unsigned int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while (neighbor == noNeigborIndex || neighbor == swap1){
                    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }
                // make the swap, sort and update P(r)
                auto pSwap2 = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), neighbor);
                std::iter_swap(itIndex, pSwap2);
                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                temp_violations = addToPrSymHelical(neighbor, bead_indices, workingLimit, binCount, pModel, pData);

                testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry

                tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel); // calculate volume of the subUnit
                temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + temp_subunit_volume_energy;

                if ((test_energy  < current_energy) || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.insert(neighbor);
                    updateCVX = true;
                } else {
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    beads_in_use_tree.insert(swap1);
                    eulerTour.removeNode(neighbor);
                    eulerTour.addNode(swap1, pModel);
                }
            }
        }

        if (updateCVX){
            current_energy = test_energy;
            currentKL = testKL;
            current_volume = subunit_test_volume;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            current_subunit_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

            if (current_energy < lowest_energy){
                lowest_energy = current_energy;
            }

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASATemp(high, step_limit, acceptRate, lowTempStop, inv_kb_temp);

        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts
        std::printf("*******************             %s                 ******************* \n", addRemoveText);
        std::printf("         MAXSTEPS : %.0f (%4i) TEMP : %.3E\n", step_limit, high, lowTempStop);
        std::printf("           ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
        std::printf("            GRAPH : %3i\n", currentNumberOfComponents);
        std::printf("            LIMIT : UPPER %i LOWER %i => LIMIT : %5i \n", upperN, lowerN, workingLimit);
        std::printf("           VOLUME : %.0f (%.0f)  MU => %.4E  MU*VOL => %.3E\n", current_volume, targetSubUnitVolume, muConstant, current_subunit_volume_energy);
        std::printf("            HELIX : RISE => %.1f THETA %.1f \n", priorRise, priorTheta);
        std::printf("       VIOLATIONS : %d  BETA => %.4E  \n", totalViolations, beta);
        std::printf("             D_KL : %.4E ENRGY: %.4E (%.4E) \n", currentKL, current_energy, lowest_energy);
        std::cout << "*******************                                        *******************" << std::endl;
    }

    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_"+nameTo, this, pData, 0, current_volume, 0);
    pModel->writeHelicalSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "hsym_annealed_"+nameTo, this, pData, 0, current_volume, 0);

}


inline unsigned int Anneal::addToPrSymHelical(const unsigned int addMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){

    const unsigned int totalRotationalIndices = pModel->getNumberOfSubUnits();
    const unsigned int totalSubUnits = pModel->getTotalNumberOfHelicalSubUnits();
    const unsigned int totalSubUnitsInSingleFilament = pModel->getTotalNumberOfSubUnitsInSingleFiber();

    unsigned int findIt=0;
    const unsigned int * ptr = &beadsInUse.front();
    for(unsigned int i=0; i<workingLimit; i++){
        if (ptr[i] == addMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const int stopAt = findIt;

    /*
     * totalsubUnits and totalSubUnitsInSingleFilament can be equal but in general
     * totalSubUnits >= totalSubUnitsInSingleFilament
     */
    const unsigned int totalCoordinates = totalSubUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < workingLimit; i++){
        Bead * tempBead = pModel->getBead(beadsInUse[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }


    unsigned int count = workingLimit;
    for (unsigned int s=1; s < totalSubUnitsInSingleFilament; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesByHelicalSymmetry(s, workingLimit, count, coordinates);
    }


    if (totalRotationalIndices > 1){ // perform rotational symmetry of all elements up to totalSubUnitsInSingleFilament
        for(unsigned int i=0; i<totalRotationalIndices; i++){

        }
    }



    // adjust Pr
    float distance_to;
    for (unsigned int s=0; s < totalSubUnits; s++){  //
        const unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

    /*
     * correct for double counting
     */
    for (unsigned int s=0; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector
        const unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int ss=s+1; ss < totalSubUnits; ss++){
            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }


    /*
     * calculate violations only on the immediate neighbors
     *
     */
    unsigned int violation=0;
    unsigned int startAt = workingLimit, stopHere = 2*workingLimit;

    for (unsigned int i=0; i < workingLimit; i++) {

        const vector3 * ptempVec = &coordinates[i];

        for (unsigned int next_i = startAt; next_i < stopHere; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*ptempVec) - coordinates[next_i]).length();
            if (distance_to < violation_limit){
                violation++;
            }
        }
    }


    if (totalRotationalIndices > 1){ // perform rotational symmetry
        startAt = workingLimit*totalSubUnitsInSingleFilament;
        stopHere = workingLimit*totalSubUnitsInSingleFilament + workingLimit;

        for (unsigned int i=0; i < workingLimit; i++) {

            const vector3 * ptempVec = &coordinates[i];

            for (unsigned int next_i = startAt; next_i < stopHere; next_i++) {
                // calculate distance and convert to bin
                distance_to = ((*ptempVec) - coordinates[next_i]).length();
                if (distance_to < violation_limit){
                    violation++;
                }
            }
        }
    }
    return violation;
}


/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSymHelical(std::vector<unsigned int> *subUnit_indices, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, unsigned int &violation,  Model *pModel, Data *pData) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);

    const unsigned int totalRotationalIndices = pModel->getNumberOfSubUnits();
    const unsigned int totalsubUnits = pModel->getTotalNumberOfHelicalSubUnits();
    const unsigned int totalSubUnitsInSingleFilament = pModel->getTotalNumberOfSubUnitsInSingleFiber();

    const unsigned int totalCoordinates = totalsubUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;
    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }


    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < totalSubUnitsInSingleFilament; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesByHelicalSymmetry(s, indicesWorkingLimit, count, coordinates);
    }

    if (totalRotationalIndices > 1){ // perform rotational symmetry
        for(unsigned int i=0; i<totalRotationalIndices; i++){

        }
    }


    // calculate Pr, order is irrelavant
    float distance_to;
    for (unsigned int i=0; i < totalCoordinates; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (pData->convertToBin(distance_to) >= maxbin){
                std::cout << "EXCEED MAX BIN " << maxbin << " " << pData->convertToBin(distance_to) << std::endl;
                exit(0);
            }
        }
    }

    /*
     * for each lattice point within parent, determine how many violations it makes with other points in symmetry mates
     *
     * What counts as a violation?
     *
     * violation should prevent double counting
     *
     *
     */
    violation=0;
    unsigned int startAt = indicesWorkingLimit, stopAt = 2*indicesWorkingLimit;
    for (unsigned int i=0; i < indicesWorkingLimit; i++) {

        const vector3 * ptempVec = &coordinates[i];

        for (unsigned int next_i = startAt; next_i < stopAt; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*ptempVec) - coordinates[next_i]).length();
            if (distance_to < violation_limit){
                violation++;
            }
        }
    }


    if (totalRotationalIndices > 1){ // perform rotational symmetry
        startAt = indicesWorkingLimit*totalSubUnitsInSingleFilament;
        stopAt = indicesWorkingLimit*totalSubUnitsInSingleFilament + indicesWorkingLimit;

        for (unsigned int i=0; i < indicesWorkingLimit; i++) {

            const vector3 * ptempVec = &coordinates[i];

            for (unsigned int next_i = startAt; next_i < stopAt; next_i++) {
                // calculate distance and convert to bin
                distance_to = ((*ptempVec) - coordinates[next_i]).length();
                if (distance_to < violation_limit){
                    violation++;
                }
            }
        }
    }
}


/**
 *
 *
 */
unsigned int Anneal::removeLatticePositionToModelSymHelical(
        std::vector<unsigned int> & bead_indices,
        std::vector<unsigned int> & modelPrBins,
        unsigned int * pWorkingLimit,
        const unsigned int * pLatticePointToRemove, Model * pModel, Data *pData){
    auto pBeginIt = bead_indices.begin();
    auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    unsigned int violations = removeFromPrSymHelical(*pLatticePointToRemove, bead_indices, *pWorkingLimit, modelPrBins, pModel, pData);
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

/*
 * position to remove should be in beadsInUse
 *
 */
unsigned int Anneal::removeFromPrSymHelical(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){


    const unsigned int totalRotationalIndices = pModel->getNumberOfSubUnits();
    const unsigned int totalSubUnits = pModel->getTotalNumberOfHelicalSubUnits();
    const unsigned int totalSubUnitsInSingleFilament = pModel->getTotalNumberOfSubUnitsInSingleFiber();

    unsigned int findIt=0;
    // find index of bead to remove
    unsigned int * ptr = &beadsInUse.front();
    for(unsigned int i=0; i < workingLimit; i++){
        if (ptr[i] == removeMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const unsigned int stopAt = findIt;
    //
    // remove interdomain distances
    //
    float distance_to;

    const unsigned int totalCoordinates = totalSubUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;

    /*
     * build the primary filament
     * 1. add vectors of first subunit
     * 2. add additional subUnits using helical symmetry parameters
     */
    for (unsigned int i=0; i < workingLimit; i++){
        Bead * tempBead = pModel->getBead(beadsInUse[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }


    unsigned int count = workingLimit;
    for (unsigned int s=1; s < totalSubUnitsInSingleFilament; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesByHelicalSymmetry(s, workingLimit, count, coordinates);
    }

    /*
     * perform any rotational symmetry on the filament
     */
    if (totalRotationalIndices > 1){ // perform rotational symmetry
        for(unsigned int i=0; i<totalRotationalIndices; i++){

        }
    }

    // adjust Pr
    for (unsigned int s=0; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector
        const unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }


    /**
     *
     * we double remove(count) the symmetry related vectors, need to add back
     */
    for (unsigned int s=0; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector
        const unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int ss=s+1; ss < totalSubUnits; ss++){
            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }


    return 0;
}



/**
 * count the number of violations of parent subunit
 * @param bead_indices
 * @param beadIndiciesWorkingLimit
 * @param pModel
 * @param pData
 * @return
 */
unsigned int Anneal::getViolationsHelical(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, Model *pModel, Data *pData){

    const unsigned int totalRotationalIndices = pModel->getNumberOfSubUnits();

    const unsigned int totalCoordinates =  (totalRotationalIndices > 1) ? 3*indicesWorkingLimit : 2*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    Bead * tempBead;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < 2; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesByHelicalSymmetry(s, indicesWorkingLimit, count, coordinates);
    }

    /*
     * perform any rotational symmetry on the filament
     */
    if (totalRotationalIndices > 1){ // perform rotational symmetry
        for(unsigned int i=0; i<totalRotationalIndices; i++){

        }
    }


    float distance_to;
    unsigned int violations = 0;
    /*
     * do any beads collide with primary subunit?
     */
    for(unsigned int i=0; i<indicesWorkingLimit; i++) {

        const vector3 *baseVec = &coordinates[i];
        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = ((*baseVec) - coordinates[s]).length();
            if (distance_to < violation_limit){
                violations++;
            }
        }
    }


    return violations;
}