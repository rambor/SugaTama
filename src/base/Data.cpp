//
// Created by xos81802 on 30/04/2018.
//

#include <random>
#include "Data.h"


Data::Data(std::string filename) : filename(filename) {

    std::vector<Datum> values;
    //open file and extract intensities and pofr distribution

    intensities.reserve(1800);
    working_probability_per_bin.resize(23);

    // determine if dat file, it is IofQ or PofR
    std::string ext = this->getFileExt(filename);

    // check if file exists
    if (!boost::filesystem::exists(filename)){
        throw std::invalid_argument("** ERROR FILE => DATA FILE NOT FOUND : " + filename);
    }


    if (ext == "dat"){
        if (this->checkIofQFile(filename)){
            this->extractIntensityData();
        } else if (this->checkPofRFile(filename)){
            this->extractPrData();
            working_probability_per_bin.cbegin(); // what the hell is this for?
        }
    }
}


Data::Data(std::string filename, float fqmax, float dmax) : filename(filename){
    this->qmax = fqmax;
    shannon_bins = (float)ceil(dmax*qmax/M_PI)*1.0f;  // shannon bins
    ns_dmax = (float)(shannon_bins*M_PI/qmax);     // dmax corresponding to Shannon Bin

    bin_width = M_PI/fqmax;
}


bool Data::checkPofRFile(std::string file) {
// read in file

    bool returnMe = false;
    std::ifstream data (file, std::ifstream::in);
    if (data.is_open()) {

        boost::regex format("REAL SPACE");
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

            if(boost::regex_search(line, format)){
                logger("READING PR VALUES", filename);
                isIntensity = false;
                isPr = true;
                returnMe = true;
                break;
            }
        }
    }

    data.close();
    return returnMe;
}


bool Data::checkIofQFile(std::string file) {

    bool returnMe = false;
    // read in file
    std::ifstream data (file, std::ifstream::in);
    if (data.is_open()) {

        boost::regex format("I_OBS");
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

            if(boost::regex_search(line, format)){
                std::cout << " Extracting Intensities : " << filename << std::endl;
                isIntensity = true;
                isPr = false;
                returnMe = true;
                break;
            }
        }
    }


    data.close();
    return returnMe;
}

void Data::extractIntensityData() {

    std::ifstream data (this->filename, std::ifstream::in);

    if (data.is_open()) {

        boost::regex dataFormat("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
        boost::regex volFormat("(volume)|(VOLUME)");
        boost::regex rgFormat("(RECI rg)|(RECI RG)|(RECI Rg)", boost::regex::icase);
        boost::regex dmaxFormat("(dmax)|(DMAX)", boost::regex::icase);
        boost::regex errorFormat("ERROR", boost::regex::icase);
        std::string line;
        float qvalue, iofq, sigma;
        totalDataPoints = 0;
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

            if ((line.c_str()[0] != '-') && (line.length() > 0 && boost::regex_search(tempLine.at(0), dataFormat)) && boost::regex_search(tempLine.at(1), dataFormat)  ) {

                qvalue = std::stof(tempLine[0]);
                iofq = std::stof(tempLine[1]);
                sigma = std::stof(tempLine[2]);

                intensities.emplace_back(Datum(qvalue, iofq, sigma, totalDataPoints));
                totalDataPoints++;

            } else if (boost::regex_search(line, volFormat)){

                std::vector<std::string> volLine;
                boost::split(volLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                long int index = findInArray("VOLUME", &volLine); // need to make case insensitive
                if (index > 0){
                    this->volume = std::stof(volLine[index + 2]);
                    std::cout << " => VOLUME " << this->volume << std::endl;
                }
                //cout << " find in vector: => " << findInArray("VOLUME", &volLine) << endl;
                //this->volume = stof(volLine[4]);
            } else if (boost::regex_search(line, rgFormat) && !(boost::regex_search(line, errorFormat))){

                std::vector<std::string> rgLine;
                boost::split(rgLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

                long int index = findInArray("RECI", &rgLine); // need to make case insensitive
                if (index > 0){
                    this->rg = std::stof(rgLine[index + 3]);
                }

                std::cout << " => RECI Rg " << this->rg << std::endl;
            } else if (boost::regex_search(line, dmaxFormat)) {
                std::vector<std::string> dmaxLine;
                boost::split(dmaxLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->dmax = std::stof(dmaxLine[4]);
                std::cout << " => DMAX " << this->dmax << std::endl;
            }
        }
        qmax = intensities[totalDataPoints-1].getQ();
        /*
         * bin_width =>   PI*n_s = q_max*d_max
         *           => PI/q_max = d_max/n_s
         */
        bin_width = (float)M_PI/qmax;
        /*
         * need a n_s estimate from user for the chi-free cross-validation
         *
         */
        shannon_bins = (float)ceil(dmax*qmax/M_PI)*1.0f;  // shannon bins

        if (shannon_bins < 4){
            throw std::invalid_argument("  ERROR => too few shannon bins");
        }

        // create workingset for analysis
        setDataBinSize(1);
        createCrossValidationDataset();

        /*
         * set target function
         */
    }
}


void Data::extractPrData() {

    std::ifstream data (this->filename, std::ifstream::in);

    unsigned int tempCount = 0;
    if (data.is_open()) {

        boost::regex dataFormat("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
        boost::regex dmaxFormat("(dmax)|(DMAX)", boost::regex::icase);
        boost::regex qmaxFormat("(qmax)|(QMAX)", boost::regex::icase);
        boost::regex volFormat("(volume)|(VOLUME)", boost::regex::icase);
        boost::regex porodFormat("POROD", boost::regex::icase);
        boost::regex raveFormat("<r>");
        //boost::regex rgFormat("(RECI rg)|(RECI RG)|(RECI Rg)", boost::regex::icase);
        boost::regex rgFormat("REAL Rg :", boost::regex::icase);

        std::string line;
        float rvalue, pvalue, sigma;

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
            if ((line.c_str()[0] != '-') && (line.length() > 0 && boost::regex_search(tempLine.at(0), dataFormat)) && boost::regex_search(tempLine.at(1), dataFormat)  ) {

                rvalue = std::strtof(tempLine[0].c_str(), nullptr);
                pvalue = std::strtof(tempLine[1].c_str(), nullptr);
                sigma = std::strtof(tempLine[2].c_str(), nullptr);

                pofr.emplace_back(RealSpace(rvalue, pvalue, sigma));
                tempCount++;

            } else if (boost::regex_search(line, dmaxFormat)){
                std::vector<std::string> dmaxLine;
                boost::split(dmaxLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->dmax = std::strtof(dmaxLine[4].c_str(), nullptr);
            } else if (boost::regex_search(line, qmaxFormat)){
                std::vector<std::string> qmaxLine;
                boost::split(qmaxLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->qmax =std::strtof(qmaxLine[6].c_str(), nullptr);
            } else if (boost::regex_search(line, volFormat)){
                std::vector<std::string> volLine;
                boost::split(volLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                if (volLine[3] == ":"){
                    this->volume = std::strtof(volLine[4].c_str(), nullptr);
                }
            } else if (boost::regex_search(line, raveFormat)){
                std::vector<std::string> raveLine;
                boost::split(raveLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->rave = std::strtof(raveLine[5].c_str(), nullptr);
            } else if (boost::regex_search(line, rgFormat)){
                std::vector<std::string> rgLine;
                boost::split(rgLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

                this->rg = std::strtof(rgLine[5].c_str(), nullptr);
            }
        }

        data.close();

        if (dmax <= 0 || qmax < 0.1){
            // throw exception
            throw std::invalid_argument("QMAX/DMAX NOT FOUND OR EQUAL TO ZERO");
        }

        //
        // shannon_bins = round(dmax*qmax/M_PI); // based on the data
        // qmax is specified in pr.dat file
        // shannon bins can be <= Moore Coefficents, e.g. downsampling the distribution
        //

        shannon_bins = (float)ceil(dmax*qmax/M_PI)*1.0f;  // shannon bins
        ns_dmax = (float)(shannon_bins*M_PI/qmax);     // dmax corresponding to Shannon Bin

        //bin_width = M_PI/qmax;

        if (shannon_bins < 3){// throw exception
            throw std::invalid_argument("Too few shannon bins, USE ANOTHER PROGRAM");
        }

        // use experimental Shannon number before scaling in setBinSize

        this->parseBins();


        std::cout << " => EXTRACTED PARAMS " << std::endl;
        logger("VOLUME", std::to_string((int)volume));
        logger("DMAX", formatNumber(dmax,1));
        logger("SHANNON DMAX", formatNumber(ns_dmax,1));
        logger("QMAX", formatNumber(qmax,5));
        logger("RG", formatNumber(rg,2));


        if (bin_coefficients.size() > 2){

            if (shannon_bins > bin_coefficients.size()){
                std::cout << "-----------------------     WARNING     -----------------------" << std::endl;
                std::cout << "-- qmax NOT SUPPORTED by Shannon Number and bin size in file          --" << std::endl;
                std::cout << "-- Validate input files                                         --" << std::endl;
                std::cout << "-- Shannon_BINS : " << shannon_bins << " != " << bin_coefficients.size() << std::endl;
                exit(0);
            }

            //bin_width = (ns_dmax/(float)bin_coefficients.size());
            bin_width = (dmax/(float)bin_coefficients.size());
            shannon_bins = bin_coefficients.size();
            normalizePrBins();

        } else { // use this for data that is not already binned by scatter
            this->normalizePofR(tempCount);
        }

        invBinSize = 1.0f/shannon_bins;
    }
}

/**
 * creates Shannon limited binning of input P(r) distribution
 * use trapezoid rule to integerate for normalization
 * not sure if I should just interpolate the middle point of the bin using interpolation theory rather than average
 * param int count is the number of r-values in P(r) file
 */
void Data::normalizePofR(unsigned int count) {
    // float sum = pofr[0].pr + pofr[count-1].pr; first and last are 0 by default
    float sum;
    // integrate per bin
    probability_per_bin.reserve((unsigned int)shannon_bins);
    // round up when calculating shannon number and use the d_max from the round up.
    float r_value;
    float totalSum=0.0, lower, upper;

    for(unsigned int i=0; i<shannon_bins; i++){ // integrate between lower and upper
        lower = i*bin_width;     //q-value
        upper = (i+1)*bin_width; //q-value
        // cc = 0; // counts ticks
        // integrate bin (area of bin)
        // interpolate P(r) values for lower and upper
        unsigned int lowerIndex=0, upperIndex=count;
        float lowerRValue, upperRValue;

        for(unsigned int j=0; j<count; j++){
            if(pofr[j].r > lower){ // equal when r is zero
                lowerIndex = j-1;  // index of first value greater than lower bound on bin
                break;
            }
        }

        for(unsigned int j=0; j<count; j++){
            upperIndex = j;
            if(pofr[j].r >= upper){ // equal when r is dmax
                //upperIndex = j; // index of first value greater than upper bound on bin
                break;
            }
        } // if upper exceeds dmax, upperIndex is last element

        // perform linear interpolation of missing value
        lowerRValue = 0;
        if (lowerIndex > 0){
            lowerRValue = pofr[lowerIndex].pr + (pofr[lowerIndex+1].pr - pofr[lowerIndex].pr)*(lower - pofr[lowerIndex].r)/(pofr[lowerIndex+1].r- pofr[lowerIndex].r);
        }

        upperRValue = 0;
        if (upperIndex < count){
            upperRValue = pofr[upperIndex-1].pr + (pofr[upperIndex].pr - pofr[upperIndex-1].pr)*(upper - pofr[upperIndex-1].r)/(pofr[upperIndex].r- pofr[upperIndex-1].r);
        }

        // do integration of bin
        sum = 0.0;
        for (unsigned int j=lowerIndex; j<count; j++){
            r_value = pofr[j].r;
            unsigned int lastIndex = 0;
            if (r_value > lower && r_value < upper){
                float height = std::min(lowerRValue, pofr[j].pr);
                sum += height*(r_value-lower) + abs(lowerRValue - pofr[j].pr)*(r_value-lower)*0.5;

                lower = r_value;
                lowerRValue = pofr[j].pr;
                lastIndex = j;
            }

            if (r_value >= upper) { // add last trapezoid to area
                // area of rectable + area of triangle
                // sum += (upper - lower)*min(upperRValue, pofr[lastIndex].pr) + abs(upperRValue - pofr[lastIndex].pr)*(r_value-lower)*0.5;
                sum += (upper - lower)*std::min(upperRValue, pofr[lastIndex].pr) + abs(upperRValue - pofr[lastIndex].pr)*(upper-lower)*0.5;
                break;
            }
        }

        probability_per_bin.push_back(sum);
        totalSum += sum;
    }

    // trapezoid rule (integrate Pr-distribution to calculate normalization constant)
    // assume last r-value is 0 at d_max
    float partialSum=0.0f;
    for(unsigned int i=0; i< count; i++){
        partialSum += 2*pofr[i].pr;
    }

    float invPartialSum = 1.0f/totalSum;

    normalize(invPartialSum);
}

long int Data::findInArray(const std::string & value, std::vector<std::string> * strings) {

    long int index = -1;
    auto beginIt = strings->begin();
    auto it = std::find(beginIt, strings->end(), value);
    if (it != strings->end()){
        index = std::distance(beginIt, it);
    }
    return index;
}


void Data::normalize(float norm){

    FILE * pFile;

    const char *outputFileName;
    std::string nameOf = "normalized_expPr.txt";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    // Normalilize
    std::cout << " NORMALIZATION CONSTANT " << norm << std::endl;
    std::cout << " NORMALIZED P(r) : " << std::endl;
    fprintf(pFile, "# NORMALIZED P(r)\n");
    fprintf(pFile, "# REMARK  BINWIDTH => %.3f \n", bin_width);
    fprintf(pFile, "# REMARK      DMAX => %.0f \n", dmax);
    fprintf(pFile, "# REMARK   NS DMAX => %.0f \n", ns_dmax);
    fprintf(pFile, "  %.3f %.8f\n", 0.0, 0.0);
    for(unsigned int i=0; i<probability_per_bin.size(); i++){
        probability_per_bin[i] = probability_per_bin[i]*norm;
        //printf("  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]);
        fprintf(pFile, "  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]); // should match bin -> P(r) mapping in input file
    }

    fclose(pFile);
}

/**
 * calculate PofR using Moore coefficients and normalize to 1
 */
void Data::normalizePrBins() {

    float norm=0;

    for (auto i : bin_coefficients){
        norm += i;
    }

    // can I down-sample distribution ?
    probability_per_bin.reserve(bin_coefficients.size());
    // bin_width = dmax/total_bins;
    // round up when calculating shannon number and use the d_max from the round up.

    // for each bin, calculate area
    float value, invNorm = 1.0f/norm;
    std::cout << "  => NORMALIZED BINNED VALUES " << std::endl;
    std::cout << "  => Q(SCATTERING VECTOR) SHOULD BE INV ANGSTROMS " << std::endl;
    std::cout << "  => RVALUE (ANGSTROMS) " << std::endl;
    logger("RVALUE (BIN)", "HEIGHT");
    for(unsigned int i=0; i<bin_coefficients.size(); i++){
        // integrate between lower and upper
        value = bin_coefficients[i]*invNorm;
        probability_per_bin.push_back(value);
        logger(formatNumber((bin_width*(0.5+i)),2), formatNumber(value,6));
    }
}

void Data::parseBins() {

    std::ifstream data (this->filename.c_str());
    // boost::regex binValueLine("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
    boost::regex coefficient("BIN_[0-9]+");
    boost::regex remarkFormat("REMARK");

    int binCount=0;
    logger("", "PARSING BINS");

    if (data.is_open()) {

        while(!data.eof()) {
            std::string line;
            std::vector<std::string> tempLine;

            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if(boost::regex_search(tempLine.at(0), remarkFormat) && boost::regex_search(line, coefficient) ){
                std::vector<std::string> mline;
                boost::trim(line);
                boost::split(mline, line, boost::is_any_of("\t  "), boost::token_compress_on);
                unsigned long int total = mline.size();

                try {
                    float value = std::abs(stof(mline.back()));

                    if(std::strcmp(mline[total-2].c_str(), ":")==0 && value > 0){
                        bin_coefficients.push_back(stof(mline.back()));
                        binCount++;
                    } else {
                        throw std::invalid_argument( "Improper Bin Value at : \n\t" + line  + " \n Can not be zero or negative");
                    }

                } catch (std::exception &err) {
                    std::cerr<<"Caught "<<err.what()<<std::endl;
                    std::cerr<<"Type "<<typeid(err).name()<<std::endl;
                    exit(0);
                };
            }
        }

    }

    logger("TOTAL BINS READ", std::to_string(bin_coefficients.size()));
    data.close();
}



std::string Data::getFileExt(const std::string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
        return(s.substr(i+1, s.length() - i));
    }

    return("");
}


/**
 *  Generate Randomly selected points from input intensity data set
 *  Based on pointsPerBin
 */
void Data::createCrossValidationDataset(){

    std::cout << "--------------  CREATING CROSS-VALIDATION SET     -------------" << std::endl;
    std::random_device rd;
    std::mt19937 gen(rd());

    auto total = (unsigned int) pointsPerBin.size();
    unsigned int toGet, index=0;

    std::vector<unsigned int> testSelection;

    qvalues.clear();
    workingSet.clear();

    std::srand ( unsigned ( std::time(0) ) );
    float lower, upper, deltaQ = qmax/(float)total;
    for(unsigned int i=0; i < total; i++){

        testSelection.clear();
        toGet = pointsPerBin[i];

        lower = lowerBounds[i];
        upper = lower + deltaQ;

        // need indices within the lower and upper bound

        for(unsigned int j=index; j < totalDataPoints; j++){
            Datum temp = intensities[j];

            if ((temp.getQ() >= lower) && (temp.getQ() < upper)) {
                testSelection.push_back(j);
            } else if ((temp.getQ() >= upper)){
                index=j;
                break;
            }
        }

        // if no points in bin move to next bin and see if there are points
        std::shuffle ( testSelection.begin(), testSelection.end(), gen );
        std::sort (testSelection.begin(), testSelection.begin()+toGet);

        // using built-in random generator:
        // randomly sample from testSelection based on toGet
        for(unsigned int m=0; m < toGet; m++){
            workingSet.push_back(intensities[testSelection[m]]);
        }
        std::printf("   Q-RANGE : %.5f < q <= %.5f N %d (%d) \n", lower, upper, toGet, (int)testSelection.size());
    }

    auto workingSetSize = (unsigned int)workingSet.size();
    std::cout << "        Total I(q) : " << totalDataPoints << std::endl;
    std::cout << "  working set size : " << workingSetSize << std::endl;
    // create vector of qvalues
    qvalues.resize(workingSetSize);
    invVarianceWorkingSet.resize(workingSetSize);

    for(unsigned int m=0; m < workingSetSize; m++){
        qvalues[m] = workingSet[m].getQ();
        invVarianceWorkingSet[m] = workingSet[m].getInvVar();
    } // write workingSet to File

}

/**
 * bins can be 1 or more times the number of shannon bins
 */
void Data::setDataBinSize(unsigned int bins) {

    lowerBounds.clear();
    auto data_bins = (unsigned int)round(bins*shannon_bins);

    float binwidth = qmax/(float)data_bins;
    std::vector<float> shValues(data_bins);

    unsigned int index = 0;
    float sumSN, averageNoise, invCount;

    // determine points per bin
    float lower, upper;

    for (unsigned int i=0; i < data_bins; i++){
        lower = binwidth*i;
        upper = binwidth*(i+1);
        float shannonHartley;
        unsigned int count = 0;
        sumSN = 0;
        // go through data and count and determine average noise in bin

        while( (intensities[index].getQ() >= lower && intensities[index].getQ() < upper) && index < totalDataPoints ){
            Datum temp = intensities[index];
            sumSN += temp.getI()/temp.getSigma();
            count++; // number of points within range
            index++;
        }

        invCount = 1.0f/(float)count;
        averageNoise = sumSN*invCount;

        shannonHartley = (float)(2.0f*M_PI/(float)dmax*log(1.0f + averageNoise));
        shValues[i] = shannonHartley;

        if (count > 0){
            pointsPerBin.push_back(getPointsPerBin(averageNoise, count));
            lowerBounds.push_back(lower);
        }
    }
}

unsigned int Data::getPointsPerBin(float ioversigma, unsigned int totalInBin) {

    float value = 2.0f;
    if (ioversigma < 2){
        value = 0.80f*totalInBin;
    } else if (ioversigma > 100) {
        value = 0.1f*totalInBin;
    } else if (ioversigma >= 2 && ioversigma < 10){
        value = 0.47f*totalInBin;
    } else if (ioversigma >= 10 && ioversigma < 30) {
        value = 0.31f*totalInBin;
    } else if (ioversigma >= 30 && ioversigma < 60){
        value = 0.23f*totalInBin;
    } else if (ioversigma >= 60 && ioversigma < 100){
        value = 0.13f*totalInBin;
    }

    return (unsigned int)value;
}


bool Data::getIsIntensity() const {
    return isIntensity;
}

bool Data::getIsPr() const {
    return isPr;
}


double Data::getScore(std::vector<unsigned int> &modelPR){
    return (*score)(modelPR);
}


void Data::setScoringFunction(unsigned int maxBin){

    if (isIntensity){
        /*
         * calculate sin(qr)/qr for all bins
         */
        sinc_qr.clear();
        unsigned int workingSetSize = qvalues.size();

        std::cout << "  SCORING FUNCTION => CHI-SQUARE " << std::endl;
        for (unsigned int q = 0; q < workingSetSize; q++) {
            float qValue = qvalues[q];
            for (unsigned int i=0; i < maxBin; i++) { // for fixed q value, iterate over all possible bins in P(r)
                float rvalue = bin_width*i + 0.5f*bin_width;
                float qr = qValue*rvalue;
                sinc_qr.push_back(sin(qr)/qr);
                // std::cout<< rvalue << " " << q << " " << qr << " " << sin(qr)/qr << std::endl;
            }
        }

        score = new ReciSpaceScore(maxBin, workingSet, sinc_qr, invVarianceWorkingSet);

        probability_per_bin.resize(maxBin);
        working_probability_per_bin.resize(maxBin);
        std::fill(working_probability_per_bin.begin(), working_probability_per_bin.end(), 0);
        std::fill(probability_per_bin.begin(), probability_per_bin.end(), 0);

    } else if (isPr){ // effectively creates working distribution

        working_probability_per_bin.resize(maxBin);
        std::fill(working_probability_per_bin.begin(), working_probability_per_bin.end(), 0);
        /*
         * probability_per_bin set by bin_coefficients.size()
         *
         * maxBin must be > bin_coefficients.size()
         */
        if (probability_per_bin.size() > working_probability_per_bin.size()){
            std::cout << "ERROR => PR DATA BINS : " << probability_per_bin.size() << "> WORKING BINS : " << maxBin << std::endl;
            exit(0);
        }

        std::copy(probability_per_bin.begin(), probability_per_bin.end(), working_probability_per_bin.begin());

        //last nonzero bin
        zeroBin=0;
        for (unsigned int i=0; i<maxBin; i++){ // if maxbin is greater than probability_per_bin.size, we have an empty bin

            if (working_probability_per_bin[zeroBin] <= 0){ // maybe this should be flt_epsilon
                break;
            }
            zeroBin++;
        }

        score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
        //Score realSpaceScore = RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
        //score = &realSpaceScore;
    }
}



void Data::printKLDivergence(std::vector<unsigned int> &modelPR){

    float totalCounts = 0.0;
    float prob;
    int totalm = modelPR.size();

    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    for (int i=0; i<totalm; i++){
        totalCounts += modelPR[i];
    }

    float tempPR, r;
    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    std::cout << "FINAL MODEL " << std::endl;
    std::cout << "       r       MODEL      EXP        D_KL" << std::endl;
    for (unsigned int i=0; i < shannon_bins; i++){
        prob = probability_per_bin[i];  // bounded by experimental Shannon Number
        tempPR = modelPR[i];
        r = bin_width*i+0.5*bin_width;
        printf(" %10.4f  %.6f  %.6f  % .6f\n", r, prob, tempPR/totalCounts, (prob * log(prob/tempPR*totalCounts)));
        //cout << r << " " << prob << " " << tempPR/totalCounts << " " << prob * log(prob/tempPR*totalCounts) << endl;
    }
}

void Data::printICalc(std::vector<unsigned int> &modelPR){

    // calculate chi2
    unsigned int totalInWorkingSet= workingSet.size();
    unsigned int maxBin = modelPR.size();

    std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
    double totalCounts = 0.0d;
    for (auto & pr : modelPR_float) {
        totalCounts += pr;
    }
    double invTotalCounts = 1.0d/totalCounts;

    FILE * pFile;

    const char *outputFileName;
    std::string nameOf = "icalcf.txt";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    // for each q-value, calculate I_calc from P(r)
    for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {

        float qValue = qvalues[qr];
        float iofq=0;
        unsigned int count=0;
        for (auto pr : modelPR_float) {
            iofq += sinc_qr[qr*maxBin + count]*pr*invTotalCounts;
            count++;
        }
        fprintf(pFile, "%.5E  %.5E \n", qValue, iofq);
    }
    //fprintf(pFile, "END\n");
    fclose(pFile);
}

Data::Data(const Data &dataToCopy) {
    filename = dataToCopy.filename;
    totalDataPoints = dataToCopy.totalDataPoints;
    totalPhases = dataToCopy.totalPhases;
    zeroBin = dataToCopy.zeroBin;
    volume = dataToCopy.volume;
    rg = dataToCopy.rg;
    rave = dataToCopy.rave;
    qmin = dataToCopy.qmin;
    qmax = dataToCopy.qmax;
    ns_dmax = dataToCopy.ns_dmax;
    dmax = dataToCopy.dmax;
    shannon_bins = dataToCopy.shannon_bins;
    weight = dataToCopy.weight;
    bin_width = dataToCopy.bin_width;
    invBinSize = dataToCopy.invBinSize;
    isIntensity = dataToCopy.isIntensity;
    isPr = dataToCopy.isPr;

    intensities = std::vector<Datum>(dataToCopy.intensities);
    workingSet = std::vector<Datum> (dataToCopy.workingSet);
    invVarianceWorkingSet = std::vector<float> (dataToCopy.invVarianceWorkingSet);
    pointsPerBin= std::vector<unsigned int> (dataToCopy.pointsPerBin);
    lowerBounds= std::vector<float> (dataToCopy.lowerBounds);
    qvalues= std::vector <float> (dataToCopy.qvalues);
    probability_per_bin= std::vector<float> (dataToCopy.probability_per_bin);
    bin_coefficients= std::vector<float> (dataToCopy.bin_coefficients);
    working_probability_per_bin= std::vector<float> (dataToCopy.working_probability_per_bin);
    sinc_qr= std::vector<float> (dataToCopy.sinc_qr);
    pofr= std::vector<RealSpace> (dataToCopy.pofr);


    if (isIntensity){
        score = new ReciSpaceScore(probability_per_bin.size(), workingSet, sinc_qr, invVarianceWorkingSet);
    } else if (isPr){
        score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
    }

}
