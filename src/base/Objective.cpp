//
// Created by xos81802 on 30/04/2018.
//

#include "Objective.h"

using namespace boost::filesystem;

Objective::Objective(){

}

/**
 * returns a pointer to Data Object
 * holds the memory address to element in datasets vector
 *
 */
Data * Objective::getMainDataset() {

    Data * preturnMe = &(datasets.begin()->second);

    int phaseCount = 0;

    if (datasets.size() > 1){
        // iterate over each dataset
        // determine total number of phases in each dataset
        // identify dataset with largest number of phases
        for(auto & iterator : datasets) {
            // getData object
            Data * temp = &iterator.second;
//            if (temp.getTotalPhases() > phaseCount){
//                preturnMe = &(iterator->second);
//                phaseCount = temp.getTotalPhases();
//            }
        }

    }

    return preturnMe;
}


/**
 *
 */
//void Objective::createWorkingSets(Model & model) {
//
//    int totalDatasets = datasets.size();
//
//    // generates different sets of q-values per dataset
//    for(int i=0; i<totalDatasets; i++){
//        this->getDataObjectByIndex(i)->creatingWorkingSet(model);
//    }
//
//}


/**
 * Create dataset object using both Intensity (Reciprocal Space) and Real-Space datasets
 * Modeling is based on real-space data.
 *
 * Why keep Reciprocal?
 *
 * @param iofqfile
 * @param pofrfile
 */
void Objective::addDataObject(std::string datFile) {
    // load file and create Data object
    // get name of iofqfile, strip away

    path p = datFile;
    boost::regex slashes("/|\\\\");

    std::string identifierKey;
    std::vector<std::string> tempLine;

    try {
        if (exists(p)){
            //split name if forward or backslash is present
            if(boost::regex_search(datFile, slashes)){

                boost::split(tempLine, datFile, boost::is_any_of("/|\\\\"), boost::token_compress_on);

                //int last = tempLine.size()-1;
                identifierKey = tempLine[tempLine.size()-1];

            } else {
                identifierKey = datFile;
            }

            // create IofQ object
            datasets.emplace(identifierKey,Data(datFile));
            //datasets.insert(std::pair<std::string, Data> (identifierKey, Data(datFile)));
            keys.push_back(identifierKey);

            // Add PofR Dataset
            //Data * pdata;
            //pdata = &datasets[identifierKey];
            //(*pdata).addPofRData(pofrfile);
            //datasets[identifierKey].addPofRData(pofrfile);
        }

    } catch (const filesystem_error& ex) {
        std::cout << ex.what() << '\n';
    }

    /*
    cout  <<  "  root_name()----------: " << p.root_name() << '\n';
    cout  <<  "  root_directory()-----: " << p.root_directory() << '\n';
    cout  <<  "  root_path()----------: " << p.root_path() << '\n';
    cout  <<  "  relative_path()------: " << p.relative_path() << '\n';
    cout  <<  "  parent_path()--------: " << p.parent_path() << '\n';
    cout  <<  "  filename()-----------: " << p.filename() << '\n';
    */
}
