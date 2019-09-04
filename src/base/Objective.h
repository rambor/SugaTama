//
// Created by xos81802 on 30/04/2018.
//

#ifndef SUGATAMA_OBJECTIVE_H
#define SUGATAMA_OBJECTIVE_H

#include "Data.h"

class Objective {
    // HASH KEY => DATA
    std::map <std::string, Data> datasets;
    std::vector<std::string> keys;
    std::vector<float> weights;

public:
    Objective();
    ~Objective(){
//        for(auto & it :datasets){
//            datasets.erase(it.first);
//        }
        datasets.clear();
    }

    void addDataObject(std::string datfile);
    Data * getMainDataset();
};


#endif //SUGATAMA_OBJECTIVE_H
