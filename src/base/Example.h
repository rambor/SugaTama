//
// Created by xos81802 on 21/08/2018.
//

#ifndef SUGATAMA_EXAMPLE_H
#define SUGATAMA_EXAMPLE_H

#include <string>

class Example {

    const std::string type;

    void printAnchor();
    void printHelical();
    void printManual();

public:
    Example(std::string what);
};


#endif //SUGATAMA_EXAMPLE_H
