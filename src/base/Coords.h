//
// Created by xos81802 on 13/07/2018.
//

#ifndef SUGATAMA_COORDS_H
#define SUGATAMA_COORDS_H


#include <iostream>
#include <cstdio>
#include <string>


class Coords {

public:
    Coords();
    Coords(float x, float y, float z, std::string atomType, float occ);
    Coords(const Coords &e2);

    Coords & operator=(const Coords & tourtomove) {
        Coords tmp(tourtomove);
        return *this;
    }


    float x;
    float y;
    float z;
    std::string type;
    float occ;

};


#endif //SUGATAMA_COORDS_H
