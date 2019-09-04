//
// Created by xos81802 on 13/07/2018.
//

#include "Coords.h"

using namespace std;

Coords::Coords(){
    this->x = 0;
    this->y = 0;
    this->z = 0;
    this->type = "";
    this->occ = 0;
}

Coords::Coords(float x, float y, float z, std::string atomType, float occ) : x(x), y(y), z(z), type(atomType), occ(occ) {
}

Coords::Coords(const Coords &e2) : Coords(e2.x, e2.y, e2.z, e2.type, e2.occ) {

}
