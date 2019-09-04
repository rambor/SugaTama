//
// Created by xos81802 on 10/07/2018.
//

#ifndef SUGATAMA_BEAD_H
#define SUGATAMA_BEAD_H


#include "vector3.h"

class Bead {
private:
    float x;
    float y;
    float z;
    vector3 vec;
    float contrast;

    /**
      * octant 1 (q1) (x+y+z)
      * octant 2 (q2) (x-y-z)
      * octant 3 (q3) (x+y-z)
      * octant 4 (q4) (x-y+z)
      *
      */

    float q1, q2, q3, q4;

public:
    Bead(){};
    Bead(float x, float y, float z, float contrast): x(x), y(y), z(z), vec(x,y,z), contrast(contrast){
        /**
         * (x+y+z)
         * (x-y-z)
         * (x+y-z)
         * (x-y+z)
         * octant 1 (q1) (x+y+z)
         * octant 2 (q2) (x-y-z)
         * octant 3 (q3) (x+y-z)
         * octant 4 (q4) (x-y+z)
         */
        q1 = x+y+z;
        q2 = x-y-z;
        q3 = x+y-z;
        q4 = x-y+z;
    }

    Bead(const Bead & other) : x(other.x), y(other.y), z(other.z), vec(other.x,other.y,other.z), contrast(other.contrast) {
        q1 = other.q1;
        q2 = other.q2;
        q3 = other.q3;
        q4 = other.q4;
    }

    ~Bead(){

    }

    float const getContrast() const {return contrast;} // likely relative to protein ?
    void setContrast(float contrast);

    float const & getX() const {return vec.x;}
    float const & getY() const {return vec.y;}
    float const & getZ() const {return vec.z;}

    /**
     * x+y+z
     */
    float const & getQ1() const {return q1;}

    /**
     * x-y-z
     */
    float const & getQ2() const {return q2;}

    /**
    * x+y-z
    */
    float const & getQ3() const {return q3;}

    /**
     * x-y+z
     */
    float const & getQ4() const {return q4;}

    vector3 const & getVec() const {return vec;}

    const vector3 & getVec() {return vec;}

    void translateTo(const vector3 & tvec){
        vec += tvec;
        x = vec.x;
        y = vec.y;
        z = vec.z;
    }

    void printCoordinates(){ std::cout << "COORDINATE " << this->x << " " << this->y << " " << this->z << std::endl;}
};


#endif //SUGATAMA_BEAD_H
