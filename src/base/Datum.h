//
// Created by xos81802 on 30/04/2018.
//

#ifndef SUGATAMA_DATUM_H
#define SUGATAMA_DATUM_H


class Datum {

    float q, iofq, sigma;
    unsigned int index;
    float var, invvar;

public:

    Datum(float q, float iofq, float sigma, unsigned int index) : q(q), iofq(iofq), sigma(sigma), index(index) {
        var = sigma*sigma;
        invvar = 1.0/var;
    }

    // copy constructor
    Datum(const Datum & dat) : q(dat.q), iofq(dat.iofq), sigma(dat.sigma), index(dat.index) {
        this->var = sigma*sigma;
        this->invvar = 1.0d/this->var;
    }

    ~Datum(){

    }

    const float getQ() const {return q;}
    const float getI() const {return iofq;}
    float getSigma(){return sigma;}
    float getVar(){return var;}
    float getInvVar(){return invvar;}

    const unsigned int getIndex() const {
        return index;
    }

};


#endif //SUGATAMA_DATUM_H
