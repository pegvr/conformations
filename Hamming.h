#ifndef HAMMING_H
#define	HAMMING_H

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "randomfunc.h"
#include "time.h"

using namespace std;

class Hamming {
public:
    Hamming(string temp, string temp1);
    Hamming(const Hamming& orig);
    virtual ~Hamming();
    string ConstructGFunction(int k);
    string getId();
    string getName();
private:
    string id, name;
    int length;
};


#endif	/* HAMMING_H */
