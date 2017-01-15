#ifndef COSINESIM_H
#define COSINESIM_H
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "randomfunc.h"
#include "time.h"


using namespace std;

class CosineSim {
public:
    CosineSim(string temp,string temp1);
    CosineSim(const CosineSim& orig);
    virtual ~CosineSim();
    string ConstructGFunctionC(int k);
    string getId();
    string getName();
private:
    string id;
    int length;
    int pos;
    string name;
};

#endif /* COSINESIM_H */


