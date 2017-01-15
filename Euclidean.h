#ifndef EUCLIDEAN_H
#define	EUCLIDEAN_H
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>

using namespace std;

class Euclidean {
public:
    Euclidean(string temp, string temp1);
    Euclidean(const Euclidean& orig);
    virtual ~Euclidean();
    int ConstructFiFunctionC(int L, int k);
    string getId();
    string getName();
private:
    string id, name;
    int length;
    int euclideanid;
};

#endif	/* EUCLIDEAN_H */

