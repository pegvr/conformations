#ifndef CONFORMATIONS_H
#define	CONFORMATIONS_H
#include <iostream>
#include <vector>

using namespace std;

class Conformations {
public:
    Conformations();
    Conformations(const Conformations& orig);
    void InsertPoint(string tmp, float dist);
    void PrintConformations();
    float GetDistance(int i){return distances[i];};
    string GetPair(int i){return pairs[i];};
    virtual ~Conformations();
private:
    vector <string> pairs;
    vector <float> distances;
};

#endif	/* CONFORMATIONS_H */

