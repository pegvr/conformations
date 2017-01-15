#include "Conformations.h"

Conformations::Conformations() 
{
    //cout << "here" << endl;
}


Conformations::Conformations(const Conformations& orig) {
    vector <string> pairs(orig.pairs);
    vector <float> distances(orig.distances);
    cout << "haaaa" << endl;
}

Conformations::~Conformations() {
}

void Conformations::InsertPoint(string tmp, float dist)
{
    pairs.push_back(tmp);
    distances.push_back(dist);
}

void Conformations:: PrintConformations()
{
    for(int i = 0; i < pairs.size(); i++)
    {
        cout << "point  " << i << " =    " << pairs[i] << endl;
    }
}