#include "Euclidean.h"
#include "randomfunc.h"



Euclidean::Euclidean(string temp, string temp1) 
{
    name = temp1;       //Name = itemY
    id = temp;          //id = distances from other items
    length = id.length();
}

Euclidean::Euclidean(const Euclidean& orig) {
}

Euclidean::~Euclidean() 
{
    //cout << "Euclidean is destroyed";
}

string Euclidean :: getId()
{
    return id;
}

string Euclidean :: getName()
{
    return name;
}

int Euclidean::ConstructFiFunctionC(int L, int k)
{

    int sum = 0 , p = 0, r, w = 4;
    double res , v;
    string st;
    float t;
    long long M = pow(2.0,32) - 5;
    for( int i = 0; i < length; i++)            //Count how many distances from other items the point has
    {
        st = id[i];
        if(st == "\t")
        {
            p++;
        }
    }

    double array[p + 1];
    istringstream iss(id);
    for (auto& i : array)                   //Array to store each dimension of itemY
    {
        iss>> i;
    }
    srand(time(NULL));    

    int h[k];   
        
    for(int i = 0; i < k; i++)              
    {
        int n = 0; 
        t = float (rand() / (RAND_MAX + 1.0)) * ( w + 1);
        v = marsagliarandom(i);
        
        for (int j = 0;j < (p + 1); j++)
        {
            res += array[j]*v;        //inner product           
        }
        h[i] = int((res + t) / w);
    }
    for (int i = 0; i < k;i++)
    {   
        r = (rand() / (RAND_MAX + 1.0)) * ( (pow(2.0,31)-1) + 1);       
        sum = mod((h[i] * r), M) + sum;        //(a+b)modm = ((amodm)+(bmodm))modm
    }
    
    euclideanid =mod(sum,M);

    return mod(euclideanid, L);     //mod
}

