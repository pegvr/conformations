#include "Hamming.h"

Hamming::Hamming(string temp, string temp1) 
{
    name = temp1;               //Name = itemY
    id = temp;                  //id = string of '0' and '1'
    length = id.size() - 1;
}

Hamming::Hamming(const Hamming& orig) {
}

Hamming::~Hamming() 
{
    //cout << "hamming destructor" << endl;
}

string Hamming:: getId()
{
    return id;
}

string Hamming:: getName()
{
    return name;
}

string Hamming::ConstructGFunction(int k)
{
    int i; 
    string g;
    //string h;
    string f="goo";
    //long long int a;
    

    /*int array[length];
    istringstream iss(id);
    for (auto& i : array)                   //Array to store each dimension of itemY
    {
        iss>> i;
    }
    
    for (i=0;i<length;i++)
    {
        cout << array[i] << "\n" ; 
    }
    
    srand(time(NULL)); */
    for (i = 0; i < k; i++)     //Pick random k bits of item
    {
        //pick uniformly k random numbers from 0 to 63
        int x1 = (rand() / (RAND_MAX + 1.0))*(id[id.size()-1] + 1);
        g = g + id[x1];
        //a = combine(a,array[x1]);
        //cout << "a: " << a << "\n" ;
        
    }
    //cout << "k: " << k << "\n";
   // cout << a << "\n"; 
    return g;
}



