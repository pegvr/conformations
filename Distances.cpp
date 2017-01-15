#include "Distances.h"


int DistanceHamming(string point1, string point2)
{
    int dist = 0;
    for (int z = 0; z < point1.size() ; z++)          //Find how many bits are different
        if (point1[z] != point2[z]) dist++; 
    return dist;
}


int DistanceCosine(string point1, string point2)
{
    int num = 0, norma1 = 0, norma2 = 0, inner_product = 0, dist, i;
    string  h;
    for(int t = 0; t < point1.length(); t++)                     //Count how many dimensions the point has
    {
        h = point1[t];
        if(h == "\t")   num++;
    }
    double array1[num + 1];                         //Array to store each dimension of bucket point
    istringstream iss(point1);
    for (auto& t : array1)                          //For every dimension, store it
    {
        iss>> t;
    }
    for(int t = 0; t < num + 1; t++)
        norma1 += array1[t] * array1[t];           //norma of bucket point
    norma1 = sqrt(norma1);
    double array2[num + 1];                         //Array to store each dimension of query point
    istringstream iss2(point2);
    for (auto& i : array2)                          //For every dimension, store it
    {
        iss2>> i;
    }
    for(i = 0; i < num + 1; i++)
        norma2 += array2[i] * array2[i];        //norma of query point
    norma2 = sqrt(norma2); 
    for(i = 0; i < num + 1; i++)
        inner_product += array1[i] * array2[i];          //inner product of query and bucket point
    dist = 1 - (inner_product / (norma1 * norma2));    
    return dist;
}


float DistanceEuclidean(float point1, float point2)
{
    int sum = point1 - point2;
    sum = sum * sum;
    return sqrt(sum);
}

float DistanceEuclideanA(string point1, string point2)
{
    int num = 0, sum = 0;
    float distf;
    string h;
    for(int t = 0; t < point1.length(); t++)                     //Count how many dimensions the point has
    {
        h = point1[t];
        if(h == "\t")   num++;
    }

    double array1[num + 1];                         //Array to store each distance of bucket point
    istringstream iss(point1);

    for (auto& t : array1)                          //For every distance, store it
    {
        iss >> t;
    }

    double array2[num + 1];                         //Array to store each distance of query point
    istringstream iss2(point2);

    for (auto& t : array2)                          //For every distance, store it
    {
        iss2 >> t;
    }
    for (int i = 0; i < num + 1 ; i++)              //Calculate "almost" euclidean distance 
        sum += (array1[i] - array2[i]) * ((array1[i] - array2[i]));
    return sqrt(sum);
}


/*int MinDistance(string temp1, string temp2,string temp3, int *row , string method)
{
    int dist1,dist2, min_dist; 
    float distf1, distf2, min_distf;

    if (method == "@metric_space hamming")
    {
        dist1 = DistanceHamming(temp1, temp2);
        dist2 = DistanceHamming(temp1, temp3);                                            
    }                
    else if (method == "@metric euclidean")
    {
        distf1 = DistanceEuclidean(temp1, temp2);
        distf2 = DistanceEuclidean(temp1, temp3);

    }
    else if (method == "@metric cosine")
    {

        dist1 = DistanceEuclidean(temp1, temp2);
        dist2 = DistanceEuclidean(temp1, temp3);
    }
    else
    {   
        dist1 = row[0];
    }  
 
    if ((dist1 < dist2 && dist1 != 0) || (distf1 < distf2 && distf1 != 0))
        return 1;
    else
        return 2;
}*/