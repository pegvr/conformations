#ifndef DISTANCES_H
#define	DISTANCES_H
#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>
#include <sstream>

using namespace std;

int DistanceHamming(string point1, string point2);
int DistanceCosine(string point1, string point2);
float DistanceEuclidean(float point1, float point2);
float DistanceEuclideanA(string point1, string point2);
//int MinDistance(string temp1, string temp2,string temp3, int *row , string method);

#endif	/* DISTANCES_H */

