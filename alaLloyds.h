#ifndef ALALLOYDS_H
#define	ALALLOYDS_H
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "Cluster.h"
#include "Distances.h"
#include "PAM.h"

void alalloyds(Cluster **cluster, int k, int N, int numConform, string method, double points[][3],  float * ObjectiveFunctionF);

void alalloydsb(Cluster **cluster, int k, int N, int numConform, string method, Conformations conformations[], float * ObjectiveFunctionF);


#endif	/* ALALLOYDS_H */

