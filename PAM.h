#ifndef PAM_H
#define	PAM_H
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "Cluster.h"
#include "Distances.h"
#include "cRMSD.h"
#include "Conformations.h"

using namespace std;

void PAM(Cluster **cluster, int k, int N, int numConform, string method,double points[][3], float * ObjectiveFunctionF);
void PAM_Update(Cluster **cluster, int k, int N, int numConform, string method,double points[][3], float * ObjectiveFunctionF);

void PAMB(Cluster **cluster, int k, int N, int numConform, string method, Conformations conformations[], float * ObjectiveFunctionF);
void PAM_UpdateB(Cluster **cluster, int k, int N, int numConform, string method, Conformations conformations[], float * ObjectiveFunctionF);

#endif	/* PAM_H */

