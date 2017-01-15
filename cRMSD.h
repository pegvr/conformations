#ifndef CRMSD_H
#define	CRMSD_H
#include <iostream>
#include <string>
#include <time.h>
#include <ctime>
#include <math.h>
#include "../Eigen/Dense"
#include "../Eigen/Eigen"
#include "../Eigen/SVD"
#include "Conformations.h"

using namespace std;
using namespace Eigen;

float CRMSD(int N, int numConform, double points[][3], string centroid, int conf2);


#endif	/* CRMSD_H */

