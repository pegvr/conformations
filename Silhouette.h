#ifndef SILHOUETTE_H
#define	SILHOUETTE_H
#include "Cluster.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "Cluster.h"
#include "Distances.h"

void Silhouette(Cluster **cluster, int k, int N, string method, float *s);

#endif	/* SILHOUETTE_H */

