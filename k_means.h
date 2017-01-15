#ifndef K_MEANS_H
#define	K_MEANS_H
#include <cstdlib>
#include "Cluster.h"
#include "Conformations.h"

void k_means(Cluster **cluster, int numConform, int N, int k);

void k_meansb(Cluster **cluster, int numConform, int N, int k, Conformations conformations[]);

#endif	/* K_MEANS_H */

