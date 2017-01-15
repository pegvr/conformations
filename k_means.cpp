#include "k_means.h"
#include "Conformations.h"
#include <time.h>

void k_means(Cluster **cluster, int numConform, int N, int k)
{
    for(int i = 0; i < k; i++)
    {
        int rand_num = (rand() / (RAND_MAX + 1.0))* numConform; 
        cluster[i] = new Cluster(to_string(rand_num));
    }
}

void k_meansb(Cluster **cluster, int numConform, int N, int k, Conformations conformations[])
{
    for(int i = 0; i < k; i++)
    {
        int rand_num = (rand() / (RAND_MAX + 1.0))* numConform; 
        int rand_num2 = (rand() / (RAND_MAX + 1.0))* N;
        cluster[i] = new Cluster(to_string(conformations[rand_num].GetDistance(rand_num2)));
    }
}
