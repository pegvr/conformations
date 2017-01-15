#include "alaLloyds.h"

void alalloyds(Cluster **cluster, int k, int N, int numConform, string method, double points[][3], float * ObjectiveFunctionF)
{
    int i, j, t, z, dist , medoid;
    float OFF, min_OFF,distf;
    
    for ( i = 0; i < k; i++)        //For every cluster
    {
        min_OFF = ObjectiveFunctionF[i];
        for (t = 0; t < cluster[i]->GetSize(); t++)         //Find if t is a medoid
        {
            float newdistancesf[cluster[i]->GetSize()];

            OFF = 0;
            for (j = 0; j < cluster[i]->GetSize(); j++)     //Compute distance from possible medoid t to j points in same cluster
            {
                string tmp = to_string(cluster[i]->GetPoint(j));
                distf = CRMSD(N, numConform,  points, tmp,  cluster[i]->GetPoint(t));
                newdistancesf[j] = distf; 

                OFF += distf;
            }
            if (OFF < min_OFF) //If objective function of point is better then we have a possible medoid
            {
                min_OFF = OFF;
                medoid = t;   
            }
     
        }
        /*If objective fynction of medoid is better than the one we already have, update centroid*/
        if (min_OFF < ObjectiveFunctionF[i])
        {
            cluster[i]->UpdateCentroid(to_string(cluster[i]->GetPoint(medoid)));
            ObjectiveFunctionF[i] = min_OFF;
            //cout << "8" << endl;
            PAM_Update(cluster, k, N, numConform, method, points, ObjectiveFunctionF);
            //cout << "after pam" << endl;
        }   
    }
}

void alalloydsb(Cluster **cluster, int k, int N, int numConform, string method, Conformations conformations[], float * ObjectiveFunctionF)
{
    int i, j, t, z, dist, medoid;
    float OFF, min_OFF,distf;
    for ( i = 0; i < k; i++)        //For every cluster
    {
        min_OFF = ObjectiveFunctionF[i];

        for (t = 0; t < cluster[i]->GetSize(); t++)         //Find if t is a medoid
        {
            float newdistancesf[cluster[i]->GetSize()];
            OFF = 0;
            for (j = 0; j < cluster[i]->GetSize(); j++)     //Compute distance from possible medoid t to j points in same cluster
            {               
                float tmp = stof(cluster[i]->getCentroid());
                distf = DistanceEuclidean(cluster[i]->GetDistance(j), tmp);
                newdistancesf[j] = distf; 

                OFF += distf;
            }
            if (OFF < min_OFF) //If objective function of point is better then we have a possible medoid
            {
                min_OFF = OFF;
                medoid = t;   
            }
     
        }
        /*If objective fynction of medoid is better than the one we already have, update centroid*/
        if (min_OFF < ObjectiveFunctionF[i])
        {           
            cluster[i]->UpdateCentroid(to_string(cluster[i]->GetPoint(medoid)));
            ObjectiveFunctionF[i] = min_OFF;
            //cout << "8" << endl;
            PAM_UpdateB(cluster, k, N, numConform, method, conformations, ObjectiveFunctionF);
            //cout << "after pam" << endl;
        }   
    }
}
