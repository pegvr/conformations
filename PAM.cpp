#include "PAM.h"
#include "Cluster.h"
#include "cRMSD.h"

 
void PAM(Cluster **cluster, int k, int N, int numConform, string method,double points[][3], float * ObjectiveFunctionF)
{
    int i, j, min_k;
    float distf, min_distf;
    for (i = 1; i < numConform+1; i++)             //For every point
    {

        min_distf = 1000000000000.0;
        for (j = 0; j < k; j++)         //For every cluster compute nearest centroid from point
        {         
            distf = CRMSD( N,  numConform,  points, cluster[j]->getCentroid(), i);
            if (distf < min_distf && distf != 0)
            {
                min_distf = distf;
                min_k = j;
            }  
        }
        ObjectiveFunctionF[min_k] += distf;
        
        cluster[min_k]->InsertPoint(i, min_distf); //Insert point into nearest cluster
        min_distf = 10000;
        int sc;
        for (int p = 0; p < k; p++)     //Find second best cluster
        {
            if (p != min_k)
            {
                distf = CRMSD(N,numConform,  points, cluster[p]->getCentroid(), i);
                if (distf < min_distf && distf != 0)
                {
                    min_distf = distf;
                    sc = p;
                }  
            }
        }
        cluster[min_k]->SecondCentroid(sc, min_distf);
    }

}


void PAM_Update(Cluster **cluster, int k, int N,int numConform, string method, double points[][3], float * ObjectiveFunctionF)
{
    int i, j, min_k = 0, flag = 0, position = 0;
    float distf, min_distf;
    for (i = 0; i < numConform; i++)   //For every point
    {
        min_distf = 100000.0;
        for (j = 0; j < k; j++)         //For every cluster
        {
            distf = CRMSD(N,numConform,  points, cluster[j]->getCentroid(), i);
            for (int t = 0; t < cluster[j]->GetSize(); t++)  //Check where point already belongs
            {
                if (i == cluster[j]->GetPoint(t))
                {
                    flag = j;           //Cluster where point belongs
                    position = t;           //Position in cluster where point belongs
                    break;
                }                
            }
            if (distf < min_distf && distf != 0)       
            {
                min_distf = distf;
                min_k = j;
            }           
        }
        ObjectiveFunctionF[min_k] += distf;
        /*Insert now point into the right cluster, after the update of the centroid*/
        if (min_k != flag && flag < k)                         //If point needs to be assigned again to another cluster
        { 
            cluster[flag]->PopPoint(position);     //Pop point from cluster
            cluster[min_k]->InsertPoint(i, min_distf);
            min_distf = 10000.0;
            int sc;
            for (int p = 0; p < k; p++)
            {
                if (p != min_k)
                {
                    distf = CRMSD(N,numConform,  points, cluster[p]->getCentroid(), i);
                    if (distf < min_distf && distf != 0)
                    {
                        min_distf = distf;
                        sc = p;
                    }  
                }
                cluster[min_k]->SecondCentroid(sc, min_distf);
            }
           
        }
    }
}


void PAMB(Cluster **cluster, int k, int N, int numConform, string method, Conformations conformations[], float * ObjectiveFunctionF)
{
    int i, j, min_k;
    float distf, min_distf;
    for(int c = 0; c < numConform; c++)
    {
        for (i = 0; i < N; i++)             //For every point
        {
            min_distf = 1000000000000.0;
            for (j = 0; j < k; j++)         //For every cluster compute nearest centroid from point
            {
               //cout << "pam 2   " <<  conformations[c].GetPair(i)<< endl;
                float tmp = stof(cluster[j]->getCentroid());
                distf = DistanceEuclidean(conformations[c].GetDistance(i), tmp);
                if (distf < min_distf && distf != 0)
                {
                    min_distf = distf;
                    min_k = j;
                }  
            }
            ObjectiveFunctionF[min_k] += distf;

            int tmp = stoi(conformations[c].GetPair(i), nullptr, 10);
            cluster[min_k]->InsertPoint(tmp, min_distf); //Insert point into nearest cluster
            min_distf = 10000;
            int sc;
            for (int p = 0; p < k; p++)     //Find second best cluster
            {
                if (p != min_k)
                {
                    float tmp = stof(cluster[p]->getCentroid());
                    distf = DistanceEuclidean(conformations[c].GetDistance(i), tmp);
                    if (distf < min_distf && distf != 0)
                    {
                        min_distf = distf;
                        sc = p;
                    }  
                }
            }
                cluster[min_k]->SecondCentroid(sc, min_distf);
        }
    }
}


void PAM_UpdateB(Cluster **cluster, int k, int N,int numConform, string method, Conformations conformations[], float * ObjectiveFunctionF)
{
    int i, j, min_k = 0, flag = 0, position = 0;
    float distf, min_distf;
    for(int c = 0; c < numConform; c++)
    {
        for (i = 0; i < N; i++)   //For every point
        {
            min_distf = 100000.0;
            for (j = 0; j < k; j++)         //For every cluster
            {
                float tmp = stof(cluster[j]->getCentroid());
                distf = DistanceEuclidean(conformations[c].GetDistance(i), tmp);
                for (int t = 0; t < cluster[j]->GetSize(); t++)  //Check where point already belongs
                {
                    if (i == cluster[j]->GetPoint(t))
                    {
                        flag = j;           //Cluster where point belongs
                        position = t;           //Position in cluster where point belongs
                        break;
                    }                
                }
                if  (distf < min_distf && distf != 0)     
                {
                    min_distf = distf;
                    min_k = j;
                }

            }

            ObjectiveFunctionF[min_k] += distf;
            /*Insert now point into the right cluster, after the update of the centroid*/
            if (min_k != flag && flag < k)                         //If point needs to be assigned again to another cluster
            { 
                cluster[flag]->PopPoint(position);     //Pop point from cluster
                cluster[min_k]->InsertPoint(i, min_distf);
                    min_distf = 10000.0;
                    int sc;
                    for (int p = 0; p < k; p++)
                    {
                        if (p != min_k)
                        {
                            float tmp = stof(cluster[p]->getCentroid());
                            distf = DistanceEuclidean(conformations[c].GetDistance(i), tmp);
                            if (distf < min_distf && distf != 0)
                            {
                                min_distf = distf;
                                sc = p;
                            }  
                        }
                    }
                    cluster[min_k]->SecondCentroid(sc, min_distf);

            }
        }
    }
}
