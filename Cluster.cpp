#include "Cluster.h"


Cluster::Cluster(string point) 
{
    centroid = point;
}

Cluster::Cluster(const Cluster& orig) 
{
    return;
}

Cluster::~Cluster() 
{
    
}

/*ClusterH :: ClusterH(Hamming *hamming)
{
    centroid = hamming->getId();
    
}

ClusterE :: ClusterE(Euclidean *euclidean)
{
    centroid = euclidean->getId();
}

ClusterC :: ClusterC(CosineSim *cosine)
{
    centroid = cosine->getId();
}

ClusterD::ClusterD(int* row)
{
    centroid = row;
}


void ClusterH :: InsertPointH(Hamming *point, int dist)
{
    points.push_back(point);
    distance.push_back(dist);
}

void ClusterE :: InsertPointE(Euclidean *point, int dist)
{
    points.push_back(point);
    distance.push_back(dist);
}

void ClusterC :: InsertPointC(CosineSim *point, int dist)
{
    points.push_back(point);
    distance.push_back(dist);
}

void ClusterD::InsertPointD(int* point, int dist)
{
   points.push_back(point);
}*/

void Cluster:: PrintCluster(ofstream &outputFile)
{
    for(int i = 0; i < conformations.size(); i++)
    {
        outputFile << conformations[i] << "\t";
        outputFile.flush();
    }
}

/*void ClusterH:: PrintCluster()
{
    //cout << "size =  " << points.size() << endl;
    for(int i = 0; i < points.size(); i++)
    {
        cout << "point  " << i << " =    " << points[i]->getId() << endl;
    }
}

void ClusterC:: PrintCluster()
{
    for(int i = 0; i < points.size(); i++)
    {
        cout << "point  " << i << " =    " << points[i]->getId() << endl;
    }
}

void ClusterD:: PrintCluster()
{
    for(int i = 0; i < points.size(); i++)
    {
        cout << "point  " << i << " =    " << points[i] << endl;
    }
}

void ClusterE:: PrintCluster()
{
    for(int i = 0; i < points.size(); i++)
    {
        cout << "point  " << i << " =    " << points[i]->getId() << endl;
    }
}*/

void Cluster::UpdateDistanceM(int* matrix)
{
    int j = 0;
    for (int i = 0; i < secondcentroiddist.size(); i++)
    {
        if (matrix[i] != 0)
        {
            secondcentroiddist[j] = matrix[i];
            j++;
        }
        
    }
}

/*void ClusterH::UpdateDistanceCentroidH(int* matrix, int size)
{
    int j = 0;
    for (int i = 0; i < size; i++)
    {
        if (matrix[i] != 0 && j < distance.size())
        {
            distance[j] = matrix[i];
            j++;
        }
        
    }
}*/

void Cluster::PopSecondCentroid(int i)
{
    secondcentroid.erase(secondcentroid.begin()+(i-1));
    secondcentroiddist.erase(secondcentroiddist.begin()+(i-1));
}

/*void ClusterH::PopPoint(int i)
{
    points.erase(points.begin()+i);
    distance.erase(distance.begin()+i);
}

void ClusterE::PopPoint(int i)
{
    points.erase(points.begin()+i);
    distance.erase(distance.begin()+i);
    
}

void ClusterC::PopPoint(int i)
{
    points.erase(points.begin()+i);
    distance.erase(distance.begin()+i);
}
*/
