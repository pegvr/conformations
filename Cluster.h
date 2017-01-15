#ifndef CLUSTER_H
#define	CLUSTER_H
#include <vector>
#include <cstdlib>
#include <string>
#include "Hamming.h"
#include "CosineSim.h"
#include "DistanceMatrix.h"
#include "Euclidean.h"



class Cluster {
public:
    Cluster(string point);
    Cluster(const Cluster& orig);
    virtual string getCentroid() {return centroid;} ;
    virtual int * getCentroidD() {} ;
    void InsertPoint(int i, float min_dist){ conformations.push_back(i); distance.push_back(min_dist);};
    virtual void InsertPointE(Euclidean *point, int dist){};
    virtual void InsertPointC(CosineSim *point, int dist){};
    virtual void InsertPointD(int *point, int dist){};
    virtual void InsertPointH(Hamming *point, int dist){};
    void SecondCentroid(int p, float dist){secondcentroid.push_back(p); secondcentroiddist.push_back(dist);};
    int GetPoint(int i){return conformations[i];};
    float GetDistance(int i){return distance[i];};
    virtual Hamming * GetPointH(int i){};
    virtual Euclidean * GetPointE(int i){};
    virtual CosineSim * GetPointC(int i){};
    virtual int * GetPointD(int i){};
    void UpdateCentroid(string ce){centroid = ce;};
    int GetSize(){return conformations.size();};
    int GetSecondCentroidDistance(int i){return secondcentroiddist[i];};
    int GetSecondCentroid(int i){return secondcentroid[i];};
    virtual void PrintCluster(ofstream &outputFile);
    virtual void PopPoint(int i){};
    void PopSecondCentroid(int i);
    void UpdateDistanceM(int *matrix);
    virtual void UpdateDistanceCentroidH(int * matrix, int size){};
    virtual ~Cluster();
private:   
    string centroid;
    vector <int> secondcentroid;
    vector <float> secondcentroiddist;
    vector <int> conformations;
    vector <float> distance;
};

/*class ClusterH : public Cluster{
public:
    ClusterH(Hamming *hamming);
    ClusterH(const ClusterH& orig);
    string getCentroid() {return centroid;} ;
    void InsertPointH(Hamming *point, int dist);
    void PrintCluster();
    Hamming * GetPointH(int i){return points[i];};
    void PopPoint(int i);
    void UpdateDistanceCentroidH(int * matrix, int size);
    virtual ~ClusterH(){        
        for (auto it = points.begin(); it != points.end(); ++it){
        delete *it;
        }
    };
private:
    string centroid;
    vector <Hamming *> points;
    vector <int> distance;
    
};

class ClusterE : public Cluster{
public:
    ClusterE(Euclidean *euclidean);
    ClusterE(const ClusterE& orig);
    string getCentroid() {return centroid;} ;
    void InsertPointE(Euclidean *point, int dist);
    void PrintCluster();
    Euclidean * GetPointE(int i){return points[i];};
    void PopPoint(int i);
    virtual ~ClusterE(){ 
        for (auto it = points.begin(); it != points.end(); ++it){
        delete *it;
    }};
private:
    string centroid;
    vector <Euclidean*> points;
    vector <int> distance;
};

class ClusterC : public Cluster{
public:
    ClusterC(CosineSim *cosine);
    ClusterC(const ClusterC& orig);
    string getCentroid() {return centroid;} ;
    void InsertPointC(CosineSim *point, int dist);
    void PrintCluster();
    CosineSim * GetPointC(int i){return points[i];};
    void PopPoint(int i);
    virtual ~ClusterC(){ 
        for (auto it = points.begin(); it != points.end(); ++it){
        delete *it;
    }};
private:
    string centroid;
    vector <CosineSim*> points;
    vector <int> distance;
};

class ClusterD : public Cluster{
public:
    ClusterD(int *row);
    ClusterD(const Cluster& orig);
    int * getCentroidD() {return centroid;} ;
    void InsertPointD(int *point, int dist);
    void PrintCluster();
    int * GetPointD(int i){return points[i];};
    void PopPoint(int *row);
    virtual ~ClusterD(){ 
        for (auto it = points.begin(); it != points.end(); ++it){
        delete *it;
    }};
private:
    int * centroid;
    vector <int *> points;
};*/

#endif	/* CLUSTER_H */

