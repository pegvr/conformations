#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <ctime>
#include <math.h>
#include "../Eigen/Dense"
#include "../Eigen/Eigen"
#include "../Eigen/SVD"
#include "Conformations.h"
#include "Cluster.h"
#include "PAM.h"
#include "alaLloyds.h"
#include "cRMSD.h"
#include "k_means.h"
#include "Silhouette.h"

using namespace std;
using namespace Eigen;

/*
 * 
 */
int main(int argc, char** argv) 
{
    int N, numConform, i, j;
    string line;
    ifstream inputfile(argv[2]);
    srand(time(NULL));
    if (inputfile.is_open())
    {
        getline (inputfile,line);
        numConform = stoi(line, nullptr, 10);
        getline (inputfile,line);
        N = stoi(line, nullptr, 10);
    }
    double points[N*numConform][3];
    
    if (inputfile.is_open())
    {
        for(i = 0; i < N * numConform; i++)
        {
            getline (inputfile,line);
            size_t pos = line.find("\t");
            string temp1 = line.substr(0, pos);
            string temp3 = line.substr(pos+1, line.size());
            pos = temp3.find("\t");
            string temp2 = temp3.substr(0, pos);
            temp3 = temp3.substr(pos+1, temp3.size());
            points[i][0] = stod(temp1);
            points[i][1] = stod(temp2);
            points[i][2] = stod(temp3);
        }     
        inputfile.close();        
    }
    else    cout << "Unable to open file" << endl;
    int k = 5;
    Cluster **cluster = new Cluster*[k];
    float ObjectiveFunctionF[k];
    
    k_means(cluster, numConform, N, k);
    //for (i = 0; i < k; i++) cout << cluster[i]->getCentroid() << endl;    
    
    for (i = 0; i < k; i++) ObjectiveFunctionF[i] = 0;

    PAM(cluster, k, N,numConform, line, points, ObjectiveFunctionF);
    //for (i = 0; i < k; i++) cout << ObjectiveFunction[i] << endl;
    //for (i = 0; i < k; i++) cout << cluster[i]->GetSize() << endl;
    
    alalloyds(cluster, k, N, numConform, line, points, ObjectiveFunctionF);
    
    ofstream outputFile("conform.txt");
    outputFile << "k: "<< k << endl;
    for (i = 0; i < k; i++)
    {
        cluster[i]->PrintCluster(outputFile); 
        outputFile << "\n" << endl;
    }
    outputFile.close();
    for (i = 0; i < k; i++) delete cluster[i];

    
    /******************************************PART B*****************************************************/
    
    Conformations conformations[numConform];
    int r = (rand() / (RAND_MAX + 1.0))* ((N * (N - 1)) / 2);
    cout << r << endl;
    float x, y , z, x2, y2, z2, dist;
    string tmp, tmp2;
    r = N;
    for( i = 1; i <= numConform; i++)
    {
        for(j = 0; j < r; j++)
        {
            int rand_num1 = (rand() / (RAND_MAX + 1.0))* (N); 
            int rand_num2 = (rand() / (RAND_MAX + 1.0))* (N);
            x = points[(i * rand_num1) - 1][0];
            y = points[(i * rand_num1) - 1][1];
            z = points[(i * rand_num1) - 1][2];
            x2 = points[(i * rand_num2) - 1][0];
            y2 = points[(i * rand_num2) - 1][1];
            z2 = points[(i * rand_num2) - 1][2];
            tmp = to_string(x) + "\t" + to_string(y) + "\t" + to_string(z);
            tmp2 = to_string(x2) + "\t" + to_string(y2) + "\t" + to_string(z2);
            dist = DistanceEuclideanA(tmp, tmp2);
            tmp = to_string(rand_num1) + "\t" + to_string(rand_num2);
            conformations[i - 1].InsertPoint(tmp, dist);
        }
    }
    
    cluster = new Cluster*[numConform];
    k_meansb(cluster, numConform, N, k, conformations);
    for (i = 0; i < k; i++) cout << cluster[i]->getCentroid() << endl;    
        
    for (i = 0; i < k; i++) ObjectiveFunctionF[i] = 0;

    PAMB(cluster, k, N, numConform, line, conformations, ObjectiveFunctionF);
    //for (i = 0; i < k; i++) cout << ObjectiveFunction[i] << endl;
    for (i = 0; i < k; i++) cout << cluster[i]->GetSize() << endl;
    
    alalloydsb(cluster, k, N, numConform, line, conformations, ObjectiveFunctionF);
    
    float s[k], sum = 0;
    Silhouette(cluster,  k,  N,  "", s);
    cout << "Silhouette: [ " ;
    for (i = 0; i < k; i++)
    {
        sum += s[i];
        cout << s[i] << " , ";
    }       
    cout <<  " ] "<< endl;
    cout <<  " Sum: "<< sum << endl;
    
    for (i = 0; i < k; i++) delete cluster[i];
    return 0;
}

