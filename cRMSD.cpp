#include "cRMSD.h"


float CRMSD(int N, int numConform, double points[][3], string centroid, int conf2)
{
    int conf1 = stoi(centroid, nullptr, 10);
    float norm;
    
    MatrixXd X(N, 3);
    MatrixXd Y(N, 3);

    for( int i = 0; i < N; i++)
    {
        X(i,0) = points[N * (conf1 - 1)+ i][0];
        X(i,1) = points[N * (conf1 - 1)+ i][1];
        X(i,2) = points[N * (conf1 - 1)+ i][2];
        Y(i,0) = points[N * (conf2 - 1)+ i][0];
        Y(i,1) = points[N * (conf2 - 1)+ i][1];
        Y(i,2) = points[N * (conf2 - 1)+ i][2];
    }

    double x_c [3] = {0.0 , 0.0, 0.0}, y_c[3] = {0.0 , 0.0, 0.0} ;
    
    for (int i = 0; i < N; i++)
    {
        x_c[0] = x_c[0] + X(i,0);
        x_c[1] = x_c[1] + X(i,1);
        x_c[2] = x_c[2] + X(i,2);
        y_c[0] = y_c[0] + Y(i,0);
        y_c[1] = y_c[1] + Y(i,1);
        y_c[2] = y_c[2] + Y(i,2);
    }

    for( int i = 0; i < N; i++)
    {
        X(i,0) = X(i,0) - x_c[0];
        X(i,1) = X(i,1) - x_c[1];
        X(i,2) = X(i,2) - x_c[2];
        Y(i,0) = Y(i,0) - y_c[0];
        Y(i,1) = Y(i,1) - y_c[1];
        Y(i,2) = Y(i,2) - y_c[2];
    }
    
   
    
    MatrixXd X_T(3, N);
    
    X_T = X.adjoint();
 
    
    MatrixXd mul = X_T * Y;
  
    
    JacobiSVD<MatrixXd> svd( mul, ComputeFullV | ComputeFullU );
    //cout << svd.computeU() << endl;
    //cout << svd.computeV() << endl;

    
    JacobiSVD<MatrixXd>::SingularValuesType singular = svd.singularValues();
    MatrixXd Cp = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().adjoint();
    //cout << Cp << endl;
    
    if (singular(2) > 0)
    {       
        MatrixXd Q = svd.matrixU() * svd.matrixV().adjoint();
        if (Q.determinant() < 0)
        {
            MatrixXd T = svd.matrixU();
            T(0,2) = T(0,2) * (-1);
            T(1,2) = T(1,2) * (-1);
            T(2,2) = T(2,2) * (-1);
            Q = T * svd.matrixV().adjoint();
        }
        
        MatrixXd B = (X * Q) ;
        Q = B - Y; 
        norm = Q.norm()/sqrt(N);
        //cout << "norm:  " << norm << endl; 
        
    }
    return norm;
}
