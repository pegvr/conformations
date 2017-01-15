#include "randomfunc.h"


double marsagliarandom(int j) // pseudo-random function for normal uniform
{
    
    static double x1 , x2 , r, y1, y2, z;
    int t;
    
    do
    {
        x1 = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        x2 = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        //std:: cout<< "\n\nto x1 einai :" << x1;
        //std:: cout << "\n\n to x2 einai :" << x2;
        r = x1*x1 + x2*x2;
    }
    while ((r>=1.0)|| (r==0.0));
    y1 = sqrt(-2.0 *( log(r) / (r)));
    y1 = y1 * x1;
    y2 = sqrt(-2.0 *log(r) / (r));
    y2 = y2 *x2;
    for (int i = 0; i < 100; i++) // choose between y1 and y2
    {t = mod(j,2);}
    if (t==0)
    {
        return y1;
    }
    return y2;
}


int mod (int a, int b)
{
    if(b < 0)
        return mod(a, -b);   
    int ret = a % b;
    if(ret < 0)
        ret+=b;
    return ret;
}

/*long long int combine(long long int a, int b) {
   int times = 1;
   while (times <= b)
      times *= 10;
   return a*times + b;
} */