#include "RasRandomNumber.h"

#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds
#include <Eigen/Cholesky>
using namespace Eigen;


unsigned RasRandomNumber::ras_now_nanoseconds(void)
{
    std::chrono::nanoseconds ns = std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::system_clock::now().time_since_epoch());
    return unsigned(ns.count() % 100000000)+std::rand();
}


std::vector<std::vector<double> > RasRandomNumber::ras_mvnorm(int n, std::vector<double> mu, std::vector<std::vector<double> > corr)
{
    unsigned seed=ras_now_nanoseconds();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);
    
    // chol is cholesky decomposition, should be square
    int ncol=(int)corr.size();
    
    MatrixXd A(ncol,ncol);
    for (int i=0; i<ncol; i++)
        for (int j=0; j<ncol; j++)
            A(i,j)=corr[i][j];
        
    MatrixXd U = A.llt().matrixU();
    std::vector<std::vector<double> > chol(ncol,std::vector<double> (ncol,0) );
    for (int i=0; i<ncol; i++)
        for (int j=0; j<ncol; j++)
            chol[i][j]=U(i,j);
    
    
    // initialization
    std::vector<std::vector<double> > ret(n,std::vector<double> (ncol,0));
    std::vector<std::vector<double> > z  (n,std::vector<double> (ncol,0));
    
    
    for (int i=0; i<n; i++){
        for (int j=0; j<ncol; j++){
            double number = distribution(generator);
            z[i][j]=number;
        }
    }
    std::vector<std::vector<double> > C= RasMatrix::ras_prod_mat(z,chol);
    for (int i=0; i<n; i++){
        for (int j=0; j<ncol; j++){
            ret[i][j]=mu[j]+C[i][j];
        }
    }
    return ret;
}



std::vector<int> RasRandomNumber::ras_rpois(int n, double lam)
{
//    std::default_random_engine generator(time(0));
    unsigned seed=ras_now_nanoseconds();
    std::default_random_engine generator(seed);
    std::poisson_distribution<int> distribution(lam);
    
    std::vector<int> number(n);
    for (int i=0; i<n; i++)
        number[i] = distribution(generator);
    return number;
}

int RasRandomNumber::ras_rpois(double lam)
{
    //    std::default_random_engine generator(time(0));
    unsigned seed=ras_now_nanoseconds();
    std::default_random_engine generator(seed);
    std::poisson_distribution<int> distribution(lam);
    
    int number = distribution(generator);
    return number;
}



double RasRandomNumber::ras_uniform()
{
    static std::default_random_engine re;
    static std::uniform_real_distribution<double> Dist(0,1);
    return Dist(re);
}

std::vector<int> RasRandomNumber::ras_SampleWithoutReplacement(int populationSize, int sampleSize)
{
    // Use Knuth's variable names
    std::vector<int> samples(sampleSize);
    int n = sampleSize;
    int N = populationSize;
    
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;
    
    while (m < n)
    {
        u = ras_uniform(); // call a uniform(0,1) random number generator
        
        if ( (N - t)*u >= n - m )
        {
            t++;
        }
        else
        {
            samples[m] = t;
            t++; m++;
        }
    }
    return samples;
}


