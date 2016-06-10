#include "RasRandomNumber.h"

#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds
#include <Eigen/Cholesky>
using namespace Eigen;


unsigned RasRandomNumber::ras_now_nanoseconds(void)
{
    std::chrono::nanoseconds ns = std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::system_clock::now().time_since_epoch());
    return unsigned(ns.count() % 100000000)+std::rand();
}


std::vector<std::vector<double> > RasRandomNumber::ras_mvnorm(unsigned long int n, std::vector<double> &mu, std::vector<std::vector<double> > &corr, unsigned seed)
{
    if (seed==0) seed=ras_now_nanoseconds();
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
    
    
    for (unsigned long int i=0; i<n; i++){
        for (int j=0; j<ncol; j++){
            z[i][j] = distribution(generator); // standard normal rv
        }
    }
    std::vector<std::vector<double> > C= RasMatrix::ras_prod_mat(z,chol);
    for (unsigned long int i=0; i<n; i++){
        for (int j=0; j<ncol; j++){
            ret[i][j]=mu[j]+C[i][j];
        }
    }
    return ret;
}



std::vector<int> RasRandomNumber::ras_rpois(unsigned long int n, double lam, unsigned seed)
{
    if (seed==0) seed=ras_now_nanoseconds();
    std::default_random_engine generator(seed);
    std::poisson_distribution<int> distribution(lam);
    
    std::vector<int> number(n);
    for (unsigned long int i=0; i<n; i++)
        number[i] = distribution(generator);
    return number;
}

int RasRandomNumber::ras_rpois(double lam, unsigned seed)
{
    //    std::default_random_engine generator(time(0));
    if (seed==0) seed=ras_now_nanoseconds();
    std::default_random_engine generator(seed);
    std::poisson_distribution<int> distribution(lam);
    
    int number = distribution(generator);
    return number;
}



double RasRandomNumber::ras_uniform(unsigned seed)
{
    if (seed==0) seed=ras_now_nanoseconds();
    static std::default_random_engine re(seed);
    static std::uniform_real_distribution<double> Dist(0,1);
    return Dist(re);
}

std::vector<unsigned long int> RasRandomNumber::ras_SampleWithoutReplacement(unsigned long int populationSize, unsigned long int sampleSize, unsigned seed)
{
    if (seed==0) seed=ras_now_nanoseconds();
    static std::default_random_engine re(seed);
    static std::uniform_real_distribution<double> Dist(0,1);

    // Use Knuth's variable names
    std::vector<unsigned long int> samples(sampleSize);
    int n = sampleSize;
    int N = populationSize;
    
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;
    
    while (m < n)
    {
        u = Dist(re); // call a uniform(0,1) random number generator
        
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


