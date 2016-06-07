#include "RasMatrix.h"


std::vector<std::vector<double> > RasMatrix::ras_prod_mat(std::vector<std::vector<double> > &A, std::vector<std::vector<double> > &B)
{
    int r=A.size();
    int m1=A[0].size();
    int m2=B.size();
    int c=B[0].size();
    
    // initial
    std::vector<std::vector<double> > ret(r,std::vector<double> (c,0));
    
    if(m1!=m2) return ret;
    
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            double s=0;
            for(int k=0;k<m1;k++){
                s=s+(A[i][k] * B[k][j]);
            }
            ret[i][j]=s;
        }
    }
    return ret;
}


