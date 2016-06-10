#include <cstdio>
#include <vector>
#include <string>

/*
 HAP
 each row of a hap file is one snp and each column is a haplotype: 0 0 1 1 0 1
 where 0 means al0 and 1 means al1
 // .hap file: nrow=nsnps,  ncols=2*nind
 // Hap_SNP  : nrow=2*nind, ncol=nsnps

 LEGEND
 .legend file has a header. Columns are: ID pos allele0 allele1
 */


// for IMPUTE2 file format
class Legend
{
public:
    std::vector<std::string> id;
    std::vector<unsigned long int> pos;
    std::vector<std::string> al0;
    std::vector<std::string> al1;
};

class Hap_SNP
{
public:
    std::vector<std::vector<bool> > hap; //each row is nhaps=2*ind, ncols=nSNPs
};

namespace format_hap
{
    unsigned long int read_legend(Legend &legend, std::string f_name);
    bool read_hap(Hap_SNP &hap_snp, std::string f_name, unsigned long int nhap, unsigned long int nsnp, bool show_iterations);
    bool write_hap(Hap_SNP &hap_snp, std::string file_out_name);
    bool write_indv(std::vector<unsigned long int> &indv_id, std::string file_out_name);
};


