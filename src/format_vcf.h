#include <cstdio>
#include <vector>
#include <string>
#include <iostream> // for cout
#include <fstream> // for ofstream
#include <cmath> // for round, floor, ceil, trunc

#include "VcfFile.h"
#include "VcfFileReader.h"
#include "VcfFileWriter.h"


// this is for biallelic SNPs
class vcf_structure
{
    public:
    std::vector<std::string> meta_lines;
    std::vector<std::string> SAMPLES;
    std::vector<std::string> CHROM;
    std::vector<int>         POS;
    std::vector<std::string> ID;
    std::vector<std::string> REF;
    std::vector<std::string> ALT;
    std::vector<float> QUAL;
    std::vector<std::string> FILTER;
    std::vector<std::string> INFO;
    std::vector<std::string> FORMAT;
    std::vector<std::vector<bool> > data; // nhap*nsnp with 0=REF=false, 1=true=ALT. In C, 0 means false
};


namespace format_vcf
{
    bool write_vcf_file(std::string file_out_name, vcf_structure &vcf_st);
    bool read_vcf_file(std::string filename, vcf_structure &vcf_st);
    bool read_vcf_header_sample(std::string filename, std::vector<std::string> &sample);
};
