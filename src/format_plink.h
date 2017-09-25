#include <cstdio>
#include <vector>
#include <string>
#include <iostream> // for cout
#include <fstream> // for ofstream

/*
 .ped file has no header with columns: FID IID PID MID sex pheno snp1 snp1 snp2 snp2 ...
 .map file has no header with exactly 4 columns
 
 */



// For the first 6 column of .ped file
class plink_PED_ids
{
public:
    std::vector<std::string> FID;
    std::vector<std::string> IID;
    std::vector<std::string> PID;
    std::vector<std::string> MID;
    std::vector<int> sex;
    std::vector<double> phen;

    plink_PED_ids() { }
    plink_PED_ids(unsigned long int n)
    {
        FID.resize(n);
        IID.resize(n);
        PID.resize(n);
        MID.resize(n);
        sex.resize(n);
        phen.resize(n);
    }
    void alloc(unsigned long int n)
    {
        FID.resize(n);
        IID.resize(n);
        PID.resize(n);
        MID.resize(n);
        sex.resize(n);
        phen.resize(n);
    }
};


// For the 4 column of .map file (chr, rs, cM, pos) and two more columns
class plink_MAP
{
public:
    std::vector<std::string> chr;
    std::vector<std::string> rs;
    std::vector<double> cM;
    std::vector<unsigned long int> pos;
    std::vector<std::string> al0;
    std::vector<std::string> al1;

    plink_MAP() { }
    plink_MAP(unsigned long int n)
    {
        chr.resize(n);
        rs.resize(n);
        cM.resize(n);
        pos.resize(n);
        al0.resize(n);
        al1.resize(n);
    }
    void alloc(unsigned long int n)
    {
        chr.resize(n);
        rs.resize(n);
        cM.resize(n);
        pos.resize(n);
        al0.resize(n);
        al1.resize(n);
    }
};


namespace format_plink
{
    bool write_ped_map(std::string file_out_name, std::vector<std::vector<bool> > &matrix_plink_ped, plink_PED_ids &plink_ped_ids, plink_MAP &plink_map);
    bool write_ped01_map(std::string file_out_name, std::vector<std::vector<bool> > &matrix_plink_ped, plink_PED_ids &plink_ped_ids, plink_MAP &plink_map);
}

