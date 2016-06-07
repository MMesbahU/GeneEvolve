#include "format_plink.h"



bool format_plink::write_ped_map(std::string file_out_name, std::vector<std::vector<bool> > &matrix_plink_ped, plink_PED_ids &plink_ped_ids, plink_MAP &plink_map)
{
    // create .ped and .map files
    // create a file for each chromosome
    std::string sep=" ";
    std::string name_ped=file_out_name + ".ped";
    
    unsigned long int nind=matrix_plink_ped.size();
    if (nind<1)
    {
        std::cout << "Error: no individula to write." << std::endl;
        return false;
    }
    
    unsigned long int nsnp=plink_map.pos.size();
    if (nsnp*2 != matrix_plink_ped[0].size())
    {
        std::cout << "Error: mismatch in SNPs number." << std::endl;
        return false;
    }
    
    
    ////////////////////
    // ped file
    std::ofstream file_ped;
    file_ped.open(name_ped.c_str());
    
    for (unsigned long int i=0; i<nind; i++)
    {
        file_ped << plink_ped_ids.FID[i] << sep;
        file_ped << plink_ped_ids.IID[i] << sep;
        file_ped << plink_ped_ids.PID[i] << sep;
        file_ped << plink_ped_ids.MID[i] << sep;
        file_ped << plink_ped_ids.sex[i] << sep;
        file_ped << plink_ped_ids.phen[i];
        if (i%1000==0) std::cout << "\r      " << i << " of " << nind << " wrote ..." << std::flush;
        for (unsigned long int j=0; j<nsnp; j++)
        {
            // true (1) means al1, false (0) means al0
            std::string s1=matrix_plink_ped[i][2*j] ? plink_map.al1[j] : plink_map.al0[j];
            std::string s2=matrix_plink_ped[i][2*j+1] ? plink_map.al1[j] : plink_map.al0[j];
            file_ped << sep << s1 << sep << s2;
        }
        file_ped << std::endl;
    }
    std::cout << std::endl;
    file_ped.close();
    
 
    ////////////////////
    // map file
    std::string name_map=file_out_name + ".map";
    std::ofstream file_map;
    file_map.open(name_map.c_str());
    for (unsigned long int j=0; j<nsnp; j++)
    {
        file_map << plink_map.chr[j] << sep;
        file_map << plink_map.rs[j]  << sep;
        file_map << plink_map.cM[j]  << sep;
        file_map << plink_map.pos[j] << std::endl;
    }
    file_map.close();

    return true;
}





