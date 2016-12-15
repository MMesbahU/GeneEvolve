#include "format_hap.h"
#include "CommFunc.h"



bool format_hap::write_hap(Hap_SNP &hap_snp, std::string file_out_name)
{
    // create .hap files
    // create a file for each chromosome
    std::string f_out=file_out_name + ".hap";
    unsigned long int nhap=hap_snp.hap.size();
    unsigned long int nsnp=hap_snp.hap[0].size();
    
    std::ofstream outfile;
    outfile.open(f_out.c_str());
    
    for (unsigned long int k=0; k<nsnp; k++) // for snps
    {
        if (k%1000==0) std::cout << "\r      " << k << " of " << nsnp << " wrote ..." << std::flush;
        for (unsigned long int ih=0; ih<nhap; ih++) // for hap
        {
            outfile << (int)hap_snp.hap[ih][k] << " ";
        }
        outfile << std::endl;
    }
    std::cout << std::endl;

    outfile.close();
    return true;
}




// create .indv files
// create a file for each chromosome
bool format_hap::write_indv(std::vector<unsigned long int> &indv_id, std::string file_out_name)
{
    std::string f_out=file_out_name + ".indv";
    unsigned long int nind=indv_id.size();
    
    std::ofstream outfile;
    outfile.open(f_out.c_str());
    
    for (unsigned long int ih=0; ih<nind; ih++)
    {
        outfile << indv_id[ih] << std::endl;
    }
    
    outfile.close();
    
    return true;
}




// no header: nrow=nsnps, ncols=2*nind
// we convert it to nrow=2*nind, ncol=nsnps
// input : std::string f_name, long int nind, long int nsnp, bool show_iterations
// output: hap_snp
bool format_hap::read_hap(Hap_SNP &hap_snp, std::string f_name, unsigned long int nind, unsigned long int nsnp, bool show_iterations)
{
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return false;
    }
    if(!ifile)
    {
        throw("Error: can not open the file ["+file_name+"] to read.");
        return false;
    }
    
    unsigned long int nhap=nind*2;
    
    
    if (show_iterations) std::cout << "      Allocating memory ..." << std::flush;
    hap_snp.hap.resize(nhap, std::vector<bool>(nsnp) );
    if (show_iterations) std::cout << "      done." << std::endl;
    
    
    std::vector<std::string> st1;
    
    std::string line;
    unsigned long int j=0,i=0;
    
    //cout << endl;
    // each line is one SNP
    while (std::getline(ifile, line)){
        if (show_iterations && j%1000==0) std::cout << "\r      " << j << " of " << nsnp << " read ..." << std::flush;
        for(i=0; i<nhap; i++)
        {
            if (line[2*i]=='0')
                hap_snp.hap[i][j]=false;
            else if (line[2*i]=='1')
                hap_snp.hap[i][j]=true;
            else
            {
                std::cout << "Error: undefined character [" << line[2*i] << "] in file ["+ file_name +"], line number:" << j << std::endl;
                return false;
            }
        }
        j++;
    }
    if (show_iterations) std::cout << std::endl;
    
    if (j!=nsnp)
    {
        throw("Error: in file ["+file_name+"].");
        return false;
    }
    
    ifile.clear();
    ifile.close();
    
    return true;
}

// .legend file has header
unsigned long int format_hap::read_legend(Legend &legend, std::string f_name)
{
    std::string sep=" ";
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    if(!ifile) throw("Error: can not open the file ["+file_name+"] to read.");
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+file_name+"] to read." << std::endl;
        return 0;
    }
    
    std::string id, al0,al1;
    unsigned long int pos;
    
    // discard the first line which is header
    ifile >> id >> id >> id >> id;
    
    unsigned long int j=0;
    while (ifile >> id >> pos >> al0 >> al1){
        legend.id.push_back(id);
        legend.pos.push_back(pos);
        legend.al0.push_back(al0);
        legend.al1.push_back(al1);
        j++;
    }
    
    ifile.clear();
    ifile.close();
    
    return j; // number of snps
}


/*
 bool format_hap::read_hap_chr_old(std::vector<std::vector<bool> > &haps0_chr, std::string f_name, long int &nhap, long int &nsnp)
 {
 std::string sep=" ";
 
 std::string file_name=f_name + ".hap";
 std::ifstream ifile(file_name.c_str());
 
 if(!ifile)
 {
 cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
 return false;
 }
 if(!ifile)
 {
 throw("Error: can not open the file ["+file_name+"] to read.");
 return false;
 }
 
 long int n_r=CommFunc::ras_FileLineNumber(file_name);
 long int n_c=CommFunc::ras_FileColNumber(file_name, sep);
 
 nhap=n_c;
 nsnp=n_r;
 
 
 std::cout << "      Allocating memory ..." << std::flush;
 haps0_chr.resize(nhap, std::vector<bool>(nsnp) );
 std::cout << "      done." << std::endl;
 
 
 std::vector<string> st1;
 
 std::string line;
 long int j=0,i=0;
 
 //cout << endl;
 while (std::getline(ifile, line)){
 if (j%1000==0) cout << "\r      " << j << " of " << n_r << " read ..." << std::flush;
 st1=CommFunc::split(line,sep);
 for(i=0; i<(long int)st1.size(); i++)
 {
 haps0_chr[i][j]= (st1[i]=="1" ? true : false);
 }
 j++;
 }
 cout << endl;
 
 ifile.clear();
 ifile.close();
 
 return true;
 }
 */


