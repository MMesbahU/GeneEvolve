#include <random>
#include "Population.h"
#include "CommFunc.h"
#include <cmath> // for round, floor, ceil, trunc
#include <exception>
#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds


// this file has header
// file_gen_info is a space delimited file with 3 columns: _pop_size _mat_cor _offspring_dist
// each row is a generation
// output: npop
int Population::ras_read_generation_info_file(std::string file_gen_info)
{
    char sep=' ';
    
    std::string file_name=file_gen_info;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }
    
    if (CommFunc::ras_FileColNumber(file_name," ") !=6 )
    {
        std::cout << "Error: file ["+ file_name +"] must have 6 columns: pop_size, mat_cor, offspring_dist, selection_func, selection_func_par1 and selection_func_par2." << std::endl;
        return 0;
    }
    
    std::vector<std::string> st1;
    
    std::string line;
    int ngen=0;
    
    // discard the first line which is header
    std::getline(ifile, line);

    while (std::getline(ifile, line)){
        std::istringstream iss(line);
        std::string token;
        
        int ps;
        double mc;
        std::string offspring_d;
        std::string selection_func="";
        double selection_func_p1;
        double selection_func_p2;

        std::getline(iss, token, sep); // pop_size
        ps=std::stod(token); // it should be stod for 3e+05
        std::getline(iss, token, sep); //mate_corr
        mc=std::stod(token);
        std::getline(iss, token, sep); //offspring_dist
        offspring_d=token;
        std::getline(iss, token, sep); //selection_func
        selection_func=token;
        std::getline(iss, token, sep); //selection_func_p1
        selection_func_p1=std::stod(token);
        std::getline(iss, token, sep); //selection_func_p2
        selection_func_p2=std::stod(token);
        
        if(mc>1 || mc<-1)
        {
            mc=0;
            std::cout << " Warning in file ["+ file_name +"]: mate_corr should be in range [-1,1]. We set it to 0." << std::endl;
        }
        if (!(offspring_d=="p" || offspring_d=="f"))
        {
            offspring_d="p";
            std::cout << " Warning in file ["+ file_name +"]: offspring_dist should be [p] or [f]. We set it to [p]." << std::endl;
        }
        if (!(selection_func=="logit" || selection_func=="probit" || selection_func=="stab" || selection_func=="thr"))
        {
            selection_func="logit";
            selection_func_p1=0;
            selection_func_p2=1;
            std::cout << " Warning in file ["+ file_name +"]: selection_func should be [logit,probit,stab,thr]. We set it to [logit 0 1]." << std::endl;
        }
        _pop_size.push_back(ps);
        _mat_cor.push_back(mc);
        _offspring_dist.push_back(offspring_d);
        _selection_func.push_back(selection_func);
        _selection_func_par1.push_back(selection_func_p1);
        _selection_func_par2.push_back(selection_func_p2);
        ngen++;
    }
    
    return ngen;
}


// this file has header
// f_name is a space dilimited file with 4 columns: chr file.hap file.legend fiel.sample
// each row is a chromosome
// output: npop
int Population::ras_read_hap_legend_sample_address_name(std::string f_name)
{
    std::string sep=" ";
    
    _hap_legend_sample_name.clear();
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }
    
    std::vector<std::string> st1;
    
    std::string line;
    int i=0;
    
    // discard first line which is heder
    std::getline(ifile, line);
    //cout << endl;
    _all_active_chrs.clear();
    while (std::getline(ifile, line)){
        std::istringstream iss(line);
        int chr;
        iss >> chr;
        _all_active_chrs.push_back(chr);
        std::string c_hap,c_legend,c_sample;
        iss >> c_hap >> c_legend >> c_sample;
        std::vector<std::string> v1(3);
        v1[0]=c_hap;
        v1[1]=c_legend;
        v1[2]=c_sample;
        _hap_legend_sample_name.push_back(v1);
        i++;
    }
    
    return i; // number of lines read
}





// THIS IS FOR Additive and Dominance

// processing the '--file_cv_info f_name'
// just for '_all_active_chrs' chrs, saves the cvs info in 'cv_info'
// f_name has header
// f_name is a space dilimited file with 4 columns: chr pos a d
// each row is a chromosome
// note to read ras_read_hap_legend_sample_address_name before running this func
//output: ncv
int Population::ras_read_cv_info_dominace_model_file(std::string f_name, int iphen)
{
    unsigned long int file_ncol=4;
    _pheno_scheme[iphen]._cv_info.resize(_nchr);
    
    char sep=' ';
    
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std:: cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }
    
    if(CommFunc::ras_FileColNumber(file_name," ") != file_ncol)
    {
       std:: cout << "Error: file ["+ file_name +"] should have " << file_ncol << " columns." << std::endl;
       return 0;
    }
    
    std::string line;
    
    int i=0;
    // this file has header
    // discard the first line which is heder
    std::getline(ifile, line);
    //std::cout << std::endl;
    while (std::getline(ifile, line)){
        std::istringstream iss(line);
        std::string token;
        
        int chr;
        unsigned long int bp;
        double genetic_value_a, genetic_value_d;
        //iss >> chr >> bp >> a >> d;
        std::getline(iss, token, sep);
        chr=std::stoi(token);
        int j=ras_get_ind_active_chr(chr);
        if (j>=0)
        {
            std::getline(iss, token, sep);
            bp=std::stod(token); // it shopuld be stod and not stol, because in the input file we may have 1.2e+4
            std::getline(iss, token, sep);
            genetic_value_a=std::stod(token);
            std::getline(iss, token, sep);
            genetic_value_d=std::stod(token);
            
            _pheno_scheme[iphen]._cv_info[j].bp.push_back(bp);
            _pheno_scheme[iphen]._cv_info[j].genetic_value_a.push_back(genetic_value_a);
            _pheno_scheme[iphen]._cv_info[j].genetic_value_d.push_back(genetic_value_d);
        }
        
        i++;
    }
    
    return i;
}





int Population::ras_get_ind_active_chr(int chr)
{
    for (int j=0; j<(int)_all_active_chrs.size(); j++)
        if (_all_active_chrs[j]==chr) return j;
    
    return -1;
}


// this function reads '--file_cvs file' and just for '_all_active_chrs' chrs, saves the address of .hap files in name_cv_hap
// f_name has NO header
// f_name is a space delimited file with 2 columns: chr cv.hap
// each row is a chromosome
// output: nchr
int Population::ras_read_cvs_address_name(std::string f_name, int iphen)
{
    _pheno_scheme[iphen]._name_cv_hap.resize(_nchr,"");
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }
    
    std::string line;
    
    int i=0;
    //cout << endl;
    while (std::getline(ifile, line)){
        std::istringstream iss(line);
        int chr;
        iss >> chr;
        std::string s1;
        iss >> s1;
        int j=ras_get_ind_active_chr(chr);
        if (j>=0)
            _pheno_scheme[iphen]._name_cv_hap[j]=s1;
        i++;
    }
    
    return i;
}


// this function loads the cvs for all chrs which their address is in 'name_cv_hap'
// and saves them in 'cvs'
bool Population::ras_load_cvs(int iphen)
{
    int ncv_hapfile=_pheno_scheme[iphen]._name_cv_hap.size();
    _pheno_scheme[iphen]._cvs.resize(ncv_hapfile);
    std::string sep=" ";
    unsigned long int nhap=CommFunc::ras_FileColNumber(_pheno_scheme[iphen]._name_cv_hap[0],sep);
    unsigned long int nind=nhap/2;
    for (int ichr=0; ichr<ncv_hapfile; ichr++)
    {
        Hap_SNP hap_snp;
        if(_pheno_scheme[iphen]._name_cv_hap[ichr].size()>0)
        {
            unsigned long int ncv_chr=CommFunc::ras_FileLineNumber(_pheno_scheme[iphen]._name_cv_hap[ichr]);
            if (ncv_chr<1)
            {
                std::cout << "Error reading file [" << _pheno_scheme[iphen]._name_cv_hap[ichr] << "]. No CV." << std::endl;
                return false;
            }
            bool show_iterations=false;
            bool t=format_hap::read_hap(hap_snp, _pheno_scheme[iphen]._name_cv_hap[ichr], nind, ncv_chr, show_iterations);
            if (!t)
            {
                std::cout << "Error reading file [" << _pheno_scheme[iphen]._name_cv_hap[ichr] << "]" << std::endl;
                return false;
            }
            _pheno_scheme[iphen]._cvs[ichr].val=hap_snp.hap;
        }
    }
    return true;
}


// this file has header
// f_name is a space delimited file with 3 columns: chr bp cM
// each row is a snp
unsigned long int Population::ras_read_rmap(std::string f_name)
{
    _rmap.resize(_nchr);
    char sep=' ';
    
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }
    
    std::string line;
    
    unsigned long int i=0;
    // discard first line which is heder
    std::getline(ifile, line);
    
    // we just read rmaps with chr in '_all_active_chrs'
    while (std::getline(ifile, line)){
        std::istringstream iss(line);
        std::string token;
        int chr;
        unsigned long int bp;
        double cM;
        //iss >> chr >> bp >> cM;
        std::getline(iss, token, sep);
        chr=std::stoi(token);
        int j=ras_get_ind_active_chr(chr);
        if (j>=0)
        {
            std::getline(iss, token, sep);
            bp=std::stod(token); // it should be stod and not stol, because in the input file we may have 1.2e+4
            std::getline(iss, token, sep);
            cM=std::stod(token);
            
            //4 1.91e+08 220.003969423109
            //if(chr==4) cout << chr << " qqqqqqq " << bp  << " " << cM<< endl;
            _rmap[j].bp.push_back(bp);
            _rmap[j].cM.push_back(cM);
        }
        i++;
    }
    
    // compute bp_dist_in_rmap
    for (int ichr=0; ichr<_nchr; ichr++)
        _rmap[ichr].bp_dist_in_rmap=_rmap[ichr].bp[1]-_rmap[ichr].bp[0];
    
    
    if (_debug)
    {
        for (int ichr=0; ichr<_nchr; ichr++)
        {
            std::cout << "  rmap bp distance in chr " << _all_active_chrs[ichr] << "=" << _rmap[ichr].bp_dist_in_rmap << std::endl;
            unsigned long int ll=_rmap[ichr].cM.size();
            std::cout << "  rmap: ";
            for (unsigned long int kj=20; kj>0; kj--)
                std::cout << _rmap[ichr].cM[ll-kj] << " ";
            std::cout << std::endl;
        }
    }
    
    return i;
}


// this file has header
// f_name is a space delimited file with 3 columns: chr bp mutation_rate
// each row is a snp
unsigned long int Population::ras_read_file_mutation(std::string f_name)
{
    _mutation_map.resize(_nchr);
    char sep=' ';
    
    std::string file_name=f_name;
    std::ifstream ifile(file_name.c_str());
    
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }
    
    std::string line;
    
    unsigned long int i=0;
    // discard first line which is heder
    std::getline(ifile, line);
    
    // we just read _mutation_map with chr in '_all_active_chrs'
    while (std::getline(ifile, line)){
        std::istringstream iss(line);
        std::string token;
        int chr;
        unsigned long int bp;
        double mutation_rate;
        //iss >> chr >> bp >> mutation_rate;
        std::getline(iss, token, sep);
        chr=std::stoi(token);
        int j=ras_get_ind_active_chr(chr);
        if (j>=0)
        {
            std::getline(iss, token, sep);
            bp=std::stod(token); // it shopuld be stod and not stol, because in the input file we may have 1.2e+4
            std::getline(iss, token, sep);
            mutation_rate=std::stod(token);
            // check its range
            if (mutation_rate<0 || mutation_rate>1) mutation_rate=0;
            
            //4 1.91e+08 220.003969423109
            _mutation_map[j].bp.push_back(bp);
            _mutation_map[j].mutation_rate.push_back(mutation_rate);
        }
        i++;
    }
    
    return i; // number of read lines (zero based)
}


bool Population::ras_compute_recom_prob(void)
{
    _recom_prob.resize(_nchr);
    for (int chr=0; chr<_nchr; chr++)
    {
        _recom_prob[chr].resize(_rmap[chr].cM.size());
        _recom_prob[chr][0]=0;
        for (unsigned long int k=1; k<_rmap[chr].cM.size(); k++)
            _recom_prob[chr][k]=(_rmap[chr].cM[k]-_rmap[chr].cM[k-1])*.01;
    }
    
    
    /*
     This part tests the accuracy of this function
     for (int chr=0; chr<_nchr; chr++)
     {
     std::ofstream fo;
     fo.open((_out_prefix+".computed_recom_prob.chr"+to_string(_all_active_chrs[chr])+".txt").c_str());
     for (int k=0; k<(int)_rmap[chr].cM.size(); k++)
     {
     fo << _recom_prob[chr][k] << std::endl;
     }
     fo.close();
     }
     */
    
    if (_debug)
    {
        for (int chr=0; chr<_nchr; chr++)
        {
            //for (int k=0; k<60; k++)
            //    cout << recom_prob[chr][k] << " ";
            std::cout << "  mean(recom_prob)=" << CommFunc::mean(_recom_prob[chr]) << ", recom_prob[end]=" << _recom_prob[chr][_recom_prob[chr].size()-1] << std::endl;
        }
    }
    return true;
}


bool Population::ras_save_human_info(int gen_num)
{
    std::string sep=" ";
    std::ofstream file_human;
    file_human.open((_out_prefix+".info.pop"+std::to_string(_pop_num+1)+".gen"+ std::to_string(gen_num)+".txt").c_str());
    int npheno=_pheno_scheme.size();
    
    //create header
    file_human << "ID" << sep;
    file_human << "ID_Father" << sep;
    file_human << "ID_Mother" << sep;
    file_human << "ID_Fathers_Father" << sep;
    file_human << "ID_Fathers_Mother" << sep;
    file_human << "ID_Mothers_Father" << sep;
    file_human << "ID_Mothers_Mother" << sep;
    file_human << "sex" << sep;
    for (int j=0; j<npheno; j++)
    {
        file_human << "ph" << j+1 << "_A" << sep;
        file_human << "ph" << j+1 << "_D" << sep;
        file_human << "ph" << j+1 << "_G" << sep;
        file_human << "ph" << j+1 << "_C" << sep;
        file_human << "ph" << j+1 << "_E" << sep;
        file_human << "ph" << j+1 << "_F" << sep;
        file_human << "ph" << j+1 << "_P" << sep;
    }
    file_human << "MV" << sep;
    file_human << "SV" << sep;
    file_human << "SV_f" << std::endl;
    
    unsigned long int n_human=h.size();
    for (unsigned long int i=0; i<n_human; i++)
    {
        // IDs are added 1 in order to start from 1, not zero
        file_human << h[i].ID+1 << sep;
        file_human << h[i].ID_Father+1 << sep;
        file_human << h[i].ID_Mother+1 << sep;
        file_human << h[i].ID_Fathers_Father+1 << sep;
        file_human << h[i].ID_Fathers_Mother+1 << sep;
        file_human << h[i].ID_Mothers_Father+1 << sep;
        file_human << h[i].ID_Mothers_Mother+1 << sep;
        file_human << h[i].sex << sep;
        for (int j=0; j<npheno; j++)
        {
            file_human << h[i].additive[j] << sep;
            file_human << h[i].dominance[j] << sep;
            file_human << h[i].bv[j] << sep;
            file_human << h[i].common_sibling[j] << sep;
            file_human << h[i].e_noise[j] << sep;
            file_human << h[i].parental_effect[j] << sep;
            file_human << h[i].phen[j] << sep;
        }
        file_human << h[i].mating_value << sep;
        file_human << h[i].selection_value << sep;
        file_human << h[i].selection_value_func << std::endl;
    }
    file_human.close();
    return true;
}




std::vector<double> Population::get_additive(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].additive[iphen];
    }
    return x;
}


std::vector<double> Population::get_dominance(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].dominance[iphen];
    }
    return x;
}


std::vector<double> Population::get_common(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].common_sibling[iphen];
    }
    return x;
}


std::vector<double> Population::get_bv(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].bv[iphen];
    }
    return x;
}


std::vector<double> Population::get_e_noise(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].e_noise[iphen];
    }
    return x;
}


std::vector<double> Population::get_parental_effect(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].parental_effect[iphen];
    }
    return x;
}


std::vector<double> Population::get_phen(int iphen)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].phen[iphen];
    }
    return x;
}


std::vector<double> Population::get_mating_value(void)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].mating_value;
    }
    return x;
}


std::vector<double> Population::get_selection_value(void)
{
    unsigned long int nh=h.size();
    std::vector<double> x(nh,0);
    for (unsigned long int i=0; i<nh; i++)
    {
        x[i]=h[i].selection_value;
    }
    return x;
}


// these function needs corrections, A and phen not bv
/*
double Population::compute_couple_var_bv(int sex, int iphen)
{
    unsigned long int n_couples=_couples_info.size();
    if (sex==1)//male
    {
        std::vector<double> bv_male;
        bv_male.reserve(n_couples);
        for (unsigned long int i=0; i<n_couples; i++)
        {
            if(!_couples_info[i].inbreed)
            {
                unsigned long int pos_m=_couples_info[i].pos_male;
                bv_male.push_back(h[pos_m].bv[iphen]);
            }
        }
        return CommFunc::var(bv_male);
    }
    else if (sex==2)//female
    {
        std::vector<double> bv_female;
        bv_female.reserve(n_couples);
        for (unsigned long int i=0; i<n_couples; i++)
        {
            if(!_couples_info[i].inbreed)
            {
                unsigned long int pos_f=_couples_info[i].pos_female;
                bv_female.push_back(h[pos_f].bv[iphen]);
            }
        }
        return CommFunc::var(bv_female);
    }
    else return 0;
}

 */

double Population::compute_couple_cor_bv(int iphen)
{
    unsigned long int n_couples=_couples_info.size();
    std::vector<double> bv_male(n_couples);
    std::vector<double> bv_female(n_couples);
    
    for (unsigned long int i=0; i<n_couples; i++)
    {
        unsigned long int pos_m=_couples_info[i].pos_male;
        unsigned long int pos_f=_couples_info[i].pos_female;
        bv_male[i]=h[pos_m].bv[iphen];
        bv_female[i]=h[pos_f].bv[iphen];
    }
    return CommFunc::cor(bv_male,bv_female);
}


double Population::compute_couple_cor_phen(int iphen)
{
    unsigned long int n_couples = _couples_info.size();
    std::vector<double> phen_male(n_couples);
    std::vector<double> phen_female(n_couples);
    
    for (unsigned long int i=0; i<n_couples; i++)
    {
        unsigned long int pos_m=_couples_info[i].pos_male;
        unsigned long int pos_f=_couples_info[i].pos_female;
        phen_male[i]=h[pos_m].phen[iphen];
        phen_female[i]=h[pos_f].phen[iphen];
    }
    return CommFunc::cor(phen_male,phen_female);
}


double Population::compute_couple_cor_mating_value(void)
{
    unsigned long int n_couples = _couples_info.size();
    std::vector<double> b_male(n_couples);
    std::vector<double> b_female(n_couples);
    
    for (unsigned long int i=0; i<n_couples; i++)
    {
        unsigned long int pos_m=_couples_info[i].pos_male;
        unsigned long int pos_f=_couples_info[i].pos_female;
        b_male[i]=h[pos_m].mating_value;
        b_female[i]=h[pos_f].mating_value;
    }
    return CommFunc::cor(b_male,b_female);
}



