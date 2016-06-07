#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdlib>      // std::rand, std::srand

#include "RasRandomNumber.h"
#include "format_hap.h" // to convert to .hap, .sample .legend file
#include "format_plink.h" // to convert to plink file format


#ifdef _OPENMP
#include <omp.h>
#endif



class part
{
public:
    unsigned long int st;
    unsigned long int en;
    unsigned long int hap_index;
    int root_population; // when we have several populations, the genetic informtion comes from which pupolation
    std::vector<unsigned long int> mutation_pos;
    
    part(void)
    {
        st=0;
        en=0;
        hap_index=0;
        root_population=0;
    }
    part(unsigned long int s,unsigned long int e,unsigned long int i, int r)
    {
        st=s;
        en=e;
        hap_index=i;
        root_population=r;
    }
    bool check_interval(unsigned long int t)
    {
        if (st<=t && t<en) return true;
        else return false;
    }
};

/*
class Chromatid
{
    std::vector<unsigned long int> mutation; // Mutation positions
    //to do: Hap should come here
};
*/


class chromosome
{
public:
    // to do: move Hap to Chromatid
    std::vector<std::vector<part> > Hap; //Hap.size() for humans is 2; for pat and mat
    //std::vector<Chromatid> chromatid; //Chromatid.size() for humans is 2; for pat and mat
    std::vector<double> bv_chr; // each human can have several phenotypes
    std::vector<double> additive_chr; // each human can have several phenotypes
    std::vector<double> dominance_chr; // each human can have several phenotypes
    
/*
    ~chromosome()
    {
        Hap[0].resize(0);
        Hap[1].resize(0);
        Hap[0].shrink_to_fit();
        Hap[1].shrink_to_fit();
        Hap.clear();
    }
*/
};


class Human_CV
{
public:
    struct Chromatid
    {
        std::vector<bool> cv; // size=ncv
        std::vector<double> genetic_value_a; // size=ncv
        std::vector<double> genetic_value_d; // size=ncv
    };
    std::vector<Chromatid> chromatid;
    Human_CV()
    {
    
    }
    Human_CV(unsigned ncv)
    {
        chromatid.resize(2);
        for (int i=0; i<2; i++)
        {
            chromatid[i].cv.resize(ncv);
            chromatid[i].genetic_value_a.resize(ncv);
            chromatid[i].genetic_value_d.resize(ncv);
        }
    }
};



class Human
{
public:
    std::vector<chromosome> chr;
    std::vector<double> additive; // for each pheno
    std::vector<double> dominance; // for each pheno
    std::vector<double> bv; // for each pheno
    std::vector<double> e_noise; // for each pheno
    std::vector<double> common_sibling; // for each pheno
    std::vector<double> parental_effect; // for each pheno
    std::vector<double> phen; // each human can have several phenotypes
    double mating_value;
    double selection_value;
    int gen_num; // the i_th generation
    int sex; // 1=male, 2=female
    unsigned long int ID;
    unsigned long int ID_Father;
    unsigned long int ID_Mother;
    unsigned long int ID_Fathers_Father;
    unsigned long int ID_Fathers_Mother;
    unsigned long int ID_Mothers_Father;
    unsigned long int ID_Mothers_Mother;
/*
    ~Human()
    {
        int nchr=chr.size();
        for (int ichr=0; ichr<nchr; ichr++)
        {
            chr[ichr].Hap[0].resize(0);
            chr[ichr].Hap[1].resize(0);
            chr[ichr].Hap[0].shrink_to_fit();
            chr[ichr].Hap[1].shrink_to_fit();
        }
        //cout << "Human cleared. nchr=" << nchr << endl;
        chr.clear();
    }
 */
};




class Ind_MatingValue
{
public:
    unsigned long int ind;
    double mating_value;
};



class Couples_Info
{
public:
    unsigned long int pos_male; // pos human, not pos hap
    unsigned long int pos_female;
    bool inbreed; // true means can not marry, false mean can marry
    int num_offspring;
    
    Couples_Info()
    {
        pos_male=0;
        pos_female=0;
        inbreed=true;
        num_offspring=0;
    }
};


class rMap // for each chr; recombination map
{
public:
    std::vector<unsigned long int> bp;
    std::vector<double> cM;
    unsigned long int bp_dist_in_rmap;
};


class MutationMap // for each chr; Mutation map
{
public:
    std::vector<unsigned long int> bp;
    std::vector<double> mutation_rate;
};



class CV_INFO // for each chr
{
public:
    std::vector<unsigned long int> bp; // for each cv
    std::vector<double> maf; // for each cv
    std::vector<double> alpha; // for each cv
    std::vector<double> genetic_value_a; // for additive and dominace model
    std::vector<double> genetic_value_d;
};



class CV // for each chr
{
public:
    std::vector<std::vector<bool> > val; // for each hap and cv; i.e., CV.val[10][8] hap10 and cv 8
};



class Phenotype_scheme  // for one phenotype
{
public:
    std::vector<CV_INFO> _cv_info; // for each chr
    std::vector<CV> _cvs; // for each chr
    std::vector<std::string> _name_cv_hap; // for each chr
    double _va;
    double _vd;
    double _ve;
    double _vc;
    double _vf;
    double _omega; // mating value coefficient
    double _beta;  // transmission_of_environmental_effects_from_parents_to_offspring
    double _lambda; // selection value coefficient
};


class Parameters
{
public:
    int _n_pop;
    // pupulations information
    std::vector<std::string> file_gen_info;               // for each population
    std::vector<std::string> file_hap_name;               // for each population
    std::vector<std::string> file_recom_map;              // for each population
    std::vector<std::string> file_mutation_map;           // for each population
    std::vector<std::vector<std::string> > file_cv_info;  // for each population and each phenotype
    std::vector<std::vector<std::string> > file_cvs;      // for each population and each phenotype
    std::vector<std::vector<double> > _va;                // for each population and each phenotype
    std::vector<std::vector<double> > _vd;                // for each population and each phenotype
    std::vector<std::vector<double> > _ve;                // for each population and each phenotype
    std::vector<std::vector<double> > _vc;                // for each population and each phenotype
    std::vector<std::vector<double> > _vf;                // for each population and each phenotype
    std::vector<std::vector<double> > _omega;             // coeff for computing mating value=1
    std::vector<std::vector<double> > _beta;              //  _transmission_of_environmental_effects_from_parents_to_offspring=0;
    std::vector<std::vector<double> > _lambda;            //  selection value coefficient=1
    std::vector<double> _MM_percent;                      // for each population
    std::vector<bool> _RM;                                // for each population
    
    // migration matrix
    std::vector<double> _gamma;                           // for each phenotype and all populations,  _environmental_effects_specific_to_each_population=0
    std::string file_migration;                          // just one for all population
    std::vector<std::vector<double> > migration_mat_gen; //nrow: the number of generations, ncol: _pop_num^2
    
    
    // options
    bool avoid_inbreeding; // if you use, it becomes true and then we have no inbreeding
    bool help;
    bool debug;
    
    // output
    std::string prefix;
    std::string format_output;
    bool no_output;
    bool _interval;
    bool output_all_generations;
    
    // default values
    Parameters(void)
    {
        avoid_inbreeding=false;
        no_output=false;
        _interval=false;
        output_all_generations=false;
        prefix="out";
        help=false;
        debug=false;
        format_output="hap";
        file_migration="";
    }
    void init(int npop)
    {
        _n_pop=npop;
        avoid_inbreeding=false;
        no_output=false;
        _interval=false;
        output_all_generations=false;
        prefix="out";
        help=false;
        debug=false;
        format_output="hap";
        file_migration="";

        file_gen_info.resize(npop);
        file_hap_name.resize(npop);
        file_recom_map.resize(npop);
        file_mutation_map.resize(npop);
        file_cv_info.resize(npop);
        file_cvs.resize(npop);
        _va.resize(npop); // for each population and not phenotype
        _vd.resize(npop);
        _ve.resize(npop);
        _vc.resize(npop);
        _vf.resize(npop);
        _omega.resize(npop);
        _beta.resize(npop);
        _lambda.resize(npop);
        _MM_percent.resize(npop,0);
        _RM.resize(npop,false);
        
        // we can not set _gamma.resize(npheno,0), because _gamma.size=npheno
        //_gamma.resize(npheno,0); // for each phenotype and all populations
    }
    
};


class Population
{
public:
    int _pop_num;
    int _nchr;
    std::vector<Phenotype_scheme> _pheno_scheme; // each population has several pheno, so several Phenotype_scheme
    std::vector<Human> h;
    std::vector<Couples_Info> _couples_info;
    std::vector<unsigned long int> _pop_size; // for each generation
    std::vector<double> _mat_cor;             // for each generation
    std::vector<std::string> _offspring_dist; // for each generation
    std::vector<std::string> _selection_func; // for each generation
    std::vector<double> _selection_func_par1; // for each generation
    std::vector<double> _selection_func_par2; // for each generation
    std::vector<std::vector<std::string> > _hap_legend_sample_name; // for each chr, with 3 columns
    std::vector<rMap> _rmap;                                        // for each chr
    std::vector<MutationMap> _mutation_map;                         // for each chr
    std::vector<std::vector<double> > _recom_prob;                  // for each chr
    std::vector<double> _var_bv_gen0;                               // for each pheno
    std::vector<double> _var_a_gen0;                               // for each pheno
    std::vector<double> _var_d_gen0;                               // for each pheno
    std::vector<int> _all_active_chrs;                              // for each chr
    bool _avoid_inbreeding;
    bool _no_output;
    bool _interval; // for generating interval output
    bool _output_all_generations;
    bool _debug;
    std::string _out_prefix;
    std::string _format_output;
    double _MM_percent; // Random mating percent (inds who have 2 spouses)
    bool _RM; // Random mating
    

    // return
    std::vector<double> ret_var_mating_value;
    std::vector<std::vector<double> > ret_var_parental_effect;
    std::vector<std::vector<double> > ret_var_phen;
    std::vector<std::vector<double> > ret_var_A;
    std::vector<std::vector<double> > ret_var_D;
    std::vector<std::vector<double> > ret_var_C;
    std::vector<std::vector<double> > ret_var_G;
    std::vector<std::vector<double> > ret_var_E;
    std::vector<std::vector<double> > ret_h2;
    
private:
    int ras_get_ind_active_chr(int chr);
    unsigned ras_now_nanoseconds(void);
    
public:
    int ras_read_generation_info_file(std::string file_gen_info);
    int ras_read_hap_legend_sample_address_name(std::string f_name);
    int ras_read_cv_info(std::string f_name, int iphen); // for each phenotype
    int ras_read_cv_info_dominace_model_file(std::string f_name, int iphen);
    int ras_read_cvs_address_name(std::string f_name, int iphen); // for each phenotype
    bool ras_load_cvs(int iphen); // for each phenotype
    unsigned long int ras_read_rmap(std::string f_name);
    unsigned long int ras_read_file_mutation(std::string f_name);
    bool ras_compute_recom_prob(void);
    bool ras_save_human_info(int gen_num);
    std::vector<double> get_bv(int iphen);
    std::vector<double> get_additive(int iphen);
    std::vector<double> get_dominance(int iphen);
    std::vector<double> get_common(int iphen);
    std::vector<double> get_e_noise(int iphen);
    std::vector<double> get_parental_effect(int iphen);
    std::vector<double> get_phen(int iphen);
    std::vector<double> get_mating_value(void);
    std::vector<double> get_selection_value(void);

};



// for migration
class Camp
{
public:
    std::vector<Human> humans;
};



