#include <iostream>
#include <string> // str1.compare(str2)
#include <vector>
#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds



class Parameters
{
    public:
    int _n_pop;
    unsigned _seed;
    // pupulations information
    std::vector<std::string> _file_gen_info;               // for each population
    std::vector<std::string> _file_hap_name;               // for each population
    std::vector<std::string> _file_ref_vcf;                // for each population
    std::vector<std::string> _file_recom_map;              // for each population
    std::vector<std::string> _file_mutation_map;           // for each population
    std::vector<std::vector<std::string> > _file_cv_info;  // for each population and each phenotype
    std::vector<std::vector<std::string> > _file_cvs;      // for each population and each phenotype
    std::vector<std::vector<double> > _va;                // for each population and each phenotype
    std::vector<std::vector<double> > _vd;                // for each population and each phenotype
    std::vector<std::vector<double> > _ve;                // for each population and each phenotype
    std::vector<std::vector<double> > _vc;                // for each population and each phenotype
    std::vector<std::vector<double> > _vf;                // for each population and each phenotype
    std::vector<std::vector<double> > _omega;             // coeff for computing mating value=1
    std::vector<std::vector<double> > _beta;              // _transmission_of_environmental_effects_from_parents_to_offspring=0;
    std::vector<std::vector<double> > _lambda;            //  selection value coefficient=1
    std::vector<double> _MM_percent;                      // for each population
    std::vector<bool> _RM;                                // for each population
    int _vt_type;                                         // type of vertical transition (familial effect)
    bool _ref_is_hap;
    bool _ref_is_vcf;
    
    // migration matrix
    std::vector<double> _gamma;                           // for each phenotype and all populations,  _environmental_effects_specific_to_each_population=0
    std::string _file_migration;                          // just one for all population
    std::vector<std::vector<double> > _migration_mat_gen; //nrow: the number of generations, ncol: _pop_num^2
    
    
    // options
    bool _avoid_inbreeding; // if you use, it becomes true and then we have no inbreeding
    bool _help;
    bool _debug;
    
    // output
    std::string _prefix;
    bool _out_hap;
    bool _out_plink;
    bool _out_plink01;
    bool _out_vcf;
    bool _out_interval;
    bool _output_all_generations;
    std::string _file_output_generations;
    
    // default values
    Parameters(void)
    {
        _avoid_inbreeding=false;
        _out_hap=false;
        _out_plink=false;
        _out_plink01=false;
        _out_interval=false;
        _out_vcf=false;
        _output_all_generations=false;
        _prefix="out";
        _help=false;
        _debug=false;
        _file_migration="";
        _file_output_generations="";
        _ref_is_vcf = false;
        _ref_is_hap = false;
    }
    void init(int npop)
    {
        _n_pop=npop;
        _seed=0;
        _avoid_inbreeding=false;
        _out_hap=false;
        _out_plink=false;
        _out_plink01=false;
        _out_interval=false;
        _out_vcf=false;
        _output_all_generations=false;
        _prefix="out";
        _help=false;
        _debug=false;
        _file_migration="";
        
        _file_gen_info.resize(npop);
        _file_hap_name.resize(npop);
        _file_ref_vcf.resize(npop);
        _file_recom_map.resize(npop);
        _file_mutation_map.resize(npop);
        _file_cv_info.resize(npop);
        _file_cvs.resize(npop);
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
        _file_output_generations="";
        _vt_type = 1;
        _ref_is_vcf = false;
        _ref_is_hap = false;

        // we can not set _gamma.resize(npheno,0), because _gamma.size=npheno
        //_gamma.resize(npheno,0); // for each phenotype and all populations
    }
    
    
    public:
    bool read(std::vector<std::string> &vec_arg);
    bool check(void);
    bool print(void);

    
};


