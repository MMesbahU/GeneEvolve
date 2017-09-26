
#include "parameters.h"



unsigned ras_now_nanoseconds(void)
{
    std::chrono::nanoseconds ns = std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::system_clock::now().time_since_epoch());
    return unsigned(ns.count() % 100000000)+std::rand();
}




bool Parameters::read(std::vector<std::string> &vec_arg)
{
    unsigned ipop=0, npop=1;
    for (unsigned i=1; i<vec_arg.size(); i++)
    {
        if(vec_arg[i]=="--next_population"){
            npop++;
        }
    }
    
    init(npop);
    
    for (unsigned i=1; i<vec_arg.size(); i++)
    {
        if(vec_arg[i]=="--next_population"){
            ipop++;
        }
        //=====================================================================
        // population information
        else if(vec_arg[i]=="--file_gen_info"){
            _file_gen_info[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--file_hap_name"){
            _file_hap_name[ipop]=vec_arg[++i];
            _ref_is_hap = true;
        }
        else if(vec_arg[i]=="--file_ref_vcf"){
            _file_ref_vcf[ipop]=vec_arg[++i];
            _ref_is_vcf = true;
        }
        else if(vec_arg[i]=="--file_recom_map"){
            _file_recom_map[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--file_mutation_map"){
            _file_mutation_map[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--MM"){
            _MM_percent[ipop]=std::stod(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--RM"){
            _RM[ipop]=true;
        }
        else if(vec_arg[i]=="--vt_type"){
            _vt_type = std::stoi(vec_arg[++i]);
        }
        
        /////////////////////////////////////////////////////////
        // for each phenotype
        else if(vec_arg[i]=="--file_cv_info"){ // for each pheno
            _file_cv_info[ipop].push_back(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--file_cvs"){ // for each pheno
            _file_cvs[ipop].push_back(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--va"){ // for each pheno
            _va[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--vd"){ // for each pheno
            _vd[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--vc"){ // for each pheno
            _vc[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--ve"){ // for each pheno
            _ve[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--vf"){ // for each pheno
            _vf[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--omega"){
            _omega[ipop].push_back(std::stod(vec_arg[++i])); // mating value
        }
        else if(vec_arg[i]=="--beta"){
            _beta[ipop].push_back(std::stod(vec_arg[++i])); // transmission_of_environmental_effects_from_parents_to_offspring
        }
        else if(vec_arg[i]=="--lambda"){
            _lambda[ipop].push_back(std::stod(vec_arg[++i])); // selection value
        }
        //=====================================================================
        // for each population
        // _gamma
        else if(vec_arg[i]=="--gamma"){
            _gamma.push_back(std::stod(vec_arg[++i])); // environmental_effects_specific_to_each_population
        }
        //=====================================================================
        // migration matrix
        else if(vec_arg[i]=="--file_migration"){
            _file_migration=vec_arg[++i];
        }
        //=====================================================================
        // other options
        else if(vec_arg[i]=="--avoid_inbreeding"){
            _avoid_inbreeding=true;
        }
        else if(vec_arg[i]=="--seed"){
            _seed=std::stod(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--debug"){
            _debug=true;
        }
        //=====================================================================
        // output
        else if(vec_arg[i]=="--prefix"){
            _prefix=vec_arg[++i];
        }
        else if(vec_arg[i]=="--out_hap"){
            _out_hap=true;
        }
        else if(vec_arg[i]=="--out_plink"){
            _out_plink=true;
        }
        else if(vec_arg[i]=="--out_plink01"){
            _out_plink01=true;
        }
        else if(vec_arg[i]=="--out_vcf"){
            _out_vcf=true;
        }
        else if(vec_arg[i]=="--out_interval"){
            _out_interval=true;
        }
        else if(vec_arg[i]=="--output_all_generations"){
            _output_all_generations=true;
        }
        else if(vec_arg[i]=="--file_output_generations"){
            _file_output_generations=vec_arg[++i];
        }
        //=====================================================================
        // --help
        else if(vec_arg[i]=="--help" || vec_arg[i]=="-h" || vec_arg[i]=="?"){
            _help=true;
        }
        else if(vec_arg[i]=="nothing"){
            // do nothing
        }
        else
        {
            std::cout << " Error: unknown parameter [" << vec_arg[i] << "]" << std::endl;
            return false;
        }
    }
    
    // set default value for optinal vectors:
    for (ipop=0; ipop<npop; ipop++)
    {
        unsigned pop_npheno=_file_cv_info[ipop].size();
        
        // setting va=-1 and vd=-1 means that program will use the variance of real a and d in cv_info file
        // set _va=-1
        if(_va[ipop].size()==0)
        {
            _va[ipop].resize(pop_npheno,-1);
        }
        // set _vd=-1
        if(_vd[ipop].size()==0)
        {
            _vd[ipop].resize(pop_npheno,-1);
        }
        // set _vc=0
        if(_vc[ipop].size()==0)
        {
            _vc[ipop].resize(pop_npheno,0);
        }
        // set _ve=1
        if(_ve[ipop].size()==0)
        {
            _ve[ipop].resize(pop_npheno,1);
        }
        // set _vf=0
        if(_vf[ipop].size()==0)
        {
            _vf[ipop].resize(pop_npheno,0);
        }
        // set _omega=1 // for mating value
        if(_omega[ipop].size()==0)
        {
            _omega[ipop].resize(pop_npheno,1);
        }
        // set _beta=1 // for parental effect to offspring; it is 1, but user can define ve as zero
        if(_beta[ipop].size()==0)
        {
            _beta[ipop].resize(pop_npheno,1);
        }
        // set _lambda=1 // for selection value
        if(_lambda[ipop].size()==0)
        {
            _lambda[ipop].resize(pop_npheno,1);
        }
    }
    // set _gamma=0
    if(_gamma.size()==0)
    {
        _gamma.resize(_file_cv_info[0].size(),0); // set 0 for all phenotypes
    }
    //_seed
    if(_seed==0)
    {
        _seed=ras_now_nanoseconds();
    }
    
    
    return true;
}

bool Parameters::check(void)
{
    if (_file_gen_info[0].size()==0 && !_help)
    {
        return false;
    }
    
    unsigned nphen=_file_cv_info[0].size();
    
    //cehck the number of populations in the inputted options
    for (unsigned ipop=0; ipop< _n_pop; ipop++)
    {
        if (_file_gen_info[ipop].size()==0)
        {
            std::cout << "Error: missing parameter [--file_gen_info] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_file_hap_name[ipop].size()==0 && _file_ref_vcf[ipop].size()==0)
        {
            std::cout << "Error: missing the reference file. Check the parameter [--file_hap_name] or [--file_ref_vcf] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_file_recom_map[ipop].size()==0)
        {
            std::cout << "Error: missing parameter [--file_recom_map] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        
        // check phenotypes prameters
        unsigned pop_npheno=_file_cv_info[ipop].size();
        if (pop_npheno==0)
        {
            std::cout << "Error: missing parameter [--file_cv_info] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_file_cvs[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--file_cvs]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_file_cvs[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--file_cvs]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_va[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--va]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_vd[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--vd]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_vc[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--vc]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_ve[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--ve]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_vf[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--vf]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_omega[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--omega]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_beta[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--beta]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (_lambda[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--lambda]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (pop_npheno!=nphen)
        {
            std::cout << "Error: The number of phenotypes should be the same for each population." << std::endl;
            return false;
        }
        
        // check range for va, vd, ve, vf
        for (unsigned iphen=0; iphen<_va[ipop].size(); iphen++)
        {
            if (!(_va[ipop][iphen]>0 || _va[ipop][iphen]==-1))
            {
                std::cout << "Error: The parameter [--va] should be positive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (unsigned iphen=0; iphen<_vd[ipop].size(); iphen++)
        {
            if (!(_vd[ipop][iphen]>=0 || _vd[ipop][iphen]==-1))
            {
                std::cout << "Error: The parameter [--vd] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (unsigned iphen=0; iphen<_vc[ipop].size(); iphen++)
        {
            if (_vc[ipop][iphen]<0)
            {
                std::cout << "Error: The parameter [--vc] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (unsigned iphen=0; iphen<_ve[ipop].size(); iphen++)
        {
            if (_ve[ipop][iphen]<0)
            {
                std::cout << "Error: The parameter [--ve] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (unsigned iphen=0; iphen<_vf[ipop].size(); iphen++)
        {
            if (_vf[ipop][iphen]<0)
            {
                std::cout << "Error: The parameter [--vf] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        if (_MM_percent[ipop]<0 || _MM_percent[ipop]>1)
        {
            std::cout << "Error: The parameter [--MM] should be between 0 and 1. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        
    } // for ipop
    
    
    //_gamma
    if (_gamma.size()!=nphen)
    {
        std::cout << "Error: the number of [--gamma] must be equal to the number of phenotypes (" << nphen << ")." << std::endl;
        return false;
    }
    
    //_seed
    if (_seed<=0)
    {
        std::cout << "Eroor: the parameter [--seed] should be positive." << std::endl;
    }
    
    //check migration
    if (_file_gen_info.size()>1 && _file_migration.size()==0)
    {
        std::cout << "Error: When you have more than one populations, you must specify the [--file_migration] option." << std::endl;
        return false;
    }
    if (_file_gen_info.size()==1 && _file_migration.size()>0)
    {
        std::cout << " Warning: when there is one population, the parameter [--file_migration] is redundant." << std::endl;
    }
    
    
    return true;
}

bool Parameters::print(void)
{
    if (_help) return true; // no need to pritnt options
    
    std::cout << std::endl;
    std::cout << " Options:" << std::endl;
    std::cout << std::endl;
    unsigned npop=_file_gen_info.size();
    
    for (unsigned ipop=0; ipop<npop; ipop++){
        std::cout << "  Population " << ipop+1 << ":" << std::endl;
        std::cout << "      --file_gen_info          : [" << _file_gen_info[ipop] << "]" << std::endl;
        std::cout << "      --file_hap_name          : [" << _file_hap_name[ipop] << "]" << std::endl;
        std::cout << "      --file_ref_vcf           : [" << _file_ref_vcf[ipop] << "]" << std::endl;
        std::cout << "      --file_recom_map         : [" << _file_recom_map[ipop] << "]" << std::endl;
        std::cout << "      --file_mutation_map      : [" << _file_mutation_map[ipop] << "]" << std::endl;
        std::cout << "      --MM                     : [" << _MM_percent[ipop] << "]" << std::endl;
        std::cout << "      --RM                     : [" << (_RM[ipop] ? "On" : "Off") << "]" << std::endl;
        std::cout << "      --vt_type                : [" << _vt_type << "]" << std::endl;
        
        unsigned pop_npheno=_file_cv_info[ipop].size();
        for (unsigned j=0; j<pop_npheno; j++)
        {
            std::cout << "      phenotype: " << j+1 << std::endl;
            std::cout << "        --file_cv_info         : [" << _file_cv_info[ipop][j] << "]" << std::endl;
            std::cout << "        --file_cvs             : [" << _file_cvs[ipop][j] << "]" << std::endl;
            std::cout << "        --va                   : [" << _va[ipop][j] << "]" << std::endl;
            std::cout << "        --vd                   : [" << _vd[ipop][j] << "]" << std::endl;
            std::cout << "        --vc                   : [" << _vc[ipop][j] << "]" << std::endl;
            std::cout << "        --ve                   : [" << _ve[ipop][j] << "]" << std::endl;
            std::cout << "        --vf                   : [" << _vf[ipop][j] << "]" << std::endl;
            std::cout << "        --omega                : [" << _omega[ipop][j] << "]" << std::endl;
            std::cout << "        --lambda               : [" << _lambda[ipop][j] << "]" << std::endl;
            std::cout << "        --beta                 : [" << _beta[ipop][j] << "]" << std::endl;
        }
    }
    
    std::cout << "  Immigration parameters" << std::endl;
    std::cout << "      --file_migration         : [" << _file_migration << "]" << std::endl;
    
    std::cout << "  Environmental effects specific to each population (for each phenotype)" << std::endl;
    unsigned pop_npheno=_file_cv_info[0].size();
    for (unsigned j=0; j<pop_npheno; j++)
    {
        std::cout << "      --gamma                  : [" << _gamma[j] << "]" << std::endl;
    }
    
    std::cout << "  Other parameters" << std::endl;
    std::cout << "      --seed                   : [" << _seed << "]" << std::endl;
    std::cout << "      --avoid_inbreeding       : [" << (_avoid_inbreeding ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --debug                  : [" << (_debug ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --prefix                 : [" << _prefix << "]" << std::endl;
    std::cout << "      --out_hap                : [" << (_out_hap ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --out_plink              : [" << (_out_plink ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --out_plink01            : [" << (_out_plink01 ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --out_vcf                : [" << (_out_vcf ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --out_interval           : [" << (_out_interval ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --output_all_generations : [" << (_output_all_generations ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --file_output_generations: [" << _file_output_generations << "]" << std::endl;
    std::cout << std::endl;
    return true;
}





