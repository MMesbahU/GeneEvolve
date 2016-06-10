#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <ostream>
#include <iostream>
#include <ctime>
#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds
#include <string> // str1.compare(str2)
//#include "Parameters.h"
//#include "StringBasics.h"
//#include "HaplotypeSet.h"
#include "Simulation.h"
//#include "hap.h" it is in Simulation.h
#include "CommFunc.h"
//#include "Unique.h"

//using Eigen::MatrixXd;
//using namespace Eigen;

int transFactor = 3;
int cisFactor = 2;


void ProgramVersion();
void helpFile();
bool parameter_proc(std::vector<std::string> &vec_arg, Parameters &par);
bool parameter_check(Parameters &par);
bool parameter_print(Parameters &par);
unsigned ras_now_nanoseconds(void);


int main(int argc, char ** argv)
{
    int cpus;
    
    int start_time = time(0);

#ifdef _OPENMP
    omp_set_num_threads(cpus);
#else
    cpus=1;
#endif
    
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    ProgramVersion();

    
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    // parameters
    Parameters par;
    std::vector<std::string> vec_arg(argc+2); // we add 2 nothing
    for(int i=0; i<argc; i++) vec_arg[i]=(std::string) argv[i];
    vec_arg[argc]="nothing"; // for checking
    vec_arg[argc+1]="nothing";

    if (!parameter_proc(vec_arg, par))
    {
        std::cout << std::endl;
        std::cout << " For more information type [GeneEvolve --help]." << std::endl;
        std::cout << std::endl;
        return -1;
    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    if (par.help)
    {
        helpFile();
        return -1;
    }
    
    if (!parameter_check(par))
    {
        std::cout << std::endl;
        std::cout << " For more information type [GeneEvolve --help]." << std::endl;
        std::cout << std::endl;
        return -1;
    }
    
    
    parameter_print(par);
    
    
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    // main simulation
    Simulation simulation;
    simulation.par=par;

    if(!simulation.run())
    {
        std::cout << "Program exits with an error." << std::endl;
        return -1;
    }
    
    //if(log)
    //    fclose(LogFile);

    
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    int time_tot = time(0) - start_time;
    std::cout << "\n Program Successfully finished.\n ";
    printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
           time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);
    std::cout << "\n Thank You for using GeneEvolve !!! " << std::endl << std::endl;

	return 0;
}




void ProgramVersion()
{
    printf("\n");
    printf(" ------------------------------------------------------------------------------ \n");
	printf("                                  GeneEvolve                                    \n");
	printf(" ------------------------------------------------------------------------------ \n");
    printf(" (c) 2015 - Rasool Tahmasbi and Matthew C. Keller.\n");
    std::cout << "\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
    std::cout << std::endl;
    printf(" URL = https://github.com/rtahmasbi/GeneEvolve\n");
    printf(" URL = http://matthewckeller.com/html/GeneEvolve.html\n");
}

void helpFile()
{
    printf("\n\n About:\n");
    printf("\n\t A fast and memory efficient forward-time simulator of whole-genome data.\n");
    printf("\n\n");
    printf(" GeneEvolve is a user-friendly and efficient population genetics simulator that handles\n");
    printf(" complex evolutionary scenarios and generates individual-level phenotypes and realistic\n");
    printf(" whole-genome sequence or SNP data. GeneEvolve runs forward-in-time, which allows it to\n");
    printf(" provide a wide range of scenarios for mating systems, selection, population size and\n");
    printf(" structure, migration, recombination, and environmental effects.\n");
   
    
    printf("\n\n\n");
    printf(" ----------------------------------------------------------------------------------------- \n");
	printf("                            GeneEvolve - List of Usage Options \n");
	printf(" -----------------------------------------------------------------------------------------\n\n");

    std::cout << std::endl;
    printf(" --------- Population information\n");
    printf("        --file_gen_info          : [filename]\n");
    printf("          File containing information about generations, each line represent one generation.\n");
    printf("        --file_hap_name          : [filename]\n");
    printf("        --file_recom_map         : [filename]\n");
    printf("        --file_mutation_map      : [filename]\n");
    printf("        --selection_func         : [logit 0 1]\n");
    printf("        --RM                     : [off]  ->  Random Mating\n");
    printf("        --MM                     : [0]\n");
    printf("          Percentage of individuals who have more than 1 spouse (0<=MM<=1).\n");
    printf("        --next_population        :\n");
    printf("          This keyword can be used to distinguish between populations.\n");

    std::cout << std::endl;
    printf(" --------- Phenotypes\n");
    printf("        --file_cv_info           : [filename]\n");
    printf("        --file_cvs               : [filename]\n");
    printf("        --va                     : [-1]   ->  variance of additive effect (--va -1 means the program will not transform the variance)\n");
    printf("        --vd                     : [-1]   ->  variance of dominance effect (--vd -1 means the program will not transform the variance)\n");
    printf("        --vc                     : [0]    ->  variance of sibling effect (common effect)\n");
    printf("        --ve                     : [1]    ->  variance of envirorment effect\n");
    printf("        --vf                     : [0]    ->  variance of familial effect\n");
    printf("        --omega                  : [1]    ->  coefficient for mating value\n");
    printf("        --lambda                 : [1]    ->  coefficient for selection value\n");
    printf("        --beta                   : [1]    ->  coefficient for familial effect\n");

    std::cout << std::endl;
    printf(" --------- Immigration parameters\n");
    printf("        --file_migration         : [filename]\n");
    
    std::cout << std::endl;
    printf(" --------- Environmental effects specific to each population (for each phenotype)\n");
    printf("        --gamma                  : [0]\n");
    
    std::cout << std::endl;
    printf(" --------- Other parameters\n");
    printf("        --seed                   : [0]  ->  can be any passitive number.\n");
    printf("        --avoid_inbreeding       : [Off]  ->  [On] means no inbreeding.\n");
    printf("        --debug                  : [Off]\n");
    printf("        --prefix                 : [out]\n");
    printf("        --format_output          : [hap]  ->  can be [hap], [plink] or [interval].\n");
    printf("        --interval               : [Off]  ->  [On] means also output the [interval] format.\n");
    printf("        --no_output              : [Off]\n");
    printf("        --output_all_generations : [Off]\n");
    

    printf("\n Please visit <https://github.com/rtahmasbi/GeneEvolve> for detailed documentation.\n\n");
    std::cout << std::endl;



	return;
}


bool parameter_proc(std::vector<std::string> &vec_arg, Parameters &par)
{
    int ipop=0, npop=1;
    for (int i=1; i<(int)vec_arg.size(); i++)
    {
        if(vec_arg[i]=="--next_population"){
            npop++;
        }
    }
    
    par.init(npop);
    
    for (int i=1; i<(int)vec_arg.size(); i++)
    {
        if(vec_arg[i]=="--next_population"){
            ipop++;
        }
        //=====================================================================
        // population information
        else if(vec_arg[i]=="--file_gen_info"){
            par.file_gen_info[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--file_hap_name"){
            par.file_hap_name[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--file_recom_map"){
            par.file_recom_map[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--file_mutation_map"){
            par.file_mutation_map[ipop]=vec_arg[++i];
        }
        else if(vec_arg[i]=="--MM"){
            par._MM_percent[ipop]=std::stod(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--RM"){
            par._RM[ipop]=true;
        }
        /////////////////////////////////////////////////////////
        // for each phenotype
        else if(vec_arg[i]=="--file_cv_info"){ // for each pheno
            par.file_cv_info[ipop].push_back(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--file_cvs"){ // for each pheno
            par.file_cvs[ipop].push_back(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--va"){ // for each pheno
            par._va[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--vd"){ // for each pheno
            par._vd[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--vc"){ // for each pheno
            par._vc[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--ve"){ // for each pheno
            par._ve[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--vf"){ // for each pheno
            par._vf[ipop].push_back(std::stod(vec_arg[++i]));
        }
        else if(vec_arg[i]=="--omega"){
            par._omega[ipop].push_back(std::stod(vec_arg[++i])); // mating value
        }
        else if(vec_arg[i]=="--beta"){
            par._beta[ipop].push_back(std::stod(vec_arg[++i])); // transmission_of_environmental_effects_from_parents_to_offspring
        }
        else if(vec_arg[i]=="--lambda"){
            par._lambda[ipop].push_back(std::stod(vec_arg[++i])); // selection value
        }
        //=====================================================================
        // for each population
        // _gamma
        else if(vec_arg[i]=="--gamma"){
            par._gamma.push_back(std::stod(vec_arg[++i])); // environmental_effects_specific_to_each_population
        }
        //=====================================================================
        // migration matrix
        else if(vec_arg[i]=="--file_migration"){
            par.file_migration=vec_arg[++i];
        }
        //=====================================================================
        // other options
        else if(vec_arg[i]=="--avoid_inbreeding"){
            par.avoid_inbreeding=true;
        }
        else if(vec_arg[i]=="--seed"){
            par._seed=std::stod(vec_arg[++i]);
        }
        else if(vec_arg[i]=="--debug"){
            par.debug=true;
        }
        //=====================================================================
        // output
        else if(vec_arg[i]=="--prefix"){
            par.prefix=vec_arg[++i];
        }
        else if(vec_arg[i]=="--format_output"){
            par.format_output=vec_arg[++i];
        }
        else if(vec_arg[i]=="--no_output"){
            par.no_output=true;
        }
        else if(vec_arg[i]=="--interval"){
            par._interval=true;
        }
        else if(vec_arg[i]=="--output_all_generations"){
            par.output_all_generations=true;
        }
        //=====================================================================
        // --help
        else if(vec_arg[i]=="--help" || vec_arg[i]=="-h" || vec_arg[i]=="?"){
            par.help=true;
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
        unsigned pop_npheno=par.file_cv_info[ipop].size();

        // setting va=-1 and vd=-1 means that program will use the variance of real a and d in cv_info file
        // set _va=-1
        if(par._va[ipop].size()==0)
        {
            par._va[ipop].resize(pop_npheno,-1);
        }
        // set _vd=-1
        if(par._vd[ipop].size()==0)
        {
            par._vd[ipop].resize(pop_npheno,-1);
        }
        // set _vc=0
        if(par._vc[ipop].size()==0)
        {
            par._vc[ipop].resize(pop_npheno,0);
        }
        // set _ve=1
        if(par._ve[ipop].size()==0)
        {
            par._ve[ipop].resize(pop_npheno,1);
        }
        // set _vf=0
        if(par._vf[ipop].size()==0)
        {
            par._vf[ipop].resize(pop_npheno,0);
        }
        // set _omega=1 // for mating value
        if(par._omega[ipop].size()==0)
        {
            par._omega[ipop].resize(pop_npheno,1);
        }
        // set _beta=1 // for parental effect to offspring; it is 1, but user can define ve as zero
        if(par._beta[ipop].size()==0)
        {
            par._beta[ipop].resize(pop_npheno,1);
        }
        // set _lambda=1 // for selection value
        if(par._lambda[ipop].size()==0)
        {
            par._lambda[ipop].resize(pop_npheno,1);
        }
    }
    // set _gamma=0
    if(par._gamma.size()==0)
    {
        par._gamma.resize(par.file_cv_info[0].size(),0); // set 0 for all phenotypes
    }
    //_seed
    if(par._seed==0)
    {
        par._seed=ras_now_nanoseconds();
    }
    
    
    return true;
}

bool parameter_check(Parameters &par)
{
    if (par.file_gen_info[0].size()==0 && !par.help)
    {
        return false;
    }
    
    unsigned nphen=par.file_cv_info[0].size();
    
    //cehck the number of populations in the inputted options
    for (int ipop=0; ipop< par._n_pop; ipop++)
    {
        if (par.file_gen_info[ipop].size()==0)
        {
            std::cout << "Error: missing parameter [--file_gen_info] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par.file_hap_name[ipop].size()==0)
        {
            std::cout << "Error: missing parameter [--file_hap_name] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par.file_recom_map[ipop].size()==0)
        {
            std::cout << "Error: missing parameter [--file_recom_map] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        
        // check phenotypes prameters
        unsigned pop_npheno=par.file_cv_info[ipop].size();
        if (pop_npheno==0)
        {
            std::cout << "Error: missing parameter [--file_cv_info] in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par.file_cvs[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--file_cvs]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par.file_cvs[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--file_cvs]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._va[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--va]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._vd[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--vd]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._vc[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--vc]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._ve[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--ve]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._vf[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--vf]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._omega[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--omega]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._beta[ipop].size()!=pop_npheno)
        {
            std::cout << "Error: each phenotype needs one [--beta]. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        if (par._lambda[ipop].size()!=pop_npheno)
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
        for (int iphen=0; iphen<(int)par._va[ipop].size(); iphen++)
        {
            if (!(par._va[ipop][iphen]>0 || par._va[ipop][iphen]==-1))
            {
                std::cout << "Error: The parameter [--va] should be positive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (int iphen=0; iphen<(int)par._vd[ipop].size(); iphen++)
        {
            if (!(par._vd[ipop][iphen]>=0 || par._vd[ipop][iphen]==-1))
            {
                std::cout << "Error: The parameter [--vd] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (int iphen=0; iphen<(int)par._vc[ipop].size(); iphen++)
        {
            if (par._vc[ipop][iphen]<0)
            {
                std::cout << "Error: The parameter [--vc] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (int iphen=0; iphen<(int)par._ve[ipop].size(); iphen++)
        {
            if (par._ve[ipop][iphen]<0)
            {
                std::cout << "Error: The parameter [--ve] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        for (int iphen=0; iphen<(int)par._vf[ipop].size(); iphen++)
        {
            if (par._vf[ipop][iphen]<0)
            {
                std::cout << "Error: The parameter [--vf] should not be negetive. Error in population " << ipop+1 << "." << std::endl;
                return false;
            }
        }
        if (par._MM_percent[ipop]<0 || par._MM_percent[ipop]>1)
        {
            std::cout << "Error: The parameter [--MM] should be between 0 and 1. Error in population " << ipop+1 << "." << std::endl;
            return false;
        }
        
    } // for ipop

    
    //_gamma
    if (par._gamma.size()!=nphen)
    {
        std::cout << "Error: the number of [--gamma] must be equal to the number of phenotypes (" << nphen << ")." << std::endl;
        return false;
    }

    //_seed
    if (par._seed<0)
    {
        std::cout << "Eroor: the parameter [--seed] can't be negative." << std::endl;
    }
    
    //check migration
    if (par.file_gen_info.size()>1 && par.file_migration.size()==0)
    {
        std::cout << "Error: When you have more than one populations, you must specify the [--file_migration] option." << std::endl;
        return false;
    }
    if (par.file_gen_info.size()==1 && par.file_migration.size()>0)
    {
        std::cout << " Warning: when there is one population, the parameter [--file_migration] is redundant." << std::endl;
    }
    
    //check par.format_output
    if (!(par.format_output=="hap" || par.format_output=="plink" || par.format_output=="interval"))
    {
        std::cout << "Error: unknown output format [" +par.format_output+ "]." << std::endl;
        return false;
    }
    
    return true;
}

bool parameter_print(Parameters &par)
{
    if (par.help) return true; // no need to pritnt options
    
    std::cout << std::endl;
    std::cout << " Options:" << std::endl;
    std::cout << std::endl;
    int npop=par.file_gen_info.size();
    
    for (int ipop=0; ipop<npop; ipop++){
        std::cout << "  Population " << ipop+1 << ":" << std::endl;
        std::cout << "      --file_gen_info          : [" << par.file_gen_info[ipop] << "]" << std::endl;
        std::cout << "      --file_hap_name          : [" << par.file_hap_name[ipop] << "]" << std::endl;
        std::cout << "      --file_recom_map         : [" << par.file_recom_map[ipop] << "]" << std::endl;
        std::cout << "      --file_mutation_map      : [" << par.file_mutation_map[ipop] << "]" << std::endl;
        std::cout << "      --MM                     : [" << par._MM_percent[ipop] << "]" << std::endl;
        std::cout << "      --RM                     : [" << (par._RM[ipop] ? "On" : "Off") << "]" << std::endl;
        
        int pop_npheno=par.file_cv_info[ipop].size();
        for (int j=0; j<pop_npheno; j++)
        {
            std::cout << "      phenotype: " << j+1 << std::endl;
            std::cout << "        --file_cv_info         : [" << par.file_cv_info[ipop][j] << "]" << std::endl;
            std::cout << "        --file_cvs             : [" << par.file_cvs[ipop][j] << "]" << std::endl;
            std::cout << "        --va                   : [" << par._va[ipop][j] << "]" << std::endl;
            std::cout << "        --vd                   : [" << par._vd[ipop][j] << "]" << std::endl;
            std::cout << "        --vc                   : [" << par._vc[ipop][j] << "]" << std::endl;
            std::cout << "        --ve                   : [" << par._ve[ipop][j] << "]" << std::endl;
            std::cout << "        --vf                   : [" << par._vf[ipop][j] << "]" << std::endl;
            std::cout << "        --omega                : [" << par._omega[ipop][j] << "]" << std::endl;
            std::cout << "        --lambda               : [" << par._lambda[ipop][j] << "]" << std::endl;
            std::cout << "        --beta                 : [" << par._beta[ipop][j] << "]" << std::endl;
        }
    }

    std::cout << "  Immigration parameters" << std::endl;
    std::cout << "      --file_migration        : [" << par.file_migration << "]" << std::endl;

    std::cout << "  Environmental effects specific to each population (for each phenotype)" << std::endl;
    int pop_npheno=par.file_cv_info[0].size();
    for (int j=0; j<pop_npheno; j++)
    {
        std::cout << "      --gamma                  : [" << par._gamma[j] << "]" << std::endl;
    }
    
    std::cout << "  Other parameters" << std::endl;
    std::cout << "      --seed                   : [" << par._seed << "]" << std::endl;
    std::cout << "      --avoid_inbreeding       : [" << (par.avoid_inbreeding ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --debug                  : [" << (par.debug ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --prefix                 : [" << par.prefix << "]" << std::endl;
    std::cout << "      --format_output          : [" << par.format_output << "]" << std::endl;
    std::cout << "      --no_output              : [" << (par.no_output ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --interval               : [" << (par._interval ? "On" : "Off") << "]" << std::endl;
    std::cout << "      --output_all_generations : [" << (par.output_all_generations ? "On" : "Off") << "]" << std::endl;
    std::cout << std::endl;
    return true;
}



unsigned ras_now_nanoseconds(void)
{
    std::chrono::nanoseconds ns = std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::system_clock::now().time_since_epoch());
    return unsigned(ns.count() % 100000000)+std::rand();
}


