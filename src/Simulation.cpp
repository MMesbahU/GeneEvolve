#include <random>
#include "Simulation.h"
#include <exception>
#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds
#include <sstream>      // std::istringstream
#include "CommFunc.h"
#include "RasRandomNumber.h"
#include <unistd.h> // for sysconf



// todo: change bv to G


bool ValueCmp(Ind_MatingValue const &a, Ind_MatingValue const &b)
{
    return (a.mating_value < b.mating_value);
};


// Returns f'(x) given a pointer to f, and value x.
double Simulation::Fprime(int iphen, double x)
{
    const double dx = 0.001; // x/1000.0; //small_delta( x, 5 );
    return ( ras_combined_variance(iphen, x+dx) - ras_combined_variance(iphen, x-dx) ) / (2*dx);
}


// Solves f(x) = 0 for x,
// Provide: function pointer f, a guess (x0), and the required precision.
double Simulation::NewtonRaphson(int iphen, double x0, double precision)
{
    double fx0 = ras_combined_variance(iphen, x0);
    double x1 = x0 - fx0 / Fprime(iphen, x0);
    double fx1 = ras_combined_variance(iphen, x1);
    
    //std::cout << "fx0=" << fx0 << ", x0=" << x0 << std::endl;
    //std::cout << "fx1=" << fx1 << ", x1=" << x1 << std::endl;
    //std::cout << "precision=" << precision << std::endl;

    //if ( std::abs(x1 - x0) < precision )
    //   return x1;
    
    if ( std::abs(fx1) < precision )
    {
        return x1;
    }
    
    return NewtonRaphson(iphen, x1 , precision); // Recursion baby!
}




bool Simulation::run(void)
{
    int start_time = time(0);
    int time_prev = start_time;
    int time_load;
    
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    std::cout << "                                INITIALIZATION                                 " << std::endl;
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    
    if( !ras_init_parameters() )
        return false;
    
    std::cout << std::endl;
    
    if( !ras_init_generation0() )
        return false;
    
    time_load = time(0) - time_prev;
    std::cout << std::endl << " Time taken for initialization = " << time_load << " seconds." << std::endl;
    
    
    
    
    
    time_prev = time(0);
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    std::cout << "                                MAIN PROCEDURE                                 " << std::endl;
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    
    if ( !ras_main_sim() )
    {
        std::cout << "Error in ras_main_sim() function!" << std::endl;
        return false;
    }
    time_load = time(0) - time_prev;
    std::cout << std::endl <<  " Time taken for simulation = " << time_load << " seconds." << std::endl;
    
    

    
    
    // give some info
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    std::cout << "                                    RESULTS                                    " << std::endl;
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    
    ras_show_res();
    ras_save_res();
    
    std::cout << std::endl;
    
    
    

    // read the hap file and write the output for the last generation
    if (!par.no_output && !par.output_all_generations)
    {
        time_prev = time(0);
        std::cout << " ------------------------------------------------------------------------------" << std::endl;
        std::cout << "                                    OUTPUTS                                    " << std::endl;
        std::cout << " ------------------------------------------------------------------------------" << std::endl;
        
        if(!ras_save_genotypes(_tot_gen))
        {
            return false;
        }
        
        time_load = time(0) - time_prev;
        std::cout << " Time taken for reading and writing = " << time_load << " seconds." << std::endl;
    }

   
    
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    std::cout << "                                END OF PROGRAM                                 " << std::endl;
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
   
    return true;

}


bool Simulation::ras_init_parameters(void)
{
    _n_pop=par._n_pop;
    population.resize(_n_pop);
    _tot_gen=0;
    _output_all_generations=par.output_all_generations;
    _format_output=par.format_output;
    _out_prefix=par.prefix;
    _debug=par.debug;
    _interval=par._interval;
    
    
    _Pop_info_prev_gen.resize(_n_pop); // to save mating_value for the previous generation, for each population
    _gen0_SV_var.resize(_n_pop);
    _gen0_SV_mean.resize(_n_pop);
    //it will be used for transmition of parental pheno to offsprings
    

    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        std::cout << " Population " << ipop+1 << std::endl;
        
        population[ipop]._pop_num=ipop;
        population[ipop]._avoid_inbreeding=par.avoid_inbreeding;
        population[ipop]._no_output=par.no_output;
        population[ipop]._interval=par._interval; // for interval output
        population[ipop]._output_all_generations=par.output_all_generations;
        population[ipop]._debug=par.debug;
        population[ipop]._out_prefix=par.prefix;
        population[ipop]._format_output=par.format_output;
        population[ipop]._MM_percent=par._MM_percent[ipop];
        population[ipop]._RM=par._RM[ipop];

        // phenotypes
        int nphen = par._va[ipop].size();
        population[ipop]._pheno_scheme.resize(nphen);
        for (int iphen=0; iphen<nphen; iphen++)
        {
            population[ipop]._pheno_scheme[iphen]._va = par._va[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._vd = par._vd[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._vc = par._vc[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._ve = par._ve[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._vf = par._vf[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._omega = par._omega[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._beta = par._beta[ipop][iphen];
            population[ipop]._pheno_scheme[iphen]._lambda = par._lambda[ipop][iphen];
        }

        ///////////////////////////
        //ras_read_gene_info ; genertion
        int tot_gen_temp=population[ipop].ras_read_generation_info_file(par.file_gen_info[ipop]);
        if (tot_gen_temp==0)
        {
            std::cout << "Error: The number of generations should be > 0." << std::endl;
            return false;
        }
        if (_tot_gen>0 && _tot_gen!=tot_gen_temp)
        {
            std::cout << "Error: The number of generations in each population differ.";
            return false;
        }
        _tot_gen=tot_gen_temp;
        std::cout << "     Number of generations            = " << _tot_gen << std::endl;
        
        
        ///////////////////////////
        //ras_read_hap_legend_sample_address_name
        int nl=population[ipop].ras_read_hap_legend_sample_address_name(par.file_hap_name[ipop]);
        int nchr=population[ipop]._hap_legend_sample_name.size();
        population[ipop]._nchr=nchr;
        if(nl!=nchr)
        {
            std::cout << "Error: Number of chromosomes is not equal to the input file." << std::endl;
            return false;
        }
        std::cout << "     Number of chromosomes            = " << nchr << std::endl;
        
        

        // phenotypes
        int pop_npheno = par.file_cv_info[ipop].size();
        
        // alloc returns
        population[ipop].ret_var_phen.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_A.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_D.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_C.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_G.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_E.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_h2.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_parental_effect.resize(pop_npheno, std::vector<double> (_tot_gen+1,0));
        population[ipop].ret_var_mating_value.resize(_tot_gen+1,0);

        population[ipop]._var_bv_gen0.resize(pop_npheno);
        population[ipop]._var_a_gen0.resize(pop_npheno);
        population[ipop]._var_d_gen0.resize(pop_npheno);
        std::cout << "     Number of phenotypes             = " << pop_npheno << std::endl;
        population[ipop]._pheno_scheme.resize(pop_npheno);
        for (int iphen=0; iphen<pop_npheno; iphen++)
        {
            std::cout << "     phenotype: " << iphen+1 << std::endl;

            ///////////////////////////
            //ras_read_cv_info
            //int ncv = population[ipop].ras_read_cv_info(par.file_cv_info[ipop][iphen],iphen);
            int ncv = population[ipop].ras_read_cv_info_dominace_model_file(par.file_cv_info[ipop][iphen],iphen);
            std::cout << "       Number of CVs                  = " << ncv << std::endl;
            if (ncv==0) return false;
            
            ///////////////////////////
            //ras_read_cvs_address_name
            int ncv_hapfile = population[ipop].ras_read_cvs_address_name(par.file_cvs[ipop][iphen],iphen);
            std::cout << "       Number of cv hap files         = " << ncv_hapfile << std::endl;
            if (ncv_hapfile==0) return false;

            ///////////////////////////
            //load cv values
            bool ret=population[ipop].ras_load_cvs(iphen); // name of cv hap file is extracted in func ras_read_cvs_address_name
            if(!ret) return false;
            
            if(_debug)
            {
                for (int ideb=0; ideb<40; ideb++)
                {
                    std::cout << population[ipop]._pheno_scheme[iphen]._cvs[0].val[0][ideb] << " ";
                }
                std::cout << std::endl;
            }
        }

        ///////////////////////////
        //read ras_read_rmap
        long int tot_recom_map = population[ipop].ras_read_rmap(par.file_recom_map[ipop]); // we inputed just one rec map file
        std::cout << "     Number of recombination map snps = " << tot_recom_map << std::endl;
        if (tot_recom_map==0) return false;

        
        ///////////////////////////
        //recom_prob
        population[ipop].ras_compute_recom_prob();
        
        
        
        ///////////////////////////
        //check CVs are in genomic map or not
        for (int iphen=0; iphen<pop_npheno; iphen++)
        {
            int nchr = population[ipop]._pheno_scheme[iphen]._cv_info.size();
            for (int ichr=0; ichr<nchr; ichr++)
            {
                unsigned long int rmap_st = population[ipop]._rmap[ichr].bp[0];
                int nrmap_chr=population[ipop]._rmap[ichr].bp.size();
                unsigned long int rmap_en = population[ipop]._rmap[ichr].bp[nrmap_chr-1];
                int ncv_chr=population[ipop]._pheno_scheme[iphen]._cv_info.size();
                for (int icv=0; icv<ncv_chr; icv++)
                {
                    unsigned long int cv_bp = population[ipop]._pheno_scheme[iphen]._cv_info[ichr].bp[icv];
                    if (cv_bp<rmap_st || cv_bp>rmap_en)
                    {
                        std::cout << "Error: CVs should be in range of genomic map. Error in ichr " << ichr << std::endl;
                        return false;
                    }
                }
            }

        }
        
        ///////////////////////////
        //read ras_read_file_mutation
        if (par.file_mutation_map[ipop].size() >0)
        {
            // user inputed the mutation_map_file
            long int tot_mutation_map = population[ipop].ras_read_file_mutation(par.file_mutation_map[ipop]);
            std::cout << "     Number of intervals in the mutation map file = " << tot_mutation_map << std::endl;
            if (tot_mutation_map==0) return false;
        }
        
    }
    
    int nphen = par._va[0].size();
    _gamma.resize(nphen);
    for (int iphen=0; iphen<nphen; iphen++)
    {
        _gamma[iphen] = par._gamma[iphen];
    }
    
    
    if (_n_pop>1)
    {
        if (!read_migration_file())
        {
            return false;
        }
    }
    
    _all_active_chrs=population[0]._all_active_chrs;
    // todo: check they are equal for all the populations

    return true;
}



bool Simulation::ras_init_generation0(void)
{
    int nphen = par._va[0].size();
    // create genotype and phenotypes for gen0
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        ///////////////////////////
        //human
        ras_initial_human_gen0(ipop);
        
        ///////////////////////////////////////////////////////
        //computing breeding value
        std::cout << "      computing breeding value" << std::endl;
        if(!ras_compute_breeding_value(ipop) ) return false;

        
        ///////////////////////////
        // in order to use the parental effect, here we should
        // fill _Pop_info_prev_gen with some random numbers
        // this will be used in ras_create_pheno()
        ras_fill_Pop_info_prev_gen_for_gen0_prev(ipop);

        
        ///////////////////////////
        // create pheno for each iphen
        for (int iphen=0; iphen<nphen; iphen++)
        {
            double var_G_gen0=CommFunc::var(population[ipop].get_bv(iphen));
            double var_a_gen0=CommFunc::var(population[ipop].get_additive(iphen));
            double var_d_gen0=CommFunc::var(population[ipop].get_dominance(iphen));
            population[ipop]._var_bv_gen0[iphen]=var_G_gen0;
            population[ipop]._var_a_gen0[iphen]=var_a_gen0;
            population[ipop]._var_d_gen0[iphen]=var_d_gen0;
            std::cout << "        var(A) before transformation for phenotype " << iphen+1 << " = " << var_a_gen0 << std::endl;
            std::cout << "        var(D) before transformation for phenotype " << iphen+1 << " = " << var_d_gen0 << std::endl;
            std::cout << "        var(G) before transformation for phenotype " << iphen+1 << " = " << var_G_gen0 << std::endl;
            ras_create_pheno(ipop, iphen, population[ipop]._var_a_gen0[iphen], population[ipop]._var_d_gen0[iphen]);
        }
    } // for each pop
    
   
    ///////////////////////////////////////////////////////
    // environmental effects specific to each population
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      environmental effects specific to each population" << std::endl;
    for (int iphen=0; iphen<nphen; iphen++)
    {
        std::cout << "        phenotype: " << iphen+1 << std::endl;
        sim_environmental_effects_specific_to_each_population(iphen);
    }
    
    
    ///////////////////////////////////////////////////////
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      computing mating value and selection value" << std::endl;
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        std::cout << "      Population: " << ipop+1 << std::endl;
        ///////////////////////////////////////////////////////
        // compute_mating_value and selection value
        ras_compute_mating_value_selection_value(0, ipop);
    }
    

    ///////////////////////////////////////////////////////
    //save pheno to class _Pop_info_prev_gen
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        ras_save_human_info_to_Pop_info_prev_gen(ipop);
    }
    
    ///////////////////////////////////////////////////////
    //save the human info to file
    // give some info
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      saving human info" << std::endl;
    int gen_num=0;
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        std::cout << "      Population: " << ipop+1 << std::endl;
        //save the human info to file
        population[ipop].ras_save_human_info(gen_num);
        
        // give some info
        double var_mv = CommFunc::var(population[ipop].get_mating_value());
        double var_sv = CommFunc::var(population[ipop].get_selection_value());
        population[ipop].ret_var_mating_value[gen_num]=var_mv;
        std::cout << "        var(mv)           = " << var_mv << std::endl;
        std::cout << "        var(sv)           = " << var_sv << std::endl;

        int nphen=population[ipop]._pheno_scheme.size();
        for (int iphen=0; iphen<nphen; iphen++)
        {
            double var_phen=CommFunc::var(population[ipop].get_phen(iphen));
            double var_A = CommFunc::var(population[ipop].get_additive(iphen));
            double var_D = CommFunc::var(population[ipop].get_dominance(iphen));
            double var_C = CommFunc::var(population[ipop].get_common(iphen));
            double var_G=CommFunc::var(population[ipop].get_bv(iphen));
            double var_e=CommFunc::var(population[ipop].get_e_noise(iphen));
            double var_parental_effect=CommFunc::var(population[ipop].get_parental_effect(iphen));
            population[ipop].ret_var_phen[iphen][gen_num]=var_phen;
            population[ipop].ret_var_A[iphen][gen_num]=var_A;
            population[ipop].ret_var_D[iphen][gen_num]=var_D;
            population[ipop].ret_var_C[iphen][gen_num]=var_C;
            population[ipop].ret_var_G[iphen][gen_num]=var_G;
            population[ipop].ret_var_E[iphen][gen_num]=var_e;
            population[ipop].ret_h2[iphen][gen_num]=var_G/var_phen;
            population[ipop].ret_var_parental_effect[iphen][gen_num]=var_parental_effect;
            std::cout << "        phenotype: " << iphen+1 << std::endl;
            std::cout << "          var(A)          = " << var_A << std::endl;
            std::cout << "          var(D)          = " << var_D << std::endl;
            std::cout << "          var(G)          = " << var_G << std::endl;
            std::cout << "          var(C)          = " << var_C << std::endl;
            std::cout << "          var(E)          = " << var_e << std::endl;
            std::cout << "          var(F)          = " << var_parental_effect << std::endl;
            std::cout << "          var(P)          = " << var_phen << std::endl;
            std::cout << "          h2              = " << var_G/var_phen << std::endl;
        }
    }

    
    ///////////////////////////////////////////////////////
    //memory info
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      memory used" << std::endl;
    double vm, rss;
    process_mem_usage(vm, rss);
    std::cout << "        VM                = " << vm  << " (Mb)" << std::endl;
    std::cout << "        RSS               = " << rss << " (Mb)" << std::endl;

    if (_debug)
        std::cout << " end of [ras_init]" << std::endl;
    
    return true;
}


//note: std::vector<Human> &h already has phenotype
// and population[ipop].ret_var_phen and population[ipop].ret_var_G are computed
bool Simulation::ras_main_sim(void)
{
    std::cout << " Starting the main function ..." << std::endl;
    
    unsigned seed=ras_now_nanoseconds();
    std::srand(seed);

    for (int gen_num=1; gen_num<=_tot_gen; gen_num++) // for each generation
    {
        std::cout << "    -------------------------------------------------------------------" << std::endl;
        std::cout << "    Start generation " << gen_num << std::endl;
        
        if(!sim_next_generation(gen_num, rand())) return false;
        
        
        std::cout << "    done." << std::endl; // end of one generation
    }
    return true;
}

void Simulation::ras_show_res(void)
{
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        std::cout << " ---------- " << "Population " << ipop+1 << std::endl;
        // phenotypes
        unsigned nphen=population[ipop]._pheno_scheme.size();
        for (unsigned iphen=0; iphen<nphen; iphen++)
        {
            std::cout << " phenotype: " << iphen+1 << std::endl;
            std::cout << "   ret_var_A:";
            for (unsigned i=0; i<population[ipop].ret_var_A[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_A[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_D:";
            for (unsigned i=0; i<population[ipop].ret_var_D[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_D[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_G:";
            for (unsigned i=0; i<population[ipop].ret_var_G[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_G[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_C:";
            for (unsigned i=0; i<population[ipop].ret_var_C[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_C[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_E:";
            for (unsigned i=0; i<population[ipop].ret_var_E[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_E[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_parental_effect:";
            for (unsigned i=0; i<population[ipop].ret_var_parental_effect[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_parental_effect[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_phen:";
            for (unsigned i=0; i<population[ipop].ret_var_phen[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_phen[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_h2:";
            for (unsigned i=0; i<population[ipop].ret_h2[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_h2[iphen][i];
            std::cout << std::endl;
            
            std::cout << "   ret_var_G/ret_var_G[0]:";
            for (unsigned i=0; i<population[ipop].ret_var_G[iphen].size(); i++)
                std::cout << " " << population[ipop].ret_var_G[iphen][i]/population[ipop].ret_var_G[iphen][0];
            std::cout << std::endl;
        }
        // mating_value
        std::cout << " ret_var_mating_value:";
        for (int i=0; i<(int)population[ipop].ret_var_mating_value.size(); i++)
            std::cout << " " << population[ipop].ret_var_mating_value[i];
        std::cout << std::endl;
    }
}


void Simulation::ras_save_res(void)
{
    std::string sep=" ";
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        int nphen=population[ipop]._pheno_scheme.size();

        std::string outfile_name=_out_prefix +".pop"+std::to_string(ipop+1) + ".summary";
        std::ofstream file_summary;
        file_summary.open(outfile_name.c_str());
        
        //create header
        file_summary << "gen" << sep;
        for (int iphen=0; iphen<nphen; iphen++)
        {
            file_summary << "ph" << iphen+1 << "_var_A" << sep;
            file_summary << "ph" << iphen+1 << "_var_D" << sep;
            file_summary << "ph" << iphen+1 << "_var_G" << sep;
            file_summary << "ph" << iphen+1 << "_var_C" << sep;
            file_summary << "ph" << iphen+1 << "_var_E" << sep;
            file_summary << "ph" << iphen+1 << "_var_parental_effect" << sep;
            file_summary << "ph" << iphen+1 << "_var_phen" << sep;
            file_summary << "ph" << iphen+1 << "_h2" << sep;
            file_summary << "ph" << iphen+1 << "_var_G_std" << sep;
        }
        file_summary << "var_mating_value" << std::endl;

        
        // phenotypes
        for (int igen=0; igen<(int)population[ipop].ret_var_G[0].size(); igen++)
        {
            file_summary << igen << sep;
            for (int iphen=0; iphen<nphen; iphen++)
            {
                file_summary << population[ipop].ret_var_A[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_D[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_G[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_C[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_E[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_parental_effect[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_phen[iphen][igen] << sep;
                file_summary << population[ipop].ret_h2[iphen][igen] << sep;
                file_summary << population[ipop].ret_var_G[iphen][igen]/population[ipop].ret_var_G[iphen][0] << sep;
            }
            // mating_value
            file_summary << population[ipop].ret_var_mating_value[igen] << std::endl;
        }
        
        file_summary.close();
    }// for ipop
}


// this file has no header, with n lines equal to the genertions
// each column is the the elements of the transition matrix in format: [a11 a12 a13 ... a21 a22 a23 ...]
bool Simulation::read_migration_file(void)
{
    char sep=' ';
    std::string file_name=par.file_migration;
    std::ifstream ifile(file_name.c_str());
    if(!ifile)
    {
        std::cout << "Error: can not open the migration file [" + file_name + "] to read." << std::endl;
        return false;
    }
    
    int n2=_n_pop*_n_pop;
    std::string line;
    std::string token;
    std::vector<double> temp(n2,0);
    while (getline(ifile, line)){
        std::istringstream iss(line);
        for (int icol=0; icol<n2; icol++)
        {
            std::getline(iss, token, sep);
            if (token=="")
            {
                std::cout << "Error: The file [" + file_name + "] must have n^2 columns, where n is the number of populations." << std::endl;
                return false;
            }
            temp[icol]=std::stod(token);
        }
        migration_mat_gen.push_back(temp);
    }
    if ((int)migration_mat_gen.size() != _tot_gen)
    {
        std::cout << "Error: The file [" + file_name + "] must have "  << _tot_gen << " lines, equal to the number of generations." << std::endl;
        return false;
    }
    return true;
}


bool Simulation::ras_do_migration(int gen_ind)
{
    std::vector<std::vector<double> > migration_mat(_n_pop,std::vector<double>(_n_pop,0));
    int k=0;
    // create transition matrix
    for (int i=0; i<_n_pop; i++)
    {
        double s=0;
        for (int j=0; j<_n_pop; j++)
        {
            migration_mat[i][j]=migration_mat_gen[gen_ind][k];
            k++;
            s+=migration_mat[i][j];
        }
        if (s<0.99999 || s>1.00001)
        {
            std::cout << "Error: The sum of columns in transition matrix in [--file_migration] must be 1." << std::endl;
            return false;
        }
    }
    
    
    // alloc an empty _n_pop*_n_pop Camp
    std::vector<std::vector<Camp> > camp(_n_pop , std::vector<Camp>(_n_pop));
    
    std::vector<std::vector<int> > num_move(_n_pop, std::vector<int>(_n_pop,0));
    // number of migrant inds
    for (int i=0; i<_n_pop; i++)
    {
        for (int j=0; j<_n_pop; j++)
        {
            if (i==j) continue;
            int n_migrant=round(migration_mat[i][j]*(double)population[i].h.size());
            std::cout << "        migration from population " << i+1 << " to " << j+1 << " = " << n_migrant << std::endl;
            num_move[i][j]=n_migrant;
        }
    }
    
    // in order to moev unique inds from i to other pops, we first sample at size equal to the all migrants
    std::vector<std::vector<int> > move_sample_pop(_n_pop);
    for (int i=0; i<_n_pop; i++)
    {
        //std::vector<double> temp(num_move[i].begin(), num_move[i].end());
        int s=CommFunc::sum(num_move[i]);
        std::vector<int> sample = RasRandomNumber::ras_SampleWithoutReplacement((int)population[i].h.size(), s);
        std::sort(sample.begin(),sample.end(),std::greater<int>()); // sort descending
        move_sample_pop[i]=sample;
        int k=0;
        for (int j=0; j<_n_pop; j++)
        {
            if (i==j) continue;
            camp[i][j].humans.resize(num_move[i][j]);
            int it=0;
            while (k<num_move[i][j])
            {
                camp[i][j].humans[it]=population[i].h[sample[k]];
                k++;
                it++;
            }
        }
    }
    
    /*
    for (int i=0; i<_n_pop; i++)
    {
        for (int j=0; j<_n_pop; j++)
        {
            if (i==j) continue;
            std::cout << "        some of the indexes:";
            for (int to=0; to<20; to++) std::cout << " " <<  camp[i][j].humans[to].ID;
            std::cout << std::endl;
        }
    }
    */

     
    for (int i=0; i<_n_pop; i++)
    {
        std::cout << "        size pop "  << i+1 << " before immigration    = " << population[i].h.size() << std::endl;
    }

    
    // remove migrants from the origin population
    for (int i=0; i<_n_pop; i++)
    {
        for (unsigned k=0; k<move_sample_pop[i].size(); k++)
        {
            int ind=move_sample_pop[i][k];
            population[i].h.erase(population[i].h.begin()+ind);
        }
        //std::cout << "        size pop "  << i+1 << " after remove          = " << population[i].h.size() << std::endl;
    }
    
    // send migrants from the camp to the destination population
    for (int i=0; i<_n_pop; i++)
    {
        for (int j=0; j<_n_pop; j++)
        {
            if (i==j) continue;
            for (unsigned k=0; k<camp[i][j].humans.size(); k++)
            {
                population[j].h.push_back(camp[i][j].humans[k]);
            }
        }
    }
    
    for (int i=0; i<_n_pop; i++)
    {
        std::cout << "        size pop "  << i+1 << " after immigration     = " << population[i].h.size() << std::endl;
    }
    
    return true;
}


// this file is for generarting random seed
unsigned Simulation::ras_now_nanoseconds(void)
{
    std::chrono::nanoseconds ns = std::chrono::duration_cast< std::chrono::nanoseconds >(std::chrono::system_clock::now().time_since_epoch());
    return unsigned(ns.count() % 100000000)+std::rand();
}




bool Simulation::ras_save_genotypes(int gen_num)
{
    std::cout << " saving genotypes" << std::endl;
    
    // hap
    if (_format_output=="hap")
    {
        if (!ras_write_hap_legend_sample(gen_num))
        {
            std::cout << "Error in reading and writing haplotypes!" << std::endl;
            return false;
        }
    }
    // plink
    else if (_format_output=="plink")
    {
        if (!ras_write_hap_to_plink_format(gen_num))
        {
            std::cout << "Error in reading and writing haplotypes!" << std::endl;
            return false;
        }
    }
    // interval
    else if (_format_output=="interval")
    {
        if (!ras_write_hap_to_interval_format(gen_num))
        {
            std::cout << "Error in reading and writing haplotypes!" << std::endl;
            return false;
        }
    }
    else
    {
        std::cout << "Error: The output format [" << _format_output << "] is not defined." << std::endl;
        return false;
    }
    
    
    // interval combined by other parameters
    if (_interval)
    {
        if (!ras_write_hap_to_interval_format(gen_num))
        {
            std::cout << "Error in reading and writing haplotypes!" << std::endl;
            return false;
        }
    }

    
    return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// read hap files for each chr and pop

bool Simulation::ras_read_hap_legend_sample_chr(std::vector<Legend> &pops_legend, std::vector<Hap_SNP> &pops_hap, int ichr)
{
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        std::cout << "     --Population: " << (ipop+1) << std::endl;
        // reading sample file;
        // note that .sample file has header file, so -1
        // but .impute.hap.indv has no header!!!!
        //long int nind=CommFunc::ras_FileLineNumber(file_in_name[ichr][2])-1; // sample_file
        long int nind=CommFunc::ras_FileLineNumber(population[ipop]._hap_legend_sample_name[ichr][2]); // .impute.hap.indv
        
        // reading legend file;
        Legend legend;
        long int nsnp=format_hap::read_legend(legend, population[ipop]._hap_legend_sample_name[ichr][1]); //legend file
        std::cout << "       nsnp=" << nsnp << ", nind=" << nind << std::endl;
        
        // reading hap file
        Hap_SNP hap_snp;
        if(!format_hap::read_hap(hap_snp, population[ipop]._hap_legend_sample_name[ichr][0], nind, nsnp, true) ) // hap file
        {
            std::cout << "Error in reading the hap file." << std::endl;
            return false;
        }
        pops_legend[ipop]=legend;
        pops_hap[ipop]=hap_snp;
    }
    return true;
}



///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// impute2: .hap .legend .sample

bool Simulation::ras_write_hap_legend_sample(int gen_num)
{
    int n_chr=population[0].h[0].chr.size();
    for (int ichr=0; ichr<n_chr; ichr++) // for chr
    {
        // reading hap_legend_sample for all populations
        std::cout << "    Start chr " << _all_active_chrs[ichr] << std::endl << std::flush;
        std::cout << "       reading hap files" << std::endl << std::flush;
        std::vector<Legend> pops_legend(_n_pop); // the population's legend file
        std::vector<Hap_SNP> pops_hap(_n_pop); // the population's hap file
        ras_read_hap_legend_sample_chr(pops_legend, pops_hap, ichr); // for all populations
        std::cout << "       done." << std::endl << std::flush;
        
        
        //convert hap matrix to matrix_hap according to human interval information
        for (int ipop=0; ipop<_n_pop; ipop++) // for pops
        {
            std::cout << "    --Population: " << (ipop+1) << std::endl;
            
            std::cout << "      converting to hap file matrix" << std::endl << std::flush;
            Hap_SNP hap_snp;
            ras_convert_interval_to_hap_matrix(ipop, pops_hap, pops_legend, ichr, hap_snp);
            
            
            // writing the hap
            std::cout << "      writing" << std::endl << std::flush;
            std::string outfile_name=_out_prefix +".pop"+std::to_string(ipop+1) + ".gen" + std::to_string(gen_num)+".chr" + std::to_string(_all_active_chrs[ichr]);
            format_hap::write_hap(hap_snp, outfile_name);
        }
        std::cout << "    --------------------------------------------------------------" << std::endl;

    }
    return true;
}


// return: hap_snp
bool Simulation::ras_convert_interval_to_hap_matrix(int ipop, std::vector<Hap_SNP> &pops_hap, std::vector<Legend> &pops_legend, int ichr, Hap_SNP &hap_snp)
{
    // convert hap0 matrix to matrix_hap according to Human interval information
    unsigned long int n_human=population[ipop].h.size();
    unsigned long int nsnp=pops_legend[ipop].id.size();
    std::cout << "      n_human=" << n_human << std::endl;
    int n_chromatid=population[ipop].h[0].chr[ichr].Hap.size(); // =2
    
    std::cout << "      Allocating memory ..." << std::flush;
    hap_snp.hap.resize(n_human*2, std::vector<bool> (nsnp) );
    std::cout << "      done." << std::endl;
    
    
    for (unsigned long int ih=0; ih<n_human; ih++) // for humans
    {
        for (int ihaps=0; ihaps<n_chromatid; ihaps++) // for haps
        {
            int n_parts=(int)population[ipop].h[ih].chr[ichr].Hap[ihaps].size();
            for (int ip=0; ip<n_parts; ip++) // for parts
            {
                part p=population[ipop].h[ih].chr[ichr].Hap[ihaps][ip];
                int root_pop=p.root_population;
                for (unsigned long int ii=0; ii<pops_legend[root_pop].pos.size(); ii++)
                {
                    if (p.check_interval(pops_legend[root_pop].pos[ii]))
                    {
                        if (p.hap_index< 0 || p.hap_index >= pops_hap[root_pop].hap.size())
                        {
                            std::cout << "Error: p.hap_index=" << p.hap_index << " is not in range, ih=" << ih << std::endl;
                            return false;
                        }
                        // check mutation
                        std::vector<unsigned long int>::iterator it = find(p.mutation_pos.begin(), p.mutation_pos.end(), pops_legend[root_pop].pos[ii]);
                        if (it != p.mutation_pos.end()) // mutation
                            hap_snp.hap[2*ih+ihaps][ii] = !pops_hap[root_pop].hap[p.hap_index][ii];
                        else // no mutation
                            hap_snp.hap[2*ih+ihaps][ii] = pops_hap[root_pop].hap[p.hap_index][ii];

                    }
                    // for the last part and last snps
                    if (ii==pops_legend[root_pop].pos.size()-1 && ip==(int)population[ipop].h[ih].chr[ichr].Hap[ihaps].size()-1)
                    {
                        hap_snp.hap[2*ih+ihaps][ii]=pops_hap[root_pop].hap[p.hap_index][ii];
                    }
                }
            }
        }
    }
    return true;
}





///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// plink

bool Simulation::ras_write_hap_to_plink_format(int gen_num)
{
    int n_chr=population[0].h[0].chr.size();
    for (int ichr=0; ichr<n_chr; ichr++) // for chr
    {
        // reading hap_legend_sample for all populations
        std::cout << "    Start chr " << _all_active_chrs[ichr] << std::endl << std::flush;
        std::cout << "       reading hap files" << std::endl << std::flush;
        std::vector<Legend> pops_legend(_n_pop); // the population's legend file
        std::vector<Hap_SNP> pops_hap(_n_pop); // the population's hap file
        ras_read_hap_legend_sample_chr(pops_legend, pops_hap, ichr); // for all populations
        std::cout << "       done." << std::endl << std::flush;
        
        //convert hap matrix to matrix_hap according to human interval information
        for (int ipop=0; ipop<_n_pop; ipop++) // for pops
        {
            std::cout << "    --Population: " << (ipop+1) << std::endl;

            std::cout << "      converting to plink file matrix" << std::endl << std::flush;
            std::vector<std::vector<bool> > matrix_plink_ped;
            plink_PED_ids plink_ped_ids;
            plink_MAP plink_map;
            ras_convert_interval_to_format_plink(ipop, pops_hap, pops_legend, ichr, matrix_plink_ped, plink_ped_ids, plink_map);
            
            
            // writing the hap
            std::cout << "      writing" << std::endl << std::flush;
            std::string outfile_name=_out_prefix +".pop"+ std::to_string(ipop+1) + ".gen" + std::to_string(gen_num)+".chr"+ std::to_string(_all_active_chrs[ichr]);
            format_plink::write_ped_map(outfile_name, matrix_plink_ped, plink_ped_ids, plink_map);
        }
        std::cout << "    --------------------------------------------------------------" << std::endl;
    }
    return true;
}


//inputs: ipop, pops_hap, pops_legend, ichr
//output: matrix_plink_ped, plink_ped_ids, plink_map
bool Simulation::ras_convert_interval_to_format_plink(int ipop, std::vector<Hap_SNP> &pops_hap, std::vector<Legend> &pops_legend, int ichr, std::vector<std::vector<bool> > &matrix_plink_ped, plink_PED_ids &plink_ped_ids, plink_MAP &plink_map)
{
    // convert hap0 matrix to matrix_hap according to Human interval information
    
    unsigned long int n_human=population[ipop].h.size();
    unsigned long int nsnp=pops_legend[ipop].id.size();
    std::cout << "      n_human=" << n_human << std::endl;
    int n_chromatid=population[ipop].h[0].chr[ichr].Hap.size(); // =2
    
    
    std::cout << "      Allocating memory ..." << std::flush;
    try
    {
        matrix_plink_ped.clear();
        
        matrix_plink_ped.resize(n_human, std::vector<bool> (nsnp*n_chromatid, false) );
        plink_ped_ids.alloc(n_human); // FID IID PID MID sex phen
        plink_map.alloc(nsnp);
        std::cout << "      done." << std::endl;
    }
    catch (std::exception& e)
    {
        std::cout << "Standard exception: " << e.what() << std::endl;
        return false;
    }
    
    
    
    for (unsigned long int ih=0; ih<n_human; ih++) // for humans
    {
        //if(ih % 1000==0) cout << "ih=" << ih << endl << flush;
        for (int ihaps=0; ihaps<n_chromatid; ihaps++) // for 2 haps
        {
            int n_parts=(int)population[ipop].h[ih].chr[ichr].Hap[ihaps].size();
            for (int ip=0; ip<n_parts; ip++) // for parts
            {
                part p=population[ipop].h[ih].chr[ichr].Hap[ihaps][ip];
                int root_pop=p.root_population;
                for (unsigned long int ii=0; ii<nsnp; ii++)
                {
                    if (p.check_interval(pops_legend[root_pop].pos[ii]))
                    {
                        if (p.hap_index< 0 || p.hap_index >= pops_hap[root_pop].hap.size())
                        {
                            std::cout << "Error: p.hap_index=" << p.hap_index << " is not in range, ih=" << ih << std::endl;
                            return false;
                        }
                        // check mutation
                        std::vector<unsigned long int>::iterator it = find(p.mutation_pos.begin(), p.mutation_pos.end(), pops_legend[root_pop].pos[ii]);
                        if (it != p.mutation_pos.end()) // mutation
                            matrix_plink_ped[ih][2*ii+ihaps]= !pops_hap[root_pop].hap[p.hap_index][ii];
                        else // no mutation
                            matrix_plink_ped[ih][2*ii+ihaps]= pops_hap[root_pop].hap[p.hap_index][ii];
                    }
                    // for the last part and last snps
                    if (ii==pops_legend[root_pop].pos.size()-1 && ip==(int)population[ipop].h[ih].chr[ichr].Hap[ihaps].size()-1)
                    {
                        matrix_plink_ped[ih][2*ii+ihaps]=pops_hap[root_pop].hap[p.hap_index][ii];
                    }
                }
            }
        }
    }
    
    
    if (population[ipop]._debug)
    {
        //check MAF for the last SNPs
        std::cout << "The last allele frequencies" << std::endl;
        for (unsigned long int k=10; k>0; k--)
        {
            std::vector<double> temp(2*n_human);
            unsigned long int ii_ind=nsnp-k;
            for (unsigned long int ih=0; ih<n_human; ih++) // for humans
            {
                for (int ihaps=0; ihaps<n_chromatid; ihaps++) // for 2 haps
                {
                    temp[2*ih+ihaps]=(double)matrix_plink_ped[ih][2*ii_ind+ihaps];
                }
            }
            std::cout << "AF = " << CommFunc::mean(temp) << std::endl;
        }
    }
    
    
    // fill plink_PED_ids class
    for (unsigned long int ih=0; ih<n_human; ih++) // for humans
    {
        // FID IID PID MID
        // all IDs are +1, since we want IDs>=1 (in c++, indexes are zero-based)
        // in plink zero means missing value
        plink_ped_ids.FID[ih] = std::to_string(population[ipop].h[ih].ID_Father+1); // letting ID_Father as FID
        plink_ped_ids.IID[ih] = std::to_string(population[ipop].h[ih].ID+1);
        plink_ped_ids.PID[ih] = std::to_string(population[ipop].h[ih].ID_Father+1);
        plink_ped_ids.MID[ih] = std::to_string(population[ipop].h[ih].ID_Mother+1);
        plink_ped_ids.sex[ih] = population[ipop].h[ih].sex; // 1=male, 2=female
        plink_ped_ids.phen[ih] = -9; // we can also let them -9
    }
    
    // fill plink_MAP class
    for (unsigned long int i=0; i<nsnp; i++) // for SNPs
    {
        plink_map.chr[i]=std::to_string(_all_active_chrs[ichr]);
        plink_map.rs[i]=pops_legend[ipop].id[i];
        plink_map.cM[i]=0;
        plink_map.pos[i]=pops_legend[ipop].pos[i];
        plink_map.al0[i]=pops_legend[ipop].al0[i];
        plink_map.al1[i]=pops_legend[ipop].al1[i];
    }
    
    return true;
}



///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// interval file format (for checking the IBDs)

// the ID is zero-based
bool Simulation::ras_write_hap_to_interval_format(int gen_num)
{
    std::string sep=" ";
    int n_chr=population[0].h[0].chr.size();
    
    
    for (int ipop=0; ipop<_n_pop; ipop++) // for pops
    {
        std::cout << "      --Population: " << (ipop+1) << std::endl;
        unsigned long int nind = population[ipop].h.size();

        for (int ichr=0; ichr<n_chr; ichr++) // for chr
        {
            std::cout << "      Start chr " << _all_active_chrs[ichr] << std::endl << std::flush;
            std::string outfile_name=_out_prefix +".pop"+ std::to_string(ipop+1) + ".gen" + std::to_string(gen_num)+".chr"+ std::to_string(_all_active_chrs[ichr])+".int";
            std::ofstream file_out;
            file_out.open(outfile_name.c_str());
            
            //header
            file_out << "h_ID" << sep;
            file_out << "chr" << sep;
            file_out << "hap" << sep;
            file_out << "st" << sep;
            file_out << "en" << sep;
            file_out << "hap_index" << sep;
            file_out << "root_pop" << std::endl;


            for (unsigned long int ih=0; ih<nind; ih++)
            {
                for (int ihap=0; ihap<2; ihap++)
                {
                    int npart=population[ipop].h[ih].chr[ichr].Hap[ihap].size();
                    for (int ip=0; ip<npart; ip++)
                    {
                        part p=population[ipop].h[ih].chr[ichr].Hap[ihap][ip];
                        file_out << population[ipop].h[ih].ID << sep;
                        file_out << _all_active_chrs[ichr] << sep;
                        file_out << ihap << sep;
                        file_out << p.st << sep;
                        file_out << p.en << sep;
                        file_out << p.hap_index << sep;
                        file_out << p.root_population << std::endl;
                    }// for parts
                }// for ihap

            }// for ih
            
            file_out.close();
        }// for chr
    
    }// for pops

    return true;
}


// main body for each generation
bool Simulation::sim_next_generation(int gen_num, unsigned seed)
{
    
    int time_start_this_generation = time(0);
    int nphen=population[0]._pheno_scheme.size();
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        std::cout << "      Population: " << ipop+1 << std::endl;
        std::cout << "      user specified parameters" << std::endl;
        std::cout << "        pop_size          = " << population[ipop]._pop_size[gen_num-1] << std::endl;
        std::cout << "        mat_cor           = " << population[ipop]._mat_cor[gen_num-1] << std::endl;
        std::cout << "        offspring_dist    = " << population[ipop]._offspring_dist[gen_num-1] << std::endl;
        std::cout << "        selection_func    = " << population[ipop]._selection_func[gen_num-1] << " " << population[ipop]._selection_func_par1[gen_num-1] << " " << population[ipop]._selection_func_par2[gen_num-1] << std::endl;

        
        ///////////////////////////////////////////////////////
        //mate the people
        if (population[ipop]._RM)
        {
            std::cout << "      random mating" << std::endl;
            if(!random_mate(ipop, gen_num-1, seed)) // gen_num-1, because it starts from 1
                return false;
        }
        else
        {
            std::cout << "      assortative mating" << std::endl;
            if(!assort_mate(ipop, gen_num-1)) // gen_num-1, because it starts from 1
                return false;
        }
        
        // give some info
        double couple_cor_mv=compute_couple_cor_mating_value(ipop);
        std::cout << "        couples_info.size = " << population[ipop]._couples_info.size() << std::endl;
        std::cout << "        couple_cor_mv     = " << couple_cor_mv << std::endl;
        std::cout << std::flush;
        
        ///////////////////////////////////////////////////////
        //A8 Reproduce: create new humans
        std::cout << "      reproducing" << std::endl;
        population[ipop].h=reproduce(ipop, gen_num);
        std::cout << "        human size        = " << population[ipop].h.size() << std::endl;
        
        ///////////////////////////////////////////////////////
        //computing breeding value
        std::cout << "      computing breeding value" << std::endl;
        if(!ras_compute_breeding_value(ipop) ) return false;

        
        ///////////////////////////////////////////////////////
        //A6 Create Phenotypes (e_noise, parental_effect, phen)
        std::cout << "      creating phenotypes" << std::endl;
        // create pheno for each iphen
        //unsigned long int n=population[ipop].h.size();
        for (int iphen=0; iphen<nphen; iphen++)
        {
            if(!ras_create_pheno(ipop, iphen, population[ipop]._var_a_gen0[iphen], population[ipop]._var_d_gen0[iphen]))
                return false;
            // give some info
            double var_phen=CommFunc::var(population[ipop].get_phen(iphen));
            double var_A = CommFunc::var(population[ipop].get_additive(iphen));
            double var_D = CommFunc::var(population[ipop].get_dominance(iphen));
            double var_C = CommFunc::var(population[ipop].get_common(iphen));
            double var_G=CommFunc::var(population[ipop].get_bv(iphen));
            double var_e=CommFunc::var(population[ipop].get_e_noise(iphen));
            double var_parental_effect=CommFunc::var(population[ipop].get_parental_effect(iphen));
            std::cout << "        phenotype: " << iphen+1 << std::endl;
            std::cout << "          var(A)          = " << var_A << std::endl;
            std::cout << "          var(D)          = " << var_D << std::endl;
            std::cout << "          var(G)          = " << var_G << std::endl;
            std::cout << "          var(C)          = " << var_C << std::endl;
            std::cout << "          var(E)          = " << var_e << std::endl;
            std::cout << "          var(F)          = " << var_parental_effect << std::endl;
            std::cout << "          var(P)          = " << var_phen << std::endl;
            std::cout << "          h2              = " << var_G/var_phen << std::endl;
        }
        
    } // for each population
    
    
    ///////////////////////////////////////////////////////
    // environmental effects specific to each population
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      environmental effects specific to each population" << std::endl;
    for (int iphen=0; iphen<nphen; iphen++)
    {
        std::cout << "        phenotype: " << iphen+1 << std::endl;
        sim_environmental_effects_specific_to_each_population(iphen);
    }
    
    
    ///////////////////////////////////////////////////////
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      computing mating value and selection value" << std::endl;
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        std::cout << "      Population: " << ipop+1 << std::endl;
        ///////////////////////////////////////////////////////
        // compute_mating_value and selection value
        ras_compute_mating_value_selection_value(gen_num, ipop);
    }
    
    
    ///////////////////////////////////////////////////////
    //migration
    if (_n_pop>1)
    {
        std::cout << "      -------------------------" << std::endl;
        std::cout << "      migration" << std::endl;
        if(!ras_do_migration(gen_num-1)) return false;
    }
    
    ///////////////////////////////////////////////////////
    //save pheno to class _Pop_info_prev_gen
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        ras_save_human_info_to_Pop_info_prev_gen(ipop);
    }
    
    ///////////////////////////////////////////////////////
    //save the human info into the file
    // give some info
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      saving human info" << std::endl;
    
    for (int ipop=0; ipop<_n_pop; ipop++) // for each population
    {
        std::cout << "      Population: " << ipop+1 << std::endl;
        //save the human info to file
        population[ipop].ras_save_human_info(gen_num);
        
        // give some info
        double var_mv = CommFunc::var(population[ipop].get_mating_value());
        population[ipop].ret_var_mating_value[gen_num]=var_mv;
        std::cout << "        var(mv)           = " << var_mv << std::endl;
        std::cout << "        var(sv)           = " << CommFunc::var(population[ipop].get_selection_value()) << std::endl;
        
        for (int iphen=0; iphen<nphen; iphen++)
        {
            double var_phen = CommFunc::var(population[ipop].get_phen(iphen));
            double var_A = CommFunc::var(population[ipop].get_additive(iphen));
            double var_D = CommFunc::var(population[ipop].get_dominance(iphen));
            double var_C = CommFunc::var(population[ipop].get_common(iphen));
            double var_G = CommFunc::var(population[ipop].get_bv(iphen));
            double var_e = CommFunc::var(population[ipop].get_e_noise(iphen));
            double var_parental_effect = CommFunc::var(population[ipop].get_parental_effect(iphen));
            population[ipop].ret_var_phen[iphen][gen_num]=var_phen;
            population[ipop].ret_var_A[iphen][gen_num]=var_A;
            population[ipop].ret_var_D[iphen][gen_num]=var_D;
            population[ipop].ret_var_C[iphen][gen_num]=var_C;
            population[ipop].ret_var_G[iphen][gen_num]=var_G;
            population[ipop].ret_var_E[iphen][gen_num]=var_e;
            population[ipop].ret_h2[iphen][gen_num]=var_G/var_phen;
            population[ipop].ret_var_parental_effect[iphen][gen_num]=var_parental_effect;
            std::cout << "        phenotype: " << iphen+1 << std::endl;
            std::cout << "          var(A)          = " << var_A << std::endl;
            std::cout << "          var(D)          = " << var_D << std::endl;
            std::cout << "          var(G)          = " << var_G << std::endl;
            std::cout << "          var(C)          = " << var_C << std::endl;
            std::cout << "          var(E)          = " << var_e << std::endl;
            std::cout << "          var(F)          = " << var_parental_effect << std::endl;
            std::cout << "          var(P)          = " << var_phen << std::endl;
            std::cout << "          h2              = " << var_G/var_phen << std::endl;
        }
    }
    
    ///////////////////////////////////////////////////////
    //save genotypes
    if(_output_all_generations) // for all generations
    {
        ras_save_genotypes(gen_num);
    }
    

    ///////////////////////////////////////////////////////
    //memory info
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      memory used" << std::endl;
    double vm, rss;
    process_mem_usage(vm, rss);
    std::cout << "        VM                = " << vm  << " (Mb)" << std::endl;
    std::cout << "        RSS               = " << rss << " (Mb)" << std::endl;
    
    int time_this_generation = time(0)-time_start_this_generation;
    ///////////////////////////////////////////////////////
    //memory info
    std::cout << "      -------------------------" << std::endl;
    std::cout << "      time used for this generation: " << time_this_generation << " seconds" << std::endl;

    return true;
}



// _offspring_dist has no meaning here
// _avoid_inbreeding has no meaning here
// _mat_cor has no meaning here

bool Simulation::random_mate(int ipop, int gen_ind, unsigned seed)
{
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    
    unsigned long int n_h =(int)population[ipop].h.size();
    
    // MARRIAGEABLE PEOPLE - this section creates an equal number of males and females who will be paired off below
    
    // check selection_value
    std::vector<unsigned long int> pos_male;
    std::vector<unsigned long int> pos_female;
    for (unsigned long int i=0; i<n_h; i++)
    {
        // check selection_value
        double r=distribution(generator);
        if (r < population[ipop].h[i].selection_value)
        {
            if (population[ipop].h[i].sex==1) pos_male.push_back(i); // this male can mary
            else if (population[ipop].h[i].sex==2) pos_female.push_back(i); //this female can mary
        }
    }
    unsigned long int num_males_mate=pos_male.size();
    unsigned long int num_females_mate=pos_female.size();
    
    std::cout << "        num_males_mate    = " << num_males_mate << std::endl;
    std::cout << "        num_females_mate  = " << num_females_mate << std::endl;
    
    if (num_males_mate==0 || num_females_mate==0)
    {
        std::cout << "Error: No one can marry, num_males_mate=" << num_males_mate << ", num_females_mate=" << num_females_mate << std::endl;
        return false;
    }
    

    std::default_random_engine g_uint_f(seed+1);
    std::default_random_engine g_uint_m(seed+2);
    std::uniform_int_distribution<unsigned long int> d_uint_f(0,pos_male.size()-1);
    std::uniform_int_distribution<unsigned long int> d_uint_m(0,pos_female.size()-1);
    
    //marry the males and females according to the assortative mating coefficient earlier inputed
    std::vector<Couples_Info> couples_info(population[ipop]._pop_size[gen_ind]);

    // randomly select father and mother form pos_male_marriageable and pos_female_marriageable
    for (unsigned long int i=0; i< population[ipop]._pop_size[gen_ind]; i++)
    {
        //generate random index for mother and father
        unsigned long int i_f=d_uint_f(g_uint_f);
        unsigned long int i_m=d_uint_m(g_uint_m);
        couples_info[i].pos_male=pos_male[i_f];
        couples_info[i].pos_female=pos_female[i_m];
        couples_info[i].inbreed=false; // no inbreed, so can marry
        couples_info[i].num_offspring=1;
    }
    

    population[ipop]._couples_info=couples_info;
    
    return true;
    
}






// mate based on the phenotype
// it is possible taht an individual has 2 spouses, based on _MM_percent
bool Simulation::assort_mate(int ipop, int gen_ind)
{
    //std::srand(unsigned(std::time(0)));
    unsigned seed=ras_now_nanoseconds();
    std::srand(seed);
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    unsigned long int n_h =(int)population[ipop].h.size();
    
    // MARRIAGEABLE PEOPLE - this section creates an equal number of males and females who will be paired off below
    std::vector<Ind_MatingValue> pos_male_marriageable;
    std::vector<Ind_MatingValue> pos_female_marriageable;
    for (unsigned long int i=0; i<n_h; i++)
    {
        //double r=(double)std::rand()/RAND_MAX;
        double r=distribution(generator);
        if (r < population[ipop].h[i].selection_value)
        {
            if (population[ipop].h[i].sex==1) //this male can marry
            {
                Ind_MatingValue ind_mv;
                ind_mv.ind=i;
                ind_mv.mating_value=population[ipop].h[i].mating_value;
                pos_male_marriageable.push_back(ind_mv);
                // cjeck if this ind has 2 spouses
                //double r=(double)std::rand()/RAND_MAX;
                double r=distribution(generator);
                if (r < population[ipop]._MM_percent)
                    pos_male_marriageable.push_back(ind_mv);
            }
            else if (population[ipop].h[i].sex==2) //this female can marry
            {
                Ind_MatingValue ind_mv;
                ind_mv.ind=i;
                ind_mv.mating_value=population[ipop].h[i].mating_value;
                pos_female_marriageable.push_back(ind_mv);
                // cjeck if this ind has 2 spouses
                //double r=(double)std::rand()/RAND_MAX;
                double r=distribution(generator);
                if (r < population[ipop]._MM_percent)
                    pos_female_marriageable.push_back(ind_mv);
            }
        }
    }
    unsigned long int num_males_mate=pos_male_marriageable.size();
    unsigned long int num_females_mate=pos_female_marriageable.size();
    
    unsigned long int couples=std::min(num_males_mate,num_females_mate);
    
    std::cout << "        num_males_mate    = " << num_males_mate << std::endl;
    std::cout << "        num_females_mate  = " << num_females_mate << std::endl;
    std::cout << "        couples           = " << couples << std::endl;

    if (couples==0)
    {
        std::cout << "Error: couples=0, num_males_mate=" << num_males_mate << ", num_females_mate=" << num_females_mate << std::endl;
        return false;
    }
    
    //order the males & females by their mating phenotypic value, lowest to highest
    std::sort(pos_male_marriageable.begin(),pos_male_marriageable.end(),ValueCmp);
    std::sort(pos_female_marriageable.begin(),pos_female_marriageable.end(),ValueCmp);
    
    //making the template correlation matrix and ordering the two variables
    //as it is (empirical=TRUE), the correlation is *exactly* the AM coefficient each generation; may change to FALSE to make it more realistic
    //nevertheless, there is a stochastic element to it due to the sorting
    std::vector<double> mu(2,0);
    std::vector<std::vector<double> > corr(2,std::vector<double>(2,0));
    corr[0][0]=1;
    corr[0][1]=population[ipop]._mat_cor[gen_ind];
    corr[1][0]=population[ipop]._mat_cor[gen_ind];
    corr[1][1]=1;
    
    unsigned long int n_couples2=std::min(pos_male_marriageable.size(), pos_female_marriageable.size());
    std::vector<std::vector<double> > template_AM_dist = RasRandomNumber::ras_mvnorm(n_couples2,mu,corr);
    
    std::vector<double> t1(template_AM_dist.size());
    for (int i=0; i<(int)template_AM_dist.size(); i++) t1[i]=template_AM_dist[i][0];
    std::vector<double> t2(template_AM_dist.size());
    for (int i=0; i<(int)template_AM_dist.size(); i++) t2[i]=template_AM_dist[i][1];
    
    std::cout << "        cor(mates)        = " << CommFunc::cor(t1,t2) << std::endl;
    
    std::vector<unsigned long int> rank_template_males = CommFunc::ras_rank(t1);
    std::vector<unsigned long int> rank_template_females = CommFunc::ras_rank(t2);
    //std::cout << "debug: rank" << std::endl;
    
    // to check the corr rank
    // convert to double
    //std::vector<double> v1(rank_template_males.begin(), rank_template_males.end());
    //std::vector<double> v2(rank_template_females.begin(), rank_template_females.end());
    //std::cout << "        cor(rank(t1),rank(t2))=" << CommFunc::cor(v1,v2) << std::endl;
    
    
    
    //marry the males and females according to the assortative mating coefficient earlier inputed
    std::vector<Couples_Info> couples_info(n_couples2);
    
    unsigned long int n_inbreed=0;
    for (unsigned long int i=0; i<n_couples2; i++)
    {
        //if (i%1000==0)
        //    std::cout << "\r        couples: " << i << std::flush;
        
        unsigned long int ind_male=rank_template_males[i];
        unsigned long int pos_m=pos_male_marriageable[ind_male].ind;
        couples_info[i].pos_male=pos_m;
        unsigned long int ind_female=rank_template_females[i];
        unsigned long int pos_f=pos_female_marriageable[ind_female].ind;
        couples_info[i].pos_female=pos_f;
        
        //AVOID INBREEDING
        if (population[ipop]._avoid_inbreeding)
        {
            //for now, we just drop the inbreeding pairs; an alternative is to remate them
            bool inbreeding_sib = (population[ipop].h[pos_m].ID_Father==population[ipop].h[pos_f].ID_Father);
            bool inbreeding_cousin=(population[ipop].h[pos_m].ID_Fathers_Father==population[ipop].h[pos_f].ID_Fathers_Father ||
                                    population[ipop].h[pos_m].ID_Fathers_Father==population[ipop].h[pos_f].ID_Mothers_Father ||
                                    population[ipop].h[pos_m].ID_Mothers_Father==population[ipop].h[pos_f].ID_Fathers_Father ||
                                    population[ipop].h[pos_m].ID_Mothers_Father==population[ipop].h[pos_f].ID_Mothers_Father ||
                                    population[ipop].h[pos_m].ID_Fathers_Mother==population[ipop].h[pos_f].ID_Fathers_Mother ||
                                    population[ipop].h[pos_m].ID_Fathers_Mother==population[ipop].h[pos_f].ID_Mothers_Mother ||
                                    population[ipop].h[pos_m].ID_Mothers_Mother==population[ipop].h[pos_f].ID_Fathers_Mother ||
                                    population[ipop].h[pos_m].ID_Mothers_Mother==population[ipop].h[pos_f].ID_Mothers_Mother);
            
            bool inbreeding =(inbreeding_sib || inbreeding_cousin);
            couples_info[i].inbreed=inbreeding;
            if (inbreeding) n_inbreed++;
        }
        else
            couples_info[i].inbreed=false; // can marry
    }
    
    //num_offspring distribution
    if (population[ipop]._offspring_dist[gen_ind]=="p" || population[ipop]._offspring_dist[gen_ind]=="P") //Poisson
    {
        double lam=(double)population[ipop]._pop_size[gen_ind]/(n_couples2-n_inbreed);
        std::vector<int> num_offspring = RasRandomNumber::ras_rpois(n_couples2,lam);
        for (unsigned long int i=0; i<n_couples2; i++)
        {
            couples_info[i].num_offspring=num_offspring[i]; // we want just one Poisson rv
        }
    }
    else if (population[ipop]._offspring_dist[gen_ind]=="f" || population[ipop]._offspring_dist[gen_ind]=="F") // fixed distribution
    {
        int nf=round((double)population[ipop]._pop_size[gen_ind]/(n_couples2-n_inbreed));
        for (unsigned long int i=0; i<n_couples2; i++)
        {
            couples_info[i].num_offspring=nf;
        }
    }
    
    population[ipop]._couples_info=couples_info;
    
    return true;
}





bool Simulation::ras_allocate_memory_for_humans(std::vector<Human> &h, unsigned long int  n_people, int nchr, int nphen)
{
    std::cout << "        Allocating memory: " << std::flush;
    h.resize(n_people);
    for(unsigned long int it=0; it<n_people; it++)
    {
        h[it].chr.resize(nchr);
        for (int ichr=0; ichr<nchr; ichr++)
        {
            h[it].chr[ichr].Hap.resize(2);
            //h[it].chr[ichr].chromatid.resize(2);
            h[it].chr[ichr].bv_chr.resize(nphen);
            h[it].chr[ichr].additive_chr.resize(nphen);
            h[it].chr[ichr].dominance_chr.resize(nphen);
        }
        h[it].additive.resize(nphen);
        h[it].dominance.resize(nphen);
        h[it].common_sibling.resize(nphen);
        h[it].bv.resize(nphen);
        h[it].e_noise.resize(nphen);
        h[it].parental_effect.resize(nphen);
        h[it].phen.resize(nphen);
    }
    std::cout << "done." << std::endl;
    
    return true;
}

std::vector<Human> Simulation::reproduce(int ipop, int gen_num)
{
    //std::srand(unsigned(std::time(0)));
    unsigned seed=ras_now_nanoseconds();
    std::srand(seed);
    
    unsigned long int n_couples=population[ipop]._couples_info.size();
    unsigned long int n_people = 0;
    for (unsigned long int i=0; i<n_couples; i++)
        if (!population[ipop]._couples_info[i].inbreed)
            n_people+=population[ipop]._couples_info[i].num_offspring;
    
    int nchr=(int)population[ipop].h[0].chr.size();
    int nphen=population[ipop]._pheno_scheme.size();
    
    // alloc memory for humans
    std::vector<Human> h_ret;
    ras_allocate_memory_for_humans(h_ret, n_people, nchr, nphen);
    
    
    // create random normal for common effect
    std::default_random_engine generator(seed);
    std::vector<std::vector <double> > val_common(nphen,std::vector <double>(n_couples,0)); // default is zero
    for (int iphen=0; iphen<nphen; iphen++) // for pheno
    {
        if(population[ipop]._pheno_scheme[iphen]._vc>0)
        {
            std::normal_distribution<double> distribution(0.0, sqrt(population[ipop]._pheno_scheme[iphen]._vc));
            for(unsigned long int it=0; it<n_couples; it++) // for families
            {
                val_common[iphen][it]=distribution(generator);
            }
        }
    }

    // main loop
    unsigned long int i_people=0;
    for(unsigned long int it=0; it<n_couples; it++) // for each new born ind
    {
        if (it%1000==0)
            std::cout << "\r        couples: " << it << std::flush;

        unsigned long int pos_male=population[ipop]._couples_info[it].pos_male;
        unsigned long int pos_female=population[ipop]._couples_info[it].pos_female;
        if (!population[ipop]._couples_info[it].inbreed){
            Human h_pat=population[ipop].h[pos_male];
            Human h_mat=population[ipop].h[pos_female];
            for (int ns=0; ns<population[ipop]._couples_info[it].num_offspring; ns++){ // for each offspring
                for (int ichr=0; ichr<nchr; ichr++) // for each chromosome
                {
                    //get new paternal haplotype
                    unsigned seed_loc=std::rand();
                    std::vector<unsigned long int> rec_bp_loc_pat=ras_sim_loc_rec(population[ipop]._recom_prob[ichr], population[ipop]._rmap[ichr], seed_loc);
                    int starting_haplotype_pat=(std::rand() % 2);
                    std::vector<part> hap_rec_pat=recombine(h_pat.chr[ichr], starting_haplotype_pat, rec_bp_loc_pat);
                    
                    //get new maternal haplotype
                    std::vector<unsigned long int> rec_bp_loc_mat=ras_sim_loc_rec(population[ipop]._recom_prob[ichr], population[ipop]._rmap[ichr], seed_loc);
                    int starting_haplotype_mat=(std::rand() % 2);
                    std::vector<part> hap_rec_mat=recombine(h_mat.chr[ichr], starting_haplotype_mat, rec_bp_loc_mat);
                    
                    //add mutation
                    if(population[ipop]._mutation_map.size()>0)
                    {
                        ras_add_mutation(ipop, ichr, hap_rec_pat, hap_rec_mat);
                    }
                    
                    
                    // add pat and mat to chromosome
                    h_ret[i_people].chr[ichr].Hap[0]=hap_rec_pat;
                    h_ret[i_people].chr[ichr].Hap[1]=hap_rec_mat;
                }// for ichr
                
                // add more info
                h_ret[i_people].gen_num=gen_num;
                h_ret[i_people].sex=(std::rand() % 2)+1; // 1=male, 2=female
                h_ret[i_people].ID=i_people;
                h_ret[i_people].ID_Father=h_pat.ID;
                h_ret[i_people].ID_Fathers_Mother=h_pat.ID_Mother;
                h_ret[i_people].ID_Fathers_Father=h_pat.ID_Father;
                h_ret[i_people].ID_Mother=h_mat.ID;
                h_ret[i_people].ID_Mothers_Mother=h_mat.ID_Mother;
                h_ret[i_people].ID_Mothers_Father=h_mat.ID_Father;
                // common sibling effect
                for (int iphen=0; iphen<nphen; iphen++) // for pheno
                {
                    h_ret[i_people].common_sibling[iphen]=val_common[iphen][it];
                }
                i_people++;
            } // for ns
        }// if
    }// for it
    std::cout << "\r        couples: " << n_couples << std::flush;
    std::cout << std::endl;

    return h_ret;
}


// generate some mutation sites and add them to hap[0] or hap[1]
bool Simulation::ras_add_mutation(int ipop, int ichr, std::vector<part> &hap_rec_pat, std::vector<part> &hap_rec_mat)
{
    unsigned seed=ras_now_nanoseconds();
    std::srand(unsigned(seed));
    
    std::default_random_engine generator(seed);
    std::default_random_engine generator_u(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    
    // it shpould start from 1 not 0
    for (unsigned int i=1; i<population[ipop]._mutation_map[ichr].mutation_rate.size(); i++)
    {
        //double r=((double)std::rand())/RAND_MAX;
        double r=distribution(generator_u);
        if (r<population[ipop]._mutation_map[ichr].mutation_rate[i])
        {
            // here we have a mutation in the interval [_mutation_map[ichr].bp[i-1],_mutation_map[ichr].bp[i])
            unsigned long int st = population[ipop]._mutation_map[ichr].bp[i-1];
            unsigned long int en = population[ipop]._mutation_map[ichr].bp[i];
            // find its exact bp position
            std::uniform_int_distribution<unsigned long int> distribution(st,en);
            unsigned long int bp_mut = distribution(generator);
            // is the mutation in pat or mat?
            int h01=std::rand() % 2;
            // add the mutation in the appropriate part
            if(h01==0) // add to pat
            {
                for (unsigned j=0; j<hap_rec_pat.size(); j++) // for parts
                {
                    if (hap_rec_pat[j].check_interval(bp_mut)) // if the mutaion occures in this part
                    {
                        hap_rec_pat[j].mutation_pos.push_back(bp_mut);
                        continue;
                    }
                }
                
            }
            else // add to mat
            {
                for (unsigned j=0; j<hap_rec_mat.size(); j++) // for parts
                {
                    if (hap_rec_mat[j].check_interval(bp_mut)) // if the mutaion occures in this part
                    {
                        hap_rec_mat[j].mutation_pos.push_back(bp_mut);
                        continue;
                    }
                }
            }

        }
        
    } // for i
    return true;
}


/*
bool Simulation::ras_add_mutation(std::vector<part> &hap_rec, double lam, int ichr)
{
    // to do: if new mutation exist, then we should remove the old one (although its prob is small)
    int num_mutation_chr = RasRandomNumber::ras_rpois(lam);
    unsigned seed=ras_now_nanoseconds();

    std::default_random_engine generator(seed);
    // extract info from pop0, h0
    int last_part = population[0].h[0].chr[ichr].Hap[0].size();
    unsigned long int st=population[0].h[0].chr[ichr].Hap[0][0].st;
    unsigned long int en=population[0].h[0].chr[ichr].Hap[0][last_part-1].en;
    std::uniform_int_distribution<unsigned long int> distribution(st,en);
    

    for (int i=0; i<num_mutation_chr; ++i)
    {
        unsigned long int bp_mut = distribution(generator);
        // add the mutations in the appropriate part
        for (unsigned j=0; j<hap_rec.size(); j++) // for parts
        {
            if (hap_rec[j].check_interval(bp_mut)) // if the mutaion occures in this part
            {
                hap_rec[j].mutation_pos.push_back(bp_mut);
                continue;
            }
        }
    }

    return true;
}
*/

// for alpha model
// this function also checks the mutations
double Simulation::ras_compute_bv_part(std::vector<part> &p, int ichr, int iphen)
{
    // if there is a mutation in a cv, then its value (cv_val) will change to (cv_val=1-cv_val)
    double ret=0;
    for (unsigned j=0; j<p.size(); j++) // for parts
    {
        int root_population=p[j].root_population;
        unsigned hap_index=p[j].hap_index;
        unsigned ncv=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].bp.size();
        for (unsigned icv=0; icv<ncv; icv++) // for all CVs
        {
            unsigned long int bp=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].bp[icv];
            if (p[j].check_interval(bp)) // if there is a cv in this part, extrcat its root_population and hap_index
            {
                double alpha = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].alpha[icv];
                double cv_val = (double)population[root_population]._pheno_scheme[iphen]._cvs[ichr].val[hap_index][icv];
                
                // if cv is mutated
                std::vector<unsigned long int>::iterator it = find(p[j].mutation_pos.begin(), p[j].mutation_pos.end(), bp);
                if (it != p[j].mutation_pos.end())
                {
                    cv_val=1-cv_val;
                    if (_debug)
                        std::cout << "mutation in cv" << std::endl;
                }
                ret += cv_val * alpha;
                // we can't use coutinue, since it is possible that there are more than 1 cv in this interval
            }
        }
    }
    return ret;
}


bool Simulation::ras_compute_breeding_value(int ipop)
{
    unsigned nchr = population[ipop].h[0].chr.size();
    unsigned nphen = population[ipop]._pheno_scheme.size();
    unsigned long int n_human = population[ipop].h.size();
    
    
    // computing the breeding value for each phenotype
    for (unsigned iphen=0; iphen<nphen; iphen++)
    {
        //if (_debug) std::cout << "iphen" << iphen << std::endl;
        for (unsigned ichr=0; ichr<nchr; ichr++) // for each chromosome
        {
            //if (_debug) std::cout << "ichr" << ichr << std::endl;
            unsigned ncv=population[ipop]._pheno_scheme[iphen]._cv_info[ichr].bp.size();
            std::vector<Human_CV> h_cv(n_human);
            for(unsigned long int ih=0; ih<n_human; ih++) // for humans
            {
                h_cv[ih]=ras_find_cv(population[ipop].h[ih], ichr, iphen, ncv);
            }
            // we have all cvs
            //if (_debug) std::cout << "computing frq" << std::endl;
            std::vector<double> frq(ncv);
            for (unsigned icv=0; icv<ncv; icv++) // for all CVs
            {
                double f=0;
                for(unsigned long int ih=0; ih<n_human; ih++) // for humans
                {
                    f += h_cv[ih].chromatid[0].cv[icv] + h_cv[ih].chromatid[1].cv[icv];
                }
                frq[icv] = f/(2*n_human);
                if (frq[icv]==0 || frq[icv]==1)
                {
                    std::cout << "Error: frq[icv]=" << frq[icv] << "; valid range 0<frq[icv]<1; in ichr=" << ichr << ", icv=" << icv << std::endl;
                    return false;
                }
            }
            
            //if (_debug) std::cout << "A,D,bv per chr" << std::endl;
            for(unsigned long int ih=0; ih<n_human; ih++) // for humans
            {
                //if (_debug) std::cout << "h" << ih << std::endl;
                double A=0;
                double D=0;
                for (unsigned icv=0; icv<ncv; icv++) // for all CVs
                {
                    //if (_debug) std::cout << "icv" << icv << std::endl;
                    // avarage, since they may be from different populations
                    double a = (h_cv[ih].chromatid[0].genetic_value_a[icv]+h_cv[ih].chromatid[1].genetic_value_a[icv])/2;
                    double d = (h_cv[ih].chromatid[0].genetic_value_d[icv]+h_cv[ih].chromatid[1].genetic_value_d[icv])/2;
                    if (population[ipop]._pheno_scheme[iphen]._vd==0)
                        d=0;
                    int t = h_cv[ih].chromatid[0].cv[icv] + h_cv[ih].chromatid[1].cv[icv]; // 0,1,2
                    double p=frq[icv];
                    double q=1-p;
                    // Additive
                    double alpha = a+d*(q-p);
                    A += ((double)t-2*p)*alpha/std::sqrt(2*p*q);
                    // Dominance
                    std::vector<double> c_t(3);
                    c_t[0] = -2*p*p;
                    c_t[1] =  2*p*q;
                    c_t[2] = -2*q*q;
                    D += c_t[t]*d/(2*p*q);
                }
                //population[ipop].h[ih].chr[ichr].additive_chr[iphen]=G_hat;
                //population[ipop].h[ih].chr[ichr].dominance_chr[iphen]=G-G_hat;
                //population[ipop].h[ih].chr[ichr].bv_chr[iphen]=G;
                population[ipop].h[ih].chr[ichr].additive_chr[iphen]=A;
                population[ipop].h[ih].chr[ichr].dominance_chr[iphen]=D;
                population[ipop].h[ih].chr[ichr].bv_chr[iphen]=A+D;
                if (isnan(A))
                {
                    std::cout << "Error: A is nan for human " << ih << ", and chr index " << ichr << std::endl;
                    return false;
                }
            }
        }
        
    }

    
    //if (_debug) std::cout << "A,D,bv" << std::endl;
    // for all pheno, add all chrs
    for(unsigned long int ih=0; ih<n_human; ih++) // for humans
    {
        for (unsigned iphen=0; iphen<nphen; iphen++)
        {
            double bv=0;
            double add=0;
            double dom=0;
            for (unsigned ichr=0; ichr<nchr; ichr++)
            {
                add+=population[ipop].h[ih].chr[ichr].additive_chr[iphen];
                dom+=population[ipop].h[ih].chr[ichr].dominance_chr[iphen];
                bv+=population[ipop].h[ih].chr[ichr].bv_chr[iphen];
            }
            population[ipop].h[ih].additive[iphen]=add;
            population[ipop].h[ih].dominance[iphen]=dom;
            population[ipop].h[ih].bv[iphen]=bv;
        }
    }

    return true;
}


Human_CV Simulation::ras_find_cv(Human &h, unsigned ichr, unsigned iphen, unsigned ncv)
{
    Human_CV h_cv(ncv);
    
    // for pat
    for (unsigned j=0; j<h.chr[ichr].Hap[0].size(); j++) // for parts
    {
        int root_population=h.chr[ichr].Hap[0][j].root_population;
        unsigned hap_index=h.chr[ichr].Hap[0][j].hap_index;
        for (unsigned icv=0; icv<ncv; icv++) // for all CVs
        {
            unsigned long int bp = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].bp[icv];
            if (h.chr[ichr].Hap[0][j].check_interval(bp)) // if there is a cv in this part, extraxt its value from its root_population and hap_index
            {
                bool cv_val = population[root_population]._pheno_scheme[iphen]._cvs[ichr].val[hap_index][icv];
                // if cv is mutated
                std::vector<unsigned long int>::iterator it = find(h.chr[ichr].Hap[0][j].mutation_pos.begin(), h.chr[ichr].Hap[0][j].mutation_pos.end(), bp);
                if (it != h.chr[ichr].Hap[0][j].mutation_pos.end())
                {
                    cv_val = !cv_val;
                    if (_debug) std::cout << "mutation in cv" << std::endl;
                }
                h_cv.chromatid[0].cv[icv] = cv_val;
                h_cv.chromatid[0].genetic_value_a[icv] = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].genetic_value_a[icv];
                h_cv.chromatid[0].genetic_value_d[icv] = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].genetic_value_d[icv];
                // we can't use coutinue, since it is possible that there are more than 1 cv in this interval
            }//if
        }//for cv
    }// for parts
    
    // for mat
    for (unsigned j=0; j<h.chr[ichr].Hap[1].size(); j++) // for parts
    {
        int root_population=h.chr[ichr].Hap[1][j].root_population;
        unsigned hap_index=h.chr[ichr].Hap[1][j].hap_index;
        for (unsigned icv=0; icv<ncv; icv++) // for all CVs
        {
            unsigned long int bp=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].bp[icv];
            if (h.chr[ichr].Hap[1][j].check_interval(bp)) // if there is a cv in this part, extrcat its root_population and hap_index
            {
                bool cv_val = population[root_population]._pheno_scheme[iphen]._cvs[ichr].val[hap_index][icv];
                // if cv is mutated
                std::vector<unsigned long int>::iterator it = find(h.chr[ichr].Hap[1][j].mutation_pos.begin(), h.chr[ichr].Hap[1][j].mutation_pos.end(), bp);
                if (it != h.chr[ichr].Hap[1][j].mutation_pos.end())
                {
                    cv_val = !cv_val;
                    if (_debug) std::cout << "mutation in cv" << std::endl;
                }
                h_cv.chromatid[1].cv[icv] = cv_val;
                h_cv.chromatid[1].genetic_value_a[icv] = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].genetic_value_a[icv];
                h_cv.chromatid[1].genetic_value_d[icv] = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].genetic_value_d[icv];
                // we can't use coutinue, since it is possible that there are more than 1 cv in this interval
            }//if
        }//for cv
    }// for parts

    return h_cv;
}

// for additive-dominance model1

/*
// this function also check s the mutations
// if there is a mutation in a cv, then its value (cv_val) will change to (cv_val=1-cv_val)
double Simulation::ras_compute_bv_additive_dominance(int ichr, int iphen, std::vector<part> &hap_rec_pat, std::vector<part> &hap_rec_mat)
{
    unsigned ncv=population[0]._pheno_scheme[iphen]._cv_info[ichr].bp.size();
    std::vector<int> cv_ind(ncv,0); // 0, 1 or 2
    std::vector<double> genotypic_value_hap0(ncv);
    std::vector<double> dominance_value_hap0(ncv);
    std::vector<double> genotypic_value_hap1(ncv);
    std::vector<double> dominance_value_hap1(ncv);
    
    // for pat
    for (unsigned j=0; j<hap_rec_pat.size(); j++) // for parts
    {
        int root_population=hap_rec_pat[j].root_population;
        unsigned hap_index=hap_rec_pat[j].hap_index;
        for (unsigned icv=0; icv<ncv; icv++) // for all CVs
        {
            unsigned long int bp = population[root_population]._pheno_scheme[iphen]._cv_info[ichr].bp[icv];
            if (hap_rec_pat[j].check_interval(bp)) // if there is a cv in this part, extrcat its root_population and hap_index
            {
                int cv_val = (int)population[root_population]._pheno_scheme[iphen]._cvs[ichr].val[hap_index][icv];
                // if cv is mutated
                std::vector<unsigned long int>::iterator it = find(hap_rec_pat[j].mutation_pos.begin(), hap_rec_pat[j].mutation_pos.end(), bp);
                if (it != hap_rec_pat[j].mutation_pos.end())
                {
                    cv_val=1-cv_val;
                    if (_debug) std::cout << "mutation in cv" << std::endl;
                }
                genotypic_value_hap0[icv]=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].genotypic_value[icv];
                dominance_value_hap0[icv]=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].dominance_value[icv];
                cv_ind[icv] += cv_val;
                // we can't use coutinue, since it is possible that there are more than 1 cv in this interval
            }//if
        }//for cv
    }// for parts
    
    // for mat
    for (unsigned j=0; j<hap_rec_mat.size(); j++) // for parts
    {
        int root_population=hap_rec_mat[j].root_population;
        unsigned hap_index=hap_rec_mat[j].hap_index;
        for (unsigned icv=0; icv<ncv; icv++) // for all CVs
        {
            unsigned long int bp=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].bp[icv];
            if (hap_rec_mat[j].check_interval(bp)) // if there is a cv in this part, extrcat its root_population and hap_index
            {
                int cv_val = (int)population[root_population]._pheno_scheme[iphen]._cvs[ichr].val[hap_index][icv];
                // if cv is mutated
                std::vector<unsigned long int>::iterator it = find(hap_rec_mat[j].mutation_pos.begin(), hap_rec_mat[j].mutation_pos.end(), bp);
                if (it != hap_rec_mat[j].mutation_pos.end())
                {
                    cv_val=1-cv_val;
                    if (_debug) std::cout << "mutation in cv" << std::endl;
                }
                genotypic_value_hap1[icv]=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].genotypic_value[icv];
                dominance_value_hap1[icv]=population[root_population]._pheno_scheme[iphen]._cv_info[ichr].dominance_value[icv];
                cv_ind[icv] += cv_val;
                // we can't use coutinue, since it is possible that there are more than 1 cv in this interval
            }//if
        }//for cv
    }// for parts

    double ret=0;
    for (unsigned icv=0; icv<ncv; icv++) // for all CVs
    {
        // it is possible that genotypic_value be different for different populations
        double a = (genotypic_value_hap0[icv]+genotypic_value_hap1[icv])/2;
        double d = (dominance_value_hap0[icv]+dominance_value_hap1[icv])/2;
        if (cv_ind[icv]==0)
            ret += -a;
        else if(cv_ind[icv]==1)
            ret += d;
        else if(cv_ind[icv]==2)
            ret += a;
    }

    return ret;
}
*/



std::vector<part> Simulation::recombine(chromosome &d1, int starting_haplotype, std::vector<unsigned long int> &recombination_locs)
{
    // recombination_locs is c(start,...,end)
    std::vector<part> ret;
    
    int hap_index=starting_haplotype;
    
    if (recombination_locs.size()<3) return d1.Hap[hap_index];
    
    //ret.reserve(max(l1,l2)+(long int)recombination_locs.size());
    
    for (unsigned long int i1=1; i1<recombination_locs.size(); i1++)
    {
        unsigned long int i2=0;
        
        while(d1.Hap[hap_index].size()>i2 && d1.Hap[hap_index][i2].en <= recombination_locs[i1-1])
        {
            i2++;
        }
        if(d1.Hap[hap_index].size()>i2 && d1.Hap[hap_index][i2].st < recombination_locs[i1-1] && recombination_locs[i1-1] < d1.Hap[hap_index][i2].en && recombination_locs[i1]<d1.Hap[hap_index][i2].en)
        {
            part p=d1.Hap[hap_index][i2];
            p.st=recombination_locs[i1-1];
            p.en=recombination_locs[i1];
            modify_part_for_mutation_pos(p);
            ret.push_back(p);
            i2++;
        }
        if(d1.Hap[hap_index].size()>i2 && d1.Hap[hap_index][i2].st < recombination_locs[i1-1] && recombination_locs[i1-1] < d1.Hap[hap_index][i2].en && recombination_locs[i1] >= d1.Hap[hap_index][i2].en)
        {
            part p=d1.Hap[hap_index][i2];
            p.st=recombination_locs[i1-1];
            modify_part_for_mutation_pos(p);
            ret.push_back(p);
            i2++;
        }
        while(d1.Hap[hap_index].size()>i2 && d1.Hap[hap_index][i2].en <= recombination_locs[i1] &&
              recombination_locs[i1-1]<= d1.Hap[hap_index][i2].st)
        {
            part p=d1.Hap[hap_index][i2];
            modify_part_for_mutation_pos(p);
            ret.push_back(p);
            i2++;
        }
        if(d1.Hap[hap_index].size()>i2 && d1.Hap[hap_index][i2].st < recombination_locs[i1] &&
           recombination_locs[i1] < d1.Hap[hap_index][i2].en)
        {
            part p=d1.Hap[hap_index][i2];
            p.en=recombination_locs[i1];
            modify_part_for_mutation_pos(p);
            ret.push_back(p);
        }
        hap_index=(hap_index+1)%2;
    }
    return ret;
}


void Simulation::modify_part_for_mutation_pos(part &p)
{
    std::vector<unsigned long int> ret;
    for(unsigned i=0; i<p.mutation_pos.size(); i++)
    {
        if(p.check_interval(p.mutation_pos[i]))
            ret.push_back(p.mutation_pos[i]);
    }
    p.mutation_pos = ret;
}


std::vector<unsigned long int> Simulation::ras_sim_loc_rec(std::vector<double> &recom_prob, rMap &rmap, unsigned seed)
{
    //if(seed==0) seed=unsigned(std::time(0));
    //unsigned seed=ras_now_nanoseconds();
    std::srand(seed);
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    std::vector<unsigned long int> rec_bp_loc_mat;
    rec_bp_loc_mat.clear();
    // start of rec_bp_loc_mat should be rmap.bp[0]
    rec_bp_loc_mat.push_back(rmap.bp[0]);
    for (unsigned long int j=0; j< recom_prob.size(); j++)
    {
        //double r=((double)std::rand())/RAND_MAX;
        double r=distribution(generator);
        if (r<recom_prob[j])
            rec_bp_loc_mat.push_back(rmap.bp[j] + std::rand() % rmap.bp_dist_in_rmap);
    }
    // end of rec_bp_loc_mat should be rmap.bp[end]
    rec_bp_loc_mat.push_back(rmap.bp[rmap.bp.size()-1]);
    return rec_bp_loc_mat;
}


double Simulation::compute_couple_var_bv(int ipop, int sex, int iphen)
{
    unsigned long int n_couples=population[ipop]._couples_info.size();
    if (sex==1)//male
    {
        std::vector<double> bv_male;
        bv_male.reserve(n_couples);
        for (unsigned long int i=0; i<n_couples; i++)
        {
            if(!population[ipop]._couples_info[i].inbreed)
            {
                unsigned long int pos_m=population[ipop]._couples_info[i].pos_male;
                bv_male.push_back(population[ipop].h[pos_m].bv[iphen]);
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
            if(!population[ipop]._couples_info[i].inbreed)
            {
                unsigned long int pos_f=population[ipop]._couples_info[i].pos_female;
                bv_female.push_back(population[ipop].h[pos_f].bv[iphen]);
            }
        }
        return CommFunc::var(bv_female);
    }
    else return 0;
}


double Simulation::compute_couple_cor_bv(int ipop, int iphen)
{
    unsigned long int n_couples=population[ipop]._couples_info.size();
    std::vector<double> bv_male(n_couples);
    std::vector<double> bv_female(n_couples);
    
    for (unsigned long int i=0; i<n_couples; i++)
    {
        unsigned long int pos_m=population[ipop]._couples_info[i].pos_male;
        unsigned long int pos_f=population[ipop]._couples_info[i].pos_female;
        bv_male[i]=population[ipop].h[pos_m].bv[iphen];
        bv_female[i]=population[ipop].h[pos_f].bv[iphen];
    }
    return CommFunc::cor(bv_male,bv_female);
}

double Simulation::compute_couple_cor_phen(int ipop, int iphen)
{
    unsigned long int n_couples = population[ipop]._couples_info.size();
    std::vector<double> bv_male(n_couples);
    std::vector<double> bv_female(n_couples);
    
    for (unsigned long int i=0; i<n_couples; i++)
    {
        unsigned long int pos_m=population[ipop]._couples_info[i].pos_male;
        unsigned long int pos_f=population[ipop]._couples_info[i].pos_female;
        bv_male[i]=population[ipop].h[pos_m].phen[iphen];
        bv_female[i]=population[ipop].h[pos_f].phen[iphen];
    }
    return CommFunc::cor(bv_male,bv_female);
}


double Simulation::compute_couple_cor_mating_value(int ipop)
{
    unsigned long int n_couples = population[ipop]._couples_info.size();
    std::vector<double> b_male(n_couples);
    std::vector<double> b_female(n_couples);
    
    for (unsigned long int i=0; i<n_couples; i++)
    {
        unsigned long int pos_m=population[ipop]._couples_info[i].pos_male;
        unsigned long int pos_f=population[ipop]._couples_info[i].pos_female;
        b_male[i]=population[ipop].h[pos_m].mating_value;
        b_female[i]=population[ipop].h[pos_f].mating_value;
    }
    return CommFunc::cor(b_male,b_female);
}



// this func finally creates bv
bool Simulation::ras_initial_human_gen0(int ipop)
{
    unsigned seed=ras_now_nanoseconds();
    std::srand(seed);
 
    std::cout << "     Initializing humans at generation zero for population " << ipop+1 << std::endl;
    int nphen=population[ipop]._pheno_scheme.size();
    unsigned long int nhaps=population[ipop]._pheno_scheme[0]._cvs[0].val.size(); // for gen0, h.size() is 0
    unsigned long int n_people=nhaps/2;
    int nchr=population[ipop]._nchr;
    
    // allocate memory for humans
    ras_allocate_memory_for_humans(population[ipop].h, n_people, nchr, nphen);

    //population[ipop].h.resize(n_people);
    //_Pop_info_prev_gen[ipop].mating_value.resize(n_people,0); // for gen0, set them to 0
    std::cout << "        nhaps = " << nhaps << std::endl;
    
    for (unsigned long int i=0; i<n_people; i++)
    {
        for (int ichr=0; ichr<nchr; ichr++)
        {
            // initialize population[ipop].h[i].chr[ichr].Hap
            long int st=population[ipop]._rmap[ichr].bp[0];
            long int en=population[ipop]._rmap[ichr].bp[population[ipop]._rmap[ichr].bp.size()-1];
            part p1(st,en,2*i,population[ipop]._pop_num);
            population[ipop].h[i].chr[ichr].Hap[0].push_back(p1);
            part p2(st,en,2*i+1,population[ipop]._pop_num);
            population[ipop].h[i].chr[ichr].Hap[1].push_back(p2);
        }
        population[ipop].h[i].sex=(std::rand() % 2)+1; // 1=male, 2=female
        population[ipop].h[i].ID=i;
        population[ipop].h[i].ID_Father=i;// all should be i, in order to check inbreed
        population[ipop].h[i].ID_Mother=i;
        population[ipop].h[i].ID_Fathers_Father=i;
        population[ipop].h[i].ID_Fathers_Mother=i;
        population[ipop].h[i].ID_Mothers_Father=i;
        population[ipop].h[i].ID_Mothers_Mother=i;
    }
    
    
    if (_debug)
        std::cout << " end of [ras_initial_human_gen0]" << std::endl;

    return true;
}

bool Simulation::ras_create_pheno(int ipop, int iphen, double s2_a_gen0, double s2_d_gen0)
{
    unsigned seed=ras_now_nanoseconds();
    std::default_random_engine generator(seed);
    
    std::normal_distribution<double> distribution(0.0, 1);
    
    double beta=population[ipop]._pheno_scheme[iphen]._beta;

    // Get breeding values and parental_effect
    unsigned long int n_human=population[ipop].h.size();
    std::vector<double> a(n_human);
    std::vector<double> d(n_human);
    //std::vector<double> bv(n_human);
    std::vector<double> e(n_human);
    std::vector<double> par_eff(n_human);
    for (unsigned long int i=0; i<n_human; i++)
    {
        // e_noise
        e[i] = distribution(generator); // create noise with N(0,1)

        a[i]=population[ipop].h[i].additive[iphen];
        d[i]=population[ipop].h[i].dominance[iphen];
        //bv[i]=population[ipop].h[i].bv[iphen];
        // parental_effect: transmission_of_environmental_effects_from_parents_to_offspring
        unsigned long int ind_f = population[ipop].h[i].ID_Father;
        unsigned long int ind_m = population[ipop].h[i].ID_Mother;
        double f_father = _Pop_info_prev_gen[ipop].phen[iphen][ind_f];
        double f_mother = _Pop_info_prev_gen[ipop].phen[iphen][ind_m];
        par_eff[i] = beta*(f_father+f_mother)/std::sqrt(2);
    }
    
    // additive
    double s_a=1;
    if (population[ipop]._pheno_scheme[iphen]._va>0)
        s_a=std::sqrt(s2_a_gen0 / population[ipop]._pheno_scheme[iphen]._va);
    else if (population[ipop]._pheno_scheme[iphen]._va==-1) // for using the a and d columns in cv_info file
        s_a=1;
    
    
    //dominance
    double s_d=0;
    if (population[ipop]._pheno_scheme[iphen]._vd>0)
        s_d=std::sqrt(s2_d_gen0 / population[ipop]._pheno_scheme[iphen]._vd);
    else if (population[ipop]._pheno_scheme[iphen]._vd==-1) // for using the a and d columns in cv_info file
        s_d=1;
    
    //double s_bv=std::sqrt(s2_gen0 / population[ipop]._pheno_scheme[iphen]._va);
    //double m_bv=CommFunc::mean(bv);
    double s_par_eff=0;
    if (population[ipop]._pheno_scheme[iphen]._vf>0)
        s_par_eff=std::sqrt(CommFunc::var(par_eff) / population[ipop]._pheno_scheme[iphen]._vf);
    
    // almost no need
    double s_ev=0;
    if (population[ipop]._pheno_scheme[iphen]._ve>0)
        s_ev=std::sqrt(CommFunc::var(e) / population[ipop]._pheno_scheme[iphen]._ve);

    
    //Create phenotypes
    // bv shold be standardized
    for (unsigned long int i=0; i<n_human; i++)
    {
        // e_noise
        if (s_ev>0)
            population[ipop].h[i].e_noise[iphen] = e[i]/s_ev;
        else
            population[ipop].h[i].e_noise[iphen] = 0; // if s_ev==0, then e_noise=0
        // additive
        population[ipop].h[i].additive[iphen] = a[i]/s_a; // scale to gen0
        // dominance
        if (s_d>0)
            population[ipop].h[i].dominance[iphen] = d[i]/s_d; // scale to gen0
        else
            population[ipop].h[i].dominance[iphen] = 0;
        // bv
        population[ipop].h[i].bv[iphen] = population[ipop].h[i].additive[iphen]+population[ipop].h[i].dominance[iphen];
        // parental_effect
        if (s_par_eff>0)
            population[ipop].h[i].parental_effect[iphen] = par_eff[i]; // /s_par_eff; // for gen0 _Pop_info_prev_gen[ipop] is 0
        else
            population[ipop].h[i].parental_effect[iphen] = 0;
        // phen
        population[ipop].h[i].phen[iphen]= population[ipop].h[i].bv[iphen] + population[ipop].h[i].common_sibling[iphen] + population[ipop].h[i].e_noise[iphen] + population[ipop].h[i].parental_effect[iphen];
    }
    
    return true;
}

bool Simulation::ras_save_human_info_to_Pop_info_prev_gen(int ipop)
{
    unsigned long int n_human=population[ipop].h.size();
    int nphen=population[ipop]._pheno_scheme.size();
    std::vector<double> mating_value(n_human); // for each ind
    std::vector<double> selection_value(n_human); // for each ind
    std::vector<std::vector<double> > phen(nphen, std::vector<double>(n_human)); // for each phen and ind

    for (unsigned long int i=0; i<n_human; i++)
    {
        mating_value[i] = population[ipop].h[i].mating_value;
        selection_value[i] = population[ipop].h[i].selection_value;
        for (int iphen=0; iphen<nphen; iphen++)
            phen[iphen][i] = population[ipop].h[i].phen[iphen];
    }
    _Pop_info_prev_gen[ipop].mating_value=mating_value;
    _Pop_info_prev_gen[ipop].selection_value=selection_value;
    _Pop_info_prev_gen[ipop].phen=phen;

    return true;
}


bool Simulation::ras_fill_Pop_info_prev_gen_for_gen0_prev(int ipop)
{

    unsigned long int n=population[ipop].h.size();
    int nphen=population[ipop]._pheno_scheme.size();
    _Pop_info_prev_gen[ipop].mating_value.resize(n);
    _Pop_info_prev_gen[ipop].selection_value.resize(n);
    _Pop_info_prev_gen[ipop].phen.resize(nphen, std::vector<double>(n));
    
    for (int iphen=0; iphen<nphen; iphen++)
    {
        double vf=population[ipop]._pheno_scheme[iphen]._vf;
        if (vf != 0)
        {
            unsigned seed=ras_now_nanoseconds()+iphen;
            std::default_random_engine generator(seed);
            std::normal_distribution<double> distribution(0.0, std::sqrt(vf));

            for (unsigned long int i=0; i<n; i++)
            {
                _Pop_info_prev_gen[ipop].mating_value[i] = 0;
                _Pop_info_prev_gen[ipop].selection_value[i] = 0;
                _Pop_info_prev_gen[ipop].phen[iphen][i] = distribution(generator); // N(0,_vf)
            }
        }//if
    }
    
    return true;
}


double Simulation::ras_combined_variance(int iphen, double a)
{
    double s2x=0, s2y=0;
    unsigned long int N=0;
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        unsigned long int ni=population[ipop].h.size();
        N += ni;
    }
    std::vector<double> x(N);
    std::vector<double> y(N);
    
    unsigned long int k=0;
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        double bi=a*(2*ipop/(_n_pop-1)-1);
        for (unsigned long int j=0; j<population[ipop].h.size(); j++)
        {
            x[k]=population[ipop].h[j].phen[iphen];
            y[k]=population[ipop].h[j].phen[iphen] + bi;
            k++;
        }
    }
    
    s2x=CommFunc::var(x);
    s2y=CommFunc::var(y);
    
    return s2y-(1+_gamma[iphen])*s2x;
}


bool Simulation::ras_add_environmental_effects_specific_to_each_population(int iphen, double a)
{
    for (int ipop=0; ipop<_n_pop; ipop++)
    {
        double gamma_i=a*(2*ipop/(_n_pop-1)-1);
        for (unsigned long int j=0; j<population[ipop].h.size(); j++)
        {
            population[ipop].h[j].phen[iphen] += gamma_i;
        }
    }

    return true;
}


bool Simulation::ras_compute_mating_value_selection_value(int gen_num, int ipop)
{
    unsigned long int nh=population[ipop].h.size();
    int nphen=population[ipop]._pheno_scheme.size();
    //if (_debug)
    //    std::cout << " Population " << ipop+1 << " has " << nphen << " phenotypes." << std::endl;
    
    std::vector<double> x_sv(nh);
    for (unsigned long int i=0; i<nh; i++)
    {
        double mv=0, sv=0;
        for (int iphen=0; iphen<nphen; iphen++)
        {
            double omega=population[ipop]._pheno_scheme[iphen]._omega;
            double lambda=population[ipop]._pheno_scheme[iphen]._lambda;
            mv += omega*population[ipop].h[i].phen[iphen];
            sv += lambda*population[ipop].h[i].phen[iphen];
            //if (i<30) std::cout << population[ipop].h[i].phen[iphen] << std::endl;
        }
        population[ipop].h[i].mating_value=mv;
        population[ipop].h[i].selection_value=sv;
        x_sv[i]=sv;
    }
    
    // make selection_value standardized to generation 0
    if (gen_num==0)
    {
        _gen0_SV_var[ipop]=CommFunc::var(x_sv);
        _gen0_SV_mean[ipop]=CommFunc::mean(x_sv);
    }

    for (unsigned long int i=0; i<nh; i++)
    {
        double sv_standardized=(x_sv[i]-_gen0_SV_mean[ipop])/std::sqrt(_gen0_SV_var[ipop]);
        population[ipop].h[i].selection_value=ras_selection_func(gen_num, ipop, sv_standardized);
    }
    
    return true;
}


bool Simulation::sim_environmental_effects_specific_to_each_population(int iphen)
{
    //environmental_effects_specific_to_each_population
    if(_gamma[iphen]!=0)
    {
        double x0=10;
        double prc=1e-4;
        //double (Simulation::*function_1d_ras)(double)= &Simulation::ras_combined_variance;
        double ah = NewtonRaphson(iphen,x0, prc);
        std::cout << "        ah                = " << ah << std::endl;
        std::cout << "        s2y-(1+Gamma)*s2x = " << ras_combined_variance(iphen,ah) << std::endl;
        /*
        std::cout << "      before environmental effect:" << std::endl;
        for (int ipop=0; ipop<_n_pop; ipop++) // for each population
        {
            std::cout << "       Pop " << ipop+1 << ": ";
            std::cout << "mean(mv) = " << CommFunc::mean(population[ipop].get_mating_value()) << ",\t";
            std::cout << "var(mv) = " << CommFunc::var(population[ipop].get_mating_value()) << std::endl;
        }
         */
        if(!ras_add_environmental_effects_specific_to_each_population(iphen,ah))
        {
            std::cout << "Error in environmental_effects_specific_to_each_population." << std::endl;
            return false;
        }
        /*
        std::cout << "      after environmental effect:" << std::endl;
        for (int ipop=0; ipop<_n_pop; ipop++) // for each population
        {
            std::cout << "       Pop " << ipop+1 << ": ";
            std::cout << "mean(mv) = " << CommFunc::mean(population[ipop].get_mating_value()) << ",\t";
            std::cout << "var(mv) = " << CommFunc::var(population[ipop].get_mating_value()) << std::endl;
        }
        */
    }
    return true;
}


// z is standardized selection value
// return value is in interval (0,1)
double Simulation::ras_selection_func(int gen_num, int ipop, double z)
{
    if (gen_num==0) // for generation 0, all people can marry
        return 1;
    
    int index_generartion=gen_num-1;
    
    if (population[ipop]._selection_func[index_generartion]=="")    // for default value use (logit 0 1)
    {
        double b0=0;
        double b1=1;
        double y=exp(b0 + b1 * z);
        return y/(1+y);
    }
    else if (population[ipop]._selection_func[index_generartion]=="logit") // inverse logit function, for directional selection
    {
        double b0=population[ipop]._selection_func_par1[index_generartion];
        double b1=population[ipop]._selection_func_par2[index_generartion];
        double y=exp(b0 + b1 * z);
        return y/(1+y);
    }
    else if (population[ipop]._selection_func[index_generartion]=="probit") // inverse probit function=normal CDF, for directional selection
    {
        double mu = population[ipop]._selection_func_par1[index_generartion];
        double sigma = population[ipop]._selection_func_par2[index_generartion];
        return CommFunc::NormalCDF(z, mu, sigma);
    }
    else if (population[ipop]._selection_func[index_generartion]=="stab") // for stabilizing selection
    {
        double mu = population[ipop]._selection_func_par1[index_generartion];
        double sigma = population[ipop]._selection_func_par2[index_generartion];
        return CommFunc::NormalPDF(z, mu, sigma);
    }
    else if (population[ipop]._selection_func[index_generartion]=="thr") // for threshold selection
    {
        double p1 = population[ipop]._selection_func_par1[index_generartion];
        double thr = population[ipop]._selection_func_par2[index_generartion];
        double p2 = 1;
        return (z<=thr ? p1 : p2);
    }
    
    return 1;
}



//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in MB.
//
// On failure, returns 0.0, 0.0

void Simulation::process_mem_usage(double& vm_usage, double& resident_set)
{
    using std::ios_base;
    using std::ifstream;
    using std::string;
    
    vm_usage     = 0.0;
    resident_set = 0.0;
    
    // 'file' stat seems to give the most reliable results
    //
    std::ifstream stat_stream("/proc/self/stat",std::ios_base::in);
    
    // dummy vars for leading entries in stat that we don't care about
    //
    std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;
    
    // the two fields we want
    //
    unsigned long vsize;
    long rss;
    
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
    >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
    >> utime >> stime >> cutime >> cstime >> priority >> nice
    >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
    
    stat_stream.close();
    
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0 / 1024.0;
    resident_set = rss * page_size_kb / 1024.0;
}


