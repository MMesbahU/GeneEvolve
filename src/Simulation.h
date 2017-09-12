#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath> // for round, floor, ceil, trunc
#include <vector>
#include <cstdlib>      // std::rand, std::srand
#include <random>

#include "Population.h"


#ifdef _OPENMP
#include <omp.h>
#endif




class Pop_phen_info
{
public:
    std::vector<double> mating_value; // for each ind
    std::vector<double> selection_value; // for each ind
    std::vector<std::vector<double> > phen; // for each phen and ind
};



class Simulation
{
public:
    int _n_pop;
    int _tot_gen;
    bool _debug;
    std::vector<Population> population;
    Parameters par;
    //_tot_gen rows, _n_pop^2 cols: each column is the the elements of the transition matrix in format: [a11 a12 a13 ... a21 a22 a23 ...]
    std::vector<std::vector<double> > migration_mat_gen;
    std::string _out_prefix;
    bool _output_all_generations;
    std::vector<int> _file_output_generations_list; // 0-based list of generation for genotype outputs
    bool _out_hap;
    bool _out_plink;
    bool _out_interval;
    bool _out_vcf;
    std::vector<int> _all_active_chrs;
    // for Fprime NewtonRaphson
    //typedef double (Simulation::*function_1d)(double); // function_1d is a pointer to a double f(double) function
    std::vector<Pop_phen_info> _Pop_info_prev_gen; // Saving mating_value for the next generation
    
    // environmental effects specific to each_population
    std::vector<double> _gamma; // for each phenotype and all the populations
    std::vector<double> _gen0_SV_var;
    std::vector<double> _gen0_SV_mean;
    
    
private:
    bool ras_init_parameters(void);
    bool ras_init_generation0(void);
    bool ras_main_sim(void);
    void ras_show_res(void);
    void ras_save_res(void);
    bool read_migration_file(void);
    bool ras_do_migration(int gen_ind);
    bool ras_save_genotypes(int gen_num);
    
    std::vector<part> recombine(chromosome &d1, int starting_haplotype, std::vector<unsigned long int> &recombination_locs);
    void modify_part_for_mutation_pos(part &p);
    std::vector<unsigned long int> ras_sim_loc_rec(std::vector<double> &recom_prob, rMap &rmap, unsigned seed);
    double ras_compute_bv_part(std::vector<part> &p, int ichr, int iphen);
    //double ras_compute_bv_additive_dominance(int ichr, int iphen, std::vector<part> &hap_rec_pat, std::vector<part> &hap_rec_mat);
    Human_CV ras_find_cv(Human &h, unsigned ichr, unsigned int iphen, unsigned ncv);
    bool ras_compute_AD(int ipop, int gen_num);
    std::vector<Human> reproduce(int ipop, int gen_num);
    bool assort_mate(int ipop, int gen_ind);
    bool random_mate(int ipop, int gen_ind);
    bool sim_next_generation(int gen_num);
    bool ras_scale_AD_compute_GEF(int ipop, int iphen, double s2_a_gen0, double s2_d_gen0);
    bool ras_initial_human_gen0(int ipop);
    bool ras_allocate_memory_for_humans(std::vector<Human> &h, unsigned long int  n_people, int nchr, int nphen);

    //environmental_effects_specific_to_each_population
    double Fprime(int iphen, double x);
    double NewtonRaphson(int iphen, double x0, double precision);
    double ras_combined_variance(int iphen, double a);
    bool ras_add_environmental_effects_specific_to_each_population(int iphen, double a);

    //save transmission_of_environmental_effects_from_parents_to_offspring
    bool ras_save_human_info_to_Pop_info_prev_gen(int ipop);
    bool ras_fill_Pop_info_prev_gen_for_gen0_prev(int ipop);

    
    // reading the inputed hap files
    bool ras_read_hap_legend_sample_chr(std::vector<Legend> &pops_legend, std::vector<Hap_SNP> &pops_hap, int ichr);

    // dealing with .hap .legend .sample files
    bool ras_write_hap_legend_sample(int gen_num);
    bool ras_convert_interval_to_hap_matrix(int ipop, std::vector<Hap_SNP> &pops_hap, std::vector<Legend> &pops_legend, int ichr, Hap_SNP &hap_snp);
    bool ras_convert_pop_to_indv(int ipop, std::vector<unsigned long int> &indv_id);

    
    // dealing with plink .ped .map files
    bool ras_write_hap_to_plink_format(int gen_num);
    bool ras_convert_interval_to_format_plink(int ipop, std::vector<Hap_SNP> &pops_hap, std::vector<Legend> &pops_legend, int ichr, std::vector<std::vector<bool> > &matrix_plink_ped, plink_PED_ids &plink_ped_ids, plink_MAP &plink_map);

    // dealing with interval (.int) file
    bool ras_write_hap_to_interval_format(int gen_num);
    
    
    // dealing with vcf files
    bool ras_write_vcf_to_vcf_format(int gen_num);
    bool ras_convert_interval_to_vcf_structure(vcf_structure &vcf_out, int ipop, int ichr, std::vector<vcf_structure> &vcf_structure_allpops_chr);
    bool ras_read_vcf_pops_chr(std::vector<vcf_structure> &vcf_structure_pops, int ichr);

    

    bool ras_compute_mating_value_selection_value(int gen_num, int ipop);
    bool sim_environmental_effects_specific_to_each_population(int iphen);

    // selection func
    double ras_selection_func(int gen_num, int ipop, double x);
    void process_mem_usage(double& vm_usage, double& resident_set);

    // mutation
    //bool ras_add_mutation(std::vector<part> &hap_rec, double lam, int ichr);
    bool ras_add_mutation(int ipop, int ichr, std::vector<part> &hap_rec_pat, std::vector<part> &hap_rec_mat);

    // seed
    std::default_random_engine glob_generator;
    unsigned ras_glob_seed(void);

    // output generations
    bool read_file_output_generation_list(void);

public:
    bool run(void);
    
};


