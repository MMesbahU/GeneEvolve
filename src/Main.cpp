#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <ostream>
#include <iostream>
#include <ctime>
//#include <chrono> // for seconds, milliseconds, nanoseconds, picoseconds
#include <string> // str1.compare(str2)
#include "Simulation.h"
//#include "parameters.h" is Simulation.h
//#include "hap.h" it is in Simulation.h
#include "CommFunc.h"
//#include "Unique.h"

//using Eigen::MatrixXd;
//using namespace Eigen;

//int transFactor = 3;
//int cisFactor = 2;


void ProgramVersion();
void helpFile();


int main(int argc, char ** argv)
{
    int start_time = time(0);

    int cpus;

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

    if (!par.read(vec_arg))
    {
        std::cout << std::endl;
        std::cout << " For more information type [GeneEvolve --help]." << std::endl;
        std::cout << std::endl;
        return -1;
    }

    if (par._help)
    {
        helpFile();
        return -1;
    }
    
    if (!par.check())
    {
        std::cout << std::endl;
        std::cout << " For more information type [GeneEvolve --help]." << std::endl;
        std::cout << std::endl;
        return -1;
    }
    
    
    par.print();
    
    
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
    std::cout << std::endl;
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    std::cout << "                                  GeneEvolve                                   " << std::endl;
    std::cout << " ------------------------------------------------------------------------------" << std::endl;
    std::cout << " (c) 2015-2017 - Rasool Tahmasbi and Matthew C. Keller." << std::endl;
    std::cout << std::endl;
    std::cout << " Version: " << VERSION << std::endl;
    std::cout << " Built: " << DATE << " by " << USER << std::endl;
    std::cout << std::endl;
    std::cout << " URL = https://github.com/rtahmasbi/GeneEvolve" << std::endl;
    //std::cout << " URL = http://matthewckeller.com/html/GeneEvolve.html" << std::endl;
    std::cout << std::endl;
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
    printf("        --file_ref_vcf           : [filename]\n");
    printf("        --file_recom_map         : [filename]\n");
    printf("        --file_mutation_map      : [filename]\n");
    printf("        --RM                     : [off]  ->  Random Mating\n");
    printf("        --MM                     : [0]\n");
    printf("        --vt_type                : [1]\n");
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
    //printf("        --beta                   : [1]    ->  coefficient for familial effect\n");

    std::cout << std::endl;
    printf(" --------- Immigration parameters\n");
    printf("        --file_migration         : [filename]\n");
    
    std::cout << std::endl;
    printf(" --------- Environmental effects specific to each population (for each phenotype)\n");
    printf("        --gamma                  : [0]\n");
    
    std::cout << std::endl;
    printf(" --------- Output parameters\n");
    printf("        --out_hap                : [Off]  ->  [On] means output in the [hap] format.\n");
    printf("        --out_plink              : [Off]  ->  [On] means output in the [plink] format.\n");
    printf("        --out_plink01            : [Off]  ->  [On] means output in the [plink01] format.\n");
    printf("        --out_vcf                : [Off]  ->  [On] means output in the [vcf] format.\n");
    printf("        --out_interval           : [Off]  ->  [On] means output in the [interval] format.\n");
    printf("        --file_output_generations: [filename]\n");
    printf("          A file contating the list of the generations for genotype output. Each generation number should be in a line.\n");
    
    std::cout << std::endl;
    printf(" --------- Other parameters\n");
    printf("        --prefix                 : [out]\n");
    printf("        --avoid_inbreeding       : [Off]  ->  [On] means no inbreeding.\n");
    printf("        --seed                   : [0]    ->  can be any passitive number.\n");
    printf("        --debug                  : [Off]\n");
    

    printf("\n Please visit <https://github.com/rtahmasbi/GeneEvolve> for detailed documentation.\n\n");
    std::cout << std::endl;


}



