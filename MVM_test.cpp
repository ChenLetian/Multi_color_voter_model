#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <algorithm>
#include "MVM_v2.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;				
				
int main(int argc, char *argv[])
{
	int dimension=atoi(argv[1]);
	int linear_size=atoi(argv[2]);
	int num_species=atoi(argv[3]);
	double dif_specie=atoi(argv[4]);
	int tfin=atoi(argv[5]);
	int trials=atoi(argv[6]);
	unsigned int refresh_time=(atoi(argv[7]))*60; //min
	int thread=atoi(argv[8]);
	
	char FILE_NAME [100];
	//sprintf(FILE_NAME,"Collision_%d_%d_%g_%d_%d.txt",linear_size,num_species,dif_specie,tfin,thread);
	sprintf(FILE_NAME,"Collision_detailed_%d_%d_%g_%d_%d.txt",linear_size,num_species,dif_specie,tfin,thread);
	vector <double> freq_specie;
	double slice;
	
	for(int i=0;i<num_species;i++) freq_specie.push_back(1/double(num_species));
	slice=0;
	for(int i=0;i<num_species-1;i++){
	
		freq_specie[i]+=slice;
		slice=freq_specie[i]*(1-1./dif_specie);
		freq_specie[i]-=slice;
	}
	freq_specie[num_species-1]+=slice;
	
	for(int i=0;i<num_species;i++) cout<<i<<" "<<freq_specie[i]<<endl;

	Dynamics dynamics(dimension,linear_size,num_species,freq_specie);
	//Simulation_collisions simulation(dynamics,tfin,trials,refresh_time,thread,FILE_NAME);
	Simulation_collisions_detailed simulation(dynamics,tfin,trials,refresh_time,thread,FILE_NAME,num_species);

	simulation.promediate();

	return EXIT_SUCCESS;

}
