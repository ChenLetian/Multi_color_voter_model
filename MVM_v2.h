#ifndef MVM_H
#define MVM_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>

using std::vector;
using namespace std;

class Site{
	public:
		int specie;
		vector<int> vneighbors;
};

//---------------------- Lattice class-------------
class Lattice{

	private:
		const int linear_size;
		const int dimension;
		vector <double> freq_specie;
	public:
		vector <Site> list_sites;
		int num_sites;
		int num_neighbors;
		const int num_species;
		Lattice(int dim, int L, int Ns, vector<double>& vfreq_specie): dimension(dim),linear_size(L),num_species(Ns) 
		{
			num_sites=pow(linear_size,dimension);
			num_neighbors=2*dimension;
			Site site_aux;
			list_sites.clear();
			for(int i=0;i<num_sites;i++) list_sites.push_back(site_aux);
			freq_specie.clear();
			for(int i=0;i<num_species;i++) freq_specie.push_back(vfreq_specie[i]);
		};
		~Lattice();
		void set_neighbors();
		int set_specie(gsl_rng *r);
		int get_specie(int index);
		int get_neighbor0(int index);
		int get_sites();
};

Lattice::~Lattice(){ cout<<"destroyed"<<endl; }
int Lattice::get_sites(){ return num_sites; }
int Lattice::get_specie(int index){ return list_sites[index].specie;	 }
int Lattice::get_neighbor0(int index){ return list_sites[index].vneighbors[0]; }

int Lattice::set_specie(gsl_rng *r){

	double sum;
	double sum_ant;
	double rand_aux;
	
	for(int i=0;i<num_sites;i++){
		sum=0;
		sum_ant=0;
		rand_aux=gsl_rng_uniform(r);
		for(int j=0;j<num_species;j++){
			sum+=freq_specie[j];
			if((rand_aux>=sum_ant)&&(rand_aux<sum)) list_sites[i].specie=j;
			sum_ant=sum;
		}
	}
}

void Lattice::set_neighbors(){
	
	for(int i=0;i<num_sites;i++) list_sites[i].vneighbors.clear();
	//Network 1D
	if(dimension==1){
			
		for(int i=1;i<num_sites-1;i++){
		
			list_sites[i].vneighbors.push_back(i-1);
			list_sites[i].vneighbors.push_back(i+1);
		}		
			
		list_sites[0].vneighbors.push_back(num_sites-1);
		list_sites[0].vneighbors.push_back(1);
	
		list_sites[num_sites-1].vneighbors.push_back(num_sites-2);
		list_sites[num_sites-1].vneighbors.push_back(0);
	}
	
	//Network 2D		
	else if(dimension==2){
	
		int site;
			
		for(int isite=0;isite<linear_size;isite++){
			for(int jsite=0;jsite<linear_size;jsite++){
						
				site=isite*linear_size+jsite;
	
				list_sites[site].vneighbors.push_back(((isite+1+linear_size)%linear_size)*linear_size+jsite);
				list_sites[site].vneighbors.push_back(((isite-1+linear_size)%linear_size)*linear_size+jsite);
				list_sites[site].vneighbors.push_back(isite*linear_size+(jsite+1+linear_size)%linear_size);
				list_sites[site].vneighbors.push_back(isite*linear_size+(jsite-1+linear_size)%linear_size);
				
				//list_sites[site].row=isite;
				//list_sites[site].column=jsite;
			}
		}
	}	
}

//------------------Class Dynamics----

class Dynamics: public Lattice{

	private:
		vector<int> vinterfases;
		vector<int> vinterfase_positions;
		int check_interfase(int index);
		void add_interfase(int index);
		void remove_interfase(int index);
	public:
		Dynamics(int dim, int L, int Ns,vector<double>& vfreq_specie):Lattice(dim,L,Ns,vfreq_specie){};
		int num_interfases;
		int num_vinterfases;
		void get_interfases();
		int get_die(gsl_rng *r);
		int get_reproduce(gsl_rng *r, int site_die);
		void reproduce(int sdie,int srepr);
		void clean_variables();
};
void Dynamics::clean_variables(){

	vinterfases.clear();
	vinterfase_positions.clear();
}
void Dynamics::get_interfases(){
	
	vinterfases.clear();
	vinterfase_positions.clear();
	vinterfase_positions.assign(num_sites,0);

	for(int i=0;i<num_sites;i++){
		if(check_interfase(i)==1) add_interfase(i);
		else vinterfase_positions[i]=-1;
	}
	
	num_vinterfases=vinterfases.size();
	num_interfases=round(num_vinterfases/2.);
	if((num_interfases%2!=0)&&(num_interfases!=1)) num_interfases++;
}

int Dynamics::check_interfase(int index){

	int specie1,specie2;
	int neighbor_aux;
	
	for(int j=0;j<num_neighbors;j++){
		neighbor_aux=list_sites[index].vneighbors[j];
		specie1=list_sites[index].specie;
		specie2=list_sites[neighbor_aux].specie;
		if(specie1!=specie2) return 1;
	}
	return 0;
}
void Dynamics::add_interfase(int index){

	vinterfase_positions[index]=vinterfases.size();
	vinterfases.push_back(index);
}
void Dynamics::remove_interfase(int site){

	int last_index=vinterfases.size()-1;
	int last_element=vinterfases[last_index];
	int current_index=vinterfase_positions[site];
	
	vinterfases[current_index]=last_element;
	vinterfase_positions[last_element]=vinterfase_positions[site];
	vinterfase_positions[site]=-1;
	
	vinterfases.pop_back();
}

int Dynamics::get_die(gsl_rng *r){ return vinterfases[gsl_rng_uniform_int(r, num_vinterfases)]; }
int Dynamics::get_reproduce(gsl_rng *r, int site_die){ return list_sites[site_die].vneighbors[gsl_rng_uniform_int(r, num_neighbors)]; }
void Dynamics::reproduce(int sdie,int srepr){
	
	int neighbor;
	if(list_sites[sdie].specie!=list_sites[srepr].specie){
	
		list_sites[sdie].specie=list_sites[srepr].specie;
		if(check_interfase(sdie)==0) remove_interfase(sdie);
		
		for(int i=0;i<num_neighbors;i++){
			neighbor=list_sites[sdie].vneighbors[i];
			if(vinterfase_positions[neighbor]>=0){
			
				if(check_interfase(neighbor)==0) remove_interfase(neighbor);
			}
			else{
			
				if(check_interfase(neighbor)==1) add_interfase(neighbor);
			}
		}
	}
	num_vinterfases=vinterfases.size();
	num_interfases=round(num_vinterfases/2.);
	if((num_interfases%2!=0)&&(num_interfases!=1)) num_interfases++;
}

//-------------------- Class simulation

class Simulation{

	protected:
		const int tfin;
		const int trials;
		const int refresh_time;
		const int rand_seed;
		
		Dynamics dynamics;
		
		int time_increment;
		time_t itime, ftime;
		fstream FILE;
		const char* CURRENT_FILE_NAME;
		
		Simulation(Simulation&);
		void operator=(Simulation&);
		
	public:
		const int num_points;
		Simulation(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME): 
		dynamics(dynamics_in),tfin(T),trials(itrials),refresh_time(irefresh_time),rand_seed(iseed),num_points(int(log10(double(tfin))*2e2)+2),CURRENT_FILE_NAME(FILE_NAME){}; 
		
		//General simulation
		void sequential_trial(gsl_rng *r);
		void promediate();
		
		//Determined at the concret simulation
		virtual void initialice_time_variables(gsl_rng *r){};
		virtual void initialice_time_measures() {};
		virtual void take_initial_measures() {};
		virtual void take_time_measures(int ipoint) {};
		virtual bool conditions_run(double t){ return false; };
		virtual double get_time_increment() { return 1.; };
		virtual void sequential_time_step(gsl_rng *r){};
		virtual void initiate_counter(){};
		virtual void initialice_global_variables() {};
		virtual void initialice_global_measures() {};
		virtual void write_measures(int trial) {};
};

void Simulation::promediate(){ //Maybe here the random generator!?

	//// ------ GENERADOR DE NUMEROS ALEATORIOS ------//
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	//if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	gsl_rng_default_seed=rand_seed; //pongo semilla a mano
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	////----------------------------------------------//
	
	initialice_global_variables();
	//cout<<" variable global"<<endl;
	initialice_global_measures();
	//cout<<" meas glob"<<endl;
	itime=time(NULL);
	
	for(int itrial=0;itrial<trials;itrial++){
	
		sequential_trial(r);
		//cout<<" seq trial"<<endl;
		write_measures(itrial);
	}
	
	gsl_rng_free (r);
}

void Simulation::sequential_trial(gsl_rng *r){

	double t=1;
	int ipoint=0;
	
	initialice_time_variables(r);
	initialice_time_measures();
	take_initial_measures();
	while(conditions_run(t)){
	
		t+=get_time_increment();
		ipoint=int(log10(t)*2e2)+1;
		sequential_time_step(r);
		take_time_measures(ipoint);
		initiate_counter();
	}
	
	dynamics.clean_variables();
}

//--------- Measure of pairs



//----------

class Simulation_collisions: public Simulation{

	private:
	 	int num_collisions;
		int num_anh;
		int num_coal;
		int num_triplets_coal;
		int num_triplets_anh;
		vector<double> p_anh;
		vector<double> p_coal;
		vector<double> prom_num_interfases;
		vector<double> p_triplets_anh;
		vector<double> p_triplets_coal;
	protected:
		vector<double> points;
	public:
		Simulation_collisions(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME)
		: Simulation(dynamics_in,T,itrials,irefresh_time,iseed,FILE_NAME){};
		
		virtual void sequential_time_step(gsl_rng *r);
		virtual void initialice_time_variables(gsl_rng *r);
		virtual void initialice_time_measures();
		//void take_initial_measures() {}; //In this concrete case there aren't initial measures!
		virtual void take_time_measures(int ipoint);
		virtual bool conditions_run(double t);
		virtual double get_time_increment();
		virtual void get_collisions(int index_die);
		virtual void initiate_counter();
		virtual void get_triplets();
		virtual void initialice_global_variables();
		virtual void initialice_global_measures();
		virtual void write_measures(int itrial);
};
	
void Simulation_collisions::initialice_global_variables(){ 

	dynamics.set_neighbors();
}

void Simulation_collisions::initialice_global_measures(){ 

	p_anh.assign(num_points,0);
	p_coal.assign(num_points,0);
	p_triplets_anh.assign(num_points,0);
	p_triplets_coal.assign(num_points,0);
	points.assign(num_points,0);
	prom_num_interfases.assign(num_points,0);
}

void Simulation_collisions::initialice_time_variables(gsl_rng *r){
 
	dynamics.set_specie(r);
	dynamics.get_interfases();
}

void Simulation_collisions::initialice_time_measures(){

	num_coal=0;
	num_anh=0;
	num_collisions=0;
	num_triplets_coal=0;
	num_triplets_anh=0;
	
}

void Simulation_collisions::initiate_counter(){

	num_coal=0;
	num_anh=0;
	num_collisions=0;
	num_triplets_coal=0;
	num_triplets_anh=0;
}

bool Simulation_collisions::conditions_run(double t){

	if((t<tfin)&&(dynamics.num_vinterfases>0)) return true;
	else return false;
}

double Simulation_collisions::get_time_increment(){ return 1./double(dynamics.num_vinterfases); }
void Simulation_collisions::sequential_time_step(gsl_rng *r){
	
	int index_die,index_repr;

	get_triplets();
	index_die=dynamics.get_die(r);
	get_collisions(index_die);
	index_repr=dynamics.get_reproduce(r,index_die);
	dynamics.reproduce(index_die,index_repr);
}

void Simulation_collisions::get_collisions(int index_die){
	
	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	neighbor1=dynamics.list_sites[index_die].vneighbors[0];
	neighbor2=dynamics.list_sites[index_die].vneighbors[1];
	
	specie1=dynamics.list_sites[neighbor1].specie;
	specie2=dynamics.list_sites[index_die].specie;
	specie3=dynamics.list_sites[neighbor2].specie;
	
	comp1=specie1-specie2;
	comp2=specie2-specie3;
	comp3=specie1-specie3;
	
	if((comp1!=0)&&(comp2!=0)){
	
		num_collisions++;
		if(comp3!=0) num_coal++;
		else num_anh++;
	}
}

void Simulation_collisions::get_triplets(){

	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	for(int i=0;i<dynamics.num_vinterfases;i++){
	
		neighbor1=dynamics.list_sites[i].vneighbors[0];
		neighbor2=dynamics.list_sites[i].vneighbors[1];
	
		specie1=dynamics.list_sites[neighbor1].specie;
		specie2=dynamics.list_sites[i].specie;
		specie3=dynamics.list_sites[neighbor2].specie;
	
		comp1=specie1-specie2;
		comp2=specie2-specie3;
		comp3=specie1-specie3;
	
		if((comp1!=0)&&(comp2!=0)){
		
			if(comp3!=0) num_triplets_coal++;
			else num_triplets_anh++;
		}	
	}
}

void Simulation_collisions::take_time_measures(int ipoint){

	p_triplets_anh[ipoint]+=num_triplets_anh;///double(dynamics.num_vinterfases);
	p_triplets_coal[ipoint]+=num_triplets_coal;///double(dynamics.num_vinterfases);
	p_anh[ipoint]+=num_anh;
	p_coal[ipoint]+=num_coal;
	points[ipoint]++;
	prom_num_interfases[ipoint]+=dynamics.num_vinterfases;
}

void Simulation_collisions::write_measures(int itrial){

	int tlogaux,tlog;
	
	ftime=time(NULL);
	if(difftime(ftime, itime)>refresh_time){
		
		FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc);
		tlogaux=0.;
		for(int i=0;i<num_points;i++){
			tlog=int(pow(10,((i-1)/(2e2))));
			if((points[i]>1)&&(tlog!=tlogaux)){ 
				FILE<<tlog<<" ";
				FILE<<p_anh[i]/points[i]<<" ";
				FILE<<p_coal[i]/points[i]<<" ";
				if((p_anh[i]/points[i]+p_coal[i]/points[i])!=0) FILE<<(p_anh[i]+p_coal[i])/points[i]<<" ";
				else FILE<<1<<" ";
				FILE<<p_triplets_anh[i]/points[i]<<" ";
				FILE<<p_triplets_coal[i]/points[i]<<" ";
				if((p_triplets_anh[i]/points[i]+p_triplets_coal[i]/points[i])!=0) FILE<<(p_triplets_anh[i]+p_triplets_coal[i])/points[i]<<" ";
				else FILE<<1<<" ";
				FILE<<prom_num_interfases[i]/points[i]<<" ";
				FILE<<itrial<<endl;
			}
			tlogaux=tlog;
		}
		FILE.close();
		itime=time(NULL);
	}
}



//------------ Concrete simulations classes

class Simulation_collisions_detailed: public Simulation_collisions{

	private:
		double ****p_triplet;
		int ***num_triplet;
		int num_collisions;

	public:
		Simulation_collisions_detailed(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME,int Ns): 
		Simulation_collisions(dynamics_in,T,itrials,irefresh_time,iseed,FILE_NAME){

  			p_triplet = new double***[Ns];
  			for (int i = 0; i < Ns; ++i) {
    				p_triplet[i] = new double**[Ns];

    				for (int j = 0; j < Ns; ++j){
     				 p_triplet[i][j] = new double*[Ns];
     				 for (int k = 0; k < Ns; ++k) p_triplet[i][j][k] = new double[int(log10(double(tfin))*2e2)+2];
     			}
  			}
  
			num_triplet = new int**[Ns];
  			for (int i = 0; i < Ns; ++i) {
    				num_triplet[i] = new int*[Ns];
    				for (int j = 0; j < Ns; ++j) num_triplet[i][j] = new int[Ns];	
 			}
		};
		
		void initialice_global_measures();
		void initialice_time_measures();
		void get_collisions(int index_die);
		void take_time_measures(int ipoint);
		void write_measures(int itrial);
		void initiate_counter();
};
  
void Simulation_collisions_detailed::initialice_global_measures(){

	for(int i=0;i<dynamics.num_species;i++){
		for(int j=0;j<dynamics.num_species;j++){
			for(int k=0;k<dynamics.num_species;k++){
				for(int l=0;l<num_points;l++) p_triplet[i][j][k][l]=0;
			}
		}
	}
	
	points.assign(num_points,0);
}

void Simulation_collisions_detailed::initialice_time_measures(){

	num_collisions=0;
	
	for(int i=0;i<dynamics.num_species;i++){
		for(int j=0;j<dynamics.num_species;j++){
			for(int k=0;k<dynamics.num_species;k++) num_triplet[i][j][k]=0;
		}
	}
}

void Simulation_collisions_detailed::get_collisions(int index_die){
	
	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	neighbor1=dynamics.list_sites[index_die].vneighbors[0];
	neighbor2=dynamics.list_sites[index_die].vneighbors[1];
	
	specie1=dynamics.list_sites[neighbor1].specie;
	specie2=dynamics.list_sites[index_die].specie;
	specie3=dynamics.list_sites[neighbor2].specie;
	
	comp1=specie1-specie2;
	comp2=specie2-specie3;
	//comp3=specie1-specie3;
	
	if((comp1!=0)&&(comp2!=0)){
	
		num_collisions++;
		num_triplet[specie1][specie2][specie3]++;
	}
}

void Simulation_collisions_detailed::initiate_counter(){

	num_collisions=0;
	
	for(int i=0;i<dynamics.num_species;i++){
		for(int j=0;j<dynamics.num_species;j++){
			for(int k=0;k<dynamics.num_species;k++) num_triplet[i][j][k]=0;
		}
	}
}

void Simulation_collisions_detailed::take_time_measures(int ipoint){
	
	for(int i=0;i<dynamics.num_species;i++){
		for(int j=0;j<dynamics.num_species;j++){
			for(int k=0;k<dynamics.num_species;k++) p_triplet[i][j][k][ipoint]+=num_triplet[i][j][k];//double(num_collisions);
		}
	}
	points[ipoint]++;
}

void Simulation_collisions_detailed::write_measures(int itrial){

	int tlogaux,tlog;
	double sum_probabilities;
	
	ftime=time(NULL);
	if(difftime(ftime, itime)>refresh_time){
		
		FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc );
		tlogaux=0.;
		
		for(int i=0;i<num_points;i++){
			tlog=int(pow(10,((i-1)/(2e2))));
			if((points[i]>1)&&(tlog!=tlogaux)){
				FILE<<tlog<<" ";
				sum_probabilities=0;
				for(int is=0;is<dynamics.num_species;is++){
					for(int js=0;js<dynamics.num_species;js++){
					
						if(is!=js){
							FILE<<p_triplet[is][js][is][i]/points[i]<<" ";
							sum_probabilities+=p_triplet[is][js][is][i];
						}	
					}
				}
				
				for(int is=0;is<dynamics.num_species;is++){
					for(int js=0;js<dynamics.num_species;js++){
						for(int ks=0;ks<dynamics.num_species;ks++){
							if((is!=js)&&(js!=ks)&&(is!=ks)){
								FILE<<p_triplet[is][js][ks][i]/points[i]<<" ";
								sum_probabilities+=p_triplet[is][js][ks][i];
							}
						}
					}
				}
				FILE<<sum_probabilities/points[i]<<" ";
				FILE<<itrial<<endl;
			}
			tlogaux=tlog;
		}
		FILE.close();
		itime=time(NULL);
	}
}

//-----------

class Simulation_collisions_detailed_pairs: public Simulation_collisions{

	private:
		double ***p_pair;
		int **num_pair;
		int num_collisions;

	public:
		Simulation_collisions_detailed_pairs(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME,int Ns): 
		Simulation_collisions(dynamics_in,T,itrials,irefresh_time,iseed,FILE_NAME){

  			p_pair = new double**[Ns];
  			for (int i = 0; i < Ns; ++i) {
    				p_pair[i] = new double*[Ns];
    				for (int k = 0; k < Ns; ++k) p_pair[i][k] = new double[int(log10(double(tfin))*2e2)+2];
  			}
  
			num_pair = new int*[Ns];
  			for (int i = 0; i < Ns; ++i) {
    				num_pair[i] = new int[Ns];	
 			}
		};
		
		void initialice_global_measures();
		void initialice_time_measures();
		void get_collisions(int index_die);
		void take_time_measures(int ipoint);
		void write_measures(int itrial);
		void initiate_counter();
};
  
void Simulation_collisions_detailed_pairs::initialice_global_measures(){

	for(int i=0;i<dynamics.num_species;i++){
		for(int j=0;j<dynamics.num_species;j++){
			for(int l=0;l<num_points;l++) p_pair[i][j][l]=0;
		}
	}
	
	points.assign(num_points,0);
}

void Simulation_collisions_detailed_pairs::initialice_time_measures(){

	num_collisions=0;
	
	for(int i=0;i<dynamics.num_species;i++){
		for(int k=0;k<dynamics.num_species;k++) num_pair[i][k]=0;
	}
}

void Simulation_collisions_detailed_pairs::get_collisions(int index_die){
	
	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	neighbor1=dynamics.list_sites[index_die].vneighbors[0];
	neighbor2=dynamics.list_sites[index_die].vneighbors[1];
	
	specie1=dynamics.list_sites[neighbor1].specie;
	specie2=dynamics.list_sites[index_die].specie;
	specie3=dynamics.list_sites[neighbor2].specie;
	
	comp1=specie1-specie2;
	comp2=specie2-specie3;
	//comp3=specie1-specie3;
	
	if((comp1!=0)&&(comp2!=0)){
	
		num_collisions++;
		num_pair[specie1][specie2]++;
		num_pair[specie2][specie3]++;
	}
}

void Simulation_collisions_detailed_pairs::initiate_counter(){

	num_collisions=0;
	
	for(int i=0;i<dynamics.num_species;i++){
		for(int k=0;k<dynamics.num_species;k++) num_pair[i][k]=0;
	}
}

void Simulation_collisions_detailed_pairs::take_time_measures(int ipoint){
	
	for(int i=0;i<dynamics.num_species;i++){
		for(int k=0;k<dynamics.num_species;k++) p_pair[i][k][ipoint]+=num_pair[i][k];
	}
	points[ipoint]++;
}

void Simulation_collisions_detailed_pairs::write_measures(int itrial){

	int tlogaux,tlog;
	double sum_probabilities;
	
	ftime=time(NULL);
	if(difftime(ftime, itime)>refresh_time){
		
		FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc );
		tlogaux=0.;
		
		for(int i=0;i<num_points;i++){
			tlog=int(pow(10,((i-1)/(2e2))));
			if((points[i]>1)&&(tlog!=tlogaux)){
				FILE<<tlog<<" ";
				sum_probabilities=0;
				for(int is=0;is<dynamics.num_species;is++){
					for(int js=0;js<dynamics.num_species;js++){
						FILE<<(p_pair[is][js][i]+p_pair[js][is][i])/points[i]<<" "; //Tendré columnas repetidas!!!!
						sum_probabilities+=p_pair[is][js][i];	
					}
				}
				FILE<<sum_probabilities/points[i]<<" ";
				FILE<<itrial<<endl;
			}
			tlogaux=tlog;
		}
		FILE.close();
		itime=time(NULL);
	}
}

//------------ Num pair en un tiempo dado

class Simulation_collisions_pairs: public Simulation{

	private:
	 	int num_collisions;
		int num_anh;
		int num_coal;
		int num_pairs_coal;
		int num_pairs_anh;
		vector<double> p_anh;
		vector<double> p_coal;
		vector<double> prom_num_interfases;
		vector<double> p_pairs_anh;
		vector<double> p_pairs_coal;
	protected:
		vector<double> points;
	public:
		Simulation_collisions_pairs(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME)
		: Simulation(dynamics_in,T,itrials,irefresh_time,iseed,FILE_NAME){};
		
		virtual void sequential_time_step(gsl_rng *r);
		virtual void initialice_time_variables(gsl_rng *r);
		virtual void initialice_time_measures();
		//void take_initial_measures() {}; //In this concrete case there aren't initial measures!
		virtual void take_time_measures(int ipoint);
		virtual bool conditions_run(double t);
		virtual double get_time_increment();
		virtual void get_collisions(int index_die);
		virtual void initiate_counter();
		virtual void get_pairs();
		virtual void initialice_global_variables();
		virtual void initialice_global_measures();
		virtual void write_measures(int itrial);
};
	
void Simulation_collisions_pairs::initialice_global_variables(){ 

	dynamics.set_neighbors();
}

void Simulation_collisions_pairs::initialice_global_measures(){ 

	p_anh.assign(num_points,0);
	p_coal.assign(num_points,0);
	p_pairs_anh.assign(num_points,0);
	p_pairs_coal.assign(num_points,0);
	points.assign(num_points,0);
	prom_num_interfases.assign(num_points,0);
}

void Simulation_collisions_pairs::initialice_time_variables(gsl_rng *r){
 
	dynamics.set_specie(r);
	dynamics.get_interfases();
}

void Simulation_collisions_pairs::initialice_time_measures(){

	num_coal=0;
	num_anh=0;
	num_collisions=0;
	num_pairs_coal=0;
	num_pairs_anh=0;
	
}

void Simulation_collisions_pairs::initiate_counter(){

	num_coal=0;
	num_anh=0;
	num_collisions=0;
	num_pairs_coal=0;
	num_pairs_anh=0;
}

bool Simulation_collisions_pairs::conditions_run(double t){

	if((t<tfin)&&(dynamics.num_vinterfases>0)) return true;
	else return false;
}

double Simulation_collisions_pairs::get_time_increment(){ return 1./double(dynamics.num_vinterfases); }
void Simulation_collisions_pairs::sequential_time_step(gsl_rng *r){
	
	int index_die,index_repr;

	get_pairs();
	index_die=dynamics.get_die(r);
	get_collisions(index_die);
	index_repr=dynamics.get_reproduce(r,index_die);
	dynamics.reproduce(index_die,index_repr);
}

void Simulation_collisions_pairs::get_collisions(int index_die){
	
	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	neighbor1=dynamics.list_sites[index_die].vneighbors[0];
	neighbor2=dynamics.list_sites[index_die].vneighbors[1];
	
	specie1=dynamics.list_sites[neighbor1].specie;
	specie2=dynamics.list_sites[index_die].specie;
	specie3=dynamics.list_sites[neighbor2].specie;
	
	comp1=specie1-specie2;
	comp2=specie2-specie3;
	comp3=specie1-specie3;
	
	if((comp1!=0)&&(comp2!=0)){
	
		num_collisions++;
		if(comp3!=0) num_coal++;
		else num_anh++;
	}
}

void Simulation_collisions_pairs::get_pairs(){

	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	for(int i=0;i<dynamics.num_vinterfases;i++){
	
		neighbor1=dynamics.list_sites[i].vneighbors[0];
		neighbor2=dynamics.list_sites[i].vneighbors[1];
	
		specie1=dynamics.list_sites[neighbor1].specie;
		specie2=dynamics.list_sites[i].specie;
		specie3=dynamics.list_sites[neighbor2].specie;
	
		comp1=specie1-specie2;
		comp2=specie2-specie3;
		comp3=specie1-specie3;
	
		if((comp1!=0)&&(comp2!=0)){
		
			if(comp3!=0) num_pairs_coal+=2;
			else num_pairs_anh+=2;
		}	
	}
}

void Simulation_collisions_pairs::take_time_measures(int ipoint){

	p_pairs_anh[ipoint]+=num_pairs_anh;///double(dynamics.num_vinterfases);
	p_pairs_coal[ipoint]+=num_pairs_coal;///double(dynamics.num_vinterfases);
	p_anh[ipoint]+=num_anh;
	p_coal[ipoint]+=num_coal;
	points[ipoint]++;
	prom_num_interfases[ipoint]+=dynamics.num_vinterfases;
}

void Simulation_collisions_pairs::write_measures(int itrial){

	int tlogaux,tlog;
	
	ftime=time(NULL);
	if(difftime(ftime, itime)>refresh_time){
		
		FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc);
		tlogaux=0.;
		for(int i=0;i<num_points;i++){
			tlog=int(pow(10,((i-1)/(2e2))));
			if((points[i]>1)&&(tlog!=tlogaux)){ 
				FILE<<tlog<<" ";
				FILE<<p_anh[i]/points[i]<<" ";
				FILE<<p_coal[i]/points[i]<<" ";
				if((p_anh[i]/points[i]+p_coal[i]/points[i])!=0) FILE<<(p_anh[i]+p_coal[i])/points[i]<<" ";
				else FILE<<1<<" ";
				FILE<<p_pairs_anh[i]/points[i]<<" ";
				FILE<<p_pairs_coal[i]/points[i]<<" ";
				if((p_pairs_anh[i]/points[i]+p_pairs_coal[i]/points[i])!=0) FILE<<(p_pairs_anh[i]+p_pairs_coal[i])/points[i]<<" ";
				else FILE<<1<<" ";
				FILE<<prom_num_interfases[i]/points[i]<<" ";
				FILE<<itrial<<endl;
			}
			tlogaux=tlog;
		}
		FILE.close();
		itime=time(NULL);
	}
}




#endif
