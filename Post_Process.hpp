#pragma once
#include <armadillo>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "Point_Vortex.hpp"


using namespace arma;

const bool flag_remove_hamil = true;
const bool flag_remove_momentum = true;
const double cluster_tol=intervortex*0.5;      //should be intervortex length/2 or something 
const double dipole_tol=intervortex*0.5;      //should be intervortex length/2 or something 

class Post_Process {
public:
  Post_Process(int i) {
	  if (i != 0) file_num = i;
	  else; 
  }
  void load_file(int file);
  void update_4n(double inc_angle, double Length);
  void count_files(void); 
  void find_circs(void);
 

  //calculation of key quantities
  void calc_impact(void);
  void calc_impact_4r(void);
  void calc_impact_gen(void);
  void calc_theta(int vor);
  void dipole_cluster_populate(int file);
  void dipole_cluster_populate_fix(int file);
  void calc_cluster_seps(void); 
  void calc_mat_separations(int file);
  void calc_mat_separations_periodic(int file);
  vec min_lengths(void); 
  mat critical_lengths(int vor1, int vor2);  
  void exchange_find(void);

  //rountines to easily print to files 
  void print_result(void);
  void print_result_4n(void);
  void print_min_lengths(void); 
  void print_critical_lengths(int vor1, int vor2);
  void print_critical_lengths_4n(int vor1, int vor2);
  void print_dipole_extremes(int vor1, int vor2);      //assumes dipoles start with vortex 1 and 2, with vortex 1 being antivortex 
  //helper functions to save time or make things easier  
  vec mat_2_vec(mat input);
  void dipole_cluster_refresh(void){clusters=eye(N,N); dipoles=eye(N,N);}
  void file_delete_routine(void); 

  //access functions 
  int get_file_num(void) {return file_num;};
  vec get_circs(void){return circs;}
  double get_di_length(void){return di_length;}
  double get_impact(void){return rho;}
  double get_theta(void){return theta;}
  double get_psi(void){return psi;}
  mat get_clusters(void){return clusters;}
  vec get_cluster_separations(void){return clusterSeps;}
  vec get_separations(void){return mat_2_vec(mat_separations);}
  double get_mat_separation(int vor1, int vor2){return mat_separations(vor1,vor2);}
  mat get_mat_separations(void){return mat_separations;}
  bool get_exchange(void){return exchange;}
  int get_final_dipole_vortex(void){return final_dipole_vortex;}
  int get_exchange_num(void){return exchange_num;}
  double get_preIntermin(void){return preInterMin;}
  double get_preIntermax(void){return preInterMax;}
  double get_postIntermin(void){return postInterMin;}
  double get_postIntermax(void){return postInterMax;}
  mat get_dipoles(void){ return dipoles;}
private:
  double preInterMin, preInterMax, postInterMin, postInterMax;
  bool exchange; 
  int final_dipole_vortex;
  int exchange_num;
  mat temp_xy;
  int file_num;
  double theta;
  double rho;
  double d;
  vec circs; 
  vec clusterSeps;
  mat mat_separations=mat(N, N, fill::zeros);
  mat clusters=eye(N, N);
  mat dipoles=eye(N, N);
  double di_length;
  double psi;
  std::ofstream print;
};
