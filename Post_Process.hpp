#pragma once
#include <armadillo>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "Point_Vortex.hpp"


using namespace arma;

const bool flag_remove_hamil = true;
const bool flag_remove_momentum = true;
const double cluster_tol=2.0;

class Post_Process {
public:
  Post_Process(int i) {
	  if (i != 0) file_num = i;
	  else; 
  }
  void load_file(int file);
  void update_4n(double inc_angle, double Length);
  void count_files(void);
  int get_file_num(void);
  void find_circs(void);
  vec get_circs(void);
  void calc_impact(void);
  void calc_impact_4r(void);
  void calc_impact_gen(void);
  double get_di_length(void); 
  double get_impact(void);
  void calc_theta(void);
  void cluster_populate(int file);
  void calc_cluster_seps(void);
  double get_theta(void);
  void calc_mat_separations(int file);
  vec min_lengths(void); 
  vec critical_lengths(int vor1, int vor2);  
  void print_result(void);
  void print_result_4n(void);
  void print_min_lengths(void); 
  void print_critical_lengths(int vor1, int vor2);
  void print_critical_lengths_4n(int vor1, int vor2);
  void file_delete_routine(void); 
  double get_psi(void);
  mat get_clusters(void);  
  vec get_cluster_separations(void);
  vec get_separations(void); 
  double get_mat_separation(int vor1, int vor2);
  mat get_mat_separations(void);
  vec mat_2_vec(mat input);
private:
  mat temp_xy;
  int file_num;
  double theta;
  double rho;
  double d;
  vec circs; 
  vec clusterSeps;
  mat mat_separations=mat(N, N, fill::zeros);
  mat clusters=eye(N, N);
  double di_length;
  double psi;
  std::ofstream print;
};
