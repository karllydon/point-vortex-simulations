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
  double get_theta(void);
  void calc_separations(int file); 
  void calc_mat_separations(int file);
  vec min_lengths(void); 
  vec critical_lengths(int v1, int v2);  
  //vec cluster_find(int file);
  void print_result(void);
  void print_result_4n(void);
  void print_min_lengths(void); 
  void print_critical_lengths(int v1, int v2);
  void file_delete_routine(void);
  double get_separation(int i);
  double get_psi(void);
  mat get_clusters(void);
  vec get_separations(void);
  mat get_mat_separations(void);
private:
  int file_num;
  double theta;
  double rho;
  double d;
  vec circs;
  vec separations=vec(N*(N-1)/2, fill::zeros);
  mat mat_separations=mat(N, N, fill::zeros);
  mat clusters=eye(N, N);
  double di_length;
  double psi;
  ofstream print;
};
