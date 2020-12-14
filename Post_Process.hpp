#pragma once
#include <armadillo>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "Point_Vortex.hpp"

using namespace std;
using namespace arma;

const bool flag_remove_hamil = true;
const bool flag_remove_momentum = true;

class Post_Process {
public:
  Post_Process(int i) {
	  if (i != 0) file_num = i;
		else;
	}
  void update_4n(double inc_angle, double Length);
  void count_files(void);
  int get_file_num(void);
  void calc_impact(void);
  void calc_impact_4r(void);
  void calc_impact_gen(void);
  double get_di_length(void); 
  double get_impact(void);
  void calc_theta(void);
  double get_theta(void);
  void calc_separations(int file); 
  vec min_lengths(void); 
  vec critical_lengths(void);  
  void print_result(void);
  void print_result_4n(void);
  void print_min_lengths(void); 
  void print_critical_lengths(void);
  void file_delete_routine(void);
  double get_separation(int i);
  vec get_separations(void);
private:
  int file_num;
  double theta;
  double rho;
  double d;
  vec separations=vec(N*(N-1)/2, fill::zeros);
  double di_length;
  double phi;
  ofstream print;
};
