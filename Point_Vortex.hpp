#pragma once

#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <unistd.h>

using namespace std;
using namespace arma;



//------------------------------------------------------------------ Constants + Flags -----------------------------------------------------------------------------------------------------------

const int N = 4;			//number of vortices in syste
const int File_Max = pow(10,5);		//max files allowed to print
const int write_limit=pow(10,0);
const double pi = 3.14159265358979;

const bool Flag_is_Periodic = false; //periodic BCs
const bool Flag_Print_Hamil = true; //print hamiltonian 
const bool Flag_Print_Momen = true; //print momentums

const double tol = pow(10, -12);//tolerance of rk45 error
const double Rtol = pow(10, -12); // relative tolerqance of rk45 error
static double h=0.1;   //step length
const double hmin = pow(10,-7); //min step allowed
const double hmax = 1; //max step allowed
const double hfac=0.9;
const double hfacmax=2;
const double hfacmin=0.5;

					   //-------------------------------------------------------------------- Misc. Functions --------------------------------------------------------------------------------------------------------

double square(double x);

vec timesteps(double a, double b, int m);


double Per_Sum(double x, double y);
double Per_Hamil_Sum(double x, double y);
//----------------------------------------------------------------- Lengths + Velocities --------------------------------------------------------------------------------------------------------
//Functions that we will use to solve the system

double xij(int i, int j, mat &xy);
double yij(int i, int j, mat &xy);

double length(int i, int j, mat &xy);

double X_Vel(int i, vec &circ, mat &xy);

double Per_X_Vel(int i, vec &circ, mat &xy);
double Y_Vel(int i, vec &circ, mat &xy);

double Per_Y_Vel(int i, vec &circ, mat &xy);

mat Vel_Matrix(mat &xy, vec &circ);
mat Per_Vel_Matrix(mat &xy, vec &circ);

//---------------------------------------------------------------------------- Conserved Quantities------------------------------------------------------------------------------------------------------------

double P(vec &circ, mat &xy);
double Q(vec &circ, mat &xy);
double M(vec &circ, mat &xy);
double Hamil(vec &circ, mat &xy);

//-----------------------------------------------------------------------------Timestepping procedures--------------------------------------------------------------------------------------------------------------
double error_func(mat &RK1, mat &RK2); 
void Dormand_Price(vec &t_steps, vec &circ, mat &xy);
void Bound_Dormand_Price(vec &circ, mat &xy);
