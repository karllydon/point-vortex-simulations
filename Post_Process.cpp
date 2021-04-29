#include "Post_Process.hpp"


void::Post_Process::load_file(int file){
	temp_xy.load("Vortex_xy_"+std::to_string(file)+".txt");
}


void Post_Process::update_4n(double inc_angle, double len) {    //Gives the angle and L to post process in the non integrable dipole/dipole case
	load_file(0);
	psi = inc_angle;
	di_length = len;
	d = temp_xy(1,1)-temp_xy(0,1);
}


void Post_Process::count_files() {	//Counts files printed through the simulation
	std::ifstream file_num_read;
	file_num_read.open("python.txt");
	file_num_read >> file_num;
	file_num_read.close();
}


int Post_Process::get_file_num() {   //prints the private int File_num if necessary
	return file_num;
}

void Post_Process::find_circs(){
	vec init;
	init.load("python.txt");
	circs=init.tail(N);
}

vec Post_Process::get_circs(){
	return circs;
}

void Post_Process::calc_impact() {      //calculates impact explicitly from 3rd vortex, useful in dipole/voprtex case, assumes vortices 1 and 2 form dipole 
	load_file(0);	
	d=temp_xy(1,1)-temp_xy(0,1);
	rho=(temp_xy(0,1)+temp_xy(1,1))/2.0-temp_xy(2,1);
}


void Post_Process::calc_impact_4r() {   //calculated impact excplicitly from the centre of a 2cluster
	load_file(0);	
	di_length=temp_xy(0,0)-temp_xy(2,0);
	d = hypot(temp_xy(2,0)-temp_xy(3,0), temp_xy(2,1)-temp_xy(3,1));
	rho=((temp_xy(0,1)+temp_xy(1,1))/2)-((temp_xy(2,1)+temp_xy(3,1))/2);
}


void Post_Process::calc_impact_gen() {   //calculates impact assuming vortices 1 and 2 form dipole and impact is calculated to (0,0)
	load_file(0);
	di_length=hypot(temp_xy(0,0)-temp_xy(2,0),0);
	d = hypot(temp_xy(1,1)-temp_xy(0,1),0);
	rho=(temp_xy(1,1)+temp_xy(0,1))/2;
}


double Post_Process::get_di_length() {    // Gives private variables di_length...
	return di_length;
}

double Post_Process::get_impact() {      // ...and impact
	return rho;
}

double Post_Process::get_psi() {      // ...and psi
	return psi;
}


mat Post_Process::get_clusters(){
	return clusters;
}

vec Post_Process::get_cluster_separations(){
	return clusterSeps;
}


void Post_Process::calc_theta() {       //calculates angle of vortex 1 from the horizontal 
	double delta_y;
	double delta_x;
	mat A;
	mat B;
	A.load("Vortex_xy_" + std:: to_string(file_num - 3) + ".txt");
	B.load("Vortex_xy_" + std:: to_string(file_num - 2) + ".txt");
	delta_y = B(0, 1) - A(0, 1);
	delta_x = B(0, 0) - A(0, 0);
	theta = atan2(delta_y, delta_x);
}

void Post_Process::cluster_populate(int file){  // populate the cluster matrix of pairs closer than cluster_tol
		calc_mat_separations(file);		
		for (int i=0; i<N-1; i++){
			for (int j=i+1; j<N; j++){
				if (sign(circs(i))==sign(circs(j)) && mat_separations(i,j)<cluster_tol){
					clusters(i,j)=1.0;
					clusters(j,i)=1.0;
				}
			}
		}
		vec unique_rows=vec(N, fill::ones);	
		for (int i=0; i<N-1; i++){
			for (int j=i+1; j<N; j++){
				if ((unique_rows(i)==1) && (approx_equal(clusters.row(i), clusters.row(j), "absdiff", 0.01))) unique_rows(j)=0;
			}
		}		
		clusters.shed_rows(find(unique_rows==0));	
}

void  Post_Process::calc_cluster_seps(void){   //return average separations of clusters from their midpoints 
	clusterSeps = vec(clusters.n_rows, fill::zeros); //define cluster separation vector 
	double cluster_mid_x;
	double cluster_mid_y;
	double clusterN; //number of vortices in cluster 
	for (int i=0; i<(int)clusters.n_rows; i++){
		clusterN=sum(clusters.row(i));
		if (clusterN==1) ;  //dont bother calculating anything for 1clusters, leave length = to zero 
		else {
			cluster_mid_x=dot(clusters.row(i), temp_xy.col(0))/clusterN;  //calculate average xy points of each cluster, temp_xy wiill already be loaded from cluster_populate 
			cluster_mid_y=dot(clusters.row(i), temp_xy.col(1))/clusterN;
			for (int j=0; j<(int)clusters.n_cols; j++){  //loop to add each separation from midpoint to each cluster vortex  
				if (clusters(i,j)==1) clusterSeps(i)+= hypot(temp_xy(j,0)-cluster_mid_x, temp_xy(j,1)-cluster_mid_y);
			}
			clusterSeps(i) *= (2/clusterN); //divide by cluster N and double to take the average radius of the circle circumscribed by the cluster 
		}
	}
}

double Post_Process::get_theta() {                 //returns the calculated scattering angle 
	return theta;
}
/*
vec Post_Process::cluster_find(int file){         //find the index of cluster vortices in a given file
	//first need to give a file(arugment of function, then need to find circs, then need to see how close other vortices of the same circulation are 
	mat data;
	data.load("Vortex_xy_" + to_string(file) + ".txt");
	find_circs();
	uvec posvort=find(circs>0);
	uvec negvort=find(circs<0);
	mat clusvort(N,N, fill::zeros);
	int clustcount=0;
	for (int i=0; i<posvort.elem();i++){
		for (int j=i+1, )
	}
	return clustvort;
}
*/


void Post_Process::print_result() {                            //Prints the scattering angle to a file, with normalized impact parameter
	print.open("Scatter_Angles.txt", std:: ofstream::app);
	print << rho / d << "\t" << theta << endl;
	print.close();
}

void Post_Process::print_result_4n() {                  //prints the scattering angle to a file in the 4N case, with L and angle of incidence with L=L2/L1
	print.open("Scatter_Angles.txt", std:: ofstream::app);
	print << di_length << "\t" << psi << "\t" << theta << endl;
	print.close();
}

void Post_Process::calc_mat_separations(int file){
	load_file(file); 	
	for (int i=0; i<N-1; i++){
		for (int j=1+i; j<N; j++){
			mat_separations(i,j)=hypot(temp_xy(i,0)-temp_xy(j,0), temp_xy(i,1)-temp_xy(j,1));	
			mat_separations(j,i)=mat_separations(i,j);
		}
	}
}

mat Post_Process::get_mat_separations(){
	return mat_separations;
}


double Post_Process::get_mat_separation(int vor1, int vor2){
	return mat_separations(vor1, vor2);
}

vec Post_Process::get_separations(){
	return mat_2_vec(mat_separations);
}

vec Post_Process::min_lengths() {     //return min lengths
	calc_mat_separations(0);
	mat min_length_matrix=mat_separations;	
	uvec minimums;
	for (int i = 1; i < file_num; ++i) {
		calc_mat_separations(i);
		minimums=find(mat_separations<min_length_matrix);
		min_length_matrix(minimums)=mat_separations(minimums);
		
	}	
	return mat_2_vec(min_length_matrix);
}

vec Post_Process::critical_lengths(int vor1, int vor2) {    //gives critical lengths when lij is at minimum 
	calc_mat_separations(0);	
	mat crit_length_matrix=mat_separations;
	for (int i = 1; i < file_num; ++i) {
		calc_mat_separations(i); 
		if (mat_separations(vor1,vor2 ) < crit_length_matrix(vor1,vor2)) crit_length_matrix = mat_separations;
	}
	return mat_2_vec(crit_length_matrix);
}

void Post_Process::print_min_lengths() {    //print minimum lengths to file 
	vec Lengths = min_lengths();
	print.open("Minimum_Lengths.txt", std:: ofstream::app);
	print << rho / d;
	for (int i=0; i<(int)Lengths.n_elem; i++) print << "\t" << Lengths(i) / d;
	print <<endl;
	print.close();
}

void Post_Process::print_critical_lengths(int vor1, int vor2) {     //print the 3N critical lengths to file with the normed impact parameter 
	vec Lengths = critical_lengths(vor1,vor2);
	print.open("Scatter_Lengths.txt", std:: ofstream::app);	
	print << rho/d;
	for (int i=0; i<(int)Lengths.n_elem; i++) print << "\t" << Lengths(i) / d;
	print << endl;
	print.close();
}

void Post_Process::print_critical_lengths_4n(int vor1, int vor2) {     //print the 4N critical lengths to file with the delta L and angle of incidence  
	vec Lengths = critical_lengths(vor1,vor2);
	print.open("critical_lengths.txt", std:: ofstream::app);	
	print << di_length << "\t" << psi << "\t" ;
	for (int i=0; i<(int)Lengths.n_elem; i++) print << "\t" << Lengths(i) / d;
	print << endl;
	print.close();
}

void Post_Process::file_delete_routine() {    //deletes all position files 
	int success = 0;
	for (int i = 0; i < file_num; ++i) {
		if (remove(("Vortex_xy_" + std:: to_string(i) + ".txt").c_str()) != 0) success = 1;
	}
	if (flag_remove_hamil) remove("hamiltonian.txt");
	if (flag_remove_momentum) remove("momentums.txt");
	remove("python.txt");
	if (success == 0) std:: cout << "\nFiles deleted successfully" << endl;
	else std:: cout << "\nOne or more files could not be deleted" << endl;
}



vec Post_Process:: mat_2_vec(mat input){ 				//helper function to return a vector of separations from a separation matrix
	vec outVec=vec(N*(N-1)/2, fill::zeros);
	int indexCounter=0;
	for (int i=0; i<N-1;i++){
		for (int j=i+1; j<N; j++)
		{
			outVec(indexCounter)=input(i,j);
			indexCounter++;
		}
	}
	return outVec;
}
