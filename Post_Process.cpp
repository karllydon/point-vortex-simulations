#include "Post_Process.hpp"

void Post_Process::update_4n(double inc_angle, double len) {    //Gives the angle and L to post process in the non integrable dipole/dipole case
	mat init;
	init.load("Vortex_xy_0.txt");
	phi = inc_angle;
	di_length = len;
	d = init(1,1)-init(0,1);
}


void Post_Process::count_files() {	//Counts files printed through the simulation
	ifstream file_num_read;
	file_num_read.open("python.txt");
	file_num_read >> file_num;
	file_num_read.close();
}


int Post_Process::get_file_num() {   //prints the private int File_num if necessary
	return file_num;
}



void Post_Process::calc_impact() {      //calculates impact explicitly from 3rd vortex, useful in dipole/voprtex case, assumes vortices 1 and 2 form dipole 
	mat init;
	init.load("Vortex_xy_0.txt");
	d=init(1,1)-init(0,1);
	rho=(init(0,1)+init(1,1))/2.0-init(2,1);
}


void Post_Process::calc_impact_4r() {   //calculated impact excplicitly from the centre of a 2cluster
	mat init;
	init.load("Vortex_xy_0.txt");
	di_length=init(0,0)-init(2,0);
	d = hypot(init(2,0)-init(3,0), init(2,1)-init(3,1));
	rho=((init(0,1)+init(1,1))/2)-((init(2,1)+init(3,1))/2);
}


void Post_Process::calc_impact_gen() {   //calculates impact assuming vortices 1 and 2 form dipole and impact is calculated to (0,0)
	mat init;
	init.load("Vortex_xy_0.txt");
	di_length=hypot(init(0,0)-init(2,0),0);
	d = hypot(init(1,1)-init(0,1),0);
	rho=(init(1,1)+init(0,1))/2;
}


double Post_Process::get_di_length() {    // Gives private variables di_length...
	return di_length;
}

double Post_Process::get_impact() {      // ...and impact
	return rho;
}

void Post_Process::calc_theta() {       //calculates angle of vortex 1 from the horizontal 
	double delta_y;
	double delta_x;
	mat A;
	mat B;
	A.load("Vortex_xy_" + to_string(file_num - 2) + ".txt");
	B.load("Vortex_xy_" + to_string(file_num - 1) + ".txt");
	delta_y = B(0, 1) - A(0, 1);
	delta_x = B(0, 0) - A(0, 0);
	theta = atan2(delta_y, delta_x);
}

double Post_Process::get_theta() {                 //returns the calculated scattering angle 
	return theta;
}

void Post_Process::print_result() {                            //Prints the scattering angle to a file, with normalized impact parameter
	print.open("Scatter_Angles.txt", ofstream::app);
	print << rho / d << "\t" << theta << endl;
	print.close();
}

void Post_Process::print_result_4n() {                  //prints the scattering angle to a file in the 4N case, with L and angle of incidence
	print.open("Scatter_Angles.txt", ofstream::app);
	print << di_length/5 << "\t" << phi << "\t" << theta << endl;
	print.close();
}


void Post_Process::calc_separations(int file){      //calculates the vortex separations from a given file
	int counter=0;
	mat init;
	init.load("Vortex_xy_" + to_string(file) + ".txt");
	for (int i=0; i<N-1; i++){
		for (int j=1+i; j<N; j++){
			separations(counter)=hypot(init(i,0)-init(j,0), init(i,1)-init(j,1));
			counter++;
		}
	}
}



double Post_Process::get_separation(int i){     //return a known length 
	return separations(i);
}

vec Post_Process::get_separations(){
	return separations;
}



vec Post_Process::min_lengths() {     //return min lengths
	calc_separations(0); 
	vec out=separations;	
	for (int i = 1; i < file_num; ++i) {
		calc_separations(i);
		uvec minimums=find(separations<out);
		for (int j=0; j<minimums.n_elem; j++) out(minimums(j))=separations(minimums(j));
	}
	return out;
}

vec Post_Process::critical_lengths() {    //gives critical lengths in the 3N case
	calc_separations(0);	
	vec out = separations;
	for (int i = 1; i < file_num; ++i) {
		calc_separations(i);
		if (separations(0) < out(0)) {
			out(0) = separations(0);
			out(1) = separations(1);
			out(2) = separations(2);
		}
	}
	return out;
}

void Post_Process::print_min_lengths() {    //print minimum lengths to file 
	vec Lengths = min_lengths();
	print.open("Minimum_Lengths.txt", ofstream::app);
	//print << rho / d << "\t" << Lengths(0) / d << "\t" << Lengths(1) / d << "\t" << Lengths(2) / d << endl;
	for (int i=0; i<Lengths.n_elem; i++) print << "\t" << Lengths(0) / d;
	print <<endl;
	print.close();
}

void Post_Process::print_critical_lengths() {     //print the 3N critical lengths to file with the normed impact parameter 
	vec Lengths = critical_lengths();
	print.open("Scatter_Lengths.txt", ofstream::app);
	print << rho / d << "\t" << Lengths(0) / d << "\t" << Lengths(1) / d << "\t" << Lengths(2) / d << endl;
	print.close();
}



void Post_Process::file_delete_routine() {    //deletes all position files 
	int success = 0;
	for (int i = 0; i < file_num; ++i) {
		if (remove(("Vortex_xy_" + to_string(i) + ".txt").c_str()) != 0) success = 1;
	}
	if (flag_remove_hamil) remove("hamiltonian.txt");
	if (flag_remove_momentum) remove("momentums.txt");
	remove("python.txt");
	if (success == 0) cout << "\nFiles deleted successfully" << endl;
	else cout << "\nOne or more files could not be deleted" << endl;
}
