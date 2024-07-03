#include "Post_Process.hpp"


void::Post_Process::load_file(int file){
	temp_xy.load("Vortex_xy_"+std:: to_string(file)+".txt");
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

void Post_Process::find_circs(){
	vec init;
	init.load("python.txt");
	circs=init.tail(N);
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


void Post_Process::calc_theta(int vor) {       //calculates angle of vortex v from the horizontal 
	double delta_y;
	double delta_x;
	mat A;
	mat B;
	A.load("Vortex_xy_" + std:: to_string(file_num - 3) + ".txt");
	B.load("Vortex_xy_" + std:: to_string(file_num - 2) + ".txt");
	delta_y = B(vor, 1) - A(vor, 1);
	delta_x = B(vor, 0) - A(vor, 0);
	theta = atan2(delta_y, delta_x);
}


void Post_Process::dipole_cluster_populate_fix(int file)
{
	calc_mat_separations_periodic(file);
	
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			if (i != j)
			{
				if (circs(i) == -circs(j) && mat_separations(i, j) < dipole_tol)
				{ // populate dipole matrix if circulations are opposite and equal and vortices close enough
					dipoles(i, j) = 1;
					dipoles(j, i) = 1;
				}
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		while (arma::sum(dipoles.row(i)) > 1.0)
		{
			for (int j = 0; j < N-1; j++)
			{
				for (int k = j+1; k < N; k++)
				{
					if (dipoles(i,j) == 1 && dipoles(i,k) == 1)
					{						
						if (mat_separations(i, j) < mat_separations(i, k)){
							dipoles(i, k) = 0;
							dipoles(k,i) = 0;
						}
						else{
							dipoles(i, j) = 0;
							dipoles(j,i) = 0;
						}
					}
						
				}
			}
		}
	}

	for (int i = 0; i < N ; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j &&arma::sum(dipoles.row(i))== 0 &&arma::sum(dipoles.row(j))==0)
			{
				clusters(i, i) = 1;
				if (circs(i) == circs(j) && mat_separations(i, j) < cluster_tol  )
				{ // populate dipole matrix if circulations are opposite and equal and vortices close enough
					clusters(i, j) = 1;
					clusters(j, i) = 1;
				}
			}
		}
	}

	double irowsum = 0 ;
	double jrowsum = 0;

	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			if ((i!=j) && (clusters(i,j)==1))
			{
				irowsum =arma::sum(clusters.row(i));
				jrowsum =arma::sum(clusters.row(j));
				if (irowsum != jrowsum){
					if (irowsum > jrowsum) clusters.row(j) = clusters.row(i);
					else clusters.row(i) = clusters.row(j);
				}

			}
		}
	}
	
	
	rowvec temp_row;
	mat temp_mat= zeros(N,N);
	temp_mat.row(0) = clusters.row(0);
	bool matched;
	int count = 1;
	for (int i=1; i<N; i++){
		matched = false;
		temp_row = clusters.row(i);
		for (int j=0; j<i; j++){
			if (i!=j && arma::approx_equal(temp_row, clusters.row(j), "absdiff", 0.001)){
				matched = true;
	 		}
	 	}
		if (!matched){
			temp_mat.row(count) = temp_row;
			count++;
		}
	 }
	clusters = temp_mat;


	vec nonzeros = zeros(N);
	for (int i = 0; i< N; i++){
		if (arma::sum(clusters.row(i)) != 0){
			nonzeros(i) = 1;
		}
	}
	temp_mat = zeros(arma::sum(nonzeros), N);
	count=0;
	for (int i=0; i<N; i++){
		if (nonzeros(i) == 1){
			temp_mat.row(count) = clusters.row(i);
			count++;
		}
	}

	clusters = temp_mat;
}

















void Post_Process::dipole_cluster_populate(int file){  // populate the dipole and cluster matrices of pairs closer than cluster_tol and dipole_tol using separations in posittion file 
		calc_mat_separations(file);			
		for (int i=0; i<N-1; i++){
			for (int j=i+1; j<N; j++){  
				if (sign(circs(i)) == sign(circs(j)) && mat_separations(i,j) < cluster_tol){   //populate cluster matrix if circulations are same sign and vortices close enough   
					clusters(i,j)=1;
					clusters(j,i)=1;
				} 
				else if (circs(i) == -circs(j) && mat_separations(i,j) < dipole_tol ){  //populate dipole matrix if circulations are opposite and equal and vortices close enough  
					dipoles(i,j) = 1; 
					dipoles(j,i) = 1;
				}
			}
		}		
		vec cluster_unique_rows = vec(N, fill::ones);	
		vec dipole_unique_rows = vec(N, fill::ones);		
		for (int i=0; i<N; i++) if (accu(dipoles.row(i))==1) dipole_unique_rows(i) = 0;  
		for (int i=0; i<N-1; i++){
			for (int j=i+1; j<N; j++){
				if ((cluster_unique_rows(i) == 1) && (approx_equal(clusters.row(i), clusters.row(j), "absdiff", 0.01))) cluster_unique_rows(j) = 0;    //remove rows that are equal (i.e. remove storing twice)
				if ((dipole_unique_rows(i) == 1) && (approx_equal(dipoles.row(i), dipoles.row(j), "absdiff", 0.01))) dipole_unique_rows(j) = 0;
			}
		}
		clusters.shed_rows(find(cluster_unique_rows == 0));	
		dipoles.shed_rows(find(dipole_unique_rows == 0));	 		
		vec cluster_unique_rows2 = vec(clusters.n_rows , fill::ones);	 //now to make sure no clusters "share" vortices 	  
		vec dipole_unique_rows2 = vec(dipoles.n_rows , fill::ones);	 //now to make sure no dipoles "share" vortices 	  
		for (int i=0; i<N; i++){
			if (accu(clusters.col(i))!=1){    //if algorithm has operated correctly each row should sum to 1 unless a vortex is part of two clusters 	
				uvec sharedVortices = find(clusters.col(i) == 1);	
				for (int j=1; j < (int)sharedVortices.n_elem; j++){
					if (cluster_unique_rows2(sharedVortices(j))==1 ){
						clusters.row(sharedVortices(0))= conv_to<rowvec>::from(clusters.row(sharedVortices(0)) || clusters.row(sharedVortices(j)));
						cluster_unique_rows2(sharedVortices(j))=0;
					}
				}	
			}
			if (accu(dipoles.col(i))!=1){    //if algorithm has operated correctly each row should sum to 1 unless a vortex is part of two clusters 	
				uvec sharedVortices = find(dipoles.col(i) == 1);	
				for (int j=1; j < (int)sharedVortices.n_elem; j++){
					if (dipole_unique_rows2(sharedVortices(j))==1 ){
						dipoles.row(sharedVortices(0))= conv_to<rowvec>::from(dipoles.row(sharedVortices(0)) || dipoles.row(sharedVortices(j)));
						dipole_unique_rows2(sharedVortices(j))=0;
					}
				}	
			}
		}		
		clusters.shed_rows(find(cluster_unique_rows2 == 0));
		dipoles.shed_rows(find(dipole_unique_rows2 == 0));


/*
		//for (int i=0; i<N; i++) if ( accu(dipoles.row(i)) != 2 ) dipole_unique_rows(i) = 0;    //remove rows of dipoles that are not actually dipoles (i.e. passed over by the population loop, close to two vortices etc)   	
		for (int i=0; i<N; i++){	
			for (int j=0; j<N; j++){	
				if (((clusters(i,j) == 1) && (dipoles(i,j) == 1)) && ((cluster_unique_rows(i) == 1) && (dipole_unique_rows(i) == 1)))  //remove clusters that are actually part of dipoles 
					clusters(i,j) = 0;

			}
		}	
		clusters.shed_rows(find(cluster_unique_rows == 0));	
		dipoles.shed_rows(find(dipole_unique_rows == 0));	 
		
		   for (int i=0; i<N; i++) if ( accu(clusters.row(i)) == 0) cluster_unique_rows (i) = 0;    //remove rows of dipoles that are not actually dipoles (i.e. passed over by the population loop, close to two vortices etc) 	
		*/
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

void Post_Process::print_dipole_extremes(int vor1, int vor2){      //assumes dipoles start with vortex 1 and 2, with vortex 1 being antivortex 
	calc_mat_separations(0);	
	mat crit_length_matrix = mat_separations;
	int crit_file = 0; 	
	for (int i = 1; i < file_num; ++i) {
		calc_mat_separations(i); 
		if (mat_separations(vor1,vor2 ) < crit_length_matrix(vor1,vor2)) 
		{
			crit_length_matrix = mat_separations;
			crit_file = i;
		}	
	}
	preInterMin = 10.0;	//some dummy values to initialize these variables 
	preInterMax = 0.0;
	for (int j=0; j <= crit_file; j++){
		calc_mat_separations(j);
		if (mat_separations(0,1) < preInterMin) preInterMin = mat_separations(0,1);
		if (mat_separations(0,1) > preInterMax)	preInterMax = mat_separations(0,1); 
	}
	postInterMin = 10.0;
	postInterMax = 0.0;
	calc_mat_separations(file_num - 1);
	vec finalseps = mat_2_vec(mat_separations).head(N-1);	//for checking which vortex is part of final dipole
	int ind = 0; 	
	for (int x = 1; x < (int)circs.n_elem; x++) if ((circs(0) == -circs(x)) && (finalseps(x-1) < finalseps(ind)) ) ind = x;		
	if (ind == 0) ind++;
	for (int k = crit_file + 1; k < file_num; k++){
		calc_mat_separations(k);		
		if (mat_separations(0,ind) < postInterMin) postInterMin = mat_separations(0,ind);
		if (mat_separations(0,ind) > postInterMax) postInterMax = mat_separations(0,ind); 
	}	
}

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

void Post_Process::calc_mat_separations_periodic(int file)
{
	load_file(file);
	double dx = 0;
	double dy = 0;
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 1 + i; j < N; j++)
		{
			if (std::abs(temp_xy(i, 0) - temp_xy(j, 0)) > pi) 
			{
				dx = (2 * pi) - std::abs(temp_xy(i, 0) - temp_xy(j, 0));
			}
			else {
				dx = std::abs(temp_xy(i, 0) - temp_xy(j, 0));	
			}
			if (std::abs(temp_xy(i, 1) - temp_xy(j, 1)) > pi) 
			{
				dy = (2 * pi) - std::abs(temp_xy(i, 1) - temp_xy(j, 1));
			}
			else {
				dy = std::abs(temp_xy(i, 1) - temp_xy(j, 1));
			}
			mat_separations(i, j) = hypot(dx, dy);
			mat_separations(j, i) = mat_separations(i, j);
		}
	}
}

vec Post_Process::min_lengths() {     //return min lengths
	calc_mat_separations(0); 
	mat out=mat_separations;	
	for (int i = 1; i < file_num; ++i) {
		calc_mat_separations(i);
		uvec minimums=find(mat_separations<out);
		out(minimums)=mat_separations(minimums);
	}
	return out;
}

mat Post_Process::critical_lengths(int vor1, int vor2) {    //gives critical lengths when lij is at minimum 
	calc_mat_separations(0);	
	mat crit_length_matrix = mat_separations;
	for (int i = 1; i < file_num; ++i) {
		calc_mat_separations(i); 
		if (mat_separations(vor1,vor2 ) < crit_length_matrix(vor1,vor2)) crit_length_matrix = mat_separations;
	}
	return crit_length_matrix;
}

void Post_Process::exchange_find(){        //calculates if exchange scattering and total number of exchanges through system assuming vortex 1 is anti-vortex of initial dipole and vortex 2 is other dipole vortex 
	exchange = false;
	final_dipole_vortex = 1;
	calc_mat_separations(file_num - 1);	
	double dipole_size = mat_separations(0,1);
	for (int i=2; i<N; i++ ){
		if ((circs(0) == -circs(i)) && (mat_separations(0,i)<dipole_size)) {
			exchange = true;
			dipole_size = mat_separations(0,i);
			final_dipole_vortex = i;
		}
	}
	exchange_num = 0;
	int temp_dipole_vortex=1;	
	for (int i=0; i<file_num; i++){
		calc_mat_separations(i);		
		for (int j=0; j<N; j++){
			if ((temp_dipole_vortex != j) && (mat_separations(0,j) == min(mat_separations.row(0).tail(N-1)))) {
				temp_dipole_vortex = j;
				exchange_num++;
			}
		}
	}
}
//(circs(0) == -circs(j)) && 
void Post_Process::print_min_lengths() {    //print minimum lengths to file 
	vec Lengths = min_lengths();
	print.open("Minimum_Lengths.txt", std:: ofstream::app);
	print << rho / d;
	for (int i=0; i<(int)Lengths.n_elem; i++) print << "\t" << Lengths(i) / d;
	print <<endl;
	print.close();
}
/*
void Post_Process::print_critical_lengths(int v1, int v2) {     //print the 3N critical lengths to file with the normed impact parameter 
	vec Lengths = critical_lengths(v1,v2);
	print.open("Scatter_Lengths.txt", ofstream::app);
	print << rho / d << "\t" << Lengths(0) / d << "\t" << Lengths(1) / d << "\t" << Lengths(2) / d << endl;
	print.close();
}*/


void Post_Process::print_critical_lengths(int v1, int v2) {     //print the 3N critical lengths to file with the normed impact parameter 
	vec Lengths = mat_2_vec(critical_lengths(v1,v2));
	print.open("Scatter_Lengths.txt", std:: ofstream::app);
	//print << rho / d << "\t" << Lengths(0) / d << "\t" << Lengths(1) / d << "\t" << Lengths(2) / d << endl;
	print << rho/d;
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
	if (success == 0) cout << "\nFiles deleted successfully" << endl;
	else cout << "\nOne or more files could not be deleted" << endl;
}



vec Post_Process:: mat_2_vec(mat input){ 				//helper function to return a vector of separations from a separation matrix
	vec outVec=vec(N*(N-1)/2, fill::zeros );
	int counter=0;
	for (int i=0; i<N-1; i++){
		for (int j=i+1; j<N; j++){
			outVec(counter)=input(i,j);
			counter++;
		}
	}
	return outVec;
}
