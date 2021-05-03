
#include "Point_Vortex.hpp"
#include "Post_Process.hpp"


//-------------------------------------------------------------------- Misc. Functions --------------------------------------------------------------------------------------------------------


vec timesteps(double a, double b, int m) {	//Generate timestep mesh vector + Set h
	vec t_vals(m + 1, fill::zeros);
	h = (b - a) / m;
	t_vals[0] = a;
	for (int i = 1; i < m + 1; ++i) {
		t_vals[i] = t_vals[i - 1] + h;
	}
	return t_vals;
}

double x_per_sum(double x, double y){    //infinite sum for the equation of motion in the periodic case 
	double old_val = -sin(y)/(cosh(x)- cos(y));    //initialized for m=0
	double new_val;
	for (int m=1; m<=50; m++){
		new_val=old_val- (sin(y)/(cosh(x-(2*pi*m))-cos(y)))-(sin(y)/(cosh(x+(2*pi*m))-cos(y)));
		if (std::abs(new_val-old_val)< pow(10,-12)) return new_val;
		old_val=new_val;
	}
	return new_val;
}


double y_per_sum(double x, double y){    //infinite sum for the equation of motion in the periodic case 
	double old_val = sin(x)/(cosh(y)- cos(x));    //initialized for m=0
	double new_val;
	for (int m=1; m<=50; m++){
		new_val=old_val+ (sin(x)/(cosh(y-(2*pi*m))-cos(x)))+(sin(x)/(cosh(y+(2*pi*m))-cos(x)));
		if (abs(new_val-old_val)< pow(10,-12)) return new_val;
		old_val=new_val;
	}
	return new_val;
}


double per_hamil_sum(double x, double y) {		//infinite sum approximation for the periodic hamiltonian energy function
	double old_val = log(cosh(x)-cos(y));
	double new_val;
	for (int m = 1; m < 2000; m++) {
		new_val = old_val + log((cosh(x-(2*pi*m))-cos(y))/cosh(2*pi*m))+log((cosh(x+(2*pi*m))-cos(y))/cosh(-2*pi*m));
		if (std:: abs(new_val - old_val) < pow(10, -12)) return new_val;
		old_val = new_val;
	}
	return new_val;
}

/*mat position_generate(double L, double phi) {       //generates 4 position matrix with 2nd dipole at angle phi at length L
	mat out = { {5,4.9},{5,5.1},{0,0},{0,0} };
	rowvec Vor_3 = { (0.1*cos(phi - (pi / 2))) + (L*cos(phi)),(0.1*sin(phi - (pi / 2))) + (L*sin(phi))+5 };
	rowvec Vor_4 = { (0.1*cos(phi + (pi / 2))) + (L*cos(phi)),(0.1*sin(phi + (pi / 2))) + (L*sin(phi))+5 };
	out.row(2) = Vor_3;
	out.row(3) = Vor_4;
	return out;
}*/


mat position_generate(double L, double phi){      //generates dipole/dipole collision with L1=5 and a given L and incidence phi
	mat out = {{-5, 0.1}, {-5, -0.1}, {0,0}, {0,0}};
	double L2=5/L;
	rowvec Vor_3={(-L2*cos(phi))+(0.2/2*sin(phi)),(L2*sin(phi))+(0.2/2*cos(phi))};
	rowvec Vor_4={(-L2*cos(phi))-(0.2/2*sin(phi)), (L2*sin(phi))-(0.2/2*cos(phi))};
	out.row(2)= Vor_3;
	out.row(3)= Vor_4;
	return out;
}





mat rand_per_mat(){              //generates a random batrix inside a -pi->pi box for use in periodic case
	mat pos=randu(N,2);
	pos=(2*pi*pos)-pi;
	return pos; 
}

vec per_circ(){             //generates a circulation vector with N/2=-1 and N/2=1
	vec circ=ones<vec>(N);
	circ.head(N/2)*=-1;
	return circ;
}



//----------------------------------------------------------------- Lengths + Velocities --------------------------------------------------------------------------------------------------------
//Functions that we will use to solve the system

double xij(int i, int j, mat &xy) {  //calculate x distance between two vortices
	return xy(i, 0) - xy(j, 0);
}

double yij(int i, int j, mat &xy) { //calculate y distance between two vortices
	return xy(i, 1) - xy(j, 1);
}

double x_vel(int i, vec &circ, mat &xy) {		//calculate the x velocity of the ith vortex
	double sum = 0;
	for (int j = 0; j < N; ++j) {	
		if (i != j) { sum += (circ(j)*yij(i, j, xy)) / pow(hypot(xy(i,0)-xy(j,0), xy(i,1)-xy(j,1)), 2); }
	}
	return -sum / (2 * pi);
}

double per_x_vel(int i, vec &circ, mat &xy) {	//calculate the x velocity of ith vortex (periodic case)
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) sum += circ(j)*x_per_sum(xij(i,j,xy), yij(i, j, xy)); 
	}
	return sum;
}

double y_vel(int i, vec &circ, mat &xy) {   //calculate the y velocity of the ith vortex
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) { sum += (circ(j)*xij(i, j, xy)) / pow(hypot(xy(i,0)-xy(j,0), xy(i,1)-xy(j,1)), 2); }
	}
	return sum / (2 * pi);
}

double per_y_vel(int i, vec &circ, mat &xy) {	//calculate the y velocity of ith vortex (periodic case)
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) sum += circ(j)*y_per_sum(xij(i, j, xy), yij(i, j, xy)); 
	}
	return sum;
}

mat vel_matrix(mat &xy, vec &circ) {		//Creates matrix of vortex velocities at given positions
	mat vels(N, 2, fill::zeros);
	for (int i = 0; i < N; ++i) {
		vels(i, 0) = x_vel(i, circ, xy);
		vels(i, 1) = y_vel(i, circ, xy);
	}
	return vels;
}

mat per_vel_matrix(mat &xy, vec &circ) {   //Creates matrix of vortex velocities (periodic)
	mat vels(N, 2, fill::zeros);
	for (int i = 0; i < N; ++i) {
		vels(i, 0) = per_x_vel(i, circ, xy);
		vels(i, 1) = per_y_vel(i, circ, xy);
	}
	return vels;
}

//---------------------------------------------------------------------------- Conserved Quantities------------------------------------------------------------------------------------------------------------

double p(vec &circ, mat &xy) {			//linear momentum x direction
	double sum = 0;
	for (int i = 0; i < N; ++i) sum += circ(i)*xy(i, 0);
	return sum;
}

double q(vec &circ, mat &xy) {					//linear momentum y direction
	double sum = 0;
	for (int i = 0; i < N; ++i) sum += circ(i)*xy(i, 1);
	return sum;
}

double m(vec &circ, mat &xy) {				//angular momentum
	double sum = 0;
	for (int i = 0; i < N; ++i) sum += circ(i)*(pow(xy(i, 0),2) + pow(xy(i, 1),2));
	return sum;
}


double hamil(vec &circ, mat &xy) {		//Calculate the hamiltonian(energy) of the system, can be negative, periodic and infinite domain accounted for
	if (Flag_is_Periodic) {
		double sum = 0;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (i < j) sum += circ(i)*circ(j)*(per_hamil_sum(xij(i,j,xy),yij(i,j,xy))-pow(xij(i,j,xy),2)/(2*pi));
			}
		}
		return -sum;
	}
	else {
		double sum = 0;
		for (int i = 0; i < N-1; ++i) {
			for (int j = i+1; j < N; ++j) {	
				sum += circ(i)*circ(j)*log(hypot(xy(i,0)-xy(j,0), xy(i,1)-xy(j,1)));	
			}
		}
		return -sum / (2*pi);
	}
}

//------------------------------------------------------------------Dormand Price-------------------------------------------------------------------------------------------------------------
//Solves the sytstem and outputs as files  

double error_func(mat &RK1, mat &RK2){
	double sum=0;	
	double sc_x=0;
	double sc_y;
	for (int i=0; i<N; i++){
		sc_x=tol+(Rtol*fmax(RK1(i,0), RK2(i,0)));
		sc_y=tol+(Rtol*fmax(RK1(i,1), RK2(i,1)));
		sum+=pow((RK1(i,0)-RK2(i,0))/sc_x, 2);
		sum+=pow((RK1(i,1)-RK2(i,1))/sc_y, 2);
	}
	return sqrt(sum/(2*N));
}



cube rk45(mat &xy, vec &circ) {
	cube out(N, 2, 2, fill::zeros);			//initialize cube in order to return both solution matrices from the function
	if (Flag_is_Periodic) {			//periodic case, in a doubly periodic square box ( -pi<x<pi, -pi<y<pi
		mat k1 = h*per_vel_matrix(xy, circ);
		mat xy_temp = xy + (0.2*k1);
		mat k2 = h*per_vel_matrix(xy_temp, circ);
		xy_temp = xy + (0.075*k1) + (0.225*k2);
		mat k3 = h*per_vel_matrix(xy_temp, circ);
		xy_temp = xy + ((44.0 / 45)*k1) - ((56.0 / 15)*k2) + ((32.0 / 9)*k3);
		mat k4 = h*per_vel_matrix(xy_temp, circ);
		xy_temp = xy + ((19372.0 / 6561)*k1) - ((25360.0 / 2187)*k2) + ((64448.0 / 6561)*k3) - ((212.0 / 729)*k4);
		mat k5 = h*per_vel_matrix(xy_temp, circ);
		xy_temp = xy + ((9017.0 / 3168)*k1) - ((355.0 / 33)*k2) + ((46732.0 / 5247)*k3) + ((49.0 / 176)*k4) - ((5103.0 / 18656)*k5);
		mat k6 = h*per_vel_matrix(xy_temp, circ);
		xy_temp = xy + ((35.0 / 384)*k1) + ((500.0 / 1113)*k3) + ((125.0 / 192)*k4) - ((2187.0 / 6784)*k5) + ((11.0 / 84)*k6);
		mat k7 = h*per_vel_matrix(xy_temp, circ);
		out.slice(0) = xy_temp;												//runge kutta order four
		out.slice(1) = xy + ((5179.0 / 57600)*k1) + ((7571.0 / 16695)*k3) + ((393.0 / 640)*k4) - ((92097.0 / 339200)*k5) + ((187.0 / 2100)*k6) + ((1.0 / 40)*k7);	//runge kutta order 5
		while (abs(out).max()>pi) out[abs(out).index_max()]-=sign(out[abs(out).index_max()])*2*pi;           // performing mod operation on values of x and y greater than pi/less than -pi
		return out;
	}
	else {			//infinite domain
		mat k1 = h*vel_matrix(xy, circ);
		mat xy_temp = xy + (0.2*k1);
		mat k2 = h*vel_matrix(xy_temp, circ);
		xy_temp = xy + (0.075*k1) + (0.225*k2);
		mat k3 = h*vel_matrix(xy_temp, circ);
		xy_temp = xy + ((44.0 / 45)*k1) - ((56.0 / 15)*k2) + ((32.0 / 9)*k3);
		mat k4 = h*vel_matrix(xy_temp, circ);
		xy_temp = xy + ((19372.0 / 6561)*k1) - ((25360.0 / 2187)*k2) + ((64448.0 / 6561)*k3) - ((212.0 / 729)*k4);
		mat k5 = h*vel_matrix(xy_temp, circ);
		xy_temp = xy + ((9017.0 / 3168)*k1) - ((355.0 / 33)*k2) + ((46732.0 / 5247)*k3) + ((49.0 / 176)*k4) - ((5103.0 / 18656)*k5);
		mat k6 = h*vel_matrix(xy_temp, circ);
		xy_temp = xy + ((35.0 / 384)*k1) + ((500.0 / 1113)*k3) + ((125.0 / 192)*k4) - ((2187.0 / 6784)*k5) + ((11.0 / 84)*k6);
		mat k7 = h*vel_matrix(xy_temp, circ);
		out.slice(0) = xy_temp;
		out.slice(1) = xy + ((5179.0 / 57600)*k1) + ((7571.0 / 16695)*k3) + ((393.0 / 640)*k4) - ((92097.0 / 339200)*k5) + ((187.0 / 2100)*k6) + ((1.0 / 40)*k7);
		return out;
	}
}

void dormand_price(vec &t_steps, vec &circ, mat &xy) {
	std::ofstream writer, hamil_print, momen_print;		//initialize the three file output objects
	double error=0;
	double h_old = h;
	if (Flag_Print_Hamil){	
		hamil_print.open("hamiltonian.txt");
		hamil_print.precision(12);
	}
	if (Flag_Print_Momen){	
		momen_print.open("momentums.txt");
		momen_print.precision(12);
	}
	writer.precision(12);	
	writer.open("Vortex_xy_0.txt");	//open initial output file	
	xy.raw_print(writer);	//print initial condition matrix
	writer.close();	
	double end = t_steps[t_steps.size() - 1];	//set when to end print process (when end time is reached)
	cube RK_Results(N, 2, 2);			//to store rk results
	int file_counter = 1;				//count the files that have already been made
	int run_counter = 1;
	for (double t = t_steps[1]; t <= end; t += h) {
		if (Flag_Print_Hamil) hamil_print << hamil(circ, xy) << endl;
		if (Flag_Print_Momen) momen_print << p(circ, xy) << "\t" << q(circ, xy) << "\t" << m(circ, xy) << endl;	
		RK_Results = rk45(xy, circ);
		error=error_func(RK_Results.slice(0), RK_Results.slice(1));
		while (error>1){
			h*=std:: min(hfacmax,std:: max(hfacmin,hfac*pow(1/error,0.2)));
			if (h<hmin){
				h=hmin;
				RK_Results = rk45(xy, circ);
				break;
			}
			RK_Results = rk45(xy, circ);
			error=error_func(RK_Results.slice(0), RK_Results.slice(1));
		}	
		xy = RK_Results.slice(1); //set xy
		if (run_counter%write_limit==0 && file_counter <File_Max){	
			writer.open("Vortex_xy_" + std:: to_string(file_counter) + ".txt");	
			xy.raw_print(writer);
			writer.close();
			file_counter++;
		}
		else if (file_counter==File_Max) break;
		run_counter++;	
		cout <<t <<"\t" << file_counter<< endl;
	}
	writer.open("Vortex_xy_" + std:: to_string(file_counter) + ".txt");
	xy.raw_print(writer);
	writer.close(); 
	file_counter++;
	writer.open("python.txt");		//print circulation and file number for plotting
	writer << file_counter << endl << circ << "\n\n";
	writer.close();
	if (Flag_Print_Momen) momen_print.close();
	if (Flag_Print_Hamil) hamil_print.close();	
	std::cout << "File output complete\n"<<run_counter;
	h = h_old;
}

void bound_dormand_price(vec &circ, mat &xy) {      //run dormand price until final dipole is 20d + L away
  std:: ofstream writer, hamil_print, momen_print;			//initialize the three file output objects
  double error=0;
  double h_old = h;  
  if (Flag_Print_Hamil) { 
	  hamil_print.open("hamiltonian.txt");
	  hamil_print.precision(12);
  }
  if (Flag_Print_Momen) {
	  momen_print.open("momentums.txt");
	  momen_print.precision(12);
  }
  writer.precision(12); 
  writer.open("Vortex_xy_0.txt");	//open initial output file
  xy.raw_print(writer);	//print initial condition matrix
  writer.close(); 
  cube RK_Results(N, 2, 2);			//to store rk result
  int file_counter=1;
  int run_counter = 1;
  double init_d=hypot(0,xy(1,1)-xy(0,1));
  double init_L=hypot(xy(0,0)-xy(2,0), xy(0,1)-xy(2,1));
  vec lengths={hypot(xy(0,0)-xy(2,0), xy(0,1)-xy(2,1)),hypot(xy(0,0)-xy(1,0), xy(0,1)-xy(1,1))};
  while(lengths.max()<(20*init_d)+init_L) {
        RK_Results = rk45(xy, circ); 
    error = error_func(RK_Results.slice(0), RK_Results.slice(1));
    if (error < pow(hfac,5)) h*=std::min(hfacmax,hfac*pow(1/error,0.2)); 
    if (h>hmax) h=hmax;
   while (error>1){		
	if (h<hmin){
		h=hmin;
		RK_Results = rk45(xy, circ);
		break;
	}
	RK_Results = rk45(xy, circ);
	h*=std::min(hfacmax,std::max(hfacmin,hfac*pow(1/error,0.2)));
	error=error_func(RK_Results.slice(0), RK_Results.slice(1));
    }
    xy = RK_Results.slice(1); //set xy, print xy then repeat
    if (run_counter%write_limit==0&& file_counter<File_Max){
	    writer.open("Vortex_xy_" + std:: to_string(file_counter) + ".txt");
	    xy.raw_print(writer);
	    writer.close(); 
	    file_counter++;
	    if (Flag_Print_Hamil) hamil_print <<hamil(circ, xy) << endl;
	    if (Flag_Print_Momen) momen_print << p(circ, xy) << "\t" << q(circ, xy) << "\t" << m(circ, xy) << endl;
    }
    else if (file_counter==File_Max) break;
    lengths={hypot(xy(0,0)-xy(2,0), xy(0,1)-xy(2,1)),hypot(xy(0,0)-xy(1,0), xy(0,1)-xy(1,1))};
    run_counter++;
  }
  writer.open("Vortex_xy_" + std:: to_string(file_counter) + ".txt");
  xy.raw_print(writer);
  writer.close(); 
  file_counter++;
  writer.open("python.txt");		//print circulation and file number for plotting
  writer <<file_counter<<endl <<circ << "\n\n";
  writer.close();
  if (Flag_Print_Momen) momen_print.close();
  if (Flag_Print_Hamil) hamil_print.close(); 
  cout << "File output complete\n";
  h = h_old;
}


//------------------------------------------------------------------------------------------------------------------

int main(){
	mat xy;	
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 1.2/320.0;
	double psi = 0;
	double psi_inc = 2*pi / 200;
	double L=0.4;
	for (int i = 0; i < 200; ++i) {
		psi += psi_inc;	
		for (int j = 0; j < 320; ++j) {
			L = 0.4+(L_inc*j);
			xy = position_generate(L, psi);
			bound_dormand_price(circ, xy);
			post.update_4n(psi, L);
			post.count_files();
			post.print_critical_lengths_4n(0,2);
			cout << i << "\t"<< j <<endl;
			post.file_delete_routine();
			}
	}
	return 0;
}

	/*mat xy;
	vec circ={-1,1,1,1, 1, 1};
	vec seps;
	double dipole_end_d;
	double cluster_end_d;
	double avg_x;
	double avg_y;
	vec avg_seps;
	Post_Process post(0);
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	cout<< "matric+circ initialized";
	for (int i=0; i<=320; i++){
		for (int j=0;j<=200;j++){ 	
			xy={{-7.0-(j*pi/12.0*0.2/200),-0.26+(i*0.22/320)}, {-7.0-(j*pi/12.0*0.2/200),-0.06+(i*0.22/320)}, {0.0,0.1}, {0.1,0.0}, {0.0,-0.1}, {-0.1,0.0}};
			bound_dormand_price(circ,xy);
			post.count_files();
			post.calc_impact_gen();
			post.calc_separations(post.get_file_num()-1); 	
			seps=post.get_separations();
	if (seps.head(5).index_min()==0){
			dipole_end_d=seps(0);
			avg_x=(xy(2,0)+xy(3,0)+xy(4,0)+xy(5,0))/4;  //avg of vortices 3456
			avg_y=(xy(2,1)+xy(3,1)+xy(4,1)+xy(5,1))/4;
			avg_seps={hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y), hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 3456  from avg point
			cluster_end_d=avg_seps.max()*2;
		}
		if (seps.head(5).index_min()==1){
			dipole_end_d=seps(1);
			avg_x=(xy(1,0)+xy(3,0)+xy(4,0)+xy(5,0))/4;  //avg of vortices 2456
			avg_y=(xy(1,1)+xy(3,1)+xy(4,1)+xy(5,1))/4;  //avg of vortices 2456
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y), hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 2456  from avg point
			cluster_end_d=avg_seps.max()*2;
		}
		if (seps.head(5).index_min()==2){
			dipole_end_d=seps(2);
			avg_x=(xy(1,0)+xy(2,0)+xy(4,0)+xy(5,0))/4;  //avg of vortices 2356
			avg_y=(xy(1,1)+xy(2,1)+xy(4,1)+xy(5,1))/4;  //avg of vortices 2356
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y),hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 2356  from avg point
			cluster_end_d=avg_seps.max()*2;
		}
		if (seps.head(5).index_min()==3){
			dipole_end_d=seps(3);
			avg_x=(xy(1,0)+xy(2,0)+xy(3,0)+xy(5,0))/4;  //avg of vortices 2346
			avg_y=(xy(1,1)+xy(2,1)+xy(3,1)+xy(5,1))/4;  //avg of vortices 2346
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y), hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 2346  from avg point
			cluster_end_d=avg_seps.max()*2;
		}
		if (seps.head(5).index_min()==4){
			dipole_end_d=seps(4);
			avg_x=(xy(1,0)+xy(2,0)+xy(3,0)+xy(4,0))/4;  //avg of vortices 2345
			avg_y=(xy(1,1)+xy(2,1)+xy(3,1)+xy(4,1))/4;  //avg of vortices 2345
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y), hypot(xy(4,0)-avg_x, xy(4,1)-avg_y)};                     //separations of vortices 2345  from avg point
			cluster_end_d=avg_seps.max()*2;
		}		
		output<<post.get_di_length()/0.2<<"\t"<<post.get_impact()/0.2<<"\t"<<dipole_end_d<<"\t"<<cluster_end_d<<endl;
			post.file_delete_routine();	
			cout<<i<<"\t"<<j<<"\t"<<endl; }
		} 
	output.close();			
	return 0;
	*/


/*mat xy;
	vec circ={-1,1,1, 1, 1};	
	Post_Process post(0);
	ofstream output;
	output.open("Scattering_Angles.txt",ofstream::app);	
	double phase1_ang, phase2_ang, phase3_ang;
	for (int i=0; i<=100; i++){
		xy={{-15.0,-4.1+(i*10.0/100.0)}, {-15.0,-3.9+(i*10.0/100.0)},{0.0,0.1}, {sqrt(3.0)*0.05,-0.05}, {-sqrt(3.0)*0.05, -0.05}};
		bound_dormand_price(circ,xy);
		post.count_files();
		post.calc_impact_gen();
		post.calc_theta(); 
		phase1_ang=post.get_theta();
		post.file_delete_routine();	
		xy={{-15.0-(0.2*pi/18.0),-4.1+(i*10.0/100.0)}, {-15.0-(0.2*pi/18.0),-3.9+(i*10.0/100.0)},{0.0,0.1}, {sqrt(3.0)*0.05,-0.05}, {-sqrt(3.0)*0.05, -0.05}};
		bound_dormand_price(circ,xy);
		post.count_files();
		post.calc_impact_gen();
		post.calc_theta(); 
		phase2_ang=post.get_theta();
		post.file_delete_routine();	
		xy={{-15.0-(0.2*pi/9.0),-4.1+(i*10.0/100.0)}, {-15.0-(0.2*pi/9.0),-3.9+(i*10.0/100.0)},{0.0,0.1}, {sqrt(3.0)*0.05,-0.05}, {-sqrt(3.0)*0.05, -0.05}} ;
		bound_dormand_price(circ,xy);
		post.count_files();
		post.calc_impact_gen();
		post.calc_theta(); 
		phase3_ang=post.get_theta();
		post.file_delete_routine();
		output<<post.get_impact()/0.2<<"\t"<< phase1_ang <<"\t"<< phase2_ang<< "\t"<< phase3_ang<< endl;
	} 
	output.close();			
	return 0;
*/


/*
chdir("/home/karl/Documents/source/data/2cluster/L=10.05, rho=0.1/");
	Post_Process post(0);
	post.count_files();	
	ofstream out;
	out.open("dipole-cluster growth.txt", ofstream::app);
	double dipole_size;
	double cluster_size;
	vec lengths;
	for (int i=0; i<post.get_file_num(); i++){
		post.calc_separations(i);
		lengths=post.get_separations();	
		if (lengths.head(3).index_min()==0){
			dipole_size=lengths(0);
			cluster_size=lengths(5);
		}
		if (lengths.head(3).index_min()==1){
			dipole_size=lengths(1);
			cluster_size=lengths(4);
		}
		if (lengths.head(3).index_min()==2){
			dipole_size=lengths(2);
			cluster_size=lengths(3);
		}
		out<<dipole_size<<"\t"<<cluster_size<<endl;
	}
	out.close();
	return 0;
}
/*
chdir("/home/karl/Documents/source/data/2cluster/L=10.05, rho=0.1/");
	Post_Process post(0);
	post.count_files();	
	ofstream out;
	out.open("dipole-cluster growth.txt", ofstream::app);
	double dipole_size;
	double cluster_size;
	vec lengths;
	for (int i=0; i<post.get_file_num(); i++){
		post.calc_separations(i);
		lengths=post.get_separations();	
		if (lengths.head(3).index_min()==0){
			dipole_size=lengths(0);
			cluster_size=lengths(5);
		}
		if (lengths.head(3).index_min()==1){
			dipole_size=lengths(1);
			cluster_size=lengths(4);
		}
		if (lengths.head(3).index_min()==2){
			dipole_size=lengths(2);
			cluster_size=lengths(3);
		}
		out<<dipole_size<<"\t"<<cluster_size<<endl;
	}
	out.close();
	return 0;

	}


/*
chdir("/home/karl/Documents/source/data/4cluster/L=4.06, rho=0.02/");
	Post_Process post(0);
	post.count_files();	
	ofstream out;
	out.open("dipole-cluster growth.txt", ofstream::app);
	double dipole_size;
	double cluster_size;
	double avg_x;
	double avg_y;
	vec avg_seps;
	vec lengths;
	mat xy;
	for (int i=0; i<post.get_file_num(); i++){	
		xy.load("Vortex_xy_" + to_string(i) + ".txt");
		post.calc_separations(i);
		lengths=post.get_separations();	
		if (lengths.head(5).index_min()==0){
			dipole_size=lengths(0);
			avg_x=(xy(2,0)+xy(3,0)+xy(4,0)+xy(5,0))/4;  //avg of vortices 3456
			avg_y=(xy(2,1)+xy(3,1)+xy(4,1)+xy(5,1))/4;
			avg_seps={hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y), hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 3456  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(5).index_min()==1){
			dipole_size=lengths(1);
			avg_x=(xy(1,0)+xy(3,0)+xy(4,0)+xy(5,0))/4;  //avg of vortices 2456
			avg_y=(xy(1,1)+xy(3,1)+xy(4,1)+xy(5,1))/4;  //avg of vortices 2456
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y), hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 2456  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(5).index_min()==2){
			dipole_size=lengths(2);
			avg_x=(xy(1,0)+xy(2,0)+xy(4,0)+xy(5,0))/4;  //avg of vortices 2356
			avg_y=(xy(1,1)+xy(2,1)+xy(4,1)+xy(5,1))/4;  //avg of vortices 2356
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y),hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 2356  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(5).index_min()==3){
			dipole_size=lengths(3);
			avg_x=(xy(1,0)+xy(2,0)+xy(3,0)+xy(5,0))/4;  //avg of vortices 2346
			avg_y=(xy(1,1)+xy(2,1)+xy(3,1)+xy(5,1))/4;  //avg of vortices 2346
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y), hypot(xy(5,0)-avg_x, xy(5,1)-avg_y)};                     //separations of vortices 2346  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(5).index_min()==4){
			dipole_size=lengths(4);
			avg_x=(xy(1,0)+xy(2,0)+xy(3,0)+xy(4,0))/4;  //avg of vortices 2345
			avg_y=(xy(1,1)+xy(2,1)+xy(3,1)+xy(4,1))/4;  //avg of vortices 2345
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y), hypot(xy(4,0)-avg_x, xy(4,1)-avg_y)};                     //separations of vortices 2345  from avg point
			cluster_size=avg_seps.max()*2;
		}
		out<<dipole_size<<"\t"<<cluster_size<<endl;
	}
	out.close();
	return 0;
*/

/* chdir("/home/karl/Documents/source/data/3cluster/L=4.08, rho=0.02/");
	Post_Process post(0);
	post.count_files();	
	ofstream out;
	out.open("dipole-cluster growth.txt", ofstream::app);
	double dipole_size;
	double cluster_size;
	double avg_x;
	double avg_y;
	vec avg_seps;
	vec lengths;
	mat xy;
	for (int i=0; i<post.get_file_num(); i++){	
		xy.load("Vortex_xy_" + to_string(i) + ".txt");
		post.calc_separations(i);
		lengths=post.get_separations();	
		if (lengths.head(4).index_min()==0){
			dipole_size=lengths(0);
			avg_x=(xy(2,0)+xy(3,0)+xy(4,0))/3;  //avg of vortices 345
			avg_y=(xy(2,1)+xy(3,1)+xy(4,1))/3;
			avg_seps={hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y)};                     //separations of vortices 345  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(4).index_min()==1){
			dipole_size=lengths(1);
			avg_x=(xy(1,0)+xy(3,0)+xy(4,0))/3;  //avg of vortices 245
			avg_y=(xy(1,1)+xy(3,1)+xy(4,1))/3;  //avg of vortices 245
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y)};                     //separations of vortices 245  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(4).index_min()==2){
			dipole_size=lengths(2);
			avg_x=(xy(1,0)+xy(2,0)+xy(4,0))/3;  //avg of vortices 235
			avg_y=(xy(1,1)+xy(2,1)+xy(4,1))/3;  //avg of vortices 235
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(4,0)-avg_x, xy(4,1)-avg_y)};                     //separations of vortices 235  from avg point
			cluster_size=avg_seps.max()*2;
		}
		if (lengths.head(4).index_min()==3){
			dipole_size=lengths(3);
			avg_x=(xy(1,0)+xy(2,0)+xy(3,0))/3;  //avg of vortices 234
			avg_y=(xy(1,1)+xy(2,1)+xy(3,1))/3;  //avg of vortices 234
			avg_seps={hypot(xy(1,0)-avg_x, xy(1,1)-avg_y),hypot(xy(2,0)-avg_x, xy(2,1)-avg_y),hypot(xy(3,0)-avg_x, xy(3,1)-avg_y)};                     //separations of vortices 234  from avg point
			cluster_size=avg_seps.max()*2;
		}
		out<<dipole_size<<"\t"<<cluster_size<<endl;
	}
	out.close();
	return 0;
*/
/*
chdir("/home/karl/Documents/source/data/2cluster/L=10.05, rho=0.1/");
	Post_Process post(0);
	post.count_files();	
	ofstream out;
	out.open("dipole-cluster growth.txt", ofstream::app);
	double dipole_size;
	double cluster_size;
	vec lengths;
	for (int i=0; i<post.get_file_num(); i++){
		post.calc_separations(i);
		lengths=post.get_separations();	
		if (lengths.head(3).index_min()==0){
			dipole_size=lengths(0);
			cluster_size=lengths(5);
		}
		if (lengths.head(3).index_min()==1){
			dipole_size=lengths(1);
			cluster_size=lengths(4);
		}
		if (lengths.head(3).index_min()==2){
			dipole_size=lengths(2);
			cluster_size=lengths(3);
		}
		out<<dipole_size<<"\t"<<cluster_size<<endl;
	}
	out.close();
	return 0;
*/


/*mat xy;
	vec circ={-1,1,1,1, 1,1};
	double dipole_end_d;
	Post_Process post(0);
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	cout<< "matric+circ initialized";
	for (int i=0; i<=320; i++){
		for (int j=0;j<=200;j++){ 	
			xy={{-4.0-(j*pi/4.0*0.2/200),-0.26+(i*0.0006875)}, {-4.0-(j*pi/4.0*0.2/200),-0.06+(i*0.0006875)}, {0.0,0.1}, {0.1,0.0}, {0.0,-0.1}, {-0.1,0.0}}; 
			bound_dormand_price(circ,xy);
			post.count_files();
			post.calc_impact_gen();
			post.calc_separations(post.get_file_num()-1); 
			dipole_end_d=post.get_separations().head(N-1).min(); 
			output<<post.get_di_length()/0.2<<"\t"<<post.get_impact()/0.2<<"\t"<<dipole_end_d<<"\t"<<endl;
			post.file_delete_routine();	
			cout<<i<<"\t"<<j<<"\t"<<endl;
			}
		} 
	output.close();			
	return 0;
*/

/*mat xy;
	vec circ={-1,1,1,1, 1, 1};	
	Post_Process post(0);
	ofstream output;
	output.open("Scattering_Angles.txt",ofstream::app);	
	double phase1_ang, phase2_ang, phase3_ang;
	for (int i=0; i<=100; i++){
		xy={{-15.0,-4.1+(i*10.0/100.0)}, {-15.0,-3.9+(i*10.0/100.0)}, {0.0,0.1},{0.1,0.0}, {0.0,-0.1}, {-0.1,0.0}}; 
		bound_dormand_price(circ,xy);
		post.count_files();
		post.calc_impact_gen();
		post.calc_theta(); 
		phase1_ang=post.get_theta();
		post.file_delete_routine();	
		xy={{-15.0-(0.2*pi/24.0),-4.1+(i*10.0/100.0)}, {-15.0-(0.2*pi/24.0),-3.9+(i*10.0/100.0)}, {0.0,0.1},{0.1,0.0}, {0.0,-0.1}, {-0.1,0.0}}; 
		bound_dormand_price(circ,xy);
		post.count_files();
		post.calc_impact_gen();
		post.calc_theta(); 
		phase2_ang=post.get_theta();
		post.file_delete_routine();	
		xy={{-15.0-(0.2*pi/12.0),-4.1+(i*10.0/100.0)}, {-15.0-(0.2*pi/12.0),-3.9+(i*10.0/100.0)}, {0.0,0.1},{0.1,0.0}, {0.0,-0.1}, {-0.1,0.0}}; 
		bound_dormand_price(circ,xy);
		post.count_files();
		post.calc_impact_gen();
		post.calc_theta(); 
		phase3_ang=post.get_theta();
		post.file_delete_routine();
		output<<post.get_impact()/0.2<<"\t"<< phase1_ang <<"\t"<< phase2_ang<< "\t"<< phase3_ang<< endl;
	} 
	output.close();			
	return 0;
*/
/*
mat xy;
	vec circ={-1,1,1,1, 1};
	double dipole_end_d;
	Post_Process post(0);
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	cout<< "matric+circ initialized";
	for (int i=0; i<=320; i++){
		for (int j=0;j<=200;j++){ 	
			xy={{-4.0-(j*pi/3.0*0.2/200),-0.34+(i*0.0010625)}, {-4.0-(j*pi/3.0*0.2/200),-0.14+(i*0.0010625)}, {0.0,0.1}, {sqrt(3.0)*0.05,-0.05}, {-sqrt(3.0)*0.05, -0.05}}; 
			bound_dormand_price(circ,xy);
			post.count_files();
			post.calc_impact_gen();
			post.calc_separations(post.get_file_num()-1); 
			dipole_end_d=post.get_separations().head(N-1).min(); 
			output<<post.get_di_length()/0.2<<"\t"<<post.get_impact()/0.2<<"\t"<<dipole_end_d<<"\t"<<endl;
			post.file_delete_routine();	
			cout<<i<<"\t"<<j<<"\t"<<endl;
			}
		} 
	output.close();			
	return 0;
*/






/*mat xy;
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 3.0/1600.0;
	double psi = 0;
	double psi_inc = pi/200.0;
	double L=0.4;	
	xy = position_generate(L, psi);
	bound_dormand_price(circ, xy);	
	for (int i = 1; i <= 400; ++i) {
		psi = (i*psi_inc);	
		for (int j = 0; j <= 640; ++j) {
			L = 0.4+(j*L_inc);
			xy = position_generate(L, psi);
			bound_dormand_price(circ, xy);	
			post.update_4n(psi, L);
			post.count_files();	
			post.calc_theta();
			post.print_result_4n();
			cout << i<<"\t" << j << "\t" << post.get_theta() << endl;
			post.file_delete_routine();
		}
	}	
	return 0;
*/
/*mat xy;
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 3.0/800.0;
	double psi = 0;
	double psi_inc = pi/100;
	double L;	
	for (int i = 1; i <= 200; ++i) {
		psi = (i*psi_inc);	
		for (int j = 0; j <= 320; ++j) {
			L = 0.4+(j*L_inc);
			post.update_4n(psi, L);
			xy = position_generate(L, psi);
			bound_dormand_price(circ, xy);	
			post.count_files();	
			post.calc_theta();
			post.print_result_4n();
			cout << i<<"\t" << j << endl;
			post.file_delete_routine();
		}
	}	
	return 0;
*/

/*mat xy;
	vec circ={-1,1,1,1, 1};
	double dipole_end_d;
	Post_Process post(0);
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	cout<< "matric+circ initialized";
	for (int i=0; i<=80; i++){
		for (int j=0;j<=50;j++){ 	
			xy={{-4.0-(j*pi/3.0*0.2/50),-0.34+(i*0.00425)}, {-4.0-(j*pi/3.0*0.2/50),-0.14+(i*0.00425)}, {0.0,0.1}, {sqrt(3.0)*0.05,-0.05}, {-sqrt(3.0)*0.05, -0.05}}; 
			bound_dormand_price(circ,xy);
			post.count_files();
			post.calc_impact_gen();
			post.calc_separations(post.get_file_num()-1); 
			dipole_end_d=post.get_separations().head(N-1).min(); 
			output<<post.get_di_length()/0.2<<"\t"<<post.get_impact()/0.2<<"\t"<<dipole_end_d<<"\t"<<endl;
			post.file_delete_routine();	
			cout<<i<<"\t"<<j<<"\t"<<endl;
			}
		} 
	output.close();			
	return 0;
*/



/*
	mat xy;
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 3.0/800.0;
	double psi = 0;
	double psi_inc = pi/100;
	double L;
	vec minimums;
	ofstream output;
	output.open("Minimum_Lengths.txt",ofstream::app);
	output.precision(12);
	for (int i = 1; i <= 200; ++i) {
		psi = (i*psi_inc);	
		for (int j = 0; j <= 320; ++j) {
			L = 0.4+(j*L_inc);
			xy = position_generate(L, psi);
			bound_dormand_price(circ, xy);
			post.update_4n(psi, L);
			post.count_files();	
			minimums = post.min_lengths();
			output<<post.get_di_length()<<"\t"<<post.get_psi()<<"\t" << minimums(0)<< "\t" << minimums(1)<< "\t"<< minimums(2)<< "\t"<< minimums(3)<< "\t"<< minimums(4)<< "\t"<< minimums(5)<< endl;
			cout << i<<"\t" << j << endl;
			post.file_delete_routine();
			}
	}
	output.close();
	return 0;
}
*/
/*
mat xy;
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 0.001875;
	double psi = 0;
	double psi_inc = pi/200.0;
	double L;
	double dipole1;
	double dipole2;
	ofstream output;
	output.open("region_plot.txt",ofstream::app);
	for (int i = 1; i <= 400; ++i) {
		psi = (i*psi_inc);	
		for (int j = 0; j <= 640; ++j) {
			L = 0.4+(j*L_inc);
			xy = position_generate(L, psi);
			bound_dormand_price(circ, xy);
			post.update_4n(psi, L);
			post.count_files();	
			post.calc_separations(post.get_file_num()-1);
			if (post.get_separations()(0)<post.get_separations()(2)) output<<post.get_di_length()<<"\t"<<post.get_psi()<<"\t"<< 0 << endl;
			else output<<post.get_di_length()<<"\t"<<post.get_psi()<<"\t"<< 1 << endl;	
			cout << i<<"\t" << j << endl;
			post.file_delete_routine();
			}
	}
	output.close();
	return 0;
*/
/*int main(){ 
	mat xy;
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 0.001875;
	double psi = 0;
	double psi_inc = 0.01554435;
	double L;
	double dipole1;
	double dipole2;
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	for (int i = 1; i <= 200; ++i) {
		psi = (i*psi_inc);	
		for (int j = 0; j <= 320; ++j) {
			L = 0.4+(j*L_inc);
			xy = position_generate(L, psi);
			bound_dormand_price(circ, xy);
			post.update_4n(psi, L);
			post.count_files();	
			post.calc_separations(post.get_file_num()-1);
			dipole1=post.get_separations().head(N-1).min(); 
			dipole2=post.get_separations().tail(N-1).min(); 
			output<<post.get_di_length()<<"\t"<<post.get_psi()<<"\t"<<dipole1<<"\t"<<dipole2<<endl;
			cout << i<<"\t" << j << endl;
			post.file_delete_routine();
			}
	}
	output.close();
	return 0;
}
*/
/*mat xy;
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	ofstream writer;
	double L;
	double a=0.4;
	double b=0.99;
	double di_sc;
	double ex_sc;	
	double mid;
	writer.open("interpolation_points.txt", ofstream::app);
	for (int i=0; i<20; i++){
		di_sc=0.05;
		ex_sc=pi;	
		mid=(di_sc+ex_sc)/2;
		L=(0.5*(a+b))+(0.5*(b-a)*cos(((2*i)+1)*pi/(2*20)));
		while (abs(di_sc-ex_sc)>pow(10,-5)){	
			xy = position_generate(L, mid);
			bound_dormand_price(circ, xy);
			post.count_files();	
			post.calc_separations(post.get_file_num()-1);
			if (post.get_separations()(0)<post.get_separations()(2)) di_sc=mid;
			else ex_sc=mid;
			post.file_delete_routine();
			mid=(di_sc+ex_sc)/2;
		}
		writer << L << "\t" << mid <<endl;
		cout<< "run "<<i<<" complete"<<endl;
	}
	writer.close();	
	return 0;	*/

/*
mat xy; 	
	vec circ = { 1,-1,1,-1 };	
	double phi = 2.95535;
	double L = 0.7;
	xy=position_generate(L, phi);
	system("mkdir 'phi=2.95535,L=0.7'");
	chdir("phi=2.95535,L=0.7");
	bound_dormand_price(circ,xy);			
	return 0;

*/
/*	mat xy; 
	vec t = timesteps(0, 15, 100);
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 0.625;
	double phi_inc = pi / 6;
	double phi = 0;
	double L;
	for (int i = 1; i < 12; ++i) {
		phi += phi_inc;
		L = 1.875;
		for (int j = 1; j <= 10; ++j) {
			L += L_inc;
			post.update_4N(phi, L);
			xy = position_generate(L, phi);
			bound_dormand_price(circ, xy);
			post.Count_Files();
			post.Calc_Theta();
			post.Print_Result_4N();
			post.File_Delete_Routine();
		}*/



/*      periodic random simulation 
mat xy=rand_per_mat();
	vec circ=per_circ();	
	vec t=timesteps(0,700,2000);
	dormand_price(t, circ, xy);
	return 0;
*/


/*mat xy;
	vec circ={-1,1,1,1,1, 1, 1};
	double dipole_end_d;
	Post_Process post(0);
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	cout<< "matric+circ initialized";
	for (int i=0; i<=80; i++){
		for (int j=0;j<=160;j++){ 	
			if ((i!=5)||(j!=5)){
			xy={{-10.0-(j*pi*0.00125),-0.2+(i*0.001875)}, {-10.0-(j*pi*0.00125),0.0+(i*0.001875)}, {0,0.1},{0.1, 0.0}, {0, -0.1}, {-0.1,0}, {0,0}}; 
			bound_dormand_price(circ,xy);
			post.count_files();
			post.calc_impact_gen();
			post.calc_separations(post.get_file_num()-1); 
			dipole_end_d=post.get_separations().head(N-1).min(); 
			output<<post.get_di_length()/0.2<<"\t"<<post.get_impact()/0.2<<"\t"<<dipole_end_d<<"\t"<<endl;
			post.file_delete_routine();	
			cout<<i<<"\t"<<j<<"\t"<<endl;
			}
		}
	  }
	output.close();			
	return 0;
	Post_Process post(0);
	post.count_files();
	post.file_delete_routine();
	return 0;*/

/*
double inc = 0;
	mat xy;
	vec t = timesteps(0, 100, 300);
	vec circ = { -1, 1, 3 };
	Post_Process test(0);
	for (int i = 0; i <= 100; ++i) {
		xy = { { -20.0,-4.1+inc },{ -20.0,-3.9 + inc },{ 0,0} };
		Dormand_Price(t, circ, xy);
		test.Count_Files();
		test.Calc_Impact_Gen();
		test.Calc_Theta();
		test.Print_Result();
		test.File_Delete_Routine();
		inc += 0.1;
	}
	return 0;
*/
/*
double inc = 0;
	mat xy;
	vec t = timesteps(0, 100, 300);
	vec circ = { -1, 1, 1,1,1 };
	Post_Process test(0);
	for (int i = 0; i <= 100; ++i) {
		xy = { { -20.0,-4.1+inc },{ -20.0,-3.9 + inc },{ 0,0.15}, {0.0866, -0.05}, {-0.0866, -0.05} };
		Dormand_Price(t, circ, xy);
		test.Count_Files();
		test.Calc_Impact_Gen();
		test.Calc_Theta();
		test.Print_Result();
		test.File_Delete_Routine();
		inc += 0.1;
	}
	return 0;
}
*/




/*double inc = 0;
	mat xy;
	vec t = timesteps(0, 15, 100);
	vec circ = { -1,1,1,1 };
	Post_Process test(0);
	for (int i = 0; i <= 200; ++i) {
	xy = { { -10,-2 + inc },{ -10,-1.8 + inc },{ 0,-0.1 }, {0,0.1} };
	Dormand_Price(t, circ, xy);
	test.Count_Files();
	test.Calc_Impact_4R();
	test.Calc_Theta();
	test.Print_Result();
	test.File_Delete_Routine();
	inc += 0.025;
*/

/*mat xy(N,2,fill::zeros);
  vec circ={1,1,1,-1};
  Post_Process post(0);
  ofstream cluster_output;
  ofstream dipole_output;
  double temp_length;
  cluster_output.open("Cluster_Length.txt");
  dipole_output.open("Dipole_Length.txt");
  xy={ {2,0.1},{2,-0.1},{-2.43982,0.1},{-2.43982,-0.1} };
  Bound_Dormand_Price(circ,xy);
  post.Count_Files();
  for (int i=0; i<post.get_File_num(); i++){
    if (i%write_limit==0){
    post.Calc_Lengths_4N(i);
    cluster_output<<i/write_limit<<"\t"<<((post.get_L_12()<post.get_L_13()&&post.get_L_12()<post.get_L_23())?post.get_L_12():((post.get_L_13()<post.get_L_23())?post.get_L_13():post.get_L_23()))<<endl;
    dipole_output<<i/write_limit<<"\t"<<((post.get_L_34()<post.get_L_14()&&post.get_L_34()<post.get_L_24())?post.get_L_34():((post.get_L_14()<post.get_L_24())?post.get_L_14():post.get_L_24()))<<endl;
    }
  }
  return 0;
*/


