
#include "Point_Vortex.hpp"
#include "Post_Process.hpp"


//-------------------------------------------------------------------- Misc. Functions --------------------------------------------------------------------------------------------------------

double square(double x) {
	return x*x;
}

vec timesteps(double a, double b, int m) {	//Generate timestep mesh vector + Set h
	vec t_vals(m + 1, fill::zeros);
	h = (b - a) / m;
	t_vals[0] = a;
	for (int i = 1; i < m + 1; ++i) {
		t_vals[i] = t_vals[i - 1] + h;
	}
	return t_vals;
}


double per_s(double x, double y) {		//infinite sums approximation for periodic case 
	double old_val = sin(x) / (cosh(y) - cos(x));
	for (int m = 1; m < 10; ++m) {
		double new_val = old_val + (sin(x) / (cosh(y - (2 * pi*m)) - cos(x))) + (sin(x) / (cosh(y + (2 * pi*m)) - cos(x)));
		if (abs(new_val - old_val) < pow(10, -12)) return 0.5*old_val;
		old_val = new_val;
	}
	return 0.5*old_val;
}

double per_hamil_sum(double x, double y) {		//infinite sum approximation for the periodic hamiltonian
	double old_val = log(cosh(x) - cos(y));
	for (int m = 1; m < 1000; ++m) {
		double new_val = old_val + log((cosh(x - (2 * pi*m)) - cos(y)) / cosh(2 * pi*m)) + log((cosh(x + (2 * pi*m)) - cos(y)) / cosh(-2 * pi*m));
		if (abs(new_val - old_val) < pow(10, -12)) return old_val;
		old_val = new_val;
	}
	return old_val;
}

mat position_generate(double L, double phi) {       //generates 4 position matrix with 2nd dipole at angle phi at length L
	mat out = { {5,4.9},{5,5.1},{0,0},{0,0} };
	rowvec Vor_3 = { (0.1*cos(phi - (pi / 2))) + (L*cos(phi)),(0.1*sin(phi - (pi / 2))) + (L*sin(phi))+5 };
	rowvec Vor_4 = { (0.1*cos(phi + (pi / 2))) + (L*cos(phi)),(0.1*sin(phi + (pi / 2))) + (L*sin(phi))+5 };
	out.row(2) = Vor_3;
	out.row(3) = Vor_4;
	return out;
}

//----------------------------------------------------------------- Lengths + Velocities --------------------------------------------------------------------------------------------------------
//Functions that we will use to solve the system

double xij(int i, int j, mat xy) {  //calculate x distance between two vortices
	return xy(i, 0) - xy(j, 0);
}

double yij(int i, int j, mat xy) { //calculate y distance between two vortices
	return xy(i, 1) - xy(j, 1);
}

double length(int i, int j, mat xy) {         //calculate Euclidean legnth between two vortices
	return sqrt(square(xij(i, j, xy)) + square(yij(i, j, xy)));
}

double x_vel(int i, vec circ, mat xy) {		//calculate the x velocity of the ith vortex
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) { sum += (circ(j)*yij(i, j, xy)) / square(length(i, j, xy)); }
	}
	return -sum / (2 * pi);
}

double per_x_vel(int i, vec circ, mat xy) {	//calculate the x velocity of ith vortex (periodic case)
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) { sum += circ(j)*per_s(yij(i, j, xy), xij(i, j, xy)); }
	}
	return -sum / (2 * pi);
}

double y_vel(int i, vec circ, mat xy) {   //calculate the y velocity of the ith vortex
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) { sum += (circ(j)*xij(i, j, xy)) / square(length(i, j, xy)); }
	}
	return sum / (2 * pi);
}

double per_y_vel(int i, vec circ, mat xy) {	//calculate the y velocity of ith vortex (periodic case)
	double sum = 0;
	for (int j = 0; j < N; ++j) {
		if (i != j) { sum += circ(j)*per_s(xij(i, j, xy), yij(i, j, xy)); }
	}
	return sum / (2 * pi);
}

mat vel_matrix(mat xy, vec circ) {		//Creates matrix of vortex velocities at given positions
	mat vels(N, 2, fill::zeros);
	for (int i = 0; i < N; ++i) {
		vels(i, 0) = x_vel(i, circ, xy);
		vels(i, 1) = y_vel(i, circ, xy);
	}
	return vels;
}

mat per_vel_matrix(mat xy, vec circ) {   //Creates matrix of vortex velocities (periodic)
	mat vels(N, 2, fill::zeros);
	for (int i = 0; i < N; ++i) {
		vels(i, 0) = per_x_vel(i, circ, xy);
		vels(i, 1) = per_y_vel(i, circ, xy);
	}
	return vels;
}

//---------------------------------------------------------------------------- Conserved Quantities------------------------------------------------------------------------------------------------------------

double p(vec circ, mat xy) {			//linear momentum x direction
	double sum = 0;
	for (int i = 0; i < N; ++i) sum += circ(i)*xy(i, 0);
	return sum;
}

double q(vec circ, mat xy) {					//linear momentum y direction
	double sum = 0;
	for (int i = 0; i < N; ++i) sum += circ(i)*xy(i, 1);
	return sum;
}

double m(vec circ, mat xy) {				//angular momentum
	double sum = 0;
	for (int i = 0; i < N; ++i) sum += circ(i)*(square(xy(i, 0)) + square(xy(i, 1)));
	return sum;
}


double hamil(vec circ, mat xy) {		//Calculate the hamiltonian(energy) of the system, can be negative, periodic and infinite domain accounted for

	if (Flag_is_Periodic) {
		double sum = 0;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (i < j) sum += circ(i)*circ(j)*(per_hamil_sum(xij(i, j, xy), yij(i, j, xy)) - (xij(i, j, xy) / (2 * pi)));
			}
		}
		return -sum / (2 * pi);
	}
	else {
		double sum = 0;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (i < j) sum += circ(i)*circ(j)*log(length(i, j, xy));
				else;
			}
		}
		return -sum / pi;
	}
}

//------------------------------------------------------------------Dormand Price-------------------------------------------------------------------------------------------------------------
//Solves the sytstem and outputs as files  

double error_func(mat RK1, mat RK2){
	double sum=0;	
	double sc_x=0;
	double sc_y;
	for (int i=0; i<N; i++){
		sc_x=tol+(Rtol*fmax(RK1(i,0), RK2(i,0)));
		sc_y=tol+(Rtol*fmax(RK1(i,1), RK2(i,1)));
		sum+=square((RK1(i,0)-RK2(i,0))/sc_x);
		sum+=square((RK1(i,1)-RK2(i,1))/sc_y);
	}
	return sqrt(sum/(2*N));
}





cube rk45(mat xy, vec circ) {
	cube out(N, 2, 2, fill::zeros);			//initialize cube in order to return both solution matrices from the function
	if (Flag_is_Periodic) {			//periodic case
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

void dormand_price(vec t_steps, vec circ, mat xy) {
	ofstream writer, hamil_print, momen_print;		//initialize the three file output objects
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
	int counter = 1;				//count the files that have already been made
	for (double t = t_steps[1]; t <= end; t += h) {
		if (Flag_Print_Hamil) hamil_print << hamil(circ, xy) << endl;
		if (Flag_Print_Momen) momen_print << p(circ, xy) << "\t" << q(circ, xy) << "\t" << m(circ, xy) << endl;
		writer.open("Vortex_xy_" + to_string(counter) + ".txt");
		RK_Results = rk45(xy, circ);
		error=error_func(RK_Results.slice(0), RK_Results.slice(1));
		while (error>1){
			h*=min(hfacmax,max(hfacmin,hfac*pow(1/error,0.2)));
			RK_Results = rk45(xy, circ);
			error=error_func(RK_Results.slice(0), RK_Results.slice(1));
		}
		xy = RK_Results.slice(1); //set xy, print xy then repeat
		xy.raw_print(writer);
		writer.close();	
		++counter;
	}
	writer.open("python.txt");		//print circulation and file number for plotting
	writer << (int)counter << endl << circ << "\n\n";
	writer.close();
	if (Flag_Print_Momen) momen_print.close();
	if (Flag_Print_Hamil) hamil_print.close();	
	cout << "File output complete\n";
	h = h_old;
}

void bound_dormand_price(vec circ, mat xy) {      //run dormand price until final dipole is 20d + L away
  ofstream writer, hamil_print, momen_print;			//initialize the three file output objects
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
  int counter = 1;
  double init_d=xy(1,1)-xy(0,1);
  double init_L=hypot(xy(0,0)-xy(2,0), xy(0,1)-xy(2,1));
  vec lengths={hypot(xy(0,0)-xy(2,0),xy(0,1)-xy(2,1)),hypot(xy(0,0)-xy(1,0),xy(0,1)-xy(1,1)), hypot(xy(0,0)-xy(3,0), xy(0,1)-xy(3,1))}; 
  while(lengths.max()<(20*init_d)+init_L) {
    if (Flag_Print_Hamil) hamil_print <<counter/write_limit<<"\t"<<hamil(circ, xy) << endl;
    if (Flag_Print_Momen) momen_print << p(circ, xy) << "\t" << q(circ, xy) << "\t" << m(circ, xy) << endl;
    RK_Results = rk45(xy, circ);
    error = error_func(RK_Results.slice(0), RK_Results.slice(1));
    while (error>1){
			h*=min(hfacmax,max(hfacmin,hfac*pow(1/error,0.2)));
			RK_Results = rk45(xy, circ);
			error=error_func(RK_Results.slice(0), RK_Results.slice(1));
		}
    xy = RK_Results.slice(1); //set xy, print xy then repeat
    if (counter%write_limit==0&& counter<File_Max) {
      writer.open("Vortex_xy_" + to_string(counter/write_limit) + ".txt");
      xy.raw_print(writer);
      writer.close();
    } 
    else if (counter=File_Max){
    	break;
    }
    lengths={hypot(xy(0,0)-xy(2,0),xy(0,1)-xy(2,1)),hypot(xy(0,0)-xy(1,0),xy(0,1)-xy(1,1))};
    counter++;
  }
  writer.open("python.txt");		//print circulation and file number for plotting
  writer <<(int)counter/write_limit<<endl <<circ << "\n\n";
  writer.close();
  if (Flag_Print_Momen) momen_print.close();
  if (Flag_Print_Hamil) hamil_print.close(); 
  cout << "File output complete\n";
  h = h_old;
}


//------------------------------------------------------------------------------------------------------------------



int main() {
	mat xy;
	vec circ={-1,1,1,1,1};
	double dipole_end_d;
	Post_Process post(0);
	ofstream output;
	output.open("min-end-d.txt",ofstream::app);
	cout<< "matric+circ initialized";
	for (int i=0; i<=80; i++){
		for (int j=0;j<=160;j++){ 	
			xy={{-4.0-(j*pi*0.00125),-0.26+(i*0.00275)}, {-4.0-(j*pi*0.00125),-0.06+(i*0.00275)}, {0.0,0.1},{0.1, 0.0}, {0.0, -0.1}, {-0.1,0}}; 
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
}

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
/*
  mat xy(N,2,fill::zeros);
  vec circ={1,1,1,-1};
  double clustering_end_d;
  double dipole_end_d;
  Post_Process post(0);
  ofstream output;
  output.open("Min-end-d.txt",ofstream::app);
  for (int i=0; i<=20; i++){
    for (int j=0;j<=40;j++){
      xy={ {2,0.1},{2,-0.1},{-2-(j*pi*0.005),-0.5+(i*0.06)},{-2-(j*pi*0.005),-0.7+(i*0.06)} };
      Bound_Dormand_Price(circ,xy);
      post.Count_Files();
      post.Calc_Impact_4R();
      post.Calc_Lengths_4N(post.get_File_num()-1);
      dipole_end_d=(post.get_L_34()<post.get_L_14()&&post.get_L_34()<post.get_L_24())?post.get_L_34():((post.get_L_14()<post.get_L_24())?post.get_L_14():post.get_L_24());
      clustering_end_d=(post.get_L_12()<post.get_L_13()&&post.get_L_12()<post.get_L_23())?post.get_L_12():((post.get_L_13()<post.get_L_23())?post.get_L_13():post.get_L_23());
      output<<post.get_Di_Length()<<"\t"<<post.get_Impact()/0.2<<"\t"<<post.Min_Lengths_4N()(5)<<"\t"<<dipole_end_d<<"\t"<<clustering_end_d<<endl;
      post.File_Delete_Routine();
    }
  }
  output.close();
  return 0;*/

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
  return 0;*/


