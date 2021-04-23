/*mat xy;    4n scatter angles
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
		xy = Position_Generate(L, phi);
		Dormand_Price(t, circ, xy);
		post.Count_Files();
		post.Calc_Theta();
		post.Print_Result_4N();
		post.File_Delete_Routine();
	}
}
*/

/* Plot 4N Min Lengths
   	mat xy;
	vec t = timesteps(0, 15, 100);
	vec circ = { 1,-1,1,-1 };
	Post_Process post(0);
	double L_inc = 0.078125;
	double phi = -pi/96;;
	double phi_inc = pi / 96;
	double L;
	for (int i = 1; i < 97; ++i) {
		phi += phi_inc;
		L = 1.875;
		for (int j = 1; j <= 80; ++j) {
			L += L_inc;
			xy = Position_Generate(L, phi);
			Dormand_Price(t, circ, xy);
			post.update_4N(phi, L);
			post.Count_Files();
			post.Print_Min_Lengths_4N();
			cout << "Min Lengths updated"<<endl;
			post.File_Delete_Routine();
			}
	}
*/


/* Plot scatter angles
double inc = 0;
mat xy;
vec t = timesteps(0, 15, 100);
vec circ = { -1,1,1 };
Post_Process test(0);
for (int i = 0; i <= 30; ++i) {
xy = { { -2,-0.5 + inc },{ -2,-0.3 + inc },{ 2,0 } };
Dormand_Price(t, circ, xy);
test.Count_Files();
test.Calc_Impact();
if (test.get_Impact() != 0.0 && test.get_Impact() != 0.9) {
test.Calc_Theta();
test.Print_Result();
}
test.File_Delete_Routine();
inc += 0.05;
}
return 0;



*/

/* Plot min lengths
double inc = 0;
mat xy;
vec t = timesteps(0, 15, 100);
vec circ = { -1,1,1 };
Post_Process test(0);
for (int i = 0; i <= 30; ++i) {
xy = { { -2,-0.5 + inc },{ -2,-0.3 + inc },{ 2,0 } };
Dormand_Price(t, circ, xy);
test.Count_Files();
test.Calc_Impact();
if (test.get_Impact() != 0.0 && test.get_Impact() != 0.9) {
test.Print_Min_Lengths();
}
test.File_Delete_Routine();
inc += 0.05;
}
return 0;*/

/*Plotting for min length theory
double inc = 0;
mat xy;
vec t = timesteps(0, 15, 100);
vec circ = { -1,1,1 };
ofstream writer;
int End_Flag = 0;
writer.open("min_length_theory.txt", ofstream::app);

while (End_Flag!=1){
xy = { { -2,-0.5 + inc },{ -2,-0.3 + inc },{ 2,0 } };
double rho = (xy(0, 1) + xy(1, 1)) / 2;
double d = length(0, 1, xy);
vec min_lengths = { 0,0,0 };
if (rho / d < 0) {
min_lengths(0) = ((2 * rho) - d + sqrt((4 * square(rho)) - (20 * rho*d) + (9 * square(d)))) / 4;
min_lengths(2) = (d - (2 * rho)) / 2;
min_lengths(1) = min_lengths(0) + min_lengths(2);
}
if (rho / d > 4.5) {
min_lengths(1) = ((2 * rho) - d + sqrt((4 * square(rho)) - (20 * rho*d) + (9 * square(d)))) / 4;
min_lengths(2) = ((2 * rho) - d) / 2;
min_lengths(0) = d;
}
if (rho / d > 0 && rho / d < 4.5) {
min_lengths(0) = d;
min_lengths(1) = d;
min_lengths(2) = d + sqrt(2 * rho*d);
}
writer << rho/d << "\t" << min_lengths(0)/d<<"\t"<< min_lengths(1)/d << "\t" << min_lengths(2)/d << endl;
if (rho / d > 6) End_Flag = 1;
inc += 0.001;

}
writer.close();
return 0
*/


/*plotting for scattering lengths theory
double inc = 0;
mat xy;
vec t = timesteps(0, 15, 100);
vec circ = { -1,1,1 };
ofstream writer;
int End_Flag = 0;
writer.open("scatter_length_theory.txt", ofstream::app);

while (End_Flag!=1){
xy = { { -2,-0.5 + inc },{ -2,-0.3 + inc },{ 2,0 } };
double rho = (xy(0, 1) + xy(1, 1)) / 2;
double d = length(0, 1, xy);
vec min_lengths = { 0,0,0 };
if (rho / d < 0) {
min_lengths(0) = ((2 * rho) - d + sqrt((4 * square(rho)) - (20 * rho*d) + (9 * square(d)))) / 4;
min_lengths(2) = (d - (2 * rho)) / 2;
min_lengths(1) = min_lengths(0) + min_lengths(2);
}
if (rho / d > 4.5) {
min_lengths(1) = ((2 * rho) - d + sqrt((4 * square(rho)) - (20 * rho*d) + (9 * square(d)))) / 4;
min_lengths(2) = ((2 * rho) - d) / 2;
min_lengths(0) = min_lengths(2) - min_lengths(1);
}
if (rho / d > 0 && rho / d < 4.5) {
min_lengths(0) = sqrt(square(d) + (d*sqrt(2 * rho*d)));
min_lengths(1) = sqrt(square(d) + (d*sqrt(2 * rho*d)));
min_lengths(2) = d + sqrt(2 * rho*d);
}
writer << rho/d << "\t" << min_lengths(0)/d<<"\t"<< min_lengths(1)/d << "\t" << min_lengths(2)/d << endl;
if (rho / d > 6) End_Flag = 1;
inc += 0.001;

}
writer.close();
return 0;*/
