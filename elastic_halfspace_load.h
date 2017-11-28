#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14159265358979323846

extern double r_ij(int i, int j, double x, double y, double z, double a, double b);
extern double beta_ij(int i, int j, double x, double y, double z, double a, double b);
extern double phi_ij(int i, int j, double x, double y, double z, double a, double b);
extern double J_ij(int j, double x, double y, double z, double a, double b);
extern double K_ij(int j, double x, double y, double z, double a, double b);
extern double L_ij(int j, double x, double y, double z, double a, double b);
extern double u(double x, double y, double z, double a, double b, double lame, double mu, double pressure);
extern double v(double x, double y, double z, double a, double b, double lame, double mu, double pressure);
extern double w(double x, double y, double z, double a, double b, double lame, double mu, double pressure);

//
//  These functions are meant to compute the elastic displacement from a rectangular load 
//   with half-width a and half-height b (meters) and pressure p centered at the origin. 
//   They will return the 3-component displacements at a point (x, y, z) in meters. 
//   


double r_ij(int i, int j, double x, double y, double z, double a, double b){
/*
	r10 and r20: This is the radius between the point xyz and 
			the points x=a,y=0,z=0 or x=-a,y=0,z=0 
	r01 and r02: This is the radius between the point xyz and 
			the points x=0,y=b,z=0 or x=0,y=-b,z=0 
*/
	double r = 0;
	if (i == 1 && j == 0) {
		r=sqrt(pow(a-x,2) + pow(y,2) + pow(z,2) );
	}
	else if (i == 2 && j == 0) {
		r=sqrt(pow(a+x,2) + pow(y,2) + pow(z,2) );
	}
	else if (i == 0 && j == 1) {
		r=sqrt(pow(x,2) + pow(b-y,2) + pow(z,2) );
	}
	else if (i == 0 && j == 2) {  // i == 0 and j == 2
		r=sqrt(pow(x,2) + pow(b+y,2) + pow(z,2) );
	}
	else {
		printf("ERROR in R_IJ: i = %i and j = %i", i, j);
	}
	return r;
}


double beta_ij(int i, int j, double x, double y, double z, double a, double b){
/*
	This is a projection of the radius into one of the central planes.
	beta10 and beta20 are (a-/+x)^2 + z^2
	beta01 and beta02 are (b-/+y)^2 + z^2
*/
	double beta;
	if (j==0){  // the 10 or 20 cases
		beta = sqrt( pow(r_ij(i,j,x,y,z,a,b),2) - pow(y,2) );
	}
	else{  // the 01 or 02 cases
		beta = sqrt( pow(r_ij(i,j,x,y,z,a,b),2) - pow(x,2) );
	}
	return beta;
}


double phi_ij(int i, int j, double x, double y, double z, double a, double b){
/*
	phi10 and phi20 are fractional radius with y on top.
	phi01 and phi02 are fractional radius with x on top. 
*/
	double tol = 1.0*pow(10,-9);
	double phi; 
	if ( (r_ij(i,j,x,y,z,a,b)+beta_ij(i,j,x,y,z,a,b)) < tol ) {
		phi=1;   // prevent dividing by zero; replace with the largest possible number. 
	}
	else{
		if (j==0) {  // the 10 or 20 cases
			phi=y/(r_ij(i,j,x,y,z,a,b)+beta_ij(i,j,x,y,z,a,b));
		}
		if (j>0) {  // the 01 or 02 cases
			phi=x/(r_ij(i,j,x,y,z,a,b)+beta_ij(i,j,x,y,z,a,b));
		}
	}
	return phi;
}


double J_j(int j, double x, double y, double z, double a, double b){
	double J;
	double tol= 1.0*pow(10, -9);
	double prefactor, term1, term2, term3; 
	if (j==1) {
		prefactor=fabs(a-x);
	}
	if (j==2) {
		prefactor=fabs(a+x);
	}
	if ( fabs(z)<tol ) {  // limit issues with solving for surface displacements
		if ( fabs(z+r_ij(j,0,x,y,z,a,b)) < tol){
			term1=0.0;
		}
		else {
			term1= y * ( log(z+r_ij(j,0,x,y,z,a,b)) - 1.0);
		}
		term2=0;
		term3=2*prefactor*atan(phi_ij(j,0,x,y,z,a,b));
	}
	else {
		term1=y*(log(z+r_ij(j,0,x,y,z,a,b))-1.0);
		term2=z*log((1.0+phi_ij(j,0,x,y,z,a,b))/(1.0-phi_ij(j,0,x,y,z,a,b)));
		term3=2*prefactor*atan(prefactor*phi_ij(j,0,x,y,z,a,b)/(z+beta_ij(j,0,x,y,z,a,b)));
	}
	J=term1+term2+term3;
	return J;
}


double K_j(int j, double x, double y, double z, double a, double b){
	double K;
	double tol= 1.0*pow(10, -9);
	double prefactor, term1, term2, term3; 
	if (j==1) {
		prefactor=fabs(b-y);
	}
	if (j==2) {
		prefactor=fabs(b+y);
	}
	if ( fabs(z)<tol ) {  // limit issues with solving for surface displacements. 
		if ( fabs(z+r_ij(0,j,x,y,z,a,b))<tol ){
			term1=0;
		}
		else {
			term1 = x*(log(z+r_ij(0,j,x,y,z,a,b))-1.0);
		}
		term2=0;
		term3 = 2*prefactor*atan(phi_ij(0,j,x,y,z,a,b));
	}
	else {
		term1 = x*(log(z+r_ij(0,j,x,y,z,a,b))-1.0);
		term2 = z*log((1+phi_ij(0,j,x,y,z,a,b))/(1-phi_ij(0,j,x,y,z,a,b)));
		term3 = 2*prefactor*atan(prefactor*phi_ij(0,j,x,y,z,a,b)/(z+beta_ij(0,j,x,y,z,a,b)));
	}
	K = term1+term2+term3;
	return K;
}

double L_j(int j, double x, double y, double z, double a, double b){
	double L;
	double tol= 1.0*pow(10, -9);
	double prefactor, term1, term2, term3; 

	if ( j==1 ){
		prefactor=a-x;
	}
	if ( j==2 ){
		prefactor=-a-x;
	}
	if ( fabs(z)<tol ) {
		if ( fabs(prefactor+r_ij(j,0,x,y,z,a,b))<tol ) {
			term1=0;
		}
		else {
			term1=y*(log(prefactor+r_ij(j,0,x,y,z,a,b))-1);
		}
		if ( fabs(prefactor)<tol ) {
			term2=0;
		}
		else {
			term2=prefactor*log((1+phi_ij(j,0,x,y,z,a,b))/(1-phi_ij(j,0,x,y,z,a,b)));
		}
		term3=0;
	}
	else {
		term1=y*(log(prefactor+r_ij(j,0,x,y,z,a,b))-1);
		term2=prefactor*log((1+phi_ij(j,0,x,y,z,a,b))/(1-phi_ij(j,0,x,y,z,a,b)));
		term3=2*z*atan((z*phi_ij(j,0,x,y,z,a,b))/(prefactor+beta_ij(j,0,x,y,z,a,b)));
	}
	L = term1+term2+term3;
	return L;
}



double u(double x, double y, double z, double lame, double mu, double p, double a, double b){
/*
	Displacements in the x direction. For the corners of the box, we set the exploding terms manually to 0.
*/
	double tol=1e-9;
	double u_top, u_bottom;
	double deltay_top = -b-y;
	double deltay_bottom = b-y;
	if ( (fabs(deltay_top)<tol) || (fabs(deltay_bottom)<tol) || (fabs(-a-x)<tol) || (fabs(a-x)<tol) ){
		if ( fabs(z)<tol ) {
			u_top=(-p/4/pi)*((J_j(2,x,deltay_top,z,a,b)-J_j(1,x,deltay_top,z,a,b))*(1/(lame+mu)));
			u_bottom=(-p/4/pi)*((J_j(2,x,deltay_bottom,z,a,b)-J_j(1,x,deltay_bottom,z,a,b))*(1/(lame+mu)));
		}
	}
	else {
		u_top=(-p/4/pi)*((J_j(2,x,deltay_top,z,a,b)-J_j(1,x,deltay_top,z,a,b))*(1/(lame+mu))+(z/mu)*(log((deltay_top+r_ij(2,0,x,deltay_top,z,a,b))/(deltay_top+r_ij(1,0,x,deltay_top,z,a,b)))));
		u_bottom=(-p/4/pi)*((J_j(2,x,deltay_bottom,z,a,b)-J_j(1,x,deltay_bottom,z,a,b))*(1/(lame+mu))+(z/mu)*(log((deltay_bottom+r_ij(2,0,x,deltay_bottom,z,a,b))/(deltay_bottom+r_ij(1,0,x,deltay_bottom,z,a,b)))));
	}
	return -1*(u_top - u_bottom);
}


double v(double x, double y, double z, double lame, double mu, double p, double a, double b){
/*
	Displacements in the y direction. For the corners of the box, we set the exploding terms manually to 0.
*/
	double tol=1e-9;
	double v_top, v_bottom; 
	double deltax_top = -a-x;
	double deltax_bottom = a-x;
	if ( (fabs(deltax_top)<tol) || (fabs(deltax_bottom)<tol) || (fabs(-b-y)<tol) || (fabs(b-y)<tol) ) {
		v_top=(-p/4/pi)*((K_j(2,deltax_top,y,z,a,b)-K_j(1,deltax_top,y,z,a,b))*(1/(lame+mu)));
		v_bottom=(-p/4/pi)*((K_j(2,deltax_bottom,y,z,a,b)-K_j(1,deltax_bottom,y,z,a,b))*(1/(lame+mu)));
	}
	else {
		v_top=(-p/4/pi)*((K_j(2,deltax_top,y,z,a,b)-K_j(1,deltax_top,y,z,a,b))*(1/(lame+mu))+(z/mu)*(log((deltax_top+r_ij(0,2,deltax_top,y,z,a,b))/(deltax_top+r_ij(0,1,deltax_top,y,z,a,b)))));
		v_bottom=(-p/4/pi)*((K_j(2,deltax_bottom,y,z,a,b)-K_j(1,deltax_bottom,y,z,a,b))*(1/(lame+mu))+(z/mu)*(log((deltax_bottom+r_ij(0,2,deltax_bottom,y,z,a,b))/(deltax_bottom+r_ij(0,1,deltax_bottom,y,z,a,b)))));
	}
	return -1*(v_top - v_bottom);
}


double w(double x, double y, double z, double lame, double mu, double p, double a, double b){
/*
	Displacements in the z direction. For the corners of the box, we set the exploding terms manually to 0.
*/
	double tol=1e-9;
	double w_top, w_bottom; 
	double deltay_top = -b-y;
	double deltay_bottom = b-y;
	if ( fabs(z)<tol ) {
		w_top=(p/4/pi/mu)*((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_top,z,a,b)-L_j(2,x,deltay_top,z,a,b));
		w_bottom=(p/4/pi/mu)*((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_bottom,z,a,b)-L_j(2,x,deltay_bottom,z,a,b));
	}
	else {	
		w_top=(p/4/pi/mu)*(((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_top,z,a,b)-L_j(2,x,deltay_top,z,a,b)) + z*(atan(((a-x)*deltay_top)/(z*r_ij(1,0,x,deltay_top,z,a,b)))+atan(((a+x)*deltay_top)/(z*r_ij(2,0,x,deltay_top,z,a,b)))));
		w_bottom=(p/4/pi/mu)*(((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_bottom,z,a,b)-L_j(2,x,deltay_bottom,z,a,b)) + z*(atan(((a-x)*deltay_bottom)/(z*r_ij(1,0,x,deltay_bottom,z,a,b)))+atan(((a+x)*deltay_bottom)/(z*r_ij(2,0,x,deltay_bottom,z,a,b)))));
	}
	return (w_top - w_bottom);
}




