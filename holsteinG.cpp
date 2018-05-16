#include <iostream>
#include<fstream>
#include <complex>
#include <stdio.h>
#include <ctime>
#include <sys/time.h>
#include <stdlib.h> // to allocate additional memory
#include <omp.h>
#include <algorithm>

using namespace std;
typedef complex<double> dcomp;

const dcomp I(0.0, 1.0);	// imaginary unit i
const double pi = 3.1415926535897;	// pi

// Model Parameters
const double epsimp = 1.0;
const double lambda = 0.1;
const double Omega = 6.0;
const double t = 1;
const double g = sqrt(2*t*abs(Omega)*lambda);
const double eta = 1e-3; // width of delta peaks

// Approximations
const double eps = 1e-3; // error threshold

// Discretization (k,w)
const int Nk = 1e4;
const int Nw = 2e4;
double wstart = -10.0;
double dw = 20.0/Nw;
double dk = pi/Nk;

double Ek(int k);

dcomp g0(double w);
dcomp calcA(double w);

dcomp calcG0(double w);

int main(){
	int iw;
	double w;
	dcomp A, G, G0;

	
	FILE * param;
	param = fopen("Parameters.txt","w");
	fprintf(param, "t\t%e\n",t);
	fprintf(param, "g\t%e\n",g);
	fprintf(param, "epsimp\t%e\n",g);	
	fprintf(param, "Omega\t%e\n",Omega);
	fprintf(param, "Eta\t%e\n",eta);
	fprintf(param, "eps\t%e\n",eps);
	fprintf(param, "dk\t%e\n",dk);
	fprintf(param, "dw\t%e\n",dw);
	
		
	FILE * greensOut; // Ouput Data
	greensOut = fopen("greensOut.txt","w");
	
	for (iw = 0; iw<Nw; iw++){
		w = double(iw)*dw+wstart;
		A = calcA(w);
		G0 = calcG0(w);
			G = G0/(1.0-G0*(epsimp+g*A));
			fprintf(greensOut, "%e\t%e\n", real(G),imag(G));
	}
	fclose(greensOut);
	

	return 0;
}

dcomp calcA(double w){ // calculate A1
	dcomp b0, bn, a1, an;
	dcomp x; // A
	
	int j; // phonon # 
	dcomp C, D, Delta;
	
	dcomp Gj, Gj1;
		
	b0 = 0.0;

	if (b0==0.0){
		x = 1e-30;
	}
	else x = b0;
	
	C = x;
	D = 0.0;
	j = 1;
	Delta = 0.0;
	
	Gj = calcG0(w-Omega);
	
	bn = 1.0-Gj*epsimp;
	a1 = g*Gj;
	D = bn;
	if (D==0.0){
		D = 1e-30;
	}
	C = bn+a1/C;
	if (C==0.0){
		C = 1e-30;
	}
	D = 1.0/D;
	Delta = C*D;
	x = x*Delta;
	
	while (abs(Delta-1.0)>eps){
		j++;
		
		Gj = calcG0(w-double(j)*Omega);
		Gj1 = calcG0(w-double(j-1)*Omega);
		
		bn = 1.0-Gj*epsimp;
		an = -pow(g,2.0)*double(j)*Gj*Gj1;
		D = bn+an*D;
		if (D==0.0){
			D = 1e-30;
		}
		C = bn+an/C;
		if (C==0.0){
			C = 1e-30;
		}
		D = 1.0/D;
		Delta = C*D;
		x = x*Delta;
	}
	return x;	
}

dcomp g0(double w){ // Calculate the MA
	int ik;
	dcomp x = 0.0;
	for (ik=0; ik<Nk; ik++){
		x+= 1.0/(w-Ek(ik)+I*eta);
		//x+= 1.0/(w+Ek(ik)+I*eta);
	}
	return x/double(Nk);
}


dcomp calcG0(double w){
	
	dcomp b0, bn, a1, an, a2;
	dcomp x; // A
	
	int j; // phonon # 
	dcomp C, D, Delta;
		
	b0 = 0.0;

	if (b0==0.0){
		x = 1e-30;
	}
	else x = b0;
	
	C = x;
	D = 0.0;
	Delta = 0.0;
	
	j = 1;
	
	bn = w+I*eta;
	a1 = 1.0;
	D = bn+a1*D;
	if (D==0.0){
		D = 1e-30;
	}
	C = bn+a1/C;
	if (C==0.0){
		C = 1e-30;
	}
	D = 1.0/D;
	Delta = C*D;
	x = x*Delta;
	
	j = 2;
	
	bn = w+I*eta;
	a2 = -2.0*t*t;
	D = bn+a2*D;
	if (D==0.0){
		D = 1e-30;
	}
	C = bn+a1/C;
	if (C==0.0){
		C = 1e-30;
	}
	D = 1.0/D;
	Delta = C*D;
	x = x*Delta;
	
	while (abs(Delta-1.0)>eps){
		j++;
		bn = w+I*eta;
		an = -t*t;
		D = bn+an*D;
		if (D==0.0){
			D = 1e-30;
		}
		C = bn+an/C;
		if (C==0.0){
			C = 1e-30;
		}
		D = 1.0/D;
		Delta = C*D;
		x = x*Delta;
	}
	return x;	
	
}

double Ek(int k){ // calculate the band dispersion at momentum k
	return -2.0*t*cos(double(k)*dk);
}
