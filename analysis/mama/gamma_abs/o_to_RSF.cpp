/********************************************************************/
/*	Program to calculate gamma strength function from cross section */
/*	measurement data (photoneutron GDR experiments of type (n,g)	*/
/*																	*/
/*	Version Wednesday 25 April 2007									*/
/********************************************************************/

using namespace std;

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

int main()
{
	// open files
	
	ifstream gdr("varmalow.dat");	// input file
	ofstream rsf("rsf_varmalow.dat");	// output file
	
	int l = 49;					// length of input & output file
	
	float energy[l], cross[l], crosserr[l];	// vectors for input file. Units: energy: MeV, cross section: mb
	float RSF[l], RSFerr[l];	// vectors for output file. Units: energy: MeV, RSF: MeV^(-3)
	
	const float factor = 8.674E-08;	// const. factor in mb^(-1) MeV^(-2)
	
	// read from input to vectors
	
	int i = 0;
	while(gdr){
		gdr >> energy[i] >> cross[i] >> crosserr[i];
		i++;
	}
	
	// calculate RSF from cross section and write to outfile
	
	for(i = 0; i < l; i++){
		RSF[i] = factor*cross[i]/energy[i];
		RSFerr[i] = factor*crosserr[i]/energy[i];
		rsf << energy[i] << '\t'<< RSF[i] << '\t' << RSFerr[i] << endl;
	}
	
	// close files
	
	gdr.close();
	rsf.close();
}
