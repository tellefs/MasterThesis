/****************************************************************************/
/*	Program to calculate theoretical gamma strength functions of			*/
/*	the GEDR, GMDR and GEQR.												*/
/*	Written by Cecilie.														*/
/*	Version Monday 18 February 2008.										*/
/*																			*/
/*	Useful databases: http://www-nds.iaea.org/RIPL-2/						*/
/*					  http://cdfe.sinp.msu.ru/services/unifsys/index.html	*/
/*					  ROBIN 1.3												*/
/****************************************************************************/

using namespace std;

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

// General constants
const double Pi = 3.14159;
const double factor   = 8.6737E-08;	// const. factor in mb^(-1) MeV^(-2)

// Constants for the specific nucleus 
const double A = 233.0;
const double Z = 92.0;
const double N = A - Z;
const double a = 25.395;	// level density parameter from Egidy&Bucurescu (MeV^{-1})
const double E_1      = -0.519 ;	// backshift parameter (MeV)
const double beta2    = 0.22;	// G.s. deformation of 232Th. Check RIPL
const double temp     = 0.20; // constant temperature (MeV) with pygmy
const double E2       = 1.083; // Energy of first excited VIBRATIONAL 2+ state, to be used in function GenFerLiq (see RIPL-2)

// Limits in first-generation matrix for RhoSigChi, to be used in to be used in function KMF_var
//const double Ex_min   = 4.51;	// Min excitation energy in f.g. matrix, total and lower range
//const double Ex_max   = 7.0745;	// Max excitation in f.g. matrix, lower Ex range
//const double Ex_min   = 7.3079;	// Min excitation in f.g. matrix, higher Ex range
const double Ex_min   = 4;	// Max excitation in f.g. matrix, total and higher Ex range 
const double Ex_max   = 5;	// Max excitation in f.g. matrix, total and higher Ex range 

// GEDR parameters. Two sets if the nucleus has a significant ground-state deformation, else, comment away one set. 
// (GEDR: Giant Electric Dipole Resonance)
// If experimental GEDR parameters are known, use the last sets. Parameterization is taken from RIPL-2.

// Spherical nuclei, parameterization:
const double E_0      = 31.2*pow(A,-1./3.) + 20.6*pow(A,-1./6.); //  (MeV)
const double Gamma_0  = 0.026*pow(E_0,1.91); //  (MeV)
const double sigma_0  = 1.2*120.*N*Z/(A*Pi*Gamma_0); // (mb)

// If spherical nucleus and no exp. GEDR parameters, use the following three lines
//const double E_r1      = E_0; //  (MeV)
//const double Gamma_r1  = Gamma_0; //  (MeV)
//const double sigma_r1  = sigma_0; // (mb)

// Deformed nuclei, parameterization:
const double alpha2 = beta2*sqrt(5/(4*Pi));
const double lambda3 = 1. + 0.6*pow(alpha2,2.0) + 2.*pow(alpha2,3.0)/35.;
const double lambda  = pow(lambda3, 1./3.);
const double a0      = (1. + alpha2)/lambda;
const double b0      = (1. - 0.5*alpha2)/lambda;

// If deformed nucleus and no exp. GEDR parameters, use the following six lines	  
//const double E_r2 = E_0*(1. - 0.0151*(pow(a0,2.) - pow(b0,2.)))/b0;
//const double E_r1 = E_r2/((0.911*a0/b0) + 0.089);
//const double Gamma_r1 = 0.026*pow(E_r1,1.91);
//const double Gamma_r2 = 0.026*pow(E_r2,1.91);
//const double sigma_r1 = 0.667*sigma_0; 
//const double sigma_r2 = 0.333*sigma_0;

// Experimental GEDR parameters (see atlas of Dietrich & Berman, or RIPL-2)
const double E_r1      = 11.3; //  (MeV)
const double Gamma_r1  = 3.3*0.9; //  (MeV)
const double sigma_r1  = 440.0*2.2; // (mb)

const double E_r2      = 14.6; // (MeV)
const double Gamma_r2  = 4.4; // (MeV)
const double sigma_r2  = 800.*1.9; // (mb)

//const double E_r3      = 23.2; // (MeV)
//const double Gamma_r3  = 4.0; // (MeV)
//const double sigma_r3  = 30.0; // (mb)

// Declaring all functions
double StandLor(double , double , double , double );
double EnGeLor(double , double , double , double , double );
double EnGeLorTemp(double , double , double , double );
double ModLor(double , double , double , double );
double ModLorDef(double , double , double , double , double , double , double );
double GenFerLiq(double , double , double , double );
double HybridFormulaE1(double , double , double , double );
double KadMarFur(double , double , double , double );
double SpinFlipM1(double , double , double , double );
double IsoScalarE2(double );
double KMF(double , double , double , double , double );
double KMF_var(double , double , double , double );
double PygmyT1(double );
double PygmyT2(double );
double PygmyH1(double );
double PygmyH2(double );
double PygmyL1(double );
double PygmyL2(double );
double ModGenLor(double , double , double , double , double );


int main()
{	
	// Calculate TRK sum rule and print to screen
	double TRK =  60.*N*Z/A;
	cout << endl;
	cout << "=================================================================================" << endl;
	cout << "The Thomas-Reiche-Kuhn sum rule is " << TRK << " mb MeV. " << endl;		
	cout << "GDR parameters: E_r1 = " << E_r1 << " MeV, Gamma_r1 = " << Gamma_r1 << " MeV, sigma_r1 = " << sigma_r1 << " mb. " << endl;
	cout << "                E_r2 = " << E_r2 << " MeV, Gamma_r2 = " << Gamma_r2 << " MeV, sigma_r2 = " << sigma_r2 << " mb. " << endl;
	cout << "=================================================================================" << endl;
	
	cout << endl;
	
	// open file to write results on
	ofstream outfile("GDR_models_238Np.dat");	// output file
//	ofstream outfile2("GLOvarT_Tmax.dat");	// output file with experimental energy calibration for checking consistency with f.g. spectra
	
	const double step  = 0.1;		// Steps of calculation (MeV)
	const double E_max = 40;			// Maximum gamma energy (MeV)
	
	int n = (int) E_max/step;		// Number of calculated points
	double E_g = step;			// Start value of E_gamma
	
	// Declarations of vectors
	double energy[n], SLo[n], EGLo[n], MLo[n], GFL[n], HybFor[n]; // Energy (MeV), strength functions (MeV^(-3))
	double f_KMF[n], f_KMFvar[n], f_M1[n], f_E2[n], f_tot[n];	
	double E1comp1[n], E1comp2[n], E1comp3[n];
	double ModifiedGLO[n];
	double pygmyT[n], pygmyT2[n], pygmyT1[n], pygmyH[n], pygmyH2[n], pygmyH1[n], pygmyL[n], pygmyL2[n], pygmyL1[n];
	
	//Calculation of E1, M1 and E2 strength. Sum of two functions for the E1 part due to deformation splitting of GDR
	for(int i=0; i<n; i++){		
		SLo[i]      = StandLor(E_g,E_r1,Gamma_r1,sigma_r1) + StandLor(E_g,E_r2,Gamma_r2,sigma_r2);
		E1comp1[i]  = EnGeLor(E_g,E_r1,Gamma_r1,sigma_r1,temp);
		E1comp2[i]  = EnGeLor(E_g,E_r2,Gamma_r2,sigma_r2,temp);
//		E1comp3[i]  = EnGeLor(E_g,E_r3,Gamma_r3,sigma_r3);
		EGLo[i]     = EnGeLor(E_g,E_r1,Gamma_r1,sigma_r1,temp) + EnGeLor(E_g,E_r2,Gamma_r2,sigma_r2,temp); // + EnGeLor(E_g,E_r3,Gamma_r3,sigma_r3);
		MLo[i]      = ModLorDef(E_g,E_r1,Gamma_r1,sigma_r1,E_r2,Gamma_r2,sigma_r2);
		GFL[i]      = GenFerLiq(E_g,E_r1,Gamma_r1,sigma_r1) + GenFerLiq(E_g,E_r2,Gamma_r2,sigma_r2);
		HybFor[i]   = HybridFormulaE1(E_g,E_r1,Gamma_r1,sigma_r1) + HybridFormulaE1(E_g,E_r2,Gamma_r2,sigma_r2);
		f_KMF[i]    = KMF(E_g,E_r1,Gamma_r1,sigma_r1,temp) + KMF(E_g,E_r2,Gamma_r2,sigma_r2,temp);
		f_KMFvar[i] = KMF_var(E_g,E_r1,Gamma_r1,sigma_r1) + KMF_var(E_g,E_r2,Gamma_r2,sigma_r2);		
		f_M1[i]     = SpinFlipM1(E_g,E_r1,Gamma_r1,sigma_r1);
		f_E2[i]     = IsoScalarE2(E_g);
		f_tot[i]    = EGLo[i] + f_M1[i] + f_E2[i];
        
        pygmyL1[i]   = PygmyL1(E_g);
        pygmyL2[i]   = PygmyL2(E_g);
        pygmyL[i]    = PygmyL1(E_g) + PygmyL2(E_g);		
        
        pygmyH1[i]   = PygmyH1(E_g);
        pygmyH2[i]   = PygmyH2(E_g);
        pygmyH[i]    = PygmyH1(E_g) + PygmyH2(E_g);
        
		pygmyT1[i]   = PygmyT1(E_g);
        pygmyT2[i]   = PygmyT2(E_g);
        pygmyT[i]    = PygmyT1(E_g) + PygmyT2(E_g);

		ModifiedGLO[i] = ModGenLor(E_g,E_r1,Gamma_r1,sigma_r1,temp) + ModGenLor(E_g,E_r2,Gamma_r2,sigma_r2,temp);
		
		// Write some of the functions to file (before Slo[i] was at f_E2[i])
		outfile <<  E_g << '\t' << f_E2[i] << '\t' << EGLo[i] << '\t' << 
        EGLo[i]+f_M1[i]+f_E2[i] << '\t' <<  f_M1[i] << '\t' << 
        pygmyT1[i] << '\t' <<  pygmyT2[i] << '\t' << pygmyT[i] << '\t' << 
        pygmyH1[i] << '\t' <<  pygmyH2[i] << '\t' << pygmyH[i] << '\t' <<
        pygmyL1[i] << '\t' <<  pygmyL2[i] << '\t' << pygmyL[i] << '\t' <<
        HybFor[i]+f_M1[i] << endl;
        
		E_g += step;
	}
	
/*	double a0 = 0.000;
	double a1 = 0.100;
	double e_exp[150], RSF[150];
	
	for(int i=0; i<150; i++){
		e_exp[i] = a0 + a1*(double)i;
//		RSF[i]   = KMF_var(e_exp[i],E_r1,Gamma_r1,sigma_r1) + KMF_var(e_exp[i],E_r2,Gamma_r2,sigma_r2) + SpinFlipM1(e_exp[i],E_r1,Gamma_r1,sigma_r1);
		RSF[i]   = EnGeLorTemp(e_exp[i],E_r1,Gamma_r1,sigma_r1) + EnGeLorTemp(e_exp[i],E_r2,Gamma_r2,sigma_r2) + SpinFlipM1(e_exp[i],E_r1,Gamma_r1,sigma_r1);
		if(i>8 && i<150) outfile2 << i << '\t' << e_exp[i] << '\t' << RSF[i] << endl;
		else outfile2 << i << '\t' << e_exp[i] << '\t' << 0.0000 << endl;
	}
 */
	// close files
	outfile.close();
//	outfile2.close();
	
}	// END OF main()

// Functions to be used in main(), see RIPL-1, RIPL-2 and references therein

// Standard Lorentzian, Brink-Axel approach
double StandLor(double E_g, double E_r, double Gamma_r, double sigma_r){
	double SLo = 0.;
	double denominator = (pow(E_g,2.0) - pow(E_r,2.0))*(pow(E_g,2.0) - pow(E_r,2.0)) + pow(Gamma_r,2.0)*pow(E_g,2.0);
	
	SLo = factor*sigma_r*pow(Gamma_r,2.0)*E_g/denominator;
	return SLo;
}

// Enhanced Generalized Lorentzian
double EnGeLor(double E_g, double E_r, double Gamma_r, double sigma_r, double temp){
	double EGLo = 0.;
	double Kappa = 0.;
	double Kappa_0 = 0.;
	double k = 0.;
	double Gamma_k = 0.;
	double Gamma_k0 = 0.;
	double denominator = 0.;
	double epsilon_0 = 4.5;	// (MeV)
	
	if(A<148)  k = 1.0;
	if(A>=148) k = 1 + 0.09*(A-148)*(A-148)*exp(-0.18*(A-148));
	
	Kappa = k + (1.0-k)*(E_g-epsilon_0)/(E_r-epsilon_0);
	Kappa_0 = k + (k-1.)*(epsilon_0)/(E_r-epsilon_0);
	Gamma_k = Kappa*Gamma_r*(pow(E_g,2.0) + (2.0*Pi*temp)*(2.0*Pi*temp))/pow(E_r,2.0);
	Gamma_k0 = Kappa_0*Gamma_r*((2.*Pi*temp)*(2.*Pi*temp))/pow(E_r,2.0);
	denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + pow(E_g,2.)*pow(Gamma_k,2.);
	
	EGLo = factor*sigma_r*Gamma_r*((E_g*Gamma_k)/denominator + 0.7*Gamma_k0/pow(E_r,3.));
	return EGLo;
}

// Enhanced Generalized Lorentzian, variable temperature
double EnGeLorTemp(double E_g, double E_r, double Gamma_r, double sigma_r){

	double GLo = 0.;
	double GLo_var = 0.;
	double Gamma_k = 0.;
	double Gamma_k0 = 0.;
	double denominator = 0.;
	
// VARIABLE TEMPERATURE			
	int i,j,k;
	double dEx      = 0.01;	// steps of 10 keV in excitation energy
	double Ex_final = 0.;	// final excitation energy
	double U        = 0.;	// intrinsic excitation energy
	double T        = 0.;	// variable temperature
	int fsum        = 0;	// sum of strength; to make average
	double E_low    = 0.;
	double E_high   = 0.;
	
	E_low  = Ex_min - E_g;	// Final Ex low
	E_high = Ex_max - E_g;	// Final Ex high
	
	if(E_g > Ex_max) return 0;	// no gammas bigger than max excitation energy
	
	if(E_low < 0.0) E_low = 0.0; 
	i = (int) (E_low/dEx);
	j = (int) (E_high/dEx);
	
	if(E_g <= Ex_max){
		for(k=i; k<=j; k++){
			Ex_final = dEx*(double) k;
			U = Ex_final - E_1;
			if(U < 0.0) U = 0.000001;
			T = sqrt(U/a);
			if(T < 0.0) T = 0.0;
//			T = 0.;
			
			Gamma_k = Gamma_r*(pow(E_g,2.0) + (2.0*Pi*T)*(2.0*Pi*T))/pow(E_r,2.0);
			Gamma_k0 = Gamma_r*((2.*Pi*T)*(2.*Pi*T))/pow(E_r,2.0);	// Gamma_k for E_gamma=0
			denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + pow(E_g,2.)*pow(Gamma_k,2.);
	
			GLo_var += factor*sigma_r*Gamma_r*((E_g*Gamma_k)/denominator + 0.7*Gamma_k0/pow(E_r,3.));

			fsum   += 1;
		}
	
		GLo_var = GLo_var/(double)fsum;	// average KMF strength function for a certain E_g over the whole Ex range used in RhoSigChi
		return GLo_var;
	}

/*	double EGLo = 0.;
	double Kappa = 0.;
	double Kappa_0 = 0.;
	double k = 0.;
	double Gamma_k = 0.;
	double Gamma_k0 = 0.;
	double denominator = 0.;
	double epsilon_0 = 4.5;	// (MeV)
	if(E_g > Ex_max) return 0.;	// no gammas bigger than max excitation energy
	
	if(A<148)  k = 1.0;
	if(A>=148) k = 1 + 0.09*(A-148)*(A-148)*exp(-0.18*(A-148));
	
	Kappa = k + (1.0-k)*(E_g-epsilon_0)/(E_r-epsilon_0);
	Kappa_0 = k + (k-1.)*(epsilon_0)/(E_r-epsilon_0);
	
	int i,j,n;
	double dEx      = 0.01;	// steps of 10 keV in excitation energy
	double Ex_final = 0.;	// final excitation energy
	double U        = 0.;	// intrinsic excitation energy
	double T        = 0.;	// temperature
	int fsum        = 0;	// sum of strength; to make average
	double E_low    = 0.;
	double E_high   = 0.;
	
	E_low  = Ex_min - E_g;	// Final Ex low
	E_high = Ex_max - E_g;	// Final Ex high
	
	if(E_low < 0.0) E_low = 0.0; 
	i = (int) (E_low/dEx);
	j = (int) (E_high/dEx);
	
	for(n=i; n<=j; n++){
		Ex_final = dEx*(double) n;
		U = Ex_final - E_1;
		if(U < 0.0) U = 0.000001;
		T = sqrt(U/a);
		if(T < 0.0) T = 0.0;
	
		Gamma_k = Kappa*Gamma_r*(pow(E_g,2.0) + (2.0*Pi*T)*(2.0*Pi*T))/pow(E_r,2.0);
		Gamma_k0 = Kappa_0*Gamma_r*((2.*Pi*T)*(2.*Pi*T))/pow(E_r,2.0);
		denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + pow(E_g,2.)*pow(Gamma_k,2.);
		EGLo   += factor*sigma_r*Gamma_r*((E_g*Gamma_k)/denominator + 0.7*Gamma_k0/pow(E_r,3.));
		fsum   += 1;
	}
	
	EGLo = EGLo/(double)fsum;	// average KMF strength function for a certain E_g over the whole Ex range used in RhoSigChi
	
	return EGLo;
*/	
}


// Modified Lorentzian, spherical (beta2 == 0)
double ModLor(double E_g, double E_r, double Gamma_r, double sigma_r){
	double MLo = 0.;
	double L = 0.;
	double c = 5.39E-03;	// See corrections to RIPL-2 from V. Plujko
	double C_coll = 0.;
	double F = 1.0;		// sigma(np)/sigma_f(np), see RIPL2
	double k_s = 0.;
	double k_r = 0.0;	
	double k_0 = 0.3;
	double Gamma_c = 0.;
	double Gamma_F = 0.;
	double Gamma_w = 32.85*pow(A,(-1./3.));
	double Gamma_tot = 0.;
	double denominator = 0.;
	
	C_coll = c*F;
	k_r = (Gamma_r - C_coll*E_r*E_r)/Gamma_w;
	if(E_g < 2.*E_r)  k_s = k_r + (k_0 - k_r)*fabs((E_g - E_r)/E_r);
	if(E_g >= 2.*E_r) k_s = k_0;
	
	Gamma_c = C_coll*E_r*(E_g + a*pow(temp,2));
	Gamma_F = k_s*Gamma_w;
	Gamma_tot = Gamma_c + Gamma_F;

	L = 1./(1. - exp(-E_g/temp));
	denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + pow(E_g,2.)*pow(Gamma_tot,2.);
	
	MLo = factor*L*sigma_r*Gamma_r*E_g*Gamma_tot/denominator;
	
	return MLo;
}

// Modified Lorentzian, deformed (beta2 different from 0)
double ModLorDef(double E_g, double E_r1, double Gamma_r1, double sigma_r1, double E_r2, double Gamma_r2, double sigma_r2){
	double MLo1 = 0.;
	double MLo2 = 0.;
	double MLo_tot = 0.;	
	double L = 0.;
	double c = 5.39E-03;	// See corrections to RIPL-2 from V. Plujko
	double C_coll = 0.;
	double F = 1.0;		// sigma(np)/sigma_f(np), see RIPL2
	double k_s1 = 0.;
	double k_s2 = 0.;	
	double k_r1 = 0.0;	
	double k_r2 = 0.0;	

	double k_0 = 0.3;
	double Gamma_c1 = 0.;
	double Gamma_c2 = 0.;	
	double Gamma_F1 = 0.;
	double Gamma_F2 = 0.;	
	double Gamma_w = 32.72*pow(A,(-1./3.));
	double Gamma_s1 = Gamma_w/pow(a0,1.6);
	double Gamma_s2 = Gamma_w/pow(b0,1.6);	 
	double Gamma_tot1 = 0.;
	double Gamma_tot2 = 0.;
	double denominator1 = 0.;
	double denominator2 = 0.;
	
	C_coll = c*F;
	
	k_r1 = (Gamma_r1 - C_coll*E_r1*E_r1)/Gamma_w;
	if(E_g < 2.*E_r1)  k_s1 = k_r1 + (k_0 - k_r1)*fabs((E_g - E_r1)/E_r1);
	if(E_g >= 2.*E_r1) k_s1 = k_0;
	
	k_r2 = (Gamma_r2 - C_coll*E_r2*E_r2)/Gamma_w;
	if(E_g < 2.*E_r2)  k_s2 = k_r2 + (k_0 - k_r2)*fabs((E_g - E_r2)/E_r2);
	if(E_g >= 2.*E_r2) k_s2 = k_0;	
	
	Gamma_c1 = C_coll*E_r1*(E_g + a*pow(temp,2));
	Gamma_c2 = C_coll*E_r2*(E_g + a*pow(temp,2));
	
	Gamma_F1 = k_s1*Gamma_s1;
	Gamma_F2 = k_s2*Gamma_s2;
	
	Gamma_tot1 = Gamma_c1 + Gamma_F1;
	Gamma_tot2 = Gamma_c2 + Gamma_F2;	

	L = 1./(1. - exp(-E_g/temp));
	denominator1 = (pow(E_g,2.) - pow(E_r1,2.))*(pow(E_g,2.) - pow(E_r1,2.)) + pow(E_g,2.)*pow(Gamma_tot1,2.);
	denominator2 = (pow(E_g,2.) - pow(E_r2,2.))*(pow(E_g,2.) - pow(E_r2,2.)) + pow(E_g,2.)*pow(Gamma_tot2,2.);
	
	MLo1 = factor*L*sigma_r1*Gamma_r1*E_g*Gamma_tot1/denominator1;
	MLo2 = factor*L*sigma_r2*Gamma_r2*E_g*Gamma_tot2/denominator2;
	MLo_tot = MLo1 + MLo2;
	
	return MLo_tot;
}
 

// Generalized Fermi Liquid (Mughabghab PLB 487, 155 (2000), and RIPL-2)	
double GenFerLiq(double E_g, double E_r, double Gamma_r, double sigma_r){
	double GFL = 0.;
	double K_GFL = 0.63;
	double Gamma_m = 0.;
	double Gamma_C = 0.;
	double Gamma_dq = 0.;
	double Gamma_dq_0 = 0.;
	
	double C_dq = 1.05;
	double C_f = 0;
	double s2 = 217.16/pow(A,2.);	// systematics
//	s2 = E2*beta2*beta2;			// if E2 and beta2 is known, else, use systematics
	double denominator = 0.;	

	Gamma_dq = C_dq*sqrt(pow(E_g,2.)*pow(beta2,2.) + E_g*s2);
	Gamma_dq_0 = C_dq*sqrt(pow(E_r,2.)*pow(beta2,2.) + E_r*s2);
	
	C_f = (Gamma_r - Gamma_dq_0)/pow(E_r,2.);	// From the claim that Gamma_m(E_g=E_r,temp=0) = Gamma_r (see RIPL-2, p.123)
	
	Gamma_C  = C_f*(pow(E_g,2.) + 4.*pow(Pi,2.)*pow(temp,2.));
	Gamma_m = Gamma_C + Gamma_dq;
	denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + K_GFL*pow(E_g,2.)*pow(Gamma_m,2.);
//	denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)); //original from Mughabghab, diverges as E_g -> E_r
	
	GFL = factor*sigma_r*Gamma_r*K_GFL*E_g*Gamma_m/denominator;
	
	return GFL;
}	

// Hybrid Formula, Goriely PLB 436, 10 (1998) and RIPL-2
double HybridFormulaE1(double E_g, double E_r, double Gamma_r, double sigma_r){
	double HybFormE1 = 0.0;
	double F1prime   = -0.04;
	double F0prime   = 1.49;
	double F1        = 0.12;
	double Landau    = sqrt((1+(0.667*F1prime))/(1+(2*F1)));
	double Gamma_h   = Landau*Gamma_r*(pow(E_g,2.)+4.*pow(Pi,2.)*pow(temp,2.))/E_r/E_g;
	double denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + pow(E_g,2.)*Gamma_r*Gamma_h;
	
	HybFormE1 = 1.2*factor*sigma_r*Gamma_r*Gamma_h*E_g/denominator;
	
	return HybFormE1;
}

// Kadmenskij, Markushev and Furman model, Yad. Fiz. 37, 277-283 (1983), here with constant temperature
double KMF(double E_g, double E_r, double Gamma_r, double sigma_r, double temp){
	double f_KMF = 0.;
	double Gamma_k = 0.;
	double denominator = 0.;
	
	Gamma_k = Gamma_r*Gamma_r*(pow(E_g,2.0) + (2.0*Pi*temp)*(2.0*Pi*temp))/E_r;
	denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.));
	
	f_KMF = factor*0.7*sigma_r*Gamma_k/denominator;
	return f_KMF;
}	

// Kadmenskij, Markushev and Furman model, Yad. Fiz. 37, 277-283 (1983), with variable temperature 
double KMF_var(double E_g, double E_r, double Gamma_r, double sigma_r){
	if(E_g > Ex_max) return 0.;	// no gammas bigger than max excitation energy

	double f_KMF = 0.;
	double Gamma_k = 0.;
	double denominator = 0.;
	double diff = abs(pow(E_g,2.) - pow(E_r,2.));
	denominator = diff*diff;
	
	int i,j,k;
	double dEx      = 0.01;	// steps of 10 keV in excitation energy
	double Ex_final = 0.;	// final excitation energy
	double U        = 0.;	// intrinsic excitation energy
	double T        = 0.;	// temperature
	int fsum        = 0;	// sum of strength; to make average
	double E_low    = 0.;
	double E_high   = 0.;
	
	E_low  = Ex_min - E_g;	// Final Ex low
	E_high = Ex_max - E_g;	// Final Ex high
	
	if(E_low < 0.0) E_low = 0.0; 
	i = (int) (E_low/dEx);
	j = (int) (E_high/dEx);
	
	for(k=i; k<=j; k++){
		Ex_final = dEx*(double) k;
		U = Ex_final - E_1;
		if(U < 0.0) U = 0.000001;
		T = sqrt(U/a);
		if(T < 0.0) T = 0.0;
//		cout << "Temperature T_f = " << T << " MeV at Ex_final = " << Ex_final << " MeV. " << endl;

		Gamma_k = Gamma_r*Gamma_r*(pow(E_g,2.0) + (2.0*Pi*T)*(2.0*Pi*T))/E_r;
		f_KMF  += factor*0.7*sigma_r*Gamma_k/denominator;
		fsum   += 1;
	}
	
	f_KMF = f_KMF/(double)fsum;	// average KMF strength function for a certain E_g over the whole Ex range used in RhoSigChi
	

	return f_KMF;		
}

// Giant Magnetic Dipole Resonance due to spin-flip excitations (see RIPL-1)
double SpinFlipM1(double E_g, double E_r, double Gamma_r, double sigma_r){
	double f_M1     = 0.;
	double E_M1     = 41.*pow(A,(-1./3.)) + 0.0;  //REMEMBER TO REMOVE DIRTY +0.0 MeV in energy
	double Gamma_M1 =  4.0  ; //REMEMBER TO REMOVE DURTY -0.0 TO GIVE 4.0
	double R1       = 0.0588*pow(A,0.878);
	double E1_7     = (49.0 - pow(E_r,2.0))*(49.0 - pow(E_r,2.0))*E_r;
	double E1_7MeV  = 0.7*sigma_r*pow(Gamma_r,2.)*49.0/E1_7;
	double R2       = (49.0 - pow(E_M1,2.))*(49.0 - pow(E_M1,2.))/(7.0*Gamma_M1)/(7.0*Gamma_M1);
	double M1_7     = 1.0/(R2 + 1.0);
	double M1_7MeV  = M1_7/7.0;
	double sigma_M1 = 1.0 * E1_7MeV/M1_7MeV/R1;  //REMEMBER TO REMOVE DIRTY 0.5 TO GIVE RATIO 7.05
    cout << "Ratio E1/M1 = " << R1 << " Centroid = " << E_M1 << " M1 width = " << Gamma_M1 << " sigma M1 = " << sigma_M1 <<endl;

	
//	E_M1 = 7.5;
//	sigma_M1 = 0.8;
//	Gamma_M1 = 4.0;
	
	double denominator = (pow(E_g,2.) - pow(E_M1,2.))*(pow(E_g,2.) - pow(E_M1,2.)) + pow(E_g,2.)*pow(Gamma_M1,2.);
	
	f_M1 = factor*sigma_M1*pow(Gamma_M1,2.)*E_g/denominator;
	
//	cout << "E_M1 = " << E_M1 << ", sigma_M1 = " << sigma_M1 << ", Gamma_M1 = " << Gamma_M1 << endl;
    
    f_M1= 0.;       /////////////TAKE CARE!
	
	return f_M1;
}

// Giant Electric Quadrupole Resonance, from RIPL-1 with corrections (see next comment)	
double IsoScalarE2(double E_g){

	/* Isoscalar E2 model. Values taken from RIPL1 page 103 with corrections:
	0.00014->0.00015=2*pi*alpha*r**2/(5*mc**2), alpha=fine structure constant,
	r=1.233 fm, m=nucleon mass (=1u). Value given in mb/MeV. Further, it should 
	be E_E2**2 in the formula for sigma. Finally, it is -0.012 in the slope for 
	Gamma */
	
	double f_E2        = 0.;
	double E_E2        = 63.*pow(A,(-1./3.));
	double Gamma_E2    = 6.11 - (0.012*A);
	double sigma_E2    = 0.00015*pow(Z,2.)*pow(E_E2,2.)/(pow(A,1./3.)*Gamma_E2);
    /* Using f_e2 as unknown E1 resonance to fake the Varlamox 2007 data*/
    E_E2     =    7.5;      //7.0;
    Gamma_E2 =    1.4;      //1.8;
    sigma_E2 =    60.;      //4.;
	double denominator = (pow(E_g,2.) - pow(E_E2,2.))*(pow(E_g,2.) - pow(E_E2,2.)) + pow(E_g,2.)*pow(Gamma_E2,2.);
//	f_E2 =(3./5.)*factor*sigma_E2*pow(Gamma_E2,2.)*E_g/denominator;  To be used for E2!!!
  f_E2 =(1.)*factor*sigma_E2*pow(Gamma_E2,2.)*E_g/denominator; //fake for unknown E1
//    f_E2 = sigma_E2*factor*(1./(sqrt(2.*Pi)*Gamma_E2))*exp(-pow((E_g-E_E2),2.)/(2.*pow(Gamma_E2,2.)));
    
    // ONCE MORE
    
    E_E2     =    5.5;      //7.0;
    Gamma_E2 =    0.7;      //1.8;
    sigma_E2 =    50.;      //4.;
	    denominator = (pow(E_g,2.) - pow(E_E2,2.))*(pow(E_g,2.) - pow(E_E2,2.)) + pow(E_g,2.)*pow(Gamma_E2,2.);
    //	f_E2 =(3./5.)*factor*sigma_E2*pow(Gamma_E2,2.)*E_g/denominator;  To be used for E2!!!
    f_E2 =f_E2+(1.)*factor*sigma_E2*pow(Gamma_E2,2.)*E_g/denominator; //fake for unknown E1
    //    f_E2 = sigma_E2*factor*(1./(sqrt(2.*Pi)*Gamma_E2))*exp(-pow((E_g-E_E2),2.)/(2.*pow(Gamma_E2,2.)));
    
    
    
    
    
	return f_E2;
}

double PygmyL1(double E_g){
	double pygmyStrength = 0.;
	double E_pygmy            = 6.20;	// position of pygmy resonance (MeV)
	double Gamma_pygmy		  = 2.2;	// Lorentzian width
	double peak_pygmy         = 200;	// Lorentzian peak cross section (mb)
	// assuming a Lorentzian pygmy resonance
	double denominator = (pow(E_g,2.0) - pow(E_pygmy,2.0))*(pow(E_g,2.0) - pow(E_pygmy,2.0)) + pow(Gamma_pygmy,2.0)*pow(E_g,2.0);
	pygmyStrength = factor*peak_pygmy*pow(Gamma_pygmy,2.0)*E_g/denominator;
	return pygmyStrength;
}		

double PygmyL2(double E_g){
	double pygmyStrength = 0.;
	double E_pygmy            = 2.65;	// position of pygmy resonance (MeV) 
	double Gamma_pygmy		  = 0.60;	// Lorentzian width
	double peak_pygmy         = 0.60;	// Lorentzian peak cross section (mb)
	// assuming a Lorentzian pygmy resonance
	double denominator = (pow(E_g,2.0) - pow(E_pygmy,2.0))*(pow(E_g,2.0) - pow(E_pygmy,2.0)) + pow(Gamma_pygmy,2.0)*pow(E_g,2.0);
	pygmyStrength = factor*peak_pygmy*pow(Gamma_pygmy,2.0)*E_g/denominator;
	return pygmyStrength;
}		

double PygmyH1(double E_g){
	double pygmyStrength = 0.;
	double E_pygmy            = 1.95;	// position of pygmy resonance (MeV) 
	double Gamma_pygmy		  = 0.90;	// Lorentzian width
	double peak_pygmy         = 0.60;	// Lorentzian peak cross section (mb)
	// assuming a Lorentzian pygmy resonance
	double denominator = (pow(E_g,2.0) - pow(E_pygmy,2.0))*(pow(E_g,2.0) - pow(E_pygmy,2.0)) + pow(Gamma_pygmy,2.0)*pow(E_g,2.0);
	pygmyStrength = factor*peak_pygmy*pow(Gamma_pygmy,2.0)*E_g/denominator;
	return pygmyStrength;
}		

double PygmyH2(double E_g){
	double pygmyStrength = 0.;
	double E_pygmy            = 2.65;	// position of pygmy resonance (MeV) 
	double Gamma_pygmy		  = 0.70;	// Lorentzian width
	double peak_pygmy         = 0.50;	// Lorentzian peak cross section (mb)
	// assuming a Lorentzian pygmy resonance
	double denominator = (pow(E_g,2.0) - pow(E_pygmy,2.0))*(pow(E_g,2.0) - pow(E_pygmy,2.0)) + pow(Gamma_pygmy,2.0)*pow(E_g,2.0);
	pygmyStrength = factor*peak_pygmy*pow(Gamma_pygmy,2.0)*E_g/denominator;
	return pygmyStrength;
}		
		
double PygmyT1(double E_g){
	double pygmyStrength = 0.;
	double E_pygmy            = 1.95;	// position of pygmy resonance (MeV)
	double Gamma_pygmy		  = 0.61;	// Lorentzian width
	double peak_pygmy         = 0.41;	// Lorentzian peak cross section (mb)
	// assuming a Lorentzian pygmy resonance
	double denominator = (pow(E_g,2.0) - pow(E_pygmy,2.0))*(pow(E_g,2.0) - pow(E_pygmy,2.0)) + pow(Gamma_pygmy,2.0)*pow(E_g,2.0);
	pygmyStrength = factor*peak_pygmy*pow(Gamma_pygmy,2.0)*E_g/denominator;
	return pygmyStrength;
}		
		
double PygmyT2(double E_g){
	double pygmyStrength = 0.;
	double E_pygmy            = 2.48;	// position of pygmy resonance (MeV)
	double Gamma_pygmy		  = 0.90;	// Lorentzian width
	double peak_pygmy         = 0.49;	// Lorentzian peak cross section (mb)
	// assuming a Lorentzian pygmy resonance
	double denominator = (pow(E_g,2.0) - pow(E_pygmy,2.0))*(pow(E_g,2.0) - pow(E_pygmy,2.0)) + pow(Gamma_pygmy,2.0)*pow(E_g,2.0);
	pygmyStrength = factor*peak_pygmy*pow(Gamma_pygmy,2.0)*E_g/denominator;
	return pygmyStrength;
}		
		
double ModGenLor(double E_g, double E_r, double Gamma_r, double sigma_r, double temp){	
// modified, Generalized Lorentzian. THe width is modified to fit the pygmy in 232Th
	
	double GLo = 0.;
	double Gamma_k = 0.;
	double Gamma_k0 = 0.;
	double denominator = 0.;
	double delta = 0.5;
	double front_factor = Gamma_r/pow(E_r,2.);
	double little_denominator = pow(E_g,0.9) + delta;
//	double little_denominator = pow(E_g,1.) + delta;
	
	Gamma_k = front_factor*(pow(E_g,2.0) + (2.0*Pi*temp)*(2.0*Pi*temp)*E_r/little_denominator);
	Gamma_k0 = front_factor*(2.*Pi*temp)*(2.*Pi*temp)*E_r/little_denominator;	// Gamma_k for E_gamma=0
	denominator = (pow(E_g,2.) - pow(E_r,2.))*(pow(E_g,2.) - pow(E_r,2.)) + pow(E_g,2.)*pow(Gamma_k,2.);
	
	GLo = factor*sigma_r*Gamma_r*((E_g*Gamma_k)/denominator + 0.7*Gamma_k0/pow(E_r,3.)); 
	return GLo;
	
}	

	
		
		
		
		
		
		
		
		
		
		
		
	
	