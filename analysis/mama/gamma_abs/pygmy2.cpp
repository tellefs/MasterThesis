{
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle.SetOptStat(0);
	gStyle.SetFillColor(0);
	gStyle.SetPadBorderMode(0);
    
	m = (TH1F*)gROOT->FindObject("h");
	if (m) m->Delete();
	TCanvas *c1 = new TCanvas("c1","Gamma-ray strength function",900,600);
        
        
	TH2F *h = new TH2F("h"," ", 10,0.0,  4.2,10,-0.09e-7,0.45e-07);

	const double factor = 8.674E-08;	// const. factor in mb^(-1) MeV^(-2)
	int    i, j;
	float x, y, z, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12;
	float bg, eg;
    int i1, i2;
    
    ///////////////////////////////////////////////////////	
    float teoenergy[200];
    float teoT[200], teoT1[200], teoT2[200];
    float back[200];
    float a0_teo =   0.1;
	float a1_teo =   0.1;
    
	i = 0;
    ifstream teofile("GDR_models_238Np.dat");
	while(teofile){
		teofile  >> x >> y >> z >> d1 >> d2 >>  
        d3 >> d4 >> d5 >>
        d6 >> d7 >> d8 >>
        d9 >> d10 >> d11 >>
        d12;
		if(i<200){
			teoenergy[i] = x;
			teoT1[i] = d3;
			teoT2[i] = d4;
            teoT[i]  = d5;
			back[i]  = d1;
        }
        i++;
	}
	TGraph *strengthteoT = new TGraph(169,teoenergy,teoT);
    TGraph *strengthteoT1 = new TGraph(169,teoenergy,teoT1);
    TGraph *strengthteoT2 = new TGraph(169,teoenergy,teoT2);
    
	

///////////////////////////////////////////////////////	
    float strength[60],strengtherr[60],energy[60],energyerr[60];
    float a0 =   0.117371;
	float a1 =   0.121372;
	i = 0;
    ifstream strengthfile("../strength.nrm");
	while(strengthfile){
		strengthfile >> x;
		if(i<48){
			if(x == 0.) x = -1.;
                eg = a0 + (a1*i);
                i1 = (int)(((eg-a0_teo)/a1_teo)+0.5);
                i2 = i1 + 1;
            bg = back[i1] + (back[i2]-back[i1])*((eg -(a0_teo + a1_teo*i1))/(a1_teo));
            
                if(x !=0)strength[i] = x - bg;
                    energy[i] = a0 + (a1*i);
                    energyerr[i] = 0.0;
                    }	
		else{strengtherr[i-48] = x;}
		i++;
	}
    energy[8]=20; //getting away first point
	TGraphErrors *strengthexp = new TGraphErrors(48,energy,strength,energyerr,strengtherr);
    
    ///////////////    MAKING CHI**2   /////////////////////
    double E_1          = 0;	// position of pygmy resonance (MeV)
    double G_1          = 0;	// Lorentzian width
    double S_1          = 0;
    double E_2          = 0;	// position of pygmy resonance (MeV)
    double G_2          = 0;	// Lorentzian width
    double S_2          = 0;
    double E_1p          = 1.95;//+.04;	// position of pygmy resonance (MeV)
    double G_1p          = 0.61;//+0.05;	// Lorentzian width
    double S_1p          = 0.40;
    double E_2p          = 2.47;//+0.06;	// position of pygmy resonance (MeV)
    double G_2p          = 0.91;//+0.1;	// Lorentzian width
    double S_2p          = 0.49;

    double E_1b         = 0;	// position of pygmy resonance (MeV)
    double G_1b         = 0;	// Lorentzian width
    double S_1b         = 0;
    double E_2b         = 0;	// position of pygmy resonance (MeV)
    double G_2b         = 0;	// Lorentzian width
    double S_2b         = 0;
    
    double teo;
    int k1,k2,k3,k4,k5,k6;
    int j1 = 9, j2 = 31;
    float sum = 0., chi2 = 100000., chi2b=100000;
    int ktot = 4; //10/2;
    float stepE1 = 0.01,stepE2= 0.01, stepG1=0.01, stepG2=0.01, stepS1=0.01, stepS2=0.01;
    for (k1=-ktot;k1<=ktot;k1++){
        E_1=E_1p+k1*stepE1;
        for (k2=-ktot;k2<=ktot;k2++){
            G_1=G_1p+k2*stepG1;
            for (k3=-ktot;k3<=ktot;k3++){
                S_1=S_1p+k3*stepS1;
                for (k4=-ktot;k4<=ktot;k4++){
                    E_2=E_2p+k4*stepE2;
                    for (k5=-ktot;k5<=ktot;k5++){
                        G_2=G_2p+k5*stepG2;
                        for (k6=-ktot;k6<=ktot;k6++){
                            S_2=S_2p+k6*stepS2;
                            sum=0.;
    for(i=j1; i<=j2; i++){
        eg = a0 + (a1*i);
        i1 = (int)(((eg-a0_teo)/a1_teo)+0.5);
        i2 = i1 + 1;
//        teo = teoT[i1] + (teoT[i2]-teoT[i1])*((eg -(a0_teo + a1_teo*i1))/(a1_teo));

        teo =(factor*S_1*pow(G_1,2.0)*eg/((pow(eg,2.0) - pow(E_1,2.0))*(pow(eg,2.0) - pow(E_1,2.0)) + pow(G_1,2.0)*pow(eg,2.0)))  +   (factor*S_2*pow(G_2,2.0)*eg/((pow(eg,2.0) - pow(E_2,2.0))*(pow(eg,2.0) - pow(E_2,2.0)) + pow(G_2,2.0)*pow(eg,2.0)));
        sum = sum + ((strength[i]-teo)*(strength[i]-teo)/(strengtherr[i]*strengtherr[i]));
//        printf(" cha =%d  eg = %f exp= %e teo = %e \n",i, eg, strength[i],teo);
    }
    chi2 = sum/((float)(j2-j1+1-6));
    if(chi2<chi2b){
        chi2b=chi2;
        E_1b=E_1;
        G_1b=G_1;
        S_1b=S_1;
        E_2b=E_2;
        G_2b=G_2;
        S_2b=S_2;
    }
//                            printf("chi %f  chi-best %f \n",chi2,chi2b);
                        }
                    }
                }
            }
        }
    }
    
    
    printf("Chi**2 = %12.3f\n",chi2b);
    printf("E_1b= %f G_1b= %f S_1b= %f\n",E_1b,G_1b,S_1b);
    printf("E_2b= %f G_2b= %f S_2b= %f\n",E_2b,G_2b,S_2b);

    ///////////////////////////////////////////////////////
    
    
    //////////////////////////////////////////////
    
	c1.SetRightMargin(0.01);
	c1.SetLeftMargin(0.15);
	c1.SetBottomMargin(0.18);
	c1.SetTopMargin(0.15);
    
    h.GetYaxis().CenterTitle();
	h.GetYaxis().SetTitleOffset(0.9);
	h.GetYaxis().SetTitle("#gammaSF (MeV^{-3})");
    h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetLabelSize(0.06);
//    h->GetXaxis()->SetLabelOffset(1.8);

	h.GetXaxis().CenterTitle();
	h.GetXaxis().SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitleOffset(1.2);

	h.Draw();
    strengthexp->SetMarkerStyle(21);
	strengthexp->SetMarkerSize(0.9);
	strengthexp->Draw("P");
    
    strengthteoT->SetLineColor(2);
    strengthteoT->SetLineWidth(2);
	strengthteoT->Draw("L");
    
    strengthteoT1->SetLineColor(1);
    strengthteoT1->SetLineWidth(2);
    strengthteoT1->SetLineStyle(7);
	strengthteoT1->Draw("L");
    
    strengthteoT2->SetLineColor(1);
    strengthteoT2->SetLineWidth(2);
    strengthteoT2->SetLineStyle(7);
	strengthteoT2->Draw("L");
    
    TLatex t;
    t->SetTextSize(0.06);
    t->DrawLatex(0.3,35E-09,"(d,p)^{238}Np");
    t->DrawLatex(1.7,35E-09,"E = 3.0 - 5.7 MeV");
    
//    t.SetTextSize(0.045);
//    t->DrawLatex(0.3,35E-09,"#omega = 1.97(4) MeV");
//    t->DrawLatex(0.3,29E-09,"#Gamma = 0.59(5) MeV");
//    t->DrawLatex(0.3,23E-09,"#sigma = 0.39(4) mb");
//    t->DrawLatex(3.2,35E-09,"2.51(6) MeV");
//    t->DrawLatex(3.2,29E-09,"0.79(10) MeV");
//    t->DrawLatex(3.2,23E-09,"0.43(6) mb");

	
	c1->Update();
	c1->Print("pygmy2.pdf");
    c1->Print("pygmy2.png");
}
