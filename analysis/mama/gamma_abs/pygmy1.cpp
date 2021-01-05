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
        
        
	TH2F *h = new TH2F("h"," ", 10,0.0,  4.2,10,-0.09e-7,0.599e-07);

	const float factor = 8.674E-08;	// const. factor in mb^(-1) MeV^(-2)
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
	TGraphErrors *strengthexp = new TGraphErrors(48,energy,strength,energyerr,strengtherr);
    
   
    ///////////////////////////////////////////////////////	
    
	c1.SetRightMargin(0.01);
	c1.SetLeftMargin(0.15);
	c1.SetBottomMargin(0.18);
	c1.SetTopMargin(0.15);
    
    h.GetYaxis().CenterTitle();
	h.GetYaxis().SetTitleOffset(0.9);
	h.GetYaxis().SetTitle("#gammaSF (MeV^{-3})");
    h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
//    h->GetXaxis()->SetLabelOffset(1.8);

	h.GetXaxis().CenterTitle();
	h.GetXaxis().SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
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
    t->SetTextSize(0.05);
    t->DrawLatex(0.3,50E-09,"(d,p)^{238}Np");
    t->DrawLatex(1.7,50E-09,"E = 3.0 - 5.7 MeV");
    
    t.SetTextSize(0.045);
    t->DrawLatex(0.3,35E-09,"#omega = 2.00(5) MeV");
    t->DrawLatex(0.3,29E-09,"#Gamma = 0.70(5) MeV");
    t->DrawLatex(0.3,23E-09,"#sigma = 0.70(5) mb");
    t->DrawLatex(3.2,35E-09,"2.55(5) MeV");
    t->DrawLatex(3.2,29E-09,"0.70(10) MeV");
    t->DrawLatex(3.2,23E-09,"0.50(5) mb");

	
	c1->Update();
	c1->Print("pygmy1.pdf");
    c1->Print("pygmy1.png");
}
