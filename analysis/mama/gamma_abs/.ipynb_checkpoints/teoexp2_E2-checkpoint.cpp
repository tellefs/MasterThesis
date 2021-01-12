void teoexp2_E2(){
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetFillColor(0);
    gStyle->SetPadBorderMode(0);
    /*
	m = (TH1F*)gROOT->FindObject("h");
	if (m) m->Delete();
	*/
	TCanvas *c1 = new TCanvas("c1","Gamma-ray strength function",900,700);
	TH2F *h = new TH2F("h"," ",10,0.0,  17.2,10,5.0e-9,19e-06);

	const double factor = 8.674E-08;	// const. factor in mb^(-1) MeV^(-2)
	int    i, j;
	double x, y, z, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12;	
    
    ///////////////////////////////////////////////////////	
    double teo[200], fM1[200],fE2[200], teoenergy[200];
	i = 0;
    ifstream teofile("GDR_models_233U.dat");
	while(teofile){
		teofile  >> x >> y >> z >> d1 >> d2 >>  
        d3 >> d4 >> d5 >>
        d6 >> d7 >> d8 >>
        d9 >> d10 >> d11 >>
        d12;
        
		if(i<200){
			teoenergy[i] = x;
			teo[i] = d1;
            fM1[i] = d2;
            fE2[i] = y;
        }
        i++;
	}
	TGraph *strengthteo = new TGraph(199,teoenergy,teo);
    TGraph *strengthfM1 = new TGraph(199,teoenergy,fM1);
    TGraph *strengthfE2 = new TGraph(199,teoenergy,fE2);
	
///////////////////////////////////////////////////////	
    double strength[60],strengtherr[60],energy[60],energyerr[60];
    double a0 =   0.117371;
	double a1 =   0.121372;
	i = 0;
    ifstream strengthfile("../strength.nrm");
	while(strengthfile){
		strengthfile >> x;
		if(i<48){
			strength[i] = x;
			energy[i] = a0 + (a1*i);
			energyerr[i] = 0.0;
		}	
		else{strengtherr[i-48] = x;}
		i++;
	}
    energy[7]=20;
    energy[8]=20; //getting away first points

	TGraphErrors *strengthexp = new TGraphErrors(48,energy,strength,energyerr,strengtherr);

    
///////////////////////////////////////////////////////	
 /*   
    double gdr[70],gdr_err[70],gdr_energy[70],gdr_energyerr[70];
	i =	0;
    ifstream gdrfile("237Np_berman_1986x.txt");
    
	while(gdrfile){
		gdrfile >> x >> y >> z;
		if(i<27){
            gdr_energy[i]	= x;
			gdr[i]			= y;
			gdr_err[i]		= z;
			gdr_energyerr[i]= 0;
		}
		i++;
	}
	for(i = 0; i < 27; i++){
		gdr[i]     = factor*gdr[i]    /gdr_energy[i];
        gdr_err[i] = factor*gdr_err[i]/gdr_energy[i];
	}
    gdr[0]=0; //gdr[1]=0; gdr[2]=0;
    gdr_energy[0]=20; //gdr_energy[1]=20; gdr_energy[2]=20;
    
	TGraphErrors *gdrexp = new TGraphErrors(27,gdr_energy,gdr,gdr_energyerr,gdr_err);
    
    
    */
    ///////////////////////////////////////////////////////
    /*
    
    double gdr[70],gdr_err[70],gdr_energy[70],gdr_energyerr[70];
	i =	0;
    ifstream gdrfile("238U_caldwell_1980x.txt");
    
	while(gdrfile){
		gdrfile >> x >> y >> z;
		if(i<59){
            gdr_energy[i]	= x;
			gdr[i]			= y;
			gdr_err[i]		= z;
			gdr_energyerr[i]= 0;
		}
		i++;
	}
	for(i = 0; i < 59; i++){
		gdr[i]     = factor*gdr[i]    /gdr_energy[i];
        gdr_err[i] = factor*gdr_err[i]/gdr_energy[i];
	}
//    gdr[0]=0; gdr[1]=0; gdr[2]=0;
//    gdr_energy[0]=20; gdr_energy[1]=20; gdr_energy[2]=20;
    
	TGraphErrors *gdrexp1 = new TGraphErrors(59,gdr_energy,gdr,gdr_energyerr,gdr_err);
    */
//////////////////////////////KOPECKY/////////////////	

        double fEx[1]={12.74E-08};
        double fEerr[1]={3.14E-08};
        double fEene[1]={4.1};
        double fEeneerr[1]={0.00};
    
        double fMx[1]={ 3.22E-08};
        double fMerr[1]={ 0.962E-08};
        double fMene[1]={4.2};
        double fMeneerr[1]={0.0};

    
	TGraphErrors *fE = new TGraphErrors(1,fEene,fEx,fEeneerr,fEerr);
    TGraphErrors *fM = new TGraphErrors(1,fMene,fMx,fEeneerr,fMerr);
    
/////////////////////////////////////RICK_FIRESTONE////////////
    
    double frickEx[1]={1.64E-07};
    double frickEerr[1]={0.741E-07};
    double frickEene[1]={3.79};
    double frickEeneerr[1]={0.00};
    
    double frickMx[1]={ 4.07E-08};
    double frickMerr[1]={ 1.21E-08};
    double frickMene[1]={3.87};
    double frickMeneerr[1]={0.0};
    
    
	TGraphErrors *frickE = new TGraphErrors(1,frickEene,frickEx,frickEeneerr,frickEerr);
    TGraphErrors *frickM = new TGraphErrors(1,frickMene,frickMx,frickEeneerr,frickMerr);

    
///////////////////////////////////////////////////////

    c1->SetLogy();
    c1->SetRightMargin(0.05);
	c1->SetLeftMargin(0.1);
	h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.1);
	h->GetXaxis()->SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitle("#gamma-ray strength function (MeV^{-3})");
	h->Draw();
	
	strengthexp->SetMarkerStyle(21);
	strengthexp->SetMarkerSize(0.7);
	strengthexp->Draw("P");
    
//  strengthexp2->SetMarkerStyle(22);
//	strengthexp2->SetMarkerSize(0.9);
//	strengthexp2->Draw("P");
    
    strengthteo->SetLineColor(2);
    strengthteo->SetLineWidth(2);
	strengthteo->Draw("L");
    
    strengthfM1->SetLineColor(2);
    strengthfM1->SetLineWidth(2);
    strengthfM1->SetLineStyle(7);
	strengthfM1->Draw("L");
    
    strengthfE2->SetLineColor(2);
    strengthfE2->SetLineWidth(2);
    strengthfE2->SetLineStyle(7);
	strengthfE2->Draw("L");
    
    fE->SetMarkerStyle(24);
	fE->SetMarkerSize(1.2);
//	fE->Draw("P");
    
    fM->SetMarkerStyle(26);
	fM->SetMarkerSize(1.2);
//	fM->Draw("P");
    
    frickE->SetMarkerStyle(25);
	frickE->SetMarkerSize(1.2);
//	frickE->Draw("P");
    
    frickM->SetMarkerStyle(26);
	frickM->SetMarkerSize(1.2);
//	frickM->Draw("P");
	
/*	gdrexp->SetMarkerStyle(23);
	gdrexp->SetMarkerSize(0.9);
	gdrexp->Draw("P");
    
    gdrexp1->SetMarkerStyle(32);
	gdrexp1->SetMarkerSize(0.9);
//	gdrexp1->Draw("P");
*/	
	TLegend *leg = new TLegend(0.55,0.55,0.85,0.65);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
	leg->AddEntry(strengthexp,"(a, a') ^{233}U, present exp.","P");
    leg->Draw();
//    leg->AddEntry(fE,"E1 from ARC (n,#gamma)","P");
//    leg->Draw();
//    leg->AddEntry(fM,"M1 from ARC (n,#gamma)","P");
//	leg->Draw();
//    leg->AddEntry(frickE,"E1 from thermal (n,#gamma)","P");
//    leg->Draw();
//    leg->AddEntry(frickM,"M1 from thermal (n,#gamma)","P");
//	leg->Draw();
//    leg->AddEntry(gdrexp,"(#gamma, x)  ^{237}Np, Berman (1986)","P");
//	leg->AddEntry(gdrexp1,"^{238}U(#gamma,x), Caldwell (1980)","P");
	leg->Draw();
    
    TLatex t;
	t.SetTextSize(0.03);
    t.DrawLatex(9.0,1.e-05,"GEDR");
	t.DrawLatex(1.5,7.e-08,"SR");
    t.DrawLatex(3.9,3.5e-08,"pygmy1");
	t.DrawLatex(11.,3.5e-08,"pygmy2");
    
    t.DrawLatex(4.2,4.1e-06,"S_{n}(^{ 238}Np)");
    TArrow *arrow1 = new TArrow(5.488,3.00E-06,5.488,1.70E-06,0.02,">");
    arrow1->SetLineWidth(2);
    arrow1->Draw();
	
	c1->Update();
	c1->Print("teoexp2_E2.pdf");
    c1->Print("teoexp2_E2.png");
    c1->Print("teoexp2_E2.eps");
}
