{
   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetFillColor(0);
   gStyle->SetPadBorderMode(0);
   m = (TH1F*)gROOT->FindObject("h");
   if (m) m->Delete();
   TCanvas *c1 = new TCanvas("c1","Normalization of gamma-transmission coefficient",600,600);
   TH2F *h = new TH2F("h"," ",10,-0.867078,   6.356,50,1.526e+01,3.246e+07);
   ifstream sigfile("sigpaw.cnt");
   float sig[62],sigerr[62];
   float energy[193],energyerr[193];
   float extL[194],extH[194];
   int i;
   float a0 = -0.8671;
   float a1 =  0.1037;
   for(i = 0; i < 70; i++){
   	energy[i] = a0 + (a1*i);
   	energyerr[i] = 0.0;
   	extL[i] = 0.0;
   	extH[i] = 0.0;
   }
   float x, y;
   i = 0;
   while(sigfile){
   	sigfile >> x;
   	if(i<61){
   		sig[i]=x;
   	}
   	else{sigerr[i-61]=x;}
   	i++;
   }
   ifstream extendfile("extendLH.cnt");
   i = 0;
   while(extendfile){
   	extendfile >> x >> y ;
   	extL[i]=x;
   	extH[i]=y;
   	i++;
   }
   TGraph *extLgraph = new TGraph(70,energy,extL);
   TGraph *extHgraph = new TGraph(70,energy,extH);
   TGraphErrors *sigexp = new TGraphErrors(61,energy,sig,energyerr,sigerr);
   c1->SetLogy();
   c1->SetLeftMargin(0.14);
   h->GetXaxis()->CenterTitle();
   h->GetXaxis()->SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
   h->GetYaxis()->CenterTitle();
   h->GetYaxis()->SetTitleOffset(1.4);
   h->GetYaxis()->SetTitle("Transmission coeff. (arb. units)");
   h->Draw();
   sigexp->SetMarkerStyle(21);
   sigexp->SetMarkerSize(0.8);
   sigexp->Draw("P");
   extLgraph->SetLineStyle(1);
   extLgraph->DrawGraph(23,&extLgraph->GetX()[0],&extLgraph->GetY()[0],"L");
   extHgraph->SetLineStyle(1);
   extHgraph->DrawGraph(26,&extHgraph->GetX()[44],&extHgraph->GetY()[44],"L");
   TArrow *arrow1 = new TArrow(1.104e+00,8.641e+03,1.104e+00,1.533e+03,0.02,">");
   arrow1->Draw();
   TArrow *arrow2 = new TArrow(1.415e+00,2.374e+04,1.415e+00,4.210e+03,0.02,">");
   arrow2->Draw();
   TArrow *arrow3 = new TArrow(3.697e+00,5.741e+05,3.697e+00,1.018e+05,0.02,">");
   arrow3->Draw();
   TArrow *arrow4 = new TArrow(4.941e+00,5.191e+06,4.941e+00,9.207e+05,0.02,">");
   arrow4->Draw();
   c1->Update();
   c1->Print("sigext.pdf");
   c1->Print("sigext.eps");
   c1->Print("sigext.ps");
}
