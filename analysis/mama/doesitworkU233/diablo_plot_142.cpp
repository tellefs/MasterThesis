// ROOT script to plot shape-method results with the OCL results,
// the case of 56Fe.
// Also a little chi^2 function from Dag Gillberg to get the scaling factor.
// From Magne's diablo.c, version 1.5, givin diablo_plot.cpp
// Modified script: Ann-Cecilie Larsen, a.c.larsen@fys.uio.no
// 16 June 2020
// Last update: 17 June 2020

void diablo_plot_142(){
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetFillColor(0);
    gStyle->SetPadBorderMode(0);

    // Read fg matrix
        // Pointer to input file, MAMA matrix 
    ifstream  ifile("fg");

    // Declarations of various help variables
    string line;
    string cal_dummy;
    string dim_dummy;
    int dim;
    int dim_start;
    int dim_stop;
    int dim_size;
    int position;
    int file_length;
    line.resize(200);   // need long enough line to read MAMA headers
    double x_cal[3] = {0.,1.,0.};   // calibration coeffs. on x axis: a0, a1, a2
    double y_cal[3] = {0.,1.,0.};   // calibration coeffs. on y axis: a0, a1, a2
    int dx, dy; // dimension on x and y axis
    int ix, iy;
    double value;
    double x,y;
    double number_of_counts = 0.;
    double new_y1, new_y2; 

    
    // read MAMA header (fixed format). The 10 first lines are info text 
    if(!getline(ifile,line) || line.substr(0,10) != "!FILE=Disk"){  // check correct format
        printf("\n This is not a MAMA file!!!\n ");
        exit(2);
    }   
    getline(ifile,line);    // skip !KIND=Spectrum
    getline(ifile,line);    // skip !LABORATORY=Oslo Cyclotron Laboratory (OCL) 
    getline(ifile,line);    // skip !EXPERIMENT=mama
    getline(ifile,line);    // skip !COMMENT=Sorted simulated data
    getline(ifile,line);    // skip !TIME=DATE:    19/11/09 11:47:26
    getline(ifile,line);    // get line with calibration
    cout << "\n Reading calibration coeffs.:" << endl;
    // calibration on x axis
    cal_dummy = line.substr(20,13); // position 20, length 13 characters
    if(!(istringstream(cal_dummy) >> x_cal[0])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a0 on x axis is: " << x_cal[0] << " keV." << endl;
    cal_dummy = line.substr(34,13); 
    if(!(istringstream(cal_dummy) >> x_cal[1])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a1 on x axis is: " << x_cal[1] << " keV/ch." << endl;
    cal_dummy = line.substr(48,13); 
    if(!(istringstream(cal_dummy) >> x_cal[2])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a2 on x axis is: " << x_cal[2] << " (keV/ch)^2." << endl;
    // calibration on y axis
    cal_dummy = line.substr(62,13); 
    if(!(istringstream(cal_dummy) >> y_cal[0])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a0 on y axis is: " << y_cal[0] << " keV." << endl;
    cal_dummy = line.substr(76,13); 
    if(!(istringstream(cal_dummy) >> y_cal[1])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a1 on y axis is: " << y_cal[1] << " keV/ch." << endl;
    cal_dummy = line.substr(90,13); 
    if(!(istringstream(cal_dummy) >> y_cal[2])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a2 on y axis is: " << y_cal[2] << " (keV/ch)^2.\n" << endl;
    getline(ifile,line);    // skip !PRECISION=16
    getline(ifile,line);    // get dimension
    // dimension of matrix
    dim_start = line.find_first_of("=") + 1;
    dim_dummy = line.substr(dim_start,1);
    if(!(istringstream(dim_dummy) >> dim)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension of matrix is: " << dim << endl;    
    getline(ifile,line);    // get channels
    // dimension on x axis
    dim_start = line.find_first_of(":") + 1;
    dim_stop = line.find_last_of(",");
    dim_size = dim_stop - dim_start;
    dim_dummy = line.substr(dim_start,dim_size);
    if(!(istringstream(dim_dummy) >> dx)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension on x axis is: " << dx+1 << " ch." << endl; 
    dx = dx+1;
    // dimension on y axis
    dim_start = line.find_last_of(":");
    dim_stop = line.find_last_of(")");
    dim_size = dim_stop - dim_start;
    dim_dummy = line.substr(dim_start+1,dim_size-1);
    if(!(istringstream(dim_dummy) >> dy)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension on y axis is: " << dy+1 << " ch." << endl; 
    dy = dy+1;

    
    // Make histogram 
    TH2D *matrix = new TH2D("matrix"," ",dx,(x_cal[0]/1000.),(dx*x_cal[1]+x_cal[0])/1000.,dy,(y_cal[0]/1000.),(dy*y_cal[1]+y_cal[0])/1000.);
    matrix->SetOption("colz");
    gStyle->SetPalette(1);
    
    for(iy=0;iy<dy;iy++){
        for(ix=0;ix<dx;ix++){
            ifile >> value;
            number_of_counts += value;
            matrix->SetBinContent(ix,iy,value);
        }
    }

    
    cout << " Matrix is now filled." << endl;
    cout << " Total number of counts in the matrix: " << number_of_counts << endl;
    cout << " Closing file."  << endl;
    // close file
    ifile.close();  

                   // Standard Oslo method begin
         ifstream strengthfilet("../../nd142/strength.nrm");
          float strength[90],strengtherr[90],egerr[90];
          float energy[731];
          int i = 0;
    float a0 =  -0.8588;
    float a1 =   0.1261;
    //      float x;
          while(strengthfilet){
              strengthfilet >> x;
              if(i<71){
                  strength[i] = 5.*x;
                  energy[i] = a0 + (a1*i);
                  if(energy[i]<=0.9)energy[i]=0.;
                  egerr[i] = 0.0;
              }
              else{strengtherr[i-71] = 5.*x;}
              i++;
          }
          TGraphErrors *strengthexpt = new TGraphErrors(72,energy,strength,egerr,strengtherr);
        //Standard Oslo method end

               // Standard Oslo method begin
     ifstream strengthfile("../../nd142r/strength.nrm");
        i = 0;
       while(strengthfile){
          strengthfile >> x;
          if(i<71){
              strength[i] = x;
              energy[i] = a0 + (a1*i);
              if(energy[i]<=0.9)energy[i]=0.;
              egerr[i] = 0.0;
          }
          else{strengtherr[i-71] = x;}
          i++;
      }
      TGraphErrors *strengthexp = new TGraphErrors(72,energy,strength,egerr,strengtherr);
    //Standard Oslo method end
    
    

    ifstream diafile("results.dia");
    float eg[500], energyerr[500];
    float gsf_ave[500],dgsf_ave[500],dsys_ave[500];
    float gsf1[500],gsf2[500],dgsf1[500],dgsf2[500];
    float sysL[500], sysH[500];
    float Eg, gsf, dgsf, dsys, gsf_1,gsf_2,dgsf_1,dgsf_2;
    int   idum;
    i = 0;
    float xnorm = 0.05*.8;
    while(diafile){
        diafile >>idum >> Eg >>gsf>>dgsf>>dsys>>gsf_1>>gsf_2>>dgsf_1>>dgsf_2;
        eg[i]        = Eg/1000.;
        if(Eg<=2200.)eg[i]=0.;
        gsf_ave[i]   = gsf*xnorm;
        dgsf_ave[i]  = dgsf*xnorm;
        sysL[i]      = (gsf-dsys)*xnorm;
        sysH[i]      = (gsf+dsys)*xnorm;
        gsf1[i]      = gsf_1*xnorm;
        gsf2[i]      = gsf_2*xnorm;
        dgsf1[i]     = dgsf_1*xnorm;
        dgsf2[i]     = dgsf_2*xnorm;
        energyerr[i] = 0.;
        i++;
    }
    TGraphErrors *gsfgr = new TGraphErrors(i+1,eg,gsf_ave,energyerr,dgsf_ave);
    TGraph       *sysLgr = new TGraphErrors(i,eg,sysL);
    TGraph       *sysHgr = new TGraphErrors(i,eg,sysH);
    TGraphErrors *gsf1gr = new TGraphErrors(i,eg,gsf1,energyerr,dgsf1);
    TGraphErrors *gsf2gr = new TGraphErrors(i,eg,gsf2,energyerr,dgsf2);
   
    /////****************************************************////
       ////Begin finding best factor and Chi**2 to match gsf////
       int ii_max = i;
    float EgFitLow = 4.0, EgFitHigh = 7.0;
       float shape[500], dshape[500];
       int ii=0,ii1=0,ii2=0;
       //Making shape-data the same calibration as gsf data //
       //by interpolating shape-data//
       for(i=0;i<73;i++){
           if(energy[i]>EgFitLow && energy[i]<EgFitHigh){
               for(ii=0; ii<ii_max;ii++){
                   if(eg[ii]<energy[i]){
                       ii1=ii;ii2=ii1+1;
                       if(ii1>ii_max){
                           ii1=ii-1;ii2=ii1;
                       }
                   }
               }
               shape[i]=gsf_ave[ii1]+((gsf_ave[ii2]-gsf_ave[ii1])/(eg[ii2]-eg[ii1]))*(energy[i]-eg[ii1]);
               dshape[i]=dgsf_ave[ii1]+((dgsf_ave[ii2]-dgsf_ave[ii1])/(eg[ii2]-eg[ii1]))*(energy[i]-eg[ii1]);
           }
       }
    //First estimate: Matching in the midle channel
       float EgMiddle = (EgFitLow + EgFitHigh)/2.0;
       int   ChMiddle = ((EgMiddle-a0)/a1 + 0.5);
       float factor, factor_basis= strength[ChMiddle]/shape[ChMiddle];
       float factor_best=1.;
       float alpha, alpha_min=0.3, alpha_max = 1.7, d_alpha = 0.01;
       double chix=0., chi2=0., chi2_old = 1000000000.,chi_best;
       int nch=0.;
       alpha = alpha_min;
       while (alpha < alpha_max){
           factor = alpha*factor_basis;
           nch =0, chix=0.;
           for(i=0;i<73;i++){
                  if(energy[i]>EgFitLow && energy[i]<EgFitHigh){
                      chix = chix + (factor*shape[i]-strength[i])*(factor*shape[i]-strength[i])/(factor*factor*dshape[i]*dshape[i]+(strengtherr[i])*(strengtherr[i]));
                      nch = nch + 1;
                  }
              }
           chi2 = chix;
           printf("alpha= %6.4f factor = %6.4f chi2= %e\n", alpha, factor, chi2);
           if (chi2<chi2_old){
               factor_best = factor;
               chi_best = chi2;
           }
           chi2_old = chi2;
           alpha = alpha + d_alpha;
       }
        printf(" Best normalization of shape to standard Oslo is:\n Factor= %6.4f found in fit-region (EgL,EgH) = (%6.3f,%6.3f)MeV\n Number of data points nch = %3d, Reduced Chi2 = %8.3e and total Chi2 = %8.3e \n" ,factor_best, EgFitLow, EgFitHigh, nch, chi_best/(nch-1), chi_best);
       ////End finding best factor and Chi**2 to match gsf, now making plots for normalized shape////
       /////****************************************************////

       for(i=0;i<73;i++){
           gsf1[i]  = factor_best*gsf1[i];
           dgsf1[i] = factor_best*dgsf1[i];
           gsf2[i]  = factor_best*gsf2[i];
           dgsf2[i] = factor_best*dgsf2[i];
       }
       TGraphErrors *gsf1grn = new TGraphErrors(ii_max,eg,gsf1,energyerr,dgsf1);
       TGraphErrors *gsf2grn = new TGraphErrors(ii_max,eg,gsf2,energyerr,dgsf2);
    
    TCanvas *c1 = new TCanvas("c1","FG matrix and Gamma-ray strength function",1250,600);
    c1->Divide(2,1,0,0);
    TH2F *h     = new TH2F("h"," ",100,.5, 8.1,100,1.e-9,1.e-06);

    c1->cd(1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);

    matrix->GetXaxis()->SetTitle("#gamma-ray energy #it{E}_{#gamma} (MeV)");
    matrix->GetXaxis()->CenterTitle();
    matrix->GetXaxis()->SetTitleFont(42);
    matrix->GetXaxis()->SetTitleSize(0.04);
    matrix->GetXaxis()->SetLabelFont(42);
    matrix->GetXaxis()->SetTitleOffset(1.3);
    matrix->GetYaxis()->SetTitle("Excitation energy #it{E_{i}} (MeV)");
    matrix->GetZaxis()->SetLabelFont(42);
    matrix->GetZaxis()->SetTitleFont(42);
    matrix->GetZaxis()->SetTitleOffset(1.4);
    matrix->GetZaxis()->SetTitle("Counts ");
    matrix->GetYaxis()->SetTitleOffset(1.1);
    matrix->GetYaxis()->SetTitleSize(0.04);
    matrix->GetYaxis()->CenterTitle();
    matrix->GetYaxis()->SetTitleFont(42);
    matrix->GetYaxis()->SetLabelFont(42);
    // Set user range in the matrix if you want to
    matrix->GetXaxis()->SetRangeUser(0.5,8.500);
    matrix->GetYaxis()->SetRangeUser(0.5,8.200);
    matrix->GetZaxis()->SetRangeUser(1.,2000);
    //matrix->SetContour(90); // number of divisions in the color range
    // Nice draw options: surf1, colz, contz, surf2z. The FB command suppresses the front-end axis.
    // If you want to draw in black&white, use surf or surf4, but then you might need to rebin a lot
    matrix->Draw("colz ");
    
    gPad->Update();
    
    // Play with the palette to make it look nice,
    // only works when the "z" option is in the draw command, e.g. surf2z
    TPaletteAxis *palette = (TPaletteAxis*)matrix->GetListOfFunctions()->FindObject("palette");
    if(!palette) cout << "nono." << endl;
    palette->SetX1NDC(0.901);
    palette->SetX2NDC(0.93);
    palette->SetY1NDC(0.12);
    palette->SetY2NDC(0.92);
    palette->SetLabelSize(0.04);
    palette->SetLabelFont(42);
    palette->Draw();        
    
    
    ///////////
    float x1L1 , y1L1, x2L1, y2L1;
    float x1U1 , y1U1, x2U1, y2U1;
    float x1L2 , y1L2, x2L2, y2L2;
    float x1U2 , y1U2, x2U2, y2U2;
    ifstream infile("input.dia");
    infile >>x1L1>>y1L1>>x2L1>>y2L1;
    infile >>x1U1>>y1U1>>x2U1>>y2U1;
    infile >>x1L2>>y1L2>>x2L2>>y2L2;
    infile >>x1U2>>y1U2>>x2U2>>y2U2;
    infile.close();
    TLine *line1 = new TLine(x1L1/1000., y1L1/1000., x2L1/1000., y2L1/1000.);
    line1->SetLineStyle(1);
    line1->SetLineWidth(1);
    line1->Draw();
    TLine *line2 = new TLine(x1U1/1000., y1U1/1000., x2U1/1000., y2U1/1000.);
    line2->SetLineStyle(1);
    line2->SetLineWidth(1);
    line2->Draw();
    TLine *line3 = new TLine(x1L2/1000., y1L2/1000., x2L2/1000., y2L2/1000.);
    line3->SetLineStyle(1);
    line3->Draw();
    TLine *line4 = new TLine(x1U2/1000., y1U2/1000., x2U2/1000., y2U2/1000.);
    line4->SetLineStyle(1);
    line4->Draw();
    /////////////

    TLatex text;
    text.SetTextSize(0.04);
    text.SetTextFont(42);
    text.DrawLatex(6.40,8.4,"  D_{2}             D_{1}");
    text.DrawLatex(4.5,1.5,"");

    text.SetTextSize(0.045);
    text.DrawLatex(   1.0 ,8.4,"(a)^{142}Nd(p,p' #gamma)^{142}Nd");

    c1->cd(2);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTitle("#gamma-ray energy #it{E}_{#gamma} (MeV)");
    h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetYaxis()->SetTitle("#gamma-ray strength function #it{f}(#it{E}_{#gamma}) (MeV^{-3})");
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->Draw();

    strengthexp->SetMarkerStyle(21);
    strengthexp->SetMarkerSize(0.8);
    strengthexp->SetMarkerColor(kBlack);
    strengthexp->Draw("P");
    
    strengthexpt->SetMarkerStyle(21);
    strengthexpt->SetMarkerSize(0.8);
    strengthexpt->SetMarkerColor(kGreen);
    strengthexpt->Draw("P");

    gsf1gr->SetMarkerStyle(22);
    gsf1gr->SetMarkerSize(0.9);
    gsf1gr->SetMarkerColor(kBlue);
    gsf1gr->SetLineColor(kBlue);
//    gsf1gr->Draw("P");
    
    gsf2gr->SetMarkerStyle(26);
    gsf2gr->SetMarkerSize(0.9);
    gsf2gr->SetMarkerColor(kBlue);
    gsf2gr->SetLineColor(kBlue);
 //   gsf2gr->Draw("P");
    
    gsf1grn->SetMarkerStyle(22);
        gsf1grn->SetMarkerSize(0.9);
        gsf1grn->SetMarkerColor(kBlue);
        gsf1grn->SetLineColor(kBlue);
        gsf1grn->Draw("P");
        
        gsf2grn->SetMarkerStyle(26);
        gsf2grn->SetMarkerSize(0.9);
        gsf2grn->SetMarkerColor(kBlue);
        gsf2grn->SetLineColor(kBlue);
        gsf2grn->Draw("P");

    gsfgr->SetMarkerStyle(20);
    gsfgr->SetMarkerSize(0.6);
    gsfgr->SetMarkerColor(kAzure+1);
  //  gsfgr->Draw("P");

    sysLgr->SetMarkerStyle(20);
    sysLgr->SetMarkerSize(0.6);
    sysLgr->SetMarkerColor(kBlack);
  //  sysLgr->Draw("L");
    sysHgr->SetMarkerStyle(20);
    sysHgr->SetMarkerSize(0.6);
    sysHgr->SetMarkerColor(kBlack);
  //  sysHgr->Draw("L");

    TLegend *leg = new TLegend(0.14,0.75,0.6,0.9);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      
      leg->AddEntry(strengthexpt,"Oslo method #rho(Sn)","P");
    leg->AddEntry(strengthexp,"Oslo method 1 #times #rho(Sn)","P");

      leg->AddEntry(gsf1gr,"#gamma-rays feeding D1","P");
      leg->AddEntry(gsf2gr,"#gamma-rays feeding D2","P");
//      leg->AddEntry(sysHgr,"Systematical errors","L");
      leg->Draw();
    
    TLatex t;
    t.SetTextSize(0.047);
    t.SetTextFont(42);
    t.DrawLatex(    1. ,1.3e-06,"(b)");

    c1->Update();
    c1->Print("diablo_plot_142nd.pdf");
}
