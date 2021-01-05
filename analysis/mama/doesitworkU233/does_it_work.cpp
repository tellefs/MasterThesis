//////////////////////////////////////////////////////
//	Script to check the quality of the extracted 	//
//	level density and strength function				//
//	by plotting exp. and calculated f.g. spectra 	//
//  To be run in the folder where you have fg.rsg,  //
//  fgerr.rsg and fgteo.rsg, or give                //
//  appropriate paths.                              //
//	Cecilie, June 7, 2012. Update 5 Aug 2019        //
//////////////////////////////////////////////////////
{
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(1);
	
	// declarations
	ifstream  ifile1("fg_compressed.rsg"), ifile2("fgerr_compressed.rsg"), ifile3("fgteo_compressed.rsg");
	string line;
	string cal_dummy;
	string dim_dummy;
	char pdf_filename[512];
	int dim;
	int dim_start;
	int dim_stop;
	int dim_size;
	int position;
	int file_length;
	line.resize(200);	// need long enough line to read MAMA headers
	double x_cal[3] = {0.,1.,0.};	// calibration coeffs. on x axis: a0, a1, a2
	double y_cal[3] = {0.,1.,0.};	// calibration coeffs. on y axis: a0, a1, a2
	int dx, dy;	// dimension on x and y axis
	int ix, iy;
	double value;
	double x,y;
	double number_of_counts = 0.;
	double new_y1, new_y2; 
	int sign_ycal = 0;

	
	// open file to read
	if(!ifile1){
		cout << "\n Could not open file!!!\n ";
		exit(1);
	}
	else cout << "\n Successful opening of file"  << endl;
	
	
	// read MAMA header (fixed format). The 10 first lines are info text 
	if(!getline(ifile1,line) || line.substr(0,10) != "!FILE=Disk"){	// check correct format
		printf("\n This is not a MAMA file!!!\n ");
		exit(2);
	}	
	getline(ifile1,line);	// skip !KIND=Spectrum
	getline(ifile1,line);	// skip !LABORATORY=Oslo Cyclotron Laboratory (OCL)	
	getline(ifile1,line);	// skip !EXPERIMENT=mama
	getline(ifile1,line);	// skip !COMMENT=Sorted simulated data
	getline(ifile1,line);	// skip !TIME=DATE:    19/11/09 11:47:26
	getline(ifile1,line);	// get line with calibration
	cout << "\n Reading calibration coeffs.:" << endl;
	// calibration on x axis
	cal_dummy = line.substr(20,13);	// position 20, length 13 characters
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
	getline(ifile1,line);	// skip !PRECISION=16
	getline(ifile1,line);	// get dimension
	// dimension of matrix
	dim_start = line.find_first_of("=") + 1;
	dim_dummy = line.substr(dim_start,1);
	if(!(istringstream(dim_dummy) >> dim)) cout << "Could not convert string to number." << endl;
	else cout << " Dimension of matrix is: " << dim << endl;	
	getline(ifile1,line);	// get channels
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
	dy = dy+2;

	// Test if negative calibration coeff. on Ex, then invert axis:
	if(y_cal[1] < 0.){
		sign_ycal = -1;
		new_y1 = y_cal[0] + (y_cal[1]*(double) dy);
		new_y2 = y_cal[0];
		
		y_cal[0] = new_y1;
		y_cal[1] = (-1.)*y_cal[1];
		cout << " New calibration on y axis: y_cal[0] = " << y_cal[0] << ", y_cal[1] = " << y_cal[1] << endl;

	}
	
	// Make histogram 
	TH2D *matrix = new TH2D("matrix"," ",dx,x_cal[0],dx*x_cal[1]+x_cal[0],dy,y_cal[0],dy*y_cal[1]+y_cal[0]);
	matrix->SetOption("colz");
	gStyle->SetPalette(1);
    
    // Make vectors for the spectra we want to show
    Double_t *energy = new Double_t[dx];
    Double_t *fg1    = new Double_t[dx]; 
    Double_t *fg2    = new Double_t[dx];  
    Double_t *fg3    = new Double_t[dx]; 
    Double_t *fg4    = new Double_t[dx]; 
    Double_t *fg5    = new Double_t[dx]; 
    Double_t *fg6    = new Double_t[dx]; 

    Double_t *fg1err = new Double_t[dx];
    Double_t *fg2err = new Double_t[dx];
    Double_t *fg3err = new Double_t[dx];
    Double_t *fg4err = new Double_t[dx];
    Double_t *fg5err = new Double_t[dx];
    Double_t *fg6err = new Double_t[dx];

    //choose excitation-energy bins you want to show
    int bin1 = 13; //50
    int bin2 = 14; //53
    int bin3 = 15; //56
    int bin4 = 16; //59
    int bin5 = 17; //62
    int bin6 = 18; //65
	
	if(sign_ycal < 0.){	// if negative calibration coeff. on y axis
		for(iy=dy;iy>0;iy--){
			for(ix=0;ix<dx;ix++){
				ifile1 >> value;
				number_of_counts += value;
				matrix->SetBinContent(ix,iy,value);
                if(iy==bin1) fg1[ix] = value;
                if(iy==bin2) fg2[ix] = value;
                if(iy==bin3) fg3[ix] = value;
                if(iy==bin4) fg4[ix] = value;
                if(iy==bin5) fg5[ix] = value;
                if(iy==bin6) fg6[ix] = value;
			}
            energy[iy] = iy*x_cal[1]+x_cal[0];
		}
	}
	else{	// if positive calibration coeff. on y axis
		for(iy=0;iy<dy;iy++){
			for(ix=0;ix<dx;ix++){
				ifile1 >> value;
				number_of_counts += value;
				matrix->SetBinContent(ix,iy,value);
                if(iy==bin1) fg1[ix] = value;
                if(iy==bin2) fg2[ix] = value;
                if(iy==bin3) fg3[ix] = value;
                if(iy==bin4) fg4[ix] = value;
                if(iy==bin5) fg5[ix] = value;
                if(iy==bin6) fg6[ix] = value;
			}
            energy[iy] = iy*x_cal[1]+x_cal[0] - (x_cal[1]/2.); 
            // the last term in the expression above is to get the bin center in the graphs
            // the same as the bin center in the histograms
            //cout << "Bin y: " << iy << ", Eg (or Ex): " << energy[iy] << endl;
            //cout << iy << " " << energy[iy] << endl;
            
		}
	}
	
	cout << " Matrix is now filled." << endl;
        
    // open file to read
    if(!ifile2){
        cout << "\n Could not open file!!!\n ";
        exit(1);
    }
    else cout << "\n Successful opening of file"  << endl;
        
        
    // read MAMA header (fixed format). The 10 first lines are info text 
    if(!getline(ifile2,line) || line.substr(0,10) != "!FILE=Disk"){	// check correct format
        printf("\n This is not a MAMA file!!!\n ");
        exit(2);
    }	
    getline(ifile2,line);	// skip !KIND=Spectrum
    getline(ifile2,line);	// skip !LABORATORY=Oslo Cyclotron Laboratory (OCL)	
    getline(ifile2,line);	// skip !EXPERIMENT=mama
    getline(ifile2,line);	// skip !COMMENT=Sorted simulated data
    getline(ifile2,line);	// skip !TIME=DATE:    19/11/09 11:47:26
    getline(ifile2,line);	// get line with calibration
    cout << "\n Reading calibration coeffs.:" << endl;
    // calibration on x axis
    cal_dummy = line.substr(20,13);	// position 20, length 13 characters
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
    getline(ifile2,line);	// skip !PRECISION=16
    getline(ifile2,line);	// get dimension
    // dimension of matrix
    dim_start = line.find_first_of("=") + 1;
    dim_dummy = line.substr(dim_start,1);
    if(!(istringstream(dim_dummy) >> dim)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension of matrix is: " << dim << endl;	
    getline(ifile2,line);	// get channels
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
        
    // Test if negative calibration coeff. on Ex, then invert axis:
    if(y_cal[1] < 0.){
        sign_ycal = -1;
        new_y1 = y_cal[0] + (y_cal[1]*(double) dy);
        new_y2 = y_cal[0];
            
        y_cal[0] = new_y1;
        y_cal[1] = (-1.)*y_cal[1];
        cout << " New calibration on y axis: y_cal[0] = " << y_cal[0] << ", y_cal[1] = " << y_cal[1] << endl;
            
    }
        
    // Make histogram 
    TH2D *matrix2 = new TH2D("matrix2"," ",dx,x_cal[0],dx*x_cal[1]+x_cal[0],dy,y_cal[0],dy*y_cal[1]+y_cal[0]);
    matrix2->SetOption("colz");
    gStyle->SetPalette(1);
        
    if(sign_ycal < 0.){	// if negative calibration coeff. on y axis
        for(iy=dy;iy>0;iy--){
            for(ix=0;ix<dx;ix++){
                ifile2 >> value;
                number_of_counts += value;
                matrix2->SetBinContent(ix,iy,value);
                if(iy==bin1) fg1err[ix] = value;
                if(iy==bin2) fg2err[ix] = value;
                if(iy==bin3) fg3err[ix] = value;
                if(iy==bin4) fg4err[ix] = value;
                if(iy==bin5) fg5err[ix] = value;
                if(iy==bin6) fg6err[ix] = value;
            }
        }
    }
    else{	// if positive calibration coeff. on y axis
        for(iy=0;iy<dy;iy++){
            for(ix=0;ix<dx;ix++){
                ifile2 >> value;
                number_of_counts += value;
                matrix2->SetBinContent(ix,iy,value);
                if(iy==bin1) fg1err[ix] = value;
                if(iy==bin2) fg2err[ix] = value;
                if(iy==bin3) fg3err[ix] = value;
                if(iy==bin4) fg4err[ix] = value;
                if(iy==bin5) fg5err[ix] = value;
                if(iy==bin6) fg6err[ix] = value;
            }
        }
    }
        
    cout << " Matrix 2 is now filled." << endl;

    // open file to read
    if(!ifile3){
        cout << "\n Could not open file!!!\n ";
        exit(1);
    }
    else cout << "\n Successful opening of file"  << endl;
    
    
    // read MAMA header (fixed format). The 10 first lines are info text 
    if(!getline(ifile3,line) || line.substr(0,10) != "!FILE=Disk"){	// check correct format
        printf("\n This is not a MAMA file!!!\n ");
        exit(2);
    }	
    getline(ifile3,line);	// skip !KIND=Spectrum
    getline(ifile3,line);	// skip !LABORATORY=Oslo Cyclotron Laboratory (OCL)	
    getline(ifile3,line);	// skip !EXPERIMENT=mama
    getline(ifile3,line);	// skip !COMMENT=Sorted simulated data
    getline(ifile3,line);	// skip !TIME=DATE:    19/11/09 11:47:26
    getline(ifile3,line);	// get line with calibration
    cout << "\n Reading calibration coeffs.:" << endl;
    // calibration on x axis
    cal_dummy = line.substr(20,13);	// position 20, length 13 characters
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
    getline(ifile3,line);	// skip !PRECISION=16
    getline(ifile3,line);	// get dimension
    // dimension of matrix
    dim_start = line.find_first_of("=") + 1;
    dim_dummy = line.substr(dim_start,1);
    if(!(istringstream(dim_dummy) >> dim)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension of matrix is: " << dim << endl;	
    getline(ifile3,line);	// get channels
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
    
    // Test if negative calibration coeff. on Ex, then invert axis:
    if(y_cal[1] < 0.){
        sign_ycal = -1;
        new_y1 = y_cal[0] + (y_cal[1]*(double) dy);
        new_y2 = y_cal[0];
        
        y_cal[0] = new_y1;
        y_cal[1] = (-1.)*y_cal[1];
        cout << " New calibration on y axis: y_cal[0] = " << y_cal[0] << ", y_cal[1] = " << y_cal[1] << endl;        
    }
        
    // Make histogram 
    TH2D *matrix3 = new TH2D("matrix3"," ",dx,x_cal[0],dx*x_cal[1]+x_cal[0],dy,y_cal[0],dy*y_cal[1]+y_cal[0]);
    matrix3->SetOption("colz");
    gStyle->SetPalette(1);
    
    if(sign_ycal < 0.){	// if negative calibration coeff. on y axis
        for(iy=dy;iy>0;iy--){
            for(ix=0;ix<dx;ix++){
                ifile3 >> value;
                number_of_counts += value;
                matrix3->SetBinContent(ix,iy,value);
            }
        }
    }
    else{	// if positive calibration coeff. on y axis
        for(iy=0;iy<dy;iy++){
            for(ix=0;ix<dx;ix++){
                ifile3 >> value;
                number_of_counts += value;
                matrix3->SetBinContent(ix,iy,value);                    
            }
        }
    }
    
    cout << " Matrix 3 is now filled." << endl;

    // close files
	ifile1.close();	
    ifile2.close();	
    ifile3.close();	

    // Projections for fgteo
    TH1D *px1 = matrix3->ProjectionX("px1",bin1,bin1);
    TH1D *px2 = matrix3->ProjectionX("px2",bin2,bin2);
    TH1D *px3 = matrix3->ProjectionX("px3",bin3,bin3);
    TH1D *px4 = matrix3->ProjectionX("px4",bin4,bin4);
    TH1D *px5 = matrix3->ProjectionX("px5",bin5,bin5);
    TH1D *px6 = matrix3->ProjectionX("px6",bin6,bin6);
        
    // Print Ex values to screen:
    cout << " Ex, bin 1 = " << matrix->GetYaxis()->GetBinCenter(bin1) << " keV" << endl;
    cout << " Ex, bin 2 = " << matrix->GetYaxis()->GetBinCenter(bin2) << " keV" << endl;
    cout << " Ex, bin 3 = " << matrix->GetYaxis()->GetBinCenter(bin3) << " keV" << endl;
    cout << " Ex, bin 4 = " << matrix->GetYaxis()->GetBinCenter(bin4) << " keV" << endl;
    cout << " Ex, bin 5 = " << matrix->GetYaxis()->GetBinCenter(bin5) << " keV" << endl;
    cout << " Ex, bin 6 = " << matrix->GetYaxis()->GetBinCenter(bin6) << " keV" << endl;

	// Create TCanvas
	TCanvas *c1 = new TCanvas("c1","Primary gamma spectra",1000,500);
    
    TH2D *h = new TH2D("h"," ",100,800,6400,100,0.001,0.079);
    TH2D *h_1 = new TH2D("h_1"," ",100,800,6400,100,0.001,0.079);
    TH2D *h2 = new TH2D("h2"," ",100,800,6400,100,0.001,0.079);
    TH2D *h_2 = new TH2D("h_2"," ",100,800,6400,100,0.001,0.07);
    TH2D *h_3 = new TH2D("h_3"," ",100,800,6400,100,0.001,0.079);

    double energy_errors[100] = {0.};
    TGraphErrors *graph1 = new TGraphErrors(dx,energy,fg1,energy_errors,fg1err);        
    TGraphErrors *graph2 = new TGraphErrors(dx,energy,fg2,energy_errors,fg2err);
    TGraphErrors *graph3 = new TGraphErrors(dx,energy,fg3,energy_errors,fg3err);
    TGraphErrors *graph4 = new TGraphErrors(dx,energy,fg4,energy_errors,fg4err);
    TGraphErrors *graph5 = new TGraphErrors(dx,energy,fg5,energy_errors,fg5err);
    TGraphErrors *graph6 = new TGraphErrors(dx,energy,fg6,energy_errors,fg6err);
       
    c1->Divide(3,2,0,0);
        
    c1->cd(1);
	gPad->SetLeftMargin(0.14);  
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.01);
    h->GetYaxis()->SetTitle("Probability distribution");
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetYaxis()->SetLabelSize(0.065);
    h->GetYaxis()->SetLabelFont(42);
    h->Draw();
	px1->SetLineColor(kAzure+8);
    px1->SetLineWidth(3);
    px1->Draw("same");
    graph1->SetMarkerStyle(2);
    graph1->Draw("P");
	
	TLatex text;
	text.SetTextSize(0.08);
    text.SetTextFont(42);
    text.DrawLatex(2200,0.07,"(a) E = 3.28 MeV");
	
    TLegend *leg = new TLegend(0.6,0.5,0.89,0.68);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.08);
    leg->AddEntry(graph1," first-gen. data ","P");
    leg->AddEntry(px1," #rho x T ","L");
    leg->Draw();
        
    c1->cd(2);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.01);
    matrix3->Draw();
    h_1->GetYaxis()->SetLabelOffset(100);
    h_1->Draw();
    px2->SetLineColor(kAzure+8);
    px2->SetLineWidth(3);
    px2->Draw("same");
    graph2->SetMarkerStyle(2);
    graph2->Draw("P same");

    text.DrawLatex(2200,0.07,"(b) E = 3.59 MeV");

    c1->cd(3);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.01);
    h_1->Draw();
    px3->SetLineColor(kAzure+8);
    px3->SetLineWidth(3);
    px3->Draw("same");
    graph3->SetMarkerStyle(2);
    graph3->Draw("P same");

        
    text.DrawLatex(2200,0.07,"(c) E = 3.9 MeV");
    
    c1->cd(4);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.18);
    //h2->GetXaxis()->SetTitle("#gamma energy E_{#gamma} (keV)");
    h2->GetXaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->CenterTitle();
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetTitleSize(0.065);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetLabelSize(0.065);
    h2->GetYaxis()->SetTitle("Probability distribution");
    h2->GetYaxis()->SetTitleOffset(0.9);
    h2->GetYaxis()->CenterTitle();
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetYaxis()->SetTitleSize(0.073);
    h2->GetYaxis()->SetLabelSize(0.058);
    h2->GetYaxis()->SetLabelFont(42);
    h2->Draw();
    px4->SetLineColor(kAzure+8);
    px4->SetLineWidth(3);
    px4->Draw("same");
    graph4->SetMarkerStyle(2);
    graph4->Draw("P same");
        
    text.SetTextSize(0.08);
    text.DrawLatex(2200,0.07,"(d) E = 5.0 MeV");
    
    
    gPad->Update();
        
    c1->cd(5);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.18);
    h_2->GetXaxis()->SetTitle("#gamma energy E_{#gamma} (keV)");
    h_2->GetXaxis()->SetTitleOffset(1.1);
    h_2->GetXaxis()->CenterTitle();
    h_2->GetXaxis()->SetTitleFont(42);
    h_2->GetXaxis()->SetTitleSize(0.069);
    h_2->GetXaxis()->SetLabelFont(42);
    h_2->GetXaxis()->SetLabelSize(0.065);
    h_2->GetYaxis()->SetTitle("Probability distribution");
    h_2->GetYaxis()->SetTitleOffset(0.9);
    h_2->GetYaxis()->CenterTitle();
    h_2->GetYaxis()->SetTitleFont(42);
    h_2->GetYaxis()->SetTitleSize(0.073);
    h_2->GetYaxis()->SetLabelOffset(100.058);
    h_2->GetYaxis()->SetLabelFont(42);
    h_2->Draw();
    px5->SetLineColor(kAzure+8);
    px5->SetLineWidth(3);
    px5->Draw("same");
    graph5->SetMarkerStyle(2);
    graph5->Draw("P same");
    text.DrawLatex(2200,0.07,"(e) E = 5.3 MeV");

 
    c1->cd(6);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.18);
    h_3->GetXaxis()->SetLabelFont(42);
    h_3->GetXaxis()->SetLabelSize(0.065);
    h_3->GetYaxis()->SetLabelOffset(42);
    h_3->Draw();
    px6->SetLineColor(kAzure+8);
    px6->SetLineWidth(3);
    px6->Draw("same");
    graph6->SetMarkerStyle(2);
    graph6->Draw("P same");
    text.DrawLatex(2200,0.07,"(f) E = 5.7 MeV");


    c1->Update();
	
	// Print to pdf file
	c1->Print("does_it_work.pdf");
    c1->Print("does_it_work.eps");

    delete [] energy;
    delete [] fg1;
    delete [] fg2;
    delete [] fg3;
    delete [] fg4;
    delete [] fg5;
    delete [] fg6;
    delete [] fg1err;
    delete [] fg2err;
    delete [] fg3err;
    delete [] fg4err;
    delete [] fg5err;
    delete [] fg6err;

		
}	// END of script
