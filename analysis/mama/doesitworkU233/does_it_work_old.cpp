//////////////////////////////////////////////////////
//	Script to check the quality of the extracted 	//
//	level density and strength function				//
//	by plotting exp. and calculated f.g. spectra 	//
//	written by: ACL. June 7, 2012					//
//////////////////////////////////////////////////////
void does_it_work(){
	
	// starting root stuff
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetFillColor(0);
    gStyle->SetPadBorderMode(0);
	
    /*
	m = (TH2F*)gROOT->FindObject("matrix");
	if (m) m->Delete();
    m2 = (TH2F*)gROOT->FindObject("matrix2");
    if (m2) m2->Delete();
    m3 = (TH2F*)gROOT->FindObject("matrix3");
    if (m3) m3->Delete();
    m4 = (TH2F*)gROOT->FindObject("h");
    if (m4) m4->Delete();
    m5 = (TH2F*)gROOT->FindObject("h2");
    if (m5) m5->Delete();
    m6 = (TH2F*)gROOT->FindObject("h_1");
    if (m6) m6->Delete();
    m7 = (TH2F*)gROOT->FindObject("h_2");
    if (m7) m7->Delete();
    m8 = (TH2F*)gROOT->FindObject("h_3");
    if (m8) m8->Delete();
    
    px1 = (TH1D*)gROOT->FindObject("px1");
    if (px1) px1->Delete();
    px2 = (TH1D*)gROOT->FindObject("px2");
    if (px2) px2->Delete();
    px3 = (TH1D*)gROOT->FindObject("px3");
    if (px3) px3->Delete();
    px4 = (TH1D*)gROOT->FindObject("px4");
    if (px4) px4->Delete();
    px5 = (TH1D*)gROOT->FindObject("px5");
    if (px5) px5->Delete();
    px6 = (TH1D*)gROOT->FindObject("px6");
    if (px6) px6->Delete();            
    */
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
	int sign_ycal;

	
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
        
        // Make MeV instead of keV on the x axis (E_gamma)
        x_cal[0] /= 1000.;
        x_cal[1] /= 1000.;

	
	// Make histogram 
	TH2D *matrix = new TH2D("matrix"," ",dx,x_cal[0],dx*x_cal[1]+x_cal[0],dy,y_cal[0],dy*y_cal[1]+y_cal[0]);
	matrix->SetOption("colz");
	gStyle->SetPalette(1);
    
    // 
    double energy[dx];
    double fg1[dx], fg2[dx], fg3[dx], fg4[dx], fg5[dx], fg6[dx];
    double fg1err[dx], fg2err[dx], fg3err[dx], fg4err[dx], fg5err[dx], fg6err[dx];
    double fg1teo[dx];
    //choose excitation-energy bins
    int bin1 = 13; //27
    int bin2 = 14; //31
    int bin3 = 15; //35
    int bin4 = 16; //39
    int bin5 = 17; //43
    int bin6 = 18; //47
	
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
  
        // Make MeV instead of keV on the x axis (E_gamma)
        x_cal[0] /= 1000.;
        x_cal[1] /= 1000.;

        
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
        
        // Make MeV instead of keV on the x axis (E_gamma)
        x_cal[0] /= 1000.;
        x_cal[1] /= 1000.;

 
        
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
                    if(iy==bin1) {
                        fg1teo[ix] = value;
                        //cout << "Matrix value: " << matrix3->GetBinContent(ix,iy,value) << endl;
                        //cout << "Graph value: " << fg1teo[ix] << endl;
                    }
                    

                }
            }
        }
        
        cout << " Matrix 3 is now filled." << endl;

    // close file
	ifile1.close();	
    ifile2.close();	
    ifile3.close();	

    // Projections for fgteo
    //Changed from (for all px): matrix3->ProjectionX("px1",bin1,bin1);
    //to: 
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
    
    TH2D *h = new TH2D("h"," ",100,0,7.5,100,0.001,0.179);
    TH2D *h_1 = new TH2D("h_1"," ",100,0,7.5,100,0.001,0.179);
    TH2D *h2 = new TH2D("h2"," ",100,0,7.5,100,0.001,0.179);
    TH2D *h_2 = new TH2D("h_2"," ",100,0,7.5,100,0.001,0.179);
    TH2D *h_3 = new TH2D("h_3"," ",100,0,7.5,100,0.001,0.179);

    //Changed from this:
    //double energy_errors[dx] = {0.};

    //To this:
    double energy_errors[dx];
    memset( energy_errors, 0, dx*sizeof(int));

    TGraphErrors *graph1 = new TGraphErrors(dx,energy,fg1,energy_errors,fg1err);        
    TGraphErrors *graph2 = new TGraphErrors(dx,energy,fg2,energy_errors,fg2err);
    TGraphErrors *graph3 = new TGraphErrors(dx,energy,fg3,energy_errors,fg3err);
    TGraphErrors *graph4 = new TGraphErrors(dx,energy,fg4,energy_errors,fg4err);
    TGraphErrors *graph5 = new TGraphErrors(dx,energy,fg5,energy_errors,fg5err);
    TGraphErrors *graph6 = new TGraphErrors(dx,energy,fg6,energy_errors,fg6err);

    TGraph *graph1teo = new TGraph(dx,energy,fg1teo);        
       
    c1->Divide(3,2,0,0);
        
    c1->cd(1);
	gPad->SetLeftMargin(0.14);  
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.01);
    h->GetYaxis()->SetTitle("Probability / 160 keV");
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetYaxis()->SetLabelSize(0.065);
    h->GetYaxis()->SetLabelFont(42);
    h->Draw();
	px1->SetLineColor(kBlue);
    px1->SetLineWidth(2);
    graph1->SetMarkerStyle(2);
    graph1->Draw("P");
    px1->Draw("same");
	
    gPad->Update();

    graph1teo->SetMarkerStyle(21);
    //graph1teo->Draw("P same");	
	TLine *line1 = new TLine(20,8427,11000,8427);
	line1->SetLineStyle(3);
    line1->SetLineWidth(2);
	//line1->Draw();
        
	
	TLatex text;
	text.SetTextSize(0.08);
    text.SetTextFont(42);
    text.DrawLatex(4.40,0.149,"(a) E = 5.12 MeV");
	
    TLegend *leg = new TLegend(0.18,0.70,0.40,0.89);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.08);
    leg->AddEntry(graph1," data ","P");
    leg->AddEntry(px1," #rho T ","L");
    leg->Draw();
        
    c1->cd(2);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.01);
    matrix3->Draw();
    h_1->GetYaxis()->SetLabelOffset(100);
    h_1->Draw();
    px2->SetLineColor(kBlue);
    px2->SetLineWidth(2);
    px2->Draw("same");
    graph2->SetMarkerStyle(2);
    graph2->Draw("P same");

    gPad->Update();

    text.DrawLatex(4.40,0.149,"(b) E = 5.44 MeV");

    c1->cd(3);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.01);
    h_1->Draw();
    px3->SetLineColor(kBlue);
    px3->SetLineWidth(2);
    px3->Draw("same");
    graph3->SetMarkerStyle(2);
    graph3->Draw("P same");

    gPad->Update();
        
    text.DrawLatex(4.40,0.149,"(c) E = 5.92 MeV");
    
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
    h2->GetYaxis()->SetTitle("Probability / 160 keV");
    h2->GetYaxis()->SetTitleOffset(0.9);
    h2->GetYaxis()->CenterTitle();
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetYaxis()->SetTitleSize(0.073);
    h2->GetYaxis()->SetLabelSize(0.058);
    h2->GetYaxis()->SetLabelFont(42);
    h2->Draw();
    px4->SetLineColor(kBlue);
    px4->SetLineWidth(2);
    px4->Draw("same");
    graph4->SetMarkerStyle(2);
    graph4->Draw("P same");
        
    gPad->Update();

    text.SetTextSize(0.071);
    text.DrawLatex(4.40,0.149,"(d) E = 6.40 MeV");
    
    
    gPad->Update();
        
    c1->cd(5);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.18);
    h_2->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
    h_2->GetXaxis()->SetTitleOffset(1.1);
    h_2->GetXaxis()->CenterTitle();
    h_2->GetXaxis()->SetTitleFont(42);
    h_2->GetXaxis()->SetTitleSize(0.079);
    h_2->GetXaxis()->SetLabelFont(42);
    h_2->GetXaxis()->SetLabelSize(0.065);
    h_2->GetYaxis()->SetTitle("Probability / 120 keV");
    h_2->GetYaxis()->SetTitleOffset(0.9);
    h_2->GetYaxis()->CenterTitle();
    h_2->GetYaxis()->SetTitleFont(42);
    h_2->GetYaxis()->SetTitleSize(0.073);
    h_2->GetYaxis()->SetLabelOffset(100.058);
    h_2->GetYaxis()->SetLabelFont(42);
    h_2->Draw();
    graph5->SetMarkerStyle(2);
    graph5->Draw("P same");
    px5->SetLineColor(kBlue);
    px5->SetLineWidth(2);
    px5->Draw("same");
    text.DrawLatex(4.40,0.149,"(e) E = 6.88 MeV");

    gPad->Update();
 
    c1->cd(6);
    gPad->SetLeftMargin(0.02);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.18);
    h_3->GetXaxis()->SetLabelFont(42);
    h_3->GetXaxis()->SetLabelSize(0.065);
    h_3->GetYaxis()->SetLabelOffset(42);
    h_3->Draw();
    graph6->SetMarkerStyle(2);
    graph6->Draw("P same");
    px6->SetLineColor(kBlue);
    px6->SetLineWidth(2);
    px6->Draw("same");
    text.DrawLatex(4.40,0.149,"(f) E = 7.20 MeV");

    gPad->Update();

    c1->Update();
	
	// Print to pdf file
	c1->Print("does_it_work_U233.pdf");
		
	
		
}	// END of script
