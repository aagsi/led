#define glx__sim
#include "../src/GlxHit.h"
#include "../src/GlxEvent.h"
#include "../src/GlxLutNode.h"
#include "glxtools.C"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include <iostream>

using namespace std;


Double_t pulse(Double_t x, Double_t par, Double_t width){
	Double_t r= 0;
	Double_t d= fabs(x-par);
	if(d< width) r=1;
	else if(d<width+1) {r = (1-fabs(d-width))*TMath::Tan(45*TMath::Pi()/180.0);
	//std::cout<<"1-fabs(x-2.5)= "<< d <<std::endl;
     }
	return r;
}


Double_t conv(Double_t *x, Double_t * par ){

	// Control constants
	Double_t np = 100.0;      // number of convolution steps
	Double_t sc =   5.0;  

	// Variables
	Double_t xx;
	Double_t intg= par[0];
	Double_t mean = par[1];
	Double_t sigma = 0.128;
	Double_t width = par[2];
	Double_t sum = 0.0;
	Double_t xmin = x[0]-sc*sigma;
	Double_t xmax = x[0]+sc*sigma;
	Double_t step = (xmax-xmin)/np;

	// Convolution integral of Landau and Gaussian by sum
	for(Double_t i=1.0; i<=np/2; i++) {
		xx = xmin+ (i-0.5)* step;
		sum +=  pulse(xx, mean, width)*TMath::Gaus(x[0],xx,sigma);
		xx = xmax-(i-0.5)* step;
		sum +=  pulse(xx, mean, width)*TMath::Gaus(x[0],xx,sigma);
		
		//std::cout<<i<<" gauss= "<< TMath::Gaus(x[0],xx,sigma)<<std::endl;
	}
	//std::cout<<"sum "<< sum <<std::endl;
	return sum * intg ;


}
void tree1r(Int_t ID, TString fileName, TH1F *hPTime[nmcp][npix]);
void m3(TString in= "final_20_1_opt.root" ,TString infile="../build/20m_20_laser_25_1_4.root", Int_t ID= -99, TString out_file= "out.root"){ 
	gSystem->Load("../build/libGlxDirc.so");
	gStyle->SetOptStat(1);
	//t=0;  for i in ../build/1m*.root; do  root -l -q loadlib.C test_laser.C+"(\"$i\",$((++t)),\"out_${t}_.root\")" ;done
	std::cout <<"######################################################################= "<<ID<<std::endl;
	///////////////////////////////////////////////
	///////// variables definition  ///////////////
	///////////////////////////////////////////////
	Double_t evtime;	
	GlxLutNode* node;
	Int_t mcpid, pixid, m, n;
	//CreateMap();
	GlxInit(infile,1);
	//if(infile=="") return;
	const Int_t Nfibers = 1; // 3
	//const int nmcp = 204, npix = 64;
	TRandom fRand1;
	Int_t Angle= 25;
	Double_t r, tr, val, noise;
	Int_t pro;

	/////////////////////////////////////////////
	///////// histo initialization //////////////
	/////////////////////////////////////////////
	TH1F * hPTime[nmcp][npix];
	TH1F * hTime = new TH1F("time",";time, [ns];entries, [#]",500,0,150);
	TH2F * hPMToc[nmcp];
	// for all mcp 
	for(Int_t m=0; m<nmcp; m++){
		// loop over pix
		for(Int_t p=0; p<npix; p++){
			//histogram for each pixl
			hPTime[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),50,-25,25);//50
			axisTime800x500(hPTime[m][p]);
			//gStyle->SetOptTitle(0);
			hPTime[m][p]->SetStats(1);
			hPTime[m][p]->SetLineColor(1);
			hPTime[m][p]->SetLineWidth(3);
		}
	}
	for(int ipad=0;ipad<nmcp;ipad++){
		hPMToc[ipad]  = new TH2F(Form("hPMToc_%d",ipad),Form("PMT %d",ipad),8,0,8,8,0,8); 
		hPMToc[ipad]->SetStats(0);
	}
	const Int_t num=4;
	Double_t x[]={0,1,6,7};
	Double_t y[]={0,1,1,0};
	TGraph* gr= new TGraph(num,x,y);
	//TCanvas* c_g = new TCanvas("c_g","c_g",800,600);
	//gr->Draw("AL*");
	TH2F * randfill = new TH2F("randfill", "randfill",150,-8.0,8.0,100,0,1.2);

	////////////////////////////////
	///////// clone  ///////////////
	////////////////////////////////
	TFile* f = new TFile(infile);
	TTree* t = (TTree*)f->Get("glxlut");
	TClonesArray *fLut[3];
	for(int ifib=0;ifib<Nfibers;ifib++){
		fLut[ifib] = new TClonesArray("GlxLutNode");
		t->SetBranchAddress(Form("LUT_%d",ifib),&fLut[ifib]);
	}
	
	
	
	
	/////////////////////////////////////////////
	///////// histo initialization //////////////
	/////////////////////////////////////////////
	TH1F * hPTime_mean[nmcp][npix];
	// for all mcp 
	for(Int_t m=0; m<nmcp; m++){
		// loop over pix
		for(Int_t p=0; p<npix; p++){
			hPTime_mean[m][p]  = new TH1F(Form("hPTime_mean %d",m*100+p),Form("mcp %d, pixel %d",m, p),200,4.5,7);
			axisTime800x500(hPTime_mean[m][p]);
			hPTime_mean[m][p]->SetStats(1);
			hPTime_mean[m][p]->SetLineColor(1);
			hPTime_mean[m][p]->SetLineWidth(3);
		}
	}
	
	
	Int_t ID_ok, row, col, nMCP, nchannel, nEntries;
	Double_t mean_g_fit,mean_histo, error_mean_g_fit, error_mean_histo, sigma_g_fit, sigma_histo;
	TChain chain1("tout");
	chain1.Add(in);
	chain1.SetBranchAddress("ID",&ID);
	chain1.SetBranchAddress("nEntries",&nEntries);
	chain1.SetBranchAddress("mean_g_fit",&mean_g_fit);
	chain1.SetBranchAddress("mean_histo",&mean_histo);
	chain1.SetBranchAddress("error_mean_g_fit",&error_mean_g_fit);
	chain1.SetBranchAddress("error_mean_histo",&error_mean_histo);
	chain1.SetBranchAddress("sigma_g_fit",&sigma_g_fit);
	chain1.SetBranchAddress("sigma_histo",&sigma_histo);
	chain1.SetBranchAddress("nMCP",&nMCP);
	chain1.SetBranchAddress("row",&row);
	chain1.SetBranchAddress("col",&col);
	chain1.SetBranchAddress("nchannel",&nchannel);
	Int_t n1 = chain1.GetEntries();
	cout<<"\nNb Particles on the the file = "<<n1<<endl;
	
	for (Int_t loop=0; loop<n1; loop++){
		chain1.GetEntry(loop);
		hPTime_mean[nMCP][nchannel-100*nMCP]->Fill(mean_g_fit);
		//std::cout<<" @@@@@@@@@@@@@@@@@@@@@@ mean_g_fit= "<< mean_g_fit<<"	mean_histo="<< hPTime_mean[nMCP][nchannel-100*nMCP]->GetMean()<<std::endl;
	}
	


	
	
	
	
	
	for (Int_t k=0; k< t->GetEntries(); k++) { 
		t->GetEntry(k);	
		//cout<< "K= " <<k << "t->GetEntries= "<< t->GetEntries()<<endl;
		/////////////////////////////////////////////
		///////// Accessing  info  //////////////////
		/////////////////////////////////////////////
		for (Int_t s=0; s< fLut[0]->GetEntriesFast(); s++) { // sensor Id
			GlxLutNode * node = (GlxLutNode*)fLut[0]->At(s);
			Int_t entries=node->Entries();
			//cout<< Form("number of directions in sensor # %d = ",s) <<entries <<endl;
			/////////////////////////////////////////////
			///////// 		Digi	   //////////////////
			/////////////////////////////////////////////
			if(entries > 0){
				mcpid =( s/100 ) -102;
				pixid = s%100 -1;
				m =pixid%8;
				n =7-pixid/8;
				/////////////////////////////////////////////
				///////// 	fill pos histo	   //////////////
				/////////////////////////////////////////////
				//fhDigi[mcpid]->Fill(m, n, entries);
				Double_t counter =0;
				
				


				
				for (Int_t d=0; d< entries; d++) {
					//cout<< "entries " << entries<<endl;
					Int_t pixl_entries = node->GetWeight(d);
					TVector3 dir = node->GetEntry(d);
					//cout<< "XX " <<dir.X() <<endl;
					Double_t evtime = node->GetTime(d);
					Double_t smear= 0.120;  // PMT smearing 120 ps
					evtime+=fRand1.Gaus(0,smear);
					//adding light emission profile of the LED
					do{
						r = gRandom->Rndm();
						//cout<< "r " <<r <<endl;
						tr = (gRandom->Rndm())* 7.0;
						val= gr->Eval(tr) ;
						if(r<=val) randfill->Fill(tr, r);
						if(true){ 
							if(r<=val){
								noise = 0;
								pro= gRandom->Uniform(0,200000); 
								if(pro == 333) {noise =gRandom->Uniform(0,20); hPTime[mcpid][pixid]->Fill(noise); }
								//cout<< "pro= " <<pro <<"	noise= "<<noise<<endl;
								hPTime[mcpid][pixid]->Fill(hPTime_mean[mcpid][pixid]->GetMean()-evtime-tr); // add tr
							}
						} // add tr for LED, remove tr for laser
						if(false){ // switsh ON time offset between sources
							if(k==0) hPTime[mcpid][pixid]->Fill(evtime+tr+0.0);
							if(k==1) hPTime[mcpid][pixid]->Fill(evtime+tr+5.0);
							if(k==2) hPTime[mcpid][pixid]->Fill(evtime+tr+10.0);
						}
					}while(r>val);
				}//
		
			}
		}
	} //

	//TCanvas* ok = new TCanvas("ok","ok",800,600);
	//randfill->Draw("colz");




	
	tree1r(ID, out_file, hPTime); // calling tree1r function 
	Int_t number=0;
	for(Int_t iii=0; iii<nmcp; iii++){
		Int_t number1= fhDigi[iii]->GetEntries();
		if(number1>number) number=number1;
	}
	std::cout<<"number = "<<number<<std::endl;
drawDigi("m,p,v\n",7,0.2,0.001); //sigma error LED 20M
	cDigi->cd();
	//cDigi->SetName(Form("hp_%d",Angle));
	fhDigi[0]->GetZaxis()->SetLabelSize(0.06);
	(new TPaletteAxis(0.890977,0.0620438, 0.97995, 0.952555, ((TH1 *)(fhDigi[49])->Clone())))->Draw();
	canvasAdd(cDigi);
	canvasSave(1,0);
}

void tree1r(Int_t ID, TString fileName, TH1F *hPTime[nmcp][npix]){
	//TH1F * hPTimeError= new TH1F("timeError",";Time Error [ns];Entries [#]",50,0.0,0.140);
	TH1F * hPTimeError= new TH1F("timeError",";Time Error [ns];Entries [#]",50,4.5,7);

	/////////////////////////////////////////////
	///////// Peak searching var  ///////////////
	/////////////////////////////////////////////
	Float_t range_pixle_time, xp_pixle_time=5.5, yp_pixle_time, ra1_pixle_time, rb1_pixle_time, ra2_pixle_time, rb2_pixle_time;
	TSpectrum *s_pixle_time ;
	Int_t nfound_pixle_time, bin_pixle_time ;
	TH1 *hb_pixle_time ;
	Float_t *xpeaks_pixle_time;
	Double_t par_histo[3];
	//Double_t par_histo[6];
	TF1* fit1_histo ;
	//TF1* fit2_histo;
	Double_t mean_g_fit, mean_histo, error_mean_g_fit, error_mean_histo, sigma_g_fit, sigma_histo;
	//Double_t mean2_histo ;
	Int_t row, col, nEntries, nchannel, nMCP;
	////////////////////////////////
	///////// OutPut  //////////////
	////////////////////////////////
	TFile file(fileName,"recreate");
	TTree tout("tout","stuff");
	tout.Branch("ID",&ID,"ID/I");
	tout.Branch("nEntries",&nEntries,"nEntries/I");
	tout.Branch("mean_g_fit",&mean_g_fit,"mean_g_fit/D");
	tout.Branch("mean_histo",&mean_histo,"mean_histo/D");
	tout.Branch("error_mean_g_fit",&error_mean_g_fit,"error_mean_g_fit/D");
	tout.Branch("error_mean_histo",&error_mean_histo,"error_mean_histo/D");
	tout.Branch("sigma_g_fit",&sigma_g_fit,"sigma_g_fit/D");
	tout.Branch("sigma_histo",&sigma_histo,"sigma_histo/D");
	tout.Branch("nMCP",&nMCP,"nMCP/I");
	tout.Branch("row",&row,"row/I");
	tout.Branch("col",&col,"col/I");
	tout.Branch("nchannel",&nchannel,"nchannel/I");

	//gStyle->SetOptStat("e"); // e = 1; number of entries printed
	gStyle->SetOptFit(1); // The type of information about fit parameters printed in the histogram statistics box can be selected via the parameter mode.
	SetRootPalette(1);
	/////////////////////////////////////////////
	///////// Pixel Time processing   ///////////
	/////////////////////////////////////////////
	range_pixle_time= 20;
	TCanvas* c_histo = new TCanvas("c_histo","c_histo",800,600);
	for (Int_t m=0; m< 102; m++) {
		//for (Int_t m=50; m< 51; m++) {
		for (Int_t p=0; p< npix; p++) {
			//for (Int_t p=0; p< 1; p++) {
			// Peak finders Pixle Time plots
			s_pixle_time = new TSpectrum(1);
			nfound_pixle_time = s_pixle_time->Search(hPTime[m][p],2,"",0.10);
			printf("Found %d candidate peaks to fit\n",nfound_pixle_time);
			//Estimate background using TSpectrum::Background
			//hb_pixle_time = s_pixle_time->Background(hPTime[m][p],200,"same");
			//Loop on all found peaks. Eliminate peaks at the background level
			xpeaks_pixle_time =(Float_t*)(s_pixle_time->GetPositionX());
			for (Int_t p_pixle_time=0;p_pixle_time<nfound_pixle_time;p_pixle_time++) {
				xp_pixle_time = xpeaks_pixle_time[p_pixle_time];
				bin_pixle_time = hPTime[m][p]->GetXaxis()->FindBin(xp_pixle_time);
				yp_pixle_time = hPTime[m][p]->GetBinContent(bin_pixle_time);
				//std::cout <<"############xp (pixle_time)= "<<xp_pixle_time<<"         bin (pixle_time)= "<<bin_pixle_time<<"              yp (pixle_time)= "<<yp_pixle_time<<std::endl;
			}
			ra1_pixle_time= xpeaks_pixle_time[1]- range_pixle_time;
			rb1_pixle_time= xpeaks_pixle_time[1]+ range_pixle_time;
			ra2_pixle_time= xpeaks_pixle_time[0]- range_pixle_time;
			rb2_pixle_time= xpeaks_pixle_time[0]+ range_pixle_time;
			//fit2_histo = new TF1("fit2_histo","gaus",ra2_pixle_time,rb2_pixle_time);
			//hPTime[m][p]->Fit(fit2_histo,"R+");
			//fit2_histo->GetParameters(&par_histo[3]);
			// mean2_histo = fit2_histo->GetParameter(1);


			//fit1_histo = new TF1("fit1_histo","gaus",ra1_pixle_time,rb1_pixle_time);
			//hPTime[m][p]->Fit(fit1_histo,"R");
			
			fit1_histo = new TF1("fit1_histo",conv,0,15,3); //0, 20 the range , 2 number of parameters on the function 
			fit1_histo->SetNpx(500);
			fit1_histo->SetParameters(1,xp_pixle_time, 2.5);
			fit1_histo->SetParLimits(2,2,3);
			fit1_histo->SetParNames ("Constant","T0","1/2 Width");

			
			hPTime[m][p]->Fit(fit1_histo);
			fit1_histo->GetParameters(&par_histo[0]);
			mean_g_fit = fit1_histo->GetParameter(1);
			error_mean_g_fit = fit1_histo->GetParError(1);
			sigma_g_fit = fit1_histo->GetParameter(2);
			nEntries= hPTime[m][p]->GetEntries();
			sigma_histo = hPTime[m][p]->GetRMS();
			mean_histo= hPTime[m][p]->GetMean();
			error_mean_histo= hPTime[m][p]->GetMeanError();

			row =p%8;
			col =7-p/8;
			nchannel= m*100+p;
			nMCP=m;
			/////////////////////////////////////////////
			///////// 	fill time histos   //////////////
			/////////////////////////////////////////////
			if(hPTime[m][p]->GetEntries() > 60 ){
				//fhDigi[m]->Fill(row, col,(mean_g_fit+ mean2_histo)/2); // test
				//fhDigi[m]->Fill(row, col,hPTime[m][p]->GetMean()); // mean W/O fitting
				//fhDigi[m]->Fill(row, col,hPTime[m][p]->GetMeanError()); //Error mean W/O fitting
				fhDigi[m]->Fill(row, col,sigma_g_fit); // mean With fitting 
				//fhDigi[m]->Fill(row, col,error_mean_g_fit); // Error mean With fitting 
				
				hPTime[m][p]->Draw();
				c_histo->Update();
				c_histo->WaitPrimitive();
				
				//hPTimeError->Fill(mean_histo);
				tout.Fill();
			}
		}
	}
	tout.Write();
	//TCanvas* c_hPTimeError = new TCanvas("c_hPTimeError","c_hPTimeError",800,600);
	//hPTimeError->Draw();
}

