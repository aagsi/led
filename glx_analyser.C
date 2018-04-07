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

// define prototype 
void tree1r(TString fileName, TH1F *hPTime[nmcp][npix]);
void glx_analyser(TString in= "final_20.root" , TString out_file= "out.root"){
	gSystem->Load("../build/libGlxDirc.so");
	gStyle->SetOptStat(1);
	///////////////////////////////////////////////
	///////// variables definition  ///////////////
	///////////////////////////////////////////////
	Double_t evtime;
	GlxLutNode* node;
	Int_t mcpid, pixid, m, n;
	//CreateMap();
	GlxInit(in,1);
	//if(infile=="") return;
	const Int_t Nfibers = 1; // 3
	//const int nmcp = 204, npix = 64;
	TRandom fRand1;
	Int_t Angle= 25;
	Double_t r, tr, val;

	/////////////////////////////////////////////
	///////// histo initialization //////////////
	/////////////////////////////////////////////
	TH1F * hPTime[nmcp][npix];
	TH1F * hPTimeError= new TH1F("timeError",";Time Error [ns];Entries [#]",150,0,0.200);
	TH1F * hTime = new TH1F("time",";time, [ns];entries, [#]",500,0,150);
	// for all mcp 
	for(Int_t m=0; m<nmcp; m++){
		// loop over pix
		for(Int_t p=0; p<npix; p++){
			//histogram for each pixl
			//hPTime[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),200,0,0.200);//50
			hPTime[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),4,-1,3);//50
			axisTime800x500(hPTime[m][p]);
			// gStyle->SetOptTitle(0);
			hPTime[m][p]->SetStats(1);
			hPTime[m][p]->SetLineColor(1);
			hPTime[m][p]->SetLineWidth(3);
		}
	}

	///////////////////////////////////////
	///////// tree variables //////////////
	///////////////////////////////////////

	Int_t ID, row, col, nMCP, nchannel, nEntries;
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

	TCanvas* c_histo = new TCanvas("c_histo","c_histo",800,600);
	Int_t n1 = chain1.GetEntries();
	cout<<"\nNb Particles on the the file = "<<n1<<endl;

//////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

	for (Int_t t=0; t<n1; t++){
		chain1.GetEntry(t);
	if(nEntries > 60 ){
		hPTime[nMCP][nchannel-100*nMCP]->Fill(error_mean_g_fit);
		hPTimeError->Fill(error_mean_g_fit);
		}
	//	hPTime[nMCP][nchannel-100*nMCP]->Fill(error_mean_histo);

		//fhDigi[nMCP]->Fill(row, col, mean_histo);
	}
TCanvas* c_hPTimeError = new TCanvas("c_hPTimeError","c_hPTimeError",800,600);
hPTimeError->Draw();


	tree1r(out_file, hPTime); // calling tree1r function 
	Int_t number=0;
	for(Int_t iii=0; iii<nmcp; iii++){
		Int_t number1= fhDigi[iii]->GetEntries();
		if(number1>number) number=number1;
	}
	std::cout<<"number = "<<number<<std::endl;
	//drawDigi("m,p,v\n",7,0.0,0.); //pos
	//drawDigi("m,p,v\n",7,6.4,5.); //time LED
	//drawDigi("m,p,v\n",7,0.04,0.01); //error LED  //0.04
	//drawDigi("m,p,v\n",7,0.,1.5); //time Laser
	//drawDigi("m,p,v\n",7,0.05,0.001); //error Laser all
	
	
	
	
	//drawDigi("m,p,v\n",7,0.120,0.016); //error LED all
	drawDigi("m,p,v\n",7,0.060,0.016); //error LED 20M
	
	
	
	
	cDigi->cd();
	//cDigi->SetName(Form("hp_%d",Angle));
	fhDigi[0]->GetZaxis()->SetLabelSize(0.06);
	(new TPaletteAxis(0.890977,0.0620438, 0.97995, 0.952555, ((TH1 *)(fhDigi[49])->Clone())))->Draw();
	canvasAdd(cDigi);
	canvasSave(1,0);


}


void tree1r(TString fileName, TH1F *hPTime[nmcp][npix]){
	/////////////////////////////////////////////
	///////// Peak searching var  ///////////////
	/////////////////////////////////////////////
	Float_t range_pixle_time, xp_pixle_time, yp_pixle_time, ra1_pixle_time, rb1_pixle_time, ra2_pixle_time, rb2_pixle_time;
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

	gStyle->SetOptStat("e"); // e = 1; number of entries printed
	gStyle->SetOptFit(1); // The type of information about fit parameters printed in the histogram statistics box can be selected via the parameter mode.
	SetRootPalette(1);
	/////////////////////////////////////////////
	///////// Pixel Time processing   ///////////
	/////////////////////////////////////////////
	range_pixle_time= 0.1;
	TCanvas* c_histo = new TCanvas("c_histo","c_histo",800,600);
	for (Int_t m=0; m< 102; m++) {
	//for (Int_t m=41; m< 42; m++) {
		for (Int_t p=0; p< npix; p++) {
		//for (Int_t p=25; p< 26; p++) {
			// Peak finders Pixle Time plots
			s_pixle_time = new TSpectrum(1);
			nfound_pixle_time = s_pixle_time->Search(hPTime[m][p],2,"",0.10);
			printf("Found %d candidate peaks to fit\n",nfound_pixle_time);
			//Estimate background using TSpectrum::Background
			hb_pixle_time = s_pixle_time->Background(hPTime[m][p],200,"same");
			//Loop on all found peaks. Eliminate peaks at the background level
			xpeaks_pixle_time =(Float_t*)(s_pixle_time->GetPositionX());
			for (Int_t p_pixle_time=0;p_pixle_time<nfound_pixle_time;p_pixle_time++) {
				xp_pixle_time = xpeaks_pixle_time[p_pixle_time];
				bin_pixle_time = hPTime[m][p]->GetXaxis()->FindBin(xp_pixle_time);
				yp_pixle_time = hPTime[m][p]->GetBinContent(bin_pixle_time);
				//std::cout <<"xp (pixle_time)= "<<xp_pixle_time<<"         bin (pixle_time)= "<<bin_pixle_time<<"              yp (pixle_time)= "<<yp_pixle_time<<std::endl;
			}
			ra1_pixle_time= xpeaks_pixle_time[1]- range_pixle_time;
			rb1_pixle_time= xpeaks_pixle_time[1]+ range_pixle_time;
			ra2_pixle_time= xpeaks_pixle_time[0]- range_pixle_time;
			rb2_pixle_time= xpeaks_pixle_time[0]+ range_pixle_time;
			//fit2_histo = new TF1("fit2_histo","gaus",ra2_pixle_time,rb2_pixle_time);
			//hPTime[m][p]->Fit(fit2_histo,"R+");
			//fit2_histo->GetParameters(&par_histo[3]);
			// mean2_histo = fit2_histo->GetParameter(1);
			fit1_histo = new TF1("fit1_histo","gaus",ra1_pixle_time,rb1_pixle_time);
			hPTime[m][p]->Fit(fit1_histo,"R");

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
			if(hPTime[m][p]->GetEntries() == 1000 ){
				//fhDigi[m]->Fill(row, col,hPTime[m][p]->GetRMS()); // RMS W/O fitting
                                //fhDigi[m]->Fill(row, col,sigma_g_fit); // SIGMA 
				//fhDigi[m]->Fill(row, col,(mean_g_fit+ mean2_histo)/2); // test
				fhDigi[m]->Fill(row, col,hPTime[m][p]->GetMean()); // mean W/O fitting
				//fhDigi[m]->Fill(row, col,hPTime[m][p]->GetMeanError()); //Error mean W/O fitting
				//fhDigi[m]->Fill(row, col,mean_g_fit); // mean With fitting 
				//fhDigi[m]->Fill(row, col,error_mean_g_fit); // Error mean With fitting 
				/*
				hPTime[m][p]->Draw();
				c_histo->Update();
				c_histo->WaitPrimitive();
			*/
				tout.Fill();
			}
		}
	}
	tout.Write();
}
// mean sigma
// mean RMS
// error mean
