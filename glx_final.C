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
void tree1r(TString fileName, TH1F *hPTime_0[nmcp][npix], TH1F *hPTime_1[nmcp][npix], TH1F *hPTime_2[nmcp][npix]);
//void glx_final(TString in0= "final_0.root", TString in1= "final_1.root", TString in2= "final_2.root", TString out_file= "out.root"){
void glx_final(TString in0= "led_0.root", TString in1= "led_1.root", TString in2= "led_2.root", TString out_file= "out.root"){
	gSystem->Load("../build/libGlxDirc.so");
	gStyle->SetOptStat(1);
	///////////////////////////////////////////////
	///////// variables definition  ///////////////
	///////////////////////////////////////////////
	Double_t evtime;
	GlxLutNode* node;
	Int_t mcpid, pixid, m, n;
	//CreateMap();
	GlxInit(in0,1);
	//if(infile=="") return;
	const Int_t Nfibers = 1; // 3
	//const int nmcp = 204, npix = 64;
	TRandom fRand1;
	Int_t Angle= 25;
	Double_t r, tr, val;

	/////////////////////////////////////////////
	///////// histo initialization //////////////
	/////////////////////////////////////////////

	TH1F * hPTime_0[nmcp][npix];
	TH1F * hPTime_1[nmcp][npix];
	TH1F * hPTime_2[nmcp][npix];
	
	TH1F * hTime = new TH1F("time",";time, [ns];entries, [#]",500,0,150);
	// for all mcp 
	for(Int_t m=0; m<nmcp; m++){
		// loop over pix
		for(Int_t p=0; p<npix; p++){
			//histogram for each pixl
			//hPTime_1[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),200,4,8);//50
			hPTime_1[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),5000,4.5,20);
			axisTime800x500(hPTime_1[m][p]);
			// gStyle->SetOptTitle(0);
			hPTime_1[m][p]->SetStats(1);
			hPTime_1[m][p]->SetLineColor(1);
			hPTime_1[m][p]->SetLineWidth(3);
			
			hPTime_0[m][p]  = new TH1F(Form("hPTime_0_%d",m*100+p),Form("mcp_0 %d, pixel_0 %d",m, p),5000,4.5,20);
			axisTime800x500(hPTime_0[m][p]);
			hPTime_0[m][p]->SetStats(1);
			hPTime_0[m][p]->SetLineColor(1);
			hPTime_0[m][p]->SetLineWidth(3);
			
			hPTime_2[m][p]  = new TH1F(Form("hPTime_2_%d",m*100+p),Form("mcp_2 %d, pixel_2 %d",m, p),5000,4.5,20);
			axisTime800x500(hPTime_2[m][p]);
			hPTime_2[m][p]->SetStats(1);
			hPTime_2[m][p]->SetLineColor(1);
			hPTime_2[m][p]->SetLineWidth(3);
	
		}
	}
	///////////////////////////////////////
	///////// tree variables //////////////
	///////////////////////////////////////

	Int_t ID, row, col, nMCP, nchannel, nEntries;
	Double_t mean_g_fit,mean_histo, error_mean_g_fit, error_mean_histo, sigma_g_fit, sigma_histo;

	TChain chain0("tout");
	chain0.Add(in0);
	chain0.SetBranchAddress("ID",&ID);
	chain0.SetBranchAddress("nEntries",&nEntries);
	chain0.SetBranchAddress("mean_g_fit",&mean_g_fit);
	chain0.SetBranchAddress("mean_histo",&mean_histo);
	chain0.SetBranchAddress("error_mean_g_fit",&error_mean_g_fit);
	chain0.SetBranchAddress("error_mean_histo",&error_mean_histo);
	chain0.SetBranchAddress("sigma_g_fit",&sigma_g_fit);
	chain0.SetBranchAddress("sigma_histo",&sigma_histo);
	chain0.SetBranchAddress("nMCP",&nMCP);
	chain0.SetBranchAddress("row",&row);
	chain0.SetBranchAddress("col",&col);
	chain0.SetBranchAddress("nchannel",&nchannel);
	
	
	TChain chain1("tout");
	chain1.Add(in1);
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
	
	TChain chain2("tout");
	chain2.Add(in2);
	chain2.SetBranchAddress("ID",&ID);
	chain2.SetBranchAddress("nEntries",&nEntries);
	chain2.SetBranchAddress("mean_g_fit",&mean_g_fit);
	chain2.SetBranchAddress("mean_histo",&mean_histo);
	chain2.SetBranchAddress("error_mean_g_fit",&error_mean_g_fit);
	chain2.SetBranchAddress("error_mean_histo",&error_mean_histo);
	chain2.SetBranchAddress("sigma_g_fit",&sigma_g_fit);
	chain2.SetBranchAddress("sigma_histo",&sigma_histo);
	chain2.SetBranchAddress("nMCP",&nMCP);
	chain2.SetBranchAddress("row",&row);
	chain2.SetBranchAddress("col",&col);
	chain2.SetBranchAddress("nchannel",&nchannel);
	

	Int_t n0 = chain0.GetEntries();
	//cout<<"\nNb Particles on the the file = "<<n0<<endl;
	
	Int_t n1 = chain1.GetEntries();
	//cout<<"\nNb Particles on the the file = "<<n1<<endl;
	
	
	Int_t n2 = chain2.GetEntries();
	//cout<<"\nNb Particles on the the file = "<<n2<<endl;

//////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

	for (Int_t t=0; t<n0; t++){
		chain0.GetEntry(t);
	if(nEntries > 60 ){
		hPTime_0[nMCP][nchannel-100*nMCP]->Fill(mean_g_fit);
		}
	}
	
	
	for (Int_t e=0; e<n1; e++){
		chain1.GetEntry(e);
	if(nEntries > 60 ){
		hPTime_1[nMCP][nchannel-100*nMCP]->Fill(mean_g_fit);
		}
	}
	
	
	for (Int_t d=0; d<n2; d++){
		chain2.GetEntry(d);
	if(nEntries > 60 ){
		hPTime_2[nMCP][nchannel-100*nMCP]->Fill(mean_g_fit);
		}
	}

	tree1r(out_file, hPTime_0, hPTime_1, hPTime_2); // calling tree1r function 
	

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
//drawDigi("m,p,v\n",7,0.120,0.005); //error LED 20M



//drawDigi("m,p,v\n",7,0.2,0.001); //sigma error LED 20M  // 0.001
drawDigi("m,p,v\n",7,10.0,8.4); //time LED

//drawDigi("m,p,v\n",7,6.4,5); //time LED mix all sources 6.4 4.5
	

	cDigi->cd();
	//cDigi->SetName(Form("hp_%d",Angle));
	fhDigi[0]->GetZaxis()->SetLabelSize(0.06);
	(new TPaletteAxis(0.890977,0.0620438, 0.97995, 0.952555, ((TH1 *)(fhDigi[49])->Clone())))->Draw();
	canvasAdd(cDigi);
	canvasSave(1,0);


}


void tree1r(TString fileName, TH1F *hPTime_0[nmcp][npix], TH1F *hPTime_1[nmcp][npix], TH1F *hPTime_2[nmcp][npix]){
	TH1F * hPTime_rms= new TH1F("hPTime_rms",";T0 rms [ns];Entries [#]",50,0.0,0.200);
	TH1F * hPTime_t0= new TH1F("hPTime_t0",";T0 [ns];Entries [#]",120,4.5,20);
	
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
	Double_t mean_histo, sigma_histo;
	Double_t mean_histo_0, sigma_histo_0;
	Double_t mean_histo_1, sigma_histo_1;
	Double_t mean_histo_2, sigma_histo_2;
	//Double_t mean2_histo ;
	Int_t row, col, nchannel, nMCP, nEntries_0, nEntries_1, nEntries_2;
	////////////////////////////////
	///////// OutPut  //////////////
	////////////////////////////////

	gStyle->SetOptStat(1); // e = 1; number of entries printed
	gStyle->SetOptFit(1); // The type of information about fit parameters printed in the histogram statistics box can be selected via the parameter mode.
	SetRootPalette(1);
	/////////////////////////////////////////////
	///////// Pixel Time processing   ///////////
	/////////////////////////////////////////////

	TCanvas* c_histo = new TCanvas("c_histo","c_histo",800,600);
	for (Int_t m=0; m< 102; m++) {
	//for (Int_t m=33; m< 34; m++) {
		for (Int_t p=0; p< npix; p++) {
		//for (Int_t p=2; p< 3; p++) {
		
	
			nEntries_0= hPTime_0[m][p]->GetBinContent(hPTime_0[m][p]->GetMaximumBin());
			nEntries_1= hPTime_1[m][p]->GetBinContent(hPTime_1[m][p]->GetMaximumBin());
			nEntries_2= hPTime_2[m][p]->GetBinContent(hPTime_2[m][p]->GetMaximumBin());
			
			
			
			sigma_histo_0 = hPTime_0[m][p]->GetRMS();
			mean_histo_0  = hPTime_0[m][p]->GetMean();
		
			sigma_histo_1 = hPTime_1[m][p]->GetRMS();
			mean_histo_1  = hPTime_1[m][p]->GetMean();
		
			sigma_histo_2 = hPTime_2[m][p]->GetRMS();
			mean_histo_2  = hPTime_2[m][p]->GetMean();
			
			/*
			
			if ((sigma_histo_0>0&&sigma_histo_1>0) ||(sigma_histo_2>0&&sigma_histo_1>0))
			{
				
				
				
	if(sigma_histo_0 <= sigma_histo_1 && sigma_histo_0<= sigma_histo_2)
    {
		sigma_histo = sigma_histo_0;
		mean_histo= mean_histo_0;
    }

    if(sigma_histo_1 <= sigma_histo_0 && sigma_histo_1 <= sigma_histo_2)
    {
		sigma_histo = sigma_histo_1;
		mean_histo= mean_histo_1;
    }

    if(sigma_histo_2 <= sigma_histo_0 && sigma_histo_2 <= sigma_histo_1) {
		
		sigma_histo = sigma_histo_2;
		mean_histo= mean_histo_2;
    }
				
	}			
			
				
		else{		// remove frist if condition 
			
		*/	
			
			
	
    if(nEntries_0 >= nEntries_1 && nEntries_0>= nEntries_2)
    {
		sigma_histo = sigma_histo_0;
		mean_histo= mean_histo_0;
    }

    if(nEntries_1 >= nEntries_0 && nEntries_1 >= nEntries_2)
    {
		sigma_histo = sigma_histo_1;
		mean_histo= mean_histo_1;
    }

    if(nEntries_2 >= nEntries_0 && nEntries_2 >= nEntries_1) {
		
		sigma_histo = sigma_histo_2;
		mean_histo= mean_histo_2;
    }
    
		
	//}
		
		cout<<"nEntries_0 = "<<nEntries_0 << "	nEntries_1 ="<< nEntries_1<< "	nEntries_2 =" << nEntries_2 <<endl;
		cout<<"mean_histo_0 = "<<mean_histo_0 << "	mean_histo_1 ="<< mean_histo_1<< "	mean_histo_2 =" << mean_histo_2 <<endl;
		cout<<"mean_histo = "<< mean_histo <<endl;

			row =p%8;
			col =7-p/8;
			nchannel= m*100+p;
			nMCP=m;
			
			
			/////////////////////////////////////////////
			///////// 	fill time histos   //////////////
			/////////////////////////////////////////////
			//if(hPTime_0[m][p]->GetEntries() == 1000 && hPTime_1[m][p]->GetEntries() == 1000 && hPTime_2[m][p]->GetEntries() == 1000){
				
				//fhDigi[m]->Fill(row, col,sigma_histo);
				//hPTime_rms->Fill(sigma_histo);
				fhDigi[m]->Fill(row, col,mean_histo);
				hPTime_t0->Fill(mean_histo);
				
				/*
				
				//hPTime_0[m][p]->Draw();
				//hPTime_1[m][p]->Draw();
				//hPTime_2[m][p]->Draw();
				
				c_histo->Update();
				c_histo->WaitPrimitive();
				*/
				
			

		//	}
		}
	}
	TCanvas* c_hPTime_rms = new TCanvas("c_hPTime_rms","c_hPTime_rms",800,600);
	hPTime_rms->Draw();
	TCanvas* c_hPTime_t0 = new TCanvas("c_hPTime_t0","c_hPTime_t0",800,600);
	hPTime_t0->Draw();
		
}
