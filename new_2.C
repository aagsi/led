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
#include "TMultiGraph.h"
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
	Double_t sigma = 0.120;
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

void new_2(TString infile="../build/20m_20_laser_25_1_4.root", Int_t ID= -99, 
																				bool Canvmg=true, 
																				bool CanvFitBool=false, 
																				bool CanvSigmaBool=false , 
																				bool CanvSigmaLEDBool=false )

{	
	gSystem->Load("../build/libGlxDirc.so");
	gStyle->SetOptStat(1);
	//t=0;  for i in ../build/1m*.root; do  root -l -q loadlib.C test_laser.C+"(\"$i\",$((++t)),\"out_${t}_.root\")" ;done
	std::cout <<"######################################################################= "<<ID<<std::endl;
	///////////////////////////////////////////////
	///////// Fitting function   ///////////////
	///////////////////////////////////////////////
	TF1* fit1_histo ;
	//TF1* fit_gaus;
	fit1_histo = new TF1("fit1_histo",conv,0,15,3); //0, 20 the range , 2 number of parameters on the function 
	fit1_histo->SetNpx(500);
	fit1_histo->SetParLimits(2,2,3);
	fit1_histo->SetParNames ("Constant","Mean LED","1/2FWHM");
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
	Double_t sigma=0;
	Double_t sigmaLED=0;
	//////////////////////////////////////////////////////
	///////// histo & graphs initialization //////////////
	//////////////////////////////////////////////////////
	TGraph * SigmaGraph[nmcp][npix], * SigmaGraphLED[nmcp][npix];
	TMultiGraph *mg[nmcp][npix];
	TH1F * hPTime[nmcp][npix], * hPTimeSigma[nmcp][npix], * hPTimeSigmaLED[nmcp][npix],* hPTimeLED[nmcp][npix];
	// for all mcp 
	for(Int_t m=0; m<nmcp; m++){
		// loop over pix
		for(Int_t p=0; p<npix; p++){
			//histogram for each pixl
			hPTime[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),50,-25,25);//50   25,-10,15
			hPTimeLED[m][p]  = new TH1F(Form("hPTimeLED_%d",m*100+p),Form("mcp %d, pixel %d",m, p),50,-25,25);//50 
			hPTimeSigma[m][p]  = new TH1F(Form("hPTimeSigma_%d",m*100+p),Form("mcp %d, pixel %d",m, p),400,-0.4,0.4); // laser
			hPTimeSigmaLED[m][p]  = new TH1F(Form("hPTimeSigmaLED_%d",m*100+p),Form("mcp %d, pixel %d",m, p),400,-0.4,0.4); //LED
			//
			axisTime800x500(hPTime[m][p]);
			axisTime800x500(hPTimeLED[m][p]);
			//gStyle->SetOptTitle(0);
			hPTime[m][p]->SetStats(1);
			hPTime[m][p]->SetLineColor(1);
			hPTime[m][p]->SetLineWidth(3);
			//
			hPTimeLED[m][p]->SetStats(1);
			hPTimeLED[m][p]->SetLineColor(1);
			hPTimeLED[m][p]->SetLineWidth(3);
			// graph declaration
			SigmaGraph[m][p]= new TGraph();
			SigmaGraph[m][p]->SetTitle(Form("Laser mcp %d, pixel %d",m, p));
			//
			SigmaGraphLED[m][p]= new TGraph();
			SigmaGraphLED[m][p]->SetTitle(Form("LED mcp %d, pixel %d",m, p));
			//
			mg[m][p]= new TMultiGraph();
			mg[m][p]->SetTitle(Form("mcp %d, pixel %d",m, p));
		}
	}
	const Int_t num=4;
	//Double_t x[]={0,1,6,7};
	Double_t x[]={0,2,12,14};
	Double_t y[]={0,1,1,0};
	TGraph* gr= new TGraph(num,x,y);
	TCanvas* c_g = new TCanvas("c_g","c_g",800,600);
	gr->Draw("AL*");
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
	for (Int_t k=0; k< t->GetEntries(); k++) { 
		t->GetEntry(k);	
		//cout<< "K= " <<k << "t->GetEntries= "<< t->GetEntries()<<endl;
		/////////////////////////////////////////////
		///////// Accessing  info  //////////////////
		/////////////////////////////////////////////
		for (Int_t s=0; s< fLut[0]->GetEntriesFast(); s++) { // sensor Id
			//		if(s!= 12556 ) continue; // 12524   // 12556
			GlxLutNode * node = (GlxLutNode*)fLut[0]->At(s);
			Int_t entries=node->Entries();
			//cout<< Form("number of directions in sensor # %d = ",s) <<entries <<endl;
			/////////////////////////////////////////////
			///////// 		Digi	   //////////////////
			/////////////////////////////////////////////
			if(entries > 1000){
				mcpid =( s/100 ) -102;
				pixid = s%100 -1;
				m =pixid%8;
				n =7-pixid/8;
				if(mcpid != 50) continue; //54 //48
				//cout<< "#####################################################entries " << entries<< "	s	"<<s<<endl;
				/////////////////////////////////////////////
				///////// 	fill pos histo	   //////////////
				/////////////////////////////////////////////
				Int_t npoint=0;
				for (Int_t step=0; step< 20; step++) {
					entries= step * 100;
					for (Int_t d=0; d< entries; d++) {
						//cout<< "entries " << entries<<endl;
						//TVector3 dir = node->GetEntry(d);
						//cout<< "XX " <<dir.X() <<endl;
						Double_t evtime = node->GetTime(d);
						Double_t smear= 0.120;  // PMT smearing 120 ps
						evtime+=fRand1.Gaus(0,smear);
						//adding light emission profile of the LED
						do{
							r = gRandom->Rndm();
							//cout<< "r " <<r <<endl;
							tr = (gRandom->Rndm())* 14.0; //14
							val= gr->Eval(tr);
							if(r<=val) randfill->Fill(tr, r);
							if(true){ 
								if(r<=val){
									noise = 0;
									pro= gRandom->Uniform(0,200000); 
									if(pro == 333) {noise =gRandom->Uniform(0,20); hPTimeLED[mcpid][pixid]->Fill(noise); }
									//cout<< "pro= " <<pro <<"	noise= "<<noise<<endl;
									hPTimeLED[mcpid][pixid]->Fill(evtime+tr-7.0); // add 3.5 7 
								}
							} // add tr for LED, remove tr for laser
						}while(r>val);
						hPTime[mcpid][pixid]->Fill(evtime); // laser
					}
					TCanvas* canvFit;
					if (CanvFitBool) canvFit = new TCanvas("canvFit","canvFit",800,600);
					if(hPTime[mcpid][pixid]->GetEntries() < 10 ) continue;
					if(hPTimeLED[mcpid][pixid]->GetEntries() < 10 ) continue;
					cout<< "number hPTime= " <<hPTime[mcpid][pixid]->GetEntries()<<endl;
					Double_t mean_g_fit;
					Double_t mean_g_fit_LED;
					////////////////////////////////
					///////// Fitting  /////////////
					////////////////////////////////
					// Laser stuff
					TF1* fit_gaus = new TF1("fit_gaus","[0]*exp(-0.5*((x-[1])/[2])**2)", -10, 15); 
					fit_gaus->SetNpx(500);
					//fit_gaus->FixParameter(1,2.3);
					//fit_gaus->FixParameter(2,0.5);
					//fit_gaus->FixParameter(0,25000);
					fit_gaus->SetParLimits(2,0.5, 0.8);
					fit_gaus->SetParLimits(1,2.4, 2.6);
					fit_gaus->SetParameter(2, 0.5);
					fit_gaus->SetParameter(1, 2.5);
					cout<< "*****************************************************************************************" <<endl;
					hPTime[mcpid][pixid]->Fit("fit_gaus");
					mean_g_fit = fit_gaus->GetParameter(1);
					// LED Stuff
					fit1_histo->SetParameters(1,hPTimeLED[mcpid][pixid]->GetMean(), 5.5); // 2.5 5.5
					fit1_histo->SetParLimits(2,5.0,6.0); // commint the line
					hPTimeLED[mcpid][pixid]->Fit(fit1_histo);
					mean_g_fit_LED = fit1_histo->GetParameter(1);
					hPTime[mcpid][pixid]->Draw();
					hPTimeLED[mcpid][pixid]->Draw("same");
					fit_gaus->Draw("same");
					if (CanvFitBool) canvFit->Update();
					if (CanvFitBool) canvFit->WaitPrimitive();
					cout<< "################ mean_g_fit= " <<mean_g_fit<<endl;
					cout<< "################ mean_g_fit_LED= " <<mean_g_fit_LED<<endl;
					////////////////////////////////
					///////// Filling Delta  ///////
					////////////////////////////////           
					for (Int_t d=0; d< entries; d++) {
						Double_t evtime = node->GetTime(d);
						hPTimeSigma[mcpid][pixid]->Fill(mean_g_fit-evtime);
						hPTimeSigmaLED[mcpid][pixid]->Fill(mean_g_fit_LED-evtime);
					}
					////////////////////////////////
					///////// Fitting Delta  ///////
					////////////////////////////////
					TCanvas* CanvashPTimeSigma;
					if (CanvSigmaBool) CanvashPTimeSigma = new TCanvas("CanvashPTimeSigma","CanvashPTimeSigma",800,600);
					cout<< "integral hPTimeSigma= " <<hPTimeSigma[mcpid][pixid]->Integral(-5,5)<<endl;
					hPTimeSigma[mcpid][pixid]->Fit("gaus");
					hPTimeSigma[mcpid][pixid]->Draw();
					if (CanvSigmaBool) CanvashPTimeSigma->Update();
					if (CanvSigmaBool) CanvashPTimeSigma->WaitPrimitive();
					//
					TCanvas* CanvashPTimeSigmaLED;
					if (CanvSigmaLEDBool)CanvashPTimeSigmaLED = new TCanvas("CanvashPTimeSigmaLED","CanvashPTimeSigmaLED",800,600);
					cout<< "number hPTimeSigmaLED= " <<hPTimeSigmaLED[mcpid][pixid]->GetEntries()<<endl;
					cout<< "integral hPTimeSigmaLED= " <<hPTimeSigmaLED[mcpid][pixid]->Integral(-5,5)<<endl;
					hPTimeSigmaLED[mcpid][pixid]->Fit("gaus");
					hPTimeSigmaLED[mcpid][pixid]->Draw("same");
					if (CanvSigmaLEDBool) CanvashPTimeSigmaLED->Update();
					if (CanvSigmaLEDBool) CanvashPTimeSigmaLED->WaitPrimitive();
					//
					if(hPTimeSigma[mcpid][pixid]->Integral(5,195)== 0) continue;
					if(hPTimeSigmaLED[mcpid][pixid]->Integral(5,195)== 0) continue;
					////////////////////////////////
					///////// Filling Graphs  //////
					////////////////////////////////
					sigma = hPTimeSigma[mcpid][pixid]->GetFunction("gaus")->GetParameter(2);
					sigmaLED = hPTimeSigmaLED[mcpid][pixid]->GetFunction("gaus")->GetParameter(2);
					//sigma = hPTimeSigma[mcpid][pixid]->GetRMS();
					//sigmaLED = hPTimeSigmaLED[mcpid][pixid]->GetRMS();
					cout<< "npoint = " <<npoint<<endl;
					SigmaGraph[mcpid][pixid]->SetPoint(npoint,entries,sigma);
					SigmaGraphLED[mcpid][pixid]->SetPoint(npoint,entries,sigmaLED);
					hPTimeSigma[mcpid][pixid]->Reset();
					hPTimeSigmaLED[mcpid][pixid]->Reset();
					npoint++;	              
				}// end of steps loop
				//
				SigmaGraph[mcpid][pixid]->GetXaxis()->SetTitle("Stat");
				SigmaGraph[mcpid][pixid]->GetYaxis()->SetTitle("Sigma");
				SigmaGraph[mcpid][pixid]->SetMarkerSize(0.9);
				SigmaGraph[mcpid][pixid]->SetMarkerStyle(22);
				SigmaGraph[mcpid][pixid]->SetMarkerColor(2);
				SigmaGraph[mcpid][pixid]->SetLineColor(2);
				//
				SigmaGraphLED[mcpid][pixid]->GetXaxis()->SetTitle("Stat");
				SigmaGraphLED[mcpid][pixid]->GetYaxis()->SetTitle("Sigma");
				SigmaGraphLED[mcpid][pixid]->SetMarkerSize(0.9);
				SigmaGraphLED[mcpid][pixid]->SetMarkerStyle(21);
				SigmaGraphLED[mcpid][pixid]->SetMarkerColor(4);
				SigmaGraphLED[mcpid][pixid]->SetLineColor(4);
				//
				/*
				TCanvas *canvasgraph = new TCanvas("canvasgraph","A Simple Graph laser",200,10,700,500);
				SigmaGraph[mcpid][pixid]->Draw("APL");
				canvasgraph->Update();
				canvasgraph->WaitPrimitive();
				//
				TCanvas *canvasgraphLED = new TCanvas("canvasgraphLED","A Simple Graph LED",200,10,700,500);
				SigmaGraphLED[mcpid][pixid]->Draw("APL");
				canvasgraphLED->Update();
				canvasgraphLED->WaitPrimitive();
				*/
				mg[mcpid][pixid]->Add(SigmaGraph[mcpid][pixid],"LP");
				mg[mcpid][pixid]->Add(SigmaGraphLED[mcpid][pixid],"LP");
				////////////////////////////////
				///////// Draw Multi Graph  ////
				////////////////////////////////
				if (Canvmg){
					TCanvas *canvasMultiGraph = new TCanvas("canvasMultiGraph"," Multi Graph",200,10,700,500);
					mg[mcpid][pixid]->Draw("a");
					canvasMultiGraph->BuildLegend();
					canvasMultiGraph->Update();
					canvasMultiGraph->WaitPrimitive();
				}
				
			} // end number of entries loop
		}// end of sensor Id loop
	}// files loop
	//TCanvas* ok = new TCanvas("ok","ok",800,600);
	//randfill->Draw("colz");
}// end of the function


