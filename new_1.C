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

void new_1(TString infile="../build/20m_20_laser_25_1_4.root", Int_t ID= -99, bool laser=true){ 
	gSystem->Load("../build/libGlxDirc.so");
	gStyle->SetOptStat(1);
	//t=0;  for i in ../build/1m*.root; do  root -l -q loadlib.C test_laser.C+"(\"$i\",$((++t)),\"out_${t}_.root\")" ;done
	std::cout <<"######################################################################= "<<ID<<std::endl;
	TF1* fit1_histo ;
	fit1_histo = new TF1("fit1_histo",conv,0,15,3); //0, 20 the range , 2 number of parameters on the function 
	fit1_histo->SetNpx(500);
	fit1_histo->SetParLimits(2,2,3);
	fit1_histo->SetParNames ("Constant","Mean","Width");

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
	/////////////////////////////////////////////
	///////// histo initialization //////////////
	/////////////////////////////////////////////
	TGraph * SigmaGraph[nmcp][npix];
	TH1F * hPTime[nmcp][npix], * hPTimeSigma[nmcp][npix];
	TH1F * hTime = new TH1F("time",";time, [ns];entries, [#]",500,0,150);
	TH2F * hPMToc[nmcp];
	// for all mcp 
	for(Int_t m=0; m<nmcp; m++){
		// loop over pix
		for(Int_t p=0; p<npix; p++){
			//histogram for each pixl
			hPTime[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),500,-25,25);//50
			if(laser)hPTimeSigma[m][p]  = new TH1F(Form("hPTimeSigma%d",m*100+p),Form("mcp %d, pixel %d",m, p),200,-0.1,0.1); // laser
			else hPTimeSigma[m][p]  = new TH1F(Form("hPTimeSigma%d",m*100+p),Form("mcp %d, pixel %d",m, p),200,-0.1,0.1);
			axisTime800x500(hPTime[m][p]);
			//gStyle->SetOptTitle(0);
			hPTime[m][p]->SetStats(1);
			hPTime[m][p]->SetLineColor(1);
			hPTime[m][p]->SetLineWidth(3);
			// graph declaration
			SigmaGraph[m][p]= new TGraph();
			SigmaGraph[m][p]->SetTitle(Form("mcp %d, pixel %d",m, p));

		}
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
			if(entries > 2500){
				mcpid =( s/100 ) -102;
				pixid = s%100 -1;
				m =pixid%8;
				n =7-pixid/8;
				if(mcpid != 54) continue;
				//cout<< "#####################################################entries " << entries<< "	s	"<<s<<endl;
				/////////////////////////////////////////////
				///////// 	fill pos histo	   //////////////
				/////////////////////////////////////////////
				Int_t npoint=0;
				for (Int_t step=0; step< 50; step++) {
					entries= step * 20+10;
					for (Int_t d=0; d< entries; d++) {
						
						//cout<< "entries " << entries<<endl;
						//TVector3 dir = node->GetEntry(d);
						//cout<< "XX " <<dir.X() <<endl;
						Double_t evtime = node->GetTime(d);
						Double_t smear= 0.120;  // PMT smearing 120 ps
						evtime+=fRand1.Gaus(0,smear);
						//adding light emission profile of the LED
						
						if(!laser){
						do{
							r = gRandom->Rndm();
							//cout<< "r " <<r <<endl;
							tr = (gRandom->Rndm())* 7.0;
							val= gr->Eval(tr);
							if(r<=val) randfill->Fill(tr, r);
							if(true){ 
								if(r<=val){
									noise = 0;
									pro= gRandom->Uniform(0,200000); 
									if(pro == 333) {noise =gRandom->Uniform(0,20); hPTime[mcpid][pixid]->Fill(noise); }
									//cout<< "pro= " <<pro <<"	noise= "<<noise<<endl;
									hPTime[mcpid][pixid]->Fill(evtime+tr-3.5); // add tr-3.5

								}
							} // add tr for LED, remove tr for laser
							if(false){ // switsh ON time offset between sources
								if(k==0) hPTime[mcpid][pixid]->Fill(evtime+tr+0.0);
								if(k==1) hPTime[mcpid][pixid]->Fill(evtime+tr+5.0);
								if(k==2) hPTime[mcpid][pixid]->Fill(evtime+tr+10.0);
							}
						}while(r>val);
					}
						else{
						hPTime[mcpid][pixid]->Fill(evtime); // laser
					}

						
					}
					//TCanvas* canvGaus = new TCanvas("canvGaus","canvGaus",800,600);
					if(hPTime[mcpid][pixid]->GetEntries() < 5 ) continue;
					cout<< "number hPTime= " <<hPTime[mcpid][pixid]->GetEntries()<<endl;
					
					Double_t mean_g_fit;
					if (laser){
					// Laser stuff
					hPTime[mcpid][pixid]->Fit("gaus");
					 mean_g_fit = hPTime[mcpid][pixid]->GetFunction("gaus")->GetParameter(1);
					
					}else
					{
					// LED Stuff
					fit1_histo->SetParameters(1,hPTime[mcpid][pixid]->GetMean(), 2.5);
					hPTime[mcpid][pixid]->Fit(fit1_histo);
					 mean_g_fit = fit1_histo->GetParameter(1);
				
					}
					
					//cout<< "################ mean_g_fit= " <<mean_g_fit<<endl;
					//hPTime[mcpid][pixid]->Draw();
					//canvGaus->Update();
					//canvGaus->WaitPrimitive();
					
					
					for (Int_t d=0; d< entries; d++) {
						Double_t evtime = node->GetTime(d);
						hPTimeSigma[mcpid][pixid]->Fill(mean_g_fit-evtime);
						
					}
					//TCanvas* test = new TCanvas("test","test",800,600);
					cout<< "number hPTimeSigma= " <<hPTimeSigma[mcpid][pixid]->GetEntries()<<endl;
					cout<< "integral hPTimeSigma= " <<hPTimeSigma[mcpid][pixid]->Integral(-5,5)<<endl;

					
					
					hPTimeSigma[mcpid][pixid]->Fit("gaus");
					hPTimeSigma[mcpid][pixid]->Draw();
					//test->Update();
					//test->WaitPrimitive();
					if(hPTimeSigma[mcpid][pixid]->Integral(5,195)== 0) continue;
					sigma = hPTimeSigma[mcpid][pixid]->GetFunction("gaus")->GetParameter(2);
					cout<< "npoint = " <<npoint<<endl;
					SigmaGraph[mcpid][pixid]->SetPoint(npoint,entries,sigma);
					hPTimeSigma[mcpid][pixid]->Reset();
					npoint++;	              
				}// end of steps loop
				
			
				
				TCanvas *canvasgraph = new TCanvas("canvasgraph","A Simple Graph Example",200,10,700,500);
				SigmaGraph[mcpid][pixid]->GetXaxis()->SetTitle("Stat");
				SigmaGraph[mcpid][pixid]->GetYaxis()->SetTitle("Sigma");
				SigmaGraph[mcpid][pixid]->SetMarkerSize(0.9);
				SigmaGraph[mcpid][pixid]->SetMarkerStyle(7);	
				SigmaGraph[mcpid][pixid]->Draw("APL");
				canvasgraph->Update();
				canvasgraph->WaitPrimitive();
				
			} // end number of entries loop
		}// end of sensor Id loop
	}// files loop
	//TCanvas* ok = new TCanvas("ok","ok",800,600);
	//randfill->Draw("colz");
}// end of the function


