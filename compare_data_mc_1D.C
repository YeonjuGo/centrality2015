// comparison centrality objects between DATA and MC 
// Author : Yeonju Go
// 24 Nov 2015

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
#include <TCut.h>
//canvas, legend, latex ... //cosmetic
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>
//private setup
#include "../yjUtility.h"

void compareTwo(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, TCut cut1="(1)", TCut cut2="(1)", const string cap = "");
void compare_data_mc_1D(){
    const char* fname_data="/home/goyeonju/CMS/Files/centrality/EventTree_PbPb_data_HIMinimumBias2_run262921_15Feb.root";
	const char* fname_mc="/home/goyeonju/CMS/Files/centrality/HiForestAOD_HYDJET_Drum5Nv75_5020GeV_cmssw758patch3_genPtThr500GeV_merged.root";

    TFile* f1 = TFile::Open(fname_data);
    TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
    TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
    TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
    t1->AddFriend(t1_skim);
    t1->AddFriend(t1_hlt);
    TFile* f2 = TFile::Open(fname_mc);
    TTree* t2 = (TTree*) f2 -> Get("hiEvtAnalyzer/HiTree");
    TTree* t2_skim = (TTree*) f2 -> Get("skimanalysis/HltTree");
    TTree* t2_hlt = (TTree*) f2 -> Get("hltanalysis/HltTree");
    t2->AddFriend(t2_skim);
    t2->AddFriend(t2_hlt);

    /*
    //from DATA // this sample has no track info. 
    double hiHFMax = 5000;
    double hiBinMax = 200;
    double hiHFhitMax = 160000; 
    double hiETMax = 1600;
    double hiEEMax = 2200;
    double hiEBMax = 4000;
    double hiNpixMax = 0;
    double hiNpixelTracksMax = 40000;
    double hiNtracksMax = 3500;
    double hiNtracksPtCutMax = 2200;
    double hiNtracksEtaCutMax = 3000;
    double hiNtracksEtaPtCutMax = 1800;
    double NPixClus= 1300;
*/
   //total
    double hiHFMax = 9000;
    double hiBinMax = 200;
    double hiHFhitMax = 250000; 
    double hiETMax = 1800;
    double hiEEMax = 3500;
    double hiEBMax = 4000;
    double hiNpixMax = 50000;
    double hiNpixelTracksMax = 5000;
    double hiNtracksMax = 3500;
    double hiZDCMax = 1;
 
    int nBin = 50;
    const char* trigCut = "HLT_HIL1MinimumBiasHF1AND_v1 && pcollisionEventSelection";
    const char* mcCut = "pcollisionEventSelection";
    const string cap = "HIMinimumBias2_run262921";
    compareTwo(t1, t2, "hiHF",nBin,0,hiHFMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiBin",nBin,0,hiBinMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiHFhit",nBin,0,hiHFhitMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiET",nBin,0,hiETMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiEB",nBin,0,hiEBMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiEE",nBin,0,hiEEMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiNpix",nBin,0,hiNpixMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiNtracks",nBin,0,hiNtracksMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiNpixelTracks",nBin,0,hiNpixelTracksMax,trigCut,mcCut, cap);
    compareTwo(t1, t2, "hiZDC",nBin,0,hiZDCMax,trigCut,mcCut, cap);

} // main function

void compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, TCut cut1, TCut cut2, const string cap)  {
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0); 
    static int j = 0;
    TCanvas* c=  new TCanvas(Form("c_%s_%d",var.Data(),j),"", 400,800);
	c->Divide(1,2);
	c->cd(1);
	gPad->SetLogy();
	TH1D* h1 = new TH1D(Form("h1_%s_%d",var.Data(),j), Form(";%s;",var.Data()), nBins,xMin,xMax);
	TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%d",var.Data(),j));
	h1->Sumw2();
	h2->Sumw2();
	t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), cut1);
	t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), cut2);	
	h1->Scale( 1. / t1->GetEntries(cut1));
	h2->Scale( 1. / t2->GetEntries(cut2));
    SetHistColor(h1,2);
    SetHistColor(h2,1);
    h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.8);
    h1->SetNdivisions(505);
    h2->SetNdivisions(505);
    TLegend* l1 = new TLegend(0.5,0.6,0.90,0.8);
    legStyle(l1);
    l1->AddEntry(h1, "Data", "p");
    l1->AddEntry(h2, "Hydjet Drum5Nv75", "l");

    double range = cleverRange(h1,h2,1.5,1.e-4);
	h1->DrawCopy("L");
	h2->DrawCopy("hist same");
    l1->Draw("same");
    
    c->cd(2);
	h1->Divide(h2);
	h1->SetYTitle("DATA / HYDJET(Drum5Nv75)");
	h1->GetYaxis()->SetRangeUser(0.0,2.0);
    SetHistColor(h1,2);
    h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.8);
	h1->DrawCopy("");
	//h2->DrawCopy("le1");
	jumSun(xMin,1,xMax,1);
    drawText(cap.data(),0.2,0.2+0.7);
    //drawText(cut1.GetTitle(),0.2,0.2);
    drawText("HLT_HIL1MinimumBiasHF1AND_v1",0.2,0.2+0.06);
    drawText("pcollisionEventSelection",0.2,0.2);
	c-> SaveAs(Form("pdf/compare_%s_%s.pdf",var.Data(),cap.data()));
    j++;
}
