// Author Yeonju Go
// last modification : 2015/01/27 
// 
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "TLatex.h"
#include "stdio.h"
#include "../yjUtility.h"


void noiseLevelCheck(){

    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);
    TH1::SetDefaultSumw2();

    const char* fname = "/d3/scratch/goyeonju/files/forest/centrality/HiForestAOD_data_PbPb_withTowers_MinimumBias1_TriggerFilters_run262921_v1.root";
    TFile *f = TFile::Open(fname);
    TTree *t = (TTree*) f->Get("rechitanalyzer/tower");
    TTree *thlt = (TTree*) f->Get("hltanalysis/HltTree");
    t->AddFriend(thlt);

    const int nTrig = 3;
    TCut trigCut[nTrig] = {"","HLT_HIL1MinimumBiasHF1AND_v1","HLT_HIL1Tech5_BPTX_PlusOnly_v1 || HLT_HIL1Tech6_BPTX_MinusOnly_v1 || HLT_HIL1Tech7_NoBPTX_v1"};
    TCut hfCut = "abs(eta)>3 && abs(eta)<5";

    double eThr = 20.0;
    TH1D* h[nTrig];
    for(int i=0;i<nTrig;i++){
        h[i] = new TH1D(Form("h%d",i), ";e;Events", 100, 0, eThr);
        t->Draw(Form("e>>h%d",i), trigCut[i] && hfCut);
    }
 
    TLegend* l1 = new TLegend(0.6,0.6,0.9,0.9);
    TCanvas* c = new TCanvas();
    
    for(int i=0;i<nTrig;i++){
        SetHistColor(h[i],i+1);
        l1->AddEntry(h[i],trigCut[i]->GetTitle());
        if(i==0) h[0]->Draw("hist");
        else h[i]->Draw("same hist");
    }
    l1->Draw();
    c->SaveAs("figures/noiseLevel_e.png");

    TCanvas* c_logy = (TCanvas*) c->Copy();
    c_logy->SetLogy();
    c_logy->SaveAs("figures/noiseLevel_e_logy.png");
} 
