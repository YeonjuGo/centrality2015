// to draw the distribution of centrality variables
// with the different event selection filters in one canvas.
// should modify evtfilter definition part
//
// Author Yeonju Go
// 23 Nov 2015
//
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
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
#include "stdio.h"
#include <iostream>
#include "yjUtility.h"

//////////// modify here as you want!! ////////////////
const int nfilter = 7;
const char* evtfilter[] = {"","pBeamScrapingFilter","pPAprimaryVertexFilter","pHBHENoiseFilterResultProducer","phfCoincFilter","phfCoincFilter3","PAcollisionEventSelection"};
const double dy= 0.5;

void Get1DEffPlots(TTree* t_evt=0, TString v1="hiHF",int xbin=200, double xmin=0, double xmax=4500, TCut cut="",const char* cap="", bool isPassed=0,bool isAOD=0);

void filterDep_1D_Eff(const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163-262172_1.4M.root",
        const char* type="Run262163-262172", 
        TCut trig = "HLT_L1MinimumBiasHF1OR_part1_v1",
        bool isAOD=0)
{
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);
    SetHistTitleStyle(0.06,0.04);
    SetyjPadStyle();

    TFile *fin = TFile::Open(fname);
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
    t_evt -> AddFriend(t_hlt);
    t_evt -> AddFriend(t_skim);

    double Nevt_t = t_evt -> GetEntries(trig);
    cout << "Trigger is '" << trig.GetTitle() << "' and # of events = " << Nevt_t << endl;

    Get1DEffPlots(t_evt, "hiHF",100,0,5000,trig,type,1,isAOD);
    Get1DEffPlots(t_evt, "hiHFhit",100,0,4000,trig,type,1,isAOD);
    Get1DEffPlots(t_evt, "hiNpix",100,0,10000,trig,type,1,isAOD);
   // Get1DEffPlots(t_evt, "hiBin",105,0,210,trig,type,1,isAOD);
   // Get1DEffPlots(t_evt, "hiZDC",100,0,50000,trig,type,1,isAOD);
   // Get1DEffPlots(t_evt, "hiET",100,0,2000,trig,type,1,isAOD);
   // Get1DEffPlots(t_evt, "hiEE",100,0,4000,trig,type,1,isAOD);
   // Get1DEffPlots(t_evt, "hiEB",100,0,5000,trig,type,1,isAOD);

    Get1DEffPlots(t_evt, "hiHF",100,0,100,trig,Form("%s_shortRange",type),1,isAOD);
    Get1DEffPlots(t_evt, "hiHFhit",100,0,500,trig,Form("%s_shortRange",type),1,isAOD);
    Get1DEffPlots(t_evt, "hiNpix",100,0,300,trig,Form("%s_shortRange",type),1,isAOD);
   // Get1DEffPlots(t_evt, "hiBin",105,100,210,trig,Form("%s_shortRange",type),1,isAOD);
   // Get1DEffPlots(t_evt, "hiZDC",100,0,5000,trig,Form("%s_shortRange",type),1,isAOD);
   // Get1DEffPlots(t_evt, "hiET",20,0,20,trig,Form("%s_shortRange",type),1,isAOD);
   // Get1DEffPlots(t_evt, "hiEE",20,0,20,trig,Form("%s_shortRange",type),1,isAOD);
   // Get1DEffPlots(t_evt, "hiEB",80,0,80,trig,Form("%s_shortRange",type),1,isAOD);
}

void Get1DEffPlots(TTree* t_evt, TString v1, int xbin, double xmin, double xmax, TCut cut,const char* cap, bool isPassed, bool isAOD)
{
    TCanvas *c_tot = new TCanvas(Form("c_tot_%s_%s",v1.Data(),cap), "c_tot", 300,600);
    c_tot->Divide(1,2);

    TCut totcut[nfilter];
    for(int i=0; i<nfilter;i++){
        if(i==0) totcut[i] = cut; 
        else totcut[i] = cut && Form("%s==%d",evtfilter[i],(int)isPassed);
    }

    TCanvas *c_temp= new TCanvas(Form("c_temp_%s",cap), "", 300,300);
    TH1D *h1D[nfilter];
    TH1D *h1D_eff[nfilter];
    for(int i=0; i<nfilter; i++){
        h1D[i] = new TH1D(Form("h1D_%d",i), Form(";%s;Events", v1.Data()), xbin, xmin,xmax );
        if(i!=0) h1D[i]->SetMarkerStyle(24+i); 
        h1D[i]->SetMarkerSize(0.8);
        h1D_eff[i] = (TH1D*)h1D[i]->Clone(Form("h1D_eff_%d",i));
        h1D_eff[i] -> SetTitle(Form(";%s;Filter Rate",v1.Data()));
        SetHistColor(h1D[i],col[i]);
        SetHistColor(h1D_eff[i],col[i]);
    }
    for(int i=0; i<nfilter; i++){
        t_evt->Draw(Form("%s>>+%s",v1.Data(), h1D[i]->GetName()), totcut[i]);
        h1D[i]=(TH1D*)gDirectory->Get(h1D[i]->GetName());
        if(i!=0) h1D_eff[i] -> Divide(h1D[i],h1D[0],1,1,"B");
    }
    TLegend* l1 = new TLegend(0.35, 0.7, 0.9, 0.95, Form("%s",cap));
    legStyle(l1);
    for(int i=0; i<nfilter;i++){
        if(i==0) l1->AddEntry(h1D[0], "NO FILTER", "l"); 
        else l1->AddEntry(h1D[i], evtfilter[i]); 
    }

    for(int i=0;i<nfilter;i++){
        c_tot->cd(1);
        if(i==0) h1D[i] -> DrawCopy("hist");
        else h1D[i] -> DrawCopy("ep same");
        l1 -> Draw();
        if(i!=0){
            c_tot->cd(2);
            h1D_eff[i]->GetYaxis()-> SetRangeUser(h1D_eff[nfilter-1]->GetMinimum()*0.9,1.0);
            if(i==1) {
                h1D_eff[i] -> DrawCopy("ep"); 
                jumSun(xmin,1,xmax,1);
            }
            else h1D_eff[i] -> DrawCopy("ep same"); 
        }
    }
    c_tot->SaveAs(Form("pdf/filterDep_1D_Eff_%s_%s.pdf",v1.Data(),cap));
/*
    TFile* outf = new TFile(Form("histfiles/filterEff_1D_%s_%s.root",v1.Data(),cap), "RECREATE");
    outf->cd();
    for(int i=0; i<nfilter; i++){
        h1D[i]->Write();
        h1D_eff[i]->Write();
    }
    c_tot->Write();
    outf->Close();
*/
    for(int i=0; i<nfilter; i++){
        delete h1D[i];
        delete h1D_eff[i];
    }
    delete c_temp;
}
