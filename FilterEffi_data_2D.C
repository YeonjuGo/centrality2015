// Author Yeonju Go
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
#include "../HiForestAnalysis/hiForest.h"

void FilterEffi_data_2D()
{
    /*    const int onlineFilter = 0;//HLT_HIMinBiasHfOrBSC_v1
          const int collCut = 0; //pcollisionEventSelection
          const int vertexCut = 0; //pprimaryVertexFilter
          const int pixelCut = 0; //phltPixelClusterShapeFilter
          const int hfCoincCut = 0; //phfCoincFilter3
          */
    const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
    const TCut eventCut = runCut && lumiCut;
    const int Ncut = 5;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);

    TFile *dataf = new TFile("files/HiForest_pbpb_data_9xx_merged.root");
    TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
    TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
    TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
    datat_evt -> AddFriend(datat_hlt);
    datat_evt -> AddFriend(datat_skim);


    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);

    TLine* t1 = new TLine(0,1,1000,1);
    t1->SetLineWidth(1);
    t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
    t1->SetLineColor(1); // 2 is red
/*
    //==========================================================
    // HF pixel hits plus vs. HF pixel hits minus 2D plot
    //==========================================================


    TH2D *HFPixelPlMa_data[5]; 
    for(int i=0; i<Ncut; i++){
        HFPixelPlMa_data[i] = new TH2D(Form("HFPixelPlMa_data%d",i), ";hiHFhitMinus; hiHFhitPlus}",800,0,10000, 800,0,10000);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus>>+ HFPixelPlMa_data0",eventCut && "HLT_HIMinBiasHfOrBSC_v1==1");
    HFPixelPlMa_data[0] = (TH2D*)gDirectory->Get("HFPixelPlMa_data0");
    datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus>>+ HFPixelPlMa_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    HFPixelPlMa_data[1] = (TH2D*)gDirectory->Get("HFPixelPlMa_data1");
    datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus>>+ HFPixelPlMa_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    HFPixelPlMa_data[2] = (TH2D*)gDirectory->Get("HFPixelPlMa_data2");
    datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus>>+ HFPixelPlMa_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    HFPixelPlMa_data[3] = (TH2D*)gDirectory->Get("HFPixelPlMa_data3");
    datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus>>+ HFPixelPlMa_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    HFPixelPlMa_data[4] = (TH2D*)gDirectory->Get("HFPixelPlMa_data4");

    TCanvas *c_HFPixelPlMa = new TCanvas("c_HFPixelPlMa", "c_HFPixelPlMa", 900,600);
    c_HFPixelPlMa -> Divide(3,2);
    for(int i=0;i<Ncut;i++){
        c_HFPixelPlMa -> cd(i+1);
        HFPixelPlMa_data[i] -> Draw("colz");
    }

*/
    //==========================================================
    // HF pixel hits plus vs. HF pixel hits minus 2D plot
    //==========================================================

    TH2D *hiHF_hiHFhit_data[5]; 
    for(int i=0; i<Ncut; i++){
        hiHF_hiHFhit_data[i] = new TH2D(Form("hiHF_hiHFhit_data%d",i), ";hiHF; hiHFhit",200,0,4500, 200,0,120000);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiHFhit:hiHF>>+ hiHF_hiHFhit_data0",eventCut && "HLT_HIMinBiasHfOrBSC_v1==1");
    hiHF_hiHFhit_data[0] = (TH2D*)gDirectory->Get("hiHF_hiHFhit_data0");
    datat_evt -> Draw("hiHFhit:hiHF>>+ hiHF_hiHFhit_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiHF_hiHFhit_data[1] = (TH2D*)gDirectory->Get("hiHF_hiHFhit_data1");
    datat_evt -> Draw("hiHFhit:hiHF>>+ hiHF_hiHFhit_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiHF_hiHFhit_data[2] = (TH2D*)gDirectory->Get("hiHF_hiHFhit_data2");
    datat_evt -> Draw("hiHFhit:hiHF>>+ hiHF_hiHFhit_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiHF_hiHFhit_data[3] = (TH2D*)gDirectory->Get("hiHF_hiHFhit_data3");
    datat_evt -> Draw("hiHFhit:hiHF>>+ hiHF_hiHFhit_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiHF_hiHFhit_data[4] = (TH2D*)gDirectory->Get("hiHF_hiHFhit_data4");

    TCanvas *c_hiHF_hiHFhit = new TCanvas("c_hiHF_hiHFhit", "c_hiHF_hiHFhit", 900,600);
    c_hiHF_hiHFhit -> Divide(3,2);
    for(int i=0;i<Ncut;i++){
        c_hiHF_hiHFhit -> cd(i+1);
        hiHF_hiHFhit_data[i] -> Draw("colz");
    }
    c_temp -> cd();
    hiHF_hiHFhit_data[4] -> Draw("colz");

 
    //==========================================================
    // hiHF vs. hiNpix 2D plot
    //==========================================================

    TH2D *hiHF_hiNpix_data[5]; 
    for(int i=0; i<Ncut; i++){
        hiHF_hiNpix_data[i] = new TH2D(Form("hiHF_hiNpix_data%d",i), ";hiHF; hiNpix",100,0,40, 100,0,400);
        hiHF_hiNpix_data[i] -> SetLabelSize(0.03,"Y"); 
        hiHF_hiNpix_data[i] -> SetLabelOffset(-0.01,"Z"); 
    }

    c_temp -> cd();
    datat_evt -> Draw("hiNpix:hiHF>>+ hiHF_hiNpix_data0",eventCut && "HLT_HIMinBiasHfOrBSC_v1==1");
    hiHF_hiNpix_data[0] = (TH2D*)gDirectory->Get("hiHF_hiNpix_data0");
    datat_evt -> Draw("hiNpix:hiHF>>+ hiHF_hiNpix_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==0");
    hiHF_hiNpix_data[1] = (TH2D*)gDirectory->Get("hiHF_hiNpix_data1");
    datat_evt -> Draw("hiNpix:hiHF>>+ hiHF_hiNpix_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==0");
    hiHF_hiNpix_data[2] = (TH2D*)gDirectory->Get("hiHF_hiNpix_data2");
    datat_evt -> Draw("hiNpix:hiHF>>+ hiHF_hiNpix_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==0");
    hiHF_hiNpix_data[3] = (TH2D*)gDirectory->Get("hiHF_hiNpix_data3");
    datat_evt -> Draw("hiNpix:hiHF>>+ hiHF_hiNpix_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==0");
    hiHF_hiNpix_data[4] = (TH2D*)gDirectory->Get("hiHF_hiNpix_data4");

    TCanvas *c_hiHF_hiNpix = new TCanvas("c_hiHF_hiNpix", "c_hiHF_hiNpix", 1080,720);
    c_hiHF_hiNpix -> Divide(3,2);
    for(int i=0;i<Ncut;i++){
        c_hiHF_hiNpix -> cd(i+1);
        hiHF_hiNpix_data[i] -> Draw("colz");
    }
    c_temp -> cd();
    hiHF_hiNpix_data[4] -> Draw("colz");
} 
