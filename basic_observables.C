// Author Yeonju Go
#include "TFile.h"
#include "TTree.h"
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
//#include "TPad.h"
#include "stdio.h"
#include "../HiForestAnalysis/hiForest.h"

void basic_observables()
{
    const int onlineFilter = 0;//HLT_HIMinBiasHfOrBSC_v1
    const int collCut = 0; //pcollisionEventSelection
    const int vertexCut = 0; //pprimaryVertexFilter
    const int pixelCut = 0; //phltPixelClusterShapeFilter
    const int hfCoincCut = 0; //phfCoincFilter3
 
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);
    
    TFile *dataf = new TFile("PbPb_Centrality_data.root");
    TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
    TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
    TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
    datat_evt -> AddFriend(datat_hlt);
    datat_evt -> AddFriend(datat_skim);
    
    TFile *mcf = new TFile("PbPb_Centrality_mc.root");
    TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
    TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
    TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
    mct_evt -> AddFriend(mct_hlt);
    mct_evt -> AddFriend(mct_skim);
    cout << mcf << endl;

//
//HF sum!!
//

    TH1D *HFsum_data = new TH1D("HFsum_data", ";#Sigma E_{T}^{HF};Events",100,0,5000);
    TH1D *HFsum_mc = new TH1D("HFsum_mc", ";#Sigma E_{T}^{HF};Events",100,0,5000);
    
    if(onlineFilter==1){
        datat_evt -> Draw("hiHF >>+ HFsum_data", "1.0*(HLT_HIMinBiasHfOrBSC_v1==1)");
        HFsum_data = (TH1D*)gDirectory->Get("HFsum_data");
        mct_evt -> Draw("hiHF >>+ HFsum_mc", "1.0*(HLT_HIMinBiasHfOrBSC_v1==1)");
        HFsum_mc = (TH1D*)gDirectory->Get("HFsum_mc");
    } else {
        datat_evt -> Draw("hiHF >>+ HFsum_data");
        HFsum_data = (TH1D*)gDirectory->Get("HFsum_data");
        mct_evt -> Draw("hiHF >>+ HFsum_mc");
        HFsum_mc = (TH1D*)gDirectory->Get("HFsum_mc");
    }
    HFsum_data -> SetMarkerStyle(20);
    HFsum_mc -> SetMarkerStyle(20);
    HFsum_data -> SetMarkerColor(2);
 
    TLegend* l1 = new TLegend(0.6, 0.65, 0.85, 0.80);
//    l1 -> AddEntry((TObject*)0, "min"); 
    l1 -> AddEntry(HFsum_data, "PbPb DATA");
    l1 -> AddEntry(HFsum_mc, "PbPb MC");
    l1 -> Draw();

/*
    TCanvas *c_HFsum = new TCanvas("c_HFsum", "c_HFsum", 400,400);
    c_HFsum -> SetLogy();
    HFsum_data -> Draw("e1p");
    HFsum_mc -> Draw("same&&e1p");
*/


//
// HF pixel hits
//
/*
    TH1D *HFpixelhit_data = new TH1D("HFpixelhit_data", ";#Sigma E_{T}^{HF pixel hit};Events",100,0,5000);
    TH1D *HFpixelhit_mc = new TH1D("HFpixelhit_mc", ";#Sigma E_{T}^{HF pixel hit};Events",100,0,5000);
    
    if(onlineFilter==1){
        datat_evt -> Draw("hiHF >>+ HFpixelhit_data", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFpixelhit_data = (TH1D*)gDirectory->Get("HFpixelhit_data");
        mct_evt -> Draw("hiHF >>+ HFpixelhit_mc", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFpixelhit_mc = (TH1D*)gDirectory->Get("HFpixelhit_mc");
    } else {
        datat_evt -> Draw("hiHF >>+ HFpixelhit_data");
        HFpixelhit_data = (TH1D*)gDirectory->Get("HFpixelhit_data");
        mct_evt -> Draw("hiHF >>+ HFpixelhit_mc");
        HFpixelhit_mc = (TH1D*)gDirectory->Get("HFpixelhit_mc");
    }

    TCanvas *c_HFpixelhit = new TCanvas("c_HFpixelhit", "c_HFpixelhit", 400,400);
    c_HFpixelhit -> SetLogy();
    HFpixelhit_data -> SetMarkerStyle(20);
    HFpixelhit_mc -> SetMarkerStyle(20);
    HFpixelhit_data -> SetMarkerColor(2);
    HFpixelhit_data -> Draw("e1p");
    HFpixelhit_mc -> Draw("same&&e1p");
    l1 -> Draw();
*/
    
//
// HF pixel hits vs. HFsum 2D plot
//
/*
    //gPad-> SetLogy(1);
    TH2D *HFPixelSum_data = new TH2D("HFPixelSum_data", ";#Sigma E_{T}^{HF}; #Sigma E_{HF pixel hits}",100,0,5000,100,0,160000);
    TH2D *HFPixelSum_mc = new TH2D("HFPixelSum_mc", ";#Sigma E_{T}^{HF}; #Sigma E_{HF pixel hits}",100,0,5000,100,0,160000);
    
    if(onlineFilter==1){
        datat_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_data", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFPixelSum_data = (TH2D*)gDirectory->Get("HFPixelSum_data");
        mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFPixelSum_mc = (TH2D*)gDirectory->Get("HFPixelSum_mc");
    } else {
        datat_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_data");
        HFPixelSum_data = (TH2D*)gDirectory->Get("HFPixelSum_data");
        mct_evt -> Draw("hiHFhit:hiHF >>+ HFPixelSum_mc");
        HFPixelSum_mc = (TH2D*)gDirectory->Get("HFPixelSum_mc");
    }

    TCanvas *c_HFPixelSum = new TCanvas("c_HFPixelSum", "c_HFPixelSum", 800,400);
    c_HFPixelSum -> Divide(2,1);
    c_HFPixelSum -> cd(1);
    HFPixelSum_data -> Draw("colz");
    c_HFPixelSum -> cd(2);
    HFPixelSum_mc -> Draw("colz");
*/

//
// HF E plus vs. HF E minus 2D plot
//
/*
    TH2D *HFPlMa_data = new TH2D("HFPlMa_data", "; #Sigma E_{HF minus  hits}; #Sigma E_{HF plus  hits}",200,0,60000,200,0,60000);
    TH2D *HFPlMa_mc = new TH2D("HFPlMa_mc", "; #Sigma E_{HF minus  hits}; #Sigma E_{HF plus  hits}",200,0,60000,200,0,60000);
    
    if(onlineFilter==1){
        datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_data", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFPlMa_data = (TH2D*)gDirectory->Get("HFPlMa_data");
        mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_mc", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFPlMa_mc = (TH2D*)gDirectory->Get("HFPlMa_mc");
    } else {
        datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_data");
        HFPlMa_data = (TH2D*)gDirectory->Get("HFPlMa_data");
        mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPlMa_mc");
        HFPlMa_mc = (TH2D*)gDirectory->Get("HFPlMa_mc");
    }

    TCanvas *c_HFPlMa = new TCanvas("c_HFPlMa", "c_HFPlMa", 800,400);
    c_HFPlMa -> SetLogy();
    c_HFPlMa -> Divide(2,1);
    c_HFPlMa -> cd(1);
    HFPlMa_data -> Draw("colz");
    c_HFPlMa -> cd(2);
    HFPlMa_mc -> Draw("colz");
*/

//
// HF pixel hits plus vs. HF pixel hits minus 2D plot
//
    TH2D *HFPixelPlMa_data = new TH2D("HFPixelPlMa_data", "; #Sigma E_{HF minus pixel hits}; #Sigma E_{HF plus pixel hits}",200,0,60000,600,0,60000);
    TH2D *HFPixelPlMa_mc = new TH2D("HFPixelPlMa_mc", "; #Sigma E_{HF minus pixel hits}; #Sigma E_{HF plus pixel hits}",200,0,60000,600,0,60000);
    TH2D *ghost = (TH2D*) HFPixelPlMa_data->Clone("ghost"); 
    
    if(onlineFilter==1){
        datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_data", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFPixelPlMa_data = (TH2D*)gDirectory->Get("HFPixelPlMa_data");
        mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_mc", "HLT_HIMinBiasHfOrBSC_v1==1");
        HFPixelPlMa_mc = (TH2D*)gDirectory->Get("HFPixelPlMa_mc");
    } else {
        datat_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_data");
        HFPixelPlMa_data = (TH2D*)gDirectory->Get("HFPixelPlMa_data");
        mct_evt -> Draw("hiHFhitPlus:hiHFhitMinus >>+ HFPixelPlMa_mc");
        HFPixelPlMa_mc = (TH2D*)gDirectory->Get("HFPixelPlMa_mc");
    }

    TCanvas *c_HFPixelPlMa = new TCanvas("c_HFPixelPlMa", "c_HFPixelPlMa", 800,400);
    c_HFPixelPlMa -> Divide(2,1);
    c_HFPixelPlMa -> SetLogy(); 
    c_HFPixelPlMa -> cd(1);
    HFPixelPlMa_data -> Draw("colz");
    c_HFPixelPlMa -> cd(2);
    ghost -> Draw();
    HFPixelPlMa_mc -> Draw("same && colz");







    /*
 // TH1D *pt0 = new TH1D("pt0","pt0;p_{T}^{#gamma}",50,0,150);
   TH1D *pt0 = new TH1D("pt0","pt0;p_{T}^{#gamma}",50,0,150);
   TH1D *pt1 = (TH1D*)pt0->Clone("pt1");
   TH1D *eta0 = new TH1D("eta0","eta0;#eta^{#gamma}",50, -5, 5);
   TH1D *eta1 = (TH1D*)eta0->Clone("eta1");
   TH1D *phi0 = new TH1D("phi0","phi0;#phi^{#gamma}",50, -TMath::Pi(), TMath::Pi());
   TH1D *phi1 = (TH1D*)phi0->Clone("phi1");
   TH2D *map0 = new TH2D("map0","map0;#eta;#phi",50, -5, 5,150, -TMath::Pi(), TMath::Pi());
   TH2D *map1 = (TH2D*)map0->Clone("map1");

  // TH2D *ptrat0 = new TH2D("ptrat0","ptrat0;p_{T}^{gen #gamma};p_{T}^{reco #gamma}",50,0,150,50,0,150);
  // TH2D *ptrat1 = (TH2D*)ptrat0->Clone("ptrat1");
  const TCut noSpike = "( 1 - (eRight+eLeft+eTop+eBottom)/eMax ) <= 0.9 && abs(seedTime) <=3 && sigmaIetaIeta>=0.002 && sigmaIphiIphi>=0.002 ";
  const TCut ptCut = noSpike && "pt>20";
  const TCut etaCut = "abs(eta)>1.479";
  const TCut hltCut = "HLT_HISinglePhoton40_v2==1";
  int rerecoEntries = rereco->GetEntries(ptCut&&etaCut);
  int oldEntries = old->GetEntries(ptCut&&etaCut&&hltCut);
  cout << "rerecoEntries : " << rerecoEntries << "  ,  oldEntries : " << oldEntries << endl; 
    //rereco->Project(pt0->GetName(), "pt",etaCut);
    //old->Project(pt1->GetName(),"pt",etaCut);
 
    rereco->Project(pt0->GetName(), "pt",ptCut&&etaCut);
    old->Project(pt1->GetName(),"pt",ptCut&&etaCut&&hltCut);
    rereco->Project(eta0->GetName(),"eta",ptCut&&etaCut);
    old->Project(eta1->GetName(),"eta",ptCut&&etaCut&&hltCut);
    rereco->Project(phi0->GetName(),"phi",ptCut&&etaCut);
    old->Project(phi1->GetName(),"phi",ptCut&&etaCut&&hltCut);
    rereco->Project(map0->GetName(),"phi:eta",ptCut&&etaCut);
    old->Project(map1->GetName(),"phi:eta",ptCut&&etaCut&&hltCut);

  // rereco->Project(ptrat0->GetName(),"pt:genMatchedPt");
  // old->Project(ptrat1->GetName(),"pt:genMatchedPt");
    
    TLegend* l1 = new TLegend(0.8, 0.85, 0.95, 0.95);
    l1 -> AddEntry((TObject*)0, "HIHighPt_Photon40_Skim"); 
    l1 -> AddEntry(pt0, "53X");
    l1 -> AddEntry(pt1, "44X");

    pt0->Scale(1./rerecoEntries);
    pt1->Scale(1./oldEntries);
    eta0->Scale(1./rerecoEntries);
    eta1->Scale(1./oldEntries);
    phi0->Scale(1./rerecoEntries);
    phi1->Scale(1./oldEntries);
    map0->Scale(1./rerecoEntries);
    map1->Scale(1./oldEntries);

    TCanvas *c1 = new TCanvas();
    pt0->Draw();//rereco 0
    pt0->SetMarkerColor(kRed);
    pt0->SetLineColor(kRed);
    pt1->Draw("same");//old 1
    c1->SetLogy();
    c1->SaveAs("rereco/endcap/pt.pdf");
    //l1 -> Draw();

    TCanvas *c2 = new TCanvas();
    eta0->Draw();
    eta0->SetMarkerColor(kRed);
    eta0->SetLineColor(kRed);
    eta1->Draw("same");
    c2->SaveAs("rereco/endcap/eta_pt20.pdf");

    TCanvas *c3 = new TCanvas();
    phi0->Draw();
    phi0->SetMarkerColor(kRed);
    phi0->SetLineColor(kRed);
    phi1->Draw("same");
    c3->SaveAs("rereco/endcap/phi_pt20.pdf");

    TCanvas *c4 = new TCanvas();
    map0->Draw("colz");
    c4->SaveAs("rereco/endcap/map_new_pt20.pdf");

    TCanvas *c5 = new TCanvas();
    map1->Draw("colz");
    c5->SaveAs("rereco/endcap/map_old_pt20.pdf");

  //  TCanvas *c10 = new TCanvas();//for MC
  // ptrat0->Draw("colz");
  // ptrat0->ProfileX("prof0")->Draw("same");
  // c10->SetLogz();
  // c10->SaveAs("rereco/endcap/corr_new_pt20.pdf");

  // TCanvas *c11 = new TCanvas();//for MC
  // ptrat1->Draw("colz");
  // ptrat1->ProfileX("prof1")->Draw("same");
  // c11->SetLogz();
  // c11->SaveAs("rereco/endcap/corr_old_pt20.pdf");

  // //isolation vars
   TH1D *cc40 = new TH1D("cc40","cc40;cc4",50,-50, 150);
   TH1D *cc41 = (TH1D*)cc40->Clone("cc41");
   TH1D *cr40 = new TH1D("cr40","cr40;cr4",50,-50,350);
   TH1D *cr41 = (TH1D*)cr40->Clone("cr41");
   TH1D *ct4PtCut200 = new TH1D("ct4PtCut200","ct4PtCut200;ct4PtCut20",50,-50,350);
   TH1D *ct4PtCut201 = (TH1D*)ct4PtCut200->Clone("ct4PtCut201");
  // TH1D *genCalIsoDR040 = new TH1D("genCalIsoDR040","genCalIsoDR040;genCalIsoDR04",50,0,500);
  // TH1D *genCalIsoDR041 = (TH1D*)genCalIsoDR040->Clone("genCalIsoDR041");
   TH1D *sigmaIetaIeta0 = new TH1D("sigmaIetaIeta0",";#sigma_{#eta#eta}",50, 0, 0.05);
   TH1D *sigmaIetaIeta1 = (TH1D*)sigmaIetaIeta0->Clone("sigmaIetaIeta1");
  TH1D *ecalRecHitSumEtConeDR040 =  new TH1D("ecalRecHitSumEtConeDR040",
					     "ecalRecHitSumEtConeDR040;ecalRecHitSumEtConeDR04",50,-50, 350);
  TH1D *ecalRecHitSumEtConeDR041 = (TH1D*)ecalRecHitSumEtConeDR040->Clone("ecalRecHitSumEtConeDR041");
  TH1D *hcalTowerSumEtConeDR040 =  new TH1D("hcalTowerSumEtConeDR040",
					     "hcalTowerSumEtConeDR040;hcalTowerSumEtConeDR04",50,-50, 350);
  TH1D *hcalTowerSumEtConeDR041 = (TH1D*)hcalTowerSumEtConeDR040->Clone("hcalTowerSumEtConeDR041");
  TH1D *trkSumPtHollowConeDR040 = new TH1D("trkSumPtHollowConeDR040",
					   "trkSumPtHollowConeDR040;trkSumPtHollowConeDR04",50,-50, 350);
  TH1D *trkSumPtHollowConeDR041 = (TH1D*)trkSumPtHollowConeDR040->Clone("trkSumPtHollowConeDR041");

  TH1D *isEle0 = new TH1D("isEle0",
					   "isEle0;isEle",2,0, 2);
  TH1D *isEle1 = (TH1D*)isEle0->Clone("isEle1");

   rereco->Project(cc40->GetName(),"cc4",ptCut&&etaCut);
   old->Project(cc41->GetName(),"cc4",ptCut&&etaCut&&hltCut);
   rereco->Project(cr40->GetName(),"cr4",ptCut&&etaCut);
   old->Project(cr41->GetName(),"cr4",ptCut&&etaCut&&hltCut);
   rereco->Project(ct4PtCut200->GetName(),"ct4PtCut20",ptCut&&etaCut);
   old->Project(ct4PtCut201->GetName(),"ct4PtCut20",ptCut&&etaCut&&hltCut);
  // rereco->Project(genCalIsoDR040->GetName(),"genCalIsoDR04",ptCut&&etaCut);
  // old->Project(genCalIsoDR041->GetName(),"genCalIsoDR04",ptCut&&etaCut);
   rereco->Project(sigmaIetaIeta0->GetName(),"sigmaIetaIeta",ptCut&&etaCut);
   old->Project(sigmaIetaIeta1->GetName(),"sigmaIetaIeta",ptCut&&etaCut&&hltCut);
  rereco->Project(ecalRecHitSumEtConeDR040->GetName(),"ecalRecHitSumEtConeDR04",ptCut&&etaCut);
  old->Project(ecalRecHitSumEtConeDR041->GetName(),"ecalRecHitSumEtConeDR04",ptCut&&etaCut&&hltCut);
  rereco->Project(hcalTowerSumEtConeDR040->GetName(),"hcalTowerSumEtConeDR04",ptCut&&etaCut);
  old->Project(hcalTowerSumEtConeDR041->GetName(),"hcalTowerSumEtConeDR04",ptCut&&etaCut&&hltCut);
  rereco->Project(trkSumPtHollowConeDR040->GetName(),"trkSumPtHollowConeDR04",ptCut&&etaCut);
  old->Project(trkSumPtHollowConeDR041->GetName(),"trkSumPtHollowConeDR04",ptCut&&etaCut&&hltCut);
  rereco->Project(isEle0->GetName(),"isEle",ptCut&&etaCut);
  old->Project(isEle1->GetName(),"isEle",ptCut&&etaCut&&hltCut);
  
  cc40->Scale(1./rerecoEntries);
   cc41->Scale(1./oldEntries);
   cr40->Scale(1./rerecoEntries);
   cr41->Scale(1./oldEntries);
   ct4PtCut200->Scale(1./rerecoEntries);
   ct4PtCut201->Scale(1./oldEntries);
  // genCalIsoDR040->Scale(1./rerecoEntries);
  // genCalIsoDR041->Scale(1./oldEntries);
   sigmaIetaIeta0->Scale(1./rerecoEntries);
   sigmaIetaIeta1->Scale(1./oldEntries);
  ecalRecHitSumEtConeDR040->Scale(1./rerecoEntries);
  ecalRecHitSumEtConeDR041->Scale(1./oldEntries);
  hcalTowerSumEtConeDR040->Scale(1./rerecoEntries);
  hcalTowerSumEtConeDR041->Scale(1./oldEntries);
  trkSumPtHollowConeDR040->Scale(1./rerecoEntries);
  trkSumPtHollowConeDR041->Scale(1./oldEntries);

   TCanvas *c6 = new TCanvas();
   cc40->Draw();
   cc40->SetMarkerColor(kRed);
   cc40->SetLineColor(kRed);
   cc41->Draw("same");
   c6->SetLogy();
   c6->SaveAs("rereco/endcap/cc4_pt20.pdf");

    TCanvas *c7 = new TCanvas();
    cr40->Draw();
    cr40->SetMarkerColor(kRed);
    cr40->SetLineColor(kRed);
    cr41->Draw("same");
    c7->SetLogy();
    c7->SaveAs("rereco/endcap/cr4_pt20.pdf");
  
    TCanvas *c8 = new TCanvas();
    ct4PtCut200->Draw();
    ct4PtCut200->SetMarkerColor(kRed);
    ct4PtCut200->SetLineColor(kRed);
    ct4PtCut201->Draw("same");
    c8->SetLogy();
    c8->SaveAs("rereco/endcap/ct4PtCut20_pt20.pdf");

  // TCanvas *c9 = new TCanvas();
  // genCalIsoDR040->Draw();
  // genCalIsoDR041->SetMarkerColor(kRed);
  // genCalIsoDR041->SetLineColor(kRed);
  // genCalIsoDR041->Draw("same");
  // c9->SetLogy();
  // c9->SaveAs("rereco/endcap/genCalIsoDR04_pt20.pdf");

   TCanvas *c12 = new TCanvas();
   sigmaIetaIeta0->Draw();
   sigmaIetaIeta0->SetMarkerColor(kRed);
   sigmaIetaIeta0->SetLineColor(kRed);
   sigmaIetaIeta1->Draw("same");
   c12->SaveAs("rereco/endcap/sigmaIetaIeta_pt20.pdf");

  TCanvas *c13 = new TCanvas();
  ecalRecHitSumEtConeDR040->Draw();
  ecalRecHitSumEtConeDR040->SetMarkerColor(kRed);
  ecalRecHitSumEtConeDR040->SetLineColor(kRed);
  ecalRecHitSumEtConeDR041->Draw("same");
  c13->SetLogy();
  c13->SaveAs("rereco/endcap/ecalRecHitSumEtConeDR04_pt20.pdf");

  TCanvas *c14 = new TCanvas();
  hcalTowerSumEtConeDR040->Draw();
  hcalTowerSumEtConeDR040->SetMarkerColor(kRed);
  hcalTowerSumEtConeDR040->SetLineColor(kRed);
  hcalTowerSumEtConeDR041->Draw("same");
  c14->SetLogy();
  c14->SaveAs("rereco/endcap/hcalTowerSumEtConeDR04_pt20.pdf");

  TCanvas *c15 = new TCanvas();
  trkSumPtHollowConeDR040->Draw();
  trkSumPtHollowConeDR040->SetMarkerColor(kRed);
  trkSumPtHollowConeDR040->SetLineColor(kRed);
  trkSumPtHollowConeDR041->Draw("same");
  c15->SetLogy();
  c15->SaveAs("rereco/endcap/trkSumPtHollowConeDR04_pt20.pdf");

  TCanvas *c16 = new TCanvas();
  isEle0->Draw();
  isEle0->SetMarkerColor(kRed);
  isEle0->SetLineColor(kRed);
  isEle1->Draw("same");
  //c16->SetLogy();
  c16->SaveAs("rereco/endcap/isEle_pt20.pdf");


}
 */

 } 
