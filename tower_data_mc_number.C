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

void tower_data_mc_number()
{
    const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
   // const TCut eventCut = "";
    const TCut eventCut = runCut && lumiCut;
    const TCut threCut = "abs(eta)>3.0 && abs(eta)<5.0 && et>1.4";
    const int Ncut = 5;
    TH1::SetDefaultSumw2();
    gStyle -> SetOptStat(0);


// ===================================================================================
// Get Trees from data & mc files.
// ===================================================================================

    TFile *dataf = new TFile("files/HiForest_20141011.root");
    TTree *datat_evt = (TTree*) dataf -> Get("hiEvtAnalyzer/HiTree");
    TTree *datat_skim = (TTree*) dataf -> Get("skimanalysis/HltTree");
    TTree *datat_hlt = (TTree*) dataf -> Get("hltanalysis/HltTree");
    TTree *datat_recTower = (TTree*) dataf -> Get("rechitanalyzer/tower");
    TTree *datat_pfTower = (TTree*) dataf -> Get("pfTowers/tower");
    datat_recTower -> AddFriend(datat_hlt);
    datat_recTower -> AddFriend(datat_evt);
    datat_recTower -> AddFriend(datat_skim);
    datat_pfTower -> AddFriend(datat_hlt);
    datat_pfTower -> AddFriend(datat_evt);
    datat_pfTower -> AddFriend(datat_skim);
    cout << datat_evt << ", " << datat_skim << ", " << datat_hlt << ", " << datat_recTower << ", " << datat_pfTower << ", " << endl;
    cout << "data pfTower tree Friends : " << datat_pfTower -> GetListOfFriends() << endl;
    cout << "data rechitTower tree Friends : " << datat_recTower -> GetListOfFriends() << endl;
    int Nrec_dataEntries = datat_recTower -> GetEntries();
    int Npf_dataEntries = datat_pfTower -> GetEntries();
    cout << "# of rec DATA entries = "  << Nrec_dataEntries << ", pf DATA entries = " << Npf_dataEntries<< endl;

    TFile *mcf = new TFile("files/centrality_PbPb_minbias_MC.root");
    TTree *mct_evt = (TTree*) mcf -> Get("hiEvtAnalyzer/HiTree");
    TTree *mct_skim = (TTree*) mcf -> Get("skimanalysis/HltTree");
    TTree *mct_hlt = (TTree*) mcf -> Get("hltanalysis/HltTree");
    TTree *mct_recTower = (TTree*) mcf -> Get("rechitanalyzer/tower");
    TTree *mct_pfTower = (TTree*) mcf -> Get("pfTowers/tower");
    mct_recTower -> AddFriend(mct_hlt);
    mct_recTower -> AddFriend(mct_evt);
    mct_recTower -> AddFriend(mct_skim);
    mct_pfTower -> AddFriend(mct_hlt);
    mct_pfTower -> AddFriend(mct_evt);
    mct_pfTower -> AddFriend(mct_skim);
    int Nrec_mcEntries = mct_recTower -> GetEntries();
    int Npf_mcEntries = mct_pfTower -> GetEntries();
    cout << "# of rec MC entries = "  << Nrec_mcEntries << ", pf MC entries = " << Npf_mcEntries<< endl;
   
    // int Nevt_mct = mct_evt -> GetEntries();
    // cout << "# of MC events = " << Nevt_mct << endl;

// ***********************************************************************************************************    

// ===============================================================================================
// [et] Define et histograms (data/mc , pf tree/rechit tree) 
// ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;
    
    TH1D* data_pf_n[5];
    TH1D* data_rec_n[5];
    TH1D* mc_pf_n[5];
    TH1D* mc_rec_n[5];
    TH1D* effi_pf_n[5];// effi = data/mc
    TH1D* effi_rec_n[5];// effi = data/mc
    for(int i=0; i<Ncut; i++)
    {
        data_pf_n[i] = new TH1D(Form("data_pf_n%d",i), ";pfTowers/n;Normalized Events",200,0,4500);
        data_pf_n[i] -> SetMarkerStyle(20);
        data_pf_n[i] -> SetMarkerSize(0.7);
        data_pf_n[i] -> SetMarkerColor(21+i);
        data_pf_n[i] -> SetLabelSize(0.03);

        mc_pf_n[i] = new TH1D(Form("mc_pf_n%d",i), ";pfTowers/n;Normalized Events",200,0,4500);
       // mc_pf_n[i] -> SetMarkerStyle(24+i);
       // mc_pf_n[i] -> SetMarkerSize(0.7);
        mc_pf_n[i] -> SetMarkerColor(21+i); //marker color
        mc_pf_n[i] -> SetLineColor(21+i); //line color
        mc_pf_n[i] -> SetLabelSize(0.03);
   
        data_rec_n[i] = new TH1D(Form("data_rec_n%d",i), ";rechitTowers/n;Normalized Events",200,0,4500);
        data_rec_n[i] -> SetMarkerStyle(20);
        data_rec_n[i] -> SetMarkerSize(0.7);
        data_rec_n[i] -> SetMarkerColor(31+i);
        data_rec_n[i] -> SetLabelSize(0.03);

        mc_rec_n[i] = new TH1D(Form("mc_rec_n%d",i), ";rechitTowers/n;Normalized Events",200,0,4500);
       // mc_rec_n[i] -> SetMarkerStyle(24+i);
       // mc_rec_n[i] -> SetMarkerSize(0.7);
        mc_rec_n[i] -> SetMarkerColor(31+i); //marker color
        mc_rec_n[i] -> SetLineColor(31+i); //line color
        mc_rec_n[i] -> SetLabelSize(0.03);
        
        effi_pf_n[i] = new TH1D(Form("effi_pf_n%d",i), ";pfTowers/n;DATA/MC",200,0,4500);
        effi_pf_n[i] -> SetMarkerStyle(20);
        effi_pf_n[i] -> SetMarkerSize(0.7);
        effi_pf_n[i] -> SetMarkerColor(21+i);
        effi_pf_n[i] -> SetLabelSize(0.03);
        effi_pf_n[i] -> SetAxisRange(0.5,1.5,"Y");

        effi_rec_n[i] = new TH1D(Form("effi_rec_n%d",i), ";rechitTowers/n;DATA/MC",200,0,4500);
        effi_rec_n[i] -> SetMarkerStyle(20);
        effi_rec_n[i] -> SetMarkerSize(0.7);
        effi_rec_n[i] -> SetMarkerColor(31+i);
        effi_rec_n[i] -> SetLabelSize(0.03);
        effi_rec_n[i] -> SetAxisRange(0.5,1.5,"Y");

    }

   
// ===============================================================================================
// [et] Fill histogram by using <Draw> function in TTree.
// ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;

//    datat_pfTower -> Draw("et >>+ data_pf_n0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
//    data_pf_n[0] = (TH1D*)gDirectory->Get("data_pf_n0");
/*    datat_pfTower -> Draw("et >>+ data_pf_n1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    data_pf_n[1] = (TH1D*)gDirectory->Get("data_pf_n1");
    datat_pfTower -> Draw("et >>+ data_pf_n2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    data_pf_n[2] = (TH1D*)gDirectory->Get("data_pf_n2");
    datat_pfTower -> Draw("et >>+ data_pf_n3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    data_pf_n[3] = (TH1D*)gDirectory->Get("data_pf_n3");
  */  datat_pfTower -> Draw("n >>+ data_pf_n4",eventCut && threCut && "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    data_pf_n[4] = (TH1D*)gDirectory->Get("data_pf_n4");
    cout << "data_pf_n finished"<< endl;

//    mct_pfTower -> Draw("et >>+ mc_pf_n0");
//    mc_pf_n[0] = (TH1D*)gDirectory->Get("mc_pf_n0");
/*    mct_pfTower -> Draw("et >>+ mc_pf_n1","pprimaryVertexFilter==1");
    mc_pf_n[1] = (TH1D*)gDirectory->Get("mc_pf_n1");
    mct_pfTower -> Draw("et >>+ mc_pf_n2","phltPixelClusterShapeFilter==1");
    mc_pf_n[2] = (TH1D*)gDirectory->Get("mc_pf_n2");
    mct_pfTower -> Draw("et >>+ mc_pf_n3","phfCoincFilter3==1");
    mc_pf_n[3] = (TH1D*)gDirectory->Get("mc_pf_n3");
 */   mct_pfTower -> Draw("n >>+ mc_pf_n4",threCut && "pcollisionEventSelection==1");
    mc_pf_n[4] = (TH1D*)gDirectory->Get("mc_pf_n4");
    cout << "mc_pf_n"<< endl;

//    datat_recTower -> Draw("et >>+ data_rec_n0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
//    data_rec_n[0] = (TH1D*)gDirectory->Get("data_rec_n0");
/*    datat_recTower -> Draw("et >>+ data_rec_n1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    data_rec_n[1] = (TH1D*)gDirectory->Get("data_rec_n1");
    datat_recTower -> Draw("et >>+ data_rec_n2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    data_rec_n[2] = (TH1D*)gDirectory->Get("data_rec_n2");
    datat_recTower -> Draw("et >>+ data_rec_n3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    data_rec_n[3] = (TH1D*)gDirectory->Get("data_rec_n3");
  */  datat_recTower -> Draw("n >>+ data_rec_n4",eventCut &&threCut &&   "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    data_rec_n[4] = (TH1D*)gDirectory->Get("data_rec_n4");

//    mct_recTower -> Draw("et >>+ mc_rec_n0");
//    mc_rec_n[0] = (TH1D*)gDirectory->Get("mc_rec_n0");
 /*   mct_recTower -> Draw("et >>+ mc_rec_n1","pprimaryVertexFilter==1");
    mc_rec_n[1] = (TH1D*)gDirectory->Get("mc_rec_n1");
    mct_recTower -> Draw("et >>+ mc_rec_n2","phltPixelClusterShapeFilter==1");
    mc_rec_n[2] = (TH1D*)gDirectory->Get("mc_rec_n2");
    mct_recTower -> Draw("et >>+ mc_rec_n3","phfCoincFilter3==1");
    mc_rec_n[3] = (TH1D*)gDirectory->Get("mc_rec_n3");
  */  mct_recTower -> Draw("n >>+ mc_rec_n4",threCut && "pcollisionEventSelection==1");
    mc_rec_n[4] = (TH1D*)gDirectory->Get("mc_rec_n4");

    //for(int i=0; i<Ncut; i++){
    int i=4;    
    data_pf_n[i] -> Scale(1./data_pf_n[i]->Integral("width"));
        data_rec_n[i] -> Scale(1./data_rec_n[i]->Integral("width"));
        mc_pf_n[i] -> Scale(1./mc_pf_n[i]->Integral("width"));
        mc_rec_n[i] -> Scale(1./mc_rec_n[i]->Integral("width"));
  /*    data_pf_n[i] -> Scale(1./Npf_dataEntries);
        data_rec_n[i] -> Scale(1./Nrec_dataEntries);
        mc_pf_n[i] -> Scale(1./Npf_mcEntries);
        mc_rec_n[i] -> Scale(1./Nrec_mcEntries);
*/  //  }
    effi_pf_n[i] -> Divide(data_pf_n[i],mc_pf_n[i]);
    effi_rec_n[i] -> Divide(data_rec_n[i],mc_rec_n[i]);
 

 
// ===============================================================================================
// [et] Draw histograms in Canvas. 
// ===============================================================================================
    cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
/*
    TCanvas *c_data_n = new TCanvas("c_data_n", "c_data_n", 800,400);
    c_data_n -> Divide(2,1,0);
    c_data_n -> cd(1);
    data_pf_n[0] -> Draw("ep");
    data_pf_n[1] -> Draw("same ep");
    data_pf_n[2] -> Draw("same ep");
    data_pf_n[3] -> Draw("same ep");
    data_pf_n[4] -> Draw("same ep");
    
    c_data_n -> cd(2);
    data_rec_n[0] -> Draw("ep");
    data_rec_n[1] -> Draw("same ep");
    data_rec_n[2] -> Draw("same ep");
    data_rec_n[3] -> Draw("same ep");
    data_rec_n[4] -> Draw("same ep");
 */
    TCanvas *c_mc_n = new TCanvas("c_mc_n", "c_mc_n", 800,400);
    c_mc_n -> Divide(2,1,0);
    c_mc_n -> cd(1);
//    mc_pf_n[0] -> Draw("ep");
/*    mc_pf_n[1] -> Draw("same ep");
    mc_pf_n[2] -> Draw("same ep");
    mc_pf_n[3] -> Draw("same ep");
 */ //  mc_pf_n[4] -> Draw("same ep");
    
    c_mc_n -> cd(2);
 //   mc_rec_n[0] -> Draw("ep");
/*    mc_rec_n[1] -> Draw("same ep");
    mc_rec_n[2] -> Draw("same ep");
    mc_rec_n[3] -> Draw("same ep");
  *///  mc_rec_n[4] -> Draw("same ep");
    
    TCanvas *c_compare_n = new TCanvas("c_compare_n", "c_compare_n", 400,400);
    c_compare_n -> SetLogy();

    data_pf_n[4] -> Draw("ep");
    mc_pf_n[4] -> Draw("same ehist");
    data_rec_n[4] -> Draw("same ep");
    mc_rec_n[4] -> Draw("same ehist");
   
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetTextFont(43);
    leg->SetTextSize(15);
    leg->AddEntry(data_pf_n[4],"pfTower DATA","p");
    leg->AddEntry(mc_pf_n[4],"pfTower MC","l");
    leg->AddEntry(data_rec_n[4],"rechitTower DATA","p");
    leg->AddEntry(mc_rec_n[4],"rechitTower MC","l");
    leg->Draw();
 
    TCanvas *c_compare_n_effi = new TCanvas("c_compare_n_effi", "c_compare_n_effi", 400,200);
    effi_pf_n[4] -> Draw("ep");
    effi_rec_n[4] -> Draw("same ep");
    c_compare_n_effi->SaveAs("pdf/tower_n_effi.pdf");
   // leg->Draw();

    c_compare_n -> SaveAs("pdf/tower_n.pdf");
    
    TCanvas *c_compare_n_seperately = new TCanvas("c_compare_n_seperately", "c_compare_n_seperately", 1000,400);
    c_compare_n_seperately -> Divide(2,1);
    c_compare_n_seperately -> cd(1);
  //  c_compare_n_seperately1 -> SetLogy();
    data_pf_n[4] -> Draw("ep");
    mc_pf_n[4] -> Draw("same ehist");
    //leg -> Draw();

    c_compare_n_seperately -> cd(2);
   // c_compare_n_seperately2 -> SetLogy();
    data_rec_n[4] -> Draw("ep");
    mc_rec_n[4] -> Draw("same ehist");

    c_compare_n_seperately -> SaveAs("pdf/tower_n_sep.pdf");





 } 
