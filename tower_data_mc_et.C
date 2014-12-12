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

void tower_data_mc_et()
{
    const TCut runCut = "run==181611";
    const TCut lumiCut = "lumi>=1 && lumi<=895";
   // const TCut eventCut = "";
    const TCut eventCut = runCut && lumiCut;
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
    
    TH1D* data_pf_et[5];
    TH1D* data_rec_et[5];
    TH1D* mc_pf_et[5];
    TH1D* mc_rec_et[5];
    TH1D* effi_pf_et[5];// effi = data/mc
    TH1D* effi_rec_et[5];// effi = data/mc
    for(int i=0; i<Ncut; i++)
    {
        data_pf_et[i] = new TH1D(Form("data_pf_et%d",i), ";pfTowers/et;Normalized Events",100,0,100);
        data_pf_et[i] -> SetMarkerStyle(20);
        data_pf_et[i] -> SetMarkerSize(0.7);
        data_pf_et[i] -> SetMarkerColor(21+i);
        data_pf_et[i] -> SetLabelSize(0.03);

        mc_pf_et[i] = new TH1D(Form("mc_pf_et%d",i), ";pfTowers/et;Normalized Events",100,0,100);
       // mc_pf_et[i] -> SetMarkerStyle(24+i);
       // mc_pf_et[i] -> SetMarkerSize(0.7);
        mc_pf_et[i] -> SetMarkerColor(21+i); //marker color
        mc_pf_et[i] -> SetLineColor(21+i); //line color
        mc_pf_et[i] -> SetLabelSize(0.03);
   
        data_rec_et[i] = new TH1D(Form("data_rec_et%d",i), ";rechitTowers/et;Normalized Events",100,0,100);
        data_rec_et[i] -> SetMarkerStyle(20);
        data_rec_et[i] -> SetMarkerSize(0.7);
        data_rec_et[i] -> SetMarkerColor(31+i);
        data_rec_et[i] -> SetLabelSize(0.03);

        mc_rec_et[i] = new TH1D(Form("mc_rec_et%d",i), ";rechitTowers/et;Normalized Events",100,0,100);
       // mc_rec_et[i] -> SetMarkerStyle(24+i);
       // mc_rec_et[i] -> SetMarkerSize(0.7);
        mc_rec_et[i] -> SetMarkerColor(31+i); //marker color
        mc_rec_et[i] -> SetLineColor(31+i); //line color
        mc_rec_et[i] -> SetLabelSize(0.03);
        
        effi_pf_et[i] = new TH1D(Form("effi_pf_et%d",i), ";pfTowers/et;DATA/MC",100,0,100);
        effi_pf_et[i] -> SetMarkerStyle(20);
        effi_pf_et[i] -> SetMarkerSize(0.7);
        effi_pf_et[i] -> SetMarkerColor(21+i);
        effi_pf_et[i] -> SetLabelSize(0.03);
        effi_pf_et[i] -> SetAxisRange(0.5,1.5,"Y");

        effi_rec_et[i] = new TH1D(Form("effi_rec_et%d",i), ";rechitTowers/et;DATA/MC",100,0,100);
        effi_rec_et[i] -> SetMarkerStyle(20);
        effi_rec_et[i] -> SetMarkerSize(0.7);
        effi_rec_et[i] -> SetMarkerColor(31+i);
        effi_rec_et[i] -> SetLabelSize(0.03);
        effi_rec_et[i] -> SetAxisRange(0.5,1.5,"Y");

    }

   
// ===============================================================================================
// [et] Fill histogram by using <Draw> function in TTree.
// ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;

//    datat_pfTower -> Draw("et >>+ data_pf_et0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
//    data_pf_et[0] = (TH1D*)gDirectory->Get("data_pf_et0");
/*    datat_pfTower -> Draw("et >>+ data_pf_et1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    data_pf_et[1] = (TH1D*)gDirectory->Get("data_pf_et1");
    datat_pfTower -> Draw("et >>+ data_pf_et2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    data_pf_et[2] = (TH1D*)gDirectory->Get("data_pf_et2");
    datat_pfTower -> Draw("et >>+ data_pf_et3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    data_pf_et[3] = (TH1D*)gDirectory->Get("data_pf_et3");
  */  datat_pfTower -> Draw("et >>+ data_pf_et4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    data_pf_et[4] = (TH1D*)gDirectory->Get("data_pf_et4");
    cout << "data_pf_et finished"<< endl;

//    mct_pfTower -> Draw("et >>+ mc_pf_et0");
//    mc_pf_et[0] = (TH1D*)gDirectory->Get("mc_pf_et0");
/*    mct_pfTower -> Draw("et >>+ mc_pf_et1","pprimaryVertexFilter==1");
    mc_pf_et[1] = (TH1D*)gDirectory->Get("mc_pf_et1");
    mct_pfTower -> Draw("et >>+ mc_pf_et2","phltPixelClusterShapeFilter==1");
    mc_pf_et[2] = (TH1D*)gDirectory->Get("mc_pf_et2");
    mct_pfTower -> Draw("et >>+ mc_pf_et3","phfCoincFilter3==1");
    mc_pf_et[3] = (TH1D*)gDirectory->Get("mc_pf_et3");
 */   mct_pfTower -> Draw("et >>+ mc_pf_et4","pcollisionEventSelection==1");
    mc_pf_et[4] = (TH1D*)gDirectory->Get("mc_pf_et4");
    cout << "mc_pf_et"<< endl;

//    datat_recTower -> Draw("et >>+ data_rec_et0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
//    data_rec_et[0] = (TH1D*)gDirectory->Get("data_rec_et0");
/*    datat_recTower -> Draw("et >>+ data_rec_et1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    data_rec_et[1] = (TH1D*)gDirectory->Get("data_rec_et1");
    datat_recTower -> Draw("et >>+ data_rec_et2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    data_rec_et[2] = (TH1D*)gDirectory->Get("data_rec_et2");
    datat_recTower -> Draw("et >>+ data_rec_et3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    data_rec_et[3] = (TH1D*)gDirectory->Get("data_rec_et3");
  */  datat_recTower -> Draw("et >>+ data_rec_et4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    data_rec_et[4] = (TH1D*)gDirectory->Get("data_rec_et4");

//    mct_recTower -> Draw("et >>+ mc_rec_et0");
//    mc_rec_et[0] = (TH1D*)gDirectory->Get("mc_rec_et0");
 /*   mct_recTower -> Draw("et >>+ mc_rec_et1","pprimaryVertexFilter==1");
    mc_rec_et[1] = (TH1D*)gDirectory->Get("mc_rec_et1");
    mct_recTower -> Draw("et >>+ mc_rec_et2","phltPixelClusterShapeFilter==1");
    mc_rec_et[2] = (TH1D*)gDirectory->Get("mc_rec_et2");
    mct_recTower -> Draw("et >>+ mc_rec_et3","phfCoincFilter3==1");
    mc_rec_et[3] = (TH1D*)gDirectory->Get("mc_rec_et3");
  */  mct_recTower -> Draw("et >>+ mc_rec_et4","pcollisionEventSelection==1");
    mc_rec_et[4] = (TH1D*)gDirectory->Get("mc_rec_et4");

    //for(int i=0; i<Ncut; i++){
    int i=4;    
    data_pf_et[i] -> Scale(1./data_pf_et[i]->Integral("width"));
        data_rec_et[i] -> Scale(1./data_rec_et[i]->Integral("width"));
        mc_pf_et[i] -> Scale(1./mc_pf_et[i]->Integral("width"));
        mc_rec_et[i] -> Scale(1./mc_rec_et[i]->Integral("width"));
  /*    data_pf_et[i] -> Scale(1./Npf_dataEntries);
        data_rec_et[i] -> Scale(1./Nrec_dataEntries);
        mc_pf_et[i] -> Scale(1./Npf_mcEntries);
        mc_rec_et[i] -> Scale(1./Nrec_mcEntries);
*/  //  }
    effi_pf_et[i] -> Divide(data_pf_et[i],mc_pf_et[i]);
    effi_rec_et[i] -> Divide(data_rec_et[i],mc_rec_et[i]);
 

 
// ===============================================================================================
// [et] Draw histograms in Canvas. 
// ===============================================================================================
    cout << "LET'S DRAW HISTOGRAMS IN CANVAS" << endl;
/*
    TCanvas *c_data_et = new TCanvas("c_data_et", "c_data_et", 800,400);
    c_data_et -> Divide(2,1,0);
    c_data_et -> cd(1);
    data_pf_et[0] -> Draw("ep");
    data_pf_et[1] -> Draw("same ep");
    data_pf_et[2] -> Draw("same ep");
    data_pf_et[3] -> Draw("same ep");
    data_pf_et[4] -> Draw("same ep");
    
    c_data_et -> cd(2);
    data_rec_et[0] -> Draw("ep");
    data_rec_et[1] -> Draw("same ep");
    data_rec_et[2] -> Draw("same ep");
    data_rec_et[3] -> Draw("same ep");
    data_rec_et[4] -> Draw("same ep");
 */
    TCanvas *c_mc_et = new TCanvas("c_mc_et", "c_mc_et", 800,400);
    c_mc_et -> Divide(2,1,0);
    c_mc_et -> cd(1);
//    mc_pf_et[0] -> Draw("ep");
/*    mc_pf_et[1] -> Draw("same ep");
    mc_pf_et[2] -> Draw("same ep");
    mc_pf_et[3] -> Draw("same ep");
 */ //  mc_pf_et[4] -> Draw("same ep");
    
    c_mc_et -> cd(2);
 //   mc_rec_et[0] -> Draw("ep");
/*    mc_rec_et[1] -> Draw("same ep");
    mc_rec_et[2] -> Draw("same ep");
    mc_rec_et[3] -> Draw("same ep");
  *///  mc_rec_et[4] -> Draw("same ep");
    
    TCanvas *c_compare_et = new TCanvas("c_compare_et", "c_compare_et", 400,400);
    c_compare_et -> SetLogy();

    data_pf_et[4] -> Draw("ep");
    mc_pf_et[4] -> Draw("same ehist");
    data_rec_et[4] -> Draw("same ep");
    mc_rec_et[4] -> Draw("same ehist");
    
    TCanvas *c_compare_et_effi = new TCanvas("c_compare_et_effi", "c_compare_et_effi", 400,200);
    effi_pf_et[4] -> Draw("ep");
    effi_rec_et[4] -> Draw("same ep");
    c_compare_et_effi->SaveAs("pdf/tower_et_effi.pdf");

    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetTextFont(43);
    leg->SetTextSize(15);
    leg->AddEntry(data_pf_et[4],"pfTower DATA","p");
    leg->AddEntry(mc_pf_et[4],"pfTower MC","l");
    leg->AddEntry(data_rec_et[4],"rechitTower DATA","p");
    leg->AddEntry(mc_rec_et[4],"rechitTower MC","l");
    leg->Draw();

    c_compare_et -> SaveAs("pdf/tower_et.pdf");
    
    TCanvas *c_compare_et_seperately = new TCanvas("c_compare_et_seperately", "c_compare_et_seperately", 1000,400);
    c_compare_et_seperately -> Divide(2,1);
    c_compare_et_seperately -> cd(1);
  //  c_compare_et_seperately1 -> SetLogy();
    data_pf_et[4] -> Draw("ep");
    mc_pf_et[4] -> Draw("same ehist");
    
    c_compare_et_seperately -> cd(2);
   // c_compare_et_seperately2 -> SetLogy();
    data_rec_et[4] -> Draw("ep");
    mc_rec_et[4] -> Draw("same ehist");

    c_compare_et_seperately -> SaveAs("pdf/tower_et_sep.pdf");
    












 /*   double norm_data = HFsum_data[0]->Integral("width");
    double norm_mc = HFsum_mc[0]->Integral("width");
    for(int i=0; i<Ncut; i++){
        if(i==0) HFsum_data[i] -> Draw("ep");
        else HFsum_data[i] -> Draw("same&&ep");
        HFsum_mc[i] -> Draw("same&&ehist");

        //HFsum_data[i] -> Scale(1./Nevt_datat);
        //HFsum_mc[i] -> Scale(1./Nevt_mct);
        HFsum_data[i] -> Scale(1./norm_data);
        HFsum_mc[i] -> Scale(1./norm_mc);

        cout << "The integral of HFsum data " << i << " : " << HFsum_data[i]->Integral("width") << endl;
        cout << "The integral of HFsum mc " << i << " : " << HFsum_mc[i]->Integral("width") << endl;
    }
   
    l1 -> Draw();
    l1_mc -> Draw();
    
    c_HFsum -> SaveAs("pdf/HFsum.pdf");

 
*/


//======================================
//HF sum!!
//======================================

/*
    TLine* t1 = new TLine(0,1,1000,1);
    t1->SetLineWidth(1);
    t1->SetLineStyle(7); // 7 is jumSun , 1 is onSun
    t1->SetLineColor(1); // 2 is red

    TH1D *HFsum_data[5];
    TH1D *HFsum_mc[5];
    for(int i=0; i<Ncut; i++)
    {
        HFsum_data[i] = new TH1D(Form("HFsum_data%d",i), ";hiHF;Normalized Events",100,0,5000);
        HFsum_data[i] -> SetMarkerStyle(20+i);
        HFsum_data[i] -> SetMarkerSize(0.7);
        HFsum_data[i] -> SetMarkerColor(kRed+i);
        HFsum_data[i] -> SetLabelSize(0.03);
   	
	    HFsum_mc[i] = new TH1D(Form("HFsum_mc%d",i), ";hiHF;Normalized Events",100,0,5000);
        //HFsum_mc[i] -> SetMarkerStyle(20);
        //HFsum_mc[i] -> SetMarkerSize(1.0);
        HFsum_mc[i] -> SetLineColor(2+i); //line color
        HFsum_mc[i] -> SetLabelSize(0.03);
    }
    TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 300,300);
	
    c_temp -> cd();
    datat_evt -> Draw("hiHF >>+ HFsum_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    HFsum_data[0] = (TH1D*)gDirectory->Get("HFsum_data0");
    datat_evt -> Draw("hiHF >>+ HFsum_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    HFsum_data[1] = (TH1D*)gDirectory->Get("HFsum_data1");
    datat_evt -> Draw("hiHF >>+ HFsum_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    HFsum_data[2] = (TH1D*)gDirectory->Get("HFsum_data2");
    datat_evt -> Draw("hiHF >>+ HFsum_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    HFsum_data[3] = (TH1D*)gDirectory->Get("HFsum_data3");
    datat_evt -> Draw("hiHF >>+ HFsum_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    HFsum_data[4] = (TH1D*)gDirectory->Get("HFsum_data4");
    mct_evt -> Draw("hiHF >>+ HFsum_mc0");
    HFsum_mc[0] = (TH1D*)gDirectory->Get("HFsum_mc0");
    mct_evt -> Draw("hiHF >>+ HFsum_mc1","pprimaryVertexFilter==1");
    HFsum_mc[1] = (TH1D*)gDirectory->Get("HFsum_mc1");
    mct_evt -> Draw("hiHF >>+ HFsum_mc2","phltPixelClusterShapeFilter==1");
    HFsum_mc[2] = (TH1D*)gDirectory->Get("HFsum_mc2");
    mct_evt -> Draw("hiHF >>+ HFsum_mc3","phfCoincFilter3==1");
    HFsum_mc[3] = (TH1D*)gDirectory->Get("HFsum_mc3");
    mct_evt -> Draw("hiHF >>+ HFsum_mc4","pcollisionEventSelection==1");
    HFsum_mc[4] = (TH1D*)gDirectory->Get("HFsum_mc4");

    TLegend* l1 = new TLegend(0.3, 0.65, 0.6, 0.80, "PbPb Minbias rereco DATA");
   // l1 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    l1 -> AddEntry(HFsum_data[0], "No filters");
    l1 -> AddEntry(HFsum_data[1], "primay vertex filter");
    l1 -> AddEntry(HFsum_data[2], "pixel cluster shape filter");
    l1 -> AddEntry(HFsum_data[3], "HF coinc. 3 filter");
    l1 -> AddEntry(HFsum_data[4], "collision event filter");
  
    TLegend* l1_mc = new TLegend(0.6, 0.65, 0.9, 0.80, "PbPb Minbias MC");
    //l1 -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l1_mc -> AddEntry(HFsum_mc[0], "No filters");
    l1_mc -> AddEntry(HFsum_mc[1], "primay vertex filter");
    l1_mc -> AddEntry(HFsum_mc[2], "pixel cluster shape filter");
    l1_mc -> AddEntry(HFsum_mc[3], "HF coinc. 3 filter");
    l1_mc -> AddEntry(HFsum_mc[4], "collision event filter");
   
    TLegend* l2 = new TLegend(0.4, 0.65, 0.85, 0.80,"PbPb Minbias rereco DATA");
    //l2 -> AddEntry((TObject*)0, "PbPb Minbias rereco DATA"); 
    l2 -> AddEntry(HFsum_data[1], "primay vertex filter");
    l2 -> AddEntry(HFsum_data[2], "pixel cluster shape filter");
    l2 -> AddEntry(HFsum_data[3], "HF coinc. 3 filter");
    l2 -> AddEntry(HFsum_data[4], "collision event filter");

    TLegend* l2_mc = new TLegend(0.4, 0.65, 0.85, 0.80, "PbPb Minbias MC");
    //l2_mc -> AddEntry((TObject*)0, "PbPb Minbias MC"); 
    l2_mc -> AddEntry(HFsum_mc[1], "primay vertex filter");
    l2_mc -> AddEntry(HFsum_mc[2], "pixel cluster shape filter");
    l2_mc -> AddEntry(HFsum_mc[3], "HF coinc. 3 filter");
    l2_mc -> AddEntry(HFsum_mc[4], "collision event filter");

    TCanvas *c_HFsum = new TCanvas("c_HFsum", "c_HFsum", 400,400);
    c_HFsum -> SetLogy();

    double norm_data = HFsum_data[0]->Integral("width");
    double norm_mc = HFsum_mc[0]->Integral("width");
    for(int i=0; i<Ncut; i++){
        if(i==0) HFsum_data[i] -> Draw("ep");
        else HFsum_data[i] -> Draw("same&&ep");
        HFsum_mc[i] -> Draw("same&&ehist");

        //HFsum_data[i] -> Scale(1./Nevt_datat);
        //HFsum_mc[i] -> Scale(1./Nevt_mct);
        HFsum_data[i] -> Scale(1./norm_data);
        HFsum_mc[i] -> Scale(1./norm_mc);

        cout << "The integral of HFsum data " << i << " : " << HFsum_data[i]->Integral("width") << endl;
        cout << "The integral of HFsum mc " << i << " : " << HFsum_mc[i]->Integral("width") << endl;
    }
   
    l1 -> Draw();
    l1_mc -> Draw();
    
    c_HFsum -> SaveAs("pdf/HFsum.pdf");

    TH1D *HFsum_data_effi[5];
    for(int i=0; i<Ncut; i++){
     //   HFsum_data_effi[i] = new TH1D(Form("HFsum_data_effi%d",i), "", 500,0, 5000);
        HFsum_data_effi[i] = (TH1D*)HFsum_data[i]->Clone(Form("HFsum_data_effi%d",i));
        HFsum_data_effi[i] -> SetTitle(";hiHF;Filter Efficiency");
        if(i!=0)
            HFsum_data_effi[i] -> Divide(HFsum_data[i],HFsum_data[0]);
        HFsum_data_effi[i] -> SetMarkerStyle(20+i);
        HFsum_data_effi[i] -> SetMarkerSize(0.7);
        HFsum_data_effi[i] -> SetMarkerColor(kRed+i);
        HFsum_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        //HFsum_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_HFsum_data_effi = new TCanvas("c_HFsum_data_effi", "c_HFsum_data_effi", 400, 400);
    HFsum_data_effi[1] -> Draw("elp"); 
    HFsum_data_effi[2] -> Draw("same ep"); 
    HFsum_data_effi[3] -> Draw("same ep"); 
    HFsum_data_effi[4] -> Draw("same ep"); 
    t1 -> Draw();
    l2 -> Draw();
    
    c_HFsum_data_effi -> SetLogx();
    c_HFsum_data_effi -> SaveAs("pdf/HFsum_data_effi.pdf");

   TH1D *HFsum_mc_effi[5];
    for(int i=0; i<Ncut; i++){
     //   HFsum_mc_effi[i] = new TH1D(Form("HFsum_mc_effi%d",i), "", 500,0, 5000);
        HFsum_mc_effi[i] = (TH1D*)HFsum_mc[i]->Clone(Form("HFsum_mc_effi%d",i));
        HFsum_mc_effi[i] -> SetTitle(";hiHF;Filter Efficiency");
        if(i!=0)
            HFsum_mc_effi[i] -> Divide(HFsum_mc[i],HFsum_mc[0]);
        HFsum_mc_effi[i] -> SetMarkerStyle(20+i);
        HFsum_mc_effi[i] -> SetMarkerSize(0.7);
        HFsum_mc_effi[i] -> SetMarkerColor(2+i);
        HFsum_mc_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        //HFsum_mc_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_HFsum_mc_effi = new TCanvas("c_HFsum_mc_effi", "c_HFsum_mc_effi", 400, 400);
    HFsum_mc_effi[1] -> Draw("elp"); 
    HFsum_mc_effi[2] -> Draw("same ep"); 
    HFsum_mc_effi[3] -> Draw("same ep"); 
    HFsum_mc_effi[4] -> Draw("same ep"); 
    t1 -> Draw();
    l2_mc -> Draw();
    
    c_HFsum_mc_effi -> SetLogx();
    c_HFsum_mc_effi -> SaveAs("pdf/HFsum_mc_effi.pdf");

    TCanvas *c_HFsum_effi = new TCanvas("c_HFsum_effi","c_HFsum_effi", 400,400);
    c_HFsum_effi -> Divide(1,2,0);
    c_HFsum_effi -> SetLogx(); 
    c_HFsum_effi -> cd(1);
    HFsum_data_effi[0] -> Draw("ep");
    HFsum_data_effi[1] -> Draw("same ep");
    HFsum_data_effi[2] -> Draw("same ep");
    HFsum_data_effi[3] -> Draw("same ep");
    HFsum_data_effi[4] -> Draw("same ep");
    
    c_HFsum_effi -> cd(2);
    HFsum_mc_effi[0] -> Draw("ep");
    HFsum_mc_effi[1] -> Draw("same ep");
    HFsum_mc_effi[2] -> Draw("same ep");
    HFsum_mc_effi[3] -> Draw("same ep");
    HFsum_mc_effi[4] -> Draw("same ep");

*/
   /* 
//======================================
// hiNpix
//======================================

    TH1D *HFNpix_data[5];
    for(int i=0; i<Ncut; i++)
    {
        HFNpix_data[i] = new TH1D(Form("HFNpix_data%d",i), ";hiNpix;Events",500,0,10000);
        HFNpix_data[i] -> SetMarkerStyle(20);
        HFNpix_data[i] -> SetMarkerSize(1.0);
        HFNpix_data[i] -> SetMarkerColor(i+1);
    }
    c_temp -> cd();
    datat_evt -> Draw("hiNpix >>+ HFNpix_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    HFNpix_data[0] = (TH1D*)gDirectory->Get("HFNpix_data0");
    datat_evt -> Draw("hiNpix >>+ HFNpix_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    HFNpix_data[1] = (TH1D*)gDirectory->Get("HFNpix_data1");
    datat_evt -> Draw("hiNpix >>+ HFNpix_data2", eventCut && "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    HFNpix_data[2] = (TH1D*)gDirectory->Get("HFNpix_data2");
    datat_evt -> Draw("hiNpix >>+ HFNpix_data3", eventCut && "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    HFNpix_data[3] = (TH1D*)gDirectory->Get("HFNpix_data3");
    datat_evt -> Draw("hiNpix >>+ HFNpix_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    HFNpix_data[4] = (TH1D*)gDirectory->Get("HFNpix_data4");

    TCanvas *c_HFNpix = new TCanvas("c_HFNpix", "c_HFNpix", 400,400);
    c_HFNpix -> SetLogy();
    c_HFNpix -> SetLogx();
    HFNpix_data[0] -> Draw("e1p");
    HFNpix_data[1] -> Draw("same&&e1p");
    HFNpix_data[2] -> Draw("same&&e1p");
    HFNpix_data[3] -> Draw("same&&e1p");
    HFNpix_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *HFNpix_data_effi[5];
    for(int i=0; i<Ncut; i++){
     //   HFNpix_data_effi[i] = new TH1D(Form("HFNpix_data_effi%d",i), "", 500,0, 5000);
        HFNpix_data_effi[i] = (TH1D*)HFNpix_data[i]->Clone(Form("HFNpix_data_effi%d",i));
        HFNpix_data_effi[i] -> SetTitle(";hiNpix;Filter Efficiency");
        if(i!=0)
            HFNpix_data_effi[i] -> Divide(HFNpix_data[i],HFNpix_data[0]);
        HFNpix_data_effi[i] -> SetMarkerStyle(20);
        HFNpix_data_effi[i] -> SetMarkerSize(1.0);
        HFNpix_data_effi[i] -> SetMarkerColor(i+1);
        HFNpix_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        HFNpix_data_effi[i] -> SetAxisRange(0.0,1000.0,"X");
    }

    TCanvas *c_HFNpix_effi = new TCanvas("c_HFNpix_effi", "c_HFNpix_effi", 400, 400);
    HFNpix_data_effi[1] -> Draw("elp"); 
    HFNpix_data_effi[2] -> Draw("same elp"); 
    HFNpix_data_effi[3] -> Draw("same elp"); 
    HFNpix_data_effi[4] -> Draw("same elp"); 
    l2 -> Draw();
    t1 -> Draw();



//=========================
// hiBin
//=========================

    TH1D *HFBin_data[5];
    for(int i=0; i<Ncut; i++)
    {
        HFBin_data[i] = new TH1D(Form("HFBin_data%d",i), ";hiBin;Events",105,0,210);
        HFBin_data[i] -> SetMarkerStyle(20);
        HFBin_data[i] -> SetMarkerSize(0.8);
        HFBin_data[i] -> SetMarkerColor(i+1);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiBin >>+ HFBin_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    HFBin_data[0] = (TH1D*)gDirectory->Get("HFBin_data0");
    datat_evt -> Draw("hiBin >>+ HFBin_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    HFBin_data[1] = (TH1D*)gDirectory->Get("HFBin_data1");
    datat_evt -> Draw("hiBin >>+ HFBin_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    HFBin_data[2] = (TH1D*)gDirectory->Get("HFBin_data2");
    datat_evt -> Draw("hiBin >>+ HFBin_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    HFBin_data[3] = (TH1D*)gDirectory->Get("HFBin_data3");
    datat_evt -> Draw("hiBin >>+ HFBin_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    HFBin_data[4] = (TH1D*)gDirectory->Get("HFBin_data4");

    TCanvas *c_HFBin = new TCanvas("c_HFBin", "c_HFBin", 400,400);
    //c_HFBin -> SetLogy();
    HFBin_data[0] -> Draw("e1p");
    HFBin_data[1] -> Draw("same&&e1p");
    HFBin_data[2] -> Draw("same&&e1p");
    HFBin_data[3] -> Draw("same&&e1p");
    HFBin_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    c_HFBin -> SetLogy();
    TH1D *HFBin_data_effi[5];
    for(int i=0; i<Ncut; i++){
     //   HFBin_data_effi[i] = new TH1D(Form("HFBin_data_effi%d",i), "", 500,0, 5000);
        HFBin_data_effi[i] = (TH1D*)HFBin_data[i]->Clone(Form("HFBin_data_effi%d",i));
        HFBin_data_effi[i] -> SetTitle(";hiBin;Filter Efficiency");
        if(i!=0)
            HFBin_data_effi[i] -> Divide(HFBin_data[i],HFBin_data[0]);
        HFBin_data_effi[i] -> SetMarkerStyle(20);
        HFBin_data_effi[i] -> SetMarkerSize(0.8);
        HFBin_data_effi[i] -> SetMarkerColor(i+1);
        HFBin_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        HFBin_data_effi[i] -> SetAxisRange(0.0,210.0,"X");
    }

    TCanvas *c_HFBin_effi = new TCanvas("c_HFBin_effi", "c_HFBin_effi", 400, 400);
    HFBin_data_effi[1] -> Draw("elp"); 
    HFBin_data_effi[2] -> Draw("same elp"); 
    HFBin_data_effi[3] -> Draw("same elp"); 
    HFBin_data_effi[4] -> Draw("same elp"); 
    l2 -> Draw();
    t1 -> Draw();
 


// ======================================
// hiZDC
// ======================================

    TH1D *hiZDC_data[5];
    for(int i=0; i<Ncut; i++){
        hiZDC_data[i] = new TH1D(Form("hiZDC_data%d",i), ";hiZDC;Events",100,0,80000);
        hiZDC_data[i] -> SetMarkerStyle(20);
        hiZDC_data[i] -> SetMarkerSize(0.8);
        hiZDC_data[i] -> SetMarkerColor(i+1);
        hiZDC_data[i] -> SetLabelSize(0.03);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiZDC>>+ hiZDC_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    hiZDC_data[0] = (TH1D*)gDirectory->Get("hiZDC_data0");
    datat_evt -> Draw("hiZDC>>+ hiZDC_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiZDC_data[1] = (TH1D*)gDirectory->Get("hiZDC_data1");
    datat_evt -> Draw("hiZDC>>+ hiZDC_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiZDC_data[2] = (TH1D*)gDirectory->Get("hiZDC_data2");
    datat_evt -> Draw("hiZDC>>+ hiZDC_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiZDC_data[3] = (TH1D*)gDirectory->Get("hiZDC_data3");
    datat_evt -> Draw("hiZDC>>+ hiZDC_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiZDC_data[4] = (TH1D*)gDirectory->Get("hiZDC_data4");

    TCanvas *c_hiZDC = new TCanvas("c_hiZDC", "c_hiZDC", 400,400);
    c_hiZDC -> SetLogy();
    hiZDC_data[0] -> Draw("e1p");
    hiZDC_data[1] -> Draw("same&&e1p");
    hiZDC_data[2] -> Draw("same&&e1p");
    hiZDC_data[3] -> Draw("same&&e1p");
    hiZDC_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *hiZDC_data_effi[5];
    for(int i=0; i<Ncut; i++){
        hiZDC_data_effi[i] = (TH1D*)hiZDC_data[i]->Clone(Form("hiZDC_data_effi%d",i));
        hiZDC_data_effi[i] -> SetTitle(";hiZDC;Filter Efficiency");
        if(i!=0)
            hiZDC_data_effi[i] -> Divide(hiZDC_data[i],hiZDC_data[0]);
        hiZDC_data_effi[i] -> SetMarkerStyle(20);
        hiZDC_data_effi[i] -> SetMarkerSize(0.8);
        hiZDC_data_effi[i] -> SetMarkerColor(i+1);
        hiZDC_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        hiZDC_data_effi[i] -> SetAxisRange(0.0,80000.0,"X");
    }

    TCanvas *c_hiZDC_effi = new TCanvas("c_hiZDC_effi", "c_hiZDC_effi", 400, 400);
    hiZDC_data_effi[1] -> Draw("elp"); 
    hiZDC_data_effi[2] -> Draw("same elp"); 
    hiZDC_data_effi[3] -> Draw("same elp"); 
    hiZDC_data_effi[4] -> Draw("same elp"); 
    t1 -> Draw();
    l1 -> Draw();



// ======================================
// hiEE
// ======================================

    TH1D *hiEE_data[5];
    for(int i=0; i<Ncut; i++){
        hiEE_data[i] = new TH1D(Form("hiEE_data%d",i), ";hiEE;Events",200,0,4000);
        hiEE_data[i] -> SetMarkerStyle(20);
        hiEE_data[i] -> SetMarkerSize(0.8);
        hiEE_data[i] -> SetMarkerColor(i+1);
        hiEE_data[i] -> SetLabelSize(0.03);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiEE>>+ hiEE_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    hiEE_data[0] = (TH1D*)gDirectory->Get("hiEE_data0");
    datat_evt -> Draw("hiEE>>+ hiEE_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiEE_data[1] = (TH1D*)gDirectory->Get("hiEE_data1");
    datat_evt -> Draw("hiEE>>+ hiEE_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiEE_data[2] = (TH1D*)gDirectory->Get("hiEE_data2");
    datat_evt -> Draw("hiEE>>+ hiEE_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiEE_data[3] = (TH1D*)gDirectory->Get("hiEE_data3");
    datat_evt -> Draw("hiEE>>+ hiEE_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiEE_data[4] = (TH1D*)gDirectory->Get("hiEE_data4");

    TCanvas *c_hiEE = new TCanvas("c_hiEE", "c_hiEE", 400,400);
    c_hiEE -> SetLogy();
    hiEE_data[0] -> Draw("e1p");
    hiEE_data[1] -> Draw("same&&e1p");
    hiEE_data[2] -> Draw("same&&e1p");
    hiEE_data[3] -> Draw("same&&e1p");
    hiEE_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *hiEE_data_effi[5];
    for(int i=0; i<Ncut; i++){
        hiEE_data_effi[i] = (TH1D*)hiEE_data[i]->Clone(Form("hiEE_data_effi%d",i));
        hiEE_data_effi[i] -> SetTitle(";hiEE;Filter Efficiency");
        if(i!=0)
            hiEE_data_effi[i] -> Divide(hiEE_data[i],hiEE_data[0]);
        hiEE_data_effi[i] -> SetMarkerStyle(20);
        hiEE_data_effi[i] -> SetMarkerSize(0.8);
        hiEE_data_effi[i] -> SetMarkerColor(i+1);
        hiEE_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        hiEE_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_hiEE_effi = new TCanvas("c_hiEE_effi", "c_hiEE_effi", 400, 400);
    hiEE_data_effi[1] -> Draw("elp"); 
    hiEE_data_effi[2] -> Draw("same elp"); 
    hiEE_data_effi[3] -> Draw("same elp"); 
    hiEE_data_effi[4] -> Draw("same elp"); 
    t1 -> Draw();
    l1 -> Draw();


// ======================================
// hiET
// ======================================

    TH1D *hiET_data[5];
    for(int i=0; i<Ncut; i++){
        hiET_data[i] = new TH1D(Form("hiET_data%d",i), ";hiET;Events",100,0,2000);
        hiET_data[i] -> SetMarkerStyle(20);
        hiET_data[i] -> SetMarkerSize(0.8);
        hiET_data[i] -> SetMarkerColor(i+1);
        hiET_data[i] -> SetLabelSize(0.03);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiET>>+ hiET_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    hiET_data[0] = (TH1D*)gDirectory->Get("hiET_data0");
    datat_evt -> Draw("hiET>>+ hiET_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiET_data[1] = (TH1D*)gDirectory->Get("hiET_data1");
    datat_evt -> Draw("hiET>>+ hiET_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiET_data[2] = (TH1D*)gDirectory->Get("hiET_data2");
    datat_evt -> Draw("hiET>>+ hiET_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiET_data[3] = (TH1D*)gDirectory->Get("hiET_data3");
    datat_evt -> Draw("hiET>>+ hiET_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiET_data[4] = (TH1D*)gDirectory->Get("hiET_data4");

    TCanvas *c_hiET = new TCanvas("c_hiET", "c_hiET", 400,400);
    c_hiET -> SetLogy();
    hiET_data[0] -> Draw("e1p");
    hiET_data[1] -> Draw("same&&e1p");
    hiET_data[2] -> Draw("same&&e1p");
    hiET_data[3] -> Draw("same&&e1p");
    hiET_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *hiET_data_effi[5];
    for(int i=0; i<Ncut; i++){
        hiET_data_effi[i] = (TH1D*)hiET_data[i]->Clone(Form("hiET_data_effi%d",i));
        hiET_data_effi[i] -> SetTitle(";hiET;Filter Efficiency");
        if(i!=0)
            hiET_data_effi[i] -> Divide(hiET_data[i],hiET_data[0]);
        hiET_data_effi[i] -> SetMarkerStyle(20);
        hiET_data_effi[i] -> SetMarkerSize(0.8);
        hiET_data_effi[i] -> SetMarkerColor(i+1);
        hiET_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        hiET_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_hiET_effi = new TCanvas("c_hiET_effi", "c_hiET_effi", 400, 400);
    hiET_data_effi[1] -> Draw("elp"); 
    hiET_data_effi[2] -> Draw("same elp"); 
    hiET_data_effi[3] -> Draw("same elp"); 
    hiET_data_effi[4] -> Draw("same elp"); 
    t1 -> Draw();
    l1 -> Draw();



// ======================================
// hiEB
// ======================================

    TH1D *hiEB_data[5];
    for(int i=0; i<Ncut; i++){
        hiEB_data[i] = new TH1D(Form("hiEB_data%d",i), ";hiEB;Events",200,0,5000);
        hiEB_data[i] -> SetMarkerStyle(20);
        hiEB_data[i] -> SetMarkerSize(0.8);
        hiEB_data[i] -> SetMarkerColor(i+1);
        hiEB_data[i] -> SetLabelSize(0.03);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiEB>>+ hiEB_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    hiEB_data[0] = (TH1D*)gDirectory->Get("hiEB_data0");
    datat_evt -> Draw("hiEB>>+ hiEB_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiEB_data[1] = (TH1D*)gDirectory->Get("hiEB_data1");
    datat_evt -> Draw("hiEB>>+ hiEB_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiEB_data[2] = (TH1D*)gDirectory->Get("hiEB_data2");
    datat_evt -> Draw("hiEB>>+ hiEB_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiEB_data[3] = (TH1D*)gDirectory->Get("hiEB_data3");
    datat_evt -> Draw("hiEB>>+ hiEB_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiEB_data[4] = (TH1D*)gDirectory->Get("hiEB_data4");

    TCanvas *c_hiEB = new TCanvas("c_hiEB", "c_hiEB", 400,400);
    c_hiEB -> SetLogy();
    hiEB_data[0] -> Draw("e1p");
    hiEB_data[1] -> Draw("same&&e1p");
    hiEB_data[2] -> Draw("same&&e1p");
    hiEB_data[3] -> Draw("same&&e1p");
    hiEB_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *hiEB_data_effi[5];
    for(int i=0; i<Ncut; i++){
        hiEB_data_effi[i] = (TH1D*)hiEB_data[i]->Clone(Form("hiEB_data_effi%d",i));
        hiEB_data_effi[i] -> SetTitle(";hiEB;Filter Efficiency");
        if(i!=0)
            hiEB_data_effi[i] -> Divide(hiEB_data[i],hiEB_data[0]);
        hiEB_data_effi[i] -> SetMarkerStyle(20);
        hiEB_data_effi[i] -> SetMarkerSize(0.8);
        hiEB_data_effi[i] -> SetMarkerColor(i+1);
        hiEB_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        hiEB_data_effi[i] -> SetAxisRange(0.0,300.0,"X");
    }

    TCanvas *c_hiEB_effi = new TCanvas("c_hiEB_effi", "c_hiEB_effi", 400, 400);
    hiEB_data_effi[1] -> Draw("elp"); 
    hiEB_data_effi[2] -> Draw("same elp"); 
    hiEB_data_effi[3] -> Draw("same elp"); 
    hiEB_data_effi[4] -> Draw("same elp"); 
    t1 -> Draw();
    l1 -> Draw();



// ======================================
// hiEvtPlanes
// ======================================

    TH1D *hiEvtPlanes_data[5];
    for(int i=0; i<Ncut; i++){
        hiEvtPlanes_data[i] = new TH1D(Form("hiEvtPlanes_data%d",i), ";hiEvtPlanes;Events",50,-2,2);
        hiEvtPlanes_data[i] -> SetMarkerStyle(20);
        hiEvtPlanes_data[i] -> SetMarkerSize(0.8);
        hiEvtPlanes_data[i] -> SetMarkerColor(i+1);
        hiEvtPlanes_data[i] -> SetLabelSize(0.03);
    }

    c_temp -> cd();
    datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data0",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1");
    hiEvtPlanes_data[0] = (TH1D*)gDirectory->Get("hiEvtPlanes_data0");
    datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data1",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pprimaryVertexFilter==1");
    hiEvtPlanes_data[1] = (TH1D*)gDirectory->Get("hiEvtPlanes_data1");
    datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data2",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phltPixelClusterShapeFilter==1");
    hiEvtPlanes_data[2] = (TH1D*)gDirectory->Get("hiEvtPlanes_data2");
    datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data3",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && phfCoincFilter3==1");
    hiEvtPlanes_data[3] = (TH1D*)gDirectory->Get("hiEvtPlanes_data3");
    datat_evt -> Draw("hiEvtPlanes>>+ hiEvtPlanes_data4",eventCut &&  "HLT_HIMinBiasHfOrBSC_v1==1 && pcollisionEventSelection==1");
    hiEvtPlanes_data[4] = (TH1D*)gDirectory->Get("hiEvtPlanes_data4");

    TCanvas *c_hiEvtPlanes = new TCanvas("c_hiEvtPlanes", "c_hiEvtPlanes", 400,400);
    c_hiEvtPlanes -> SetLogy();
    hiEvtPlanes_data[0] -> Draw("e1p");
    hiEvtPlanes_data[1] -> Draw("same&&e1p");
    hiEvtPlanes_data[2] -> Draw("same&&e1p");
    hiEvtPlanes_data[3] -> Draw("same&&e1p");
    hiEvtPlanes_data[4] -> Draw("same&&e1p");
    l1 -> Draw();

    TH1D *hiEvtPlanes_data_effi[5];
    for(int i=0; i<Ncut; i++){
        hiEvtPlanes_data_effi[i] = (TH1D*)hiEvtPlanes_data[i]->Clone(Form("hiEvtPlanes_data_effi%d",i));
        hiEvtPlanes_data_effi[i] -> SetTitle(";hiEvtPlanes;Filter Efficiency");
        if(i!=0)
            hiEvtPlanes_data_effi[i] -> Divide(hiEvtPlanes_data[i],hiEvtPlanes_data[0]);
        hiEvtPlanes_data_effi[i] -> SetMarkerStyle(20);
        hiEvtPlanes_data_effi[i] -> SetMarkerSize(0.8);
        hiEvtPlanes_data_effi[i] -> SetMarkerColor(i+1);
        hiEvtPlanes_data_effi[i] -> SetAxisRange(0.0,1.1,"Y");
        hiEvtPlanes_data_effi[i] -> SetAxisRange(-2.0,2.0,"X");
    }

    TCanvas *c_hiEvtPlanes_effi = new TCanvas("c_hiEvtPlanes_effi", "c_hiEvtPlanes_effi", 400, 400);
    hiEvtPlanes_data_effi[1] -> Draw("elp"); 
    hiEvtPlanes_data_effi[2] -> Draw("same elp"); 
    hiEvtPlanes_data_effi[3] -> Draw("same elp"); 
    hiEvtPlanes_data_effi[4] -> Draw("same elp"); 
    t1 -> Draw();
    l1 -> Draw();

*/
 } 
