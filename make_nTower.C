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

void makeHist_nTower(int run = 262694,float etThr=0.5, float eThr=0);

void make_nTower(int run = 262694){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);

    double thret[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    //double thr[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
    double thre[] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
    int nbin = sizeof(thre)/sizeof(double);
   
    for(int i=0;i<nbin;i++){
        cout << "START transverse energy threshold : " << thret[i] << endl;
        makeHist_nTower(run,thret[i],0); 
        cout << "START energy threshold : " << thre[i] << endl;
        makeHist_nTower(run,0,thre[i]); 
    }
}


void makeHist_nTower(int run, float etThr, float eThr)
{
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);
    TH1::SetDefaultSumw2();

    TCut dataCut("HLT_HIL1MinimumBiasHF1AND_v1 && phfCoincFilter3 && pprimaryVertexFilter");
    TCut mcCut("");

    // ===================================================================================
    // Get Trees from data & mc files.
    // ===================================================================================
    TFile *f[2]; // 0:data, 1:MC

    TTree *t_skim[2];
    TTree *t_hlt[2];
    TTree *t[2];
    int Entries[2];   
    const char* fname_data="";
    if(run>=262693 && run<=262735) fname_data = Form("root://eoscms.cern.ch//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HiForestPromptReco_%d.root",run);
    if(run>=262811 && run<=262816) fname_data =Form("root://eoscms.cern.ch//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias2_run%d.root",run); 

    //const char* fname_data="root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HiForest_Streamer_run262315.root";
    const char* fname_mc="/afs/cern.ch/work/y/ygo/public/centrality/Hydjet_Quenched_MinBias_5020GeV_755p4.root";
    //const char* fname_mc="/afs/cern.ch/work/y/ygo/public/centrality/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_merged_forest_6.root";

    for(int i=0;i<2;i++){ 
        if(i==0) f[0] = TFile::Open(fname_data);
        else f[1] = TFile::Open(fname_mc);
        t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
        t_hlt[i] = (TTree*) f[i] -> Get("hltanalysis/HltTree");
        t[i] = (TTree*) f[i] -> Get("rechitanalyzer/tower");
        t[i] -> AddFriend(t_hlt[i]);
        t[i] -> AddFriend(t_skim[i]);
        if(i==0) Entries[i] = t[i] -> GetEntries(dataCut);
        else Entries[i] = t[i] -> GetEntries(mcCut);
    }

    // ===============================================================================================
    // [et] Define et histograms (data/mc , pf tree/rechit tree) 
    // ===============================================================================================
    cout << "LET'S DEFINE HISTOGRAMS" << endl;
    TH1F* h1F[2]; //[0:data, 1:mc]
    for(int i=0; i<2; i++){
        h1F[i] = new TH1F(Form("h1F_sample%d",i), ";# of HF towers;Events",100,0,1000);
        h1F[i] -> SetMarkerStyle(marker[i]);
        h1F[i] -> SetMarkerSize(0.9);
        if(i==0) SetHistColor(h1F[i],2);
        else SetHistColor(h1F[i],1);
    }

    // ===============================================================================================
    // [rechit] Fill histogram by using <Draw> function in TTree.
    // ===============================================================================================
    cout << "LET'S FILL HISTOGRAMS FROM TREE" << endl;
    //********************************************************
    const char* towerCut = Form("abs(tower.eta) > 3.0 && abs(tower.eta)<5.0 && tower.e > %.1f && tower.et > %.1f",eThr, etThr);
    //const char* towerCut = Form("abs(tower.eta) > 2.87 && abs(tower.eta)<5.205 && tower.et*cosh(eta) > %.1f && tower.et > %.1f",eThr, etThr);
    
    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();

    for(int i=0; i<2; i++){
        std::clock_t    startT_draw, endT_draw;
        startT_draw = std::clock();
        if(i==0) t[i]->Draw(Form("Sum$(%s)>>+%s",towerCut,h1F[i]->GetName()), dataCut);
        else if(i==1) t[i]->Draw(Form("Sum$(%s)>>+%s",towerCut,h1F[i]->GetName()), mcCut);
        h1F[i] = (TH1F*)gDirectory->Get(h1F[i]->GetName());
        endT_draw = std::clock();
        std::cout.precision(6);      // get back to default precision
        std::cout << "DRAW finished in             : " << (endT_draw - startT_draw) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    }

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited loop" << std::endl;

    TFile* outFile = new TFile(Form("histfiles/nTower_run%d_eThr%.1f_etThr%.1f.root",run,eThr,etThr),"recreate");
    outFile->cd();
    h1F[0]->Write();
    h1F[1]->Write();
    outFile->Close();

} 
