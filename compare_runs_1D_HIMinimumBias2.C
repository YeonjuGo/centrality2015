// Author : Yeonju Go
#include "../yjUtility.h"
#include <string>
#include <iostream>

const int colhere[] = {1,2,4,kYellow+2,6,kGreen+2,46,kOrange,kViolet,kOrange+10};
void draw_compare_runs_1D(string var, const char* trigCut, const char* cap, const char* dataset, bool presclOn);
//void draw_compare_runs_1D(string var, const char* trigCut, const char* cap, const char* dataset, bool presclOn = 1);
void compare_runs_1D_HIMinimumBias2(
        const char* trigCut = "pcollisionEventSelection",
        const char* cap= "_rightTrigger_pcollisionEvnetSelection",
        const char* dataset="HIMinimumBias2",
        bool presclOn=1)
{
    TStopwatch timer;
    timer.Start();

    draw_compare_runs_1D("hiBin",trigCut,cap,dataset,1);
    draw_compare_runs_1D("hiHF",trigCut,cap,dataset,1);
    //draw_compare_runs_1D("hiNpix",trigCut,cap,dataset);
    //draw_compare_runs_1D("hiNtracks",trigCut,cap,dataset);
    //draw_compare_runs_1D("sumPtVtx",trigCut,cap,dataset);
    //draw_compare_runs_1D("hiET",trigCut,cap,dataset);
    //draw_compare_runs_1D("hiEE",trigCut,cap,dataset);
    //draw_compare_runs_1D("hiEB",trigCut,cap,dataset);
    //draw_compare_runs_1D("zVtx",trigCut,cap,dataset);//change many things !! 

    timer.Stop();
    cout<<"Macro finished: "<<endl;
    cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
    cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
}

void draw_compare_runs_1D(string var, const char* trigCut, const char* cap, const char* dataset, bool presclOn){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2(); 

    ////////////////////// for HIMinimumBias1
    const char* dir = "/home/goyeonju/CMS/Files/centrality/"; //KOREA UNIV.
    //const char* dir = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODPbPb5TeV/HIMinimumBias1/";//EOS

    //string run[] =
    //{"262620", "262640", "262656", "262694", "262695", "262703"};
    
       int run[] = 
       {262620, 262640, 262656, 262694, 262695, 262703, 
       262726, 262735, 262768, 262784, 262811, 262816, 
       262818, 262834, 262837, 262893, 262921, 262988, 
       263005, 263022, 263035, 263261, 263286, 263322,
       263333, 263349, 263362, 263379, 263400, 263410,
       263502, 263584, 263604, 263614};
       
    const int Nrun = sizeof(run)/sizeof(int);
    const int Ncomp = 6;
    int Nset =0;
    if(Nrun%Ncomp!=0) Nset = (int)Nrun/Ncomp + 1;
    else Nset = (int)Nrun/Ncomp;

    int nBin = 50;
    if(var=="hiBin") nBin=50;
    //if(var=="hiBin") nBin=200;
    double xmax = 6000;
    if(var == "hiNpix") xmax = 50000;
    if(var == "hiNtracks" || var == "sumPtVtx") xmax = 3500;
    if(var == "hiBin") xmax=200;
    if(var == "hiEE") xmax = 3000;
    if(var == "hiEB") xmax = 4000;

    int nEvents[Nrun];

    TH1D* h[Nrun];
    for(int i=0; i<Nrun; ++i){
        TFile* f1;
        f1 = TFile::Open(Form("%sEventTree_PbPb_data_%s_run%d_15Feb.root",dir,dataset,run[i]));
        if(f1->IsZombie()) { cout << run[i] << " doesn't exist!! " << endl; continue;}
        cout << "Open file : " << f1->GetName() << endl;
        TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
        TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
        TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
        TTree* t1_track= (TTree*) f1 -> Get("anaTrack/trackTree");
        t1->AddFriend(t1_skim);
        t1->AddFriend(t1_hlt);
        t1->AddFriend(t1_track);
        nEvents[i] = t1->GetEntries(trigCut);
        cout << run[i]<< " # of events : " << nEvents[i] << endl;
        if(var == "zVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, -25, 25); // for zVtx 
        else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
        //draw and considering Prescl
        if(run[i]<263155) t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),Form("HLT_HIL1MinimumBiasHF1AND_v1 && %s",trigCut));
        else t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),Form("HLT_HIL1MinimumBiasHF2AND_part1_v1 && %s",trigCut));
    }

    TCanvas* can[Nset];
    TLegend* leg[Nset];
    for(int j=0;j<Nset;j++){
        leg[j] = new TLegend(0.7,0.65,0.95,0.95);
        legStyle(leg[j]);
        can[j] = new TCanvas(Form("can%d",j),"",500,500);
        if(var != "hiBin") can[j]->SetLogy();
        //compare between 'Ncomp' runs
        for(int i=0;i<Ncomp;i++){
            if((Ncomp*j+i)>=Nrun) continue;
            h[Ncomp*j+i]->Scale(1./nEvents[Ncomp*j+i],"width"); // this is for all!!
            SetHistColor(h[Ncomp*j+i],colhere[i]);
            if(i==0) h[Ncomp*j+i]->DrawCopy("hist e");
            else h[Ncomp*j+i]->DrawCopy("hist same e");
            leg[j]->AddEntry(h[Ncomp*j+i], Form("Run %d",run[Ncomp*j+i]),"el"); 
        }
        leg[j]->Draw();
        drawText(dataset,0.2,0.23);
        drawText(trigCut,0.2,0.2);
        can[j]->SaveAs(Form("pdf/compareBtwRuns_%s_run%d_%s%s.pdf",dataset,run[Ncomp*j],var.data(),cap));
    }
    cout << "c"<<endl;
    TCanvas* cratio[Nset];
    for(int j=0;j<Nset;j++){
        cratio[j] = new TCanvas(Form("cratio%d",j),"",500,500);
        if(var != "hiBin") { 
            can[j]->SetLogy();
            h[Ncomp*j]->Rebin(2);
        }
        //compare between 'Ncomp' runs
        for(int i=1;i<Ncomp;i++){
            if((Ncomp*j+i)>=Nrun) continue;
            if(var != "hiBin") h[Ncomp*j+i]->Rebin(2); 
            h[Ncomp*j+i]->Divide(h[Ncomp*j]); 
            h[Ncomp*j+i]->SetMarkerStyle(20);
            h[Ncomp*j+i]->SetMarkerSize(0.85);
            h[Ncomp*j+i]->GetYaxis()->SetRangeUser(0.5,1.5);
            h[Ncomp*j+i]->SetTitle(Form(";%s;Run XXX / Run 26%d",var.data(),run[Ncomp*j]));
            if(i==1) h[Ncomp*j+i]->DrawCopy("hist e");
            else h[Ncomp*j+i]->DrawCopy("hist same e");
        }
        leg[j]->Draw();
        drawText(dataset,0.2,0.23);
        drawText(trigCut,0.2,0.2);
        if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
        else jumSun(0,1,xmax,1);
        cratio[j]->SaveAs(Form("pdf/compareBtwRuns_%s_ratio_run%d_%s%s.pdf",dataset,run[Ncomp*j],var.data(),cap));
    }
} 

