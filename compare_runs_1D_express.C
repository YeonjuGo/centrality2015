// Author : Yeonju Go
#include "../yjUtility.h"

const int colhere[] = {1,2,4,9,6,46,kOrange,kViolet,kOrange+10};
void draw_compare_runs_1D(string var="hiHF");
void compare_runs_1D_express(){
    TStopwatch timer;
    timer.Start();
    //draw_compare_runs_1D("zVtx");//change many things !! 
    draw_compare_runs_1D("hiHF");
    //draw_compare_runs_1D("hiNpix");
    //draw_compare_runs_1D("hiNtracks");
    //draw_compare_runs_1D("sumPtVtx");
    draw_compare_runs_1D("hiBin");
    //draw_compare_runs_1D("hiEE");
    //draw_compare_runs_1D("hiEB");

    timer.Stop();
    cout<<"Macro finished: "<<endl;
    cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
    cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
}

void draw_compare_runs_1D(string var){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2(); 
#if 0 
    ///////////////////// for larger sample
    const char* dir = "root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/HIMinimumBias2/Merged/";
    const int Nrun = 7;
    string run[] = {"2695","2548","2620","2656","2694","2811","2816","3031"};
    //string run[] = {"2695","2548","2620","2656","2694","2811","2816","3033"};
    //string run[] = {"695","548","620","656","694","811","816"};
    //const char* run[] = {"695","548","620","656","694"};
    //const char* run[] = {"548","620","656","694","695"};
    //const char* trigCut = "HLT_HIL1MinimumBiasHF1AND_v1";
    TCut trigCut = "HLT_HIL1MinimumBiasHF1AND_v1 && pprimaryVertexFilter && phfCoincFilter3";
    const char* cap= "";

    double trkPtBin[] = {0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0,12.0,16.0,20.0};
    int NtrkPtBin = sizeof(trkPtBin)/sizeof(double)-1;
#endif
#if 1 
    ////////////////////// for express
    const char* dir = "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/";
    const int Nrun = 6;
    string run[] = {"2620","2837","3286","3362","3379","3400"};
    //string run[] = {"2620","2837","2921","3286","3362","3379"};
    //string run[] = {"2620","2837","2921","3286","3362","3379"};
    //string run[] = {"3400","3410","3502","3584","3604","3614","3685"};
    //string run[] = {"2620","2837","2921","3286","3362","3379","3400","3410","3502","3584","3604","3614","3685"};
    const char* trigCut = "(HLT_HIL1MinimumBiasHF1ANDExpress_v1) && pprimaryVertexFilter && phfCoincFilter3";
    //const char* trigCut = "(HLT_HIL1MinimumBiasHF1ANDExpress_v1||HLT_HIL1MinimumBiasHF2ANDExpress_v1) && pprimaryVertexFilter && phfCoincFilter3";
    const char* cap= "from2620to3400";
    //const char* cap= "from3400to3685";
#endif
    int nBin = 50;
    if(var=="hiBin") nBin=200;
    double xmax = 6000;
    if(var == "hiNpix") xmax = 50000;
    if(var == "hiNtracks" || var == "sumPtVtx") xmax = 3500;
    if(var=="hiBin") xmax=200;
    if(var == "hiEE") xmax = 3000;
    if(var == "hiEB") xmax = 4000;

    int nEvents[Nrun];

    TH1D* h[Nrun];
    TLegend* l1 = new TLegend(0.7,0.65,0.95,0.95);
    legStyle(l1);
    for(int i=0; i<Nrun; ++i){
        TFile* f1;
        if(run[i]=="2837") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262837_preHBHE.root"));
        else if(run[i]=="2620") f1 = TFile::Open(Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run262620-v6.root"));
        else f1 = TFile::Open(Form("%sHIForestExpress_run26%s.root",dir,run[i].data()));
        if(f1->IsZombie()) { cout << run[i].data() << " doesn't exist!! " << endl; continue;}
        cout << "Open file : " << f1->GetName() << endl;
        TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
        TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
        TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
        TTree* t1_track= (TTree*) f1 -> Get("anaTrack/trackTree");
        t1->AddFriend(t1_skim);
        t1->AddFriend(t1_hlt);
        t1->AddFriend(t1_track);
        nEvents[i] = t1->GetEntries(trigCut);
        cout << run[i].data() << " # of events : " << nEvents[i] << endl;
        if(var == "zVtx") h[i] = new TH1D(Form("h%d",i), Form(";%s (cm);Event Fraction",var.data()), nBin, -25, 25); // for zVtx 
        else h[i] = new TH1D(Form("h%d",i), Form(";%s;Event Fraction",var.data()), nBin, 0, xmax);
        t1->Draw(Form("%s>>%s",var.data(),h[i]->GetName()),trigCut);
        /*
        delete f1;
        delete t1;
        delete t1_skim;
        delete t1_hlt;
        delete t1_track;
        */
    }

    TCanvas* c1 = new TCanvas("c2","",500,500);
    if(var != "hiBin") c1->SetLogy();
    for(int i=0; i<Nrun; ++i){
        h[i]->Scale(1./nEvents[i]); // this is for all!! 
        //h[i]->Scale(1./h[i]->Integral("width"));
        SetHistColor(h[i],colhere[i]);
        if(i==0) h[i]->DrawCopy("hist");
        else h[i]->DrawCopy("hist same");
        l1->AddEntry(h[i], Form("Run 26%s",run[i].data()),"pl"); 
    }
    l1->Draw();
    c1->SaveAs(Form("pdf/compareBtwRuns_run26%s_%s%s.pdf",run[0].data(),var.data(),cap));

    h[0]->Rebin(2);
    TCanvas* c2 = new TCanvas("c3","", 500, 500); 
    for(int i=1; i<Nrun; ++i){
        h[i]->Rebin(2);
        h[i]->Divide(h[0]);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerSize(0.85);
        h[i]->GetYaxis()->SetRangeUser(0.5,1.3);
        h[i]->SetTitle(Form(";%s;Run XXX / Run 26%s",var.data(),run[0].data()));
        if(i==1) h[i]->DrawCopy("pe");
        else h[i]->DrawCopy("pe same");
    }
    if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
    else jumSun(0,1,xmax,1);
    c2->SaveAs(Form("pdf/compareBtwRuns_ratio_run26%s_%s%s.pdf",run[0].data(),var.data(),cap));
} 

