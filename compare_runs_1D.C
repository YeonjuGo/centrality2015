// Author : Yeonju Go
#include "../yjUtility.h"

const int colhere[] = {1,2,4,9,6,46,kOrange,kViolet,kOrange+10};
void draw_compare_runs_1D(string var="hiHF");
void compare_runs_1D(){
    TStopwatch timer;
    timer.Start();

    draw_compare_runs_1D("hiBin");
    draw_compare_runs_1D("hiHF");
    //draw_compare_runs_1D("hiNpix");
    //draw_compare_runs_1D("hiNtracks");
    //draw_compare_runs_1D("sumPtVtx");
    //draw_compare_runs_1D("hiET");
    //draw_compare_runs_1D("hiEE");
    //draw_compare_runs_1D("hiEB");
    //draw_compare_runs_1D("zVtx");//change many things !! 

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

    ////////////////////// for HIMinimumBias1
    const char* dir = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODPbPb5TeV/HIMinimumBias1/";
    const int Nrun = 32;
    const int Ncomp = 6;
    const int Nset = (int)Nrun/Ncomp + 1; 
    string run[] = 
    {"262620", "262640", "262656", "262694", "262695", "262703", 
    "262726", "262735", "262768", "262784", "262811", "262816", 
    "262818", "262834", "262837", "262893", "262921", "262988", 
    "263005", "263022", "263035", "263261", "263286", "263322",
    "263333", "263349", "263362", "263379", "263400", "263410",
    "263502", "263584"};
    const char* trigCut = "(HLT_HIL1MinimumBiasHF1AND_v1) && pcollisionEventSelection";
    //const char* trigCut = "(HLT_HIL1MinimumBiasHF1ANDExpress_v1) && pprimaryVertexFilter && phfCoincFilter3";
    const char* cap= "";
    //const char* cap= "from2620to3400";

    int nBin = 50;
    if(var=="hiBin") nBin=200;
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
        f1 = TFile::Open(Form("%sEventTree_PbPb_data_HIMinimumBias1_run%s_15Feb.root",dir,run[i].data()));
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
    }

    TCanvas* can[Nset];
    TLegend* leg[Nset];
    for(int j=0;j<Nset;j++){
        leg[j] = new TLegend(0.7,0.65,0.95,0.95);
        legStyle(l1);
        can[j] = new TCanvas(Form("can%d",j),"",500,500);
        if(var != "hiBin") can[j]->SetLogy();
        //compare between 'Ncomp' runs
        for(int i=0;i<Ncomp;i++){
            h[6*j+i]->Scale(1./nEvents[6*j+i]); // this is for all!!
            SetHistColor(h[6*j+i],colhere[j]);
            if(j==0) h[6*j+i]->DrawCopy("hist");
            else h[6*j+i]->DrawCopy("hist same");
            leg[j]->AddEntry(h[6*j+i], Form("Run %s",run[i].data()),"pl"); 
        }
        leg[j]->Draw();
        can[j]->SaveAs(Form("pdf/compareBtwRuns_run%s_%s%s.pdf",run[6*j].data(),var.data(),cap));
    }

    TCanvas* cratio[Nset];
    for(int j=0;j<Nset;j++){
        cratio[j] = new TCanvas(Form("cratio%d",j),"",500,500);
        if(var != "hiBin") can[j]->SetLogy();
        //compare between 'Ncomp' runs
        h[6*j]->Rebin(2);
        for(int i=1;i<Ncomp;i++){
            h[6*j+i]->Rebin(2); 
            h[6*j+i]->Divide(h[6*j]); 
            h[6*j+i]->SetMarkerStyle(20);
            h[6*j+i]->SetMarkerSize(0.85);
            h[6*j+i]->GetYaxis()->SetRangeUser(0.5,1.5);
            h[6*j+i]->SetTitle(Form(";%s;Run XXX / Run 26%s",var.data(),run[6*j].data()));
            if(j==1) h[6*j+i]->DrawCopy("hist");
            else h[6*j+i]->DrawCopy("hist same");
        }
        leg[j]->Draw();
        if(var == "zVtx") jumSun(-25,1,25,1);//for zVtx
        else jumSun(0,1,xmax,1);
        cratio[j]->SaveAs(Form("pdf/compareBtwRuns_ratio_run%s_%s%s.pdf",run[6*j].data(),var.data(),cap));
    }
} 

