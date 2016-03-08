#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "../yjUtility.h"

using namespace std;

static double rangeMin = 200;
static string date = "20160307";
static bool doUPC = 0;
static bool express = 0;

double chi2(TH1* h1, TH1* h2){


    double c = 0;
    int ndof = 0;

    for(int i = h1->GetXaxis()->FindBin(rangeMin); i < h1->GetNbinsX()+1; ++i){
        double y1 = h1->GetBinContent(i);
        double y2 = h2->GetBinContent(i);
        double e1 = h1->GetBinError(i);
        double e2 = h2->GetBinError(i);

        double dy = y1-y2;
        double de2 = e1*e1+e2*e2;
        if(de2 > 0) { 
            c += dy*dy/de2;
            ndof += 1;
        }
    }

    return c/ndof;
}


TH1* scale(string var, TTree* nt, double s = 1, TCut cut = ""){

    int nbins = 100;

    double xmin = 0;
    double xmax = 8000;

    if(var == "hiET"){
        xmax = 2000;
        rangeMin = 100;
    }

    if(var == "hiEB"){
        xmax = 5000;
        rangeMin = 200;
    }

    if(var == "hiEE"){
        xmax = 4000;
        rangeMin = 200;
    }

    if(var == "hiHF"){
        xmax = 6000;
        rangeMin = 200;
    }

    if(var == "hiHFhit"){
        xmax = 200000;
        rangeMin = 10000;
    }

    if(var == "hiNpix"){
        xmax = 50000;
        rangeMin = 1500;
    }

    if(var == "hiNtracks"){
        xmax = 4000;
        rangeMin = 200;
    }



    TH1D* h = new TH1D("h",Form("h;%s [GeV];event fraction",var.data()),nbins,xmin,xmax);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();

    nt->Draw(Form("%f*%s>>h",s,var.data()),cut);

    double rangeIntegral = h->Integral(h->GetXaxis()->FindBin(rangeMin),h->GetNbinsX()+1);
    if(rangeIntegral > 0){
        h->Scale(1./rangeIntegral);
    }

    return h;

}


void fit(string var = "hiHFhit", TCut dataCut = "", string cutname = "", int run = 262694){

    if(cutname == "") cutname = (const char*)dataCut;

    TH1::SetDefaultSumw2();

    const char* dir = "/home/goyeonju/CMS/Files/centrality/"; //KOREA UNIV.
    TFile *infData, *infMC, *infUPC;
    TTree *tref, *t, *tupc;
    infData = TFile::Open(Form("%sEventTree_PbPb_data_HIMinimumBias2_run%d_15Feb.root",dir,run));
    //    if(run>=262693 && run<=262735) infData = TFile::Open(Form("root://eoscms.cern.ch//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HiForestPromptReco_%d.root",run));
    //    if(run>=262811 && run<=262816) infData = TFile::Open(Form("root://eoscms.cern.ch//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HIMinimumBias2_run%d.root",run));

    cout << infData->GetName() << endl;
    infMC = TFile::Open("/home/goyeonju/CMS/Files/centrality/Hydjet_Quenched_MinBias_5020GeV_758p2_FOREST-v31.root");
    //infMC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/mverweij/jetsPbPb/Run2Prep/Hydjet_Quenched_MinBias_5020GeV_750/crab_Run2_HydjetMB/151109_130226/HiForestMerged.root");

    if(doUPC) infUPC = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/chflores/Foresting_RunPrep2015/STARLIGHTProd/HiForest_Starlight_Merge19112015.root");

    tref = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
    t = (TTree*)infMC->Get("hiEvtAnalyzer/HiTree");
    if(doUPC) tupc = (TTree*)infUPC->Get("hiEvtAnalyzer/HiTree");

    tref->AddFriend("hltanalysis/HltTree");
    tref->AddFriend("skimanalysis/HltTree");
    tref->AddFriend("anaTrack/trackTree");

    t->AddFriend("hltanalysis/HltTree");
    t->AddFriend("skimanalysis/HltTree");

    if(doUPC){
        tupc->AddFriend("hltanalysis/HltTree");
        tupc->AddFriend("skimanalysis/HltTree");
    }

    TH1* hupc, *href;

    if(doUPC){
        //WARNING : No trigger for UPC
        hupc = scale(var,tupc,1.,dataCut);
        hupc->SetLineColor(4);
        hupc->SetName("hupc");
    }

    TCut trigger("HLT_HIL1MinimumBiasHF1AND_v1");
    //if(run==262694) trigger="HLT_HIL1MinimumBiasHF1AND_v1";
    //TCut trigger("HLT_HIL1MinimumBiasHF2AND_v1");
    //TCut trigger("HLT_HIL1MinimumBiasHF2AND_part2_v1");
    if(express) trigger = "HLT_HIL1MinimumBiasHF1ANDExpress_v1";

    href = scale(var,tref,1.,trigger&&dataCut);

    cout<<"DATA ENTRIES : "<<tref->GetEntries(trigger&&dataCut)<<endl;

    static const int N = 40;
    TH1* h[N];
    double chi2s[N];
    double s[N];
    int ibest = 0;
    TGraphErrors* g = new TGraphErrors(N);

    double variation = 0.2;
    double scaleMin = 0.75;

    if(var == "hiHF" || var == "hiHFhit"){
        scaleMin = 0.50;
        // scaleMin = 0.25;
    }

    if(var == "hiNpix"){
        scaleMin = 1.0;
    }

    if(var == "hiEE" || var == "hiEB" || var == "hiET"){
        scaleMin = 0.66;
        //scaleMin = 0.85;
    }



    for(int i = 0; i < N; ++i){  

        s[i] = scaleMin+(variation/N)*i;
        h[i] = scale(var,t,s[i]);
        h[i]->SetLineColor(2);
        chi2s[i] = chi2(h[i], href);
        g->SetPoint(i,s[i],chi2s[i]);
        if(chi2s[i]< chi2s[ibest]) ibest = i;
        //      cout<<"for scale "<<s[i]<<", chi2 is : "<<chi2s[i]<<endl;

    }

    double eff = href->Integral()/h[ibest]->Integral();

    cout<<"the best scale is : "<<s[ibest]<<" i = "<<ibest<<endl;   
    cout<<"efficiency+contamination is : "<<eff<<endl;

    ofstream outputf(Form("centEff_fitResults_%s.txt",cutname.data()), ios::app);
    outputf << run << "\t" << setprecision(4) <<eff*100 << "\t" << chi2s[ibest] << endl;
    outputf.close();

    //   href->Draw();
    //   h->Draw("same");

    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    g->GetXaxis()->SetTitle("scale factor");
    g->GetYaxis()->SetTitle("reduced #chi^2 (normalized by ndf)");
    //g->GetYaxis()->SetTitle("#chi^2 (unnormalized)");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->SetMarkerStyle(20);
    g->Draw("Ap");

    c1->Print(Form("figures/figure_%s_%s_%s_%s_run%d_%s.png","chi2",var.data(),(const char*)trigger,cutname.data(),run,date.data()));

    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    c2->SetLogy();

    h[ibest]->SetMaximum(10.*h[ibest]->GetMaximum());
    h[ibest]->Draw("hist");
    href->Draw("same");
    TLine* vertical = new TLine(rangeMin, h[ibest]->GetMinimum(), rangeMin, h[ibest]->GetMaximum());
    vertical->SetLineStyle(2);
    vertical->Draw("same");

    TLegend * leg1 = new TLegend(0.4,0.6,0.95,0.9);
    leg1->SetFillStyle(1);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);

    leg1->SetTextSize(0.025);
    leg1->AddEntry(href,Form("Run %d",run),"");
    leg1->AddEntry(href,"data","p");
    leg1->AddEntry(href,trigger,"");
    leg1->AddEntry(href,dataCut,"");

    leg1->AddEntry(h[ibest],"Hydjet based fit","l");
    leg1->AddEntry(h[ibest],Form("Hydjet scaling : %f",s[ibest]),"");
    leg1->AddEntry(h[ibest],Form("Eff+Contam : %f",eff),"");
    leg1->AddEntry(h[ibest],Form("#chi^{2}/ndf : %f",chi2s[ibest]),"");
    leg1->Draw();

    c2->Print(Form("figures/figure_%s_%s_%s_%s_run%d_%s.png","fit",var.data(),(const char*)trigger,cutname.data(),run,date.data()));

    TCanvas* c3 = new TCanvas("c3","c3",600,600);
    c3->SetLogx();
    TH1D* hdiff = (TH1D*)href->Clone("hdiff");
    hdiff->SetTitle(Form(";%s [GeV];data - fit",var.data()));
    hdiff->Add(h[ibest],-1);
    hdiff->Draw();
    TLine* horizontal = new TLine(0,0,8000,0);
    horizontal->SetLineStyle(2);
    horizontal->Draw("same");

    double supc = 0;
    if(doUPC){
        supc = hdiff->Integral(1,5)/hupc->Integral(1,5);
        hupc->Scale(supc);
        hupc->Draw("hist same");
    }

    TLegend * leg2 = new TLegend(0.4,0.6,0.95,0.9);
    leg2->SetFillStyle(1);
    leg2->SetFillColor(0);
    leg2->SetBorderSize(0);

    leg2->SetTextSize(0.025);
    leg2->AddEntry(hdiff,"data - fit","p");
    leg2->AddEntry(href,trigger,"");
    leg2->AddEntry(href,dataCut,"");

    if(doUPC){
        leg2->AddEntry(hupc,"Starlight","l");
        leg2->AddEntry(hupc,Form("Vertical scaling : %f",supc),"");
    }
    leg2->Draw();

    c3->Print(Form("figures/figure_%s_%s_%s_%s_run%d_%s.png","diff",var.data(),(const char*)trigger,cutname.data(),run,date.data()));

}

