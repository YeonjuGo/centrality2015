// to draw the distribution of centrality variables
// with the different event selection filters in one canvas.
//
// Author Yeonju Go
//
#include "../yjUtility.h"

//////////// modify here as you want!! ////////////////
const int ntrig =6;
//const char* trig[] = {"(L1_Centrality_ext0_5_MinimumumBiasHF1_AND_Prescl)*(L1_Centrality_ext0_5_MinimumumBiasHF1_AND==1)","(L1_Centrality_ext0_10_MinimumumBiasHF1_AND_Prescl)*(L1_Centrality_ext0_10_MinimumumBiasHF1_AND==1)","(L1_Centrality_ext30_50_MinimumumBiasHF1_AND_Prescl)*(L1_Centrality_ext30_50_MinimumumBiasHF1_AND==1)","(L1_Centrality_ext30_100_MinimumumBiasHF1_AND_Prescl)*(L1_Centrality_ext30_100_MinimumumBiasHF1_AND==1 || L1_Centrality_ext30_100_MinimumumBiasHF2_AND","L1_Centrality_ext50_100_MinimumumBiasHF1_AND || L1_Centrality_ext50_100_MinimumumBiasHF2_AND","L1_Centrality_ext70_100_MinimumumBiasHF1_AND"};
//const char* trig[] = {"L1_Centrality_ext0_5_MinimumumBiasHF1_AND","L1_Centrality_ext0_10_MinimumumBiasHF1_AND","L1_Centrality_ext30_50_MinimumumBiasHF1_AND","L1_Centrality_ext30_100_MinimumumBiasHF1_AND || L1_Centrality_ext30_100_MinimumumBiasHF2_AND","L1_Centrality_ext50_100_MinimumumBiasHF1_AND || L1_Centrality_ext50_100_MinimumumBiasHF2_AND","L1_Centrality_ext70_100_MinimumumBiasHF1_AND"};
const char* trig[] = {"L1_Centrality_ext0_5_MinimumumBiasHF1_AND","L1_Centrality_ext0_10_MinimumumBiasHF1_AND","L1_Centrality_ext30_50_MinimumumBiasHF1_AND","L1_Centrality_ext30_100_MinimumumBiasHF2_AND","L1_Centrality_ext50_100_MinimumumBiasHF2_AND","L1_Centrality_ext70_100_MinimumumBiasHF1_AND"};
const char* legSt[] = {"Cent0-5","Cent0-10","Cent30-50","Cent30-100","Cent50-100","Cent70-100"};
void Get1DPlots(TTree* t_evt=0, TString v1="hiHF",int xbin=200, double xmin=0, double xmax=4500, TCut cut="",const char* cap="", bool isPassed=1);

void centTriggerTurnOn_prescl(
        const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/HIRun2015/HIExpressPhysics/Merged/HIForestExpress_run263261.root",
        const char* type="Run263261_express" 
        )
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
    int nBin = 50;
    Get1DPlots(t_evt, "hiBin",100,0,200,"",type);

}

void Get1DPlots(TTree* t_evt, TString v1, int xbin, double xmin, double xmax, TCut cut,const char* cap, bool isPassed)
{

    TCut totcut[ntrig];
    for(int i=0; i<ntrig;i++){
        totcut[i] = cut && Form("(%s_Prescl)*(%s==%d)",trig[i],trig[i],(int)isPassed);
        //totcut[i] = cut && Form("%s==%d",trig[i],(int)isPassed);
    }

    //TCanvas *c_temp= new TCanvas(Form("c_temp_%s",cap), "", 300,300);
    TH1D *h1D[ntrig];
    TLegend* l1 = new TLegend(0.20, 0.7, 0.95, 0.95,"L1_Centrality_ext*_*_MinimumumBiasHF1(2)_AND"); 
    legStyle(l1);

    for(int i=0; i<ntrig; i++){
        h1D[i] = new TH1D(Form("h1D_%d",i), Form(";%s;Events", v1.Data()), xbin, xmin,xmax );
        //if(i!=0) h1D[i]->SetMarkerStyle(24+i); 
        //h1D[i]->SetMarkerSize(0.8);
        //SetHistColor(h1D[i],col[i]);
        t_evt->Draw(Form("%s>>+%s",v1.Data(), h1D[i]->GetName()), totcut[i]);
        h1D[i]=(TH1D*)gDirectory->Get(h1D[i]->GetName());
        h1D[i] -> SetFillStyle(3003); 
        h1D[i] -> SetFillColor(col[i]); 
        h1D[i] -> SetLineColor(col[i]); 
        h1D[i] -> GetYaxis()->SetRangeUser(0.5,100000); 
        l1->AddEntry(h1D[i], legSt[i],"f"); 
        //l1->AddEntry(h1D[i], trig[i]); 
    }

    TCanvas *c_tot = new TCanvas(Form("c_tot_%s_%s",v1.Data(),cap), "c_tot", 500,500);
    for(int i=0;i<ntrig;i++){
        c_tot->cd(1);
        if(i==0) h1D[i] -> DrawCopy("hist");
        else h1D[i] -> DrawCopy("hist same");
        l1 -> Draw();
    }
    gPad->SetLogy();
    double height = 1000; 
    //double height = h1D[2]->GetMaximum();
    drawText(Form("%s",cap),0.60,0.8,kBlack,18);
    jumSun(10,0.5,10,height);
    jumSun(20,0.5,20,height);
    jumSun(60,0.5,60,height);
    jumSun(100,0.5,100,height);
    jumSun(140,0.5,140,height);
    c_tot->SaveAs(Form("trig_prescl_pdf/centTrig_prescl_%s_%s.pdf",v1.Data(),cap));

    for(int i=0; i<ntrig; i++){
        delete h1D[i];
    }
}
