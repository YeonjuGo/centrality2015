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
const int rebinN = 5;
float normHistHere(TH1* hNom=0, TH1* hDen=0, TH1* hRatio=0, double cut_i=700, double cut_f=900);
void nTower_plots(double *thr, int nbin, const char* var = "e", int run=262694, const char* cap="pVtx_hfCoinc3", float norm_i = 150, float norm_f=700);

void draw_nTower(){

    double etthr[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    double ethr[] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0}; 
    int nbin_et = sizeof(etthr)/sizeof(double);
    int nbin_e = sizeof(ethr)/sizeof(double);

    nTower_plots(ethr, nbin_e,"e", 262694, "",150, 700);
    nTower_plots(etthr, nbin_et,"et", 262694, "", 150, 700);
}

void nTower_plots(double *thr, int nbin, const char* var, int run, const char* cap, float norm_i, float norm_f){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);

    double thrArr[nbin+1];
    for(int i=0;i<nbin+1;i++){
        double binwidth = (thr[1]-thr[0])/2.;
        if(i!=nbin) thrArr[i] = thr[i]-binwidth;
        else thrArr[i] = thr[nbin-1]+binwidth;
    }

    double ratio[nbin];
    TH1D* h1D_e = new TH1D(Form("h1D_%s_ratio",var),Form(";(# of HF tower) %s threshold (GeV);DATA/MC",var), nbin,thrArr);
    h1D_e->SetTitleOffset(1.5, "Y");
    h1D_e->SetNdivisions(505);
    for(int i=0;i<nbin;i++){
        cout << "START " << var << " threshold : " << thr[i] << endl;
        TFile* f = new TFile(Form("histfiles/nTower_run%d_%sThr%.1f%s.root",run,var,thr[i],cap));
        TH1F* hdata = (TH1F*) f->Get("h1F_sample0");
        TH1F* hmc = (TH1F*) f->Get("h1F_sample1");
        TH1F* hratio = (TH1F*) hdata->Clone("h1F_ratio");
        ratio[i] = normHistHere(hdata, hmc, hratio, norm_i, norm_f);
        cout << i << "th ratio = " << ratio[i] << endl;
        h1D_e->SetBinContent(i+1,ratio[i]);

        TCanvas *c = new TCanvas("c","c", 500,1000);
        c->Divide(1,2);
        c->cd(1);
        gPad->SetLogy();
        hdata->SetNdivisions(505);
        hmc->SetNdivisions(505);
        Double_t range = cleverRange(hdata,hmc);
        hdata->GetYaxis()->SetRangeUser(0.1,range);
        hmc->GetYaxis()->SetRangeUser(0.1,range);
        jumSun(norm_i, 0.1, norm_i, range);
        jumSun(norm_f, 0.1, norm_f, range);
        drawText(Form("tower %s > %.1f",var,thr[i]), 0.46,0.88,18);
        hdata->DrawCopy();
        hmc->DrawCopy("hist same");
        c->cd(2);    
        hratio->Rebin(rebinN);
        hratio->Scale(1./rebinN); 
        hratio->GetYaxis()->SetTitle("DATA/MC");
        hratio->GetYaxis()->SetRangeUser(0,2);
        hratio->DrawCopy();
        jumSun(0, 1, 1000, 1);
        drawText(Form("DATA/MC = %.3f", ratio[i]), 0.26,0.88,18);
        c->SaveAs(Form("pdf/nTower_run%d_%sThr%.1f_normRange%dto%d%s.pdf",run,var,thr[i],(int)norm_i,(int)norm_f,cap));
    }
    TCanvas* c_e = new TCanvas("c_e", "",300,300);
    h1D_e->Draw("hist");
    c_e->SaveAs(Form("pdf/run%d_%sThr_vs_ratio_distribution%s.pdf",run,var,cap));
}


float normHistHere(TH1* hNom, TH1* hDen, TH1* hRatio, double cut_i, double cut_f){
    int cutBinFrom = hDen->FindBin(cut_i);
    int cutBinTo = hDen->FindBin(cut_f);
    int lastBin = hDen->FindBin(990);

    TH1::SetDefaultSumw2();
    hDen->Scale(hNom->Integral(cutBinFrom,cutBinTo)/hDen->Integral(cutBinFrom,cutBinTo));
    hRatio->Divide(hNom,hDen);

    //cout << "total integral of DATA in the whole range: " << hData->Integral()<< endl;
    //cout << "total integral of MC in the whole range: " << hMC->Integral()<< endl;
    //cout << "full efficiency in the whole range: " <<  hData->Integral()/hMC->Integral()<< endl;

    double nomInt = hNom->Integral(0,cutBinTo)+hNom->Integral(cutBinTo+1,lastBin);
    double denInt = hDen->Integral(0,cutBinTo)+hNom->Integral(cutBinTo+1,lastBin);
    return nomInt/denInt;
} 
