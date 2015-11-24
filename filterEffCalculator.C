// to get the efficiency numbers for each event selection filters.
// Author Yeonju Go 
// 20 Nov 2015
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
#include <iostream>
#include "../../yjUtility.h"

void filterEffCalculator(const char* fname="root://eoscms//eos/cms//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/Merged/ExpressHiForest_run262163-262172_1.4M.root")
{
    //trig : HLT_L1MinimumBiasHF1OR_part1_v1
    const int nfilter = 6;
    const char* evtfilter[] = {"pBeamScrapingFilter","pPAprimaryVertexFilter","pHBHENoiseFilterResultProducer","phfCoincFilter","phfCoincFilter3","PAcollisionEventSelection"};
    const char* trig = "HLT_L1MinimumBiasHF1OR_part1_v1";

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    SetHistTitleStyle();
    SetyjPadStyle();
    TFile *fin = TFile::Open(fname);
    cout << "open file : " << fname << endl;
    TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
    TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
    TTree *t = (TTree*) fin -> Get("hltanalysis/HltTree");
    t-> AddFriend(t_evt);
    t-> AddFriend(t_skim);

    int minBiasN = t->GetEntries(Form("%s==1",trig));
    cout << "# of "<< trig << " = " << minBiasN << endl;
    for(int i=0;i<nfilter;i++){
        const char* tmpTrig = Form("(%s==1 && %s==1)",trig, evtfilter[i]);
        int filterN = t->GetEntries(tmpTrig);
        cout << evtfilter[i] << "::: efficiency = \t" << filterN << "\t/\t" << minBiasN << " =\t" << (double)filterN/(double)minBiasN*100 <<"\t%" << endl;
    }
}

