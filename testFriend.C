

{

   //   TFile* inf0 = TFile::Open("root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HiForestStreamer_Run262314-262315.root");
   TFile* inf0 = new TFile("/data/velicanu/store/group/phys_heavyions/velicanu/forest/Run2015E/HIExpressPhysics/Merged/HiForest_Streamer_run262315.root");
   TFile* inf1 = new TFile("friendly_centrality_HiForest_Streamer_run262315.root");

   TTree* tref = (TTree*)inf0->Get("hiEvtAnalyzer/HiTree");
   TTree* t = (TTree*)inf1->Get("anaCentrality");

   t->AddFriend(tref);

   TH2D* h = new TH2D("h",";hiBin;anaBin",200,0,200,200,0,200);
   t->Draw("anaBin:hiBin>>h","","colz");

}

