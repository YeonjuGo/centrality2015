#include "fit.C"

void fitAll(){

    static const int Nvar = 3;
    string vars[Nvar] = {"hiNtracks","hiEE","hiHFhit"};
     //   static const int Nvar = 4;
     //   string vars[Nvar] = {"hiHF","hiNtracks","hiHFhit","hiNpix"};
     //   static const int Nvar = 7;
     //   string vars[Nvar] = {"hiHF","hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};


     int run[] = {262640, 262694, 262695, 262703, 262726, 262735, 262811, 262816};

//    int run[] =
//       {262620, 262640, 262656, 262694, 262695, 262703,
//       262726, 262735, 262768, 262784, 262811, 262816,
//       262818, 262834, 262837, 262893, 262921, 262988,
//       263005, 263022, 263035, 263261, 263286, 263322,
//       263333, 263349, 263362, 263379, 263400, 263410,
//       263502, 263584, 263604, 263614};

    //int run[] = {262620};
    //cout<< run[0] << endl;
    const int Nrun = sizeof(run)/sizeof(int);
    for(int i = 0; i < Nvar; ++i){
        for(int j=0; j< Nrun; ++j){
            if(0){
                //fit(vars[i].data(),"abs(zVtx[0])<10","innerVtx",262548);
                //fit(vars[i].data(),"abs(zVtx[0])>10","outerVtx",262548);
                //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])<10","hfCoinc3_innerVtx",262548);
                //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])>10","hfCoinc3_outerVtx",262548);
            }
            //            fit(vars[i].data(),"","",run[j]);
            //fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"epos");
            fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"CymbalEv1");
            fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"CymbalEv2");
            fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"CymbalEv3");
            fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"CymbalEv4");
            fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"CymbalEv5");
            //fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"newHydjet");
            //fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection_largeRangeMin",vars[i].data()),run[j],"epos");
            //fit(vars[i].data(),"pcollisionEventSelection",Form("%s_collisionEventSelection",vars[i].data()),run[j],"hydjet");
            //fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter",Form("%s_test",vars[i].data()),run[j]);
            //fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter",Form("%s_hfCoinc3_pVtx",vars[i].data()),run[j]);
            //            fit(vars[i].data(),"phfCoincFilter3","hfCoinc3",run[j]);
            //            fit(vars[i].data(),"pprimaryVertexFilter","pVtx",run[j]);

            //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2 && pprimaryVertexFilter","hfCoinc2_pVtx",262548);
            //fit(vars[i].data(),"phfCoincFilter && pprimaryVertexFilter","hfCoinc1_pVtx",262548);
            //fit(vars[i].data(),"phfCoincFilter","hfCoinc1",262548);
            //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2","hfCoinc2",262548);
        }     

    }

}
