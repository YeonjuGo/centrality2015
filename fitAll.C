#include "fit.C"

void fitAll(){

    static const int Nvar = 1;
    string vars[Nvar] = {"hiHF"};
    //    static const int Nvar = 7;
    //    string vars[Nvar] = {"hiHF","hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};

    // static const int Nruns = 7;
    // int run[Nruns] = {262694, 262695, 262703, 262726, 262735, 262811, 262816};


    int run[] =
       {262620, 262640, 262656, 262694, 262695, 262703,
       262726, 262735, 262768, 262784, 262811, 262816,
       262818, 262834, 262837, 262893, 262921, 262988,
       263005, 263022, 263035, 263261, 263286, 263322,
       263333, 263349, 263362, 263379, 263400, 263410,
       263502, 263584, 263604, 263614};

    //int run[] = {262620};
    //cout<< run[0] << endl;
    const int Nrun = sizeof(run)/sizeof(int);
    const int Ncomp = 6;
    int Nset =0;
    if(Nrun%Ncomp!=0) Nset = (int)Nrun/Ncomp + 1;
    else Nset = (int)Nrun/Ncomp;

    for(int i = 0; i < Nvar; ++i){
        for(int j=0; j< Nrun; ++j){
            if(0){
                //fit(vars[i].data(),"abs(zVtx[0])<10","innerVtx",262548);
                //fit(vars[i].data(),"abs(zVtx[0])>10","outerVtx",262548);
                //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])<10","hfCoinc3_innerVtx",262548);
                //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])>10","hfCoinc3_outerVtx",262548);
            }
            //            fit(vars[i].data(),"","",run[j]);
            fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter",Form("%s_hfCoinc3_pVtx",vars[i].data()),run[j]);
            //            fit(vars[i].data(),"phfCoincFilter3","hfCoinc3",run[j]);
            //            fit(vars[i].data(),"pprimaryVertexFilter","pVtx",run[j]);

            //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2 && pprimaryVertexFilter","hfCoinc2_pVtx",262548);
            //fit(vars[i].data(),"phfCoincFilter && pprimaryVertexFilter","hfCoinc1_pVtx",262548);
            //fit(vars[i].data(),"phfCoincFilter","hfCoinc1",262548);
            //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2","hfCoinc2",262548);
        }     

    }

}
