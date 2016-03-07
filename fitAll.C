#include "fit.C"

void fitAll(){

    static const int Nvar = 1;
    string vars[Nvar] = {"hiEE"};
//    static const int Nvar = 7;
//    string vars[Nvar] = {"hiHF","hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};

   //static const int Nvar = 6;
   // string vars[Nvar] = {"hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};

        static const int Nruns = 2;
        int run[Nruns] = {262811,262816};
    //    static const int Nruns = 5;
    //    int run[Nruns] = {262694, 262695, 262703, 262726, 262735};
    //static const int Nruns = 7;
    //int run[Nruns] = {262694, 262695, 262703, 262726, 262735, 262811, 262816};

    cout<< run[0] << endl;
    for(int i = 0; i < Nvar; ++i){
        for(int j=0; j< Nruns; ++j){
            if(0){
                //fit(vars[i].data(),"abs(zVtx[0])<10","innerVtx",262548);
                //fit(vars[i].data(),"abs(zVtx[0])>10","outerVtx",262548);
                //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])<10","hfCoinc3_innerVtx",262548);
                //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])>10","hfCoinc3_outerVtx",262548);
            }
            //            fit(vars[i].data(),"","",run[j]);
            fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter","hfCoinc3_pVtx",run[j]);
            //            fit(vars[i].data(),"phfCoincFilter3","hfCoinc3",run[j]);
            //            fit(vars[i].data(),"pprimaryVertexFilter","pVtx",run[j]);

            //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2 && pprimaryVertexFilter","hfCoinc2_pVtx",262548);
            //fit(vars[i].data(),"phfCoincFilter && pprimaryVertexFilter","hfCoinc1_pVtx",262548);
            //fit(vars[i].data(),"phfCoincFilter","hfCoinc1",262548);
            //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2","hfCoinc2",262548);
        }     

    }

}
