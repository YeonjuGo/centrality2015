#include "fit.C"

void fitAll(){

   static const int Nvar = 1;
   string vars[Nvar] = {"hiHF"};

//   static const int Nvar = 2;
//   string vars[Nvar] = {"hiNpix","hiEE"};
//
//   static const int Nvar = 7;
//   string vars[Nvar] = {"hiHF","hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};
//


   for(int i = 0; i < Nvar; ++i){

      if(0){
      //fit(vars[i].data(),"abs(zVtx[0])<10","innerVtx",262548);
      //fit(vars[i].data(),"abs(zVtx[0])>10","outerVtx",262548);
      
      //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])<10","hfCoinc3_innerVtx",262548);
      //fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])>10","hfCoinc3_outerVtx",262548);

//      fit(vars[i].data(),"","",262548);
//      fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter","hfCoinc3_pVtx",262548);
//      fit(vars[i].data(),"pprimaryVertexFilter","pVtx",262548);
      }
      //fit(vars[i].data(),"","",262548);
      fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2 && pprimaryVertexFilter","hfCoinc2_pVtx",262548);
      fit(vars[i].data(),"phfCoincFilter && pprimaryVertexFilter","hfCoinc1_pVtx",262548);
      //fit(vars[i].data(),"pprimaryVertexFilter","pVtx",262548);
 
      fit(vars[i].data(),"phfCoincFilter","hfCoinc1",262548);
      fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2","hfCoinc2",262548);
      //fit(vars[i].data(),"","",262315);
 
      //fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter","hfCoinc3_pVtx",262784);
      //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2 && pprimaryVertexFilter","hfCoinc2_pVtx",262784);
      //fit(vars[i].data(),"phfCoincFilter && pprimaryVertexFilter","hfCoinc1_pVtx",262784);
      //fit(vars[i].data(),"pprimaryVertexFilter","pVtx",262784);
 
      //fit(vars[i].data(),"phfCoincFilter","hfCoinc1",262784);
      //fit(vars[i].data(),"phfPosFilter2&&phfNegFilter2","hfCoinc2",262784);
      //fit(vars[i].data(),"phfCoincFilter3","hfCoinc3",262784);
      
            
   }


}
