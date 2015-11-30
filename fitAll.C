#include "fit.C"

void fitAll(){

   static const int Nvar = 7;
   string vars[Nvar] = {"hiHF","hiHFhit","hiNtracks","hiNpix","hiEB","hiEE","hiET"};


   for(int i = 0; i < Nvar; ++i){

      if(0){
      fit(vars[i].data(),"abs(zVtx[0])<10","innerVtx",262548);
      fit(vars[i].data(),"abs(zVtx[0])>10","outerVtx",262548);
      
      fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])<10","hfCoinc3_innerVtx",262548);
      fit(vars[i].data(),"phfCoincFilter3 && abs(zVtx[0])>10","hfCoinc3_outerVtx",262548);

      fit(vars[i].data(),"","",262548);
      fit(vars[i].data(),"phfCoincFilter3 && pprimaryVertexFilter","hfCoinc3_pVtx",262548);
      }

      fit(vars[i].data(),"phfCoincFilter3","",262548);
      fit(vars[i].data(),"","",262315);
      fit(vars[i].data(),"phfCoincFilter3","",262315);
            
   }


}
