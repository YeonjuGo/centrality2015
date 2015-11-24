#include "fit.C"

void fitAll(){

   fit("hiHF","");
   fit("hiHF","phfCoincFilter");
   fit("hiHF","phfCoincFilter2");
   fit("hiHF","phfCoincFilter3");
   //   fit("hiHF","phfCoincFilter4");
   //   fit("hiHF","phfCoincFilter5");

   fit("hiEB","");
   fit("hiEB","phfCoincFilter");
   fit("hiEB","phfCoincFilter2");
   fit("hiEB","phfCoincFilter3");
   //   fit("hiEB","phfCoincFilter4");
   //   fit("hiEB","phfCoincFilter5");

   fit("hiET","");
   fit("hiET","phfCoincFilter");
   fit("hiET","phfCoincFilter2");
   fit("hiET","phfCoincFilter3");
   //   fit("hiET","phfCoincFilter4");
   //   fit("hiET","phfCoincFilter5");

}
