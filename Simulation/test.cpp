// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
// [[Rcpp::export]]
SEXP MYFUNC(NumericMatrix A,const int MYC,const int NOALL,NumericMatrix PAS) {
	arma::mat table1(PAS.nrow()*2, A.ncol());
	arma::mat table2(A.begin(), A.nrow(), A.ncol(), false);
	arma::mat PA(PAS.begin(), PAS.nrow(), PAS.ncol(), false);
	int a, MY, CORUN, MYCRUN, nc,i,l,k,MYPRUN;
	float x,y,NORA;
	std::vector <int> TMPXXX;
	std::vector <int> XXX;
	const float REALL=exp(-1);
	int XXX1 = 0;
	for(MYPRUN=0;MYPRUN<PAS.nrow();MYPRUN++) {
		for (k=0;k<2;k++) {
			if (k==0) MY=PA(MYPRUN,0);
			else      MY=PA(MYPRUN,1);
			CORUN=0;
			for (MYCRUN=0;MYCRUN<MYC;MYCRUN++) {
				nc=0;
				x=REALL;
				y=x;
				NORA=((float) rand() / (RAND_MAX));
				// generate the num of recombination using poission
				while (x<=NORA) {
					nc++;
					y=y/nc;
					x=x+y;  
				}
				// randomly create the unique location of the occurance of recombination
				for (i=0; i<nc; i++) TMPXXX.push_back ((rand() % (NOALL-2))+1);
				if (nc==1) XXX.push_back (TMPXXX [0]);
				if (nc>1) {
					loop_1:
					sort (TMPXXX.begin(),TMPXXX.end());
					unique_copy (TMPXXX.begin(),TMPXXX.end(),back_inserter(XXX));
					if (XXX.size()!=nc) {
						TMPXXX.push_back ((rand() % (NOALL-2))+1);
						XXX.clear();
						goto loop_1;
					}
				}
				if (nc==0) {
					NORA=((float) rand() / (RAND_MAX));
					if (NORA<=0.5) { 
						for(;CORUN<(MYCRUN+1)*NOALL;CORUN++) table1(XXX1,CORUN)= table2(MY*2,CORUN);
					} else {
						for(;CORUN<(MYCRUN+1)*NOALL;CORUN++) table1(XXX1,CORUN)= table2(MY*2+1,CORUN);
					}
				} else {
					NORA=((float) rand() / (RAND_MAX));
					if (NORA<=0.5) {
						for (l=0;l<nc;l++){
							a= l%2;
							switch (a) {
							case 0: for (;CORUN<MYCRUN*NOALL+XXX[l];CORUN++) table1(XXX1,CORUN)=table2(MY*2,CORUN);
								break;
							case 1: for (;CORUN<MYCRUN*NOALL+XXX[l];CORUN++) table1(XXX1,CORUN)=table2(MY*2+1,CORUN);
								break;}
						} switch (a) {
							case 0: for (;CORUN<(MYCRUN+1)*NOALL;CORUN++) table1(XXX1,CORUN)=table2(MY*2+1,CORUN);
								break;
							case 1: for (;CORUN<(MYCRUN+1)*NOALL;CORUN++) table1(XXX1,CORUN)=table2(MY*2,CORUN);
								break;
						}
					} else {
						for (l=0;l<nc;l++) {
							a= l%2;
							switch (a) {
								case 0: for (;CORUN<MYCRUN*NOALL+XXX[l];CORUN++) table1(XXX1,CORUN)=table2(MY*2+1,CORUN);
								break;
								case 1: for (;CORUN<MYCRUN*NOALL+XXX[l];CORUN++) table1(XXX1,CORUN)=table2(MY*2,CORUN);
								break;
							}
						} switch (a) {
							case 0: for (;CORUN<(MYCRUN+1)*NOALL;CORUN++) table1(XXX1,CORUN)=table2(MY*2,CORUN);
							break;
							case 1: for (;CORUN<(MYCRUN+1)*NOALL;CORUN++) table1(XXX1,CORUN)=table2(MY*2+1,CORUN);
							break;
						}
					}
				}
				XXX.clear();
				TMPXXX.clear();
			}
			XXX1++;
		}
	}
    return Rcpp::wrap(table1);
}

