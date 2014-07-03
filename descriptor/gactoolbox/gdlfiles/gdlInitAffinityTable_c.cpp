///////////////////////////////////////////////////////////////////////
// by Wei Zhang (wzhang009 at gmail.com), Nov., 7, 2011

#include "mex.h"
// #include <string.h>

//#define SAFEMXDESTROYARRAY(p) { if (p != NULL) { mxDestroyArray(p); p = NULL; } }
#define MYINF 1e10
double gdlComputeAffinity (double *pW, const int height, mxArray *cluster_i, mxArray *cluster_j, double *AsymAff, bool normalize);

//////////////////////////////////////////////////////////////////////////////////////////////////////
// function [affinityTab, AsymAffTab] = gdlInitAffinityTable_c (IminuszW, initClusters)
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, mxArray *prhs[]) 
{
    mxArray *graphW = prhs[0];
    mxArray *initClusters = prhs[1];
    bool *normalize;
    normalize = mxGetLogicals(prhs[2]);
    
    //mexPrintf("normalize == %d\n", *normalize);

    if (nrhs != 3) {
        mexErrMsgTxt("Wrong number of input!");
    }
    if (mxGetNumberOfDimensions(graphW) != 2 || mxGetM(graphW) != mxGetN(graphW)) {
		mexErrMsgTxt("graphW is not a square matrix!");
	}
    double *pW = mxGetPr(graphW);
    const int height = mxGetM(graphW);
 
    if (!mxIsCell(initClusters)) {
		mexErrMsgTxt("initClusters is not a cell!");
    }
    int numClusters = mxGetNumberOfElements (initClusters);
    
    // output:
	plhs[0] = mxCreateDoubleMatrix(numClusters,numClusters,mxREAL);
    mxArray *affinityTab = plhs[0]; // reference
    double *affinityTabEntry = mxGetPr(affinityTab);
    for (int i = 0; i < numClusters*numClusters; i++) { affinityTabEntry[i] = - MYINF; }
	plhs[1] = mxCreateDoubleMatrix(numClusters,numClusters,mxREAL);
    mxArray *AsymAffTab = plhs[1]; // reference
    double *AsymAffTabEntry = mxGetPr(AsymAffTab);
    for (int i = 0; i < numClusters*numClusters; i++) { AsymAffTabEntry[i] = - MYINF; }
    
    // computing
    double tmpAsymAff[2];
    double *pTable_j = affinityTabEntry;
    for (int j = 0; j < numClusters; j++) {
        mxArray *cluster_j = mxGetCell (initClusters, j);  // cluster j
        for (int i = 0; i < j; i++) {
            pTable_j[i] = gdlComputeAffinity(pW, height, mxGetCell (initClusters, i), cluster_j, tmpAsymAff, *normalize);  // affinityTabEntry[i+j*numClusters]
            AsymAffTabEntry[i+j*numClusters] = tmpAsymAff[0];
            AsymAffTabEntry[j+i*numClusters] = tmpAsymAff[1];
        }
        pTable_j += numClusters;
    }
   
}