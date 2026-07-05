/* ========================================================================== */
/* === CMG/Source/Hierarchy/forest_components.c ============================= */
/* ========================================================================== */


/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/*
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */

/* Review 2026-07-03 (see ../../../CMG-cpp-dev/dev-comments/version-comparison.md)
 *
 * Only ONE genuine defect was found and fixed in this file:
 *   - The `sizes` workspace is allocated with mxMalloc (NOT zero-initialized)
 *     and the original init loop started at j=1, leaving sizes[0] (= the size
 *     of component label 1) uninitialized.  It is now initialized from j=0.
 *
 * Two earlier "fixes" were reverted after closer analysis:
 *   - cI init from j=1 is NOT a bug: mxCreateNumericMatrix zero-initializes
 *     its data, so cI[0] is already 0.  (Kept at j=0 as harmless belt-and-
 *     suspenders.)
 *   - The `A[j]!=(n+1)` guard was RESTORED.  It skips nodes whose parent
 *     pointer is the out-of-range sentinel n+1; removing it risks an
 *     out-of-bounds read of cI[n+1] on such inputs.
 *
 * No empty-last-component adjustment is needed: ccI is incremented exactly
 * once per completed component, so the final `ccI = ccI-1` already yields the
 * correct component count.  (The MATLAB forest_components_ uses `nc = ccI` and
 * needs its own correction; that difference does not apply here.)
 */

#include "mex.h"
#include <stdlib.h>

#define A_IN prhs[0]
#define cI_OUT plhs[0]
#define nc_OUT plhs[1]
#define sizes_OUT plhs[2]


void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{


    mwIndex i,j, jwalk, ccI, *buffer, *bfr, bufferI;
    mwSize n,m,gsize,t;
    uint32_T *A, *cI, *sizes,*sizesp;;


    if (nrhs >1)
        mexErrMsgTxt("Wrong number of input arguments");
    if ((nlhs > 3))
        mexErrMsgTxt("Wrong number of output arguments");

    n = mxGetN(A_IN);
    m = mxGetM(A_IN);
    if ((m!=1) && (n!=1))
        mexErrMsgTxt("Input must be a vector");
    if (mxGetClassID(A_IN) != mxUINT32_CLASS)
        mexErrMsgTxt("Input vector must be of type uint32");
    A = (uint32_T *) mxGetPr(A_IN);
    if (m>n)
        n=m;


    /* Create output (mxCreateNumericMatrix zero-initializes cI) */
    cI_OUT = mxCreateNumericMatrix(n, 1, mxUINT32_CLASS, mxREAL);
    if (cI_OUT == NULL)
         mexErrMsgTxt("Unsufficient memory");
    cI = (uint32_T *) mxGetPr(cI_OUT);
    for (j=0 ;j<n ;j++)
        cI[j]=0;

    t = 1;
    gsize = 100;
    buffer = mxMalloc(gsize*sizeof(mwIndex));
    if (nlhs == 3){
      sizes = mxMalloc(n*sizeof(double));
      for (j=0 ;j<n ;j++)      /* FIX: was j=1; sizes[0] was left uninitialized */
          sizes[j]=0.0;
    }

    if (nlhs <3 ){  /* sizes are not recorded */
    ccI = 1;
    for (j=0 ; j<n ; j++){
        bufferI=0;
        jwalk = j;
        /* Tree walk */
        while (cI[jwalk]==0 && A[j]!=(n+1)){
            if (bufferI==t*gsize) {
                t = t+1;
                bfr = mxRealloc(buffer,t*gsize*sizeof(mwIndex));
                buffer = bfr;
            }
            cI[jwalk] = ccI;
            buffer[bufferI] = jwalk;
            bufferI = bufferI+1;
            jwalk = A[jwalk];}
        if (cI[jwalk] != ccI)
            for (i=0; i<bufferI; i++)
                cI[buffer[i]] = cI[jwalk];
        else
            ccI = ccI+1;
    }}

    if (nlhs == 3 ){  /* sizes are recorded */
    ccI = 1;
    for (j=0 ; j<n ; j++){
        bufferI=0;
        jwalk = j;
        /* Tree walk */
        while (cI[jwalk]==0 && A[j]!=(n+1)){
            if (bufferI==t*gsize) {
                t = t+1;
                bfr = mxRealloc(buffer,t*gsize*sizeof(mwIndex));
                buffer = bfr;
            }
            cI[jwalk] = ccI;
            buffer[bufferI] = jwalk;
            bufferI = bufferI+1;
            jwalk = A[jwalk];}
        if (cI[jwalk] != ccI){
            for (i=0; i<bufferI; i++)
                cI[buffer[i]] = cI[jwalk];}
        else
            ccI = ccI+1;
        sizes[(mwIndex) cI[jwalk]-1] = sizes[(mwIndex) cI[jwalk]-1] + (double) bufferI;
    }}

    ccI = ccI-1;
    nc_OUT = mxCreateDoubleScalar(ccI);

    if (nlhs == 3) {
        sizes_OUT = mxCreateDoubleMatrix(1,ccI,mxREAL);
        sizesp = (uint32_T *) mxGetPr(sizes_OUT);
        for (j=0 ; j<ccI; j++)
            sizesp[j] = sizes[j];
        mxFree(sizes);
    }

    mxFree(buffer);

}
