/* 	
 *   preconditioner.c Solve given the hierarchy of graphs
 *
 */

/*
 CMG, Copyright (c) 2008-2020  Ioannis Koutis, Gary Miller               
 
 The CMG solver is distributed under the terms of the GNU General Public  
 Lincense Version 3.0 of the Free Software Foundation.                    
 CMG is also available under other licenses; contact authors for details.
 
 /*



/* 
 * preconditioner
 * 
 * Input:  pointer to array of hierarchy levels: H
 *      :  --precision-- vector: b
 *      :  integer: l
 *
 * Output:  --precision-- vector: x 
 *
 *
 */

#include "cmg.h"
#include "mex.h"


void preconditioner(s_hlevel *H, precision *b, int level, int iter, precision *x)
{
    mSize n,m,nc;
    int i,j;
    precision *y, *Bb, *r, *z;
    precision *b_small;
    precision s;
    
            
    n = H[level].A.n;


    /* last level direct */
    if ( H[level].islast && !(H[level].iterative)) {
        ldl_solve(&(H[level].chol) , b, x); 
        x[n]= (precision) 0.0;

        return;
    }


    /* last level iterative */
    if ( H[level].islast  && H[level].iterative) {
        vvmul(H[level].invD,b,x,n);                                                 /*  x = H{level}.invD.*b */
        return; 
    }

    /* main cycle */
    nc = H[level].nc;
    
    /* temp storage initialization */
    y  = H[level].lws2;                       
    r  = H[level].lws2;
    
    /* initialize the solution vector x */
    for (i=0;i<n;i++)
        x[i] = (precision) 0.0;
    
    Bb = H[level].lws1;                           /* fixed working space */
    vvmul(H[level].invD,b,Bb,n);                  /* Bb=H{level}.invD.*b */
 
          
    for (j=1;j<=iter; j++){
        
                
        /* Jacobi pre-smooth */
        if (j==1){
            for (i=0;i<n;i++)
            x[i] = Bb[i];}
        else{
            sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, x, y);
            vmv(b,y,r,n);                                                                   
            vvmul(H[level].invD,r,y,n);               /*  y = invD.*(b-A x)   */
            vpv(x,y,x,n);                             /*  x = x + invD.*(b-Ax) */        
        }
     
        b_small = H[level].sws1;  
        z = H[level].sws2;
        
        
        sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, x, y);
        vmv(b,y,r,n);                           /* r= b - (H{level}.A*x) */
        
   
        rmvecmul(H[level].cI , r, n, b_small, nc);                       /*  b_small = Rt*r;  */
        
        preconditioner(H,b_small,level+1, H[level].repeat,z);            /*  z = preconditioner(r) */
       
        trmvecmul(H[level].cI, z, nc, y, n);                             /*  y = R*z */
        
        vpv(y,x,x,n);  /*  x = y + R*z */
        
    
        /* Jacobi post-smooth */
        sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, x, y);
        vmv(b,y,r,n);                                                                   
        vvmul(H[level].invD,r,y,n);               /*  y = invD.*(b-A x)   */
        vpv(x,y,x,n);                             /*  x = x + invD.*(b-Ax) */        
    }
    
 
    return;
 

}