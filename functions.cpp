/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. "Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity" Nature 525, 
   no. 7568 (September 10, 2015): 261-64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file 
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/


#include <math.h>
#include "params.h"
#include "classes.h"

void save_snps(char *name,int *n, int total, int mode, int *most_abund) 
{
  const float cutoff=0.01 ;
  FILE *f=fopen(name,"w") ;
  if (f==NULL) err(name) ;
  int i, j, nsnps=0, nsnpsc=0;
  for (i=0;i<L;i++) { 
    if (n[i]>(1e-4)*total) nsnps++ ;
    if (1.*n[i]/total>cutoff) nsnpsc++ ;
  }
  float *abund=new float[nsnps], tempd ;
  int *num=new int[nsnps], temp ;
  nsnps=0 ;
  for (i=0;i<L;i++) if (n[i]>(1e-4)*total) { num[nsnps]=i ; abund[nsnps]=float(1.*n[i]/total)*(1+0.000001*i/L) ; nsnps++ ; }
  quicksort2(abund,num,0,nsnps-1) ;

  if (mode) {
    for (i=0;i<nsnps;i++) fprintf(f,"%d %d %f\n",i,num[i],abund[i]) ;
  } else {
    for (i=0;i<nsnps;i++) if (abund[i]>cutoff || i<100) fprintf(f,"%d %d %f\n",i,num[i],abund[i]) ;
  }
  if (most_abund!=NULL) { // store first 100 most abundant PMs
    for (i=0;i<MIN(100,nsnps);i++) most_abund[i]=num[i] ; 
  }
  delete [] abund ; delete [] num ;
  fclose(f) ;
}

void save_snp_corr(char *name, Hist *snps)
{
  FILE *f=fopen(name,"w") ;
  int last ;
  for (last=_bins-1;snps[last].x==0 && last>0;last--) ;
  for (int i=0;i<=last;i++) if (snps[i].n>0) fprintf(f,"%d %f %f\n",i*_resol,
      float(1.*snps[i].x/snps[i].n),float(sqrt(( (1.*snps[i].x2/snps[i].n) - (1.*snps[i].x/snps[i].n)*(1.*snps[i].x/snps[i].n))/(snps[i].n-1)))) ; 
    else fprintf(f,"%d 0 0\n",i*_resol) ;
  fclose(f) ;    
}



int how_many_SNPs_identical(Genotype *a, Genotype *b)
{
  int i,j,n=0;
  int al=a->sequence.size(), bl=b->sequence.size() ;
  for (i=0;i<al;i++)
    for (j=0;j<bl;j++) 
      if (a->sequence[i]==b->sequence[j]) n++ ;
  return n ;
}

int how_many_SNPs_identical(Genotype *a, Genotype *b, float cutoff, int *snp_no)
{
  int i,j,n=0, ntot=cells.size();
  int al=a->sequence.size(), bl=b->sequence.size() ;
  for (i=0;i<al;i++)
    for (j=0;j<bl;j++) 
      if (a->sequence[i]==b->sequence[j] && snp_no[(a->sequence[i])&L_PM]>cutoff*ntot) n++ ;
  return n ;
}


void select_two_random_cells(int &i, int &j, int &n)
{
  int ntot=cells.size() ;
  double dist ;
  do {
    i=int(_drand48()*ntot) ; j=int(_drand48()*ntot) ;
    vecd ri=lesions[cells[i].lesion]->r, rj=lesions[cells[j].lesion]->r, rij=ri-rj ;
    dist=SQR(rij.x+cells[i].x-cells[j].x)+SQR(rij.y+cells[i].y-cells[j].y)+SQR(rij.z+cells[i].z-cells[j].z) ;
  } while (_drand48()<0.5 && dist>_drand48()*100) ; // choose preferentially cells which are close to one another
  dist=sqrt(dist) ;
  
  n=int(dist/_resol) ;
  if (n<0 || n>=_bins) err("n",n);
}

float average_distance_ij()
{
  int ntot=cells.size() ;
  double avdist=0 ;
  for (int n=0;n<ntot;n++) {
    int i=int(_drand48()*ntot), j=int(_drand48()*ntot) ; 
    vecd ri=lesions[cells[i].lesion]->r, rj=lesions[cells[j].lesion]->r, rij=ri-rj ;
    avdist+=sqrt(SQR(rij.x+cells[i].x-cells[j].x)+SQR(rij.y+cells[i].y-cells[j].y)+SQR(rij.z+cells[i].z-cells[j].z)) ;
  } 
  return (avdist/ntot) ;
}

void snps_corr(Hist *snps)
{
  int i,j,k,n;
  int ntot=cells.size() ;
    
  for (i=0;i<_bins;i++) snps[i].x=snps[i].x2=snps[i].n=0 ;
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    snps[n]+=how_many_SNPs_identical(genotypes[cells[i].gen],genotypes[cells[j].gen]) ;
  }  
}

void snps_corr_cutoff(Hist *snps,float cutoff,int *snp_no) 
{
  int i,j,k,n;
  int ntot=cells.size() ;
    
  for (i=0;i<_bins;i++) snps[i].x=snps[i].x2=snps[i].n=0 ;
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    snps[n]+=how_many_SNPs_identical(genotypes[cells[i].gen],genotypes[cells[j].gen],cutoff,snp_no) ;
  }  
  
}


void snps_corr_cond_driver(Hist *snps)
{
  int ntot=cells.size() ;
  int i,j,k,n,d,ai,aj;
    
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    Genotype *gi=genotypes[cells[i].gen], *gj=genotypes[cells[j].gen]  ;    
    d=0 ;
    for (ai=0;ai<gi->sequence.size();ai++) if (gi->sequence[ai]<0) {
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
    }
    if (d>0) {
      snps[n]+=how_many_SNPs_identical(genotypes[cells[i].gen],genotypes[cells[j].gen]) ;
    }
  }  
}


void find_p_driver(Hist *pr1, Hist *pr2, Hist *pr3) 
{
//Probability of finding the same driver in two cells separated by some distance x.
// 1) prob. that two cells at distance x have at least one common driver
// 2) prob. that two cells at distance x have the same last driver
// 3) prob. that two cells at distance x have the same first driver
  int i,j,k,n,ai,aj,d,di,dj;
  int ntot=cells.size() ;
    
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    Genotype *gi=genotypes[cells[i].gen], *gj=genotypes[cells[j].gen]  ;    
    
    d=0 ; 
    for (ai=0;ai<gi->sequence.size();ai++) if (gi->sequence[ai]<0) {
      di++ ;
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
    }
    pr1[n]+=(d>0?1:0) ;

    d=0 ;
    for (ai=gi->sequence.size()-1;ai>=0;ai--) if (gi->sequence[ai]<0) {
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
      break ;
    }
    pr2[n]+=(d>0?1:0) ;

    d=0 ;
    for (ai=0;ai<gi->sequence.size();ai++) if (gi->sequence[ai]<0) {
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
      break ;
    }
    pr3[n]+=(d>0?1:0) ;
        
  }  
}
