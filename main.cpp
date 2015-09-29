/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. “A Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity.” Nature 525, 
   no. 7568 (September 10, 2015): 261–64. doi:10.1038/nature14971.

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "classes.h"
#define __MAIN
#include "params.h"

#if (defined(GILLESPIE) && defined(FASTER_KMC)) || (defined(GILLESPIE) && defined(NORMAL)) || (defined(NORMAL) && defined(FASTER_KMC))
  #error too many methods defined!
#endif

#if !defined(GILLESPIE) && !defined(FASTER_KMC) && !defined(NORMAL) 
  #error no method defined!
#endif

#if (defined(VON_NEUMANN_NEIGHBOURHOOD) && defined(MOORE_NEIGHBOURHOOD))
  #error both VON_NEUMANN_NEIGHBOURHOOD and MOORE_NEIGHBOURHOOD defined!
#endif

#if (!defined(VON_NEUMANN_NEIGHBOURHOOD) && !defined(MOORE_NEIGHBOURHOOD))
  #error neither VON_NEUMANN_NEIGHBOURHOOD nor MOORE_NEIGHBOURHOOD defined!
#endif

extern char *NUM ; 
extern int RAND, sample, treatment, max_size ;
extern double tt ;
extern float time_to_treat ;

int sample=0 ;

void save_positions(char *name, float dz) 
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
    if (abs(int(cells[i].z+ll->r.z))<dz || cells.size()<1e4) fprintf(data,"%d %d %d %u\n",int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z),genotypes[cells[i].gen]->index) ;
  }        
  fclose(data) ; 
}

float save_2d_image(char *name, vecd li)
{
  int i,j,k;
  int minx=1<<20,maxx=-minx,miny=minx,maxy=-minx,minz=minx,maxz=-minx ;  
  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    if (cells[i].x+ll->r.x<minx) minx=int(cells[i].x+ll->r.x) ;
    if (cells[i].x+ll->r.x>maxx) maxx=int(cells[i].x+ll->r.x) ;
    if (cells[i].y+ll->r.y<miny) miny=int(cells[i].y+ll->r.y) ;
    if (cells[i].y+ll->r.y>maxy) maxy=int(cells[i].y+ll->r.y) ;
    if (cells[i].z+ll->r.z<minz) minz=int(cells[i].z+ll->r.z) ;
    if (cells[i].z+ll->r.z>maxz) maxz=int(cells[i].z+ll->r.z) ;
  }
  maxx++ ; maxy++ ; 
  float density=1.*cells.size()/(float(maxx-minx)*float(maxy-miny)*float(maxz-minz)) ;
  if (cells.size()<1e3) return density ;
  float diam=pow(float(maxx-minx)*float(maxy-miny)*float(maxz-minz),1./3) ;
  printf("density=%f\n",density) ;

  int nnn=(maxx-minx)*(maxy-miny) ;
  if (float(maxx-minx)*float(maxy-miny)>2e9) err("(maxx-minx)*(maxy-miny) too large",float(maxx-minx)*float(maxy-miny)) ;
  int *types=new int[nnn] ; // 2d array of cells' types (foremost ones only)
  short int *zbuf=new short int[nnn] ; // 2d array of cells' z positions (foremost ones only)
  BYTE *br=new BYTE[nnn] ; // 2d array of cells' brightness (foremost ones only)
  Sites **bit=new Sites*[nnn] ; // 3d array of bits: cell empty/occupied
  for (i=0;i<nnn;i++) { zbuf[i]=minz ; types[i]=-1 ; br[i]=255 ; bit[i]=new Sites(maxz-minz) ; }
  printf("%d %d\t %d %d\n",minx,maxx,miny,maxy) ;
  printf("%d x %d = %d\n",maxx-minx,maxy-miny,nnn) ;

  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int z=int(cells[i].z+ll->r.z) ;
    int adr=int(cells[i].y+ll->r.y-miny)*(maxx-minx)+int(cells[i].x+ll->r.x-minx) ;
    if (adr<0 || adr>=nnn) err("adr",adr) ;
    bit[adr]->set(z-minz) ;
    if (z>zbuf[adr]) { zbuf[adr]=z ; types[adr]=genotypes[cells[i].gen]->index ; }
  }
  
  normalize(li) ;
  float dmul=0.93 ; 
  float range=-(0.916291/(density*log(dmul))) ; 
  if (range>diam) { range=diam ; dmul=pow(0.4,1/(density*range)) ; }
  for (i=0;i<maxy-miny;i++) 
    for (j=0;j<maxx-minx;j++) {
      float d=1 ; 
      int adr=i*(maxx-minx)+j ;
      k=zbuf[adr]-minz ;
      for (float o=1;o<range;o++) {
        int il=int(i-li.y*o) ; int jl=int(j-li.x*o) ; int kl=int(k-li.z*o) ; 
        if (il>0 && jl>0 && il<maxy-miny && jl<maxx-minx && kl>0 && kl<maxz-minz) {          
          if (bit[il*(maxx-minx)+jl]->is_set(kl)) d*=dmul ;
        } else break ;
        if (d<0.4) { d=0.4 ; break ; }
      }
      br[adr]*=d ;
    }

  FILE *data=fopen(name,"w") ;    
  for (i=0;i<maxy-miny;i++) {
    for (j=0;j<maxx-minx;j++) fprintf(data,"%d %d ",types[i*(maxx-minx)+j],br[i*(maxx-minx)+j]) ;
    fprintf(data,"\n") ;
  }
  fclose(data) ; 
  delete [] types ; delete [] zbuf ; delete [] br ; delete [] bit ;
  return density ;
}
  

void save_genotypes(char *name)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d  %d  %d %d  %d\t",i, g->prev_gen,g->no_resistant,g->no_drivers, g->number) ;
      for (int j=0;j<g->sequence.size();j++) fprintf(data," %u",g->sequence[j]) ; 
      fprintf(data,"\n") ;
    } 
  }

  fclose(data) ;  
}

void save_most_abund_gens(char *name, int *most_abund)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *gg=genotypes[i] ;
    if (gg!=NULL && gg->number>0) {
      int r=0,g=0,b=0 ;
      for (int j=0;j<gg->sequence.size();j++) {
        if ((gg->sequence[j]&L_PM)==most_abund[0]) r=1 ;
        if ((gg->sequence[j]&L_PM)==most_abund[1]) g=1 ; 
        if ((gg->sequence[j]&L_PM)==most_abund[2]) b=1 ;
      }
      if (r || g || b) fprintf(data,"%d %d %d\t%d\n",r,g,b,gg->index) ;
    }
  }
  fclose(data) ;  

}

int main(int argc, char *argv[])
{
#if defined(GILLESPIE)   
  cout <<"method: GILLESPIE\n" ;
#endif
#if defined(FASTER_KMC)   
  cout <<"method: FASTER_KMC\n" ;
#endif
#if defined(NORMAL)   
  cout <<"method: NORMAL\n" ;
#endif

  int nsam ;
  if (argc!=4) { err(" Error:: arguments needed: name, no_samples, RAND. Program terminated. \n"); } 
  else { 
    NUM=argv[1] ;
    nsam=atoi(argv[2]) ;
    RAND=atoi(argv[3]) ;
  }
  cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<endl ;
  _srand48(RAND) ;
  init();
  char name[256] ;
  sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
  FILE *er=fopen(name,"w") ; fclose(er) ;
  for (sample=0;sample<nsam;sample++) { 
    reset() ;
#ifdef MAKE_TREATMENT_N
    int s=0 ; while (main_proc(max_size,-1,-1, 10)==1) { s++ ; reset() ; } ; // initial growth until max size is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; 
    treatment=1 ;       
    double max_time=2*tt ;
    main_proc(1.25*max_size,-1,max_time, 10) ; // treatment
#elif defined MAKE_TREATMENT_T
    int s=0 ; while (main_proc(-1,-1,time_to_treat, 10)==1) { s++ ; reset() ; } ; // initial growth until max time is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; 
    treatment=1 ;       
    double max_time=2*tt ;
    max_size=cells.size()*1.25 ;
    main_proc(max_size,-1,max_time, 10) ; // treatment

#else    
    int s=0 ; while (main_proc(max_size,2,-1, -1)==1) { s++ ; reset() ; } // initial growth until max size is reached
    if (s>0) printf("resetted %d times\n",s) ;
    fflush(stdout) ;
    save_data() ; 
    
    sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
    FILE *er=fopen(name,"a") ;
    fprintf(er,"%d\t%d %lf\n",sample,s,tt) ;
    fclose(er) ;

    // save some more data
    
    int *snp_no=new int[L], *snp_drivers=new int[L] ; // array of SNPs abundances
    for (int i=0;i<L;i++) { snp_no[i]=snp_drivers[i]=0 ; }
    for (int i=0;i<genotypes.size();i++) {
      if (genotypes[i]!=NULL && genotypes[i]->number>0) 
        for (int j=0;j<genotypes[i]->sequence.size();j++) {
          snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
          if (((genotypes[i]->sequence[j])&DRIVER_PM)) snp_drivers[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;
        }
    }

    save_spatial(snp_no) ;

    printf("saving PMs...\n") ;
    int most_abund[100] ;
    sprintf(name,"%s/all_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_no,max_size,0,most_abund) ;
    if (driver_adv>0 || driver_migr_adv>0) { printf("saving driver PMs...\n") ; sprintf(name,"%s/drv_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_drivers,max_size,0,NULL) ; }
    delete [] snp_no ; delete [] snp_drivers ;

    if (nsam==1) {  // do this only when making images of tumours
      printf("saving images...\n") ;
      int j=0 ;
      for (int i=0;i<genotypes.size();i++) {
        if (genotypes[i]!=NULL && genotypes[i]->number>0) genotypes[i]->index=j++ ; 
      }       

      vecd li1(1,-1,-0.3) ; // lighting direction
      float density ; 
      sprintf(name,"%s/2d_image1_%d.dat",NUM,max_size) ; density=save_2d_image(name,li1) ;
      // if you want to save more images in a single run, add new lines like this
      //       vecd li2(1,-1,-1) ; sprintf(name,"%s/2d_image2_%d.dat",NUM,max_size) ; density=save_2d_image(name,li2) ;
      sprintf(name,"%s/cells_%d.dat",NUM,max_size) ; save_positions(name,1./density) ; 
      sprintf(name,"%s/genotypes_%d.dat",NUM,max_size) ; save_genotypes(name) ;
      sprintf(name,"%s/most_abund_gens_%d.dat",NUM,max_size) ; save_most_abund_gens(name,most_abund) ;
    }
#endif
  } 
  end() ;
	return 0 ;
}
