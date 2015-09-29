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

// Previous versions:

// v1.2.2 - optional pushing added
// v1.2.1 - saves most_abundant
// v1.2.0 - new modes added, files reorganised
// v1.1.0 - new modes added
// v1.0.1 - drivers not saved anymore
// v1.0

// Compilation:

// under windows, compile with g++ simulation.cpp main.cpp functions.cpp -w -O3 -lpsapi -o cancer.exe
// under linux, compile with g++ simulation.cpp main.cpp functions.cpp -w -O3 -o cancer.exe



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp
#include <vector>
#include <iostream>
using namespace std;

char *NUM ; // name given as 1st argument from the command line

#include "params.h"
#include "classes.h"

#if defined __linux
typedef unsigned int DWORD ;
int memory_taken() // return memory available in MB
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available memory in MB
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#elif defined __APPLE__
typedef unsigned int DWORD ;
int memory_taken()
{
  return 0 ; // not implemented
}
#else
#include <windows.h>
#include <psapi.h>
int memory_taken() // return memory available in MB
{
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (int) (info.WorkingSetSize/(1<<20));
}
#endif

void err(char *reason)
{
  cout <<reason<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif  
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, char *a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

static long long unsigned int _x=0x000100010001LL, _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  // works only on compilers with long long int!
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}

void _srand48(int a) { _x=a ; }

void init();
void end() ;

double tt=0, tt_at_start ;
int start_clock ;

int L=0 ; // total number of SNPs
int volume ; // total volume of the tumor
vector <int> drivers ; // vector of driver mutations
FILE *drivers_file ;
int treatment=0, cells_at_start ;
FILE *times ; 
extern int sample ;
int RAND ; // random seed
char *timesbuffer ;

int poisson(void)  // generates k from P(k)=exp(-gamma) gamma^k / k!
{
  const double l=exp(-gama) ;
  double p=1. ;
  int k=0 ;
  do {
    k++ ;
    p*=_drand48() ;
  } while (p > l) ;
  return k - 1 ;
}


vector<Cell> cells ;



int Lesion::nl=0 ;
double Lesion::maxdisp=0 ;
double max_growth_rate ;

Genotype::Genotype(void) 
{ 
  death[0]=death0 ; death[1]=death1 ; growth[0]=growth0 ; growth[1]=growth1 ;
#ifdef MIGRATION_MATRIX
  m[0]=migr[0][0] ; m[1]=migr[1][0] ;
#else
  m[0]=m[1]=migr ; 
#endif
  number=1 ; no_resistant=no_drivers=0 ; sequence.clear() ; prev_gen=-1 ;
}

Genotype::Genotype(Genotype *mother, int prevg, int no_snp) { 
  death[0]=mother->death[0] ; growth[0]=mother->growth[0] ; m[0]=mother->m[0] ;
  death[1]=mother->death[1] ; growth[1]=mother->growth[1] ; m[1]=mother->m[1] ;
  prev_gen=prevg ;
  sequence=mother->sequence ; no_resistant=mother->no_resistant ; no_drivers=mother->no_drivers; 
  for (int i=0;i<no_snp;i++) {
    if ((driver_adv>0 || driver_migr_adv>0) && _drand48()<driver_prob/gama) { 
      float q=_drand48() ;
      if (driver_mode<2 || q<0.5) {
        death[0]*=1-driver_adv*driver_balance ; 
        growth[0]*=1+driver_adv*(1-driver_balance) ; if (max_growth_rate<growth[0]) max_growth_rate=growth[0] ;
      }
      if (driver_migr_adv>0 && ((q>=0.5 && driver_mode==2) || driver_mode==1)) {
        m[0]*=1+driver_migr_adv ; if (m[0]>max_migr) m[0]=max_migr ;
        m[1]*=1+driver_migr_adv ; if (m[1]>max_migr) m[1]=max_migr ;
      }
      // drivers decrease prob. of death or increase prob. of growth
      drivers.push_back(L) ; //fprintf(drivers_file,"%d ",L) ; fflush(drivers_file) ; 
      sequence.push_back((L++)|DRIVER_PM) ; no_drivers++ ;
    } else {
      if (_drand48()<gama_res/gama) {  
        sequence.push_back((L++)|RESISTANT_PM) ; no_resistant++ ; // resistant mutation
        death[1]=death0 ; growth[1]=growth0 ; 
#ifdef MIGRATION_MATRIX
        m[0]=migr[0][1] ; m[1]=migr[1][1] ;
#endif
      } 
      else sequence.push_back(L++) ;
    }
  }
  if (L>1e9) err("L too big") ;
  number=1 ;
}

vector<Genotype*> genotypes ;
vector<Lesion*> lesions ;

void Lesion::update_wx()
{
  int i,j,k;
  int nwx=int(wx*1.25) ;
  if (nwx%2==1) nwx++ ; // make sure it's even
  int dwx=(nwx-wx)/2 ;

#ifdef PUSHING
  for (i=0;i<wx;i++) {
    for (j=0;j<wx;j++) {
      Sites *np=new Sites[nwx] ;
      for (k=0;k<nwx;k++) np[k]=-1 ;
      for (k=dwx;k<wx+dwx;k++) np[k]=p[i][j][k-dwx] ;
      delete p[i][j] ;
      p[i][j]=np ;
    }
  }

  Sites ***np=new Sites**[nwx] ;
  for (i=0;i<nwx;i++) np[i]=new Sites*[nwx] ;
    
  for (i=0;i<nwx;i++) {
    for (j=0;j<nwx;j++) {
      if (i<dwx || i>=wx+dwx || j<dwx || j>=wx+dwx) {
        np[i][j]=new Sites[nwx] ;
        for (k=0;k<nwx;k++) np[i][j][k]=-1 ;
      } else {
        np[i][j]=p[i-dwx][j-dwx] ; 
      }
    }
  }

  for (i=0;i<wx;i++) delete [] p[i] ; 

#else
  for (i=0;i<wx*wx;i++) {
    Sites *np=new Sites(nwx) ;
    for (k=dwx;k<wx+dwx;k++) if (p[i]->is_set(k-dwx)) np->set(k) ;
    delete p[i] ;
    p[i]=np ;
  }

  if (float(nwx)*float(nwx)>2e9) err("nwx too large",nwx) ;
  Sites **np=new Sites*[nwx*nwx] ;
    
  for (i=0;i<nwx;i++) {
    for (j=0;j<nwx;j++) {
      if (i<dwx || i>=wx+dwx || j<dwx || j>=wx+dwx) {
        np[i*nwx+j]=new Sites(nwx) ;
      } else {
        np[i*nwx+j]=p[(i-dwx)*wx+j-dwx] ; 
      }
    }
  }
#endif

  delete [] p ;
  p=np ;
  wx=nwx ;

}

void Lesion::one_move_step() {
  int i,j;
  double mthis=this->n ;
  for (i=0;i<closest.size();i++) {
    vecd dr=lesions[closest[i]]->r - this->r ;
    double r2=squared(dr), sumrad2=SQR(this->rad+lesions[closest[i]]->rad) ;
    if (r2<sumrad2) {
      double mi=lesions[closest[i]]->n ;
      double disp=(sqrt(sumrad2/r2)-1) ;
      if (fabs(disp)>maxdisp) maxdisp=fabs(disp) ;
      dr*=disp*1.1 ;
      this->r-=dr*mi/(mi+mthis) ;
      lesions[closest[i]]->r+=dr*mthis/(mi+mthis) ;
    }    
  }
}

void Lesion::find_closest() 
{
  rold=r ;
  closest.clear() ;
  for (int i=0;i<lesions.size();i++) {
    vecd dr=this->r - lesions[i]->r ;
    double r2=squared(dr) ;
    if (r2>0 && r2<2*(SQR(this->rad+lesions[i]->rad))) {
      closest.push_back(i) ;
    }
  }  
}

void Lesion::reduce_overlap()
{
  int i,j,k,temp ;
  int *ind=new int[lesions.size()] ;
  for (j=0;j<lesions.size();j++) ind[j]=j ;
  do {
    maxdisp=0 ;
    for (j=0;j<lesions.size();j++) { k=_drand48()*lesions.size() ; SWAP(ind[j],ind[k]) ; }
    for (j=0;j<lesions.size();j++) {  // go through a random permutation
      i=ind[j] ; 
      lesions[i]->one_move_step() ; 
        
      vecd dr=lesions[i]->r - lesions[i]->rold ; 
      if (squared(dr)>SQR(lesions[i]->rad)) lesions[i]->find_closest() ; 
    }    
  } while (maxdisp>1e-2) ;  
  delete [] ind ;
}  

void reset() 
{
  tt=0 ; L=0 ; max_growth_rate=growth0 ;
  treatment=0 ; 
  for (int i=0;i<genotypes.size();i++) if (genotypes[i]!=NULL) delete genotypes[i] ;
  genotypes.clear() ; genotypes.push_back(new Genotype) ;  
  for (int i=0;i<lesions.size();i++) delete lesions[i] ;
  lesions.clear() ;
  cells.clear() ; volume=0 ;
  drivers.clear() ;
  lesions.push_back(new Lesion(0,0, 0,0,0)) ;
  
  // erase output buffer for "times"
#if defined __linux
  times->_IO_write_ptr = times->_IO_write_base ;
#elif defined __APPLE__
  // not defined yet
#else
  times->_ptr = times->_base ; // this operates on elements of _iobuf and is specific to Windows GNU C++
#endif
}

#ifdef MOORE_NEIGHBOURHOOD
const int _nonn=26 ;
const int kx[27]={0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1},
          ky[27]={0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1},
          kz[27]={0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1};
int kln[27] ; // this is filled with lengths of (kx,ky,kz)
#endif

#ifdef VON_NEUMANN_NEIGHBOURHOOD
const int _nonn=6 ;
const int kx[7]={0,1,-1,0,0,0,0},
          ky[7]={0,0,0,1,-1,0,0},
          kz[7]={0,0,0,0,0,1,-1};
int kln[7] ; // this is filled with lengths of (kx,ky,kz)
#endif



void init()
{
  int i,j,k;
  for (i=0;i<=_nonn;i++) kln[i]=sqrt(1.*SQR(kx[i])+1.*SQR(ky[i])+1.*SQR(kz[i])) ;

  char txt[256] ;
  sprintf(txt,"mkdir %s",NUM) ; system(txt) ;
  sprintf(txt,"%s/%s_%d.dat",NUM,NUM,RAND) ; times=fopen(txt,"w") ;
  timesbuffer=new char[(1<<16)] ;
  setvbuf (times , timesbuffer , _IOFBF , (1<<16));  // this is to prevent saving data if no fflush is attempted 
                                                  // (this e.g. allows one to discard N<256)
  start_clock=clock() ;
}

void end() {
  fclose(times) ; 
}

#ifdef PUSHING
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=_nonn ;
  for (int n=1;n<=_nonn;n++) nfree-=p[(wx+z+kz[n])%wx][(wx+y+ky[n])%wx][(wx+x+kx[n])%wx]==-1?0:1 ;
  return nfree ;
}
#else
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=_nonn ;
  for (int n=1;n<=_nonn;n++)
    nfree-=p[((wx+z+kz[n])%wx)*wx + (wx+y+ky[n])%wx]->is_set((wx+x+kx[n])%wx) ;
  return nfree ;
}
inline void Lesion::choose_nn(int &x, int &y, int &z)
{
  static int nns[_nonn] ;
  int no=0,n ;
  for (n=1;n<=_nonn;n++)
    if (p[((wx+z+kz[n])%wx)*wx + (wx+y+ky[n])%wx]->is_set((wx+x+kx[n])%wx)==0) nns[no++]=n ;
  if (no==0) { x=-1000000 ; return ; }
  n=nns[int(_drand48()*no)] ; 
  z=(wx+z+kz[n])%wx ; y=(wx+y+ky[n])%wx ; x=(wx+x+kx[n])%wx ;
}
#endif


inline int free_sites(int n)
{
  Lesion *ll=lesions[cells[n].lesion] ;
  int wx=ll->wx ; 
  int k=cells[n].x+wx/2, j=cells[n].y+wx/2, i=cells[n].z+wx/2 ; 
  return ll->no_free_sites(k,j,i) ;  
}


void quicksort2(float *n, int *nums, int lower, int upper)
{
	int i, m, temp ;
  float pivot, tempd;
	
	if (lower < upper)
	{
		SWAPD(n[lower], n[(upper + lower) / 2]); SWAP(nums[lower], nums[(upper + lower) / 2]);
		pivot = n[lower];
		m = lower;
		for (i = lower + 1; i <= upper; i++)
			if (n[i] > pivot)
			{
				m++;
				SWAPD(n[m], n[i]); SWAP(nums[m], nums[i]);
			}
		SWAPD(n[lower], n[m]); SWAP(nums[lower], nums[m]);
		quicksort2(n, nums, lower, m - 1);
		quicksort2(n, nums, m + 1, upper);
	}
}


void save_data()
{
  int i,j,ntot=cells.size(), nsurf=0 ;
  double raver=0, raver2=0 ;
  int no_on_surface=0 ;
  int no_resistant=0, no_resistant_surf=0 ;
  int cells_drv=0,cells_drv_surf=0 ;
  double drv_per_cell=0,drv_per_cell_surf=0, pms_per_cell=0 ;
  double aver_growth_rate=0,av_migr=0 ;

  int *snp_no=new int[L] ; // array of SNPs abundances
  for (i=0;i<L;i++) { snp_no[i]=0 ; }
  for (i=0;i<genotypes.size();i++) {
    if (genotypes[i]!=NULL && genotypes[i]->number>0) {
      for (int j=0;j<genotypes[i]->sequence.size();j++) snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
    }
  }  
  int snps_det=0 ;
  for (i=0;i<L;i++) if (snp_no[i]>cutoff*ntot) snps_det++ ;
  delete [] snp_no ;

  for (i=0;i<ntot;i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int wx=ll->wx ; 
    double rr=SQR(cells[i].x+ll->r.x)+SQR(cells[i].y+ll->r.x)+SQR(cells[i].z+ll->r.x) ;
    raver+=sqrt(rr) ; raver2+=rr ;
#ifdef PUSHING
    if (ll->p[cells[i].z+wx/2][cells[i].y+wx/2][cells[i].x+wx/2]==-1) err("p[][][]=",i) ;
#else
    if (ll->p[(cells[i].z+wx/2)*wx+cells[i].y+wx/2]->is_set(cells[i].x+wx/2)==0) err("save: is_set: ",i) ;
#endif

    Genotype *g=genotypes[cells[i].gen] ; if (g==NULL) err("g=NULL)") ;
    int free_sites=ll->no_free_sites(cells[i].x+wx/2,cells[i].y+wx/2,cells[i].z+wx/2) ;
    int is_on_surface=(free_sites>0?1:0) ;    
    pms_per_cell+=g->sequence.size() ;

    if (g->no_resistant) {
      no_resistant++ ; 
      if (is_on_surface) no_resistant_surf++ ;
    }
    if (g->no_drivers>0) {
      cells_drv++ ; drv_per_cell+=g->no_drivers ; 
      if (is_on_surface) { cells_drv_surf++ ; drv_per_cell_surf+=g->no_drivers ; }
    }
    if (is_on_surface) nsurf++ ;    
    aver_growth_rate+=g->growth[treatment]*free_sites/float(_nonn) ;
    av_migr+=g->m[treatment] ;
  }
  raver/=ntot ; raver2/=ntot ; aver_growth_rate/=timescale ;
  drv_per_cell/=ntot ; drv_per_cell_surf/=nsurf ; pms_per_cell/=ntot ;
  av_migr/=ntot ;

  // 1.ntot 2.time  3.#genotypes  4.radius      
#ifndef CORE_IS_DEAD
  fprintf(times,"%d %lf %d %lf  ",ntot,tt,genotypes.size(),raver) ; //sqrt(raver2-raver*raver)) ;
#else
  fprintf(times,"%d %lf %d %lf  ",volume,tt,genotypes.size(),raver) ; //sqrt(raver2-raver*raver)) ;
#endif
  //  5.#cells_surf    6.#metas     7.#resistant  8.#resistant_surf
  fprintf(times,"%d %d %d %d   ",nsurf,lesions.size(),no_resistant,no_resistant_surf) ;
  // 9.#drivers   10.#cells_with_drv  11.#cells_with_drv_surf    12.#drv/cell   13.#der/cell_surf
  fprintf(times,"%d %d %d  %lf %lf  ",drivers.size(),cells_drv,cells_drv_surf,drv_per_cell,drv_per_cell_surf) ;
  // 14.growth_rate(n)   15.av_distance   16.pms_per_cell   17.snps_detected  18.<migr>
  fprintf(times,"%lf\t%f\t%lf %d\t %le\t",aver_growth_rate,average_distance_ij(),pms_per_cell,snps_det,av_migr) ;
  // #MBs   time_taken
  fprintf(times,"%d %f\n",memory_taken(),float(1.*(clock()-start_clock)/CLOCKS_PER_SEC)) ;
  if (treatment>0 || ntot>512 || ntot==max_size) fflush(times) ; // flush only when size big enough, this allows us to discard runs that died out

  if (ntot>256) { printf("%d %lf   no.les.=%d  no.res=%d drv_cell=%lf max_growth=%lf\n",ntot,tt,lesions.size(),no_resistant, drv_per_cell,max_growth_rate) ; fflush(stdout) ; }
}

void snps_corr(Hist *snps) ;
void snps_corr_cutoff(Hist *snps,float cutoff,int *snp_no) ;
void snps_corr_cond_driver(Hist *snps) ;
void find_p_driver(Hist *pr1, Hist *pr2, Hist *pr3) ;
void save_snp_corr(char *name, Hist *snps);

void save_spatial(int *snp_no)
{
#ifndef NO_MECHANICS  
  printf("save spatial\n") ;
  char tmp[256] ;
  Hist *snp_corr, *snp_corr_cutoff ;  // this is for measuring correlations between PMs in different parts of the tumor
  Hist *p_driver1, *p_driver2, *p_driver3, *snp_corr_cd ;

  snp_corr=new Hist[_bins] ; 
  snp_corr_cutoff=new Hist[_bins] ; 
  snp_corr_cd=new Hist[_bins] ; 

  snps_corr(snp_corr) ;
  sprintf(tmp,"%s/corr_%d_%d.dat",NUM,RAND,sample) ; save_snp_corr(tmp, snp_corr) ;                
  snps_corr_cutoff(snp_corr,0.1,snp_no) ;
  sprintf(tmp,"%s/cutoff01corr_%d_%d.dat",NUM,RAND,sample) ; save_snp_corr(tmp, snp_corr) ;                

  if (driver_adv>0 || driver_migr_adv>0) {
    p_driver1=new Hist[_bins] ; p_driver2=new Hist[_bins]; p_driver3=new Hist[_bins]; 
    find_p_driver(p_driver1,p_driver2,p_driver3) ;
    sprintf(tmp,"%s/P_driver1_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver1) ;
    sprintf(tmp,"%s/P_driver2_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver2) ;
    sprintf(tmp,"%s/P_driver3_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver3) ;
    snps_corr_cond_driver(snp_corr_cd) ;
    sprintf(tmp,"%s/corr_cond_driver_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, snp_corr_cd) ;
    delete [] p_driver1 ;    delete [] p_driver2 ;    delete [] p_driver3 ;  
  }

  delete [] snp_corr ;  delete [] snp_corr_cutoff ; delete [] snp_corr_cd ; 

  printf("done\n") ;
#endif
}

#ifdef PUSHING
void Lesion::find_dir_min_drag(int i, int j, int k, int &in, int &jn, int &kn)       // find direction of least drag
{
  int nn, in0,jn0,kn0 ;
  vector <IVec> vis ; // vector of visited sites. it will begin with (i,j,k) and end at empty site

rep:      
  vis.clear() ;
  vis.push_back(IVec(i,j,k)) ;
  in0=i ; jn0=j ; kn0=k ; 
      
  do {      // loop goes over subsequent pushing events
    float mind=wx ;
    nn=-1 ;
    for (int nnnn=0;nnnn<10;nnnn++) {
      int nnn=1+_drand48()*_nonn ;
      in=in0 ; jn=jn0 ; kn=kn0 ; 
      for (float drag=0;drag<mind;drag+=kln[nnn]) {
        in+=kz[nnn] ; jn+=ky[nnn] ; kn+=kx[nnn] ;
        for (int vv=0;vv<vis.size();vv++) if (vis[vv]==IVec(in,jn,kn)) goto brk ; // reject if trajectory passes through prev. visited sites
        if (p[in][jn][kn]==-1) { mind=drag ; nn=nnn ; break ; } 
      }
brk:  continue ;
    }
    if (nn==-1) goto rep ; 
        // now nn gives the direction of pushing
        
    in0+=kz[nn] ; jn0+=ky[nn] ; kn0+=kx[nn] ; // update position of the cell to be pushed
    vis.push_back(IVec(in0,jn0,kn0)) ; // and remember it...
  } while (p[in0][jn0][kn0]!=-1) ; // if the next position contains an empty site then exit

  // push all remembered cells except mother to make space for a single new daughter cell
  Sites sup=-1 ;        // sup is the elevated cell that needs to be inserted into new position
  for (int vv=1;vv<vis.size();vv++) {
    in0=vis[vv].i ; jn0=vis[vv].j ; kn0=vis[vv].k ;  
    Sites snew=p[in0][jn0][kn0] ;  
    p[in0][jn0][kn0]=sup ;
    if (sup!=-1) {
      Cell *c=&cells[sup] ; 
      c->x=kn0-wx/2 ; c->y=jn0-wx/2 ; c->z=in0-wx/2 ;
    }
    sup=snew ;
  }

  in=vis[1].i ; jn=vis[1].j ; kn=vis[1].k ; // the new cell will be the second position from the list (1st is the mother cell which is not pushed)
  if (p[in][jn][kn]!=-1) err("!!!") ;
}
#endif


//-----------------------------------------------------
#if defined(NORMAL)

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  int i,j,k,n,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;

  for(;;) {      // main loop 
#ifdef PAUSE_WHEN_MEMORY_LOW
    timeout++ ; if (timeout>1000000) {
      timeout=0 ; //printf("xxx") ;
      while (freemem()<PAUSE_WHEN_MEMORY_LOW) { sleep(1) ; } //printf("%d\n",freemem()) ; }
    }    
#endif
    double tsc=0.01*cells.size() ; if (tsc>1./max_growth_rate) tsc=1./max_growth_rate ;
    tt+=tsc*timescale/cells.size() ; 
    n=_drand48()*cells.size() ;
//    if (cells[n].lesion<0 || cells[n].lesion>=cells.size()) err("l=",cells[n].lesion) ;
    Lesion *ll=lesions[cells[n].lesion] ;
    int wx=ll->wx ; 
    k=cells[n].x+wx/2 ; j=cells[n].y+wx/2 ; i=cells[n].z+wx/2 ; 
    int need_wx_update=0 ;
    if (k<2 || k>=ll->wx-3 || j<2 || j>=ll->wx-3 || i<2 || i>=ll->wx-3) need_wx_update=1 ; 
#ifdef PUSHING
    if (ll->p[i][j][k]!=n) err("ll->p[i][j][k]!=n, p=",ll->p[i][j][k]) ;
#else
    if (ll->p[i*wx+j]->is_set(k)==0) err("ll->p[i][j][k]==0, wx=",wx) ;
#endif

    if (_drand48()<tsc*genotypes[cells[n].gen]->growth[treatment]) { // reproduction
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      int nn=1+int(_drand48()*_nonn) ;
      int in=(wx+i+kz[nn])%wx, jn=(wx+j+ky[nn])%wx, kn=(wx+k+kx[nn])%wx ;
#ifdef VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC // one more trial to find an empty site
      if (ll->p[in*wx+jn]->is_set(kn)==1) {
        nn=1+int(_drand48()*_nonn) ;
        in=(wx+i+kz[nn])%wx ; jn=(wx+j+ky[nn])%wx ; kn=(wx+k+kx[nn])%wx ;
      }
#endif
      if (ll->p[in*wx+jn]->is_set(kn)==0) {
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      int in=i, jn=j, kn=k ;
      ll->choose_nn(kn,jn,in) ;
      if (kn!=-1000000) { // if there is at least one empty n.n., then.....
#elif defined(PUSHING)
      int in,jn,kn ;
      ll->find_dir_min_drag(i,j,k, in,jn,kn) ;
      {
#else
  #error inconsistent growth conditions 
#endif
        int no_SNPs=poisson() ; // newly produced cell mutants
        if (_drand48()>genotypes[cells[n].gen]->m[treatment]) { // make a new cell in the same lesion
          Cell c ; c.x=kn-wx/2 ; c.y=jn-wx/2 ; c.z=in-wx/2 ; c.lesion=cells[n].lesion ;
#ifndef PUSHING
          ll->p[in*wx+jn]->set(kn) ;
#else
          ll->p[in][jn][kn]=cells.size() ;
#endif
          if (no_SNPs>0) { 
            c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // mutate 
          } else { 
            c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
          }
          cells.push_back(c) ; volume++ ;
          ll->n++ ; 
#ifndef NO_MECHANICS
          double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
          if (ll->rad/ll->rad0>1.05) {
            ll->reduce_overlap() ;  
            ll->find_closest() ; 
            ll->rad0=ll->rad ;
            ll->n0=ll->n ; 
          }
#endif
        } else { // make a new lesion
          int x=kn-wx/2+ll->r.x, y=jn-wx/2+ll->r.y, z=in-wx/2+ll->r.z ;
          if (no_SNPs>0) { 
            genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
            lesions.push_back(new Lesion(cells.size(),genotypes.size()-1,x,y,z)) ;
          } else {
            genotypes[cells[n].gen]->number++ ; 
            lesions.push_back(new Lesion(cells.size(),cells[n].gen,x,y,z)) ;
          }        
#ifndef NO_MECHANICS
          lesions[lesions.size()-1]->find_closest() ; 
#endif
        }
// BOTH_MUTATE          
        no_SNPs=poisson() ; // old cell mutates
        if (no_SNPs>0) { 
          genotypes[cells[n].gen]->number-- ; 
          int pn=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
          cells[n].gen=genotypes.size()-1 ;
          if (genotypes[cells[n].gen]->number<=0) { 
            delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
          }
        }
      }
    }
#ifdef CORE_IS_DEAD
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ; 
    }
#endif

// now we implement death
//    if (tsc>1./max_growth_rate) tsc=1./max_growth_rate ; // this is an alternative way but it does not change anything so no need to use it
#ifdef DEATH_ON_SURFACE    
    if (genotypes[cells[n].gen]->death[treatment]>0 && _drand48()<tsc*genotypes[cells[n].gen]->death[treatment]*ll->no_free_sites(k,j,i)/float(_nonn))  { // death on the surface
#else
    if (_drand48()<tsc*genotypes[cells[n].gen]->death[treatment]) { // death in volume
#endif
#ifndef PUSHING
      ll->p[i*wx+j]->unset(k) ;
#else
      ll->p[i][j][k]=-1 ;
#endif
      ll->n-- ; 
//      if (ll->n<0) err("ll->n<0") ;
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius
        ll->rad=0 ; 
        for (i=0;i<wx;i++) for (j=0;j<wx;j++) for (k=0;k<wx;k++) {
#ifndef PUSHING
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i*wx+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#else
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i][j][k]!=-1 && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#endif
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) {
        int nn=cells[n].lesion ;
//        if (ll!=lesions[nn]) err("ll!") ;
        ll=NULL ; 
        delete lesions[nn] ; 
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (i=0;i<cells.size();i++) if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ; 
        }        
        lesions.pop_back() ;        
#ifndef NO_MECHANICS
        for (i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;  
        }
#endif 
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      if (n!=cells.size()-1) { 
        cells[n]=cells[cells.size()-1] ;
#ifdef PUSHING
        Lesion *ll2=lesions[cells[n].lesion] ;
        int ii=cells[n].z+ll2->wx/2, jj=cells[n].y+ll2->wx/2, kk=cells[n].x+ll2->wx/2 ;
        ll2->p[ii][jj][kk]=n ;
#endif
      }
      cells.pop_back() ; volume-- ;
      //if (lesions.size()==0) err("N=",int(cells.size())) ;
    }

    if (need_wx_update && ll!=NULL) ll->update_wx() ;    
      
#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif

    if (wait_time>0 && tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }

    if (cells.size()==0) return 1 ; 
    if (max_time>0 && tt>max_time) return 3 ;
    if (exit_size>0 && ntot>=exit_size) return 4 ;

  }

}
#endif // NORMAL


//-----------------------------------------------------------------------------

#if defined(FASTER_KMC) || defined(GILLESPIE)

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  int i,j,k,n,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;

  for(;;) {      // main loop 
#ifdef PAUSE_WHEN_MEMORY_LOW
    timeout++ ; if (timeout>1000000) {
      timeout=0 ; 
      while (freemem()<PAUSE_WHEN_MEMORY_LOW) { sleep(1) ; } 
    }    
#endif

#ifdef GILLESPIE
  // Gillespie
    double tot_rate=0 ;
    for (n=0;n<cells.size();n++) {
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      tot_rate+=genotypes[cells[n].gen]->growth[treatment] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      if (free_sites(n)>0) tot_rate+=genotypes[cells[n].gen]->growth[treatment] ;
#elif defined(PUSHING)
      tot_rate+=genotypes[cells[n].gen]->growth[treatment] ;
#else
      #error inconsistent growth conditions 
#endif
#ifdef DEATH_ON_SURFACE    
      tot_rate+=genotypes[cells[n].gen]->death[treatment] * free_sites(n)/float(_nonn) ; // death on the surface
#else
      tot_rate+=genotypes[cells[n].gen]->death[treatment] ;  // death in volume
#endif
    }        
    tt+=-log(1-_drand48())*timescale/tot_rate ; 
    double q=_drand48()*tot_rate,r=0 ;
    int mode=0 ;
    for (n=0;n<cells.size();n++) {
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      r+=genotypes[cells[n].gen]->growth[treatment] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      if (free_sites(n)>0) r+=genotypes[cells[n].gen]->growth[treatment] ;
#elif defined(PUSHING)
      r+=genotypes[cells[n].gen]->growth[treatment] ;
#else
      #error inconsistent growth conditions 
#endif
      if (r>q) break ;
#ifdef DEATH_ON_SURFACE    
      r+=genotypes[cells[n].gen]->death[treatment] * free_sites(n)/float(_nonn)  ; // death on the surface
#else
      r+=genotypes[cells[n].gen]->death[treatment] ;  // death in volume
#endif
      if (r>q) { mode=1 ; break ; }
    }
    if (n==cells.size()) err("n==cells.size() at t=",tt) ;
#endif

#ifdef FASTER_KMC
    double max_death_rate=1 ;
    double tot_rate=cells.size()*(max_growth_rate+max_death_rate) ;
    tt+=-log(1-_drand48())*timescale/tot_rate ; 
    n=_drand48()*cells.size() ;
    double q=_drand48()*(max_growth_rate+max_death_rate), br,dr  ;
    int mode=0 ; 
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
    br=genotypes[cells[n].gen]->growth[treatment] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
    if (free_sites(n)>0) br=genotypes[cells[n].gen]->growth[treatment] ; else br=0 ;
#elif defined(PUSHING)
    br=genotypes[cells[n].gen]->growth[treatment] ;
#else
      #error inconsistent growth conditions 
#endif
#ifdef DEATH_ON_SURFACE    
    dr=genotypes[cells[n].gen]->death[treatment] * free_sites(n)/float(_nonn) ; // death on the surface
#else
    dr=genotypes[cells[n].gen]->death[treatment] ;  // death in volume
#endif
    if (q<br) mode=0 ;
    else if (q<br+dr) mode=1 ;
    else mode=2 ;
#endif    

    Lesion *ll=lesions[cells[n].lesion] ;
    int wx=ll->wx ; 
    k=cells[n].x+wx/2 ; j=cells[n].y+wx/2 ; i=cells[n].z+wx/2 ; 
    int need_wx_update=0 ;
    if (k<2 || k>=ll->wx-3 || j<2 || j>=ll->wx-3 || i<2 || i>=ll->wx-3) need_wx_update=1 ; 

    if (mode==0) { // reproduction
#if !defined(PUSHING)
      int in=i, jn=j, kn=k ;
      ll->choose_nn(kn,jn,in) ;
#else
      int in,jn,kn ;
      ll->find_dir_min_drag(i,j,k, in,jn,kn) ;      
#endif

      int no_SNPs=poisson() ; // newly produced cell mutants
      if (_drand48()>genotypes[cells[n].gen]->m[treatment]) { // make a new cell in the same lesion
        Cell c ; c.x=kn-wx/2 ; c.y=jn-wx/2 ; c.z=in-wx/2 ; c.lesion=cells[n].lesion ;
#ifdef PUSHING
        ll->p[in][jn][kn]=cells.size() ;
#else
        ll->p[in*wx+jn]->set(kn) ;
#endif
        if (no_SNPs>0) { 
          c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // mutate 
        } else { 
          c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
        }
        cells.push_back(c) ; volume++ ;

        ll->n++ ; 
#ifndef NO_MECHANICS
        double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
        if (ll->rad/ll->rad0>1.05) {
          ll->reduce_overlap() ;  
          ll->find_closest() ; 
          ll->rad0=ll->rad ;
          ll->n0=ll->n ; 
        }
#endif
      } else { // make a new lesion
        int x=kn-wx/2+ll->r.x, y=jn-wx/2+ll->r.y, z=in-wx/2+ll->r.z ;
        if (no_SNPs>0) { 
          genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
          lesions.push_back(new Lesion(cells.size(),genotypes.size()-1,x,y,z)) ;
        } else {
          genotypes[cells[n].gen]->number++ ; 
          lesions.push_back(new Lesion(cells.size(),cells[n].gen,x,y,z)) ;
        }        
#ifndef NO_MECHANICS
        lesions[lesions.size()-1]->find_closest() ; 
#endif
      }
// BOTH_MUTATE          
      no_SNPs=poisson() ; // old cell mutates
      if (no_SNPs>0) { 
        genotypes[cells[n].gen]->number-- ; 
        int pn=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
        cells[n].gen=genotypes.size()-1 ;
        if (genotypes[cells[n].gen]->number<=0) { 
          delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
        }
      }
    }
  
#ifdef CORE_IS_DEAD
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ; 
    }
#endif

// now we implement death
    if (mode==1) {
#ifdef PUSHING
      ll->p[i][j][k]=-1 ;
#else
      ll->p[i*wx+j]->unset(k) ;
#endif
      ll->n-- ; 
//      if (ll->n<0) err("ll->n<0") ;
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius
        ll->rad=0 ; 
        for (i=0;i<wx;i++) for (j=0;j<wx;j++) for (k=0;k<wx;k++) {
#ifdef PUSHING
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i][j][k]!=-1 && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#else
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i*wx+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#endif
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) {
        int nn=cells[n].lesion ;
//        if (ll!=lesions[nn]) err("ll!") ;
        ll=NULL ; 
        delete lesions[nn] ; 
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (i=0;i<cells.size();i++) if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ; 
        }        
        lesions.pop_back() ;        
#ifndef NO_MECHANICS
        for (i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;  
        }
#endif 
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      if (n!=cells.size()-1) { 
        cells[n]=cells[cells.size()-1] ;
#ifdef PUSHING
        Lesion *ll2=lesions[cells[n].lesion] ;
        int ii=cells[n].z+ll2->wx/2, jj=cells[n].y+ll2->wx/2, kk=cells[n].x+ll2->wx/2 ;
        ll2->p[ii][jj][kk]=n ;
#endif
      }
      cells.pop_back() ; volume-- ;
      //if (lesions.size()==0) err("N=",int(cells.size())) ;
    }

    if (need_wx_update && ll!=NULL) ll->update_wx() ;    
      
#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif

    if (wait_time>0 && tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }

    if (cells.size()==0) return 1 ; 
    if (max_time>0 && tt>max_time) return 3 ;
    if (exit_size>0 && ntot>=exit_size) return 4 ;

  }

}

#endif // FASTER_KMC or GILLESPIE