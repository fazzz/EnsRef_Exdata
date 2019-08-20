/* 
 Ensemble refinment with experimental data.

 INPUT;
 { Y } : experimental data. M dimension.
 { y } : simulation data projected to the observal space. M x N dimension.

 Where y = y(x), and {x} means samples on phase space.

 OUTPUT;
 { omega }: weights for each snapshots.
            N dimension.

 Where, 
 rho_sim : distribution derived from simulation.
 rho     : optimal distribution that realize the experimental data
           within error width Delta.

 Optimal distribution rho can be represented by weights { omega }.

 rho = omega_i ... x=x_{i}.
       0       ... otherwise.

 g_i = ln omega_i.

 ALGORITHM;
 { g^{opt.} } <- arg min_{ g } D(rho || rho_sim),
 subject to Chi^2 < Delta, and \int rho dx = 1.

 Where, 
 D means Kullback-Leiblar divergence and
 Chi := \int \rho y dx - Y.
*/

#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <lbfgs.h>

#include "error.h"
//#include "func.h"

#define M 10

#define LargeConst 1.0e10

#define ON 0
#define OFF 1

struct Lex_parameters {
  int n_pnt;            // # of points
  double Delta;         // width of error
  double r;             // weight for penalty function
  double *wi,*gi;       // initial weight for each snapshot
  double *y_ex,**y_sim; // observal (experiment), observal (simulation)
};

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step);

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls);

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,ret = 0,interval=100;
  double dummy,f;

  int n_snp,/* n_pnt,*/ N;   // # of snapshots, # of data points

  double rk[M];            // weight for penalty function
  
  //  double Del=1.0;          // width of error
  //  double r;                // weight for penalty function
  double /**wi, *gi,*/ *wopt;  // initial and optimal weight for each snapshot
  //  double *y_ex,**y_sim;    // observal (experiment), observal (simulation)
  double *Chi, ChiSqu=0.0;
  struct Lex_parameters Lex_p;    // parametrs for extended Lagrangian
  
  lbfgsfloatval_t Lex;     // extended Lagrangian
  lbfgsfloatval_t *gopt;   // (optimal) weights for snapshots
  lbfgs_parameter_t param; // parametrs for LBFGS

  char *inputfileExpDataname;
  char *inputfileSimDataname;
  char *outputfilename;

  FILE *inputfileExpData;
  FILE *inputfileSimData;
  FILE *outputfile;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  Lex_p.Delta=1.0;

  rk[0]=1.0e-1; rk[1]=1.0e0; rk[2]=1.0e1; rk[3]=1.0e2; rk[4]=1.0e3;
  rk[5]=1.0e4;  rk[6]=5.0e4; rk[7]=1.0e5; rk[8]=5.0e5; rk[9]=1.0e6;
  
  struct option long_opt[] = {
    {"Del",1,NULL,'d'},
    {"r1",1,NULL,'1'},
    {"r2",1,NULL,'2'},
    {"r3",1,NULL,'3'},
    {"r4",1,NULL,'4'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hd:1:2:3:4:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      Lex_p.Delta=atof(optarg);  break;
    case '1':
      rk[0]=atof(optarg);  break;
    case '2':
      rk[1]=atof(optarg);  break;
    case '3':
      rk[2]=atof(optarg);  break;
    case '4':
      rk[3]=atof(optarg);  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfileExpDataname = *argv;
  inputfileSimDataname = *++argv;
  outputfilename       = *++argv;

  //  printf("inputfile of Experimental Data=%s\n",inputfileExpDataname);
  //  printf("inputfile of Simulation   Data=%s\n",inputfileSimDataname);
  //  printf("outputfile      =%s\n",outputfilename);
  
  Lex_p.y_ex=(double *)gcemalloc_atomic(sizeof(double)*1);

  inputfileExpData=efopen(inputfileExpDataname,"r");
  Lex_p.n_pnt=0;
  while ( (d=fscanf(inputfileExpData,"%lf",&f)) != EOF )  {
    Lex_p.y_ex[Lex_p.n_pnt]=f;
    ++(Lex_p.n_pnt);
    Lex_p.y_ex=(double *)gcerealloc(Lex_p.y_ex,sizeof(double)*(Lex_p.n_pnt));
  }
  fclose(inputfileExpData);
  Chi = (double *)gcemalloc_atomic(sizeof(double)*Lex_p.n_pnt);
  //  printf("# of data points : %4d\n",n_pnt);

  /* debug printf */
  /*********************************/
  /* for (i=0;i<n_pnt;++i) {	   */
  /*   printf("%8.3lf\n",y_ex[i]); */
  /* }				   */
  /* printf("\n");		   */
  /*********************************/

  n_snp = 1;
  Lex_p.y_sim   =(double **)gcemalloc_atomic(sizeof(double *)*n_snp);
  Lex_p.y_sim[0]=(double  *)gcemalloc_atomic(sizeof(double)*Lex_p.n_pnt);
  
  inputfileSimData=efopen(inputfileSimDataname,"r");
  i=0;
  d = 1;
  while ( (d=fscanf(inputfileSimData,"%lf",&f)) != EOF )  {
    if (i<Lex_p.n_pnt) {
      Lex_p.y_sim[n_snp-1][i] = f;

      ++i;
    }
    else {
      ++n_snp;
      i=0;

      Lex_p.y_sim=(double **)gcerealloc(Lex_p.y_sim,sizeof(double *)*(n_snp));
      Lex_p.y_sim[n_snp-1]=(double *)gcemalloc_atomic(sizeof(double)*Lex_p.n_pnt);
      Lex_p.y_sim[n_snp-1][i] = f;

      ++i;
    }
  }
  fclose(inputfileSimData);

  //  printf("# of snap shots  : %4d\n",n_snp);
  
  /* debug printf */
  /***************************************/
  /* for (i=0;i<n_snp;++i) {		 */
  /*   for (j=0;j<n_pnt;++j) {		 */
  /*     printf("%8.3lf\n",y_sim[i][j]); */
  /*   }				 */
  /*   printf("\n");			 */
  /* }					 */
  /***************************************/

  Lex_p.wi = (double *)gcemalloc_atomic(sizeof(double)*n_snp);
  Lex_p.gi = (double *)gcemalloc_atomic(sizeof(double)*n_snp);
  wopt = (double *)gcemalloc_atomic(sizeof(double)*n_snp);

  for (i=0;i<n_snp;++i) {
    Lex_p.wi[i] = 1.0/n_snp;
    Lex_p.gi[i] = log(Lex_p.wi[i]);
  }

  for (i=0;i<n_snp;++i) {
    wopt[i] = Lex_p.wi[i];
  }

  N = n_snp;
  gopt = lbfgs_malloc(N);

  for (i=0;i<n_snp;++i) {
    gopt[i] = log(wopt[i]);
  }

  /*Lex_p.n_pnt=n_pnt;
  Lex_p.Delta=Del;
  Lex_p.wi=(double *)gcemalloc_atomic(sizeof(double)*n_snp);
  for (i=0;i<n_snp;++i) Lex_p.wi[i]=wi[i];
  Lex_p.gi=(double *)gcemalloc_atomic(sizeof(double)*n_snp);
  for (i=0;i<n_snp;++i) Lex_p.gi[i]=gi[i];
  Lex_p.y_ex=(double *)gcemalloc_atomic(sizeof(double)*n_pnt);
  for (i=0;i<n_pnt;++i) Lex_p.y_ex[i]=y_ex[i];
  Lex_p.y_sim=(double **)gcemalloc_atomic(sizeof(double *)*n_snp);
  for (i=0;i<n_snp;++i) Lex_p.y_sim[i]=(double *)gcemalloc_atomic(sizeof(double)*n_pnt);
  for (i=0;i<n_snp;++i) for (j=0;j<n_pnt;++j) Lex_p.y_sim[i][j]=y_sim[i][j];*/

  /* Initialize the parameters for the L-BFGS optimization. */
  lbfgs_parameter_init(&param);

  for (i=0;i<n_snp;++i) {
    for (j=0;j<Lex_p.n_pnt;++j) {
      printf("%8.3lf\n",Lex_p.y_sim[i][j]);
    }
  }
  
  for (i=0;i<M;++i) {
    
    Lex_p.r=rk[i];
    
    /*Start the L-BFGS optimization.*/
    if ((ret = lbfgs(N, gopt, &Lex, evaluate, progress, &(Lex_p), &param)) != 0) {
      printf("lbfgs error!\n");
      exit(1);
    }
    printf("n=%3d  rho_%-3d = %8.3e Lex = %8.3lf\n", i+1, i+1, rk[i], Lex);

  }

  for (i=0;i<n_snp;++i) {
    wopt[i] = exp(gopt[i]);
  }

  for (i=0;i<n_snp;++i) {
    for (j=0;j<Lex_p.n_pnt;++j) {
      printf("%8.3lf\n", Lex_p.y_sim[i][j]);
    }
  }
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<n_snp;++i) {
    fprintf(outputfile,"%4d %10.8lf\n",i+1,wopt[i]);
  }
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[-Del <Delta Value> ] \n");
  printf("[-r1  <rho1> ] \n");
  printf("[-r2  <rho2> ] \n");
  printf("[-r3  <rho3> ] \n");
  printf("[-r4  <rhor> ] \n");
  printf("%s inputname1(exp data) inputname2(simu data) outfilename\n",progname);
}

/*
  Chi_i := (<y_sim,i> - y_ex,i).
  Chi-2 := \Sum Chi_i^2.

  Lex^n({ omega }, Lambda) := Sum_{i=1}^{n} omega_{i} ln (omega_{i} / omega_{i}^{0}) 
                             + Lambda (Sum_{i=1}^{N} omega_{i} -1).
                             + rho_n (max[0, (\Sum Chi_i^2 - Delta])^2).

  g_i({omega}, Lambda) = ln (omega_{i} / omega_{i}^{0}) + 1
                         + 4*rho_n*(Chi^2-Delta)*\Sum_j=1{ Chi_j*y_j,i }.
  Where, 3rd of rhs is active when Chi^2>Delta.
  g_Lambda({omega}) = Sum_{i=1}^{N} omega_{i} -1.
*/

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *gopt,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
				) {
  int i, j, n_snp;
  struct Lex_parameters *Lex_p=(struct Lex_parameters *)instance;
  lbfgsfloatval_t Le=0.0;

  int n_pnt=Lex_p->n_pnt;
  double Delta=Lex_p->Delta;
  double r=Lex_p->r;
  double *wi;
  double *gi;
  double *y_ex;
  double **y_sim;
  
  double *Chi, ChiSum=0.0, ChiSqu=0.0, ratio_w_wi, ln_w_wi, Accep;
  double *wopt;
  double sum_w=0.0;

  wi=(double *)gcemalloc_atomic(sizeof(double)*n_snp);
  gi=(double *)gcemalloc_atomic(sizeof(double)*n_snp);

  y_ex=(double *)gcemalloc_atomic(sizeof(double)*n_pnt);
  y_sim=(double **)gcemalloc_atomic(sizeof(double *)*n_pnt);
  for (i=0;i<n_pnt;++i) y_sim[i]=(double *)gcemalloc_atomic(sizeof(double)*n_snp);
  
  Chi = (double *)gcemalloc_atomic(sizeof(double)*n_pnt);

  n_snp = n;
  wopt = (double *)gcemalloc_atomic(sizeof(double)*n_snp);
  
  for (i=0;i<n_snp;++i) {
    wopt[i]=exp(gopt[i]);
    sum_w+=wopt[i];
  }

  Le=0.0;
  for (i=0;i<n_snp;++i) {
    g[i]=wopt[i]*(gopt[i]-gi[i]+1.0+2.0*r*(sum_w - 1.0));

    Le+=wopt[i]*log(wopt[i]/wi[i]);
  }
  Le += r*(sum_w - 1.0)*(sum_w - 1.0);

  for (i=0;i<n_pnt;++i) {
    Chi[i]=0.0;
    for (j=0;j<n_snp;++j) {
      Chi[i]+=wopt[j]*y_sim[j][i];
    }
    Chi[i]-=y_ex[i];
    
    ChiSum += Chi[i];
    ChiSqu += Chi[i]*Chi[i];
  }

  Accep = ChiSqu - Delta;
  
  /******************************************************/
  /* if (Accep > 0.0) {				        */
  /*   for (i=0;i<n/\*-1*\/;++i) {		        */
  /*     for (j=0;j<n_pnt;++j) {		        */
  /* 	g[i]+=4.0*r*Accep*ChiSum*y_sim[i][j];	        */
  /*     }					        */
  /*   }					        */
  /*   Le+=r*Accep*Accep;			        */
  /* }						        */
  /******************************************************/

  return Le;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
		    )
{
  //  printf("Iteration %3d:", k);
  //  printf("  fx = %8.3f", fx);
  //  printf("  xnorm = %8.3e, gnorm = %8.3e, step = %8.3f\n", xnorm, gnorm, step);
  //  printf("\n");
  return 0;
}
