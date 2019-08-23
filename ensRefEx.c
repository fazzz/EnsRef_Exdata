/* 
 Ensemble refinment with experimental data using LBGFS-optimization.

 INPUT;
 { Y } : Experimental data. M dimension.
 { y } : Simulation data projected to the observal space. M x N dimension.

 Where, y = y(x), and {x} means samples on phase space.

 OUTPUT;
 { omega }: weights for each snapshots that satisfies experimental data
            with the given error width.
            N dimension.
 { <y> }: Optimized observals. M dimension.

 Where, 
 rho_sim : distribution derived from simulation.
 rho     : optimal distribution that realize the experimental data
           within error width Delta.

 Optimal distribution rho can be represented by weights { omega }.

DIFINITIONS;
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
#include <time.h>
#include <getopt.h>

#include <lbfgs.h>

#include "error.h"

#define M 100 // Maximum number for the parameters of Exterior Point M.

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
    const lbfgsfloatval_t step
				);

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
		    );

void USAGE(char *progname);

void writeConditions(
		    char *progname,
		    char *version,
		    char *inputfileExpDataname,
		    char *inputfileSimDataname,
		    char *paramfilename,
		    char *outputfilename,
		    char *outputOPTOBLfilename,
		    int n_snp,
		    int n_pnt,
		    double rk[M],
		    double Delta,
		    int n_ExPortM,
		    char date[64]
		     );

int main(int argc, char *argv[]) {
  int i,j,k,l,a,d,n,ret = 0,interval=100;
  double dummy,f;

  int n_snp, N;   // # of snapshots, # of data points

  int n_ExPortM;  // # of iterations of exterior point method
  double rk[M];   // weight for penalty function
  char buf[256];
  
  double *wopt;                // initial and optimal weight for each snapshot
  double *Chi, ChiSqu=0.0;
  double *y_opt;
  struct Lex_parameters Lex_p; // parametrs for extended Lagrangian
  
  lbfgsfloatval_t Lex;           // extended Lagrangian
  lbfgsfloatval_t *gopt;         // (optimal) weights for snapshots
  lbfgs_parameter_t LBFGS_param; // parametrs for LBFGS

  char date[64];
  time_t t;
  
  char *inputfileExpDataname;
  char *inputfileSimDataname;
  char *outputfilename;
  char *outputOPTOBLfilename;
  char *paramfilename="epmparam";

  FILE *inputfileExpData;
  FILE *inputfileSimData;
  FILE *outputfile;
  FILE *outputOPTOBLfile;
  FILE *paramfile;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  char *version="XXXX";
  int opt_idx=1;
  
  Lex_p.Delta=1.0e1;

  /* STEP 0: Get options */
  struct option long_opt[] = {
    {"Del",1,NULL,'d'},
    {"Expm",1,NULL,'e'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hd:e:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      Lex_p.Delta=atof(optarg);  break;
    case 'e':
      paramfilename=optarg;  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  /* STEP 1: Get input/output filenames */
  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfileExpDataname = *argv;
  inputfileSimDataname = *++argv;
  outputfilename       = *++argv;
  outputOPTOBLfilename = *++argv;
  
  /* STEP 2: Read exterior point method parameters */
  paramfile=efopen(paramfilename,"r");
  n_ExPortM=0;
  while ( fgets(buf,sizeof(buf),paramfile)) {
    if (strncmp(buf,"//",2)==0 || strcmp(buf,"\n")==0) continue;
    n=sscanf(buf,"%lf%lf%lf%lf%lf",&rk[n_ExPortM],&rk[n_ExPortM+1],&rk[n_ExPortM+2],&rk[n_ExPortM+3],&rk[n_ExPortM+4]);
    n_ExPortM+=n;
  }
  fclose(paramfile);

  /* STEP 3: Read Experimental Data */
  Lex_p.y_ex=(double *)emalloc(sizeof(double)*1);

  inputfileExpData=efopen(inputfileExpDataname,"r");
  Lex_p.n_pnt=0;
  while ( (d=fscanf(inputfileExpData,"%lf",&f)) != EOF )  {
    Lex_p.y_ex[Lex_p.n_pnt]=f;
    ++(Lex_p.n_pnt);
    Lex_p.y_ex=(double *)realloc(Lex_p.y_ex,sizeof(double)*(Lex_p.n_pnt));
  }
  fclose(inputfileExpData);
  Chi = (double *)emalloc(sizeof(double)*Lex_p.n_pnt);

  /* STEP 4: Read simulation data */
  n_snp = 1;
  Lex_p.y_sim   =(double **)emalloc(sizeof(double *)*n_snp);
  Lex_p.y_sim[0]=(double  *)emalloc(sizeof(double)*Lex_p.n_pnt);
  
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

      Lex_p.y_sim=(double **)realloc(Lex_p.y_sim,sizeof(double *)*(n_snp));
      Lex_p.y_sim[n_snp-1]=(double *)emalloc(sizeof(double)*Lex_p.n_pnt);
      Lex_p.y_sim[n_snp-1][i] = f;

      ++i;
    }
  }
  fclose(inputfileSimData);

  time(&t);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  
  /* STEP 5: Print conditions for calculation */
  writeConditions(progname,version,
		  inputfileExpDataname,inputfileSimDataname,paramfilename,
		  outputfilename,outputOPTOBLfilename,
		  n_snp,Lex_p.n_pnt,rk,Lex_p.Delta,n_ExPortM,date);
  
  /* STEP 6: Set parameters for extended-Lagrangian */
  Lex_p.wi = (double *)emalloc(sizeof(double)*n_snp);
  Lex_p.gi = (double *)emalloc(sizeof(double)*n_snp);
  wopt = (double *)emalloc(sizeof(double)*n_snp);

  for (i=0;i<n_snp;++i) {
    Lex_p.wi[i] = 1.0/n_snp;
    Lex_p.gi[i] = log(Lex_p.wi[i]);
  }

  for (i=0;i<n_snp;++i) wopt[i] = Lex_p.wi[i];
  
  N = n_snp;
  gopt = lbfgs_malloc(N);
  
  for (i=0;i<n_snp;++i) gopt[i] = log(wopt[i]);

  ChiSqu=0.0;
  for (a=0;a<Lex_p.n_pnt;++a) {
    Chi[a]=0.0;
    for (i=0;i<n_snp;++i) {
      Chi[a]+=wopt[i]*Lex_p.y_sim[i][a];
    }
    Chi[a]-=Lex_p.y_ex[a];
  
    ChiSqu += Chi[a]*Chi[a];
  }
  
  printf("\nInitial Value of Chi^2 = %8.3lf\n\n",ChiSqu);

  printf("Start optimization\n");
  /* STEP 7-1: Start Exportal Points Method for optimization of xLagrangian */
  for (i=0;i<n_ExPortM;++i) {
    Lex_p.r=rk[i];

    /* STEP 7-2 Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&LBFGS_param);
    
    /* STEP 7-3 Start the L-BFGS optimization.*/
    if ((ret = lbfgs(N, gopt, &Lex, evaluate, progress, &(Lex_p), &LBFGS_param)) != 0) {
      printf("lbfgs error!\n");
      exit(1);
    }
    printf("n=%3d  rho_%-3d = %8.3e Lex = %8.3lf\n", i+1, i+1, rk[i], Lex);
  }
  printf("Optimization is done\n\n");

  free(Lex_p.y_ex);
  free(Lex_p.wi);
  free(Lex_p.gi);
  
  for (i=0;i<n_snp;++i) wopt[i] = exp(gopt[i]);

  y_opt=(double *)emalloc(sizeof(double)*Lex_p.n_pnt);
  
  ChiSqu=0.0;
  for (a=0;a<Lex_p.n_pnt;++a) {
    y_opt[a]=0.0;
    for (i=0;i<n_snp;++i) {
      y_opt[a]+=wopt[i]*Lex_p.y_sim[i][a];
    }
    Chi[a] = y_opt[a]-Lex_p.y_ex[a];
  
    ChiSqu += Chi[a]*Chi[a];
  }

  printf("Final Value of Chi^2 = %8.3lf\n",ChiSqu);

  outputOPTOBLfile=efopen(outputOPTOBLfilename,"w");
  for (a=0;a<Lex_p.n_pnt;++a) {
    fprintf(outputOPTOBLfile,"%10.8lf\n",y_opt[a]);
  }
  fclose(outputOPTOBLfile);
    
  free(y_opt);

  free(Chi);
  free(Lex_p.y_sim);
  lbfgs_free(gopt);
  
  /* FINAL STEP Write outputfile */
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<n_snp;++i) {
    fprintf(outputfile,"%4d %10.8lf\n",i+1,wopt[i]);
  }
  fclose(outputfile);

  free(wopt);
  
  time(&t);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("  Fin.  Data : %s\n",date);
  printf("  Calculation is done\n");
  
  return 0;
}

/* Show usage */
void USAGE(char *progname) {
  printf("USAGE: ");
  printf("%s [ options ] input1(exp data) input2(simu data) outout1(weights) output2(ave)\n",progname);
  printf("Options: \n");
  printf("[-h] help \n");
  printf("[-Del <Delta Value> ] \n");
  printf("[-e   <param file name> ] \n");
}

/* Show conditions */
void writeConditions(
		     char *progname,
		     char *version,
		     char *inputfileExpDataname,
		     char *inputfileSimDataname,
		     char *paramfilename,
		     char *outputfilename,
		     char *outputOPTOBLfilename,
		     int n_snp,
		     int n_pnt,
		     double rk[M],
		     double Delta,
		     int n_ExPortM,
		     char date[64]) {

  printf("###############################################################\n");
  printf("###############################################################\n");
  printf("### %-20s           version %-7s :         ##\n",progname,version);
  printf("###############################################################\n");
  printf("  Start Date : %s\n\n",date);

  printf(" FILE name(input1) : %s\n",inputfileExpDataname);
  printf(" FILE name(input2) : %s\n",inputfileSimDataname);
  printf(" FILE name(output1) : %s\n",outputfilename);
  printf(" FILE name(output2) : %s\n\n",outputOPTOBLfilename);

  printf(" ExPortM parameters are given in %s\n",paramfilename);
  printf(" Delta : %8.3lf\n\n",Delta);
  
  printf("   num. of data points : %4d\n",n_pnt);
  printf("   num. of snpshots : %4d\n",n_snp);
  printf("num. of iterations of exterior point method : %4d\n", n_ExPortM);

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
  int i, j, n_snp, a;
  struct Lex_parameters *Lex_p=(struct Lex_parameters *)instance;
  lbfgsfloatval_t Le=0.0;

  double *Chi, ChiSqu=0.0, Accep;
  double *wopt;
  double sum_w=0.0, f;

  Chi = (double *)emalloc(sizeof(double)*Lex_p->n_pnt);

  n_snp = n;
  wopt = (double *)emalloc(sizeof(double)*n_snp);
  
  for (i=0;i<n_snp;++i) {
    wopt[i]=exp(gopt[i]);
    sum_w+=wopt[i];
  }

  Le=0.0;
  for (i=0;i<n_snp;++i) {
    g[i]=wopt[i]*(gopt[i]-Lex_p->gi[i]+1.0+2.0*Lex_p->r*(sum_w - 1.0));

    Le+=wopt[i]*log(wopt[i]/Lex_p->wi[i]);
  }
  Le += Lex_p->r*(sum_w - 1.0)*(sum_w - 1.0);

  ChiSqu=0.0;
  for (a=0;a<Lex_p->n_pnt;++a) {
    Chi[a]=0.0;
    for (i=0;i<n_snp;++i) {
      Chi[a]+=wopt[i]*Lex_p->y_sim[i][a];
    }
    Chi[a]-=Lex_p->y_ex[a];
  
    ChiSqu += Chi[a]*Chi[a];
  }
  
  Accep = ChiSqu - Lex_p->Delta;
  
  if (Accep > 0.0) {
    for (i=0;i<n_snp;++i) {
      f=0.0;
      for (a=0;a<Lex_p->n_pnt;++a) {
  	f+=Chi[a]*Lex_p->y_sim[i][a];
      }
      f=wopt[i]*Lex_p->r*4.0*Accep*f;
      g[i]+=f;
    }
    Le+=Lex_p->r*Accep*Accep;
  }

  free(Chi);
  free(wopt);
  
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
  return 0;
}
