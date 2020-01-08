/* 
 Ensemble refinment with experimental data using LBGFS-optimization.

 INPUT;
 { Y } : Experimental data. M dimension.
 { y } : Simulation data projected to the observal space. M x N dimension.

 M : Number of the data points.
 N : Number of the snapshots.

 Where, y = y(x), and {x} means samples on phase space.

 PARAMETERS;
 Delta : the width of the error between <y> and Y.
 r_k : weight for the k-th exterior point method.

 OUTPUT;
 { omega }: weights for each snapshots that satisfies experimental data
            with the given error width. N dimension.
 { <y> }: Optimized observals. M dimension.

 Where, 
 rho_sim : distribution derived from simulation.
 rho     : optimal distribution that realize the experimental data within error width.

 Optimal distribution rho can be represented by weights { omega }.

 DIFINITIONS;
 rho = omega_i ... x=x_{i}.
       0       ... otherwise.

 omega_i == exp(g_i) / sum_i(exp(g_i)).

 ALGORITHM;
 { g^{opt.} } <- arg min_{ g } D(rho || rho_sim),
 subject to Chi^2 < Delta, and \int rho dx = 1.

 Where, 
 D means Kullback-Leiblar divergence and Chi := \int \rho(x) y(x) dx - Y.
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
  double *g0;           // initial weight for each snapshot
  double *y_ex,**y_sim; // observal (experiment), observal (simulation)
};

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *gi,
    lbfgsfloatval_t *g_Lambda, // derivative of ex Lagrangian
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

double cChiSquare(
		  int n,
		  int m,
		  double *exp_gi,
		  double **y_sim,
		  double *y_ex);

double cKLd(
    int n,
    double *exp_g0,
    double *exp_gi
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

  int n_snp;                     // # of snapshots, # of data points

  int n_ExPortM;                 // # of iterations of exterior point method
  double rk[M];                  // weight for penalty function
  char buf[256];
  
  double *Chi, ChiSqu=0.0;
  double *y_opt;
  struct Lex_parameters Lex_p;   // parametrs for extended Lagrangian

  double KLD;                    // KL divergence
  
  lbfgsfloatval_t Lex;           // extended Lagrangian
  lbfgsfloatval_t *gopt;         // (optimal) weights for snapshots
  lbfgs_parameter_t LBFGS_param; // parametrs for LBFGS

  double *exp_gi;               // expotential of g_i
  double sum_exp_gi;            // sum of exp(g_i)
  double *exp_g0;               // expotential of g_0
  
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
  Lex_p.g0 = (double *)emalloc(sizeof(double)*n_snp);
  exp_g0 = (double *)emalloc(sizeof(double)*n_snp);

  for (i=0;i<n_snp;++i) {
    exp_g0[i] = 1.0/n_snp;
    Lex_p.g0[i] = log(exp_g0[i]);
  }
  
  gopt = lbfgs_malloc(n_snp);
  exp_gi = (double *)emalloc(sizeof(double)*n_snp);
  
  for (i=0;i<n_snp;++i) {
    gopt[i] = Lex_p.g0[i];
    exp_gi[i] = exp(gopt[i]);
  }

  KLD=cKLd(n_snp,exp_g0,exp_gi);
  ChiSqu=cChiSquare(n_snp, Lex_p.n_pnt, exp_gi, Lex_p.y_sim, Lex_p.y_ex);
  
  printf("\nInitial Value of Chi^2 = %8.3lf KLD = %10.4lf\n\n",ChiSqu,KLD);

  printf("Start optimization\n");
  /* STEP 7-1: Start Exportal Points Method for optimization of xLagrangian */
  for (i=0;i<n_ExPortM;++i) {
    Lex_p.r=rk[i];

    /* STEP 7-2 Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&LBFGS_param);
    
    /* STEP 7-3 Start the L-BFGS optimization.*/
    if ((ret = lbfgs(n_snp, gopt, &Lex, evaluate, progress, &(Lex_p), &LBFGS_param)) != 0) {
      printf("lbfgs error!\n");
      exit(1);
    }

    for (j=0;j<n_snp;++j) exp_gi[j] = exp(gopt[j]);

    KLD=cKLd(n_snp,exp_g0,exp_gi);
    ChiSqu=cChiSquare(n_snp, Lex_p.n_pnt, exp_gi, Lex_p.y_sim, Lex_p.y_ex);

    printf("n=%3d  rho_%-3d = %5.2e Lex = %8.3lf Chi^2 = %8.3lf KLD = %12.4e\n", i+1, i+1, rk[i], Lex, ChiSqu, KLD);
    if ( ChiSqu < Lex_p.Delta) {
      printf("Optimization is done\n\n");
      break;
    }
  }
  printf("Optimization is finished\n\n");

  KLD=cKLd(n_snp,exp_g0,exp_gi);
  
  free(Lex_p.g0);

  sum_exp_gi = 0.0;
  for (i=0;i<n_snp;++i) {
    exp_gi[i] = exp(gopt[i]);
    sum_exp_gi += exp_gi[i];
  }
    
  y_opt=(double *)emalloc(sizeof(double)*Lex_p.n_pnt);
  for (a=0;a<Lex_p.n_pnt;++a) {
    y_opt[a]=0.0;
    for (i=0;i<n_snp;++i) {
      y_opt[a]+=exp_gi[i]*Lex_p.y_sim[i][a];
    }
    y_opt[a]=y_opt[a]/sum_exp_gi;
  }
  
  ChiSqu=cChiSquare(n_snp, Lex_p.n_pnt, exp_gi, Lex_p.y_sim, Lex_p.y_ex);

  free(Lex_p.y_ex);

  printf("Final Value of Chi^2 = %8.3lf KLD = %10.4lf\n",ChiSqu, KLD);

  /* FINAL STEP Write outputfile (optimal expectations)*/
  outputOPTOBLfile=efopen(outputOPTOBLfilename,"w");
  for (a=0;a<Lex_p.n_pnt;++a) {
    fprintf(outputOPTOBLfile,"%10.8lf\n",y_opt[a]);
  }
  fclose(outputOPTOBLfile);
    
  free(y_opt);

  free(Chi);
  lbfgs_free(gopt);

  /* FINAL STEP Write outputfile (weights)*/
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<n_snp;++i) {
    fprintf(outputfile,"%4d %10.8lf\n",i+1,exp_gi[i]/sum_exp_gi);
  }
  fclose(outputfile);
  
  free(Lex_p.y_sim);
  free(exp_gi);
  
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
  \chi_a \equiv ( \langle y_{sim, a} \rangle - y_{ex, a}).
  \chi^2 \equiv \sum_a^{M} \chi_a^2.

  The difinition of Lagrangian $\mathscr{L}^{ex}_n$ is
  \begin{eqution}
    \mathscr{L}^{ex}_n({g}) \equiv \sum_{i=1}^{n} e^{g_{i}} ( g_{i} - g_{0} - \ln (\sum_{j=1}^{n} e^{g_{j}}) ) \ \
    + rho_{n} (\max \[0, (\chi_i^2 - \Delta^{2}])^2.
  \end{equation}

  The derivative of the Lagrangian is
  \begin{eqution}
    \frac{\partial \mathscr{L}^{ex}_n({g})}{\partial g_{i}} \\
    = \frac{e^{g_{i}}}{\sum_{j=1}^{n} e^{g_{j}}} ( g_{i} - g_{0} - \ln (\sum_{j=1}^{n} e^{g_{j}}) ) \\
    - \frac{e^{g_{i}}}{(\sum_{j=1}^{n} e^{g_{j}})^{2}} (\sum_{j=1}^{N} e^{g_{j}} \\
         (g_{j} - g_{0} - \ln (\sum_{k=1}^{n} e^{g_{k}})) ) \\
    + 4\rho_{n} (\chi^2-Delta) \sum_{a=1}^{M} \chi_{a} \\
         (  \frac{e^{g_{i}} y_{a,i}}{\sum_{k=1}^{n} e^{g_{k}}} \\
          - \frac{\sum_{j=1}^{N}e^{g_{i}} e^{g_{j}} y_{a,j}}{\sum_{k=1}^{n} e^{g_{k}}}^{2} ),
  \end{equation}
  where, 3rd of rhs is active when $\chi^2 > \delta^{2}$.
*/

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *gi,
    lbfgsfloatval_t *g_Lambda,
    const int n,
    const lbfgsfloatval_t step
				) {
  int i, j, a;
  int m; // number of data points
  struct Lex_parameters *Lex_p=(struct Lex_parameters *)instance;
  lbfgsfloatval_t Le=0.0; // Lagrangian

  double *Chi, ChiSqu=0.0, Discriminat;
  double *g0, *exp_gi;
  double sum_exp_gi, sum_exp_gi_square, ln_sum_exp_gi;
  double g_Lambda_2nd_term, g_Lambda_3rd_term, f;

  m=Lex_p->n_pnt;

  Chi = (double *)emalloc(sizeof(double)*m);

  g0 = (double *)emalloc(sizeof(double)*n);
  exp_gi = (double *)emalloc(sizeof(double)*n);
  
  sum_exp_gi = 0.0;
  for (i=0;i<n;++i) {
    g0[i] = Lex_p->g0[i];
    exp_gi[i] = exp(gi[i]);
    sum_exp_gi += exp_gi[i];
  }
  sum_exp_gi_square = sum_exp_gi*sum_exp_gi;
  ln_sum_exp_gi = log(sum_exp_gi);
  
  g_Lambda_2nd_term = 0.0;
  for (i=0;i<n;++i) {
    g_Lambda_2nd_term += exp_gi[i]*(gi[i]-g0[i]-ln_sum_exp_gi);
  }
  
  Le = 0.0;
  for (i=0;i<n;++i) {
    g_Lambda_2nd_term = exp_gi[i]/sum_exp_gi_square*g_Lambda_2nd_term;
    g_Lambda[i]=exp_gi[i]/sum_exp_gi*(gi[i]-g0[i]-ln_sum_exp_gi) - g_Lambda_2nd_term;
    
    Le+=exp_gi[i]*(gi[i]-g0[i]-ln_sum_exp_gi);
  }
  Le = Le/sum_exp_gi;

  ChiSqu=0.0;
  for (a=0;a<m;++a) {
    Chi[a]=0.0;
    for (i=0;i<n;++i) {
      Chi[a]+=exp_gi[i]*Lex_p->y_sim[i][a];
    }
    Chi[a]=Chi[a]/sum_exp_gi;
    Chi[a]-=Lex_p->y_ex[a];
  
    ChiSqu += Chi[a]*Chi[a];
  }
  
  Discriminat = ChiSqu - Lex_p->Delta;
  
  if (Discriminat > 0.0) {

    for (i=0;i<n;++i) {
      g_Lambda_3rd_term=0.0;
      for (a=0;a<m;++a) {
	f = 0.0;
	for (j=0;j<n;++j) {
	  f += exp_gi[i]*exp_gi[j]*Lex_p->y_sim[j][a];
	}
	f = f/sum_exp_gi_square;
  	g_Lambda_3rd_term+=Chi[a]*(exp_gi[i]*Lex_p->y_sim[i][a]/sum_exp_gi-f);
      }
      g_Lambda_3rd_term = Lex_p->r*4.0*Discriminat*g_Lambda_3rd_term;
      g_Lambda[i]+=g_Lambda_3rd_term;
    }

    Le+=Lex_p->r*Discriminat*Discriminat;
  }
  else {
    printf("This is internal area!\n");
  }

  free(Chi);
  free(g0);
  free(exp_gi);

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

double cChiSquare(
     int n,
     int m,
     double *exp_gi,
     double **y_sim,
     double *y_ex)
{
  int i,a;
  double sum, ChiSqu;
  double *Chi;

  Chi=(double *)emalloc(sizeof(double)*m);
  ChiSqu=0.0;
  
  sum = 0.0;
  for (i=0;i<n;++i) {
    sum += exp_gi[i];
  }
  
  for (a=0;a<m;++a) {
    Chi[a]=0.0;
    for (i=0;i<n;++i) {
      Chi[a]+=exp_gi[i]*y_sim[i][a];
    }
    Chi[a]=Chi[a]/sum;
    Chi[a]-=y_ex[a];
  
    ChiSqu += Chi[a]*Chi[a];
  }

  return ChiSqu;
}

double cKLd(
    int n,
    double *exp_g0,
    double *exp_gi)
{
  int i;
  double KLD, sum, sum_0;

  sum = 0.0;
  for (i=0;i<n;++i) {
    sum += exp_gi[i];
  }

  sum_0 = 0.0;
  for (i=0;i<n;++i) {
    sum_0 += exp_g0[i];
  }

  for (i=0;i<n;++i) {
    KLD += exp_gi[i]*(log(exp_gi[i]) - log(sum) - log(exp_g0[0]) + log(sum_0));
  }
  KLD = KLD/sum;

  return KLD;
}
