/* 
 Ensemble refinment with experimental data.

 INPUT;
 { Y } : experimental data. M dimension.
 { y } : simulation data projected to the observal space. M x N dimension.

 OUTPUT;
 omega: weights for each snapshots.

 where,
     rho_sim : distribution derived from simulation.
     rho_sim : optimal distribution.

     Chi == \int \rho y - Y.

 { omega^{opt.} } <- arg min_{ omega } D(rho || rho_sim),
 subject to Chi^2 < Delta, and Sigma (omega) = 1.
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

#define M 5

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,ret = 0;
  double dummy,f;

  int n_snp, n_pnt=1, N; // # of snapshots, # of data points

  double rk[4]; // weight for penalty function
  
  double Del=10e-2;     // width of error
  double r;             // weight for penalty function
  double *wi0;          // initial weight for each snapshot
  double *y_ex,**y_sim; // observal (experiment), observal (simulation)

  Lex_parameters Lex_p; // parametrs for extended Lagrangian
  
  lbfgsfloatval_t Lex;     // extended Lagrangian
  lbfgsfloatval_t *wopt;   // (optimal) weights for snapshots
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

  rk[0]=1.0e-1; rk[1]=1.0e0; rk[2]=1.0e1; rk[3]=1.0e2;
    
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
      Del=atof(optarg);  break;
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
  inputfileExpDataname = *++argv;
  inputfileSimDataname = *++argv;
  outputfilename       = *++argv;

  y_ex=(double *)gcemalloc(sizeof(double)*1);
  
  inputfileExpData=efopen(inputfileExpDataname,"r");
  n_pnt=0;
  while ( (d=fscanf(inputfileExpData,"%lf",&f)) != EOF )  {
    y_ex[n_pnt]=f;
    ++n_pnt;
    y_ex=(double *)gcerealloc(y_ex,sizeof(double)*(n_pnt));
  }
  fclose(inputfileExpData);

  y_sim   =(double **)gcemalloc(sizeof(double *)*1);
  y_sim[0]=(double  *)gcemalloc(sizeof(double)*n_pnt);
  
  inputfileSimData=efopen(inputfileSimDataname,"r");
  i=0;
  d = 1;
  n_snp = 0;
  while ( (d=fscanf(inputfileSimData,"%lf",&f)) != EOF )  {
    if (i<n_pnt) {
      y_sim[n_snp][i] = f;

      ++i;
    }
    else {
      y_sim=(double **)gcerealloc(y_sim,sizeof(double *)*(n_snp+1));
      y_sim[n_snp]=(double *)gcemalloc(sizeof(double)*n_pnt);
      ++n_snp;

      y_sim[n_snp][0] = f;
      i=1;
    }
  }
  fclose(inputfileSimData);

  if ( i != n_pnt ) {
    n_pnt -= 1;
  }
  
  printf("# of data points : %4d\n",n_pnt);
  printf("# of snap shots  : %4d\n",n_snp);

  /* debug printf */
  for (i=0;i<n_pnt;++i) {
    printf("%8.3lf\n",y_ex[i]);
  }
  printf("\n");
  
  /* debug printf */
  for (i=0;i<n_pnt;++i) {
    for (j=0;j<n_snp;++j) {
      printf("%8.3lf\n",y_sim[i][j]);
    }
    printf("\n");
  }
  
  wi0 = (double *)gcemalloc(sizeof(double)*n_snp);

  for (i=0;i<n_snp;++i) {
    wi0[i] = 1.0/n_snp;
  }

  wopt = lbfgs_malloc(n_snp+1);

  for (i=0;i<n_snp;++i) {
    wopt[i] = wi0[i];
  }
  wopt[i] = 0.0;
  
  /* Initialize the parameters for the L-BFGS optimization. */
  lbfgs_parameter_init(&param);

  Lex_p.Del=Del;
  Lex_p.wi0=wi0;
  Lex_p.y_ex=y_ex;
  Lex_p.y_sim=y_sim;

  N = n_snp + 1;

  //  for (i=0;i<M;++i) {
  //  
  //    Lex_p.r=rk[i];
  //
  //    /*Start the L-BFGS optimization.*/
  //    ret = lbfgs(N, wopt, &Lex, evaluate, progress, &(Lex_p), &param);
  //    printf("n=%3d  Ln = %8.3f\n", i+1,Lex);
  //  }

  /*
  outfile=efopen(outfilename,"w");
  for (i=0;i<n_snp;++i) {
    fprintf(outfile,"%4d %10.8lf\n",i+1,wopt[i]);
  }
  fclose(outfile);
  */

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[-Del <Delta Value> ] \n");
  printf("[-r1  <rho1> ] \n");
  printf("[-r2  <rho2> ] \n");
  printf("[-r3  <rho3> ] \n");
  printf("[-r4  <rhor> ] \n");
  printf("%s inputname1(exp data) inputname2(simu data) outfilename\n",progname);
}

