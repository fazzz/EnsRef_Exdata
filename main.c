
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#define ON 0
#define OFF 1

double k_B=1.98723e-3;
double T=300;

int USAGE(char *progname);

double C_func(double k, double x0, double x, int periodicflag, double periodicity);

int main(int argc, char *argv[]) {
  int i,ii,j,k,l,d,num,numnum,interval=1;
  double dummy,f;

  int periodicflag=OFF;
  double periodicity,minv,maxv;
  
  int num_sim,*n,n_sum;     // # of simulation, # of snapshots
  double **x_ij;            // raw data

  double *k_umbrella,*x0; // parameters for umbrella sampling

  double ***C;

  double *f_i,**P_ij;

  double sum=0.0;

  double lnp,L,Lnew,dL=0.0,ddL=1.0,threhold=1.0e-4; // threhold for conversion etc.

  char *inputfilelistname,*metadatafilename,*outputfilename,*pmffilename,*pmftempfilename="PMFOUT";
  FILE *inputfile,*inputfilelist,*metadatafile,*outputfile,*pmffile,*pmftempfile;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"interval",1,NULL,'l'},
    {"tol",1,NULL,'t'},
    {"pi",0,NULL,'i'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hil:t:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'l':
      interval=atoi(optarg);  break;
    case 't':
      threhold=atof(optarg);  break;
    case 'i':
      periodicflag=ON; periodicity=2.0*pi;
      minv=-1.0*pi; maxv=1.0*pi;
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  num_sim           = atoi(*argv);
  inputfilelistname = *++argv;
  metadatafilename  = *++argv;
  pmffilename       = *++argv;

  x_ij=(double **)gcemalloc(sizeof(double *)*(num_sim));
  for (i=0;i<num_sim;++i) x_ij[i]=(double *)gcemalloc(sizeof(double));
  n=(int *)gcemalloc(sizeof(int)*(num_sim));

  n_sum=0;
  inputfilelist=efopen(inputfilelistname,"r");
  for (i=0;i<num_sim;++i) {
    getline(&line,&len,inputfilelist);
    l=strlen(line);
    line[l-1]='\0';

    inputfile=efopen(line,"r");
    num=0; numnum=0;
    d = 1;
    while ( d != -1  )  {
      fscanf(inputfile,"%d",&dummy);
      d=fscanf(inputfile,"%lf",&f);
      if (numnum%interval==0) {
	x_ij[i]=(double *)gcerealloc(x_ij[i],sizeof(double)*(num+1));
	if (periodicflag==ON) {
	  while (f<minv) {
	    f+=periodicity;
	  }
	  while (f>maxv) {
	    f-=periodicity;
	  }
	}
	x_ij[i][num]=f;
	++num;
      }
      ++numnum;
    }
    fclose(inputfile);
    n[i]=num-1;
    n_sum+=n[i];
  }
  fclose(inputfilelist);

  k_umbrella=(double *)gcemalloc(sizeof(double)*num_sim);
  x0=(double *)gcemalloc(sizeof(double)*num_sim);
  metadatafile=efopen(metadatafilename,"r");
  for (i=0;i</*K*/num_sim;++i) {
    fscanf(metadatafile,"%lf",&x0[i]); fscanf(metadatafile,"%lf",&k_umbrella[i]); 
    if (periodicflag==ON) {
      while (x0[i]<minv) {
	x0[i]+=periodicity;
      }
      while (x0[i]>maxv) {
	x0[i]-=periodicity;
      }
    }
  }
  fclose(metadatafile);
  
  P_ij=(double **)gcemalloc(sizeof(double)*num_sim);
  for (i=0;i<num_sim;++i) P_ij[i]=(double *)gcemalloc(sizeof(double)*n[i]);

  f_i=(double *)gcemalloc(sizeof(double)*num_sim);

  for (i=0;i<num_sim;++i) f_i[i]=1.0;
  for (i=0;i<num_sim;++i) for (j=0;j<n[i];++j) P_ij[i][j]=1.0/n_sum;

  C=(double ***)gcemalloc(sizeof(double **)*num_sim);
  for (i=0;i<num_sim;++i) C[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) for (j=0;j<num_sim;++j) C[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);

  for (i=0;i<num_sim;++i) 
    for (j=0;j<num_sim;++j) 
      for (k=0;k<n[j];++k)
	C[i][j][k]=C_func(k_umbrella[i],x0[i],x_ij[j][k],periodicflag,periodicity);

  ii=0;
  dL=threhold;
  while (ddL >= threhold || ii==0) {

    for (i=0;i<num_sim;++i) f_i[i]=0.0;
    for (i=0;i<num_sim;++i) for (j=0;j<num_sim;++j) for (k=0;k<n[j];++k) f_i[i]+=C[i][j][k]*P_ij[j][k];
    for (i=0;i<num_sim;++i) f_i[i]=1.0/f_i[i];

    for (i=0;i<num_sim;++i) for (j=0;j<n[i];++j) P_ij[i][j]=0.0;
    for (i=0;i<num_sim;++i) 
      for (j=0;j<n[i];++j) 
	for (k=0;k<num_sim;++k) 
	  P_ij[i][j]+=n[k]*f_i[k]*C[k][i][j];

    for (i=0;i<num_sim;++i) 
      for (j=0;j<n[i];++j) 
	P_ij[i][j]=1.0/P_ij[i][j];

    sum=0.0; for (i=0;i<num_sim;++i) for (j=0;j<n[i];++j) sum+=P_ij[i][j];
    for (i=0;i<num_sim;++i) for (j=0;j<n[i];++j) P_ij[i][j]=P_ij[i][j]/sum;

    Lnew=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	lnp=0.0;
	for (k=0;k<num_sim;++k)
	  for (l=0;l<n[k];++l)
	    lnp+=f_i[k]*C[k][k][k]*P_ij[k][l];
	lnp=log(lnp);
	Lnew+=lnp;
      }
    }

    if (ii>0) {
      dL=fabs(Lnew-L);
      ddL=dL/fabs(L)*100.0;
    }
    L=Lnew;

    ++ii;
    printf("%3d %5.3e\n",ii,L);

  }

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      fprintf(pmffile,"%10.8lf %10.8lf %10.8lf\n",x_ij[i][j],P_ij[i][j],-1.0*log(P_ij[i][j]));
    }
  }
  fclose(pmffile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}

double C_func(double k, double x0, double x, int periodicflag, double periodicity){
  double y,delta;

  delta=fabs(x-x0);
  if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=periodicity-delta;

  y=exp(-0.5*k*delta*delta/k_B/T);

  return y;
}
