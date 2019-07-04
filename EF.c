#include <stdio.h>
#include <stdlib.h>

#include "EF.h"
#include "gc.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

FILE *efopen(char *filename,char *flag){
  FILE *file;

  if((file=fopen(filename,flag))==NULL) {
    printf("open error about %s\n",filename);
    exit(1);
  }
  return file;
}

void *emalloc(size_t n){
  void *p;

  if ((p=malloc(n))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

void *gcemalloc(size_t n){
  void *p;

  if ((p=GC_MALLOC(n))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

void *ecalloc(size_t t,size_t size){
  void *p;

  if ((p=calloc(t,size))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

/******************************************/
/* void *gcecalloc(size_t t,size_t size){ */
/*   void *p;				  */
/* 					  */
/*   if ((p=GC_CALLOC(t,size))==NULL) {	  */
/*     printf("open about allocation\n"); */
/*     exit(1);				  */
/*   }					  */
/* 					  */
/*   return p;				  */
/* }					  */
/******************************************/

void *gcerealloc(void *ptr, size_t n){
  void *p;

  if ((p=GC_REALLOC(ptr, n))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

char *egetenv(const char *name) {
  char *env;

  if ((env=getenv(name))==NULL) {
    printf("error; there is not eiviromental variable %s\n",name);
    exit(1);
  }
  else {
    return env;
  }
}

