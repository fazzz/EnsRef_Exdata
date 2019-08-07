#ifndef INCLUDE_ef
#define INCLUDE_ef

FILE *efopen(char *filename,char *flag);
void *emalloc(size_t n);
void *ecalloc(size_t n,size_t size);
//void *egcmalloc(size_t n);
void *gcemalloc(size_t n);
void *gcemalloc_atomic(size_t n);
void *gcerealloc(void *ptr, size_t n);
char *egetenv(const char *name);

#endif
