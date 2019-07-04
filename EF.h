#ifndef INCLUDE_ef
#define INCLUDE_ef

#include <netcdf.h>

FILE *efopen(char *filename,char *flag);
void *emalloc(size_t n);
void *ecalloc(size_t n,size_t size);
//void *egcmalloc(size_t n);
void *gcemalloc(size_t n);
void *gcerealloc(void *ptr, size_t n);
char *egetenv(const char *name);
int enc_create(const char *path, int cmode, int *ncidp);
int enc_def_dim(int ncid, const char *name, size_t len, int *idp);
int enc_def_var(int ncid, const char *name, nc_type xtype, int ndims,
		const int *dimidsp, int *varidp);
int enc_inq_varid(int ncid, const char *name, int *varidp);
int enc_open(const char *path, int mode, int *ncidp);
int encopen(const char* path, int mode);
int encclose(int ncid);

#endif
