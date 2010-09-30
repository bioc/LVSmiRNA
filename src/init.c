
/* Inizialize Lapack routines */
/* See src/main/lapack.c in R source */


#include <R_ext/Rdynload.h>

static int Lapack_initialized = 0;
static int use_lapack = 0;

#ifdef Win32
# include <fcntl.h>
#endif

void Lapack_Initialize(void)
{
    /* Initializing LAPACK */
    use_lapack = 1;
    Lapack_initialized = 1;
    return;
}


/* void Lapack_Initialize(void) */
/* { */
/*   int res = R_moduleCdynload("lapack", 1, 1); */
/*   Lapack_initialized = -1; */
/*   if(!res) return; */
/*   Lapack_initialized = 1; */
/* #ifdef Win32 */
/*   /\* gfortran initialization sets these to _O_BINARY *\/ */
/*   setmode(1,_O_TEXT); /\* stdout *\/ */
/*   setmode(2,_O_TEXT); /\* stderr *\/ */
/* #endif */
/*   return; */
/* } */
