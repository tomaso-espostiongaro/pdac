#ifdef __SP4
#include<fp.h>

int float_chk( double *a ) {
   int nan;
   nan = IS_NAN( *a );
   return nan;
}

#else

int float_chk( double *a ) {
   int nan = 0;
   return nan;
}


#endif
