#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main( int argc, char ** argv )
{
  char fil1[80];
  char fil2[80];
  FILE * pf1, * pf2;
  int n1, n2, il;
  double x1, x2;
  int ncell;

  if( argc != 4 ) {
    fprintf( stderr, "usage: %s output1 output2 ncell\n",argv[0]);
    exit(0);
  }

  pf1 = fopen( argv[1], "r" );
  pf2 = fopen( argv[2], "r" );
  ncell = atoi( argv[3] );
  

  n1 = fscanf(pf1,"%lf",&x1);
  n2 = fscanf(pf2,"%lf",&x2);

  il = 0;

  while( n1 == 1 && n2 == 1 ) {
    if( x1 != x2 ) {
      if( x1 == 0.0 && x2 == 0.0 ) {
        printf("%d %e %e %e\n", il, x1, x2, fabs(x2-x1) );
      } else {
        printf("%d %e %e %e\n", il, x1, x2, fabs(x2-x1)/(fabs(x2)+fabs(x1)) );
      }
    }
    n1 = fscanf(pf1,"%lf",&x1);
    n2 = fscanf(pf2,"%lf",&x2);
    if( ( il % ncell ) == 0 ) printf("=============== ARRAY %d ================\n", il / ncell + 1 );
    il ++;
  }

  return 1;
}
