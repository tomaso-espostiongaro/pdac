/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdio.h>
#include<time.h>
#include<ctype.h>
#include<sys/types.h>
#include<sys/time.h>

#if defined __BEOWULF
#define ELAPSED_SECONDS elapsed_seconds_
#define CCLOCK cclock_wall_
#define CCLOCK2 cclock_cpu_
#else
#define ELAPSED_SECONDS elapsed_seconds
#define CCLOCK cclock_wall
#define CCLOCK2 cclock_cpu
#endif

struct tms {
        clock_t tms_utime;              /* user time */
        clock_t tms_stime;              /* system time */
        clock_t tms_cutime;             /* user time, children */
        clock_t tms_cstime;             /* system time, children */
};


double ELAPSED_SECONDS()
{
  static time_t tstart, tend;
  static int first = 1;
  double sec;
  time(&tend);
  if( first ) {
    tstart = tend;
    first = 0;
  }
  sec = difftime( tend, tstart );
  return sec;
}


double CCLOCK()
/* Restituisce i secondi (walltime) trascorsi dalla chiamata al timer rest */
{

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;

}



double CCLOCK2( void )
/* Restituisce i secondi (user time) trascorsi dalla chiamata al timer rest */
{

    static struct tms tmp;
    double s;
    times( &tmp );
    s = (double) tmp.tms_utime ;
    s = s / (double) CLK_TCK ;
    return s;

}


#ifdef __PGI
/* This wrapper is used with PGI compilers that do not have the
   F95 intrinsic subroutine "cpu_time"
*/
  int cpu_time_( double * sec ) {
    *sec = CCLOCK2();
    return 0;
  }
#endif
