#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>


int plot_indexes( int ni, int nj, int nk, int max)
{
	int i,j,k,t,l;
	int ii,jj,kk;
	int ord[]={ 0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 };
	char si,sj,sk;
	char di[2],dj[2],dk[2];
	char sdi[10],sdj[10],sdk[10];
	char str[]="  INTEGER, PARAMETER ::";
	t=0;
	for(kk=0;kk<=(nk+nk+1);kk++) {
		k = ord[kk];
		sk = ( k<0 ? 'm' : 'p' );
                for( l = 0; l < abs( k ); l ++ ) sdk[l] = sk;
                sdk[l] = '\0';

		for(jj=0;jj<(nj+nj+1);jj++) {

			j = ord[jj];
			sj = ( j<0 ? 'm' : 'p' );

            	        for( l = 0; l < abs( j ); l ++ ) sdj[l] = sj;
               	        sdj[l] = '\0';

                        switch( j ) {
				case -2:
					strcpy(dj,"bb");
					break;
				case -1:
					strcpy(dj,"b");
					break;
				case +1:
					strcpy(dj,"t");
					break;
				case +2:
					strcpy(dj,"tt");
					break;
				default:
					strcpy(dj,"");
			}

			for(ii=0;ii<(ni+ni+1);ii++) {

				i = ord[ii];
				si = ( i<0 ? 'm' : 'p' );

            	     	   	for( l = 0; l < abs( i ); l ++ ) sdi[l] = si;
               	        	sdi[l] = '\0';

                        	switch( i ) {
					case -2:
						strcpy(di,"ll");
						break;
					case -1:
						strcpy(di,"l");
						break;
					case +1:
						strcpy(di,"r");
						break;
					case +2:
						strcpy(di,"rr");
						break;
					default:
						strcpy(di,"");
				}

				if( ( abs(i)+abs(j)+abs(k) ) <= max ) {
					t++;
					printf("%s i%c%d_j%c%d_k%c%d_ = %3d  ! i%sj%sk%s ! %s%s_\n",
                                          str,si,abs(i),sj,abs(j),sk,abs(k),t,sdi,sdj,sdk,dj,di);
				}
			}
		}
	}
}

int plot_map( int ni, int nj, int nk, int max)
{
	int i,j,k,t;
	int ii,jj,kk;
	int ord[]={ 0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 };
	char si,sj,sk;
	char str[]="   ";
	t=0;
	for(kk=0;kk<=(nk+nk+1);kk++) {
		k = ord[kk];
		sk = ( k<0 ? 'm' : 'p' );
		for(jj=0;jj<(nj+nj+1);jj++) {
			j = ord[jj];
			sj = ( j<0 ? 'm' : 'p' );
			for(ii=0;ii<(ni+ni+1);ii++) {
				i = ord[ii];
				si = ( i<0 ? 'm' : 'p' );
				if( ( abs(i)+abs(j)+abs(k) ) <= max ) {
					t++;
					printf("%s indijk( %2d, %2d, %2d ) = i%c%d_j%c%d_k%c%d_\n",
							str,i,j,k, si,abs(i),sj,abs(j),sk,abs(k));
				}
			}
		}
	}
}



int filter()
{
	char str[256];
	char myarg[256];
	char myarg_new[256];
	char str_new[256];
	char trdarg[256], iaval[256], javal[256];
	char si,sj,sk;
	int i, j, ln, ival, jval, kval, ok, good;
	char *so, *sc, *vi, *vj, *scol, *sp;
	kval = 0;
	gets( str );
	while( !feof(stdin) ) {
		ln = strlen( str );
		for( i = 0; i< (ln-5); i++ ) {
			if( ! strncmp("myijk",&str[i],5) ) {
				so = index(&str[i],'(');
				sc = index(&str[i],')');
				if( so && sc ) {
					for( sp = so, j = 0; sp <= sc; sp++ ) { myarg[j] = *sp; j++;}
					myarg[j] = '\0';
					scol = index(myarg,':');
					vi = index(myarg,',');
					vj = index(vi+1,',');
					ok = 1;
					good = 1;
					if( vi && vj ) {
						for( sp = myarg+1; sp < vi; sp++ ) { if( isalpha( *sp ) ) ok = 0; }
						for( sp = vi+1   ; sp < vj; sp++ ) { if( isalpha( *sp ) ) ok = 0; }
					} else {
						good = 0;
					}
					if( !scol && good && ok ) {
						sscanf(myarg,"(%d,%d,",&ival,&jval);
						si = ( ival<0 ? 'm' : 'p' );
						sj = ( jval<0 ? 'm' : 'p' );
						sk = ( kval<0 ? 'm' : 'p' );
						sprintf( myarg_new, "(i%c%d_j%c%d_k%c%d_%s\0",
								si,abs(ival),sj,abs(jval),sk,abs(kval),vj);
						printf("! CHANGE  %s => %s \n",myarg,myarg_new);
					} else if( !scol && good && !ok ) {
						for( sp = myarg+1, j = 0; sp < vi; sp++ ) { iaval[j] = *sp; j++;}
						iaval[j] = '\0';
						for( sp = vi+1   , j = 0; sp < vj; sp++ ) { javal[j] = *sp; j++;}
						javal[j] = '\0';
						sprintf( myarg_new, "(indijk(%s,%s,0)%s\0",iaval,javal,vj);
						printf("! CHANGE  %s => %s \n",myarg,myarg_new);
					}
				}
			}
		}
		puts( str );
		gets( str );
	} while( !feof(stdin) );
}


int cmp_inds( char * out, char * in ) {
  out[0]='\0';
  if( ! strcmp( in, "r_" ) ) strcpy( out, "ip1_jp0_kp0_" );
  if( ! strcmp( in, "l_" ) ) strcpy( out, "im1_jp0_kp0_" );
  if( ! strcmp( in, "rr_" ) ) strcpy( out, "ip2_jp0_kp0_" );
  if( ! strcmp( in, "ll_" ) ) strcpy( out, "im2_jp0_kp0_" );
  if( ! strcmp( in, "t_" ) )  strcpy( out, "ip0_jp1_kp0_" );
  if( ! strcmp( in, "tr_" ) )  strcpy( out, "ip1_jp1_kp0_" );
  if( ! strcmp( in, "tl_" ) )  strcpy( out, "im1_jp1_kp0_" );
  if( ! strcmp( in, "b_" ) )  strcpy( out, "ip0_jm1_kp0_" );
  if( ! strcmp( in, "br_" ) )  strcpy( out, "ip1_jm1_kp0_" );
  if( ! strcmp( in, "bl_" ) )  strcpy( out, "im1_jm1_kp0_" );
  if( ! strcmp( in, "tt_" ) )  strcpy( out, "ip0_jp2_kp0_" );
  if( ! strcmp( in, "bb_" ) )  strcpy( out, "ip0_jm2_kp0_" );
  
}


int filter2()
{
	char str[256];
	char myarg[256];
	char myarg_new[256];
	char str_new[256];
	char trdarg[256], in[256], out[256];
	char si,sj,sk;
	int i, j, ln, ival, jval, kval, ok, good;
	char *so, *sc, *vi, *vj, *scol, *sp;
	kval = 0;
	gets( str );
	while( !feof(stdin) ) {
		ln = strlen( str );
		for( i = 0; i< (ln-5); i++ ) {
			if( ! strncmp("myinds",&str[i],5) ) {
				so = index(&str[i],'(');
				sc = index(&str[i],')');
				if( so && sc ) {
					for( sp = so, j = 0; sp <= sc; sp++ ) { myarg[j] = *sp; j++;}
					myarg[j] = '\0';
					vi = index(myarg,',');
					if( vi ) {
						for( sp = myarg+1, j = 0; sp < vi; sp++ ) { in[j] = *sp; j++;}
						in[j] = '\0';
						good = 1;
					} else {
						good = 0;
					}
					ok = 1;
					if( ( !index( in, 'b' ) ) && 
					    ( !index( in, 'r' ) ) && 
					    ( !index( in, '_' ) ) ) ok = 0 ;
					if( good && ok ) {
						cmp_inds( out, in );	
						printf("! CHANGE  %s => %s \n",in,out);
					}
				}
			}
		}
		puts( str );
		gets( str );
	} while( !feof(stdin) );
}



int main()
{
	/* plot_indexes( 2, 2, 2, 2 ); */
	filter2(); 
	/* plot_map( 2, 2, 2, 2 ); */

}
