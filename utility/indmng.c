#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>


int plot_indexes( int ni, int nj, int nk, int max)
{
	int i,j,k,t;
	int ii,jj,kk;
	int ord[]={ 0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 };
	char si,sj,sk;
	char str[]="  INTEGER, PARAMETER ::";
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
					printf("%s i%c%d_j%c%d_k%c%d_ = %3d\n",str,si,abs(i),sj,abs(j),sk,abs(k),t);
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

int main()
{
	/* plot_indexes( 2, 2, 2, 2 ); */
	filter();
	/* plot_map( 2, 2, 2, 2 ); */

}
