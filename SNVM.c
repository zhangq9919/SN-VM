
//////////////SN Virtual Machine: low-level sleptsov net virtual machine /////////
// run low-level sleptsov net                                                   //
// On single core and multicore parallel architectures using OpenMP             //
// Lsn file format:                                                             //
// >place transition arc l NST                                                  //
// >arcs                                                                        //
// >initial marking: p mu                                                       //
// Command line:                                                                //
// SNVM lsn_file.lsn output_file.txt                                            //
// -fopenmp -O3                                                                 //
//@ 2023 Qing Zhang: zhangq9919@163.com, Dmitry Zaitsev                         //
//////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#define MAXSTRLEN 16384
#define MATRIX_SIZE(d1,d2,t) ((d1)*(d2)*(sizeof(t)))
#define MOFF(i,j,d1,d2) ((d2)*(i)+(j))
#define MELT(x,i,j,d1,d2) (*((x)+MOFF(i,j,d1,d2)))
#define SKIP_COMM  do{fgets(input_buffer,MAXSTRLEN,fi);} while(input_buffer[0]==';');//skip comments  

/*Declare error code*/
#define NO_ERR 0 /* No error */
#define ERR_NEGATIVE 1 /* Wrong header(m n k l): negative numbers */
#define ERR_LESS 2 /* Wrong header(m n k l): The number is less than 5 or non-numeric*/
#define ERR_HEADER_ZERO 3 /* Wrong header(m n k l): Value zero */
#define ERR_ARC_ZERO 4 /* Wrong arc: Value zero */
#define ERR_FORMAT 5 /* Wrong arc: Arc format error */
#define ERR_INHIBITOR 6 /* Wrong arc: The weight of inhibitior arc is not -1 */
#define ERR_EXCEEDM 7 /* Exceed number: Exceed the number of m */
#define ERR_EXCEEDN 8 /* Exceed number: Exceed the number of n */
#define ERR_DUP 9 /* Duplicate arc: The arc is duplicate */
#define ERR_OVERFLOW 10 /* Numeric overflow: The value is out of range */
#define ERR_NST 11 /* Wrong header: NST is nonzero */


static char Help[] =
  "LSN Virtual Machine\n"
  "action: run low-level sleptsov net\n"
  "        on multicore parallel architectures using OpenMP.\n"
  "usage:          lsn_pro_level [-h]\n"
  "                              [-d][-nth][-pm][-smax][-rm][-r][-lsntonet]\n"
  "                              [lsn_file.lsn][output_file.txt]\n"
  "FLAGS               WHAT                                                                       DEFAULT\n"
  "-h                  this mode\n"
  "-d n                Output data level (n = 0 or 1 or 2)                                        -d 0\n"
  "-nth n              number of thread                                                           -nth 8\n"
  "-pm n               choose format for printing marking (n = 0 - usual or 1 - sparse vector)    -pm 1\n"
  "-smax n             specified steps (n = steps)\n"
  "-rm                 print raw matices and vectors\n"
  "-r                  Sleptsov net rule\n"
  "-lsntonet           save lsn as .net\n"
  "lsn_file.lsn        Low-level sleptsov net file based on LSN format (stdin if absent)           stdin\n"
  "output_file.txt     The result of the LSN running through SN Virtual Machine (stdout if absent) stdout\n";

char* errmsg[] = {        /* 0 */ "No Error",
                                  /* 1 */ "Wrong header: Negative numbers!",
                                  /* 2 */ "Wrong header: The number is less than 5 or non-numeric!",
                                  /* 3 */ "Wrong header: Value zero!",
                                  /* 4 */ "Wrong arc: The arc has value zero!",
                                  /* 5 */ "Wrong arc: Arc format error!",
                                  /* 6 */ "Wrong arc: The weight of inhibitior arc is not -1!",
                                  /* 7 */ "Exceed number: Exceed the number of m!",
                                  /* 8 */ "Exceed number: Exceed the number of n!",
                                  /* 9 */	"Duplicate arc: The arc is duplicate!",
                                  /* 10 */ "Numeric overflow: The value is out of range!",
                                  /* 11 */ "Wrong header: NST is nonzero!",
                 };

// declare global variables
long err;//error code

double magma_wtime( void )
{
	struct timeval t;
	gettimeofday( &t, NULL );
	return t.tv_sec + t.tv_usec*1e-6;
}

/**
Function: Output error information according to the error code.
**/
void perr(char* info);


/**
Function: Generate a matrix of priority arc chains.
**/
void priority_chain(int *R,int n,int nth);


/**
Function: Parse command line.
**/
int command(int numf,int snumf,int argc,int *nth,int *debug_level,int *o,int *printm,long int *smax,int *rm,int *lsntonet,char *argv[]);


/**
Function: Allocate space for matrices and vector.
 **/
int allocate(int *m,int *n,int *k,int *l,int *NST,int **mu,int **f,int **B,int **D,int **R,FILE *fi,char input_buffer[]);


/**
Function: Read LSN file, including header, arcs, initial marking.
**/
void read_lsn(int m,int n,int k,int l,int NST,int v1,int v2,int v3,int digit,int *mu,int *f,int *B,int *D,int *R,int debug_level,int nth,int printm,int rm,int lsntonet,FILE *fi,FILE *fo,char input_buffer[]);


/**
Matrix and vector initialization
**/
void zeroBDRmu(int m,int n,  int *mu,int *f,int *B,int *D,int *R);

/**
Function: Read the arcs information according to the LSN file format, and generate the priority arc matrix R,
a pair of matrices B and D are matrices of incoming and outgoing arcs of transitions.

Case 1: If v1 < 0, v2 < 0, arc from transition v1 to transition v2. It is a priority arc with w = 0.
Case 2: If v1 > 0, v2 > 0, arc from place v1 to transition v2 of multiplicity w.
Case 3: If v1 < 0, v2 > 0, arc from transition v2 to place |v1| of multiplicity w.
Case 4: When w is negative, it is inhibitor arc.
**/
void read_arc( int v1, int v2, int v3,int m,int n,int *B,int *D, int *R,int debug_level,int rm,FILE *fo);


/**
Function: Read the number of places, the number of transitions, the number of arcs, nonzero number of marking, number of substitution transition (0 for LSN)
**/
void read_mnkl(int m,int n,int k,int l,int NST,int digit,int rm,int debug_level,FILE *fo);


/**
SN Virtual Machine: Run LSN according to the SN firing rule
Include:
Algorithm1: the multiplicity of fireable transitions
Algorithm2: The Marking Rule
Algorithm3: Remove Low Priority Transition
**/
void run_lsn(int *f, int *mu, int *B,int *D,int *R,int m, int n,int debug_level, int nth,int printm,FILE *fo,long int smax);

void printBDR(int m,int n,int *B,int *D,int *R,int debug_level,int rm,int *mu,int lsntonet,FILE *fo);

void printmk(int m,  int *mu,int printm,FILE *fo);

int main(int argc,char *argv[])
{
	int m, // number of places
	    n, // number of transitions
	    k, // number of arcs
	    l,//number of nonzero markings
	    NST;//number of substitution transition (0 - LSN)
	int i,j,z;
	int v1,v2,v3;//p t w
	int *B = NULL,*D = NULL;//input matrix & output matrix
	int *R = NULL;//Priority arc matrix
	int *f = NULL;// indicator of fireable transitions 0 or any numbers)
	int *mu; // marking (should be inputted from file)
	int p;//place
	int digit = 0;
	int debug_level=0;//-d
	int nth=omp_get_max_threads();//-nth
	int printm=1;//choose format for printing marking (0-usual or 1-sparse vector) -pm
	long int smax;//-smax
	int rm = 1;//-rm
	int lsntonet = 1;//-lsntonet
	int ii = 0;//pointer
	int numf,snumf;
	double t1,t2;
	FILE *fi, *fo;//fi - lsn_file.lsn, fo - output_file.txt
	int o = 1;
	char input_buffer[MAXSTRLEN + 1];

	/* parse command line */
	command(numf,snumf,argc,&nth,&debug_level,&o,&printm,&smax,&rm,&lsntonet,argv);
	/*lsn_ile.lsn & output_file.txt*/
	if( argv[ o ]==NULL ) fi = stdin;
	else fi = fopen( argv[ o ], "r" );
	if( fi == NULL ) {
		fprintf( stderr, "*** error open file %s\n",  argv[ o ] );
		exit(2);
	}
	if( argv[ o+1 ]==NULL ) fo = stdout;
	else fo = fopen( argv[ o+1 ], "w" );
	if( fo == NULL ) {
		fprintf( stderr, "*** error open file %s\n", argv[ o+1 ] );
		exit(2);
	}
	//allocate space for matrices and vector
	digit = allocate(&m,&n,&k,&l,&NST,&mu,&f,&B,&D,&R,fi,input_buffer);
	t1=magma_wtime();
	/*read LSN file*/
	read_lsn(m,n,k,l,NST,v1,v2,v3,digit,mu,f,B,D,R,debug_level,nth,printm,rm,lsntonet,fi,fo,input_buffer);
	/*print matrices B D R*/
	printBDR(m,n,B,D,R,debug_level,rm,mu,lsntonet,fo);
	/*SN Virtual Machine: run LSN according to firing rule of SN*/
	run_lsn(f,mu,B,D,R,m,n,debug_level,nth,printm,fo,smax);
	/*print final marking*/
	printmk(m,mu,printm,fo);
	t2=magma_wtime();
	printf("####LSN information has been created!####\n ");
	printf("####The total time is %#fs\n", t2-t1 );

	free(B);
	free(D);
	free(R);
	free(mu);
	free(f);
	return(0);
} /*main*/


//Output error information according to the error code.
void perr(char* info)
{
	if( info ) {
		fprintf(stderr,"Error %s:  %s\n",info,errmsg[err]);
		return;
	}
	fprintf(stderr,"Error %s:  %s\n",info,errmsg[err]);
}

//Generate a matrix of priority arc chains.
void priority_chain(int *R,int n,int nth)
{
	int x,y,z;
	#pragma omp parallel for private(x,y,z) num_threads(nth)
	for ( x = 0; x < n; x++) {
		for ( y = 0; y < n; y++) {
			if (MELT(R,x,y,n,n) != 0) {
				for ( z = 0; z < n; z++) {
					if (MELT(R,y,z,n,n) != 0) {
						MELT(R,x,z,n,n) = 1;
					}
				}
			}
		}
	}
}

//Parse command line.
int command(int numf,int snumf,int argc,int *nth,int *debug_level,int *o,int *printm,long int *smax,int *rm,int *lsntonet,char *argv[])
{
	int i;
	numf=0;
	for( i=1; i<argc; i++ ) {
		if( strcmp( argv[i], "-h" )==0 ) {
			printf( "%s", Help );
			*o+=1;
		}
		if( strcmp( argv[i], "-d" )==0 ) {
			snumf=numf;
			numf=1;
		} else if( numf==1 ) {
			*debug_level = atoi( argv[ i ] );
			numf=snumf;
			*o+=2;
			printf("#####The debug level is %d\n",*debug_level);
		} else if( strcmp( argv[i], "-nth" )==0 ) {
			snumf=numf;
			numf=2;
		} else if( numf==2 ) {
			*nth = atoi( argv[ i ] );
			numf=snumf;
			*o+=2;
			printf("number of threads: %d\n",*nth);
		} else if( strcmp( argv[i], "-pm" )==0 ) {
			snumf=numf;
			numf=3;
		} else if( numf==3 ) {
			*printm = atoi( argv[ i ] );
			numf=snumf;
			*o+=2;
			if(*printm == 0) printf("format for printing marking: %d - usual marking\n",*printm);
			if(*printm == 1) printf("format for printing marking: %d - sparse vector\n",*printm);
		} else if( strcmp( argv[i], "-smax" )==0 ) {
			snumf=numf;
			numf=4;
		} else if( numf==4 ) {
			*smax = atoi( argv[ i ] );
			numf=snumf;
			*o+=2;
		} else if( strcmp( argv[i], "-rm" )==0 ) {
			*rm = 0;
			*o+=1;
		} else if( strcmp( argv[i], "-r" )==0 ) {
			*o+=1;
			printf("#####The processor adopts the sleptsov net rule.\n");
		} else if( strcmp( argv[i], "-lsntonet" )==0 ) {
			*lsntonet = 0;
			*o+=1;
		}
		return( 0);
	}
}

//allocate space for matrices and vector
int allocate(int *m,int *n,int *k,int *l,int *NST,int **mu,int **f,int **B,int **D,int **R,FILE *fi,char input_buffer[])
{
	SKIP_COMM
	int digit = sscanf(input_buffer,"%d %d %d %d %d",m,n,k,l,NST);

	//check the number of n and m
	if(*m < 0 || *n < 0) {
		err = ERR_NEGATIVE;
		perr("1");
		exit(1);
	}
	if(*m == 0 || *n == 0) {
		err = ERR_HEADER_ZERO;
		perr("3");
		exit(1);
	}

	*B = (int *)malloc( MATRIX_SIZE(*m,*n,int) );
	*D = (int *)malloc( MATRIX_SIZE(*m,*n,int) );
	*R = (int *)malloc( MATRIX_SIZE(*n,*n,int) );
	*mu = (int *)malloc((*m)*sizeof(int));
	*f = (int *)malloc((*n)*sizeof(int));
	if(*B==NULL) {
		printf("***error: no enough memory for B\n");
		exit(2);
	}
	if(*D==NULL) {
		printf("***error: no enough memory for D\n");
		exit(2);
	}
	if(*R==NULL) {
		printf("***error: no enough memory for R\n");
		exit(2);
	}

	if(*mu==NULL) {
		printf("***error: no enough memory for mu\n");
		exit(2);
	}

	if(*f==NULL) {
		printf("***error: no enough memory for f\n");
		exit(2);
	}
	return digit;
}

/* Matrix and vector initialization*/
void zeroBDRmu(int m,int n,int *mu,int *f,int *B,int *D,int *R)
{
	memset(B,0,MATRIX_SIZE(m,n,int));
	memset(D,0,MATRIX_SIZE(m,n,int));
	memset(R,0,MATRIX_SIZE(n,n,int));
	memset(mu,0,m*sizeof(int));
	memset(f,0,n*sizeof(int));
}

/* Read LSN file, including header, arcs, initial marking.*/
void read_lsn(int m,int n,int k,int l,int NST,int v1,int v2,int v3,int digit,int *mu,int *f,int *B,int *D,int *R,int debug_level,int nth,int printm,int rm,int lsntonet,FILE *fi,FILE *fo,char input_buffer[])
{
	int i,j,p;
	/*zeroBDRmu*/
	zeroBDRmu(m,n,mu,f,B,D,R);
	/*read the header: m n k l NST*/
	read_mnkl(m,n,k,l,NST,digit,rm,debug_level,fo);
	/*read arcs*/
	for(i=0; i<k; i++) {
		SKIP_COMM
		sscanf(input_buffer,"%d %d %d",&v1,&v2,&v3);// p t w
		read_arc(v1,v2,v3,m,n,B,D,R,debug_level,rm,fo);
	}
	/*read initial marking*/
	for(i=0; i<l; i++) {
		SKIP_COMM
		sscanf(input_buffer,"%d %d\n",&v1,&v2);//p mu
		mu[v1-1] = v2;
		if(debug_level>0&&rm) fprintf(fo,"%d %d\n",v1,v2);
	}
	if(fi != stdin ) fclose( fi );

	/* save .lsn to .net*/
	if(!lsntonet) {
		for(j=1; j<=n; j++) {
			fprintf(fo,"tr t%d ",j);
			for(i=1; i<=m; i++) {
				if(MELT(B,i-1,j-1,m,n)==1)
					fprintf(fo,"p%d ",i);
				else if(MELT(B,i-1,j-1,m,n)==-1)
					fprintf(fo,"p%d\?-1 ",i);
				else if(MELT(B,i-1,j-1,m,n)>1)
					fprintf(fo,"p%d*%d ",i,MELT(B,i-1,j-1,m,n));
			}
			fprintf(fo,"-> ");
			for(i=1; i<=m; i++) {
				if(MELT(D,i-1,j-1,m,n)==1)
					fprintf(fo,"p%d ",i);
				else if(MELT(D,i-1,j-1,m,n)>1)
					fprintf(fo,"p%d*%d ",i,MELT(D,i-1,j-1,m,n));
			}
			fprintf(fo,"\n");
		}

		for(i=1; i<=m; i++) {
			if(mu[i-1]>0) {
				fprintf(fo,"pl p%d (%d)\n",i,mu[i-1]);
			}
		}

		for(j=1; j<=n; j++) {
			for(i=1; i<=n; i++) {
				if(MELT(R,j-1,i-1,n,n)==1)
					fprintf(fo,"pr t%d > t%d\n",j,i);
			}
		}

		fprintf(fo,"net bz\n");
		printf("#####Save .lsn to .net!\n");
		exit(0);
	}


	/*call function about priority arc chain*/
	priority_chain(R,n,nth);

	if(!rm) {
		fprintf(fo,"%d %d\n\n",m,n);
		printBDR(m,n,B,D,R,debug_level,rm,mu,lsntonet,fo);
		for(p=1; p<=m; p++) {
			fprintf(fo,"%d ",mu[p-1]);
		}
		fprintf(fo,"\n");
		printf("###The raw m, n, matrices and vectors have been generated!\n");
		exit(0);
	}

}

//read header
void read_mnkl(int m,int n,int k,int l,int NST,int digit,int rm,int debug_level,FILE *fo)
{
	if(m < 0 || n < 0 || k < 0 ) {
		err = ERR_NEGATIVE;
		perr("1");
		exit(1);
	}

	if(NST!=0) {
		err = ERR_NST;
		perr("11");
		exit(1);
	}

	if(digit != 5 ) {
		err = ERR_LESS;
		perr("2");
		exit(1);
	}
	if(m == 0 || n == 0 || k == 0) {
		err = ERR_HEADER_ZERO;
		perr("3");
		exit(1);
	}
	if(debug_level > 0&&rm) fprintf(fo,"m=%d n=%d k=%d l=%d NST=%d\n",m,n,k,l,NST);
}

//read arcs information
void read_arc( int v1, int v2, int v3,int m,int n,int *B,int *D, int *R,int debug_level,int rm,FILE *fo)
{
	int a, b;
	if(debug_level > 0&&rm)  fprintf(fo,"%d %d %d\n",v1,v2,v3);
	if(v1 == 0 || v2 == 0) {
		err = ERR_ARC_ZERO;
		perr("4");
		exit(1);
	}
	if(v1 > m) {
		err = ERR_EXCEEDM;
		perr("7");
		exit(1);
	}
	if(v2 > n) {
		err = ERR_EXCEEDN;
		perr("8");
		exit(1);
	}

	if(v1 > 0) {
		if(v2 < 0 ) {
			fprintf(stderr,"Warning: arc %d %d %d\n",v1,v2,v3);
			err = ERR_FORMAT;
			perr("5");
			exit(1);
		}
		a=v1,b=v2;
		if(MELT(B,a-1,b-1,m,n) != 0) {
			fprintf(stderr,"Warning: arc %d %d\n",a,b);
			err = ERR_DUP;
			perr("9");
			exit(1);
		}
		if(v3>0) {
			MELT(B,a-1,b-1,m,n)=v3; // B[a-1][b-1]=v3
			if(debug_level==2&&rm) fprintf(fo,"arc from place %d to transition %d of myltiplicity %d\n",v1,v2,v3);//v1,v3 are positive
		}
		if(v3<0) {
			if(v3 != -1 ) {
				fprintf(stderr,"Warning: arc %d %d %d\n",v1,v2,v3);
				err = ERR_INHIBITOR;
				perr("6");
				exit(1);
			}
			MELT(B,a-1,b-1,m,n)=v3; //B[a-1][b-1]=v3
			if(debug_level==2&&rm) {
				fprintf(fo,"inhibitor arc from place %d to transition %d\n",v1,v2);        //v3 is negative,it is inhibitor
			}
		}
		if(v3==0) {
			fprintf(stderr,"Warning: arc %d %d %d\n",v1,v2,v3);
			err = ERR_FORMAT;
			perr("5");
			exit(1);
		}
	} else if(v2 < 0) {
		if(v3 != 0) {
			fprintf(stderr,"Warning: arc %d %d %d\n",v1,v2,v3);
			err = ERR_FORMAT;
			perr("5");
			exit(1);
		}
		a = abs(v1),b = abs(v2);
		if( MELT(R,a-1,b-1,n,n) != 0) {
			fprintf(stderr,"Warning: arc %d %d\n",a,b);
			err = ERR_DUP;
			perr("9");
			exit(1);
		}
		MELT(R,a-1,b-1,n,n) = 1;//R[a-1][b-1]=1
		if(debug_level==2&&rm) fprintf(fo,"priority arc from transition %d to transition %d\n",abs(v1),abs(v2));// t --> t

	} else {
		if(v3 <= 0 ) {
			fprintf(stderr,"Warning: arc %d %d %d\n",v1,v2,v3);
			err = ERR_FORMAT;
			perr("5");
			exit(1);
		}
		a = abs(v1),b = v2;
		if( MELT(D,a-1,b-1,m,n) != 0) {
			fprintf(stderr,"Warning: arc %d %d\n",a,b);
			err = ERR_DUP;
			perr("9");
			exit(1);
		}
		MELT(D,a-1,b-1,m,n)=v3; //D[a-1][b-1]=v3
		if(debug_level==2&&rm) fprintf(fo,"arc from transition %d to place %d of myltiplicity %d\n",v2,abs(v1),v3);//v1 is negative(outgoing)
	}
}

/* SN-VM*/
void run_lsn(int *f, int *mu, int *B,int *D,int *R,int m, int n,int debug_level, int nth,int printm,FILE *fo,long int smax)
{
	long int s = 0;// current step number
	int nf;// number of fireable transitions
	int p, t;// place, transition
	int firable = 0;//multiplicity of fireable transition
	//int firable_m = INT_MAX; //the minimum multiplicity of fireable transition
	int rn; //random number
	int firing_n; //firing transition number
	int l; //remainder
	int tb,td; //columns of matrix B and D
	int ct; // count of fireable transitions
	int fm; //multiplicity of firing transition
	int i,j;

	//steps
	while(1) {

		s++;
		if(s<0) fprintf(stderr,"***error: step overflow!");
		if(debug_level==2) fprintf(fo,"--- step %ld ---\n",s);
		//form vector of fireable transitions f (0 - not fireable, numbers - fireable)
		nf = 0;

		//Algorithm1: the multiplicity of fireable transitions
		#pragma omp parallel for private(p,t,firable) num_threads(nth) //collapse(2)
		for(t=1; t<=n; t++) {
			int	firable_m = INT_MAX;//Set the initial value to a larger value
			#pragma omp simd reduction(min:firable_m)
#pragma unroll
			for(p=1; p<=m; p++) {
				if( MELT(B,p-1,t-1,m,n)>0 ) {// regular arc  B[p-1][t-1]>0
					firable = mu[p-1]/MELT(B,p-1,t-1,m,n);//multiplicity of fireable arc equals to mu marking divided by B[p-1][t-1]
				} else {
					if( MELT(B,p-1,t-1,m,n)<0) {// inhibitor arc   B[p-1][t-1]<0
						if( mu[p-1]>0 )//if mu > 0, fireable=0
							firable = 0;
						else  firable = INT_MAX;// if mu < 0, fireable= int_max
					} else firable = INT_MAX;
				}

				if (firable<firable_m) {//Compare to get the minimum  multiplicity of firable transitions
					firable_m = firable;
				}
			}
			f[t-1] = firable_m;
		}

		/*Algorithm3: Remove Low Priority Transition */
		#pragma omp parallel for private(i,j) num_threads(nth)
		for(i = 1; i <= n; i++) {
			#pragma omp simd
#pragma unroll
			for(j = 1; j <= n; j++) {
				if(MELT(R,i-1,j-1,n,n) > 0 && f[i - 1] != 0) { //R[i-1][j-1] > 0
					f[j - 1] = 0;
					//MELT(R,i-1,j-1,n,n) = 0;//R[i-1][j-1] = 0
				}
			}
		}

		#pragma omp parallel for reduction(+: nf) num_threads(nth)
		for(t = 1; t <= n; t++) {
			nf +=(f[t-1]>0) ? 1 : 0;
		}

		if(debug_level==2) {
			fprintf(fo,"vector f (%d firable transitions): ", nf);
			for(t=1; t<=n; t++) fprintf(fo,"%d ",f[t-1]);
			fprintf(fo,"\n");
		}
		if(nf==0) break;
		/*choose random fireable transition t*/
		ct = 0;
		rn = rand();
		l  = rn%nf;
		for(t=1; t<=n; t++) {
			if(f[t-1]>0) {
				if(l==ct) {
					firing_n = t;
					break;
				}
				ct++;
			}
		}
		if(debug_level==2) {
			fprintf(fo,"the random number is: %d, l=%d \n",rn,l);
			fprintf(fo,"the firing transition number is: %d\n",firing_n);
		}
		/**
		Algorithm2: The Marking Rule
		tb denotes the element in matrix B, and td denote the element in matrix D
		When tb < 0, it is inhibitor arc and put td*fm tokens to each of its output places.
		Otherwise, it puts td*fm tokens to each of its output place and extracts tb*fm to each of its input place.
		**/
		int dim = 0;
		fm = f[firing_n-1];
		if(debug_level==2)
			fprintf(fo,"firing transition multiplicity %d \n",fm);
		//	#pragma omp parallel for private(p,tb,td) num_threads(nth)
		for(p=1; p<=m; p++) {
			tb = MELT(B,p-1,firing_n-1,m,n); //B[p-1][firing_n-1]
			td = MELT(D,p-1,firing_n-1,m,n); //D[p-1][firing_n-1]
			//		#pragma omp critical
			if(tb<0) {
				mu[p-1] = mu[p-1]+td*fm;
			} else {
				mu[p-1] = mu[p-1]-tb*fm+td*fm;
			}
			if(mu[p-1] < 0) {
				err = ERR_OVERFLOW;
				perr("10");
				exit(1);
			} else if(mu[p-1] != 0) dim++;
		}


		/*print current marking mu*/
		if(debug_level==2 || s==smax) {
			fprintf(fo,"the current marking mu is:\n");
			if(printm==1) {
				for(p=1; p<=m; p++) {
					if(mu[p-1]!=0) {
						dim--;
						fprintf(fo,"%dp%d",mu[p-1],p);
						if(dim!=0)fprintf(fo,"+");
					}
				}
			} else {
				for(p=1; p<=m; p++) {
					fprintf(fo,"%d ",mu[p-1]);
				}
			}
			fprintf(fo,"\n");
			if(s==smax) {
				printf("###The result has been created! Specifies that the number of steps is %ld\n",smax);
				fprintf(fo,"###The result has been created! Specifies that the number of steps is %ld\n",smax);
				exit(0);
			}
		}

	}

}

/* output matrices*/
void printBDR(int m,int n,int *B,int *D,int *R,int debug_level,int rm,int *mu,int lsntonet,FILE *fo)
{
	int i,j;
	/* output matrices*/
	if(debug_level>0 || !rm ) {
		if(rm) {
			fprintf(fo,"matrix of transition incoming arcs B:\n");
			fprintf(fo,"(P\\T)\n");
		}
		for(i=0; i<m; i++) {
			for(j=0; j<n; j++) {
				fprintf(fo,"%d\t",MELT(B,i,j,m,n)); //B[i][j]
			}
			fprintf(fo,"\n");
		}
		fprintf(fo,"\n");

		if(rm) {
			fprintf(fo,"matrix of transition outgoing arcs D:\n");
			fprintf(fo,"(P\\T)\n");
		}
		for(i=0; i<m; i++) {
			for(j=0; j<n; j++) {
				fprintf(fo,"%d\t",MELT(D,i,j,m,n)); //D[i][j]
			}
			fprintf(fo,"\n");
		}
		fprintf(fo,"\n");

		if(rm) {
			fprintf(fo,"matrix of priority arcs R:\n");
			fprintf(fo,"(T\\T)\n");
		}
		for(i=0; i<n; i++) {
			for(j=0; j<n; j++) {
				fprintf(fo,"%d\t",MELT(R,i,j,n,n)); //R[i][j]
			}
			fprintf(fo,"\n");
		}
		fprintf(fo,"\n");
	}
}

/* output marking*/
void printmk(int m, int *mu,int printm,FILE *fo)
{
	int p;
	int c = 0;
	for(p=1; p<=m; p++) {
		if(mu[p-1]!=0) c++;
	}

	fprintf(fo,"The final marking is :\n");
	if(printm==1) {
		for(p=1; p<=m; p++) {
			if(mu[p-1]!=0) {
				c--;
				fprintf(fo,"%dp%d",mu[p-1],p);
				if(c!=0)fprintf(fo,"+");
			}
		}
		fprintf(fo,"\n");
	} else {
		for(p=1; p<=m; p++) {
			fprintf(fo,"%d ",mu[p-1]);
		}
		fprintf(fo,"\n");
	}
	if(fo != stdout ) fclose( fo );
}
//@ 2023 Qing Zhang: zhangq9919@163.com, Dmitry Zaitsev
