
//////////////////Lsn processor: low-level sleptsov net processor ////////////////
// run low-level sleptsov net                                                   //
// On single core and multicore parallel architectures using OpenMP             //
// Pass of parameters: -d debug_level -nth nth                                  //
// Command line:                                                                //
// >pl tr arc                                                                   //
// >initial marking                                                             //
// >arcs                                                                        //
// Command line:                                                                //
// SNP -d debug_level -nth nth <input_file > output_file                        //
//////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#define MATRIX_SIZE(d1,d2,t) ((d1)*(d2)*(sizeof(t)))
#define MOFF(i,j,d1,d2) ((d2)*(i)+(j))
#define MELT(x,i,j,d1,d2) (*((x)+MOFF(i,j,d1,d2)))
#define SKIP_COMM  {do { gets(input_buffer); } while(input_buffer[0]==';');}//skip comments 

/*Declare error code 声明错误代码*/
#define NO_ERR 0 /* No error */
#define ERR_NEGATIVE 1 /* Wrong header(m n k): negative numbers */ 
#define ERR_LESS 2 /* Wrong header(m n k): The number is less than 3 or non-numeric*/
#define ERR_H_ZERO 3 /* Wrong header(m n k): Value zero */
#define ERR_A_ZERO 4 /* Wrong arc: Value zero */
#define ERR_FORMAT 5 /* Wrong arc: Arc format error */
#define ERR_INHIBITOR 6 /* Wrong arc: The weight of inhibitior arc is not -1 */
#define ERR_EXCEEDM 7 /* Exceed number: Exceed the number of m */ 
#define ERR_EXCEEDN 8 /* Exceed number: Exceed the number of n */ 
#define ERR_DUP 9 /* Duplicate arc: The arc is duplicate */ 

static char Help[] =
"Lsn Processor\n" 
"action: run low-level sleptsov net\n"
"        on multicore parallel architectures using OpenMP.\n"
"usage:          lsn_pro_level [-h]\n" 
"                              [-d debug_level] [-nth nth] [-r]\n"
"-d debug_level                debug level:\n"
"                             (a): 'debug_level = 0'-----only final marking\n"
"                             (b): 'debug_level = 1'-----matrices of incoming and outgoing arcs and final marking\n"
"                             (c): 'debug_level = 2'-----all details\n"
"-nth nth                      number of thread:\n"
"                             (a): 'nth = 1'----1 thread\n"  
"                             (b): 'nth = 2'----2 threads\n"  
"                             (c): 'nth = 4'----4 threads\n"  
"                             (d): 'nth = 8'----8 threads\n"        
"-r                            sleptsov net rule\n";

char* errmsg[] ={/* 0 */ "No Error",
	/* 1 */ "Wrong header: Negative numbers",
	/* 2 */ "Wrong header: The number is less than 3 or non-numeric",
	/* 3 */ "Wrong header: Value zero",
	/* 4 */ "Wrong arc: The arc has value zero",
    /* 5 */ "Wrong arc: Arc format error",
    /* 6 */ "Wrong arc: The weight of inhibitior arc is not -1",
	/* 7 */ "Exceed number: Exceed the number of m",
	/* 8 */ "Exceed number: Exceed the number of n",
    /* 9 */	"Duplicate arc: The arc is duplicate",
	
	};

// declare global variables 
long err;//error code 
void perr(char* info)
{
	if( info ){
		printf("Error %s：%s\n",info,errmsg[err]);
		return;
	}
		printf("Error %s：%s\n",info,errmsg[err]);
}

double magma_wtime( void )
{
  struct timeval t;
  gettimeofday( &t, NULL );
  return t.tv_sec + t.tv_usec*1e-6;
}


//Function: Generate a matrix of priority arc chains 
int* tr(int *R,int n){
	int x,y,z;
  
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
        return R;	
}
  

int main (int argc,char *argv[])
{
	int m, // number of places
	    n, // number of transitions
	    k, // number of arcs
	    a,//intermediate variable
	    b;//intermediate variable
	int i,j,z;
	int v1,v2,v3;//p t w
	int *B,*D;//input matrix & output matrix
	int *R;//Priority arc matrix
	int *RC;//priority arc chain matrix 
	int debug_level=0;
	int nth;//threads
	int numf,snumf;
	double t1,t2;
    char input_buffer[60000];

  
  
    /* parse command line */
    numf=0;
    for( i=1; i<argc; i++ )
    {
      if( strcmp( argv[i], "-h" )==0 )
      {
        printf( "%s", Help ); return( 0);
      }
      if( strcmp( argv[i], "-d" )==0 ) { snumf=numf; numf=1; }
      else if( numf==1 ) { debug_level = atoi( argv[ i ] ); numf=snumf; printf("#####The debug level is %d\n",debug_level);}
      else if( strcmp( argv[i], "-nth" )==0 ) { snumf=numf; numf=2; }
      else if( numf==2 ) {nth = atoi( argv[ i ] ); numf=snumf; printf("#####The number of thread is %d\n",nth);}
      else if( strcmp( argv[i], "-r" )==0 ) {printf("#####The processor adopts the sleptsov net rule.\n");}
      else
        { printf( "*** unknown option: %s\n", argv[i] ); return(4); }
    } 
  
  
     
     SKIP_COMM
 //allocate space for matrices
  int digit = sscanf(input_buffer,"%d %d %d",&m,&n,&k);

  
  t1=magma_wtime();
  
 
  
 
  if(m < 0 || n < 0 || k < 0 ){err = ERR_NEGATIVE; perr("1"); exit(1); }
  if(digit != 3 ){err = ERR_LESS; perr("2"); exit(1); }
  if(m == 0 || n == 0 || k == 0){err = ERR_H_ZERO; perr("3"); exit(1);}
 
  if(debug_level > 0)  printf("m=%d n=%d k=%d\n",m,n,k);

 
  
  B = malloc( MATRIX_SIZE(m,n,int) );
  D = malloc( MATRIX_SIZE(m,n,int) );
  R = malloc( MATRIX_SIZE(n,n,int) );
  RC = malloc( MATRIX_SIZE(n,n,int) );
  
  if(B==NULL)
  {
     printf("***error: no enough memory for B\n");
     exit(2);
  }
  if(D==NULL)
  {
     printf("***error: no enough memory for D\n");
     exit(2);
  }
   if(R==NULL)
  {
     printf("***error: no enough memory for R\n");
     exit(2);
  }
   if(RC==NULL)
  {
     printf("***error: no enough memory for RC\n");
     exit(2);
  }
  
  
     for(i=0; i<m; i++)
  {
     for(j=0; j<n; j++)
      {
	     MELT(B,i,j,m,n) = 0; // B[i][j] = 0
	     MELT(D,i,j,m,n) = 0; // D[i][j] = 0
	    
      }
  }

  for(i=0; i<n; i++)
  {
     for(j=0; j<n; j++)
      {
	     MELT(R,i,j,n,n) = 0; // R[i][j] = 0
	     MELT(RC,i,j,n,n) = 0; // RC[i][j] = 0
	    
      }
  }


/*initializing & assignment*/ 
   int *mu; // marking (should be inputted from file)
   mu = malloc(m*sizeof(int));
   if(mu==NULL)
   {
     printf("***error: no enough memory for mu\n");
     exit(2);
   }
   
  
  
  int ii = 0;//pointer
  SKIP_COMM
  for(i=0; i<m; i++)
  { 
      sscanf(input_buffer+ii,"%d",&mu[i]);
    { 
      while(input_buffer[ii]==' ') ii++;// pass blanks 遇到空格，指针加一 
      while(input_buffer[ii]>='0' && input_buffer[ii]<='9') ii++;//pass digit 遇到一位数字，指针加一 
    }
  }
 
 
  if(debug_level > 0)
  {
  printf("marking is:(  ");
  for(i=0; i<m; i++)
  { printf("%d ",mu[i]);}
  printf("  )\n");
   }
  
  

  for(i=0;i<k;i++)
 {
     SKIP_COMM
   sscanf(input_buffer,"%d %d %d",&v1,&v2,&v3);// p t w
   
   if(debug_level > 0)  printf("%d %d %d\n",v1,v2,v3);
   if(v1 == 0 || v2 == 0) {err = ERR_A_ZERO; perr("4"); exit(1);}
   if(v1 > m)  {err = ERR_EXCEEDM; perr("7"); exit(1);}
   if(v2 > n)  {err = ERR_EXCEEDN; perr("8"); exit(1);}

  
   
 if(v1 > 0) 
        { if(v2 < 0 ) {printf("Warning: arc %d %d %d\n",v1,v2,v3);err = ERR_FORMAT; perr("5"); exit(1);}
        
		  a=v1,b=v2;
		    if( MELT(B,a-1,b-1,m,n) != 0) {printf("Warning: arc %d %d\n",a,b);err = ERR_DUP; perr("9"); exit(1);}
	    	if(v3>0) {
		    MELT(B,a-1,b-1,m,n)=v3; // B[a-1][b-1]=v3
			if(debug_level==2) printf("arc from place %d to transition %d of myltiplicity %d\n",v1,v2,v3);//v1,v3 are positive
			}
           
	        if(v3<0) {
	        if(v3 != -1 ) {err = ERR_INHIBITOR; perr("6"); exit(1);}
	        MELT(B,a-1,b-1,m,n)=v3; //B[a-1][b-1]=v3
	        if(debug_level==2) {printf("inhibitor arc from place %d to transition %d\n",v1,v2);} //v3 is negative,it is inhibitor
	       } 
	    }
 else if(v2 < 0){  if(v3 != 0) {printf("Warning: arc %d %d %d\n",v1,v2,v3);err = ERR_FORMAT; perr("5"); exit(1);}
                  a = abs(v1),b = abs(v2);
                  if( MELT(R,a-1,b-1,m,n) != 0) {printf("Warning: arc %d %d\n",a,b);err = ERR_DUP; perr("9"); exit(1);}
                   MELT(R,a-1,b-1,n,n) = 1;//R[a-1][b-1]=1
                   if(debug_level==2) printf("priority arc from transition %d to transition %d\n",abs(v1),abs(v2));// t --> t 
 	     	
                } 
      else{   if(v3 <= 0 ) {printf("Warning: arc %d %d %d\n",v1,v2,v3);err = ERR_FORMAT; perr("5"); exit(1);}
	          a = abs(v1),b = v2;
	          if( MELT(D,a-1,b-1,m,n) != 0) {printf("Warning: arc %d %d\n",a,b);err = ERR_DUP; perr("9"); exit(1);}
              MELT(D,a-1,b-1,m,n)=v3; //D[a-1][b-1]=v3
              if(debug_level==2) printf("arc from transition %d to place %d of myltiplicity %d\n",v2,abs(v1),v3);//v1 is negative(outgoing)
		  }
		
 } 
 
 

/*create input & output matrix*/
if(debug_level>0)
{
printf("matrix of transition incoming arcs B:\n");
printf("（P\\T）\n");
for(i=0;i<m;i++)
  {for(j=0;j<n;j++)
    {
	printf("%d  ",MELT(B,i,j,m,n)); //B[i][j]
    }
	printf("\n");
  } 
printf("matrix of transition outgoing arcs D:\n");
printf("（P\\T）\n");
for(i=0;i<m;i++)
 {for(j=0;j<n;j++)
    {
	printf("%d  ",MELT(D,i,j,m,n)); //D[i][j]
      }
	printf("\n");
 }
 
 printf("matrix of priority arcs R:\n");
printf("（T\\T）\n");
for(i=0;i<n;i++)
 {for(j=0;j<n;j++)
    {
	printf("%d  ",MELT(R,i,j,m,n)); //R[i][j]
     }
	printf("\n");
 }
}

//t1=magma_wtime();
 /*call function about priority arc chain*/ 
RC = tr(R,n);
//t2 = magma_wtime();


//run for v=argv[1] steps
 int v = atoi(argv[1]);
 int s; 				// current step number 当前步骤数 
 //int f[3]={0};           
 int *f;                 // indicator of fireable transitions 触发变迁的指标(0 or any numbers)
 int nf;				// number of fireable transitions 可触发变迁的总数 
 int p, t;				// place, transition
 int firable; //multiplicity of fireable transition 可触发变迁的权重 
 int firable_m; //the minimum multiplicity of fireable transition 可触发变迁权重的最小值 
 int rn; //random number 
 int firing_n; //firing transition number  正在触发的变迁个数 
 int l; //remainder 余数 
 int tb,td; //columns of matrix B and D 
 int ct;  // count of fireable transitions
 int fm;  //正在触发变迁的权重 multiplicity of firing transition
 
 
/*firing rule of SN*/
 f = malloc(n*sizeof(int));
  
  if(f==NULL)
  {
     printf("***error: no enough memory for mu\n");
     exit(2);
  }


 for(s=1;;s++)
 { 
    if(debug_level==2)
 	printf("--- step %d ---\n",s);
 	// form vector of fireable transitions f (0 - not fireable, numbers - fireable)
 	nf = 0;
 	

#pragma omp parallel for private(p,t,firable) reduction(min:firable_m) num_threads(nth)	
 for(t=1;t<=n;t++) 
 	{
 		firable_m = INT_MAX;//Set the initial value to a larger value

 		for(p=1;p<=m;p++)
 		{
 			if( MELT(B,p-1,t-1,m,n)>0 )         // regular arc  B[p-1][t-1]>0
 			  { 
 			     firable = mu[p-1]/MELT(B,p-1,t-1,m,n);
				  //multiplicity of fireable arc equals to mu marking divided by B[p-1][t-1]
 			  
			   }
	        else {if( MELT(B,p-1,t-1,m,n)<0)			// inhibitor arc   B[p-1][t-1]<0
			        { if( mu[p-1]>0 )             //if mu > 0, fireable=0          
			           firable = 0; 
					  else  firable = INT_MAX;     // if mu < 0, fireable= ∞   
				     } 
                  else firable = INT_MAX;             
				  } 

		if (firable<firable_m)           //Compare to get the minimum  multiplicity of firable transitions
			  {firable_m = firable;} 
		 } 
		f[t-1] = firable_m;	
	} 


/*remove low priority 移除低优先级弧 */
 for(i = 1;i <= n;i++){
 	for(j = 1;j <= n;j++){
 		if(MELT(RC,i-1,j-1,n,n) > 0 && f[i - 1] != 0) //R[i-1][j-1] > 0 
 		   {f[j - 1] = 0;  //置0，使其不触发 
 		   //MELT(R,i-1,j-1,n,n) = 0;//令 R[i-1][j-1] = 0，下一步中不再移除低优先级的弧 
			}
	     }   
	 }
	 
	 for(t = 1; t <= n;t++)
	 {
	 nf +=(f[t-1]>0) ? 1 : 0; 
	  } 
	




	if(debug_level==2)
{
	printf("vector f (%d firable transitions): ", nf);
	for(t=1;t<=n;t++) printf("%d ",f[t-1]);
	printf("\n");
}
	if(nf==0) break; 
	



 
/*choose random fireable transition t*/
    ct = 0;
	rn = rand();
	l  = rn%nf;//余数为0，则选择触发变迁t为第一个不为0的变迁； 余数为1，则为第二个，以此类推 
	
	

	
	for(t=1;t<=n;t++)
	{
		if(f[t-1]>0) 
		{
			if(l==ct){firing_n = t;	break;}
    //ct作为指针，指代第0个...第n-1个变迁，只有当第ct个变迁(f[t-1])触发，且余数等于ct时,选择此变迁为firing_n 
			ct++;
		}
	}
	if(debug_level==2)
	{
	printf("the random number is: %d, l=%d \n",rn,l);
	printf("the firing transition number is: %d\n",firing_n);
	}
	
	
	
	
	
//marking rules     fire transition t: 1) reduce marking mu on B; 2) increase marking mu on D
    fm = f[firing_n-1];
    
     if(debug_level==2)
    printf("firing transition multiplicity %d \n",fm);
 #pragma omp parallel for private(p,tb,td) num_threads(nth)
     for(p=1;p<=m;p++)
    {
     tb = MELT(B,p-1,firing_n-1,m,n); //B[p-1][firing_n-1]
     td = MELT(D,p-1,firing_n-1,m,n); //D[p-1][firing_n-1]
     
       { if(tb<0){
	       mu[p-1] = mu[p-1]+td*fm; 
	              }
	     else
            {
	      mu[p-1] = mu[p-1]-tb*fm+td*fm; 
	         } 
	   }
       
	}
	
	
/*print current marking mu*/
 if(debug_level==2)
 	{
	 printf("the current marking mu is:\n");
    for(p=1;p<=m;p++)
      printf("%d ",mu[p-1]);
    printf("\n");
    }

  } 
  
 free(B);
 free(D);

 printf("The final marking is :\n");
    for(p=1;p<=m;p++)
      printf("%d ",mu[p-1]);
      
  t2=magma_wtime();
  //if(debug_level>0) 
  printf( "\n##########The total time is %#fs\n", t2-t1 );
} /*main*/

