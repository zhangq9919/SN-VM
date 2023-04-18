// generator of .hsn for computing polynomials  n=k
//
// USAGE: >gen_pol 5 pol5.txt
// k=2 a2x^2+a1x+a0

#include <stdio.h>
#include <stdlib.h>
#define rand_term (rand()%10+1)   //生成随机数1-10范围 

int main(int argc, char *argv[])
{
	int k=atoi(argv[1]);
	int nt = k+k*(k+1)/2;
	int np = 3*nt+2;
	int na = 6*nt;
	int l=0;
	int mu[np];
	int NST = nt;
	int r,t,off;
	int x = rand_term;
	int c = k;
	FILE *f;
	f = fopen(argv[2], "w" );
	if(k==0) {
		fprintf(f,"k cannot equals to 0!");
		exit(1);
	}


	// making
	//Start 1个
	mu[0] = 0;
	//fprintf(f,"%d ", 0); // mu(p_1)
	//Finish nt个
	for(t=1; t<=nt; t++)
		{mu[t] = 1;
		}
	//fprintf(f,"%d ", 1); // mu(p_2)...mu(p_{nt+1})
	int ntt = nt;
	while(c!=0) {
		r=rand_term;//ak
		//fprintf(f,"%d ", r); // mu(p_{nt+2})
		mu[++ntt] = r;
		for(t=1; t<=c; t++) {
			//fprintf(f,"%d %d ",x,0);//mu{p_nt+3}.....mu{p_}
			mu[++ntt] = x;
			mu[++ntt] = 0;
		}
		if(c<k)	mu[++ntt] = 0;
		//fprintf(f,"%d ",0);//ak不需要多一个0存放结果
		c--;
	}
	r=rand_term;
	mu[++ntt] = r;
	mu[++ntt] = 0;
	//fprintf(f,"%d %d ", r,0); //a1
	//fprintf(f,"\n");

	for(t = 0; t < np; t++) {
		if(mu[t]!=0) l++;
	//	printf("%d ",mu[t]);
	}

//header np nt na l NST
	fprintf(f,"%d %d %d %d %d\n",np,nt,na,l,NST);

	// arcs
	// control flow arcs
	for(t=1; t<=nt; t++) {
		fprintf(f,"%d %d %d\n", t, t, -1); // p->t inh
		fprintf(f,"%d %d %d\n", -t, t, 1); // p<-t
		fprintf(f,"%d %d %d\n", t+1, t, 1); // t<-p
	}
// dashed arcs
	off=nt+1;//4
	for(t=1; t<=nt; t++) {
		fprintf(f,"%d %d %d\n", t+off, t, 1); // p->t
		fprintf(f,"%d %d %d\n", t+off+1, t, 1); // p->t
		fprintf(f,"%d %d %d\n", -(t+off+2), t, 1); // t->p
		off+=1;
	}
	//marking line:p mu
	for(t=1; t<=np; t++) {
		if(mu[t-1]!=0) fprintf(f,"%d %d\n",t,mu[t-1]);
	}
	int t1 = 1;
	int t2;
	int p[3];
	int i = 0;
	int b = 0;
	c = k;
	while(c!=0) {
		for(t=1; t<=c; t++) {

			fprintf(f,"%d %d %s\n", t1, 5, "mul_lsn.txt"); //substitution info
			fprintf(f,"%d %d\n", nt+2*t1,1);
			fprintf(f,"%d %d\n", nt+2*t1+1,2);
			fprintf(f,"%d %d\n", nt+2*t1+2,-4);
			fprintf(f,"%d %d\n", -t1,3);
			fprintf(f,"%d %d\n", -(t1+1),-5);
			if(t==c) {
				p[i] = nt+2*t1+2;
				if(i<1)i++;
			}
			t1++;

		}
		if(c<k) {
			fprintf(f,"%d %d %s\n", t1, 5, "add_lsn.txt"); //substitution info
			fprintf(f,"%d %d\n", p[0],1);
			fprintf(f,"%d %d\n", p[1],2);
			p[2]=p[1]>p[0]? p[1]:p[0];
			fprintf(f,"%d %d\n", p[2]+1,-3);
			fprintf(f,"%d %d\n", -t1,4);
			fprintf(f,"%d %d\n", -(t1+1),-5);
			p[i] = p[2]+1;
			i = 0;
			//t2=t1;
			t1++;
		}
		if(c == k) nt++;
		c--;

	}
	//t1++;
	fprintf(f,"%d %d %s\n", t1, 5, "add_lsn.txt"); //substitution info
	fprintf(f,"%d %d\n", p[2]+1,1);
	fprintf(f,"%d %d\n", p[2]+2,2);
	fprintf(f,"%d %d\n", p[2]+3,-3);
	fprintf(f,"%d %d\n", -t1,4);
	fprintf(f,"%d %d\n", -(t1+1),-5);

}

