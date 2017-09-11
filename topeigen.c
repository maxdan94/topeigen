/*
Maximilien Danisch
Septembre 2017
http://bit.ly/danisch
maximilien.danisch@gmail.com

to compile:
gcc topeigen.c -o topeigen -O9 -lm

to execute:
./topeigen net.txt k res.txt
net.txt should contain on each line: "i j value". The matrix has to be symetric (if "i j" is here then "j i" is not here)
res.txt contains on column k the k highest eigenvalue (in absolute value) followed by the entries of its eigenvector.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>


#define NLINKS 10000000 //maximum number of nonzero entries, will increase if needed
#define REP 30 //Number of power iterations. SHOULD BE MANUALLY INCREASED IF PRECISION IS NOT ENOUGH

// graph datastructure:
typedef struct {
	unsigned s;//source node
	unsigned t;//target node
	double w;//edge weight
} edge;

//sparse matrix structure
typedef struct {
	unsigned n;//dimensions of the squared matrix
	unsigned e;//number nonzero entries symmetric (edges)
	edge *el;//edge list
} sparse;

//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the weighted sparsematrix
sparse* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	g->el=malloc(e1*sizeof(edge));
	FILE *file;
	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	while (fscanf(file,"%u %u %le", &(g->el[g->e].s), &(g->el[g->e].t), &(g->el[g->e].w))==3) {
		g->n=max3(g->n,g->el[g->e].s,g->el[g->e].t);
		if (g->e++==e1) {
			e1+=NLINKS;
			g->el=realloc(g->el,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->el=realloc(g->el,g->e*sizeof(edge));

	return g;
}


//free the graph stucture
void freegraph(sparse *g){
	free(g->el);
	free(g);
}

//initialize to random vector
void randini(unsigned n, double* v){
	unsigned i;
	srand(time(NULL));
	for (i=0;i<n;i++){
		v[i]=((float)rand()/(float)(RAND_MAX));
	}

}

//product of matrix g and vector v1 -> strored in v2
void prod(sparse* g, double* v1, double* v2){
	unsigned i;
	bzero(v2,sizeof(double)*g->n);
	for (i=0;i<g->e;i++){
		v2[g->el[i].s]+=v1[g->el[i].t]*g->el[i].w;
		v2[g->el[i].t]+=v1[g->el[i].s]*g->el[i].w;
	}
}

void normalize(unsigned n, double* vect){
	unsigned i;
	double s=0;
	for (i=0;i<n;i++){
		s+=vect[i]*vect[i];
	}
	for (i=0;i<n;i++){
		vect[i]/=sqrt(s);
	}
}

double scallarproduct(unsigned n, double *v1, double *v2){
	unsigned i;
	double s=0;
	for (i=0;i<n;i++){
		s+=v1[i]*v2[i];
		
	}
	return s;
}

//removing components of earlier eignenvectors
void project(unsigned n, double* v,unsigned k, double** vk){
	unsigned i,j;
	double s;
	for (i=0;i<k;i++){
		s=scallarproduct(n,v,vk[i]);
		for (j=0;j<n;j++){
			v[j]-=vk[i][j]*s;
		}
	}
}

//to compute eigenvalue
double ratio(unsigned n, double* v1, double* v2){
	double s1=0,s2=0;
	unsigned i;
	for (i=0;i<n;i++){
		s1+=v1[i];
		s2+=v2[i];
	}
	return s2/s1;
}


//returns k heighest eigenvectors and eigenvalues
void topeigen(sparse* g, unsigned k, double *val, double** vk){
	double *v1,*v2,*v3;
	unsigned i,j;

	v1=malloc(g->n*sizeof(double));
	v2=malloc(g->n*sizeof(double));

	for (i=0;i<k;i++){
		randini(g->n,v1);//random initialisation
		for (j=0;j<REP;j++){
			prod(g,v1,v2);//product stored in v2
			project(g->n,v2,i,vk);//removing components of earlier eignenvectors
			normalize(g->n,v2);//normalizing
			v3=v2,v2=v1,v1=v3;
		}
		prod(g,v1,v2);
		project(g->n,v2,i,vk);
		val[i]=ratio(g->n,v1,v2);//computing eigenvalue
		normalize(g->n,v2);
		for (j=0;j<g->n;j++){
			vk[i][j]=v2[j];
		}
	}

	free(v1);
	free(v2);
}


void printres(FILE* file,unsigned n, unsigned k, double** vk, double* val){
	unsigned i,j;

	for (i=0;i<k-1;i++){
		fprintf(file,"%10le ",val[i]);
	}
	fprintf(file,"%10le\n",val[i]);

	for (j=0;j<n;j++){
		for (i=0;i<k-1;i++){
			fprintf(file,"%10le ",vk[i][j]);
		}
		fprintf(file,"%10le\n",vk[i][j]);
	}
}


int main(int argc,char** argv){
	sparse* g;
	time_t t1,t2,t3;
	unsigned k=atoi(argv[2]),i;
	double** vk;
	double* val;
	FILE* file=fopen(argv[3],"w");

	t1=time(NULL);
	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Number of nodes = %u\n",g->n);
	printf("Number of edges = %u\n",g->e);

	printf("Computing %u largest eigenvalues and associated eignenvectors\n",k);

	vk=malloc(k*sizeof(double*));
	for (i=0;i<k;i++){
		vk[i]=malloc(g->n*sizeof(double));
	}
	val=malloc(k*sizeof(double));
	topeigen(g,k,val,vk);

	printres(file,g->n,k,vk,val);
	fclose(file);

	freegraph(g);
	return 0;
}




