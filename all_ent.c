#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <fftw3.h>
#include <time.h>

SEXP fftw_transform(SEXP input){
	int n = length(input);
	fftw_complex *in, *out;
	fftw_plan p;
	SEXP result;
	PROTECT(result=allocMatrix(REALSXP,n,2));
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

	for (int i =0; i<n;i++){
		in[i][0]=REAL(input)[i];
		in[i][1]=0.0;
	
	}
	p = fftw_plan_dft_1d(n,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for (int i =0; i<n;i++){
		REAL(result)[i]=out[i][0];
		REAL(result)[i+n]=out[i][1];

	}
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	UNPROTECT(1);
	return result;

}

void swap(double *a, double *b){
	double temp=*a;
	*a=*b;
	*b=temp;
}

void shuffle(double *array_mag, int* n_ptr){
	int n = n_ptr[0];
	srand(time(NULL));
	for (int i =n-1; i>0; i--){
		int j = rand()%(i+1);
		swap(&array_mag[i],&array_mag[j]);
	}
}

void sampen(double *y, int * M_int, double * r_double, int * n_int, double * SampEn ){

	int M = M_int[0];
	double r = r_double[0];
	int  n = n_int[0];

	double *p =NULL;
	long *run = NULL, *lastrun = NULL, N;
	double *A = NULL, *B =NULL;
	int M1, m;
	unsigned long i,j,jj,nj;
	double y1;
	M++;

	if ((run = (long *) calloc(n, sizeof(long))) == NULL)
		exit(1);
	if ((lastrun = (long *) calloc(n, sizeof(long))) == NULL)
	    	exit(1);
	if ((A = (double *) calloc(M, sizeof(double))) == NULL)
		exit(1);
	if ((B = (double *) calloc(M, sizeof(double))) == NULL)
	   	exit(1);
	if ((p = (double *) calloc(M, sizeof(double))) == NULL)
		exit(1);
	for (i=0; i<n-1; i++){
		nj=n-i-1;
		y1=y[i];
		for (jj=0; jj<nj; jj++){
			j=jj+i+1;
			if (((y[j]-y1)<r)&&((y1-y[j])<r)){
				run[jj]=lastrun[jj]+1;
				M1=M<run[jj]? M : run [jj];
				for (m=0;m<M1; m++){
					A[m]++;
					if(j<n-1)
						B[m]++;
				}
			}
			else 	
				run[jj]=0;
		}
		for (j=0; j<nj;j++){
			
			lastrun[j]=run[j];
		}
		}
		N=(long)(n*(n-1)/2);
		p[0]=A[0]/N;
		SampEn[0]=-log(p[0]);
		for (m=1;m<M;m++){
			p[m]=A[m]/B[m-1];
			if (p[m]==0){
				SampEn[m]=1.;
			}
			else{ 
				SampEn[m]= -log(p[m]);
				*SampEn=SampEn[m]	;
			}
		}
		free(A);
		free(B);
		free(p);
		free(run);
		free(lastrun);
	}


//Another method for SampEn 
#define M_MAX 20
double SampleEntropy(double *y, unsigned long int nlin, int m, double r, int scale)
{
    int j = scale;
    int m_max = m;
    unsigned long int i, k, l, nlin_j; 
    unsigned long int cont[M_MAX+1];
    double r_new;
    double SE = -1.0;
    
   nlin_j = (unsigned long) ((nlin/j) - m_max); 
	r_new = r;              
                                        for (i = 0; i <= M_MAX; i++)
                       cont[i]=0;
                       for (i = 0; i < nlin_j; ++i) {
                       for (l = i+1; l < nlin_j; ++l) { /*self-matches are not counted*/
                       k = 0;
                       while (k < m_max && fabs(y[i+k] - y[l+k]) <= r_new)
                       cont[++k]++;
                       if (k == m_max && fabs(y[i+m_max] - y[l+m_max]) <= r_new)
                       cont[m_max+1]++;
                       } 
                       }     
                       for (i = 1; i <= m_max; i++)
                       if (cont[i+1] == 0 || cont[i] == 0)
                       SE = -log((double)1/((nlin_j)*(nlin_j-1)));
                       else
  		SE = -log((double)cont[i+1]/cont[i]);
                if (SE < 0.0){
                printf("Error: SampEn is negative: %lf\n",SE);
                printf("Report:\n");
                printf("nlin=%lu, j=%d, nlin_j=%lu, m=%d, r=%lf\n",nlin,j,nlin_j,m,r);
		for (i = 0; i < M_MAX; i++)
 		printf("cont[%lu] = %lu, cont[%lu] = %lu\n",i,cont[i],i+1,cont[i+1]);
                }
                return SE;
                }
 		

void CoarseGraining(double *y, unsigned long int nlin, int j, double *data){
	int i, k;

        for (i = 0; i < (unsigned long)nlin/j; i++) {
		y[i] = 0;
		for (k = 0; k < j; k++)
                	y[i] += data[i*j+k];
	y[i] /= j; 
}
}
void MSE(double *data,  int *M_int, double *r_double, int *n_int, int *max_scale_int, double * ME){
	int M = M_int[0];
	double r = r_double[0];
	int  n = n_int[0];
	int max_scale = max_scale_int[0];

	int i;
	    
	double *y = (double *)malloc((size_t)n*sizeof(double));
	if(y == NULL){
		fprintf(stderr,"Error: memory allocation problem\n");
	        return; 
	}
             
	//Perform the loop over the scales
 	for(i = 1; i <= max_scale; i++){
	           CoarseGraining(y, n, i, data);
	           ME[i-1] = SampleEntropy(y, n, M, r, i);
		   *ME=ME[i-1];
	}
	free(y);
}
double wasserstein_distance(double *p, double *q , int n ){
	double *cumulative_p=(double *)malloc(n*sizeof(double)); 
	double *cumulative_q=(double *)malloc(n*sizeof(double));
	cumulative_p[0]=p[0];
	cumulative_q[0]=q[0];
	for(int i =1; i<n; i++) {
		cumulative_p[i]=cumulative_p[i-1]+p[i];
		cumulative_q[i]=cumulative_q[i-1]+q[i];
	}
	double distance =0.0;
	for(int i=0; i<n ;i++){
		distance+=fabs(cumulative_p[i]-cumulative_q[i]);
	
	}
	free(cumulative_p);
	free(cumulative_q);
	return distance; 

}

void wasserstein_distance_r(double *p, double *q, int *n, double *result){
	*result=wasserstein_distance(p,q,*n);

}
