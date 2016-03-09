#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>

using namespace std;

#define SYSTEMTIME clock_t

 
void initMatrices(double ** pha, double ** phb, double ** phc, int m_ar, int m_br){
	long int i,j;

        
	*pha = (double *)malloc((m_ar * m_br) * sizeof(double));
	*phb = (double *)malloc((m_ar * m_br) * sizeof(double));
	*phc = (double *)malloc((m_ar * m_ar) * sizeof(double));


        for(i=0; i< m_ar; i++)
		for(j=0; j< m_br; j++){
			(*pha)[i*m_br + j] = (double)1.0;
                }

                
	for(i=0; i< m_br; i++)
		for(j=0; j< m_ar; j++){
                    (*phb)[i*m_ar + j] = (double)(i+1);
                }

}

void printResultMatrix(double ** phc, int m_br) {
	int i,j;
            
	cout << endl << "Result matrix: " << endl;
	for(i=0; i<1; i++)
	{	for(j=0; j<min(10,m_br); j++)
			cout << (*phc)[j] << " ";
	}
	cout << endl;
}

void printTimeUntilNow(SYSTEMTIME time1){
	char st[100];

	sprintf(st, "Time: %3.3f seconds\n", (double)(clock() - time1) / CLOCKS_PER_SEC);
	cout << st;
}


void freeMatrices(double * pha, double * phb, double * phc){
    free(pha);
    free(phb);
    free(phc);
}


void multMatrices(double * pha, double * phb, double * phc, int m_ar, int m_br){
	double temp;
	int i, j, k;

	for(i=0; i<m_ar; i++)
	{	for( j=0; j<m_ar; j++)
		{	temp = 0;
			for( k=0; k<m_br; k++)
			{	
				temp += pha[i*m_br+k] * phb[k*m_ar+j];
			}
			phc[i*m_ar+j]=temp;
		}
	}
}



void multMatricesLineByLine(double * pha, double * phb, double * phc, int m_ar, int m_br)
{
	double temp;
	int i, j, k;

        
        //TODO: fix thiz to line by line algorithm.
	for(i=0; i<m_ar; i++)
	{	
                for( k=0; k<m_br; k++)
		{	
                        for( j=0; j<m_ar; j++)
			{	
				phc[i*m_ar+j] = k==0 
                                                ? (pha[i*m_br+k] * phb[k*m_ar+j]) 
                                                : phc[i*m_ar+j] + (pha[i*m_br+k] * phb[k*m_ar+j]);
			}
		}
	}
}

void multiplyMatrices(int m_ar, int m_br, int op) 
{
	double *pha, *phb, *phc;
	SYSTEMTIME time1;

	time1 = clock();
	initMatrices(&pha,&phb,&phc, m_ar,m_br);
	if(op==1){
            multMatrices(pha,phb,phc,m_ar,m_br);	
	}else {
            multMatricesLineByLine(pha,phb,phc,m_ar,m_br);
	}
        printTimeUntilNow(time1);
	printResultMatrix(&phc,m_br);
	freeMatrices(pha,phb,phc);
}

float produtoInterno(float *v1, float *v2, int col)
{
	int i;
	float soma=0.0;	

	for(i=0; i<col; i++)
		soma += v1[i]*v2[i];
	
	return(soma);

}

void handle_error (int retval)
{
  printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
  exit(1);
}

void init_papi() {
  int retval = PAPI_library_init(PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT && retval < 0) {
    printf("PAPI library version mismatch!\n");
    exit(1);
  }
  if (retval < 0) handle_error(retval);

  std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
            << " MINOR: " << PAPI_VERSION_MINOR(retval)
            << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}


int main (int argc, char *argv[])
{

	int lin, col;
	int op;
	
	int EventSet = PAPI_NULL;
  	long long values[2];
  	int ret;
	

	ret = PAPI_library_init( PAPI_VER_CURRENT );
	if ( ret != PAPI_VER_CURRENT )
		std::cout << "FAIL" << endl;


	ret = PAPI_create_eventset(&EventSet);
		if (ret != PAPI_OK) cout << "ERRO: create eventset" << endl;


	ret = PAPI_add_event(EventSet,PAPI_L1_DCM );
	if (ret != PAPI_OK) cout << "ERRO: PAPI_L1_DCM" << endl;


	ret = PAPI_add_event(EventSet,PAPI_L2_DCM);
	if (ret != PAPI_OK) cout << "ERRO: PAPI_L2_DCM" << endl;


	op=1;
	do {
		cout << endl << "1. Multiplication" << endl;
		cout << "2. Line Multiplication" << endl;
		cout << "Selection?: ";
		cin >>op;
		if (op == 0)
			break;	
		printf("Dimensions: lins cols ? ");
   		cin >> lin >> col;



		// Start counting
		ret = PAPI_start(EventSet);
		if (ret != PAPI_OK) cout << "ERRO: Start PAPI" << endl;


		multiplyMatrices(lin,col,op);

  		ret = PAPI_stop(EventSet, values);
  		if (ret != PAPI_OK) cout << "ERRO: Stop PAPI" << endl;
  		printf("L1 DCM: %lld \n",values[0]);
  		printf("L2 DCM: %lld \n",values[1]);

		ret = PAPI_reset( EventSet );
		if ( ret != PAPI_OK )
			std::cout << "FAIL reset" << endl; 



	}while (op != 0);

		ret = PAPI_remove_event( EventSet, PAPI_L1_DCM );
		if ( ret != PAPI_OK )
			std::cout << "FAIL remove event" << endl; 

		ret = PAPI_remove_event( EventSet, PAPI_L2_DCM );
		if ( ret != PAPI_OK )
			std::cout << "FAIL remove event" << endl; 

		ret = PAPI_destroy_eventset( &EventSet );
		if ( ret != PAPI_OK )
			std::cout << "FAIL destroy" << endl;

}
