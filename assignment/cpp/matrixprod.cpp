#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>
#include <fstream>
#include <sstream>


using namespace std;

#define SYSTEMTIME clock_t

ofstream logFile;


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
                
        for(i=0; i< m_ar; i++)
		for(j=0; j< m_ar; j++){
                    (*phc)[i*m_ar + j] = 0;
                }

}

void printResultMatrix(double ** phc, int m_br) {
	int i,j;
            
	cout << "Result matrix: " << endl;
	for(i=0; i<1; i++)
	{
            for(j=0; j<min(10,m_br); j++)
                cout << (*phc)[j] << " ";
	}
	cout << endl << endl;
}

void printLogFileHeaders(){
    if(logFile.is_open()){
        logFile << "Function ID" << ";";
        logFile << "Rows" << ";";
        logFile << "Cols" << ";";
        logFile << "Duration (seconds)" << ";";
        logFile << "Performance (MFLOPS)" << ";";
        logFile << "L1 DCM" << ";";
        logFile << "L1 DCM" << ";";
        logFile << "\n";

    }
}

void printResults(SYSTEMTIME time1, int op, int m_ar, int m_br){
	char st[100];
        double delta = (double)(clock() - time1) / CLOCKS_PER_SEC, performance;

        if(logFile.is_open()){ 
            logFile << op << ";";
            logFile << m_ar << ";";
            logFile << m_br << ";";
        }
	sprintf(st, "Time: %3.3f seconds\n", delta);
	cout << st;
        if(logFile.is_open()) logFile << delta << ";";
        
        performance = (((double)3*m_ar*m_br*m_br)*1000000)/delta;
        sprintf(st, "Prf.: %3.3f MFLOPS\n", performance);
	cout << st;
        if(logFile.is_open()) logFile << performance << ";";
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
	{	
                for( j=0; j<m_ar; j++)
		{	temp = 0;
			for( k=0; k<m_br; k++)
			{	
				temp += pha[i*m_ar+k] * phb[k*m_br+j];
			}
			phc[i*m_ar+j]=temp;
		}
	}
	
}



void multMatricesOpt(double * pha, double * phb, double * phc, int m_ar, int m_br){
	double temp;
	int i, j, k;

	for( j=0; j<m_ar; j++)
	{	
                for(i=0; i<m_ar; i++)
		{	temp = 0;
			for( k=0; k<m_br; k++)
			{	
				temp += pha[i*m_ar+k] * phb[k*m_br+j];
			}
			phc[i*m_ar+j]=temp;
		}
	}
}



void transpose(double * matrix,int nRows,int nColumns, double * transpose){
    
    int i, j;
    for(i = 0; i < nRows; i++){
        for(j = 0; j < nColumns; j++){
            transpose[j*nRows+i] = matrix[i*nColumns+j];
        }
    }
    
}

void multMatricesLineByLine(double * pha, double * phb, double * phc, int m_ar, int m_br)
{
        int i, j, k;
        
    
        for(i=0; i<m_ar; i++) 
	{	
                for( k=0; k<m_br; k++) 
		{	 
                        for( j=0; j<m_ar; j++) 
			{	
				phc[i*m_ar+j] += (pha[i*m_ar+k] * phb[k*m_br+j]);
                                
			}
			
		}
	}
	

}

void multMatricesLineByLineByTranspose(double * pha, double * phb, double * phc, int m_ar, int m_br)
{
        double temp;
        int i, j, k;
        
        double * phbT = (double *)malloc((m_ar * m_br) * sizeof(double));
        transpose(phb, m_br, m_ar, phbT);
    
        for(i=0; i<m_ar; i++) 
	{	
                for( j=0; j<m_ar; j++) 
		{	 
                        
			 for( k=0; k<m_br; k++)
			{	
                                temp = 0;
                                for( k=0; k<m_br; k++)
                                {	
                                        temp += pha[i*m_ar+k] * phbT[i*m_ar+k]; 
                                }
                                phc[i*m_ar+j]=temp; 
			}
			
		}
	}
	
	free(phbT);
}




void multiplyMatrices(int m_ar, int m_br, int op) 
{
	double *pha, *phb, *phc;
	SYSTEMTIME time1;

	initMatrices(&pha,&phb,&phc, m_ar,m_br);
	
            
        time1 = clock();
        switch(op){
            case 1:
                multMatrices(pha,phb,phc,m_ar,m_br);
                break;
            case 2:
                multMatricesOpt(pha,phb,phc,m_ar,m_br);
                break;
            case 3: 
                multMatricesLineByLine(pha,phb,phc,m_ar,m_br);
                break;
        }
        printResults(time1,op, m_ar,m_br);
	printResultMatrix(&phc,m_br);
	
        freeMatrices(pha,phb,phc);
}

bool showMenu(){ 
    int op, lin, col;
    
    
    // menu dialog
    cout << endl;
    cout << endl;
    cout << "======================================" << endl;
    cout << "1. Multiplication" << endl;
    cout << "2. Multiplication (optimized)" << endl;
    cout << "3. Line Multiplication" << endl;
    cout << "0. End" << endl << endl;
    cout << "Selection?: ";
    cin >> op;
    
    // preliminar option validation
    if( op < 0 || op > 3) {
        cout << "Bad option. " << endl;
        return true;
    } else if (op == 0){ // exit program
        return false;
    }
    
    // secondary question
    printf("Dimensions: lins cols ? ");
    cin >> lin >> col;
    cout << endl;
    
    
    // preforming operations
    cout << "Working ..." << endl;
    multiplyMatrices(lin,col,op);
    
    return true;
}



void papiHandleError (int ret)
{
  printf("PAPI error %d: %s\n", ret, PAPI_strerror(ret));
  exit(1);
}

int papiSubscribeEvents(){
    int ret, eventSet = PAPI_NULL;
    
    
    //init papi
    ret = PAPI_library_init(PAPI_VER_CURRENT);
    if (ret != PAPI_VER_CURRENT && ret < 0) {
        printf("PAPI library version mismatch!\n");
        exit(1);
    }
    if (ret < 0) papiHandleError(ret);

    std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(ret)
                << " MINOR: " << PAPI_VERSION_MINOR(ret)
                << " REVISION: " << PAPI_VERSION_REVISION(ret) << "\n";
                
    // create event set
    ret = PAPI_create_eventset(&eventSet);
    if (ret != PAPI_OK) cout << "ERRO: create eventset" << endl;


    // subscribe events
    ret = PAPI_add_event(eventSet,PAPI_L1_DCM );
    if (ret != PAPI_OK) cout << "ERRO: PAPI_L1_DCM" << endl;

    ret = PAPI_add_event(eventSet,PAPI_L2_DCM);
    if (ret != PAPI_OK) cout << "ERRO: PAPI_L2_DCM" << endl;
    
    return eventSet;
}

void papiUnsubscribeEvents(int eventSet){
    int ret;
    // unsubscribe events
    ret = PAPI_remove_event( eventSet, PAPI_L1_DCM );
    if ( ret != PAPI_OK ) cout << "FAIL remove event" << endl; 

    ret = PAPI_remove_event( eventSet, PAPI_L2_DCM );
    if ( ret != PAPI_OK ) cout << "FAIL remove event" << endl; 

    // release eventset
    ret = PAPI_destroy_eventset( &eventSet );
    if ( ret != PAPI_OK ) cout << "FAIL destroy" << endl;
}


void papiStartCounter(int eventSet){
    int ret;
    // Start counting
    ret = PAPI_start(eventSet);
    if (ret != PAPI_OK) cout << "ERRO: Start PAPI" << endl;
}

void papiShowResults(int eventSet, bool showResults){
    int ret;
    long long values[2];
    
    // show papi results
    ret = PAPI_stop(eventSet, values);
    if (ret != PAPI_OK) cout << "ERRO: Stop PAPI" << endl;
    
    if(showResults){
        cout << "PAPI results:" << endl;
        printf("L1 DCM: %lld \n",values[0]);
        printf("L2 DCM: %lld \n",values[1]);
        
        if(logFile.is_open()){
            logFile << values[0] << ";" << values[1] << ";\n";
        }
    }
    
    

    // reset counters
    ret = PAPI_reset( eventSet );
    if ( ret != PAPI_OK ) cout << "FAIL reset" << endl; 
    
}


int main (int argc, char *argv[])
{        
     
        //create log file
        logFile.open("log.txt");
        printLogFileHeaders();
        
        
        bool cntinue = false;
        int eventSet = papiSubscribeEvents();
        do{
            papiStartCounter(eventSet);
            cntinue = showMenu();
            papiShowResults(eventSet,cntinue);
	}while (cntinue);
        papiUnsubscribeEvents(eventSet); 
        
        logFile.close();
}

