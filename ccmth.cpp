#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <thread>
#include <future>
#include <chrono>
#include "ccmsrc.h"

using namespace std;

// Create a struct to pass to the threads
struct CCM_asyncData{

    int iEmbeddingDimension,
        iLagTimeStep,
        iLibraryLength,
        iMaxAutoCorrToCheck,
        iNumberOfTimeStep2Check,
        iEstimatedYLength,
	iYEmbeddingDimension,
	iYLagTimeStep;

    double dCompareTolerance,
           dWeightToleranceLevel,
	   dCCMcorr;

    bool bL2;

    double *dX,*dY,*dEstimatedY;

};

// Simple functions to help stay organized
void CmdLineHelp();
bool ReadCmdLineArgs(int,char**);
double CCMcorr_async(CCM_asyncData);

// Globals read from command line
int iEmbeddingDimension = -1,
    iLagTimeStep = -1,
    iLibraryLength = -1,
    iNumOfTimeSeries = -1,
    iNumThreads = 0,
    iThreadWaitSec = -1,
    iMaxAutoCorrToCheck = -1,
    iNumberOfTimeStep2Check = -1,
    iYEmbeddingDimension = -1,
    iYLagTimeStep = -1;

double dCompareTolerance = -1,
       dWeightToleranceLevel = -1;

char *cOutputname,
     *cTSOutputname,
     *cFilename; //this data file is assumed to be in the proper format (see "usage")

bool bVerboseFlag = false,
     bOptThreading = false,
     bSaveEstimatedTS = false,
     bFindLag = false,
     bFindE = false,
     bPAI = false,
     bPAICalcPart1 = false,
     bPAICalcPart2 = false,
     bL2 = false;

// Main
int main(int argc, char **argv){

    bool bCmdLineOK = ReadCmdLineArgs(argc,argv);
    if( !bCmdLineOK ){
        return(-1);
    }

    char cBuffer[1024]; //read in buffer
    char *cColValue, //used for read in
         *cRowValues; //used for read in
    int iter = 0;

    int iEstimatedY_Length = iLibraryLength-((iEmbeddingDimension-1)*iLagTimeStep);
    double dX[iLibraryLength],
           dY[iLibraryLength],
           dEstimatedY[iEstimatedY_Length];

    //Populate estimate time series, just in case
    for(int iPre = 0;iPre < iEstimatedY_Length;iPre++ ){
        dEstimatedY[iPre] = nan("");
    }

    //Declare things needed for threading
    double dCCMcorrs[iNumThreads][4];
    future<double> WorkerThreads[iNumThreads][4];
    future_status WorkerStatus[4];
    CCM_asyncData sThreadData;
    sThreadData.iEmbeddingDimension = iEmbeddingDimension;
    sThreadData.iLagTimeStep = iLagTimeStep;
    sThreadData.iYEmbeddingDimension = iYEmbeddingDimension;
    sThreadData.iYLagTimeStep = iYLagTimeStep;
    sThreadData.iLibraryLength = iLibraryLength;
    sThreadData.iMaxAutoCorrToCheck = iMaxAutoCorrToCheck;
    sThreadData.dCompareTolerance = dCompareTolerance;
    sThreadData.iNumberOfTimeStep2Check = iNumberOfTimeStep2Check;
    sThreadData.dWeightToleranceLevel = dWeightToleranceLevel;
    sThreadData.iEstimatedYLength = iEstimatedY_Length;
    sThreadData.dEstimatedY = dEstimatedY;
    sThreadData.bL2 = bL2;

    if( bVerboseFlag ){
            printf("Input data file: %s\nOutput file: %s\nNumber of time series in input file: %i\nLibrary length: %i\nEmbedding dimension: %i\nLag time step: %i\n",cFilename,cOutputname,iNumOfTimeSeries,iLibraryLength,iEmbeddingDimension,iLagTimeStep);
    }

    //open data input file
    FILE *ifstream = fopen(cFilename,"r");
    if(ifstream == NULL){
        fprintf(stderr, "Error: Cannot open %s\n",cFilename);
            return(-1);
    }

    //open [-o] output file
    FILE *ofstream = fopen(cOutputname,"w");
    if(ofstream == NULL){
        fprintf(stderr, "Error: Cannot open %s\n",cOutputname);
        return(-1);
    }

    //open [-eY] output file
    FILE *ofstreamTS = fopen(cTSOutputname,"w");
    if(ofstream == NULL){
        fprintf(stderr, "Error: Cannot open %s\n",cTSOutputname);
        return(-1);
    }

    if( bVerboseFlag ){ printf("Processing %s...\n",cFilename); }
    int iCurrentThreadsUsed = 0,
        iCheckThread = 0;
    bool bFoundOpenThread;
    for(int iTSiter = 0;iTSiter < iNumOfTimeSeries;iTSiter++ ){

        if( bVerboseFlag ){ printf("TS %i:",iTSiter); }
        if( bVerboseFlag ){ printf("Reading TS..."); }
        bFoundOpenThread = false;

        //read in dX and dY
        for( iter=0;iter < iLibraryLength;iter++ ){
                cRowValues = fgets(cBuffer,sizeof(cBuffer),ifstream);
                cColValue = strtok(cRowValues,",");
                if( cColValue != NULL ){
                   dX[iter] = atof(cColValue);
                }else{
                   fprintf(stderr, "Error: %s is not formatted correctly\n",cFilename);
                }

                cColValue = strtok(NULL,";");
                if( cColValue != NULL ){
                   dY[iter] = atof(cColValue);
                }else{
                   fprintf(stderr, "Error: %s is not formatted correctly\n",cFilename);
                }
        }

        if( bVerboseFlag ){ printf("Writing CCMs..."); }

        if( iNumThreads == 0 ){

              if( bPAI ){

                  sThreadData.dX = dX;
                  sThreadData.dY = dY;
                  bPAICalcPart1 = true;
                  fprintf(ofstream,"%.20f,",CCMcorr_async(sThreadData));
                  //printf("%.20f,",CCMcorr_async(sThreadData));
                  if( bSaveEstimatedTS ){
                      for(int iPrintStep = 0;iPrintStep < sThreadData.iEstimatedYLength;iPrintStep++ ){
			fprintf(ofstreamTS,"%.20f\n",sThreadData.dEstimatedY[iPrintStep]);
                      }
                  }
                  bPAICalcPart1 = false;
                  bPAICalcPart2 = true;
                  fprintf(ofstream,"%.20f,",CCMcorr_async(sThreadData));
                  //printf("%.20f,",CCMcorr_async(sThreadData));
                  if( bSaveEstimatedTS ){
                      for(int iPrintStep = 0;iPrintStep < sThreadData.iEstimatedYLength;iPrintStep++ ){
			fprintf(ofstreamTS,"%.20f\n",sThreadData.dEstimatedY[iPrintStep]);
                      }
                  }
                  bPAICalcPart2 = false;

                  bPAICalcPart1 = true;
                  sThreadData.dX = dY;
                  sThreadData.dY = dX;
                  fprintf(ofstream,"%.20f,",CCMcorr_async(sThreadData));
                  //printf("%.20f\n",CCMcorr_async(sThreadData));
                  if( bSaveEstimatedTS ){
                      for(int iPrintStep = 0;iPrintStep < sThreadData.iEstimatedYLength;iPrintStep++ ){
			fprintf(ofstreamTS,"%.20f\n",sThreadData.dEstimatedY[iPrintStep]);
                      }
                  }
                  bPAICalcPart1 = false;
                  bPAICalcPart2 = true;
                  fprintf(ofstream,"%.20f\n",CCMcorr_async(sThreadData));
                  //printf("%.20f\n",CCMcorr_async(sThreadData));
                  if( bSaveEstimatedTS ){
                      for(int iPrintStep = 0;iPrintStep < sThreadData.iEstimatedYLength;iPrintStep++ ){
			fprintf(ofstreamTS,"%.20f\n",sThreadData.dEstimatedY[iPrintStep]);
                      }
                  }
                  bPAICalcPart2 = false;

              }else{

                  sThreadData.dX = dX;
                  sThreadData.dY = dY;
                  fprintf(ofstream,"%.20f,",CCMcorr_async(sThreadData));
                  //printf("%.20f,",CCMcorr_async(sThreadData));
                  if( bSaveEstimatedTS ){
                      for(int iPrintStep = 0;iPrintStep < sThreadData.iEstimatedYLength;iPrintStep++ ){
			fprintf(ofstreamTS,"%.20f\n",sThreadData.dEstimatedY[iPrintStep]);
                      }
                  }

                  sThreadData.dX = dY;
                  sThreadData.dY = dX;
                  fprintf(ofstream,"%.20f\n",CCMcorr_async(sThreadData));
                  //printf("%.20f\n",CCMcorr_async(sThreadData));
                  if( bSaveEstimatedTS ){
                      for(int iPrintStep = 0;iPrintStep < sThreadData.iEstimatedYLength;iPrintStep++ ){
			fprintf(ofstreamTS,"%.20f\n",sThreadData.dEstimatedY[iPrintStep]);
                      }
                  }

              }

        }else{
            if( iCurrentThreadsUsed < iNumThreads ){

                sThreadData.dX = dX;
                sThreadData.dY = dY;
                WorkerThreads[iCurrentThreadsUsed][0] = async(launch::async,CCMcorr_async,sThreadData);

                sThreadData.dX = dY;
                sThreadData.dY = dX;
                WorkerThreads[iCurrentThreadsUsed][1] = async(launch::async,CCMcorr_async,sThreadData);

                iCurrentThreadsUsed++;

            }else{

                if( bOptThreading ){

                    while( !bFoundOpenThread ){
                        if( WorkerThreads[iCheckThread][0].valid() && WorkerThreads[iCheckThread][1].valid() ){

                            WorkerStatus[0] = WorkerThreads[iCheckThread][0].wait_for(chrono::seconds(iThreadWaitSec));
                            WorkerStatus[1] = WorkerThreads[iCheckThread][1].wait_for(chrono::seconds(iThreadWaitSec));

                            if( (WorkerStatus[0] == future_status::ready) && (WorkerStatus[1] == future_status::ready) ){
				
				fprintf(ofstream,"%.20f,",WorkerThreads[iCheckThread][0].get());
				fprintf(ofstream,"%.20f\n",WorkerThreads[iCheckThread][1].get());
				//printf("%.20f,",WorkerThreads[iCheckThread][0].get());
				//printf("%.20f\n",WorkerThreads[iCheckThread][1].get());

                                sThreadData.dX = dX;
                                sThreadData.dY = dY;
                                WorkerThreads[iCheckThread][0] = async(launch::async,CCMcorr_async,sThreadData);

                                sThreadData.dX = dY;
                                sThreadData.dY = dX;
                                WorkerThreads[iCheckThread][1] = async(launch::async,CCMcorr_async,sThreadData);

                                bFoundOpenThread = true;

                            }

                        }else if(!WorkerThreads[iCheckThread][0].valid()&&!WorkerThreads[iCheckThread][1].valid() ){

                            sThreadData.dX = dX;
                            sThreadData.dY = dY;
                            WorkerThreads[iCheckThread][0] = async(launch::async,CCMcorr_async,sThreadData);

                            sThreadData.dX = dY;
                            sThreadData.dY = dX;
                            WorkerThreads[iCheckThread][1] = async(launch::async,CCMcorr_async,sThreadData);

                            bFoundOpenThread = true;

                        }else{

                            iCheckThread++;
                            if( iCheckThread > iNumThreads ){
                                iCheckThread = 0;
                            }

                        }

                    }

                }else{

                    for(int iThreadIter = 0;iThreadIter < iNumThreads;iThreadIter++ ){

                        if( WorkerThreads[iThreadIter][0].valid() && WorkerThreads[iThreadIter][1].valid() ){
                            dCCMcorrs[iThreadIter][0] = WorkerThreads[iThreadIter][0].get();
                            dCCMcorrs[iThreadIter][1] = WorkerThreads[iThreadIter][1].get();
                        }

                    }

                    iCurrentThreadsUsed = 0;

                    sThreadData.dX = dX;
                    sThreadData.dY = dY;
                    WorkerThreads[iCurrentThreadsUsed][0] = async(launch::async,CCMcorr_async,sThreadData);

                    sThreadData.dX = dY;
                    sThreadData.dY = dX;
                    WorkerThreads[iCurrentThreadsUsed][1] = async(launch::async,CCMcorr_async,sThreadData);

                    for(int iThreadIter = 0;iThreadIter < iNumThreads;iThreadIter++ ){
                        fprintf(ofstream,"%.20f,",dCCMcorrs[iThreadIter][0]);
                        fprintf(ofstream,"%.20f\n",dCCMcorrs[iThreadIter][1]);
                        //printf("%.20f,",dCCMcorrs[iThreadIter][0]);
                        //printf("%.20f\n",dCCMcorrs[iThreadIter][1]);
                        dCCMcorrs[iThreadIter][0] = nan("");
                        dCCMcorrs[iThreadIter][1] = nan("");
                    }

                    iCurrentThreadsUsed++;
                }
            }
        }

        if( bVerboseFlag ){ printf(" done.\n"); }

    }

    for(int iThreadIter = 0;iThreadIter < iNumThreads;iThreadIter++ ){

        if( WorkerThreads[iThreadIter][0].valid() && WorkerThreads[iThreadIter][1].valid() ){
            fprintf(ofstream,"%.20f,",WorkerThreads[iThreadIter][0].get());
            fprintf(ofstream,"%.20f\n",WorkerThreads[iThreadIter][1].get());

        }

    }

    fclose(ifstream);
    fclose(ofstream);
    fclose(ofstreamTS);

    if( bVerboseFlag ){ printf(" Finished.\n"); }

    return(0);
}

void CmdLineHelp(char* cName){

    printf("usage: %s <command line flags>\n\n",cName);
    printf("basic command line flags:\n");
    printf("  -E [integer]  : embedding dimension (ignored if using -FindLag or -FindE)\n");
    printf("  -t [integer]  : lag time step (ignored if using -FindLag)\n");
    printf("  -L [integer]  : library length of the time series\n");
    printf("                 (assumed to be equal for X and Y)\n");
    printf("  -f [string]   : filename of a text file containing time series\n");
    printf("                 data for X and Y with columns seperated by\n");
    printf("                 commas and rows seperated by semicolons \n");
    printf("                 (it is assumed that there are only X and Y\n");
    printf("                 in the file, i.e. only two columns)\n");
    printf("  -n [integer]  : number of time series in the -f file\n");
    printf("  -p [integer]  : number of threads to use (default is zero)\n");
    printf("  -o [string]   : filename of the output text file\n");
    printf("  -Op [integer] : 'optimize' threading check time\n");
    printf("  -v            : flag that sends various bits of information to STDOUT\n");
    printf("  -eY [string]  : filename of a text file containing the estimated time series\n");
    printf("  -PAI          : this flag calculates the pairwise asymmetric inference\n");
    printf("                      (see new output format below)\n");
    printf("  -L2           : this flag calculates the L2 norm (Euclidean norm)\n");
    printf("                      instead of the CCM correlation\n");
    printf("  -Ey           : embedding dimension for Y (ignored if not using -PAI)\n");
    printf("  -ty           : lag time step for Y (ignored if not using -PAI)\n");
    printf("Flag for finding lag time to use for embedding:\n");
    printf("  -FindLag [integer1] [integer2]\n");
    printf("     [integer1] : maximum autocorrelation lag time to check\n");
    printf("     [integer2] : comparison tolerance for autocorrelation checks\n");
    printf("  NOTE: This flag will cause the [-o] output file contents to be the \n");
    printf("        lag time step leading to the maximum autocorrelation\n");
    printf("        (within the user defined bounds)\n");
    printf("Flag for finding embedding dimension to use for embedding:\n");
    printf("  -FindE [integer1] [integer2]\n");
    printf("     [integer1] : maximum number of time steps to check\n");
    printf("     [integer2] : cutoff tolerance for weight checks\n");
    printf("  NOTE: This flag will cause the [-o] output file contents to be the \n");
    printf("        mode of cutoff embedding dimensions for requested time steps\n");
    printf("        (within the user defined bounds)\n");
    printf("input data file format:\n");
    printf("   X0,Y0;\n");
    printf("   X1,Y1;\n");
    printf("   X2,Y2;\n");
    printf("   X3,Y3;\n");
    printf("   X4,Y4;\n");
    printf("   ...\n");
    printf("[-o] output data file format:\n");
    printf("   CCM(X,Y),CCM(Y,X) [time series 1]\n");
    printf("   CCM(X,Y),CCM(Y,X) [time series 2]\n");
    printf("   ...\n");
    printf("[-o] output data file format (w/ PAI flag):\n");
    printf("   CCM(X,X),CCM(X,X and Y),CCM(Y,Y),CCM(Y,Y and X) [time series 1]\n");
    printf("   CCM(X,X),CCM(X,X and Y),CCM(Y,Y),CCM(Y,Y and X) [time series 2]\n");
    printf("   ...\n");
    printf(" NOTE: The contents of the [-o] output file depends on the flags used.\n");
    printf("[-eY] output data file format:\n");
    printf("   estimated X0\n");
    printf("   estimated X1\n");
    printf("   estimated X2\n");
    printf("   estimated X3\n");
    printf("   estimated X4\n");
    printf("   ...\n");
    printf("   estimated Y0\n");
    printf("   estimated Y1\n");
    printf("   estimated Y2\n");
    printf("   estimated Y3\n");
    printf("   estimated Y4\n");
    printf("   ...\n");
    printf("[-eY] output data file format (w/ PAI flag):\n");
    printf("   estimated X0 (given X)\n");
    printf("   estimated X1 (given X)\n");
    printf("   estimated X2 (given X)\n");
    printf("   estimated X3 (given X)\n");
    printf("   estimated X4 (given X)\n");
    printf("   ...\n");
    printf("   estimated X0 (given X and Y)\n");
    printf("   estimated X1 (given X and Y)\n");
    printf("   estimated X2 (given X and Y)\n");
    printf("   estimated X3 (given X and Y)\n");
    printf("   estimated X4 (given X and Y)\n");
    printf("   ...\n");
    printf("   estimated Y0 (given Y)\n");
    printf("   estimated Y1 (given Y)\n");
    printf("   estimated Y2 (given Y)\n");
    printf("   estimated Y3 (given Y)\n");
    printf("   estimated Y4 (given Y)\n");
    printf("   ...\n");
    printf("   estimated Y0 (given Y and X)\n");
    printf("   estimated Y1 (given Y and X)\n");
    printf("   estimated Y2 (given Y and X)\n");
    printf("   estimated Y3 (given Y and X)\n");
    printf("   estimated Y4 (given Y and X)\n");
    printf("   ...\n");
    printf("\nMore details can be found at https://github.com/<githubaccount>/...\n");

}

bool ReadCmdLineArgs(int argc, char **argv){

    bool bInputSet = false,
         bOutputSet = false;

    for(int iter=1;iter < argc;iter++ ){

        if( (strcmp("-h",argv[iter]) == 0) || (strcmp("-?",argv[iter]) == 0) ){
            CmdLineHelp(argv[0]);
            return( false );

        }else if( strcmp("-E",argv[iter]) == 0 ){

            iEmbeddingDimension = atoi(argv[iter+1]);

        }else if( strcmp("-t",argv[iter]) == 0 ){

            iLagTimeStep = atoi(argv[iter+1]);

        }else if( strcmp("-Ey",argv[iter]) == 0 ){

            iYEmbeddingDimension = atoi(argv[iter+1]);

        }else if( strcmp("-ty",argv[iter]) == 0 ){

            iYLagTimeStep = atoi(argv[iter+1]);

        }else if( strcmp("-L",argv[iter]) == 0 ){

            iLibraryLength = atoi(argv[iter+1]);

        }else if( strcmp("-f",argv[iter]) == 0 ){

            cFilename = argv[iter+1];
            bInputSet = true;

        }else if( strcmp("-n",argv[iter]) == 0 ){

            iNumOfTimeSeries = atoi(argv[iter+1]);

        }else if( strcmp("-p",argv[iter]) == 0 ){

            iNumThreads = atoi(argv[iter+1]);

        }else if( strcmp("-o",argv[iter]) == 0 ){

            cOutputname = argv[iter+1];
            bOutputSet = true;

        }else if( strcmp("-eY",argv[iter]) == 0 ){

            cTSOutputname = argv[iter+1];
            bSaveEstimatedTS = true;

        }else if( strcmp("-v",argv[iter]) == 0 ){

            bVerboseFlag = true;

        }else if( strcmp("-PAI",argv[iter]) == 0 ){

            bPAI = true;

        }else if( strcmp("-L2",argv[iter]) == 0 ){

            bL2 = true;

        }else if( strcmp("-Op",argv[iter]) == 0 ){

            bOptThreading = true;
            iThreadWaitSec = atoi(argv[iter+1]);

        }else if( strcmp("-FindLag",argv[iter]) == 0 ){

            bFindLag = true;
            iMaxAutoCorrToCheck = atoi(argv[iter+1]);
            dCompareTolerance = atof(argv[iter+2]);
        }else if( strcmp("-FindE",argv[iter]) == 0 ){

            bFindE = true;
            iNumberOfTimeStep2Check = atoi(argv[iter+1]);
            dWeightToleranceLevel = atof(argv[iter+2]);

        }

    }

    //Make sure the required things were defined
    if( iEmbeddingDimension < 0 && !bFindLag ){
        fprintf(stderr, "Error: embedding dimension is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( iLagTimeStep < 0 && !bFindLag ){
        fprintf(stderr, "Error: lag time step is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( iLibraryLength < 0 ){
        fprintf(stderr, "Error: library length is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( !bInputSet ){
        fprintf(stderr, "Error: input data filename is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( iNumOfTimeSeries < 0 ){
        fprintf(stderr, "Error: number of time series in input data file not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( !bOutputSet ){
        fprintf(stderr, "Error: output data filename is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( bOptThreading && (iThreadWaitSec < 0) ){
        fprintf(stderr, "Error: Thread wait time must be defined as an integer > 0\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }
    if( !bSaveEstimatedTS ){
        fprintf(stderr, "Error: output data filename for estimated time series not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindLag && (iMaxAutoCorrToCheck < 0) ){
        fprintf(stderr, "Error: maximum autocorrelation lag time to check is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindLag && (dCompareTolerance < 0) ){
        fprintf(stderr, "Error: autocorrelation comparison tolerance is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindE && (iNumberOfTimeStep2Check < 0) ){
        fprintf(stderr, "Error: number of time steps to check is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindE && (dWeightToleranceLevel < 0) ){
        fprintf(stderr, "Error: weight tolerance is not defined\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindE && bFindLag ){
        fprintf(stderr, "Error: -FindLag and -FindE cannot be used together\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindE && bPAI ){
        fprintf(stderr, "Error: -FindE and -PAI cannot be used together\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( bFindLag && bPAI ){
        fprintf(stderr, "Error: -FindLag and -PAI cannot be used together\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }if( (iNumThreads != 0) && bSaveEstimatedTS ){
        fprintf(stderr, "Error: -eY cannot be used with threading\n\n");
        CmdLineHelp(argv[0]);
        return( false );
    }

    //Issue any warnings
    if( iNumThreads < 0 ){
        fprintf(stderr, "Warning: number of threads improperly defined; using default\n");
        iNumThreads = 0;
    }
    if( bVerboseFlag && (iNumThreads != 0) ){
        fprintf(stderr, "Warning: verbose logging is not available with multiple threads; turning off verbose output\n");
        bVerboseFlag = false;
    }


    return( true );
}

double CCMcorr_async(CCM_asyncData sThreadData){

    double dResult;
    int iResult;

    if( bFindLag ){
        iResult = FindLagTimeStep(sThreadData.dX,sThreadData.iLibraryLength,sThreadData.iMaxAutoCorrToCheck,sThreadData.dCompareTolerance,bVerboseFlag);
        return( (double) iResult );
    }else if( bFindE ){
        iResult = FindEmbeddingDimension(sThreadData.dX,sThreadData.iLibraryLength,sThreadData.iNumberOfTimeStep2Check,sThreadData.dWeightToleranceLevel,sThreadData.iLagTimeStep,bVerboseFlag);
        return( (double) iResult );
    }else if( bPAICalcPart1 ){
        CCMcorr(sThreadData.dCCMcorr,sThreadData.dX,sThreadData.iLibraryLength,sThreadData.dX,sThreadData.iLibraryLength,sThreadData.iEmbeddingDimension,sThreadData.iLagTimeStep,sThreadData.dEstimatedY,sThreadData.iEstimatedYLength,sThreadData.bL2);
	dResult = sThreadData.dCCMcorr;
        return( dResult );
    }else if( bPAICalcPart2 ){
        CCMcorr2(sThreadData.dCCMcorr,sThreadData.dX,sThreadData.iLibraryLength,sThreadData.dX,sThreadData.iLibraryLength,sThreadData.dY,sThreadData.iLibraryLength,sThreadData.iEmbeddingDimension,sThreadData.iLagTimeStep,sThreadData.iYEmbeddingDimension,sThreadData.iYLagTimeStep,sThreadData.dEstimatedY,sThreadData.iEstimatedYLength,sThreadData.bL2);
	dResult = sThreadData.dCCMcorr;
        return( dResult );
    }else{
	CCMcorr(sThreadData.dCCMcorr,sThreadData.dX,sThreadData.iLibraryLength,sThreadData.dY,sThreadData.iLibraryLength,sThreadData.iEmbeddingDimension,sThreadData.iLagTimeStep,sThreadData.dEstimatedY,sThreadData.iEstimatedYLength,sThreadData.bL2);
	dResult = sThreadData.dCCMcorr;
        return( dResult );
    }

}
