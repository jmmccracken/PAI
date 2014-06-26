#include "ccmsrc.h"

/*FindLagTimeStep -- Find a Lag Time Step
  This function calculates autocorrelations for sequential lag times
  starting with 1 and going to a specified value, and then returns the
  lag time used to calculate the maximum value in that set.  Warnings
  are returned if the autocorrelations appear to have cyclic behavior
  because it is assumed that the user wants to know.

  References: http://en.wikipedia.org/wiki/Autocorrelation

  Dependencies:   Pcorr (function call)
                  math.h (header)

  Input Parameters:
      double dX[] - vector of values
      int iX_length - length of X
      int iMaxAutoCorrToCheck - maximum lag used to calculate an autocorrelation
      double dCompareTolerance - the tolerance value below which two values should
                                 be considered "equal"
      bool bVerboseFlag - flag that controls printing the autocorrelation values to
                          STDERR

  Output Parameters:
      int - lag time used to calculate the maximum autocorrelation within the given bounds
*/
int FindLagTimeStep(double dX[],int iXlength,int iMaxAutoCorrToCheck,double dCompareTolerance,bool bVerboseFlag){

    if( iMaxAutoCorrToCheck > (iXlength/2) ){
        fprintf(stderr, "Warning in FindLagTimeStep(): the maximum autocorreltaion lag is half the library length; ignoring requested value and checking to a lag of half the library length\n");
        iMaxAutoCorrToCheck = (iXlength/2);
    }

    double dAutoCorrs[iMaxAutoCorrToCheck],
           dCurrentAutoCorr = 0;
    int iTempVecLength = 0;

    //Calculate autocorrelations
    if( bVerboseFlag ){
            fprintf(stderr,"-- FindLagTimeStep() verbose output begin --\n");
    }
    for(int iAutoCorrLagStep = 1;iAutoCorrLagStep <= iMaxAutoCorrToCheck;iAutoCorrLagStep++ ){

        iTempVecLength = iXlength-iAutoCorrLagStep;
        double dXtemp1[iTempVecLength],dXtemp2[iTempVecLength];

        for(int iTempVecBuildStep=0;iTempVecBuildStep < iTempVecLength;iTempVecBuildStep++ ){
            dXtemp1[iTempVecBuildStep] = dX[iTempVecBuildStep];
            dXtemp2[iTempVecBuildStep] = dX[iTempVecBuildStep+iAutoCorrLagStep];
        }

        dCurrentAutoCorr = Pcorr(dXtemp1,dXtemp2,iTempVecLength);
        dAutoCorrs[iAutoCorrLagStep-1] = dCurrentAutoCorr*dCurrentAutoCorr;
        if( bVerboseFlag ){
            fprintf(stderr,"%.20f\n",dAutoCorrs[iAutoCorrLagStep-1]);
        }
    }
    if( bVerboseFlag ){
            fprintf(stderr,"-- FindLagTimeStep() verbose output end --\n");
    }

    //Find max
    int iCurrentMaxAutoCorrLag = 1;
    dCurrentAutoCorr = dAutoCorrs[0];
    for(int iter = 1;iter < iMaxAutoCorrToCheck;iter++ ){

        if( fabs(dAutoCorrs[iter]-dCurrentAutoCorr) < dCompareTolerance ){
            fprintf(stderr, "Warning in FindLagTimeStep(): autocorrelations may be cyclic; run with verbose flag to see all the autocorrelations\n");
        }

    if( fabs(dCurrentAutoCorr-1) < dCompareTolerance ){
            fprintf(stderr, "Warning in FindLagTimeStep(): autocorrelations may be near unity; run with verbose flag to see all the autocorrelations\n");
        }

        if( dAutoCorrs[iter] > dCurrentAutoCorr ){
            dCurrentAutoCorr = dAutoCorrs[iter];
            iCurrentMaxAutoCorrLag = iter+1;
        }
    }

    return(iCurrentMaxAutoCorrLag);
}

/*FindEmbeddingDimension -- Find the max embedding dimension
  This function finds the embedding dimension for which the weight of
  a nearest neighbor falls below a given threshold value.  Such a "cutoff"
  embedding dimension is found for a given number of time steps and then
  the mode of those values is returned.

  References: see the directed correlation paper and notes

  Dependencies: FindWeightsFromShadow (function call)
                iMode (function call)

  Input Parameters:
        double dX[] - vector of data used to create the shadow manifold
        int iXlength - length of X
        int iNumberOfTimeStep2Check - the number of time steps to check
        double dWeightToleranceLevel - the tolerance level used to define the
                                       "cutoff" embedding dimension
        int iLagTime - the lag time step used to create the shadow manifold
        bool bVerboseFlag - flag that controls the printing of the weights to
                            STDERR

  Output Parameters:
        int - the mode of all the requested time steps' "cutoff" embedding
              dimensions
*/
int FindEmbeddingDimension(double dX[],int iXlength,int iNumberOfTimeStep2Check,double dWeightToleranceLevel,int iLagTime,bool bVerboseFlag){

    //Set the max embedding dimension to check (to save computation time)
    int iMaxEDim = (iXlength/2);

    //Find minimum size of possible shadow manifolds and report a warning
    //if the requested number of time step checks is too large
    int iCalShadManDimMin = iXlength-((iMaxEDim-1)*iLagTime);
    if( iNumberOfTimeStep2Check > iCalShadManDimMin ){
        fprintf(stderr,"Warning in FindEmbeddingDimension(): the requested number of time step checks is too large; setting the value to %i\n",iCalShadManDimMin);
        iNumberOfTimeStep2Check = iCalShadManDimMin;
    }

    int iCutoffEDimForEachTimeStep[iNumberOfTimeStep2Check];
    for(int iter=0;iter < iNumberOfTimeStep2Check;iter++ ){
        iCutoffEDimForEachTimeStep[iter] = iMaxEDim;
    }

    if( bVerboseFlag ){
        fprintf(stderr,"-- FindEmbeddingDimension() verbose output begin --\n");
    }

    //Loop through the requested number of delay vectors to find the weights
    for(int iDelayVectorN = 1; iDelayVectorN <= iNumberOfTimeStep2Check;iDelayVectorN++ ){

        for(int iEstep = 2;iEstep < iMaxEDim;iEstep++ ){

            //Find calculated shadow manifold dimension
            int iCalShadManDim = iXlength-((iEstep-1)*iLagTime);

            //Assign some memory for the shadow manifold
            double** dXShadow;
            dXShadow = new double*[iCalShadManDim];
            for(int iter = 0;iter < iCalShadManDim;iter++ ){
                dXShadow[iter] = new double[iEstep];
            }
            int iTstep4delayvector;

            //Create the shadow manifold by creating the delay vectors sequentially
            for(int iShadowStep = 1; iShadowStep <= iCalShadManDim;iShadowStep++ ){
                iTstep4delayvector = iShadowStep+((iEstep-1)*iLagTime);
                for(int iDimStep = 1; iDimStep <= iEstep; iDimStep++ ){
                    dXShadow[iShadowStep-1][iDimStep-1] = dX[(iTstep4delayvector-((iDimStep-1)*iLagTime))-1];
                }
            }

            //Storage for sorting, calculating, and organizing
            double dDelayVectorOfInterest[iEstep],  //delay vector currently being compared
                   dWeights[(iEstep+1)],  //weights to contruct Y estimate
                   dYEstimateGivenX[iCalShadManDim];  //estimated Y used the shadow manifold of X

            int iTstepOfNearestNeighborsTempRow[iEstep+1];  //time steps of nearest neighbor delay vectors

            //Populate the temp delay vector
            for(int iCopyStep = 0; iCopyStep < iEstep; iCopyStep++ ){
                dDelayVectorOfInterest[iCopyStep] = dXShadow[iDelayVectorN-1][iCopyStep];
            }

            FindWeightsFromShadow(dWeights,iTstepOfNearestNeighborsTempRow,dDelayVectorOfInterest,dXShadow,iCalShadManDim,iEstep,iLagTime,iEstep);

            //Free the dXshadow memory
            for(int iter = 0;iter < iCalShadManDim;iter++ ){
                delete dXShadow[iter];
            }
            delete dXShadow;

            //Check weights
            for(int iWstep = 0;iWstep < (iEstep+1);iWstep++ ){
                if( dWeights[iWstep] < dWeightToleranceLevel ){
                    iCutoffEDimForEachTimeStep[iDelayVectorN] = iEstep;
                    iEstep = (iXlength/2)+1; //break out of the loop over embedding dimensions
                }
                if( bVerboseFlag ){
                    fprintf(stderr,"t = %i;E = %i; weight %i = %.20f (tol = %.20f)\n",iDelayVectorN,iEstep,iWstep,dWeights[iWstep],dWeightToleranceLevel);
                }
            }
        }
    }

   if( bVerboseFlag ){
        fprintf(stderr,"Cutoff embedding dimensions for each time step:\n");
        for(int iter=0;iter < iNumberOfTimeStep2Check;iter++ ){
            fprintf(stderr,"%i\n",iCutoffEDimForEachTimeStep[iter]);
        }
   }
   if( bVerboseFlag ){
       fprintf(stderr,"-- FindEmbeddingDimension() verbose output end --\n");
   }

   bool bCheck = true;
   for(int iter=0;iter < iNumberOfTimeStep2Check;iter++ ){
        if( iCutoffEDimForEachTimeStep[iter] == iMaxEDim){
                bCheck &= true;
        }
   }
    if( bCheck ){
        fprintf(stderr,"Warning in FindEmbeddingDimension(): cutoff dimension may have maxed out\n");
    }

    //Find and return mode
    int iEDimMode = iMode(iCutoffEDimForEachTimeStep,iNumberOfTimeStep2Check);
    return(iEDimMode);
}


/*FindWeightsFromShadow -- Find the Normalized Nearest Neighbor Weights
  This function calculates the normalized weights from the nearest
  neighbors of a given delay vector on the given shadow manifold.  This
  function is designed to work seamlessly with CCMCorr, so the variable
  names might be a little strange.

  References: http://en.wikipedia.org/wiki/Convergent_cross_mapping
              https://www.sciencemag.org/content/338/6106/496.figures-only
              https://www.sciencemag.org/content/338/6106/496/suppl/DC1

  Dependencies:   eDist (function call)
                  FindAndCheckIndex (function call)
                  math.h (header)

  Input Parameters:
      double dWeights[] - array of (iEmbeddingDimension+1) doubles into
                          which the calculated weights will be placed
      int iTstepOfNearestNeighborsTempRow[] - array of (iEmbeddingDimension+1) integers
                                              into which the time steps associated to
                                              the weights in dWeights will be placed
      double dDelayVectorOfInterest[] - delay vector used to find the nearest neighbors
      double** dXShadow - the shadow manifold
      int iCalShadManDim - the calculated shadow manifold dimension (e.g. calculated
                           in CCMCorr)
      int iEmbeddingDimension - dimension of the shadow manifold
      int iLagTime - lag time step used to construct the delay vectors
                    that make up the n-dimensional points of the
                    shadow manifold

  Output Parameters:

*/
void FindWeightsFromShadow(double dWeights[],  int iTstepOfNearestNeighborsTempRow[],  double dDelayVectorOfInterest[],double** dXShadow,int iCalShadManDim, int iEmbeddingDimension,int iLagTime,int iDelayVectorOfInterest_length){

    //Storage for sorting, calculating, and organizing
    double dDelayVector2Compare[iDelayVectorOfInterest_length],  //iterated delay vector for comparison
           dXShadowNorm,  //eucleadian distance between delay vectors being compared
           dUnsortedNorms[iCalShadManDim],  //unsorted eucleadian distances between delay vectors
           dSortedNorms[iCalShadManDim],  //sorted eucleadian distances between delay vectors
           dWeightDenominator,  //denominator in weights calculation (used for zero check)
           dWeightNormalization;  //normalization factor for weights
    int iTempValue;  //temporary value for time step of nearest neightbor delay vector

    //Loop through all of the delay vector to compare to the temp delay vector
    for(int iTstep = 1; iTstep <= iCalShadManDim; iTstep++ ){

        //Populate the temp delay vector for comparison
        for(int iCopyStep = 0; iCopyStep < iDelayVectorOfInterest_length; iCopyStep++ ){
            dDelayVector2Compare[iCopyStep] = dXShadow[iTstep-1][iCopyStep];
        }

        //Find the distance between the two temp vectors and save it
        dXShadowNorm = eDist(dDelayVectorOfInterest,dDelayVector2Compare,iDelayVectorOfInterest_length);

        //Save the norms to be sorted later
        dUnsortedNorms[iTstep-1] = dXShadowNorm;
        dSortedNorms[iTstep-1] = dXShadowNorm;
    }

    //Sort the norms
    std::sort(dSortedNorms,dSortedNorms+iCalShadManDim);

    //Save the time step locations of the E+1 closest neighbors
    iTempValue = 0;
    for(int iEDimStep = 0; iEDimStep < (iEmbeddingDimension+1); iEDimStep++ ){
        iTstepOfNearestNeighborsTempRow[iEDimStep] = -1;
    }
    for(int iEDimStep = 0; iEDimStep < (iEmbeddingDimension+1); iEDimStep++ ){
        iTempValue = FindAndCheckIndex(dUnsortedNorms,iCalShadManDim,dSortedNorms[iEDimStep+1],iTstepOfNearestNeighborsTempRow,(iEmbeddingDimension+1));
        iTstepOfNearestNeighborsTempRow[iEDimStep] = iTempValue;
    }

    //Find denominator for weight calculation
    if( dSortedNorms[1] == 0 ){
        dWeightDenominator = nan("");
        fprintf(stderr, "Warning in CCMcorr(): division by zero\n");
    }else{
        dWeightDenominator = dSortedNorms[1];
    }

    //Find weights
    dWeightNormalization = 0;
    for( int iWeightStep = 0;iWeightStep < (iEmbeddingDimension+1);iWeightStep++ ){
        dWeights[iWeightStep] = exp((-1*dSortedNorms[iWeightStep+1])/dWeightDenominator);
        dWeightNormalization += dWeights[iWeightStep];
    }

    //Check normalization factor just to be safe
    if( dWeightNormalization == 0 ){
        dWeightNormalization = nan("");
        fprintf(stderr, "Warning in CCMcorr(): division by zero\n");
    }

    //Find normalized weights
    for( int iWeightStep = 0;iWeightStep < (iEmbeddingDimension+1);iWeightStep++ ){
        dWeights[iWeightStep] = dWeights[iWeightStep]/dWeightNormalization;
    }

}

/*CCMcorr -- Convergent Cross Mapped Correlation
  This function returns the (square) Pearson correlation coefficent
  between Y and the convergent cross mapped Y given X, i.e. Y
  estimated with weights calculated from the shadow manifold
  of X or Y|X.

  References: http://en.wikipedia.org/wiki/Convergent_cross_mapping
              https://www.sciencemag.org/content/338/6106/496.figures-only
              https://www.sciencemag.org/content/338/6106/496/suppl/DC1

  Dependencies:   eDist (function call)
                  FindAndCheckIndex (function call)
                  Pcorr (function call)
                  math.h (header)

  Input Parameters:
      double dY[] - time series to compared to its convergent
                   cross mapped counterpart given X
      int iY_length - library length of Y
      double dX_UsedForShadow[] - time series used to create
                                 the shadow manifold for
                                 estimating Y
      int iX_UsedForShadow_length - library length of X
      int iEmbeddingDimension - dimension of the shadow manifold
                               (notice that this values should be at
                               least 3 (i.e. 2+1) because this is a
                               pariwise CCM correlation)
      int iLagTime - lag time step used to construct the delay vectors
                    that make up the n-dimensional points of the
                    shadow manifold
      double dYestimate[] - array of doubles into which the estimated
                            values of Y (calculated using the weights
                            found from the shadow manifold) are placed
      int iYestimate_length - length of dYestimate

  Output Parameters:
      double - square of the Pearson correlation coefficent
               between Y and Y|X
*/
void CCMcorr(double &dPcorrYYX, double dY[],int iY_length, double dX_UsedForShadow[],int iX_UsedForShadow_length, int iEmbeddingDimension,int iLagTime, double dYestimate[], int iYestimate_length, bool bL2){

    //Check Min
    //if( iEmbeddingDimension < 3){
    //    fprintf(stderr, "Error in CCMCorr(): embedding dimension is %i which is less than 3\n", iEmbeddingDimension);
    //    return(nan(""));
    //}

    //Check Max
    //TODO: Figure out maximum for Embedding dimension (it depends on
    //      LagTime and the library length of X)
    //if( EmbeddingDimension  MAX_TBD){
    //    fprintf(stderr, "Error in CCMCorr(): embedding dimension is %i which is more than %i\n", EmbeddingDimension,MAX_TBD);
    //    return(-1);
    //}

    //Find calculated shadow manifold dimension
    int iCalShadManDim = iX_UsedForShadow_length-((iEmbeddingDimension-1)*iLagTime);

    //Assign some memory for the shadow manifold
    double** dXShadow;
    dXShadow = new double*[iCalShadManDim];
    for(int iter = 0;iter < iCalShadManDim;iter++ ){
        dXShadow[iter] = new double[iEmbeddingDimension];
    }
    int iTstep4delayvector;

    //Create the shadow manifold by creating the delay vectors sequentially
    for(int iShadowStep = 1; iShadowStep <= iCalShadManDim; iShadowStep++ ){
        iTstep4delayvector = iShadowStep+((iEmbeddingDimension-1)*iLagTime);
        for(int iDimStep = 1; iDimStep <= iEmbeddingDimension; iDimStep++ ){
            dXShadow[iShadowStep-1][iDimStep-1] = dX_UsedForShadow[(iTstep4delayvector-((iDimStep-1)*iLagTime))-1];
        }
    }

    //Storage for sorting, calculating, and organizing
    double dDelayVectorOfInterest[iEmbeddingDimension],  //delay vector currently being compared
           dWeights[(iEmbeddingDimension+1)],  //weights to contruct Y estimate
           dYEstimateGivenX[iCalShadManDim];  //estimated Y used the shadow manifold of X

    int iTstepOfNearestNeighborsTempRow[iEmbeddingDimension+1];//time steps of nearest neighbor delay vectors

    //Find starting point in Y for the estimate
    int iYStart = (iEmbeddingDimension-1)*iLagTime;

    //Loop through all of the delay vectors to find it distance to every other one
    for(int iDelayVectorN = 1; iDelayVectorN <= iCalShadManDim;iDelayVectorN++ ){

        //Populate the temp delay vector
        for(int iCopyStep = 0; iCopyStep < iEmbeddingDimension; iCopyStep++ ){
            dDelayVectorOfInterest[iCopyStep] = dXShadow[iDelayVectorN-1][iCopyStep];
        }

	//printf("BE -- %.20f,%.20f",dWeights[0],dWeights[1]);
        FindWeightsFromShadow(dWeights,iTstepOfNearestNeighborsTempRow,dDelayVectorOfInterest,dXShadow,iCalShadManDim,iEmbeddingDimension,iLagTime,iEmbeddingDimension);
	//printf("AF -- %.20f,%.20f",dWeights[0],dWeights[1]);

        //Find Y point estimates from X shadow manifold
        dYEstimateGivenX[iDelayVectorN-1] = 0;
        for( int iWeightStep = 0;iWeightStep < (iEmbeddingDimension+1);iWeightStep++ ){
            dYEstimateGivenX[iDelayVectorN-1] += dWeights[iWeightStep]*dY[(iYStart+iTstepOfNearestNeighborsTempRow[iWeightStep])];
        }
    }

    //Free the dXshadow memory
    for(int iter = 0;iter < iCalShadManDim;iter++ ){
        delete dXShadow[iter];
    }
    delete dXShadow;

    //Clip Y to calculate correlation
    double dYclipped[iCalShadManDim];
    for(int iYstep = 0 ;iYstep < (iY_length-iYStart);iYstep++ ){
        dYclipped[iYstep] = dY[iYStart+iYstep];
    }

    //Put the estimate into the proper array
    if( iYestimate_length != iCalShadManDim ){
        fprintf(stderr, "Error in CCMCorr(): size of estimated Y container (%i) is smaller than size of estimated Y (%i)\n", iYestimate_length,iCalShadManDim);
    }
//    dYestimate = dYclipped;
    for(int iCopyStep = 0;iCopyStep < iCalShadManDim;iCopyStep++ ){
	//printf("Compare: %.20f,%.20f\n",dYclipped[iCopyStep],dYestimate[iCopyStep]);
	dYestimate[iCopyStep] = dYEstimateGivenX[iCopyStep];
    }

    //Find correlation or L2 norm of Y and its estimate
    double dDiff,dDiffSquared;
    double dRunningSum = 0;
    if( !bL2 ){
    	dPcorrYYX = Pcorr(dYclipped,dYEstimateGivenX,iCalShadManDim);
    }else{
	for(int iCopyStep = 0;iCopyStep < iCalShadManDim;iCopyStep++ ){
	    dDiff = dYclipped[iCopyStep] - dYEstimateGivenX[iCopyStep];
	    dDiffSquared = dDiff*dDiff;
            dRunningSum += dDiffSquared;
        }
    	dPcorrYYX = sqrt(dRunningSum);
    }
}

/*CCMcorr2 -- Convergent Cross Mapped Correlation 2
  This function returns the (square) Pearson correlation coefficent
  between Y and the convergent cross mapped Y given X and Z, i.e. Y
  estimated with weights calculated from the shadow manifold
  of X and Z or Y|X,Z.  The shadow manifold is constructed with the
  requested embedding dimension used for both X and Z, so the final
  shadow manifold has an embedding dimension of twice the requested
  embedding dimension.  See the code for details.  It is assumed that
  Y, X, and Z are all the same length.

  References: http://en.wikipedia.org/wiki/Convergent_cross_mapping
              https://www.sciencemag.org/content/338/6106/496.figures-only
              https://www.sciencemag.org/content/338/6106/496/suppl/DC1

  Dependencies:   eDist (function call)
                  FindAndCheckIndex (function call)
                  Pcorr (function call)
                  math.h (header)

  Input Parameters:
      double dY[] - time series to compared to its convergent
                   cross mapped counterpart given X
      int iY_length - library length of Y
      double dX_UsedForShadow[] - time series 1 used to create
                                 the shadow manifold for
                                 estimating Y
      int iX_UsedForShadow_length - library length of X
      double dZ_UsedForShadow[] - time series 2 used to create
                                 the shadow manifold for
                                 estimating Y
      int iZ_UsedForShadow_length - library length of Z
      int iEmbeddingDimension - dimension of the shadow manifold
                               (notice that this values should be at
                               least 3 (i.e. 2+1) because this is a
                               pariwise CCM correlation)
      int iLagTime - lag time step used to construct the delay vectors
                    that make up the n-dimensional points of the
                    shadow manifold
      double dYestimate[] - array of doubles into which the estimated
                            values of Y (calculated using the weights
                            found from the shadow manifold) are placed
      int iYestimate_length - length of dYestimate

  Output Parameters:
      double - square of the Pearson correlation coefficent
               between Y and Y|X
*/
void CCMcorr2(double &dPcorrYYX, double dY[],int iY_length,double dX_UsedForShadow[],int iX_UsedForShadow_length,double dZ_UsedForShadow[],int iZ_UsedForShadow_length,int iEmbeddingDimension,int iLagTime,int iYEmbeddingDimension,int iYLagTime,double dYestimate[],int iYestimate_length,bool bL2){

    //Check Min
    //if( iEmbeddingDimension < 3){
    //    fprintf(stderr, "Error in CCMCorr(): embedding dimension is %i which is less than 3\n", iEmbeddingDimension);
    //    return(nan(""));
    //}

    //Check Max
    //TODO: Figure out maximum for Embedding dimension (it depends on
    //      LagTime and the library length of X)
    //if( EmbeddingDimension  MAX_TBD){
    //    fprintf(stderr, "Error in CCMCorr(): embedding dimension is %i which is more than %i\n", EmbeddingDimension,MAX_TBD);
    //    return(-1);
    //}

    if( iYEmbeddingDimension < 1){
        fprintf(stderr, "Error in CCMCorr2():  Y embedding dimension is %i which is not valid\n", iYEmbeddingDimension);
        dPcorrYYX = nan("");
    }
    if( iYLagTime < 1){
        fprintf(stderr, "Error in CCMCorr2():  Y lag time is %i which is not valid\n", iYLagTime);
        dPcorrYYX  = nan("");
    }

    //Find calculated shadow manifold dimension
    int iCalShadManDim = iX_UsedForShadow_length-((iEmbeddingDimension-1)*iLagTime);

    //Assign some memory for the shadow manifold
    int iEsum = (iEmbeddingDimension+iYEmbeddingDimension);
    double** dXZShadow;
    dXZShadow = new double*[iCalShadManDim];
    for(int iter = 0;iter < iCalShadManDim;iter++ ){
        dXZShadow[iter] = new double[iEsum];
    }
    int iTstep4delayvector,iYTstep4delayvector;

    //Create the shadow manifold by creating the delay vectors sequentially
    for(int iShadowStep = 1; iShadowStep <= iCalShadManDim; iShadowStep++ ){
        iTstep4delayvector = iShadowStep+((iEmbeddingDimension-1)*iLagTime);
        for(int iDimStep = 1; iDimStep <= iEmbeddingDimension; iDimStep++ ){
            dXZShadow[iShadowStep-1][iDimStep-1] = dX_UsedForShadow[(iTstep4delayvector-((iDimStep-1)*iLagTime))-1];
        }
	iYTstep4delayvector = iShadowStep+((iYEmbeddingDimension-1)*iYLagTime);
        for(int iDimStep = iEmbeddingDimension+1; iDimStep <= iYEmbeddingDimension; iDimStep++ ){
            dXZShadow[iShadowStep-1][iDimStep-1] = dZ_UsedForShadow[(iYTstep4delayvector-((iDimStep-1)*iYLagTime))-1];
        }
    }

    //Storage for sorting, calculating, and organizing
    double dDelayVectorOfInterest[iEsum],  //delay vector currently being compared
           dWeights[(iEmbeddingDimension+1)],  //weights to contruct Y estimate
           dYEstimateGivenX[iCalShadManDim];  //estimated Y used the shadow manifold of X

    int iTstepOfNearestNeighborsTempRow[iEmbeddingDimension+1];//time steps of nearest neighbor delay vectors

    //Find starting point in Y for the estimate
    int iYStart = (iEmbeddingDimension-1)*iLagTime;

    //Loop through all of the delay vectors to find it distance to every other one
    for(int iDelayVectorN = 1; iDelayVectorN <= iCalShadManDim;iDelayVectorN++ ){

        //Populate the temp delay vector
        for(int iCopyStep = 0; iCopyStep < iEsum; iCopyStep++ ){
            dDelayVectorOfInterest[iCopyStep] = dXZShadow[iDelayVectorN-1][iCopyStep];
        }

        FindWeightsFromShadow(dWeights,iTstepOfNearestNeighborsTempRow,dDelayVectorOfInterest,dXZShadow,iCalShadManDim,iEmbeddingDimension,iLagTime,iEsum);

        //Find Y point estimates from X shadow manifold
        dYEstimateGivenX[iDelayVectorN-1] = 0;
        for( int iWeightStep = 0;iWeightStep < (iEmbeddingDimension+1);iWeightStep++ ){
            dYEstimateGivenX[iDelayVectorN-1] += dWeights[iWeightStep]*dY[(iYStart+iTstepOfNearestNeighborsTempRow[iWeightStep])];
        }
    }

    //Free the dXshadow memory
    for(int iter = 0;iter < iCalShadManDim;iter++ ){
        delete dXZShadow[iter];
    }
    delete dXZShadow;

    //Clip Y to calculate correlation
    double dYclipped[iCalShadManDim];
    for(int iYstep = 0 ;iYstep < (iY_length-iYStart);iYstep++ ){
        dYclipped[iYstep] = dY[iYStart+iYstep];
    }

    //Put the estimate into the proper array
    if( iYestimate_length != iCalShadManDim ){
        fprintf(stderr, "Error in CCMCorr(): size of estimated Y container (%i) is smaller than size of estimated Y (%i)\n", iYestimate_length,iCalShadManDim);
    }
//    dYestimate = dYclipped;

    for(int iCopyStep = 0;iCopyStep < iCalShadManDim;iCopyStep++ ){
	//printf("Compare: %.20f,%.20f\n",dYclipped[iCopyStep],dYestimate[iCopyStep]);
	dYestimate[iCopyStep] = dYEstimateGivenX[iCopyStep];
    }

    //Find correlation of Y and its estimate
    double dDiff,dDiffSquared;
    double dRunningSum = 0;
    if( !bL2 ){
    	dPcorrYYX = Pcorr(dYclipped,dYEstimateGivenX,iCalShadManDim);
    }else{
	for(int iCopyStep = 0;iCopyStep < iCalShadManDim;iCopyStep++ ){
	    dDiff = dYclipped[iCopyStep] - dYEstimateGivenX[iCopyStep];
	    dDiffSquared = dDiff*dDiff;
            dRunningSum += dDiffSquared;
        }
    	dPcorrYYX = sqrt(dRunningSum);
    }
}

/*eDist -- Euclidean Distance
  This function returns the Euclidean distance between two
  vectors X and Y.  It is assumed that X and Y are the same
  length.

  References: http://en.wikipedia.org/wiki/Euclidean_distance
              http://www.mathworks.com/help/matlab/ref/norm.html#bt0y64c-1

  Dependencies:   math.h (header)

  Input Parameters:
      double dX[] - vector of values
      double dY[] - vector of values
      int iXY_length - length of X and Y

  Output Parameters:
      double - Euclidean distance between X and Y
*/
double eDist(double dX[],double dY[],int iXY_length){

    double dRunningSum = 0;

    for(int iter = 0; iter < iXY_length; iter++ ){
        dRunningSum += (dX[iter]-dY[iter])*(dX[iter]-dY[iter]);
    }

    return(sqrt(dRunningSum));
}

/*FindAndCheckIndex -- Find index location in array
  This function returns the index location of the element
  of a 1D array that is equal to a specificed value as long as
  that index value is not present in a specified 1D array of
  values.  It is assumed that a match will be made and an
  error is reported if no match is found.

  References:

  Dependencies:

  Input Parameters:
      double dValues2BCompared[] - 1D array of values in which the match
                                   is to be found
      int iValues2BCompared_length - length of Values2BCompared array
      double dValue2Compare - value used to define the match
      int iIndices2Check[] - 1D array of indices to check against
      int iIndices2Check_length - length of Indicies2Check array
  Output Parameters:
      int - index location for value in dValues2BCompared[] that is
            equal to dValue2Compare, as long as that index value is
            not in iIndices2Check[]
*/
int FindAndCheckIndex(double dValues2BCompared[],int iValues2BCompared_length,double dValue2Compare,int iIndices2Check[],int iIndices2Check_length){

    bool bUniqueIndex = true,
         bFoundMatch = false;

    for(int iter = 0;iter < iValues2BCompared_length;iter++ ){
        if( dValues2BCompared[iter] == dValue2Compare ){
            bFoundMatch = true;
        }
        if( bFoundMatch ){
            bUniqueIndex = true;
            for(int iIndexCheckStep = 0;iIndexCheckStep < iIndices2Check_length;iIndexCheckStep++ ){
                bUniqueIndex &= (iter != iIndices2Check[iIndexCheckStep]);
            }
            if( bUniqueIndex ){
                return(iter);
            }
        }
        bFoundMatch = false;
    }

    fprintf(stderr, "Error in FindAndCheckIndex(): no match was found, so an index of -1 will be returned\n");
    return(-1);
}

/*Pcorr -- Pearson's Correlation Coefficent
  This function returns Pearson's correlation coefficent
  between X and Y.  It is assumed (i.e. not checked) that
  X and Y are the same length.

  References: http://mathworld.wolfram.com/CorrelationCoefficient.html

  Dependencies:   math.h (header)

  Input Parameters:
      double dX[] - vector of values
      double dY[] - vector of values
      int iXY_length - length of X and Y

  Output Parameters:
      double - Pearson's correlation coefficent between X and Y
*/
double Pcorr(double dX[],double dY[],int iXY_length){

    double dEX = 0,  //expectation value of X
           dEY = 0,  //expectation value of Y
           dEXdelta,  //distance from mean for point in X
           dEYdelta,  //distance from mean for point in Y
           dSX = 0,  //variance of X
           dSY = 0,  //variance of Y
           dSXY = 0,  //covariance of X and Y
           dCovNorm;  //normalization factor for the covariance

    for(int iter = 0;iter < iXY_length;iter++ ){
        dEX += dX[iter];
        dEY += dY[iter];
    }

    dEX /= iXY_length;
    dEY /= iXY_length;

    for(int iter = 0;iter < iXY_length;iter++ ){
        dEXdelta = dX[iter]-dEX;
        dEYdelta = dY[iter]-dEY;
        dSX += dEXdelta*dEXdelta;
        dSY += dEYdelta*dEYdelta;
        dSXY += dEXdelta*dEYdelta;
    }

    dCovNorm = sqrt(dSX*dSY);
    if( dCovNorm == 0 ){
        dCovNorm = nan("");
        fprintf(stderr, "Warning in Pcorr(): division by zero\n");
    }

    return(dSXY/dCovNorm);
}

/*iMode -- Find the mode of a set of integers
  This function finds the mode of a set of integers

  References: http://cforbeginners.com/mode_c++.html
              http://en.wikipedia.org/wiki/Mode_(statistics)

  Dependencies:

  Input Parameters:
        int iA[] - set of integers
        int iAlength - length of A

  Output Parameters:
        int - mode of A

*/
int iMode(int iA[],int iAlength){

    int iRepetitions[iAlength];

    for(int iter = 0; iter < iAlength;iter++ ){
        iRepetitions[iter] = 0;
        for(int iter2 = 0;iter2 < iAlength;iter2++ ){
            if( iA[iter2] == iA[iter] ){
                iRepetitions[iter]++;
            }
        }
    }

    int iMaxIter = 0,
        iCurrentRepCount = iRepetitions[0];
    for (int iter = 1;iter < iAlength;iter++ ){
        if( iRepetitions[iter] > iCurrentRepCount ){
            iMaxIter = iter;
            iCurrentRepCount = iRepetitions[iter];
        }
    }

    return( iA[iMaxIter] );
}

