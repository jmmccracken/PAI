#include <stdio.h>
#include <math.h>
#include <algorithm>

void CCMcorr(double&,double[],int,double[],int,int,int,double[],int,bool);
void CCMcorr2(double&,double[],int,double[],int,double[],int,int,int,int,int,double[],int,bool);
double eDist(double[],double[],int);
int FindAndCheckIndex(double[],int,double,int[],int);
double Pcorr(double[], double[], int);
int FindLagTimeStep(double[],int,int,double,bool);
void FindWeightsFromShadow(double[],int[],double[],double**,int,int,int,int);
int FindEmbeddingDimension(double[],int,int,double,int,bool);
int iMode(int[],int);

