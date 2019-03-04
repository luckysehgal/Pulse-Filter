#include "stdio.h"
#include <iostream>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "pulse7389_2g.h"

void root_sav_golay()
{
  const unsigned int N = 16384;
  float binX[N]; 
  
  //S-G Setting
  const unsigned int polynomial_level = 2;
  const unsigned int number_sg_point = 3001;
  const int half = (int)(number_sg_point/2);
  
  TMatrixD mX(number_sg_point,polynomial_level);
  TMatrixD mXT(polynomial_level,number_sg_point);
  TMatrixD mC(polynomial_level,polynomial_level);
  TMatrixD mF(polynomial_level,number_sg_point);
  TMatrixD mY(number_sg_point,1);
  TMatrixD mA(polynomial_level,1);

  double mx[number_sg_point];
  double my[number_sg_point];
  float A0[N], A1[N];
  // end of S-G defination

  for(int i=0;i<N;i++) 
    { binX[i] = i; }
  
  for(int j=0;j<number_sg_point;j++)
    { 
      mx[j] = (j-half-1); 
    }
  
  
  for(int i=0;i<number_sg_point;i++)
    {
      for(int j=0;j<polynomial_level;j++)
	{
	  mX(i,j) = pow(mx[i],j);
	}
    }
  
  mXT.Transpose(mX);
  mC = mXT*mX;
  mF = mC.Invert()*mXT;
 
 
 for(int i=0;i<N;i++)
   {
     if((i<=half)||(i>=(N-half)))
       {
	 A0[i] = pulse[i];
	 A1[i] = 0.0;
       }
     else
       {
	 for(int j=0;j<number_sg_point;j++)
	   { 
	     mY(j,0) = pulse[i-half+j]; 
	   }
	 mA = mF*mY;
	 A0[i] = (float)mA(0,0);
	 A1[i] = (float)mA(1,0);
       }
   }
 
 TGraph *gA0, *gA1;
 
 gA0 = new TGraph(N,binX,A0);
 gA0->SetName("gA0");
 gA0->SetMarkerStyle(1);
 
 gA1 = new TGraph(N,binX,A1);
 gA1->SetName("gA1");
 gA1->SetMarkerStyle(1);
 
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 
 
 TGraph *pp = new TGraph(N,binX,pulse);
 
 
 
 TCanvas *plot1 = new TCanvas("plot1");
 pp->Draw("al");
 gA0->SetLineColor(2);
 gA0->Draw("l");
 
 TCanvas *plot2 = new TCanvas("plot2");
 gA1->SetLineColor(2);
 gA1->Draw("al");


}
