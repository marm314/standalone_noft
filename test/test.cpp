#include<iostream>
#include<cstdio>
#include"../noft.h" 
#define zero 0e0
#define one 1e0

using namespace std;

double *hCORE,*ERI;

int main(int argc, char *argv[])
{
 int INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect;
 int imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis;
 int restart,ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag,iNOTupdateOCC,iNOTupdateORB;
 int iguess;
 int iorb,jorb,korb;
 double Enof,tolE,Vnn;
 double *Occ,*NO_COEF,*Overlap;

 // Initialize variables
 INOF=8;Ista=0;NBF_tot=0;NBF_occ=0;
 Nfrozen=0;Npairs=0;Ncoupled=1;Nbeta_elect=0;Nalpha_elect=0;
 imethocc=1;imethorb=1;itermax=10000;iprintdmn=0;iprintswdmn=0;iprintints=0;
 itolLambda=5;ndiis=5;restart=0;
 tolE=1e-9;Vnn=zero;

 if(argc!=6)
 {
  cout<<" Include the arguments:"<<endl;
  cout<<endl;
  cout<<"  INOF(integer) Nbasis(integer) Nalpha_electrons(integer) iguess(0 Hcore, 1 init MO basis) restart(integer)"<<endl;
  cout<<endl;
  return 0;
 }
 INOF=atoi(argv[1]);
 NBF_tot=atoi(argv[2]);
 Nalpha_elect=atoi(argv[3]);
 iguess=atoi(argv[4]);
 restart=atoi(argv[5]);
 NBF_occ=NBF_tot;
 Nbeta_elect=Nalpha_elect;
 Npairs=(Nbeta_elect+Nalpha_elect)/2;
 Ncoupled=(NBF_tot/Npairs)-1;
 if(restart==0)
 {
  ireadGAMMAS=0;ireadOCC=0;ireadCOEF=0;ireadFdiag=0;iNOTupdateOCC=0;iNOTupdateORB=0;
 }
 else
 {
  ireadGAMMAS=1;ireadOCC=1;ireadCOEF=1;ireadFdiag=1;iNOTupdateOCC=0;iNOTupdateORB=0;
 }
 
 // Allocate arrays
 Occ=new double[NBF_tot];
 NO_COEF=new double[NBF_tot*NBF_tot];
 Overlap=new double[NBF_tot*NBF_tot];
 for(iorb=0;iorb<NBF_tot;iorb++) // The overlap matrix is the identity
 {
  Occ[iorb]=zero;
  if(iorb<Npairs) Occ[iorb]=Occ[iorb]=one;
  for(jorb=0;jorb<NBF_tot;jorb++)
  {
   Overlap[jorb+NBF_tot*iorb]=zero;
   if(iorb==jorb)
   {
    Overlap[iorb+NBF_tot*iorb]=one;
   }
  }
 }
 if(iguess==0)  // Hcore basis
 {
  
 }
 else  // The FCIDUMP basis is the identity
 {
  for(iorb=0;iorb<NBF_tot;iorb++)
  {
   for(jorb=0;jorb<NBF_tot;jorb++)
   {
    NO_COEF[jorb+NBF_tot*iorb]=zero;
    if(iorb==jorb)
    {
     NO_COEF[iorb+NBF_tot*iorb]=one;
    }
   }
  }
 }

 // Call the module
 cout<<" Hello world cpp"<<endl;
 run_noft_c(&INOF,&Ista,&NBF_tot,&NBF_occ,&Nfrozen,&Npairs,&Ncoupled,&Nbeta_elect,&Nalpha_elect,
            &imethocc,&imethorb,&itermax,&iprintdmn,&iprintswdmn,&iprintints,&itolLambda,&ndiis,
            &Enof,&tolE,&Vnn,Occ,Overlap,NO_COEF,&restart,&ireadGAMMAS,&ireadOCC,&ireadCOEF,
            &ireadFdiag,&iNOTupdateOCC,&iNOTupdateORB);

 // Deallocate arrays
 delete[] Occ;Occ=NULL;
 delete[] NO_COEF;NO_COEF=NULL;
 delete[] Overlap;Overlap=NULL;
 return 0;
}

extern "C" void coef_2_hcore(double *NO_COEF_v,int *NBF)
{
 
}

extern "C" void hcore_ij(int *iorb,int *jorb,int *NBF,double *hVal)
{
 hVal[0]=hCORE[ (iorb[0]-1) + (jorb[0]-1)*NBF[0] ];
}

extern "C" void coef_2_ERI(double *NO_COEF_v,int *NBF)
{

}

extern "C" void ERI_ijkl(int *iorb,int *jorb,int *korb,int *lorb,int *NBF,double *EVal)
{
 int NBF2,NBF3;
 NBF2=NBF[0]*NBF[0];
 NBF3=NBF2*NBF[0];
 EVal[0]=ERI[ (iorb[0]-1) + (jorb[0]-1)*NBF[0] + (korb[0]-1)*NBF2 + (lorb[0]-1)*NBF3 ];
}

