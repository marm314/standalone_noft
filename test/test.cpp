#include<iostream>
#include<cstdio>
#include"../noft.h" 

using namespace std;

int main()
{
 int INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect;
 int imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis;
 int restart,ireadGAMMAS,ireadocc,ireadCOEF,ireadFdiag,iNOTupdateocc,iNOTupdateORB;
 double Enof,tolE,Vnn;
 double *Occ,*NO_COEF;

 NBF_tot=10;
 Occ=new double[NBF_tot];
 NO_COEF=new double[NBF_tot*NBF_tot];

 cout<<"Hello world cpp"<<endl;
 run_noft_c(&INOF,&Ista,&NBF_tot,&NBF_occ,&Nfrozen,&Npairs,&Ncoupled,&Nbeta_elect,&Nalpha_elect,
            &imethocc,&imethorb,&itermax,&iprintdmn,&iprintswdmn,&iprintints,&itolLambda,&ndiis,
            &Enof,&tolE,&Vnn,Occ,NO_COEF,&restart,&ireadGAMMAS,&ireadocc,&ireadCOEF,&ireadFdiag,
            &iNOTupdateocc,&iNOTupdateORB);

 delete[] Occ;Occ=NULL;
 delete[] NO_COEF;NO_COEF=NULL;
 return 0;
}

extern "C" void coef_2_hcore(double *NO_COEF_v,int *NBF)
{

}

extern "C" void hcore_ijkl(int *iorb, int *jorb, double *hVal)
{


}

