#include<iostream>
#include<cstdio>
#include<cstring>
#include<fstream>
#include<algorithm>
#include<cmath>
#include<iomanip>
#include"../noft.h" 
#define zero 0e0
#define one 1e0

using namespace std;

void read_fcidump(int &NBF_tot,double &Vnn);
void jacobi_donof(int n, double **m, double **v);
void sort_donof(int n, double **m, double **v);
int type;
double *hCORE_init,*ERI_init; // Read from FCIDUMP
double *hCORE_tran,*ERI_tran;       // Transformed

int main(int argc, char *argv[])
{
 int INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect;
 int imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis;
 int restart,ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag,iNOTupdateOCC,iNOTupdateORB;
 int iguess,ifort_fcidump;
 int iorb,jorb,korb;
 double Enof,tolE,Vnn;
 double *Occ,*NO_COEF,*Overlap;
 double **Cguess,**TMP_MAT;

 // Initialize variables
 INOF=8;Ista=0;NBF_tot=0;NBF_occ=0;
 Nfrozen=0;Npairs=0;Ncoupled=1;Nbeta_elect=0;Nalpha_elect=0;
 imethocc=1;imethorb=1;itermax=10000;iprintdmn=0;iprintswdmn=0;iprintints=0;
 itolLambda=5;ndiis=5;restart=0;ifort_fcidump=0;
 tolE=1e-9;Vnn=zero;Enof=zero;

 if(argc!=7)
 {
  cout<<" Include the arguments:"<<endl;
  cout<<endl;
  cout<<"  INOF(integer) Nbasis(integer) Nalpha_electrons(integer) iguess(0 Hcore, 1 init MO basis) restart(integer) type_fcidump(integer)"<<endl;
  cout<<endl;
  return 0;
 }
 INOF=atoi(argv[1]);
 NBF_tot=atoi(argv[2]);
 Nalpha_elect=atoi(argv[3]);
 iguess=atoi(argv[4]);
 restart=atoi(argv[5]);
 type=atoi(argv[6]);
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
 hCORE_tran=new double[NBF_tot*NBF_tot];
 hCORE_init=new double[NBF_tot*NBF_tot];
 ERI_tran=new double[NBF_tot*NBF_tot*NBF_tot*NBF_tot];
 ERI_init=new double[NBF_tot*NBF_tot*NBF_tot*NBF_tot];
 for(iorb=0;iorb<NBF_tot;iorb++) // The overlap matrix is the identity
 {
  Occ[iorb]=zero;
  if(iorb<Npairs) Occ[iorb]=one;
  for(jorb=0;jorb<NBF_tot;jorb++)
  {
   hCORE_tran[jorb+NBF_tot*iorb]=zero;
   hCORE_init[jorb+NBF_tot*iorb]=zero;
   Overlap[jorb+NBF_tot*iorb]=zero;
   if(iorb==jorb)
   {
    Overlap[iorb+NBF_tot*iorb]=one;
   }
  }
 }
 for(iorb=0;iorb<NBF_tot*NBF_tot*NBF_tot*NBF_tot;iorb++)
 {
  ERI_tran[iorb]=zero;
  ERI_init[iorb]=zero;
 }

 // Read the FCIDUMP file
 read_fcidump(NBF_tot,Vnn);

 // Prepare the MO guess
 if(iguess==0)  // Hcore basis
 {
  Cguess=new double*[NBF_tot];
  TMP_MAT=new double*[NBF_tot];
  for(iorb=0;iorb<NBF_tot;iorb++)
  {
   Cguess[iorb]=new double[NBF_tot];
   TMP_MAT[iorb]=new double[NBF_tot];
   for(jorb=0;jorb<NBF_tot;jorb++)
   {
    TMP_MAT[iorb][jorb]=hCORE_init[iorb+jorb*NBF_tot];
    Cguess[iorb][jorb]=zero;
   }
  }
  jacobi_donof(NBF_tot,TMP_MAT,Cguess);
  sort_donof(NBF_tot,TMP_MAT,Cguess);
  for(iorb=0;iorb<NBF_tot;iorb++)
  {
   for(jorb=0;jorb<NBF_tot;jorb++)
   {
    NO_COEF[jorb+iorb*NBF_tot]=Cguess[jorb][iorb];
   }
  }
  for(iorb=0;iorb<NBF_tot;iorb++)
  {
   delete[] Cguess[iorb];Cguess[iorb]=NULL;
   delete[] TMP_MAT[iorb];TMP_MAT[iorb]=NULL;
  }
  delete[] Cguess;Cguess=NULL;
  delete[] TMP_MAT;TMP_MAT=NULL;
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
 run_noft_c(&INOF,&Ista,&NBF_tot,&NBF_occ,&Nfrozen,&Npairs,&Ncoupled,&Nbeta_elect,&Nalpha_elect,
            &imethocc,&imethorb,&itermax,&iprintdmn,&iprintswdmn,&iprintints,&itolLambda,&ndiis,
            &Enof,&tolE,&Vnn,Occ,Overlap,NO_COEF,&restart,&ireadGAMMAS,&ireadOCC,&ireadCOEF,
            &ireadFdiag,&iNOTupdateOCC,&iNOTupdateORB,&ifort_fcidump);

 // Print results
 cout<<endl;
 cout<<setprecision(8)<<fixed<<endl;
 cout<<" Energy "<<setw(20)<<Enof<<endl;
 cout<<endl;
 cout<<" Occupancies"<<endl;
 jorb=0;
 for(iorb=0;iorb<NBF_tot;iorb++)
 {
  cout<<setw(20)<<Occ[iorb];
  jorb++;
  if(jorb==5)
  {
   cout<<endl;
   jorb=0;
  }
 }
 cout<<endl;

 // Deallocate arrays
 delete[] Occ;Occ=NULL;
 delete[] NO_COEF;NO_COEF=NULL;
 delete[] Overlap;Overlap=NULL;
 delete[] hCORE_tran;hCORE_tran=NULL;
 delete[] hCORE_init;hCORE_init=NULL;
 delete[] ERI_tran;ERI_tran=NULL;
 delete[] ERI_init;ERI_init=NULL;
 return 0;
}

extern "C" void coef_2_hcore(double *NO_COEF_v,int *NBF)
{
 int iorb,jorb,korb,porb,NBF1,NBF2;
 double *hCORE_tmp;
 NBF1=NBF[0];
 NBF2=NBF1*NBF1;
 hCORE_tmp=new double[NBF2];
 //
 // Full change: h_kp -> h_ji
 //
 // h_kp -> h_jp
 for(jorb=0;jorb<NBF1;jorb++)
 {
  for(porb=0;porb<NBF1;porb++)
  {
   hCORE_tmp[jorb+porb*NBF1]=zero;
   for(korb=0;korb<NBF1;korb++)
   {
    hCORE_tmp[jorb+porb*NBF1]+=NO_COEF_v[korb+jorb*NBF1]*hCORE_init[korb+porb*NBF1];
   }
  }
 }
 // h_jp -> h_ji
 for(jorb=0;jorb<NBF1;jorb++)
 {
  for(iorb=0;iorb<NBF1;iorb++)
  {
   hCORE_tran[jorb+iorb*NBF1]=zero;
   for(porb=0;porb<NBF1;porb++)
   {
    hCORE_tran[jorb+iorb*NBF1]+=NO_COEF_v[porb+iorb*NBF1]*hCORE_tmp[jorb+porb*NBF1];
   }
  }
 }
 delete[] hCORE_tmp;hCORE_tmp=NULL;
}

extern "C" void hcore_ij(int *iorb,int *jorb,int *NBF,double *hVal)
{
 hVal[0]=hCORE_tran[ (iorb[0]-1) + (jorb[0]-1)*NBF[0] ];
}

extern "C" void coef_2_ERI(double *NO_COEF_v,int *NBF)
{
 int iorb,jorb,korb,lorb,porb,qorb,rorb,sorb,NBF1,NBF2,NBF3,NBF4;
 double *ERI_tmp;
 NBF1=NBF[0];
 NBF2=NBF1*NBF1;
 NBF3=NBF2*NBF1;
 NBF4=NBF3*NBF1;
 ERI_tmp=new double[NBF4];
 //
 // Full change: <ij|kl> -> <pq|rs>
 //
 // <ij|kl> -> <pj|kl>
 for(porb=0;porb<NBF1;porb++)
 {
  for(jorb=0;jorb<NBF1;jorb++)
  {
   for(korb=0;korb<NBF1;korb++)
   {
    for(lorb=0;lorb<NBF1;lorb++)
    {
     ERI_tmp[porb+jorb*NBF1+korb*NBF2+lorb*NBF3]=zero;
     for(iorb=0;iorb<NBF1;iorb++)
     {
      ERI_tmp[porb+jorb*NBF1+korb*NBF2+lorb*NBF3]+=NO_COEF_v[iorb+porb*NBF1]*ERI_init[iorb+jorb*NBF1+korb*NBF2+lorb*NBF3];
     }
    }
   }
  }
 } 
 // <pj|kl> -> <pq|kl>
 for(porb=0;porb<NBF1;porb++)
 {
  for(qorb=0;qorb<NBF1;qorb++)
  {
   for(korb=0;korb<NBF1;korb++)
   {
    for(lorb=0;lorb<NBF1;lorb++)
    {
     ERI_tran[porb+qorb*NBF1+korb*NBF2+lorb*NBF3]=zero;
     for(jorb=0;jorb<NBF1;jorb++)
     {
      ERI_tran[porb+qorb*NBF1+korb*NBF2+lorb*NBF3]+=NO_COEF_v[jorb+qorb*NBF1]*ERI_tmp[porb+jorb*NBF1+korb*NBF2+lorb*NBF3];
     }
    }
   }
  }
 }
 // <pq|kl> -> <pq|rl>
 for(porb=0;porb<NBF1;porb++)
 {
  for(qorb=0;qorb<NBF1;qorb++)
  {
   for(rorb=0;rorb<NBF1;rorb++)
   {
    for(lorb=0;lorb<NBF1;lorb++)
    {
     ERI_tmp[porb+qorb*NBF1+rorb*NBF2+lorb*NBF3]=zero;
     for(korb=0;korb<NBF1;korb++)
     {
      ERI_tmp[porb+qorb*NBF1+rorb*NBF2+lorb*NBF3]+=NO_COEF_v[korb+rorb*NBF1]*ERI_tran[porb+qorb*NBF1+korb*NBF2+lorb*NBF3];
     }
    }
   }
  }
 }
 // <pq|rl> -> <pq|rs>
 for(porb=0;porb<NBF1;porb++)
 {
  for(qorb=0;qorb<NBF1;qorb++)
  {
   for(rorb=0;rorb<NBF1;rorb++)
   {
    for(sorb=0;sorb<NBF1;sorb++)
    {
     ERI_tran[porb+qorb*NBF1+rorb*NBF2+sorb*NBF3]=zero;
     for(lorb=0;lorb<NBF1;lorb++)
     {
      ERI_tran[porb+qorb*NBF1+rorb*NBF2+sorb*NBF3]+=NO_COEF_v[lorb+sorb*NBF1]*ERI_tmp[porb+qorb*NBF1+rorb*NBF2+lorb*NBF3];
     }
    }
   }
  }
 }
 delete[] ERI_tmp;ERI_tmp=NULL;
}

extern "C" void ERI_ijkl(int *iorb,int *jorb,int *korb,int *lorb,int *NBF,double *EVal)
{
 int NBF2,NBF3;
 NBF2=NBF[0]*NBF[0];
 NBF3=NBF2*NBF[0];
 EVal[0]=ERI_tran[ (iorb[0]-1) + (jorb[0]-1)*NBF[0] + (korb[0]-1)*NBF2 + (lorb[0]-1)*NBF3 ];
}

void read_fcidump(int &NBF_tot,double &Vnn)
{
 int iorb,jorb,korb,lorb;
 double Val;
 string line;
 if(type==0) // Standard
 {
  iorb=system(" echo '   0.0   -1   -1  -1   -1 ' >> FCIDUMP ");
 }
 else // CMZ style
 {
  iorb=system(" echo '  -1   -1   -1   -1   0.0 ' >> FCIDUMP ");
 }
 ifstream read_dump("FCIDUMP");
 if(type==0)
 {
  do
  {
   getline(read_dump,line);
   line.erase(std::remove_if(line.begin(),line.end(),::isspace),line.end());
  }while(line!="$" && line!="&");
  do
  {
   read_dump>>Val>>iorb>>jorb>>korb>>lorb;
   if(iorb==0 && jorb==0 && korb==0 && lorb==0){Vnn=Val;}
   if(iorb!=0 && jorb!=0 && korb==0 && lorb==0)
   {
    hCORE_init[ (iorb-1) + (jorb-1)*NBF_tot ]=Val;
    hCORE_init[ (jorb-1) + (iorb-1)*NBF_tot ]=Val;
   }
   if(iorb!=-1 && jorb!=-1 && korb!=-1 && lorb!=-1 && iorb!=0 && jorb!=0 && korb!=0 && lorb!=0)
   {
    ERI_init[ (iorb-1) + (korb-1)*NBF_tot + (jorb-1)*NBF_tot*NBF_tot + (lorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (iorb-1) + (lorb-1)*NBF_tot + (jorb-1)*NBF_tot*NBF_tot + (korb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (jorb-1) + (lorb-1)*NBF_tot + (iorb-1)*NBF_tot*NBF_tot + (korb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (jorb-1) + (korb-1)*NBF_tot + (iorb-1)*NBF_tot*NBF_tot + (lorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (korb-1) + (iorb-1)*NBF_tot + (lorb-1)*NBF_tot*NBF_tot + (jorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (lorb-1) + (iorb-1)*NBF_tot + (korb-1)*NBF_tot*NBF_tot + (jorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (lorb-1) + (jorb-1)*NBF_tot + (korb-1)*NBF_tot*NBF_tot + (iorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (korb-1) + (jorb-1)*NBF_tot + (lorb-1)*NBF_tot*NBF_tot + (iorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
   }
  }while(iorb!=-1 && jorb!=-1 && korb!=-1 && lorb!=-1);
 }
 else
 {
  do
  {
   read_dump>>iorb>>jorb>>korb>>lorb>>Val;
   if(iorb==0 && jorb==0 && korb==0 && lorb==0){Vnn=Val;}
   if(iorb!=0 && jorb!=0 && korb==0 && lorb==0)
   {
    hCORE_init[ (iorb-1) + (jorb-1)*NBF_tot ]=Val;
    hCORE_init[ (jorb-1) + (iorb-1)*NBF_tot ]=Val;
   }
   if(iorb!=-1 && jorb!=-1 && korb!=-1 && lorb!=-1 && iorb!=0 && jorb!=0 && korb!=0 && lorb!=0)
   {
    ERI_init[ (iorb-1) + (korb-1)*NBF_tot + (jorb-1)*NBF_tot*NBF_tot + (lorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (iorb-1) + (lorb-1)*NBF_tot + (jorb-1)*NBF_tot*NBF_tot + (korb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (jorb-1) + (lorb-1)*NBF_tot + (iorb-1)*NBF_tot*NBF_tot + (korb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (jorb-1) + (korb-1)*NBF_tot + (iorb-1)*NBF_tot*NBF_tot + (lorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (korb-1) + (iorb-1)*NBF_tot + (lorb-1)*NBF_tot*NBF_tot + (jorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (lorb-1) + (iorb-1)*NBF_tot + (korb-1)*NBF_tot*NBF_tot + (jorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (lorb-1) + (jorb-1)*NBF_tot + (korb-1)*NBF_tot*NBF_tot + (iorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
    ERI_init[ (korb-1) + (jorb-1)*NBF_tot + (lorb-1)*NBF_tot*NBF_tot + (iorb-1)*NBF_tot*NBF_tot*NBF_tot ]=Val;
   }
  }while(iorb!=-1 && jorb!=-1 && korb!=-1 && lorb!=-1);
 }
 read_dump.close();
}

void jacobi_donof(int n, double **m, double **v)
{
 int i,j,k,ip,iq,maxiter=10000;
 double TWO=2e0,SIX=6e0,TEN=1e1;
 double tol = pow(TEN,-(TWO*SIX)),apq,t,alpha,c,s,tau,d,delta,temp1,temp2;
 //Define v(matrix eigenvectors) as identity matrix
 for(i=0;i<n;i++)
 {
  for(j=0;j<n;j++)
  {v[i][j]=zero;}
   v[i][i] = one;
 }
 //Calculate Delta
 delta=zero;
 for(i=0;i<n;i++)
 {
  for(j=i+1;j<n;j++)
  {delta=delta+m[i][j]*m[i][j];}
 }
 //Iterations cycle
 for(i=0;i<maxiter;i++)
 {
  apq=zero;
  ip=0;
  iq=1;
  for(j=0;j<n;j++)
  {
   for(k=j+1;k<n;k++)
   {
    if(abs(m[j][k])>apq)
    {
     ip = j;
     iq = k;
     apq = abs(m[ip][iq]);
    }
   }
  }
  //Determine c, s, c, tau
  t=one;
  if(m[ip][ip]!=m[iq][iq])
  {
   alpha=(m[iq][iq]-m[ip][ip])/(TWO*m[ip][iq]);
   t=-alpha+alpha/abs(alpha)*sqrt(alpha*alpha+one);
  }
  c=one/sqrt(t*t+one);
  s=t*c;
  tau=s/(one+c);
  d=m[ip][iq];
  //Update matrix m en and the upper diagonal part
  m[ip][ip]=m[ip][ip]-t*m[ip][iq];
  m[iq][iq]=m[iq][iq]+t*m[ip][iq];
  for(j=0;j<ip;j++)
  {
   temp1=m[j][ip];
   temp2=m[j][iq];
   m[j][ip]=temp1-s*(temp2+tau*temp1);
   m[j][iq]=temp2+s*(temp1-tau*temp2);
  }
  for(j=ip+1;j<iq;j++)
  {
   temp1=m[ip][j];
   temp2=m[j][iq];
   m[j][iq]=temp2+s*(temp1-tau*temp2);
   m[ip][j]=temp1-s*(temp2+tau*temp1);
  }
  for(j=iq+1;j<n;j++)
  {
   temp1=m[ip][j];
   temp2=m[iq][j];
   m[iq][j]=temp2+s*(temp1-tau*temp2);
   m[ip][j]=temp1-s*(temp2+tau*temp1);
  }
  m[ip][iq]=zero;
  //Update v
  for(j=0;j<n;j++)
  {
   temp1=v[j][ip];
   temp2=v[j][iq];
   v[j][ip]=c*temp1-s*temp2;
   v[j][iq]=s*temp1+c*temp2;
  }
  //Update delta
  delta=delta-d*d;
  //If it has converge, update the lower diagonal part of m and finish jacobi()
  if(abs(delta)<=tol)
  {
   for(j=0;j<n;j++)
   {
    for(k=j+1;k<n;k++)
    {m[k][j] = m[j][k];}
   }
   return;
  }
 }
}

void sort_donof(int n, double **m, double **v)
{
 int iorb,jorb,korb;
 int *order;
 double *TMP_VEC,val;
 TMP_VEC=new double[n];
 order=new int[n];
 for(iorb=0;iorb<n;iorb++)
 {
  order[iorb]=iorb;
  TMP_VEC[iorb]=m[iorb][iorb];
 } 
 for(iorb=0;iorb<n;iorb++)
 {
  for(jorb=iorb+1;jorb<n;jorb++)
  {
   if(TMP_VEC[iorb]>TMP_VEC[jorb])
   {
    val=TMP_VEC[iorb];
    TMP_VEC[iorb]=TMP_VEC[jorb];
    TMP_VEC[jorb]=val;
    korb=order[iorb];
    order[iorb]=order[jorb];
    order[jorb]=korb;
   }
  }
 }
 for(iorb=0;iorb<n;iorb++)
 {
  for(jorb=0;jorb<n;jorb++)
  {
   m[iorb][jorb]=v[iorb][order[jorb]];
  }
 }
 for(iorb=0;iorb<n;iorb++)
 {
  for(jorb=0;jorb<n;jorb++)
  {
   v[iorb][jorb]=m[iorb][jorb];
   m[iorb][jorb]=zero;
  }
 }
 for(iorb=0;iorb<n;iorb++)
 {
  m[iorb][iorb]=TMP_VEC[iorb];
 }
 delete[] TMP_VEC;TMP_VEC=NULL;
 delete[] order;order=NULL;
}

