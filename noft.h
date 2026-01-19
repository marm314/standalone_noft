#ifndef NOFT_H
#define NOFT_H

extern "C" {

 void run_noft_c(int *INOF,int *Ista,int *NBF_tot,int *NBF_occ,int *Nfrozen, int *Npairs,
                 int *Ncoupled,int *Nbeta_elect,int *Nalpha_elect,int *imethocc, int *imethorb,
                 int *itermax,int *iprintdmn,int *iprintswdmn,int *iprintints,int *itolLambda,
                 int *ndiis,double *Enof,double *tolE,double *Vnn,double *Occ,double *Overlap_in,
                 double *NO_COEF_in, int *restart,int *ireadGAMMAS,int *ireadocc,int *ireadCOEF,
                 int *ireadFdiag, int *iNOTupdateocc,int *iNOTupdateORB,int *ifort_fcidump,
		 int *iskip_fcidump,int *istyle_fcidump);
 void coef_2_hcore(double *NO_COEF_v,int *NBF);
 void coef_2_ERI(double *NO_COEF_v,int *NBF);
 void hcore_ij(int *iorb,int *jorb,int *NBF,double *hVal);
 void ERI_ijkl(int *iorb,int *jorb,int *korb,int *lorb,int *NBF,double *EVal);
  
}

#endif
