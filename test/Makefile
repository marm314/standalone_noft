all:
	chmod 777 ../gitversion.sh 
	../gitversion.sh 
	gfortran -O3 -c m_gitver.F90
	gfortran -O3 -c ../m_nofoutput.F90
	gfortran -O3 -c ../m_vars.F90
	gfortran -O3 -c ../m_cg.F90
	gfortran -O3 -w -ffixed-form -fd-lines-as-comments -c ../m_lbfgs_intern.F
	gfortran -O3 -c ../m_integd.F90
	gfortran -O3 -c ../m_rdmd.F90
	gfortran -O3 -c ../m_elag.F90
	gfortran -O3 -c ../m_diagf.F90
	gfortran -O3 -c ../m_gammatodm2.F90
	gfortran -O3 -c ../m_e_grad_occ.F90
	gfortran -O3 -c ../m_optocc.F90
	gfortran -O3 -c ../m_optorb.F90
	gfortran -O3 -c ../m_noft_driver.F90
	gfortran -O3 -llapack -lblas *.o test.F90 -o test.x
clean:
	/bin/rm -rf *.o *.mod m_gitver.F90 test.x

cleannoft:
	/bin/rm -rf CANON_COEF  DYSON_COEF  F_DIAG  fort.2  GAMMAS NO_COEF  NO_COEF_BIN hubbard.noft
