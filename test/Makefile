all:
	gfortran -c ../m_nofoutput.F90
	gfortran -c ../m_vars.F90
	gfortran -c ../m_cg.F90
	gfortran -w -ffixed-form -fd-lines-as-comments -c ../m_lbfgs_intern.F
	gfortran -c ../m_integd.F90
	gfortran -c ../m_rdmd.F90
	gfortran -c ../m_elag.F90
	gfortran -c ../m_diagf.F90
	gfortran -c ../m_gammatodm2.F90
	gfortran -c ../m_e_grad_occ.F90
	gfortran -c ../m_optocc.F90
	gfortran -c ../m_optorb.F90
	gfortran -c ../m_noft_driver.F90
	gfortran -llapack -lblas *.o test.F90 -o test.x
clean:
	/bin/rm -rf *.o *.mod test.x

cleannoft:
	/bin/rm -rf CANON_COEF  DYSON_COEF  F_DIAG  fort.2  GAMMAS NO_COEF  NO_COEF_BIN hubbard.noft
