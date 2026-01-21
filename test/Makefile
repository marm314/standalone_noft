all:
	make testF90
	make testCpp
testF90:
	chmod 777 ../gitversion.sh 
	../gitversion.sh 
	gfortran -O2 -c ../m_nofoutput.F90
	gfortran -O2 -c ../m_definitions.F90
	gfortran -O2 -c ../m_fcidump_nof.F90
	gfortran -O2 -c ../m_lbfgs_intern.F
	gfortran -O2 -c ../m_integd.F90
	gfortran -O2 -c ../m_rdmd.F90
	gfortran -O2 -c ../m_anti2unit.F90
	gfortran -O2 -c ../m_elag.F90
	gfortran -O2 -c ../m_diagf.F90
	gfortran -O2 -c ../m_adam.F90
	gfortran -O2 -c ../m_hessian.F90
	gfortran -O2 -c ../m_gammatodm2.F90
	gfortran -O2 -c ../m_e_grad_occ.F90
	gfortran -O2 -c ../m_e_grad_occ_cpx.F90
	gfortran -O2 -c ../m_optocc.F90
	gfortran -O2 -c ../m_optorb.F90
	gfortran -O2 -c ../m_tz_pccd_amplitudes.F90
	gfortran -O2 -c gitver.F90
	gfortran -O2 -c ../m_noft_driver.F90
	gfortran -O2  *.o test.F90 -llapack -lblas -o test.x
testCpp:
	g++ -O2 test.cpp \
	-L./../ \
	-lnoft -llapack -lblas \
	-Wl,-rpath,'../' \
	-o testc.x	

clean:
	/bin/rm -rf *.o *.mod gitver.F90 test.x testc.x 

cleannoft:
	/bin/rm -rf CANON_COEF  DYSON_COEF  F_DIAG  fort.2  GAMMAS NO_COEF  NO_COEF_BIN hubbard.noft
