all:
	make testF90
	make testCpp
testF90:
	chmod 777 ../gitversion.sh 
	../gitversion.sh 
	gfortran -O2 -fopenmp -c ../m_nofoutput.F90
	gfortran -O2 -fopenmp -c ../m_definitions.F90
	gfortran -O2 -fopenmp -c ../m_fcidump_nof.F90
	gfortran -O2 -fopenmp -c ../m_lbfgs_intern.F90
	gfortran -O2 -fopenmp -c ../m_integd.F90
	gfortran -O2 -fopenmp -c ../m_rdmd.F90
	gfortran -O2 -fopenmp -c ../m_anti2unit.F90
	gfortran -O2 -fopenmp -c ../m_elag.F90
	gfortran -O2 -fopenmp -c ../m_diagf.F90
	gfortran -O2 -fopenmp -c ../m_adam.F90
	gfortran -O2 -fopenmp -c ../m_hessian.F90
	gfortran -O2 -fopenmp -c ../m_gammatodm2.F90
	gfortran -O2 -fopenmp -c ../m_e_grad_occ.F90
	gfortran -O2 -fopenmp -c ../m_e_grad_occ_cpx.F90
	gfortran -O2 -fopenmp -c ../m_optocc.F90
	gfortran -O2 -fopenmp -c ../m_optorb.F90
	gfortran -O2 -fopenmp -c ../m_tz_pccd_amplitudes.F90
	gfortran -O2 -fopenmp -c gitver.F90
	gfortran -O2 -fopenmp -c ../m_noft_driver.F90
	gfortran -O2 -fopenmp  *.o test.F90 -llapack -lblas -o test.x
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
