all:
	chmod 777 gitversion.sh 
	./gitversion.sh 
	gfortran -Wall -c m_nofoutput.F90
	gfortran -Wall -c m_definitions.F90
	gfortran -Wall -c m_lbfgs_intern.F
	gfortran -Wall -c m_integd.F90
	gfortran -Wall -c m_rdmd.F90
	gfortran -Wall -c m_elag.F90
	gfortran -Wall -c m_diagf.F90
	gfortran -Wall -c m_gammatodm2.F90
	gfortran -Wall -c m_e_grad_occ.F90
	gfortran -Wall -c m_e_grad_occ_cpx.F90
	gfortran -Wall -c m_optocc.F90
	gfortran -Wall -c m_optorb.F90
	gfortran -Wall -c m_tz_pccd_amplitudes.F90
	gfortran -Wall -c gitver.F90
	gfortran -Wall -c m_noft_driver.F90
clean:
	/bin/rm -rf *.o *.mod gitver.F90
tar:
	tar -cvf module_noft.tar *F90 *F README test
