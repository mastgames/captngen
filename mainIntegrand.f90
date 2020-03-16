!	Capt'n OPER General
!	program for testing
!	the good stuff is in gencap.f90
!	runs the capture code one isotope at a time, letting the user sum the total later
!	to compile this version of main, switch out 'MAIN = main.o' for 'MAIN = mainOper.o' in the make file (line 8)

PROGRAM GENCAP
	use opermod
	implicit none
	character*300 :: modfile, filename
	double precision :: capped(100), result(100), u(100)
	double precision :: beginpwr, endpwr, pwr
	double precision, external :: integrand_oper
	integer :: index, cpl, indexlength
	! character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20", "Na23", "Mg24", &
	! 										"Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]
	character (len=5) :: cplConsts(14) = [character(len=5) :: "c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", &
											"c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

	! from DarkSUSY:
	modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_ags05_nohead.dat"
	! modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_agss09_nohead.dat"
	! modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_agss09ph_nohead.dat"
	! modfile = "solarmodels/Models_From_DarkSUSY/Serenelli-model_gs98_nohead.dat"
	
	indexlength=100

	cpl = 1
	call captn_init(modfile, 0.4d0, 220.d0, 220.d0, 540.d0)
	call captn_init_oper()
	if (cpl==1) then
		call populate_array(1.65d-8, cpl, 0)
	else
		call populate_array(1.65d-8, cpl+1, 0)
	endif
	
	beginpwr = -3.0
	endpwr = 6.0
	do index=0,indexlength
		pwr = beginpwr + index/indexlength * (endpwr - beginpwr)
		u(index) = 10**pwr
	end do

	mdm = 80 !mass of dm, accessed via module
    j_chi = 0.5 !spin of dm, accessed via module
    niso = 16 !number of isotopes, accessed via module
	pickIsotope = 1 !the isotope, accessed via module
	ri_for_omega = 0 !the radius entry of the star, accessed via the module
	do index=0,indexlength
		result(index) = integrand_oper(u(index))
	end do

	filename = "integrand_oper_as_fn_of_u.dat"
	open(55,file=filename)
	do index=0,indexlength
		write(55,*) u(index), result(index)
	end do
	close(55)
	print*
END PROGRAM GENCAP
