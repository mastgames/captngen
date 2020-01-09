!   GPU-CUDA-enabled captn oper 
!   For a continuum of Q-dependent capture
!   Simplified, general solar DM capture routine
!   Standalone code for q^2n, v^2n
!   Useful stuff is run at the end; beginning is the module that does the heavy lifting
!   Future plans: add form factor handling (a la Catena & Schwabe)
!   Made for GAMBIT, with marginal competence
!   Aaron Vincent 2017
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.
!   Reference q0 is 40 MeV, and v0 is 220 km/s.


module capgpu
    use cudafor
    implicit none
    double precision, parameter :: pi=3.141592653, NAvo=6.0221409d23, GMoverR=1.908e15
    double precision, parameter :: Rsun=69.57d9
    double precision, parameter :: c0=2.99792458d10, mnuc=0.938, q0=0.04, v0=220.d5, hbar=6.582d-25
    double precision, parameter :: eps=1d-10 !stops divisions by zero
    !these are now set in captn_init
    double precision :: usun , u0 ,rho0, vesc_halo
    !this goes with the Serenelli table format
    double precision, parameter :: AtomicNumber(29) = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                                                        18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                                                        39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                                                        54.93, 55.845, 58.933, 58.693/) !29 is the max niso, corresponding to Ni
    double precision, parameter :: AtomicNumber_oper(16) = (/ 1., 3., 4., 12., 14., 16., 20., 23., 24., 27., &
                                                        28., 32., 40., 40., 56., 58./) !the isotopes the catena paper uses
    character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24", &
                                                                "Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]
    double precision, parameter :: AtomicSpin_oper(16) = (/ 0.5, 0.5, 0., 0., 1., 0., 0., 1.5, 0., 2.5, &
                                                        0., 0., 0., 0., 0., 0./) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    !tab: means tabulated from file; so as not to be confused with other variables
    double precision, allocatable :: tab_mencl(:), tab_starrho(:), tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:)

    ! nq and nv can be -1, 0, 1, 2; this is set in the main program
    integer :: nq, nv, niso, ri_for_omega, nlines, pickIsotope
    double precision :: mdm, sigma_0, j_chi
    double precision :: coupling_Array(14,2)
    double precision :: W_array(8,16,2,2,7)
    double precision, allocatable :: tab_mfr_oper(:,:)
    
    contains

    !this is to make sure the integrator does what it's supposed to
    function gaussinmod(x)
        double precision :: x,gaussinmod
        gaussinmod = nq*exp(-x**2/2.d0)
    end function gaussinmod


    attributes(device) function GFFI_H_oper(w,vesc,mq)
        double precision :: p, mu,w,vesc,u,muplus,GFFI_H_oper,G
        integer mq
        p = mdm*w
        mu = mdm/mnuc
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        if (mq .ne. -1) then
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*1./(1.+mq)*((mu/muplus**2)**(mq+1)-(u**2/w**2)**(mq+1))
        else
            !eps added to make the log finite: the value of eps does not affect the answer
            G = (p/c0)**(2.d0*mq)*mdm*w**2/(2.d0*mu**mq)*log(mu/muplus**2*w**2/(u+eps)**2)
        endif
        GFFI_H_oper = G
    end function GFFI_H_oper
    
    attributes(device) function GFFI_A_oper(w,vesc,A,mq)
        double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
        double precision :: dgamic,GFFI_A_oper
        integer :: mq
        p = mdm*w
        mu = mdm/mnuc/A
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        mN = A*mnuc
        Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
        B = .5*mdm*w**2/Ei/c0**2
        if (mq .eq. 0) then
            GFFI_A_oper = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
        else
            GFFI_A_oper = ((p+eps)/c0)**(2*mq)*Ei*c0**2/(B*mu)**mq*(dgamic(1.+dble(mq),B*u**2/w**2+eps) &
                - dgamic(1.+dble(mq),B*mu/muplus**2+eps))
        end if
    end function GFFI_A_oper
    

    ! breaks the function W into each term, and sums them with the corresponding GFFI
    attributes(device) function sumW(w,vesc,iso,tau,tauprime,Wtype,qOffset)
        double precision :: w,vesc,yConverse,sumW,tally
        integer :: k,iso,tau,tauprime,Wtype,qOffset
        double precision :: G
        ! y = yConverse * q^2
        ! yConverse is the conversion factor to go from q^2 to y
        yConverse = 2/3.*((0.91*(mnuc*AtomicNumber_oper(iso))**(1./3)+0.3)*10**-13)**2/(2*hbar*c0)**2
        tally = 0
        k = threadidx%z
        if (k.lt.8) then
        if (iso.eq.1) then
            G = GFFI_H_oper(w,vesc,(k-1+qOffset))
        else
            G = GFFI_A_oper(w,vesc,AtomicNumber_oper(iso),(k-1+qOffset))
        end if
        tally = tally + W_array(Wtype,iso,tau,tauprime,k) * yConverse**(k-1) * G
        sumW = tally
    end function sumW
    
    
    ! this is eqn 3.23 in 1501.03729
    ! large sum handled through expansion of R functions
    ! many if statements used to check for terms excluded by choice of constants (c1,c3..c15) = 0
    attributes(device) function p_tot(w,vesc,i)
        double precision :: w,vesc, p_tot
        double precision :: m_T,mu_T,GF
        integer :: i,tau,taup
        integer :: c, v2, q2, v2q2, q4, v2q4
        double precision :: RD, RM, RMP2, RP1, RP2, RS1, RS1D, RS2
        
        c = 0
        v2 = 1
        q2 = 2
        v2q2 = 3
        q4 = 4
        v2q4 = 5

        ! get the target nucleus mass m_T and the nucleus-dm reduced mass mu_T
        m_T = mnuc * AtomicNumber_oper(i)
        mu_T = (m_T*mdm)/(m_T+mdm)

        ! use GFFI_H for hydrogen or y~0
        ! do a check to see if y can simplify W to a constant
        
        ! the sum in eqn 3.23 of paper 1501.03729
        ! expanded the R functions into each term of consatnt, v^2n, and q^2n
        ! check to see if choice of couping constants has zeroed out any R term before doinng any more work with it
        ! if not zerod out, is added to sum (sumW expands the W functions as polynomials as given in paper's appendix)
        ! note! coupling constants are 2d array of length 14 (c1,c3,c4...c14,c15) (note absence of c2)
        ! this results in array call of index [1] -> c1, but a call of index [2] -> c3 !
        p_tot = 0.0
        tau = threadidx%x
        taup = threadidx%y
        if ((tau.lt.3).and.(taup.lt.3)) then
            ! RM (c, v2, q2, v2q2)
            ! c1,c1
            if ((coupling_Array(1,tau).ne.0).and.(coupling_Array(1,taup).ne.0)) then
                p_tot = p_tot + RM(mnuc,c0,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,1,0)
            end if
            ! c8,c8
            if ((coupling_Array(7,tau).ne.0).and.(coupling_Array(7,taup).ne.0)) then
                p_tot = p_tot + RM(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,1,0) &
                            - RM(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,1,1) 
            end if
            ! c11,c11
            if ((coupling_Array(10,tau).ne.0).and.(coupling_Array(10,taup).ne.0)) then
                p_tot = p_tot + RM(mnuc,c0,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,1,1)
            end if
            ! c5,c5
            if ((coupling_Array(4,tau).ne.0).and.(coupling_Array(4,taup).ne.0)) then
                p_tot = p_tot + RM(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,1,1) &
                            - RM(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,1,2) 
            end if
            
            ! RS2 (c, v2, q2, v2q2, q4)
            ! c4,c4
            if ((coupling_Array(3,tau).ne.0).and.(coupling_Array(3,taup).ne.0)) then
                p_tot = p_tot + RS2(mnuc,c0,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,2,0) 
            end if
            ! c12,c12
            if ((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0)) then
                p_tot = p_tot + RS2(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,2,0) &
                            - RS2(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,2,1)
            end if
            ! c10,c10 or c4,c6 or c6,c4
            if (((coupling_Array(9,tau).ne.0).and.(coupling_Array(9,taup).ne.0)).or. &
                    ((coupling_Array(3,tau).ne.0).and.(coupling_Array(5,taup).ne.0)).or. &
                    ((coupling_Array(5,tau).ne.0).and.(coupling_Array(3,taup).ne.0))) then
                p_tot = p_tot + RS2(mnuc,c0,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,2,1)
            end if
            ! c13,c13
            if ((coupling_Array(12,tau).ne.0).and.(coupling_Array(12,taup).ne.0)) then
                p_tot = p_tot + RS2(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,2,1) &
                            - RS2(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,2,2)
            end if
            ! c6,c6
            if ((coupling_Array(5,tau).ne.0).and.(coupling_Array(5,taup).ne.0)) then
                p_tot = p_tot + RS2(mnuc,c0,tau,taup,q4,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,2,2)
            end if
            
            ! RS1 (c, v2, q2, v2q2, v2q4)
            ! c4,c4
            if ((coupling_Array(3,tau).ne.0).and.(coupling_Array(3,taup).ne.0)) then
                p_tot = p_tot + RS1(mnuc,c0,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,3,0)
            end if
            ! c7,c7 or c12,c12
            if (((coupling_Array(6,tau).ne.0).and.(coupling_Array(6,taup).ne.0)).or. &
                    ((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
                p_tot = p_tot + RS1(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,3,0) &
                            - RS1(mnuc,c0,tau,taup,v2,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,1)
            end if
            ! c9,c9
            if ((coupling_Array(8,tau).ne.0).and.(coupling_Array(8,taup).ne.0)) then
                p_tot = p_tot + RS1(mnuc,c0,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,3,1)
            end if
            ! c3,c3 or c14,c14 or c15,c12 or c12,c15
            if (((coupling_Array(2,tau).ne.0).and.(coupling_Array(2,taup).ne.0)).or. &
                    ((coupling_Array(13,tau).ne.0).and.(coupling_Array(13,taup).ne.0)).or. &
                    ((coupling_Array(14,tau).ne.0).and.(coupling_Array(11,taup).ne.0)).or. &
                    ((coupling_Array(11,tau).ne.0).and.(coupling_Array(14,taup).ne.0))) then
                p_tot = p_tot + RS1(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array)  * w**2 * sumW(w,vesc,i,tau,taup,3,1)&
                            -RS1(mnuc,c0,tau,taup,v2q2,j_chi,coupling_Array)  * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,2)
            end if
            ! c15,c15
            if ((coupling_Array(14,tau).ne.0).and.(coupling_Array(14,taup).ne.0)) then
                p_tot = p_tot + RS1(mnuc,c0,tau,taup,v2q4,j_chi,coupling_Array) * w**2 * sumW(w,vesc,i,tau,taup,3,2) &
                            - RS1(mnuc,c0,tau,taup,v2q4,j_chi,coupling_Array) * c0**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,3) 
            end if
            
            ! RP2 (c, q2, q4)
            ! c12,c12
            if (((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RP2(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,4,1)
            end if
            ! c3,c3 or c12,c15 or c15,c12
            if (((coupling_Array(2,tau).ne.0).and.(coupling_Array(2,taup).ne.0)).or. &
                    ((coupling_Array(11,tau).ne.0).and.(coupling_Array(14,taup).ne.0)).or. &
                    ((coupling_Array(14,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RP2(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,4,2)
            end if
            ! c15,c15
            if (((coupling_Array(14,tau).ne.0).and.(coupling_Array(14,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RP2(mnuc,tau,taup,q4,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,4,3)
            end if
            
            ! RMP2 (c, q2)
            ! c3,c1 or c12,c11
            if (((coupling_Array(2,tau).ne.0).and.(coupling_Array(1,taup).ne.0)).or. &
                    ((coupling_Array(11,tau).ne.0).and.(coupling_Array(10,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RMP2(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,5,1)
            end if
            ! c15,c11
            if (((coupling_Array(14,tau).ne.0).and.(coupling_Array(10,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RMP2(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,5,2)
            end if
            
            ! RP1 (c, q2)
            ! c12,c12
            if (((coupling_Array(11,tau).ne.0).and.(coupling_Array(11,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RP1(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,6,1) 
            end if
            ! c13,c13
            if (((coupling_Array(12,tau).ne.0).and.(coupling_Array(12,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RP1(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,6,2)
            end if
            
            ! RD (c, q2)
            ! c8,c8
            if ((coupling_Array(7,tau).ne.0).and.(coupling_Array(7,taup).ne.0)) then
                p_tot = p_tot + 1./mnuc**2 * RD(mnuc,tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,7,1)
            end if
            ! c5,c5
            if ((coupling_Array(4,tau).ne.0).and.(coupling_Array(4,taup).ne.0)) then
                p_tot = p_tot + 1./mnuc**2 * RD(mnuc,tau,taup,q2,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,7,2)
            end if
            
            !RS1D (c)
            ! c5,c4 or c8,c9
            if (((coupling_Array(4,tau).ne.0).and.(coupling_Array(3,taup).ne.0)).or. &
                    ((coupling_Array(7,tau).ne.0).and.(coupling_Array(8,taup).ne.0))) then
                p_tot = p_tot + 1./mnuc**2 * RS1D(tau,taup,c,j_chi,coupling_Array) * sumW(w,vesc,i,tau,taup,8,1)
            end if
        end if
        p_tot = p_tot  * hbar**2 * c0**2
    end function p_tot

! CUDA requires the CPU to call a subroutine(not a function) with global attribute
! This subroutine is used to call the p_tot function with the device attribute
! (called from GPU, run on GPU)
    attributes(global) subroutine gpu_oper(p_tot,w,vesc,i)
        double precision,value :: w,vesc, p_tot
        integer,value :: i

        p_tot = p_tot(w,vesc,i)
    end subroutine gpu_oper


    !   this is eqn 2.1 in 1501.03729
    ! This calls the gpu kernal gpu_oper to find p_tot
    function OMEGA_oper(rindex,w)
        double precision :: w, vesc,mu,muplus,u,J, OMEGA_oper
        integer rindex, i
        double precision :: p_tot
        type(dim3) :: dimGrid, dimBlock
        
        dimGrid = dim3( 1, 1, 1 )
        dimBlock = dim3( 2, 2, 7 )

        vesc = tab_vesc(rindex)
        u = sqrt(w**2-vesc**2)
        Omega_oper = 0.d0
        
        ! if only one isotope is being run
        if ((0 .lt. pickIsotope).and.(pickIsotope .lt. 17)) then
            J = AtomicSpin_oper(pickIsotope)
            mu = mdm/mnuc/AtomicNumber_oper(pickIsotope)
            muplus = (1.+mu)/2.
            if (mu/muplus**2 .gt. u**2/w**2) then
                ! note the factor [(2*mnuc*AtomicNumber(pickIsotope))/(w**2*(2*J+1))],is the product of
                ! the factors in front of eqns 3.26 and 3.23 in paper 1501.03729
                gpu_oper(p_tot,w,vesc,pickIsotope)
                Omega_oper = Omega_oper+(NAvo*tab_starrho(rindex)/AtomicNumber_oper(pickIsotope)/mnuc* &
                    tab_mfr_oper(rindex,pickIsotope))*w*((2*mnuc*AtomicNumber_oper(pickIsotope))/(w**2*(2*J+1)))* &
                    p_tot
            end if
        else ! if all the isotopes are being run
            do i = 1,niso
                J = AtomicSpin_oper(i)
                mu = mdm/mnuc/AtomicNumber_oper(i)
                muplus = (1.+mu)/2.
                if (mu/muplus**2 .gt. u**2/w**2) then
                    ! note the factor [(2*mnuc*AtomicNumber(i))/(w**2*(2*J+1))],is the product of
                    ! the factors in front of eqns 3.26 and 3.23 in paper 1501.03729
                    gpu_oper(p_tot,w,vesc,i)
                    Omega_oper = Omega_oper+(NAvo*tab_starrho(rindex)/AtomicNumber_oper(i)/mnuc*tab_mfr_oper(rindex,i))*w* &
                                            ((2*mnuc*AtomicNumber_oper(i))/(w**2*(2*J+1)))*p_tot
                end if
            end do
        end if
    end function OMEGA_oper
end module capgpu


function dummyf(x)
    double precision :: x, dummyf
    dummyf = 1.d0
end function dummyf


subroutine captn_init_oper()
    use capmod
    implicit none
    integer :: i, j, k, l, m
    character (len=2) :: terms(7) = [character(len=2) :: "y0", "y1", "y2", "y3", "y4", "y5", "y6"]
    real :: WM, WS2, WS1, WP2, WMP2, WP1, WD, WS1D
    
    ! tab_mfr_oper is allocated in the get_solar_params subroutine
    ! take the regular array tab_mfr and extract the isotopes used in the 1501.03729 paper (otherwise indices won't match on arrays)
    do i=1,nlines
        tab_mfr_oper(i,1) = tab_mfr(i,1)
        tab_mfr_oper(i,2) = tab_mfr(i,3)
        tab_mfr_oper(i,3) = tab_mfr(i,2)
        tab_mfr_oper(i,4) = tab_mfr(i,4)
        tab_mfr_oper(i,5) = tab_mfr(i,6)
        tab_mfr_oper(i,6) = tab_mfr(i,8)
        tab_mfr_oper(i,7) = tab_mfr(i,11)
        tab_mfr_oper(i,8) = tab_mfr(i,12)
        tab_mfr_oper(i,9) = tab_mfr(i,13)
        tab_mfr_oper(i,10) = tab_mfr(i,14)
        tab_mfr_oper(i,11) = tab_mfr(i,15)
        tab_mfr_oper(i,12) = tab_mfr(i,17)
        tab_mfr_oper(i,13) = tab_mfr(i,19)
        tab_mfr_oper(i,14) = tab_mfr(i,21)
        tab_mfr_oper(i,15) = tab_mfr(i,27)
        tab_mfr_oper(i,16) = tab_mfr(i,29)
    end do
    
    ! this array stores each of the constants of the W polynomials from paper 1501.03729's appendix individually
    ! array index m handles the 8 varients of the W functions in order [M, S", S', P", MP", P', Delta, S'Delta]
    ! index i handles the 16 isotopes [H, He3, He4, C12, N14, O16, Ne20, Na 23, Mg24, Al27, Si28, S32, Ar40, Ca40, Fe56, Ni58]
    ! index j & k handle the two superscripts for each W function, each taking values of 0 and 1
    ! index L determines the power of each constant ranging from y^0 to y^6
    do m=1,8
        do i=1,16
            do j=1,2
                do k=1,2
                    do l=1,7
                        if (m.eq.1) then
                            W_array(m,i,j,k,l) = WM(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.2) then
                            W_array(m,i,j,k,l) = WS2(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.3) then
                            W_array(m,i,j,k,l) = WS1(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.4) then
                            W_array(m,i,j,k,l) = WP2(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.5) then
                            W_array(m,i,j,k,l) = WMP2(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.6) then
                            W_array(m,i,j,k,l) = WP1(j-1,k-1,isotopes(i),terms(l))
                        else if (m.eq.7) then
                            W_array(m,i,j,k,l) = WD(j-1,k-1,isotopes(i),terms(l))
                        else
                            W_array(m,i,j,k,l) = WS1D(j-1,k-1,isotopes(i),terms(l))
                        end if
                    end do
                end do
            end do
        end do
    end do

    ! initiate the coupling_Array (full of the coupling constants) with all zeros
    ! populate_array will place the non-zero value into a chosen slot at runtime
    coupling_Array = reshape((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), (/14, 2/))
end subroutine captn_init_oper


!   this is the integral over R in eqn 2.3 in 1501.03729
function integrand_oper(u)
    use capmod
    implicit none
    double precision :: u, w, vesc, integrand_oper, int
    vesc = tab_vesc(ri_for_omega)
    w = sqrt(u**2+vesc**2)
    int = get_vdist(u)/u*w*Omega_oper(ri_for_omega,w)
    integrand_oper = int
end function integrand_oper


!   Need to pass all the operators into the subroutine
subroutine captn_oper(mx_in, jx_in, niso_in, isotopeChosen, capped)
    use capmod
    implicit none
    integer, intent(in):: niso_in, isotopeChosen
    integer i, ri
    double precision, intent(in) :: mx_in, jx_in
    double precision :: capped, maxcap !this is the output
     ! array of coupling constants
    double precision :: epsabs, epsrel,limit,result,abserr,neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator
    !external gausstest !this is just for testing
    external integrand_oper
    external dummyf
    epsabs=1.d-17
    epsrel=1.d-17
    limit=1000

    mdm = mx_in
    j_chi = jx_in
    niso = niso_in
    
    pickIsotope = isotopeChosen

    ! temporary, the user will want to choose their coupling constants to match a model
    !                           c1,  c3,  c4, c5,   c6,  c7,  c8,  c9, c10, c11, c12, c13, c14, c15   
    !coupling_Array = reshape((/1.65d-8, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                !0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/), (/14, 2/))
    
    if (.not. allocated(tab_r)) then 
        print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
        return
    end if
    allocate(u_int_res(nlines))

    !As far as I can tell, the second argument (fofuoveru) does nothing in this integrator.
    !I've sent it to an inert dummy just in case.
    capped = 0.d0

    ! completes integral (2.3) in paper 1501.03729 (gives dC/dV as fn of radius)
    do ri=1,nlines !loop over the star
        result = 0.d0
        ri_for_omega = ri !accessed via the module
        !call integrator
        call dsntdqagse(integrand_oper,dummyf,1.d0,vesc_halo, &
            epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        u_int_res(ri) = result
        capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)
    end do

    !completes integral (2.4) of paper 1501.03729
    capped = 4.d0*pi*Rsun**3*capped

    if (capped .gt. 1.d100) then
      print*,"Capt'n General says: Oh my, it looks like you are capturing an  &
      infinite amount of dark matter in the Sun. Best to look into that."
    end if
end subroutine captn_oper
