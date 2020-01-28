!   Capt'n Oper for GPU
!   Module to house everying specific to captn operator, CUDA Fortran
!   Neal Avis Kozar 2020
!   all units of distance: cm
!   all units of mass/energy : GeV (or GeV/c^2, don't forget)
!   all units of time: seconds
!   Sticking with notation of 1504.04378. Cite that paper. Or 1605.06502 it's even better.


module gpumod
    use cudafor
    implicit none
    double precision, device :: hbar_dev=6.582d-25
    ! !this goes with the Serenelli table format
    
    ! integer :: pickIsotope
    ! double precision :: j_chi

    double precision, device :: AtomicNumber_oper_dev(16) = (/ 1., 3., 4., 12., 14., 16., 20., 23., 24., 27., &
                                                        28., 32., 40., 40., 56., 58./) !the isotopes the catena paper uses
    ! character (len=4) :: isotopes(16) = [character(len=4) :: "H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24", &
    !                                                             "Al27", "Si28","S32","Ar40","Ca40","Fe56","Ni58"]
    ! double precision, parameter :: AtomicSpin_oper(16) = (/ 0.5, 0.5, 0., 0., 1., 0., 0., 1.5, 0., 2.5, &
    !                                                     0., 0., 0., 0., 0., 0./) !spins pulled from https://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    double precision, device :: coupling_Array_dev(14,2)
    double precision, device :: W_array_dev(8,16,2,2,7)

    double precision, device :: mdm_dev
    double precision, device :: c0_dev=2.99792458d10, mnuc_dev=0.938
    double precision, device :: eps_dev=1d-10 !stops divisions by zero
    
    contains

    attributes(device) function GFFI_H_oper(w,vesc,mq)
        double precision :: p, mu,w,vesc,u,muplus,GFFI_H_oper,G
        integer mq
        p = mdm_dev*w
        mu = mdm_dev/mnuc_dev
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        if (mq .ne. -1) then
            G = (p/c0_dev)**(2.d0*mq)*mdm_dev*w**2/(2.d0*mu**mq)*1./(1.+mq)*((mu/muplus**2)**(mq+1)-(u**2/w**2)**(mq+1))
        else
            !eps_dev added to make the log finite: the value of eps_dev does not affect the answer
            G = (p/c0_dev)**(2.d0*mq)*mdm_dev*w**2/(2.d0*mu**mq)*log(mu/muplus**2*w**2/(u+eps_dev)**2)
        endif
        GFFI_H_oper = G
    end function GFFI_H_oper
    
    attributes(device) function GFFI_A_oper(w,vesc,A,mq)
        double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
        double precision :: dgamic,GFFI_A_oper
        integer :: mq
        p = mdm_dev*w
        mu = mdm_dev/mnuc_dev/A
        muplus = (1.+mu)/2.
        u = sqrt(w**2-vesc**2)
        mN = A*mnuc_dev
        Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
        B = .5*mdm_dev*w**2/Ei/c0_dev**2
        if (mq .eq. 0) then
            GFFI_A_oper = Ei*c0_dev**2*(exp(-mdm_dev*u**2/2/Ei/c0_dev**2)-exp(-B*mu/muplus**2))
        else
            GFFI_A_oper = ((p+eps_dev)/c0_dev)**(2*mq)*Ei*c0_dev**2/(B*mu)**mq*(dgamic(1.+dble(mq),B*u**2/w**2+eps_dev) &
                - dgamic(1.+dble(mq),B*mu/muplus**2+eps_dev))
        end if
    end function GFFI_A_oper
    
    ! breaks the function W into each term, and sums them with the corresponding GFFI
    attributes(device) function sumW(w,vesc,iso,tau,tauprime,Wtype,qOffset)
        double precision :: w,vesc,yConverse,sumW,tally
        integer :: k,iso,tau,tauprime,Wtype,qOffset
        double precision :: G
        ! y = yConverse * q^2
        ! yConverse is the conversion factor to go from q^2 to y
        yConverse = 2/3.*((0.91*(mnuc_dev*AtomicNumber_oper_dev(iso))**(1./3)+0.3)*10**-13)**2/(2*hbar_dev*c0_dev)**2
        tally = 0
        ! do k=1,7
        k = threadidx%z
        if (k.lt.8) then
            if (iso.eq.1) then
                G = GFFI_H_oper(w,vesc,(k-1+qOffset))
            else
                G = GFFI_A_oper(w,vesc,AtomicNumber_oper_dev(iso),(k-1+qOffset))
            end if
            tally = tally + W_array_dev(Wtype,iso,tau,tauprime,k) * yConverse**(k-1) * G
        end if
        ! end do
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
        double precision, device :: RD, RM_dev, RMP2, RP1, RP2, RS1, RS1D, RS2
        
        c = 0
        v2 = 1
        q2 = 2
        v2q2 = 3
        q4 = 4
        v2q4 = 5

        ! get the target nucleus mass m_T and the nucleus-dm reduced mass mu_T
        m_T = mnuc_dev * AtomicNumber_oper_dev(i)
        mu_T = (m_T*mdm_dev)/(m_T+mdm_dev)

        ! use GFFI_H for hydrogen or y~0
        ! do a check to see if y can simplify W to a constant
        
        ! the sum in eqn 3.23 of paper 1501.03729
        ! expanded the R functions into each term of consatnt, v^2n, and q^2n
        ! check to see if choice of couping constants has zeroed out any R term before doinng any more work with it
        ! if not zerod out, is added to sum (sumW expands the W functions as polynomials as given in paper's appendix)
        ! note! coupling constants are 2d array of length 14 (c1,c3,c4...c14,c15) (note absence of c2)
        ! this results in array call of index [1] -> c1, but a call of index [2] -> c3 !
        p_tot = 0.0
        ! do tau=1,2
        !     do taup=1,2
        tau = threadidx%x
        taup = threadidx%y
        if ((tau.lt.3).and.(taup.lt.3)) then
                ! RM (c, v2, q2, v2q2)
                ! c1,c1
                if ((coupling_Array_dev(1,tau).ne.0).and.(coupling_Array_dev(1,taup).ne.0)) then
                    p_tot = p_tot + RM_dev(mnuc_dev,c0_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,1,0)
                end if
                ! c8,c8
                if ((coupling_Array_dev(7,tau).ne.0).and.(coupling_Array_dev(7,taup).ne.0)) then
                    p_tot = p_tot + RM(mnuc_dev,c0_dev,tau,taup,v2,j_chi,coupling_Array_dev) * w**2 * sumW(w,vesc,i,tau,taup,1,0) &
                                - RM(mnuc_dev,c0_dev,tau,taup,v2,j_chi,coupling_Array_dev) * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,1,1) 
                end if
                ! c11,c11
                if ((coupling_Array_dev(10,tau).ne.0).and.(coupling_Array_dev(10,taup).ne.0)) then
                    p_tot = p_tot + RM(mnuc_dev,c0_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,1,1)
                end if
                ! c5,c5
                if ((coupling_Array_dev(4,tau).ne.0).and.(coupling_Array_dev(4,taup).ne.0)) then
                    p_tot = p_tot + RM(mnuc_dev,c0_dev,tau,taup,v2q2,j_chi,coupling_Array_dev) * w**2 * sumW(w,vesc,i,tau,taup,1,1) &
                                - RM(mnuc_dev,c0_dev,tau,taup,v2q2,j_chi,coupling_Array_dev) * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,1,2) 
                end if
                
                ! RS2 (c, v2, q2, v2q2, q4)
                ! c4,c4
                if ((coupling_Array_dev(3,tau).ne.0).and.(coupling_Array_dev(3,taup).ne.0)) then
                    p_tot = p_tot + RS2(mnuc_dev,c0_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,2,0) 
                end if
                ! c12,c12
                if ((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(11,taup).ne.0)) then
                    p_tot = p_tot + RS2(mnuc_dev,c0_dev,tau,taup,v2,j_chi,coupling_Array_dev) * w**2 * sumW(w,vesc,i,tau,taup,2,0) &
                                - RS2(mnuc_dev,c0_dev,tau,taup,v2,j_chi,coupling_Array_dev) * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,2,1)
                end if
                ! c10,c10 or c4,c6 or c6,c4
                if (((coupling_Array_dev(9,tau).ne.0).and.(coupling_Array_dev(9,taup).ne.0)).or. &
                        ((coupling_Array_dev(3,tau).ne.0).and.(coupling_Array_dev(5,taup).ne.0)).or. &
                        ((coupling_Array_dev(5,tau).ne.0).and.(coupling_Array_dev(3,taup).ne.0))) then
                    p_tot = p_tot + RS2(mnuc_dev,c0_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,2,1)
                end if
                ! c13,c13
                if ((coupling_Array_dev(12,tau).ne.0).and.(coupling_Array_dev(12,taup).ne.0)) then
                    p_tot = p_tot + RS2(mnuc_dev,c0_dev,tau,taup,v2q2,j_chi,coupling_Array_dev) * w**2 * sumW(w,vesc,i,tau,taup,2,1) &
                                - RS2(mnuc_dev,c0_dev,tau,taup,v2q2,j_chi,coupling_Array_dev) * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,2,2)
                end if
                ! c6,c6
                if ((coupling_Array_dev(5,tau).ne.0).and.(coupling_Array_dev(5,taup).ne.0)) then
                    p_tot = p_tot + RS2(mnuc_dev,c0_dev,tau,taup,q4,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,2,2)
                end if
                
                ! RS1 (c, v2, q2, v2q2, v2q4)
                ! c4,c4
                if ((coupling_Array_dev(3,tau).ne.0).and.(coupling_Array_dev(3,taup).ne.0)) then
                    p_tot = p_tot + RS1(mnuc_dev,c0_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,3,0)
                end if
                ! c7,c7 or c12,c12
                if (((coupling_Array_dev(6,tau).ne.0).and.(coupling_Array_dev(6,taup).ne.0)).or. &
                        ((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(11,taup).ne.0))) then
                    p_tot = p_tot + RS1(mnuc_dev,c0_dev,tau,taup,v2,j_chi,coupling_Array_dev) * w**2 * sumW(w,vesc,i,tau,taup,3,0) &
                                - RS1(mnuc_dev,c0_dev,tau,taup,v2,j_chi,coupling_Array_dev) * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,1)
                end if
                ! c9,c9
                if ((coupling_Array_dev(8,tau).ne.0).and.(coupling_Array_dev(8,taup).ne.0)) then
                    p_tot = p_tot + RS1(mnuc_dev,c0_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,3,1)
                end if
                ! c3,c3 or c14,c14 or c15,c12 or c12,c15
                if (((coupling_Array_dev(2,tau).ne.0).and.(coupling_Array_dev(2,taup).ne.0)).or. &
                        ((coupling_Array_dev(13,tau).ne.0).and.(coupling_Array_dev(13,taup).ne.0)).or. &
                        ((coupling_Array_dev(14,tau).ne.0).and.(coupling_Array_dev(11,taup).ne.0)).or. &
                        ((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(14,taup).ne.0))) then
                    p_tot = p_tot + RS1(mnuc_dev,c0_dev,tau,taup,v2q2,j_chi,coupling_Array_dev)  * w**2 * sumW(w,vesc,i,tau,taup,3,1)&
                                -RS1(mnuc_dev,c0_dev,tau,taup,v2q2,j_chi,coupling_Array_dev)  * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,2)
                end if
                ! c15,c15
                if ((coupling_Array_dev(14,tau).ne.0).and.(coupling_Array_dev(14,taup).ne.0)) then
                    p_tot = p_tot + RS1(mnuc_dev,c0_dev,tau,taup,v2q4,j_chi,coupling_Array_dev) * w**2 * sumW(w,vesc,i,tau,taup,3,2) &
                                - RS1(mnuc_dev,c0_dev,tau,taup,v2q4,j_chi,coupling_Array_dev) * c0_dev**2/(4.*mu_T**2) * sumW(w,vesc,i,tau,taup,3,3) 
                end if
                
                ! RP2 (c, q2, q4)
                ! c12,c12
                if (((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(11,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RP2(mnuc_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,4,1)
                end if
                ! c3,c3 or c12,c15 or c15,c12
                if (((coupling_Array_dev(2,tau).ne.0).and.(coupling_Array_dev(2,taup).ne.0)).or. &
                        ((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(14,taup).ne.0)).or. &
                        ((coupling_Array_dev(14,tau).ne.0).and.(coupling_Array_dev(11,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RP2(mnuc_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,4,2)
                end if
                ! c15,c15
                if (((coupling_Array_dev(14,tau).ne.0).and.(coupling_Array_dev(14,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RP2(mnuc_dev,tau,taup,q4,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,4,3)
                end if
                
                ! RMP2 (c, q2)
                ! c3,c1 or c12,c11
                if (((coupling_Array_dev(2,tau).ne.0).and.(coupling_Array_dev(1,taup).ne.0)).or. &
                        ((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(10,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RMP2(mnuc_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,5,1)
                end if
                ! c15,c11
                if (((coupling_Array_dev(14,tau).ne.0).and.(coupling_Array_dev(10,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RMP2(mnuc_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,5,2)
                end if
                
                ! RP1 (c, q2)
                ! c12,c12
                if (((coupling_Array_dev(11,tau).ne.0).and.(coupling_Array_dev(11,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RP1(mnuc_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,6,1) 
                end if
                ! c13,c13
                if (((coupling_Array_dev(12,tau).ne.0).and.(coupling_Array_dev(12,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RP1(mnuc_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,6,2)
                end if
                
                ! RD (c, q2)
                ! c8,c8
                if ((coupling_Array_dev(7,tau).ne.0).and.(coupling_Array_dev(7,taup).ne.0)) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RD(mnuc_dev,tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,7,1)
                end if
                ! c5,c5
                if ((coupling_Array_dev(4,tau).ne.0).and.(coupling_Array_dev(4,taup).ne.0)) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RD(mnuc_dev,tau,taup,q2,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,7,2)
                end if
                
                !RS1D (c)
                ! c5,c4 or c8,c9
                if (((coupling_Array_dev(4,tau).ne.0).and.(coupling_Array_dev(3,taup).ne.0)).or. &
                        ((coupling_Array_dev(7,tau).ne.0).and.(coupling_Array_dev(8,taup).ne.0))) then
                    p_tot = p_tot + 1./mnuc_dev**2 * RS1D(tau,taup,c,j_chi,coupling_Array_dev) * sumW(w,vesc,i,tau,taup,8,1)
                end if
        end if
        !     end do
        ! end do
        p_tot = p_tot  * hbar_dev**2 * c0_dev**2
    end function p_tot


    !   this is eqn 2.1 in 1501.03729
    attributes(device) function OMEGA_oper(rindex,w)
        double precision :: w, vesc,mu,muplus,u,J, OMEGA_oper
        integer rindex, i
        
        vesc = tab_vesc(rindex)
        u = sqrt(w**2-vesc**2)
        Omega_oper = 0.d0
        
        ! if only one isotope is being run
        if ((0 .lt. pickIsotope).and.(pickIsotope .lt. 17)) then
            J = AtomicSpin_oper(pickIsotope)
            mu = mdm_dev/mnuc_dev/AtomicNumber_oper_dev(pickIsotope)
            muplus = (1.+mu)/2.
            if (mu/muplus**2 .gt. u**2/w**2) then
                ! note the factor [(2*mnuc_dev*AtomicNumber(pickIsotope))/(w**2*(2*J+1))],is the product of
                ! the factors in front of eqns 3.26 and 3.23 in paper 1501.03729
                Omega_oper = Omega_oper+(NAvo*tab_starrho(rindex)/AtomicNumber_oper_dev(pickIsotope)/mnuc_dev* &
                    tab_mfr_oper(rindex,pickIsotope))*w*((2*mnuc_dev*AtomicNumber_oper_dev(pickIsotope))/(w**2*(2*J+1)))* &
                    p_tot(w,vesc,pickIsotope)
            end if
        else ! if all the isotopes are being run
            do i = 1,niso
                J = AtomicSpin_oper(i)
                mu = mdm_dev/mnuc_dev/AtomicNumber_oper_dev(i)
                muplus = (1.+mu)/2.
                if (mu/muplus**2 .gt. u**2/w**2) then
                    ! note the factor [(2*mnuc_dev*AtomicNumber(i))/(w**2*(2*J+1))],is the product of
                    ! the factors in front of eqns 3.26 and 3.23 in paper 1501.03729
                    Omega_oper = Omega_oper+(NAvo*tab_starrho(rindex)/AtomicNumber_oper_dev(i)/mnuc_dev*tab_mfr_oper(rindex,i))*w* &
                                            ((2*mnuc_dev*AtomicNumber_oper_dev(i))/(w**2*(2*J+1)))*p_tot(w,vesc,i)
                end if
            end do
        end if
    end function OMEGA_oper
end module gpumod

! If integrand_oper can be overloaded to have a device and host version that just use a different module,
! then I can just use both opermod and gpumod in the captnoper integrand? Otherwise need to create device
! sepcific integrand_oper_dev
! !get_vdist might need device label
! !   this is the integral over R in eqn 2.3 in 1501.03729
! attributes(device) function integrand_oper_dev(u)
!     use gpumod
!     implicit none
!     double precision :: u, w, vesc, integrand_oper, int
!     vesc = tab_vesc(ri_for_omega)
!     w = sqrt(u**2+vesc**2)
!     int = get_vdist(u)/u*w*Omega_oper_dev(ri_for_omega,w)
!     integrand_oper = int
! end function integrand_oper_dev


!   Need to pass all the operators into the subroutine
attributes(global) subroutine captn_oper_kernel(mx_in, jx_in, niso_in, isotopeChosen, capped)
    ! use gpumod
    implicit none
    integer, intent(in), value :: niso_in, isotopeChosen
    integer i, ri
    double precision, intent(in), value :: mx_in, jx_in
    double precision, value :: capped !this is the output
     ! array of coupling constants
    double precision :: eps_devabs, eps_devrel,limit,result,abserr,neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),rlist(1000)!for integrator
    !external gausstest !this is just for testing
    external integrand_oper
    external dummyf
    eps_devabs=1.d-17
    eps_devrel=1.d-17
    limit=1000

    mdm_dev = mx_in
    j_chi = jx_in
    niso = niso_in
    
    pickIsotope = isotopeChosen

    ! temporary, the user will want to choose their coupling constants to match a model
    !                           c1,  c3,  c4, c5,   c6,  c7,  c8,  c9, c10, c11, c12, c13, c14, c15   
    !coupling_Array_dev = reshape((/1.65d-8, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
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
            eps_devabs,eps_devrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        u_int_res(ri) = result
        capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)
    end do

    !completes integral (2.4) of paper 1501.03729
    capped = 4.d0*pi*Rsun**3*capped

    if (capped .gt. 1.d100) then
      print*,"Capt'n General says: Oh my, it looks like you are capturing an  &
      infinite amount of dark matter in the Sun. Best to look into that."
    end if
end subroutine captn_oper_kernel

subroutine captn_oper_gpu(mx_in, jx_in, niso_in, isotopeChosen, capped)
    implicit none
    integer :: niso_in, isotopeChosen
    double precision :: mx_in, jx_in
    double precision :: capped !this is the output
    ! dim3 variables to define the grid and block shapes
    type(dim3) :: dimGrid, dimBlock
    integer :: r

    dimGrid = dim3( 1, 1, 1 )
    dimBlock = dim3( 2, 2, 7 )
    call captn_oper_kernel<<<dimGrid,dimBlock>>>( mx_in, jx_in, niso_in, isotopeChosen, capped )
    r = cudathreadsynchronize()
end subroutine captn_oper_gpu