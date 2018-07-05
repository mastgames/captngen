!   Capt'n General
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


    module capmod
    implicit none
    double precision, parameter :: pi=3.141592653, NAvo=6.0221409d23,GMoverR = 1.908e15
    double precision, parameter :: c0 = 2.99792458d10, mnuc = 0.938, q0 = 0.04,v0 = 220.d5
    double precision, parameter :: eps = 1d-10 !stops divisions by zero
    !these are now set in captn_init
    double precision :: usun , u0 ,rho0,vesc_halo, Rsun
    !this goes with the Serenelli table format
    double precision :: AtomicNumber(29) !29 is is the number from the Serenelli files; if you have fewer it shouldn't matter
    double precision, allocatable :: tab_mencl(:),tab_starrho(:),tab_mfr(:,:), tab_r(:), tab_vesc(:), tab_dr(:),tab_T(:),tab_g(:)

    ! nq and nv can be -1, 0, 1, 2; this is set in the main program
    integer :: nq, nv, niso, ri_for_omega, nlines
    double precision :: mdm, sigma_0

    contains

    !velocity distribution,
    function get_vdist(u)
    double precision :: u,get_vdist, f, normfact
    f = (3./2.)**(3./2.)*4.*rho0*u**2/sqrt(pi)/mdm/u0**3 &
    *exp(-3.*(usun**2+u**2)/(2.*u0**2))*sinh(3.*u*usun/u0**2)/(3.*u*usun/u0**2)
!    normfact = .5*erf(sqrt(3./2.)*(vesc_halo-usun)/u0) + &
!    .5*erf(sqrt(3./2.)*(vesc_halo+usun)/u0)+ u0/(sqrt(6.*pi)*usun) &
!    *(exp(-3.*(usun+vesc_halo)/2./u0**2)-exp(-3.*(usun-vesc_halo)/2./u0**2))
    normfact = 1.
!    print*,normfact
    f = f/normfact
    get_vdist=f
    end function get_vdist


     !generalized form factor: hydrogen
    function GFFI_H(w,vesc)
    double precision :: p, mu,w,vesc,u,muplus,GFFI_H,G
    p = mdm*w
    mu = mdm/mnuc
    muplus = (1.+mu)/2.
    u = sqrt(w**2-vesc**2)
    if (nq .ne. -1) then
    G = (p/q0/c0)**(2.d0*nq)*mdm*w**2/(2.d0*mu**nq)*1./(1.+nq)*((mu/muplus**2)**(nq+1)-(u**2/w**2)**(nq+1))
    else
      !eps added to make the log finite: the value of eps does not affect the answer
    G = ((p)/q0/c0)**(2.d0*nq)*mdm*w**2/(2.d0*mu**nq)*log(mu/muplus**2*w**2/(u+eps)**2)
    endif
    GFFI_H = G
    end function GFFI_H

    !generalized form factor: other elements
    function GFFI_A(w,vesc,A)
    double precision :: p, mu,w,vesc,u,muplus,mN,A,Ei,B
    double precision :: dgamic,GFFI_A
    p = mdm*w
    mu = mdm/mnuc/A
    muplus = (1.+mu)/2.
    u = sqrt(w**2-vesc**2)
    mN = A*mnuc
    Ei  = 5.8407d-2/(mN*(0.91*mN**(1./3.)+0.3)**2)
    B = .5*mdm*w**2/Ei/c0**2
    if (nq .eq. 0) then
    GFFI_A = Ei*c0**2*(exp(-mdm*u**2/2/Ei/c0**2)-exp(-B*mu/muplus**2))
    else
    GFFI_A = ((p+eps)/q0/c0)**(2*nq)*Ei*c0**2/(B*mu)**nq*(dgamic(1.+dble(nq),B*u**2/w**2+eps) &
    - dgamic(1.+dble(nq),B*mu/muplus**2+eps))
    end if


    end function GFFI_A


    !this is omega/sigma_0
    function OMEGA(rindex,w)
    double precision :: sigma_N, GF,vesc,Omega,mu,muplus,muminus,u,w
    integer i, rindex

    vesc = tab_vesc(rindex)
    u = sqrt(w**2-vesc**2)
    Omega = 0.d0

    do i = 1,niso
    !this is fine for SD as long as it's just hydrogen. Otherwise, spins must be added
    sigma_N = AtomicNumber(i)**4*(mdm+mnuc)**2/(mdm+AtomicNumber(i)*mnuc)**2

    !hydrogen
    mu = mdm/mnuc/AtomicNumber(i)
    muplus = (1.+mu)/2.
    muminus = (mu-1.d0)/2.
    if (mu*vesc**2/muminus**2 .gt. u**2) then
      if (i .eq. 1) then
        GF = GFFI_H(w,vesc)
        else
        GF = GFFI_A(w,vesc,AtomicNumber(i))
      end if
      Omega = Omega+sigma_N*NAvo*tab_starrho(rindex)/AtomicNumber(i)/mnuc*tab_mfr(rindex,i)*muplus**2/mu*GF
    end if
    end do
    Omega = Omega*2.d0/mdm/w

    end function omega


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !read in solar parameters from Aldo Serenelli-style files, with header removed
    subroutine get_solar_params(filename,nlines)
    character*300 :: filename
    double precision :: Temp, Pres, Lumi !these aren't used, but dummies are required
    double precision, allocatable :: phi(:) !this is used briefly
    integer :: i,j, nlines,iostatus
    !Get number of lines in the file

    Rsun = 69.57d9 !this is set here, for other stars, this sub is not called

    open(99,file=filename)
    nlines=0
    do
    read(99,*, iostat=iostatus)
    if(iostatus/=0) then ! to avoid end of file error.
    exit
    else
    nlines=nlines+1
    end if
    end do
    close(99)
    nlines = nlines -1
    !allocate the arrays
    allocate(tab_mencl(nlines))
    allocate(tab_r(nlines))
    allocate(tab_starrho(nlines))
    allocate(tab_mfr(nlines,29)) !we could just allocate niso, but this leads to problems
    allocate(tab_vesc(nlines))
    allocate(phi(nlines))
    allocate(tab_dr(nlines))
    allocate(tab_T(nlines)) !not used in capgen; used for transgen (and anngen? )
    allocate(tab_g(nlines))


    !now actually read in the file
    open(99,file=filename)
    do i=1,nlines
    read(99,*) tab_mencl(i),tab_r(i), tab_T(i), tab_starrho(i), Pres, Lumi, tab_mfr(i,:)
    end do
    close(99)

    !we calculate the escape velocity here since all the ingredients are ready
    phi(nlines) = -GMoverR
    tab_vesc(nlines) = sqrt(-2.d0*phi(nlines))
    tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)
    do i = 1,nlines-1
    j = nlines-i !trapezoid integral
    phi(j) = phi(j+1) + GMoverR*(tab_r(j)-tab_r(j+1))/2.*(tab_mencl(j)/tab_r(j)**2+tab_mencl(j+1)/tab_r(j+1)**2)
    tab_vesc(j) = sqrt(-2.d0*phi(j)) !escape velocity in cm/s
    tab_dr(j) = -tab_r(j)+tab_r(j+1) !while we're here, populate dr
    tab_g(j) = -(-phi(j)+phi(j+1))/tab_dr(j)
    end do
    tab_g(nlines) = tab_g(nlines-1)

    open(55,file = "tab_serenelli.dat")
    do i=1,nlines
      write(55,*) tab_r(i), tab_starrho(i), tab_vesc(i), tab_mfr(i,:)
      end do
      close(55)

      ! Populate the atomic number tables here (because it relies on a specific format)
    AtomicNumber  = (/ 1., 4., 3., 12., 13., 14., 15., 16., 17., &
                      18., 20.2, 22.99, 24.3, 26.97, 28.1, 30.97,32.06, 35.45, &
                      39.948, 39.098, 40.08, 44.95, 47.86, 50.94, 51.99, &
                      54.93, 55.845, 58.933, 58.693/)


    return
    end

! this is to make sure the integrator does what it's supposed to
    function gaussinmod(x)

    double precision :: x,gaussinmod
    gaussinmod = nq*exp(-x**2/2.d0)
    end function gaussinmod



!end get_solar_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FUNCTIONS FOR USE BY TRANSGEN FOR DM TRANSPORT





    end module capmod



!   Some functions that have to be external, because of the integrator.
    function gausstest(x) !just a test for the integrator. Nothing to see here
    use capmod
    double precision :: x,gausstest
    gausstest = gaussinmod(x)
    end function gausstest


    !THIS IS THE IMPORTANT FUNCTION: the integrand for the integral over u
    function integrand(u)
    use capmod
    double precision :: u, w, vesc, integrand, int
    vesc = tab_vesc(ri_for_omega)
    w = sqrt(u**2+vesc**2)
!    print*, NAvo, tab_starrho(ri_for_omega), mnuc, tab_mfr(ri_for_omega,1) !, Omega(ri_for_omega,w)
    int = get_vdist(u)/u*w*Omega(ri_for_omega,w)
    if (nv .ne. 0) then
    int = int*(w/v0)**(2*nv)
    end if
    integrand = int


    ! print*,"integrand", u, vesc,get_vdist(u), w, ri_for_omega, Omega(ri_for_omega,w),int
    ! open(55,file = "mesabullshit.dat", status="old", position="append", action="write")
    ! do i=1,nlines
    !   write(55,*) ri_for_omega, u,w, vesc, integrand
    !   end do
    ! close(55)


    end function integrand


    function dummyf(x)
    double precision :: x, dummyf
    dummyf = 1.d0
    end function dummyf


    ! This is for use with MESA arrays specifically (though could likely be used for GARSTEC)
    ! star initialization done elsewhere
    ! subroutine captn_mesa(mx_in,sigma_0_in,niso_in,nq_in,nv_in,capped)
    !   subroutine captn_mesa()
    !   use capmod
    !   implicit none
    !   ! integer, intent(in):: nq_in, niso_in, nv_in
    !   integer i, ri
    !   ! double precision, intent(in) :: mx_in, sigma_0_in
    !   double precision :: capped, maxcap !this is the output
    !   double precision :: epsabs, epsrel,limit,result,abserr,neval !for integrator
    !   double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    !   double precision, allocatable :: u_int_res(:)
    !
    !
    ! end subroutine captn_mesa


    subroutine captn_general(mx_in,sigma_0_in,niso_in,nq_in,nv_in,capped)
    use capmod
    implicit none
    integer, intent(in):: nq_in, niso_in, nv_in
    integer i, ri
    double precision, intent(in) :: mx_in, sigma_0_in
    double precision :: capped, maxcap !this is the output
    double precision :: epsabs, epsrel,limit,result,abserr,neval !for integrator
    double precision :: ier,alist,blist,rlist,elist,iord,last!for integrator
    double precision, allocatable :: u_int_res(:)

    dimension alist(1000),blist(1000),elist(1000),iord(1000),   rlist(1000)!for integrator
    external gausstest !this is just for testing
    external integrand
    external dummyf
    epsabs=1.d-17
    epsrel=1.d-17
    limit=1000


    mdm = mx_in
    sigma_0 = sigma_0_in
    niso = niso_in
    nq = nq_in
    nv = nv_in



    if (nq*nv .ne. 0) then
    print*, "Oh no! nq and nv can't both be nonzero. "
    return
    end if

    if (.not. allocated(tab_r)) then !
        print*,"Errorface of errors: you haven't called captn_init to load the solar model!"
        return
    end if
        allocate(u_int_res(nlines))

!    As far as I can tell, the second argument (fofuoveru) does nothing in this integrator. I've sent it to an inert dummy just in case.
    capped = 0.d0

    do ri=1,nlines !loop over the star
    result = 0.d0
    ri_for_omega = ri !accessed via the module
    !call integrator
    call dsntdqagse(integrand,dummyf,1.d0,vesc_halo/3., &
    epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

    u_int_res(ri) = result*sigma_0
    ! print*,ier
    capped = capped + tab_r(ri)**2*u_int_res(ri)*tab_dr(ri)

    if (isnan(capped)) then
      capped = 0.d0
      stop 'hammer time'
    end if

    !   print*, ri, tab_r(ri), result, ier
    !   stop "nanintedfsfs"
    ! end if

    ! print*,"Capgen capped ", capped, tab_r(ri), u_int_res(ri), tab_dr(ri)
    end do
    capped = 4.d0*pi*Rsun**3*capped

    if (capped .gt. 1.d100) then
      print*,"Capt'n General says: Oh my, it looks like you are capturing an  &
      infinite amount of dark matter in the Sun. Best to look into that."
    end if



    ! return
  ! end if

    !this now has its own function:
    ! maxcap = pi/3.d0*rho0/mdm*Rsun**2 &
    ! *(exp(-3./2.*usun**2/u0**2)*sqrt(6.d0/pi)*u0 &
    ! + (6.d0*GMoverR/usun + (u0**2 + 3.d0*usun**2)/usun)*erf(sqrt(3./2.)*usun/u0))

  !  print*,"sigma_0 =", sigma_0, "; m = ", mdm, "; nq = ", nq, "; Capture rate: ", capped, "max = ", maxcap

    ! if (capped .gt. maxcap) then
    !     capped = maxcap
    ! end if


    end subroutine captn_general
  !This is fine as long as the escape velocity is large enough
  subroutine captn_maxcap(mwimp_in,maxcap)
    use capmod
    implicit none
    double precision maxcap
    double precision, intent(in) :: mwimp_in
    mdm = mwimp_in
    maxcap = pi/3.d0*rho0/mdm*Rsun**2 &
    *(exp(-3./2.*usun**2/u0**2)*sqrt(6.d0/pi)*u0 &
    + (6.d0*GMoverR/usun + (u0**2 + 3.d0*usun**2)/usun)*erf(sqrt(3./2.)*usun/u0))
  end subroutine captn_maxcap



    !! captn_specific calculates the capture rate for constant cross section.
    subroutine captn_specific(mx_in,sigma_0_SD_in,sigma_0_SI_in,capped_SD,capped_SI)
      implicit none
    double precision, intent(in) :: mx_in, sigma_0_SD_in,sigma_0_SI_in
    double precision :: capped_SD,capped_SI

    call captn_general(mx_in,sigma_0_SD_in,1,0,0,capped_SD)
    call captn_general(mx_in,sigma_0_SI_in,29,0,0,capped_SI)

    end subroutine captn_specific

!------!------!------!------!------INITIALIZATION FCT

    subroutine captn_init(solarmodel,rho0_in,usun_in,u0_in,vesc_in)
      !input velocities in km/s, not cm/s!!!
    use capmod
    use iso_c_binding, only: c_ptr
    implicit none
    character (len=300) solarmodel
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in
!    common solarmodel
!    external solarmodel


    if  (.not. allocated(tab_r)) then !
        print*,"Capgen initializing from model: ",solarmodel
        call get_solar_params(solarmodel,nlines)
    end if
    ! print*,"Capgen tabulons already allocated, you might be overdoing it by calling the init function more than once."
    usun = usun_in*1.d5
    u0 =  u0_in*1.d5
    rho0 =rho0_in
    vesc_halo = vesc_in*1.d5

    end subroutine captn_init

!     subroutine captn_init_mesa(rho0_in,usun_in,u0_in,vesc_in,starrho)
!       !input velocities in km/s, not cm/s!!!
!     use capmod
!     use iso_c_binding, only: c_ptr
!     implicit none
!     double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in
! !    common solarmodel
! !    external solarmodel
!
!
!     if  (.not. allocated(tab_r)) then !
!         print*,"Capgen initializing from model: ",solarmodel
!         call get_solar_params(solarmodel,nlines)
!     end if
!     ! print*,"Capgen tabulons already allocated, you might be overdoing it by calling the init function more than once."
!     usun = usun_in*1.d5
!     u0 =  u0_in*1.d5
!     rho0 =rho0_in
!     vesc_halo = vesc_in*1.d5
!
!   end subroutine captn_init_mesa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For mesa interface only: allocate arrays.
  subroutine allocate_stellar_arrays(nlines_mesa)
    use capmod
    integer, intent(in) :: nlines_mesa
    nlines = nlines_mesa
    allocate(tab_mencl(nlines))       !M(<r)
    allocate(tab_r(nlines))           !r
    allocate(tab_starrho(nlines))     !rho
    allocate(tab_mfr(nlines,8))       !mass fraction per isotope
    allocate(tab_vesc(nlines))        !local escape velocity
    allocate(tab_T(nlines))           !temperature
    ! allocate(phi(nlines)) !! <--- not needed; computed in wimp_support.f
    allocate(tab_dr(nlines))          !dr (nice)
    allocate(tab_g(nlines))           !local gravitational acceleration, needed for transport

    RETURN
  end subroutine allocate_stellar_arrays

  subroutine deallocate_stellar_arrays()
    use capmod
    deallocate(tab_mencl)
    deallocate(tab_r)
    deallocate(tab_starrho)
    deallocate(tab_mfr) !we could just allocate niso, but this leads to problems
    deallocate(tab_vesc)
    deallocate(tab_T)
    deallocate(tab_dr)
    deallocate(tab_g)
    RETURN
  end subroutine deallocate_stellar_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This is called INSTEAD of get_solar_params, for use with MESA interface.
  subroutine get_stellar_params(rmesa,rhomesa,mfrmesa,mesavesc,Tmesa,mesag,mesamass,mesaradius,rho0_in,usun_in,u0_in,vesc_in)
    use capmod
    double precision :: mesamass, mesaradius
    double precision :: rhomesa(nlines), rmesa(nlines), mfrmesa(8,nlines),mesavesc(nlines),mesag(nlines)
    double precision :: Tmesa(nlines)
    integer i
    double precision,intent(in) :: rho0_in,usun_in,u0_in,vesc_in

    usun = usun_in*1.d5
    u0 =  u0_in*1.d5
    rho0 =rho0_in
    vesc_halo = vesc_in*1.d5

    Rsun = rmesa(nlines)
    tab_r = rmesa/Rsun
    tab_starrho = rhomesa
    tab_vesc = mesavesc
    tab_T = tmesa
    tab_g = mesag
    do i= 1,8
    tab_mfr(:,i) = mfrmesa(i,:)
  end do

  ! tab_r(1) = tab_r(2)/10. !just to avoid some awkward?

  do i = 1, nlines-1
    tab_dr(i) = -tab_r(i)+tab_r(i+1) !while we're here, populate dr
  end do
  tab_dr(nlines) = tab_r(nlines)-tab_r(nlines-1)

  AtomicNumber(1:8) = (/1., 3., 4., 12., 14., 16., 20., 28.  /)


    ! print*,tab_mfr(1,:)
    ! %this is where the tables get populated
    ! print*,tab_r(nlines-1), tab_r(2)
    ! print*,"Mass of guy ", mesamass
    ! print*,"Radius of guy ",mesaradius
    ! open(55,file = "mesatables.dat")
    ! do i=1,nlines
    !   write(55,*) tab_r(i), tab_starrho(i), tab_vesc(i), tab_mfr(i,:)
    !   end do
    !   close(55)

    RETURN
  end subroutine get_stellar_params


  subroutine getnlines(nlines_out) !a little auxiliary trick
    use capmod
      integer, intent(out) :: nlines_out
      nlines_out = nlines
      return
    end
