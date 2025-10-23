!-----------------!
 module parameters
!-----------------!

 real(8), parameter :: length=5.12d0   ! length of the box (in Bohr)
 real(8), parameter :: mass = 1822.88839d0 != 1 atomic mass unit (=1 g/mol)
 real(8), parameter :: pi=3.141592653589793d0
 real(8), parameter :: au2kcalmol=627.509d0, fs2au=41.341373336561d0, au2invcm=219474.63d0
 real(8) :: angfreq, barrier
 character(10) :: potentialtype

 end module parameters
!---------------------!

!-----------------!
 program propagate
!-----------------!
 use parameters
 implicit none

 integer :: npoints,ntime,snapshot,i
 real(8) :: alpha,dt,t,dx,x0,norm,energy,k,p,S
 complex(8)   :: C,CE
 real(8), allocatable :: pot(:),kin(:)
 complex(8), allocatable :: psi(:),psi0(:),exppot(:),expkin(:)

 open(unit=10,file='wavepacket')
   read(10,*) npoints               !Number of lattice points
   read(10,*) x0                    !Initial position
   read(10,*) alpha                 !Governs the initial width of the wave packet
   read(10,*) dt                    !Propagation time step
   read(10,*) ntime                 !Number of propagation steps
   read(10,*) snapshot              !snapshot frequency 
 close(10)

 open(unit=11,file='potential')
   read(11,*) potentialtype         !harmonic or double well potential
   read(11,*) angfreq               !Angular frequency for harmonic potential (in fs^(-1))
   read(11,*) barrier               !Height of barrier in double well potential (in kcal/mol)
 close(11)

 dt=dt*fs2au                        !convert femtoseconds to atomic units
 angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

 allocate(psi(npoints),psi0(npoints))
 allocate(pot(npoints),exppot(npoints))
 allocate(kin(npoints),expkin(npoints))

 dx=length/dble(npoints)

 call initpsi(npoints,dx,alpha,x0,psi0)              !Obtain initial wavepacket psi0
 call fourier(0,npoints,psi0)                        !Initialize the FFT
 call operators(npoints,dx,dt,pot,kin,exppot,expkin) !Calculate the kinetic and potential operators
 
 open(1,file='energy1.dat')
 open(2,file='survival.dat')
 open(3,file='spectrum.dat')

 psi=psi0                                            !Set the wavepacket psi at t=0 equal psi0
 do i=1,ntime                                        !Start propagation
    t=i*dt
    psi=psi*exppot                                   !Multiply psi with exp(-i*dt*potential)
    call fourier(1,npoints,psi)                      !Forward FFT to momentum space
    psi=psi*expkin                                   !Multiply psi with the exp(-i*dt*kinetic operator)
    call fourier(-1,npoints,psi)                     !Backward FFT to position space
    call calcnorm(npoints,dx,psi,norm)
    call calcenergy(npoints,dx,dt,pot,kin,psi,energy,k,p)
    write(1,*) i, p, k, energy

 call autocorrelation(npoints,dx,psi,psi0,C,S)
 write(2,*) i*0.01, C, S
 end do                                            !End propagation
close(1)
close(2)

 call spectrum(ntime,C,CE)
 write(3,*) t,CE

 deallocate(psi,psi0)
 deallocate(pot,exppot)
 deallocate(kin,expkin)

 end program propagate
!---------------------!

!--------------------------------------------!
 subroutine initpsi(npoints,dx,alpha,x0,psi0)
!--------------------------------------------!
 implicit none

 integer :: i,j,npoints
 real(8) :: alpha,x,x0,dx,norm
 complex(8) :: psi0(npoints)
 
 norm=0.d0
 do i=-npoints/2+1,npoints/2
    x=dble(i)*dx
    if (i>0) then
       j=i
    else     
       j=i+npoints
    endif
    psi0(j)=exp(-alpha*(x-x0)**2)
    norm=norm+abs(psi0(j))**2*dx
 end do
 norm=1.d0/dsqrt(norm)
 do i=1,npoints
    psi0(i)=psi0(i)*norm
 end do

 end subroutine initpsi
!----------------------!

!---------------------------------------------------------!
 subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!---------------------------------------------------------!
 use parameters
 implicit none

 integer :: i,j,npoints
 real(8) :: x,p,b,dt,dx,dp,pot(npoints),kin(npoints)
 complex(8) :: exppot(npoints),expkin(npoints)

 dp=2.d0*pi/length
 do i=-npoints/2+1,npoints/2
    x=dble(i)*dx
    p=dble(i-1)*dp
    if (i>0) then
       j=i
    else
       j=i+npoints
    endif
    if (potentialtype=='harmonic') then
       pot(j)=0.5d0*mass*angfreq**2*x**2
    elseif (potentialtype=='doublewell') then
       pot(j)=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)/au2kcalmol
    endif
    kin(j)=0.5d0*p**2/mass
    exppot(j)=exp(-dt*(0,1)*pot(j))
    expkin(j)=exp(-dt*(0,1)*kin(j))
 end do

 end subroutine operators
!------------------------!

!-----------------------------------!
 subroutine fourier(dir,npoints,psi)
!-----------------------------------!
 implicit none

 integer :: i,npoints,dir
 real(8) :: nr
 complex(8) :: psi(npoints)
 real(8), allocatable, save :: wsave(:)       

 if (dir==1) then
    call dcfftf(npoints,psi,wsave)
    nr=1.d0/sqrt(dble(npoints))
    do i=1,npoints
       psi(i)=psi(i)*nr
    end do
 elseif (dir==-1) then
    call dcfftb(npoints,psi,wsave)
    nr=1.d0/sqrt(dble(npoints))
    do i=1,npoints
       psi(i)=psi(i)*nr
    end do
 elseif (dir==0) then
    if (allocated(wsave)) deallocate(wsave)
    allocate(wsave(4*npoints+20))
    call dcffti(npoints,wsave)
 endif

 end subroutine fourier
!----------------------!

!---------------------!
 subroutine calcnorm(npoints,dx,psi,norm)
!---------------------!
 use parameters
 implicit none
 integer     :: i,npoints
 real(8)     :: dx
 real(8),intent(out)  :: norm
 complex(8)     :: psi(npoints)

 norm = 0.0
 do i = 1,npoints 
  norm = norm + (psi(i)*conjg(psi(i))*dx)
 end do

 end subroutine calcnorm
!-----------------------!

!-----------------------!
 subroutine calcenergy(npoints,dx,dt,pot,kin,psi,energy,k,p)
!-----------------------!
 use parameters
 implicit none
 integer     :: i,npoints
 real(8)     :: dx,pot(npoints),kin(npoints),dt,energy,k,p
 complex(8)     :: psi(npoints),l(npoints)
 
 energy = 0.0
 call fourier(1,npoints,psi)
 l = kin*psi
 call fourier(-1,npoints,l)
 call fourier(-1,npoints,psi)

 do i = 1,npoints
  energy = energy+(conjg(psi(i))*l(i)*dx)
 end do
  k = energy
  p = 0.0
  do i = 1,npoints 
   p = p+(conjg(psi(i))*pot(i)*psi(i)*dx) 
   energy = energy+(conjg(psi(i))*pot(i)*psi(i)*dx)
 end do
 end subroutine calcenergy
!-------------------------!

!----------------------------!
 subroutine autocorrelation(npoints,dx,psi,psi0,C,S)
!----------------------------!
 use parameters
 implicit none
 integer      :: i,npoints
 real(8)      :: dx, S
 complex(8)   :: psi(npoints),psi0(npoints),C

 C = 0.d0
 S = 0.d0
 do i = 1,npoints
    C = C + (conjg(psi0(i))*psi(i)*dx)
    S = abs(C)**2
 end do
 
 end subroutine autocorrelation
!------------------------------!

!---------------------!
 subroutine spectrum(ntime,C,CE)
!---------------------!
 use parameters
 implicit none
 integer     :: i,ntime
 complex(8)     :: C,CE

 CE = C
 call fourier(0,ntime,C)
 call fourier(-1,ntime,CE)

 end subroutine spectrum
!-----------------------!
