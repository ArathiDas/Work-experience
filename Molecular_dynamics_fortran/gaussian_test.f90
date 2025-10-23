program gaussian_test
    implicit none
    integer,parameter      :: N = 10000
    real*8                 :: values(N), v(N), kb, T, m 
    integer                :: i 
    real*8                 :: gaussian_distr 
    
    ! Initialize 
    kb = 1.38d-23
    T = 293.0d0
    m = 40*1822.8884850
    
    ! Generate Gaussian values 
    do i = 1,N 
       values(i) = gaussian_distr()
    end do
    open(10, file='gaussian_data.dat')
       do i = 1,N
          write(10,*) values(i)
       end do
    close(10)
    
    ! Generate Velocities
    do i = 1,N 
       v(i) = sqrt(kb*T/m) * gaussian_distr()
    end do
    open(20, file='velocity_data.dat')
    do i = 1,N
       write(20,*) v(i)
    end do
    close(10)
       
end program 

 real*8 function gaussian_distr()
 implicit none
 integer, save :: iset=0
 real*8 :: fac,rsq,var_1,var_2,x,y
 real*8, save :: gset
 call random_seed()
 if (iset==0) then
 ! first time to calculate rsq
 call random_number(x)
 call random_number(y)
 var_1 = 2.d0*x-1.d0
 var_2 = 2.d0*y-1.d0
 rsq = var_1**2+var_2**2
 do while (rsq.ge.1..or.rsq.eq.0.) ! repeat the operation until rsq is between 0 and 1
 call random_number(x)
 call random_number(y)
 var_1 = 2.d0*x-1.d0
 var_2 = 2.d0*y-1.d0
 rsq = var_1**2+var_2**2
 enddo
 fac = sqrt(-2.d0*log(rsq)/rsq)
 gset = var_1*fac
 gaussian_distr = var_2*fac
 iset = 1
 else
 endif
 gaussian_distr = gset
 iset = 0
 end function gaussian_distr   
