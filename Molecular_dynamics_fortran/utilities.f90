module utilities

public::lecture_xyz             ! Read and save positions in an .xyz file
public::get_distances           ! Calculate square of distances between atoms
public::get_forces              ! Calculate forces from the positions
public::get_pot_ener            ! Calculate energy from the positions
public::output_structure        ! Write output structures 
public::output_energy           ! Write output energies
public::get_kin_ener            ! Calculate kinetic energy from velocities
public::get_T                   ! Calculate T from kinetic energy
public:: gaussian_distr
contains

subroutine lecture_xyz (filename,tableau,nbr_atomes)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Record frames from an .xyz file   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

character*100, intent(in)  :: filename
integer,       intent(out) :: nbr_atomes

real*8, allocatable, dimension (:,:),intent(out) :: tableau

character*10 :: var_poub
character*10 :: sep
integer      :: OK
integer      :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1,file=filename,iostat=OK)

read(1,*,iostat=OK), nbr_atomes
read(1,*,iostat=OK), sep

allocate(tableau(nbr_atomes,3))

do i=1,nbr_atomes
  read(1,*,iostat=OK), var_poub, tableau(i,1), tableau(i,2), tableau(i,3)
enddo

close(1)

endsubroutine lecture_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate square of distances between atoms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_distances (positions,distances,nbr_atomes)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: positions
real*8, allocatable, dimension (:,:),intent(inout) :: distances

integer, intent(in)  :: nbr_atomes
integer              :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,nbr_atomes
  do j=1,nbr_atomes
   distances(i,j)=      (positions(i,1)-positions(j,1))**2 + &
					&   (positions(i,2)-positions(j,2))**2 + &
					&	(positions(i,3)-positions(j,3))**2
  enddo
enddo

!do i=1,nbr_atomes-1
!  do j=i+1,nbr_atomes
!   distances(j,i)=distances(i,j)
!  enddo
!enddo

endsubroutine get_distances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ forces               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_forces (forces,distances,positions,nbr_atomes,sigma,eps)

implicit none

real*8, allocatable, dimension (:,:),intent(inout) :: forces
real*8, allocatable, dimension (:,:),intent(in)    :: distances
real*8, allocatable, dimension (:,:),intent(in)    :: positions

real*8,intent(in)  :: sigma,eps
real*8             :: sigma_12,sigma_6
real*8             :: dist2
real*8             :: dr4,dr8,dr14
real*8             :: part_6,part_12
integer,intent(in) :: nbr_atomes
integer            :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,nbr_atomes
  do j=1,3
    forces(i,j)=0.0d0
  enddo
enddo

sigma_6=sigma**6
sigma_12=sigma**12

do i=1,nbr_atomes-1
  do k=i+1,nbr_atomes
	do j=1,3
	  dist2=(positions(i,j)-positions(k,j))
	  dr4=distances(i,k)*distances(i,k)
	  dr8=dr4*dr4
	  dr14=dr8*dr4*distances(i,k)
	  part_6  = sigma_6*6.0d0*(-(dist2)/dr8)
	  part_12 = sigma_12*6.d0*(-(2.0d0*dist2)/dr14) 
	  forces(i,j)=forces(i,j)-4.0d0*eps*(part_12-part_6)
	  forces(k,j)=forces(k,j)+4.0d0*eps*(part_12-part_6)
	enddo
  enddo
enddo

endsubroutine get_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ Energy               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_pot_ener (distances,nbr_atomes,sigma,eps,scp_energy,energy)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: distances

real*8,intent(in)   :: sigma,eps,scp_energy
real*8              :: sigma_6,sigma_12
real*8              :: part_6,part_12
real*8              :: dr6,dr12
real*8,intent(inout):: energy
integer,intent(in)  :: nbr_atomes
integer             :: i,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

energy=0.0d0

sigma_6=sigma**6
sigma_12=sigma**12

do i=1,nbr_atomes-1
    do k=i+1,nbr_atomes
	  dr6=distances(i,k)*distances(i,k)*distances(i,k)
	  dr12=dr6*dr6
	  part_6  = sigma_6/dr6
	  part_12 = sigma_12/dr12
      energy=energy+4*eps*(part_12-part_6)
	enddo
enddo

endsubroutine get_pot_ener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output gemetry                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_structure (index,positions,nbr_atomes,energy_pot)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: positions

real*8,intent(in)   :: energy_pot
integer,intent(in)  :: index
integer,intent(in)  :: nbr_atomes
integer             :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(index,*), nbr_atomes
write(index,*), energy_pot
do i=1,nbr_atomes
  write(index,*), "Ar", positions(i,1), positions(i,2), positions(i,3)
enddo

endsubroutine output_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ouput Energies                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_energy (index,index_loop,timestep,energy_pot,energy_kin,temp)

implicit none

real*8, parameter   :: fs_to_timeau=41.34137314d0
real*8, parameter   :: Hartree_to_kcal=627.50947d0
real*8,intent(in)   :: energy_pot,energy_kin,temp
integer,intent(in)  :: index,index_loop
real*8,intent(in)   :: timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (index_loop==1) then
  write(index,*), "             Timestep            ", "Temp         ", "    Potential Energy      ", "Kinetic Energy      ", &
                                                                                              & "   Total Energy " 
  write(index,'(f20.4,3x,f20.10,3x,2f20.10,3x,2f20.10)'), timestep*index_loop*1.0d0/fs_to_timeau, temp, &
                                             & energy_pot*Hartree_to_kcal,energy_kin*Hartree_to_kcal,   &
                                             & (energy_pot+energy_kin)*Hartree_to_kcal
                                             
else
  write(index,'(f20.4,3x,f20.10,3x,2f20.10,3x,2f20.10)'), timestep*index_loop*1.0d0/fs_to_timeau, temp,  &
                                             & energy_pot*Hartree_to_kcal,energy_kin*Hartree_to_kcal ,   &
                                             & (energy_pot + energy_kin)*Hartree_to_kcal
                                             
endif

endsubroutine output_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Kinetic Energy          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_kin_ener (velo,nbr_atomes,kin_ener,mass)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: velo

real*8,intent(inout):: kin_ener
real*8,intent(in)   :: mass
integer,intent(in)  :: nbr_atomes
integer             :: i,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

kin_ener=0.0d0

do i=1,nbr_atomes
 do k=1,3
   kin_ener=kin_ener+mass*velo(i,k)*velo(i,k)
 enddo
enddo

kin_ener=0.5d0*kin_ener

endsubroutine get_kin_ener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Instantaneous T         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_T (kin_ener,temp,nbr_atomes)

implicit none

real*8,intent(inout) :: temp
real*8,intent(in)    :: kin_ener
integer,intent(in)   :: nbr_atomes
real*8, parameter    :: kb=3.166811d-6   ! Hartree/K

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

temp=2.0d0*kin_ener/(kb*3.d0*dble(nbr_atomes))

endsubroutine get_T


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Apply Thermostat                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine thermostat(velo,thermo,temp,temp_0,val,nbr_atomes,timestep)

implicit none

real*8, allocatable, dimension (:,:),intent(inout)    :: velo

real*8               :: lambda     ! Scaling parameter of the thermostat
real*8, intent(in)   :: temp,temp_0,val,timestep
integer, intent(in)  :: thermo,nbr_atomes              
integer              :: i,j  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(thermo==1) then   ! Rescaling Thermostat

  lambda=sqrt(temp_0/temp)

  do i=1,nbr_atomes
    do j=1,3
      velo(i,j)=velo(i,j)*lambda
    enddo
  enddo

elseif (thermo==2) then ! Berendsen Thermostat


  lambda=1.0d0+(timestep/val)*((temp_0/temp)-1.0d0)
  lambda=sqrt(lambda)
  
  do i=1,nbr_atomes
    do j=1,3
      velo(i,j)=velo(i,j)*lambda
    enddo
  enddo

endif

endsubroutine thermostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Spherical Confinement Potential   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scp(nbr_atomes,positions,fscp,Rmax,scp_energy)
        implicit none
        integer,intent(in)                         :: nbr_atomes
        real*8, intent(in)                         :: Rmax
        real*8, allocatable,intent(inout)          :: positions(:,:)
        real*8, intent(out)                        :: fscp (:,:)
        real*8                      :: scp_energy
        real*8                      :: com(3)
        real*8                      :: distcom(nbr_atomes)
        integer                     :: i

        ! compute the center of mass 
        com = 0.0d0
        fscp = 0.0d0
        scp_energy = 0.0d0 

        do i = 1, nbr_atomes
          com(1:3) = com(1:3) + positions(i,1:3)
        end do
        
        com = com / nbr_atomes

        ! translate the com of the aggregate to the origin
        do i = 1,nbr_atomes
          positions(i,1:3) = positions(i,1:3) - com(1:3)
        end do

        ! calculate distcom
        do i = 1,nbr_atomes
          distcom(i) = SQRT((positions(i,1))**2 + (positions(i,2))**2 + (positions(i,3))**2)
        
           if (distcom(i) .gt. Rmax ) then
                  ! calculate force and energy 
                  fscp(i,1:3) = fscp(i,1:3) -4.0 * 0.008 * ((distcom(i)-Rmax)**(3)) * positions(i,1:3)/ distcom(i)
                  scp_energy = scp_energy + 0.008 * (distcom(i) - Rmax)**4
           end if
        end do

end subroutine scp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Assign initial velocities   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assign_initial_velocities(velocities,nbr_atomes,mass,temp_0)
        implicit none
        integer,intent(in)         :: nbr_atomes
        real*8,intent(inout)          :: temp_0
        real*8,intent(in)             :: mass
        real*8,allocatable, intent(inout)       :: velocities(:,:)
        integer                    :: i,j
        real*8,parameter                     :: kb=3.166811d-6
        
        ! calculate the velocities using the gauss_distr function
        do i = 1, nbr_atomes
          do j = 1, 3
            velocities(i,j) = velocities(i,j)+SQRT(kb*temp_0/mass)*gaussian_distr()
          end do
        end do
end subroutine assign_initial_velocities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function for Gauss Distribution   ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real *8 function gaussian_distr()
implicit none
integer, save :: iset=0
real *8 :: fac,rsq,var_1,var_2,x,y
real *8, save :: gset
call random_seed()
if (iset==0) then
! first time to calculate rsq
call random_number(x)
call random_number(y)
var_1 = 2.d0*x-1.d0
var_2 = 2.d0*y-1.d0
rsq = var_2**2+var_2**2
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
gaussian_distr = gset
iset = 0
endif
end function gaussian_distr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Periodic Boundary Conditions   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_pbc_distances(positions,distances,nbr_atomes,box_dim)
        implicit none
        real*8, intent(in)      :: positions(:,:)
        real*8                  :: distances(:,:)
        integer,intent(in)      :: nbr_atomes
        real*8, intent(in)      :: box_dim(3)
        integer                 :: i,j
        real*8                  :: dx, dy, dz

        do i = 1, nbr_atomes
          do j = 1, nbr_atomes
               ! get the distance between the atoms 
               dx = positions(i,1) - positions(j,1)
               dy = positions(i,2) - positions(j,2)
               dz = positions(i,3) - positions(j,3)
               
               ! applying the boundary conditions 
               if (dx > box_dim(1) / 2) dx = dx - box_dim(1)
               if (dx < -box_dim(1) / 2) dx = dx + box_dim(1)
               if (dy > box_dim(2) / 2) dy = dy - box_dim(2)
               if (dy < -box_dim(2) / 2) dy = dy + box_dim(2)
               if (dz > box_dim(3) / 2) dz = dz - box_dim(3)
               if (dz < -box_dim(3) / 2) dz = dz + box_dim(3)
               
               ! calculating the distances 
               distances(i,j) = dx**2 + dy**2 + dz**2
          end do
        end do

end subroutine get_pbc_distances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module utilities
