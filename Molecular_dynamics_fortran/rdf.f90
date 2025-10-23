program calculate_distances
      implicit none
      integer                  :: i, j, n_atoms
      real*8,allocatable       :: x(:), y(:),z(:)
      real*8                   :: dx, dy, dz, dist
      character(len=10)        :: atom

      open(unit=10,file='final.xyz')
         read(10,*) n_atoms
         read(10,*)

         allocate(x(n_atoms),y(n_atoms),z(n_atoms))

         do i = 1, n_atoms
            read(10,*) atom, x(i), y(i), z(i)
         end do

      close(10)

      open(unit=20, file='distance.dat')
        do i = 1, n_atoms-1
           do j = i+1, n_atoms
              dx = x(i) - x(j)
              dy = y(i) - y(j)
              dz = z(i) - z(j)

              dist = sqrt(dx*dx + dy*dy + dz*dz)

              write(20, *) dist
            end do
        end do

      close(20)

      deallocate(x, y, z)

end program calculate_distances
