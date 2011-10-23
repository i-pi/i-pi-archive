      module LJ
         use sys_vars
      implicit none

      double precision, parameter :: sigma = 0.15d0
      double precision, parameter :: rc = 2.5d0*sigma
      double precision, parameter :: eps = 0.1d0
      double precision, parameter :: correction = 
     1(sigma/rc)**12 - (sigma/rc)**6

      contains

         subroutine LJ_functions(r, pot, force)
            double precision, intent(in) :: r
            double precision, intent(out) :: pot
            double precision, intent(out) :: force

            double precision sigma_by_r

            if (r > rc) then
               pot = 0.0d0
               force = 0.0d0
               return
            end if

            sigma_by_r = sigma/r

            pot = 4*eps*((sigma_by_r)**12 - (sigma_by_r)**6 -correction)

            force = 4*eps*(12/r*(sigma_by_r)**12 - 6/r*(sigma_by_r)**6)

         end subroutine

         subroutine separation(atoms, i, j, cell, rij, r)
            type(Atom), dimension(:), intent(in) :: atoms
            integer, intent(in) :: i
            integer, intent(in) :: j
            type(Cell_vec), intent(in) :: cell
            double precision, dimension(3), intent(out) :: rij
            double precision, intent(out) :: r

            integer k
            double precision, dimension(3) :: s

            s = mat_mul(cell%ih, atoms(i)%pos - atoms(j)%pos)
            do k = 1, 3
               s(k) = s(k) - dnint(s(k))
            end do
            rij = mat_mul(cell%h, s)
            r = sqrt(dot_product(rij, rij))

         end subroutine

      end module
