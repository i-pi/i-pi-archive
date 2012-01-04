      module SG
         use sys_vars
      implicit none

      !private
      !public :: get_all
      double precision :: sigma, rc, rn, eps, correction

      double precision, parameter :: alpha = 1.713d0
      double precision, parameter :: beta = 1.5671d0
      double precision, parameter :: delta = 0.00993d0
      double precision, parameter :: delta_diff = delta*2.0d0
      double precision, parameter :: rc_exp = 8.32d0
      double precision, parameter :: C_6 = 12.14d0
      double precision, parameter :: C_8 = 215.2d0
      double precision, parameter :: C_9 = 143.1d0
      double precision, parameter :: C_10 = 4813.9d0
      double precision, parameter :: C_6_diff = -C_6*6d0
      double precision, parameter :: C_8_diff = -C_8*8d0
      double precision, parameter :: C_9_diff = -C_9*9d0
      double precision, parameter :: C_10_diff = -C_10*10d0
      
      contains

         subroutine f_c(r, on_r, long_range, long_range_diff)
            double precision, intent(in) :: on_r
            double precision, intent(in) :: r
            double precision, intent(out) :: long_range
            double precision, intent(out) :: long_range_diff

            double precision dist_frac

            if (r > rc_exp) then
               long_range = 1.0d0
               long_range_diff = 0.0d0
            else
               dist_frac = rc_exp*on_r-1.0d0
               long_range = dexp(-(dist_frac)**2)
               long_range_diff = 2.0d0*(dist_frac)*rc_exp*on_r**2*long_range
            end if

         end subroutine

         subroutine exp_func(r, pot, force)
            double precision, intent(in) :: r
            double precision, intent(out) :: pot
            double precision, intent(out) :: force

            double precision power

            pot = dexp(alpha-r*(beta+delta*r))
            force = (beta+delta_diff*r)*pot

         end subroutine

         subroutine SG_functions(r, pot, force)
            double precision, intent(in) :: r
            double precision, intent(out) :: pot
            double precision, intent(out) :: force

            double precision on_r, exp_pot, exp_force, long_range, 
     1long_range_diff, disp, disp_diff

            if (r > rc) then
               pot = 0.0d0
               force = 0.0d0
               return
            end if

            on_r = 1.0d0/r
            call exp_func(r, exp_pot, exp_force)
            call f_c(r, on_r, long_range, long_range_diff)
      
            disp = (C_6*on_r**6+C_8*on_r**8-C_9*on_r**9+C_10*on_r**10)
            disp_diff = (C_6_diff*on_r**7 + C_8_diff*on_r**9 - 
     1C_9_diff*on_r**10 + C_10_diff*on_r**11)

            pot = exp_pot - disp*long_range - correction
            force = exp_force+disp_diff*long_range+disp*long_range_diff
         
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

            s = matmul(cell%ih, atoms(i)%pos - atoms(j)%pos)
            do k = 1, 3
               s(k) = s(k) - dnint(s(k))
            end do
            rij = matmul(cell%h, s)
            r = sqrt(dot_product(rij,rij))

         end subroutine

         subroutine SG_fij(atoms, i, j, cell, fij, rij, pot)
            type(Atom), dimension(:), intent(in) :: atoms
            integer, intent(in) :: i
            integer, intent(in) :: j
            type(Cell_vec), intent(in) :: cell
            double precision, dimension(3), intent(out) :: fij
            double precision, dimension(3), intent(out) :: rij
            double precision, intent(out) :: pot

            double precision f_tot, r
            integer k

            call separation(atoms, i, j, cell, rij, r)
            call SG_functions(r, pot, f_tot)

            do k = 1, 3
               fij(k) = f_tot*rij(k)/r
            end do

         end subroutine

         subroutine nearest_neighbours(atoms, cell, n_list, index_list)
            type(Atom), dimension(:), intent(in) :: atoms
            type(Cell_vec), intent(in) :: cell
            integer, dimension(size(atoms)-1), intent(out) :: n_list
            integer, dimension(size(atoms)*(size(atoms)-1)/2),
     1                         intent(out) :: index_list

            double precision, dimension(3) :: rij
            double precision r
            integer i, j, counter

            n_list = 0
            index_list = 0
            counter = 1
            
            do i = 1, size(atoms)-1
               do j = i+1, size(atoms)
                  call separation(atoms, i, j, cell, rij, r)
                  if (r < rn) then
                     n_list(i) = n_list(i) + 1
                     index_list(counter) = j
                     counter = counter + 1
                  end if
               end do
            end do

         end subroutine

         subroutine get_all(atoms, cell, n_list, index_list, 
     1                      pot, f, vir)
            type(Atom), dimension(:), intent(in) :: atoms
            type(Cell_vec), intent(in) :: cell
            integer, dimension(size(atoms)-1), intent(in) :: n_list
            integer, dimension(size(atoms)*(size(atoms)-1)/2),
     1                         intent(in) :: index_list
            double precision, intent(out) :: pot
            double precision, intent(out), dimension(3,size(atoms))
     1 :: f
            double precision, dimension(3,3), intent(out) :: vir

            double precision volume
            integer i, j, k, l
            integer start, finish
            double precision, dimension(3) :: fij, rij
            double precision pot_ij

            vir = 0.0d0
            pot = 0.0d0
            f = 0.0d0
            volume = cell%h(1,1)*cell%h(2,2)*cell%h(3,3)

            start = 1
            do i = 1, size(atoms)-1
               finish = start + n_list(i) - 1
               !do j = i+1, size(atoms)
               do j = start, finish
                  call SG_fij(atoms, i, index_list(j), cell, fij, 
     1                        rij, pot_ij)
                  f(:,i) = f(:,i) + fij
                  f(:,index_list(j)) = f(:,index_list(j)) - fij
                  pot = pot + pot_ij

                  do k = 1, 3
                     do l = k, 3
                        vir(k,l) = vir(k,l) + fij(k)*rij(l)
                     end do
                  end do

               end do
               start = finish + 1

            end do
            !vir = vir/volume

         end subroutine

      end module
