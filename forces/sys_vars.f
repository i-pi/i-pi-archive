      module sys_vars

      type Atom
         double precision, dimension(3) :: pos
      end type

      type Atom_vec
         type(Atom), dimension(:), pointer :: pos
      end type

      type Cell_vec
         double precision, dimension(3,3) :: h
         double precision, dimension(3,3) :: ih
      end type
      end module
