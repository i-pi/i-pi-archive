      module atoms

      type Atom
         private
         integer :: atom_index
         real, dimension(3) :: pos
         real :: mass
      end type

      end module
