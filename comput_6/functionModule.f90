module functionModule
  implicit none
  contains

    function func(X)
      integer :: i, n
      real(8), dimension(:) :: X
      real(8), dimension(size(X)) :: func
        n = size(X)
        forall(i = 1:n) func(i) = sin(X(i))
    end function func

end module functionModule
