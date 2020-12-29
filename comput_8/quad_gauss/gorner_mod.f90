module gorner_mod
  implicit none
  contains
  !_________________________________________________________________
    function gorner(x0, M)
    ! Получает на вход коэф-ы мн-на M и его корень x0, реализует
    ! метод Гoрнера деления мн-н на (x - x0)
      real(8)      :: x0
      integer(8)   :: i
      real(8)      :: M(0:), gorner(0:size(M)-1)

      gorner(0) = M(0)
      do i = 1, size(M)-1
        gorner(i) = gorner(i-1)*x0 + M(i)
      end do
      write(*,*) gorner
    end function gorner
  !_________________________________________________________________

end module gorner_mod
