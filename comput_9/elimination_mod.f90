module elimination_mod
implicit none
contains
!_______________________________________________________________________
! 1.GAUSS ELIMINATION
!_______________________________________________________________________
subroutine Gauss(n,C,res)

  integer(8)  :: n, i, j
  real(8), dimension(:,:) :: C
  real(8), dimension(:)   :: res
  real(8) :: s, eps

!  write(*,*) 'You chose Gauss elimination.'
!  write(*,*) 'Calculation solution vector..'

  eps = 0.01

  do j = 1,n
    if(abs(C(j,j))<eps) then
      write(*,*) 'WARNING','j =',j,'This diagonal matrix element < eps'
    end if
    do i = j+1,n
     C(i,:) = C(i,:) - C(j,:)*C(i,j)/C(j,j)
    end do
  end do

  do i = n,1,-1
    s = C(i,n+1)
    do j = i+1,n
      s = s - C(i,j)*res(j)
    end do
    res(i) = s/C(i,i)
  end do

end subroutine Gauss
!_______________________________________________________________________
! 2.GAUSS ELIMINATION WITH PARTIAL PIVOTING
!_______________________________________________________________________
subroutine GaussPivot(n,eps,C,res)
  integer(8)  :: n, i, j, k, m
  real(8), dimension(:,:) :: C
  real(8), dimension(:)   :: res
  real(8), dimension(:), allocatable :: bar
  integer(4), dimension(:), allocatable :: t, coords
  real(8) :: eps, s

  !write(*,*) 'You chose Gauss elimination with partial pivoting.'
  !write(*,*) 'Calculation solution vector..'

  allocate(t(n), bar(n), coords(2))    ! сделаем массив из индексов и при каждой перестановке
  do i = 1, n       ! будем менять его, чтобы затем получить порядок переменных
    t(i) = i
  end do

  do j = 1, n

    if(abs(C(j,j)) < eps) then
      write(*,*) 'WARNING','j =',j,'This diagonal matrix element < eps'
    end if

    coords = maxloc(abs(C(j:n,j:n))) + j - 1
    bar = C(j,:)
    C(j,:) = C(coords(1),:)
    C(coords(1),:) = bar
    bar = C(:, j)
    C(:,j) = C(:, coords(2))
    C(:, coords(2)) = bar

    k = t(j)
    t(j) = t(coords(2))
    t(coords(2)) = k

    forall(m = j:n+1)             C(j,m) = C(j,m)/C(j,j)
    forall(i = j+1:n, m = j:n+1)  C(i,m) = C(i,m) - C(j,m) * C(i,j)

  end do

  do i = n,1,-1
    res(t(i)) = C(i,n+1) -  dot_product(C(i,i+1:), res(t(i+1:)))
  end do

end subroutine GaussPivot
!_______________________________________________________________________
! 3.GAUSS-JORDAN ELIMINATION
!_______________________________________________________________________
subroutine GaussJordan(n,eps,C,res)

  integer(8)  :: n, i, j
  real(8), dimension(:,:) :: C
  real(8), dimension(:)  :: res
  real(8) :: eps

  write(*,*) 'You chose Gauss-Jordan elimination.'
  write(*,*) 'Calculation solution vector..'

do j = 1,n
    C(j,:) = C(j,:) / C(j,j)
    do i = 1,n
      if(i/=j) then
        C(i,j+1:) = C(i,j+1:) - C(j,j+1:) * C(i,j)
      end if
    end do
  end do

  res = C(:,n+1)

!  do i=1,n
!    res(i) = C(i,n+1)
!  end do
end subroutine GaussJordan
!_______________________________________________________________________

end module elimination_mod
