program interpolationLagrange
implicit none

integer(4) :: j, i, k, numINTERV, numPOINTS, numPOINTS100

real(8)    :: pointL, pointR, find_polyINTERPOL, polyINTERPOL, pi, x, fi

real(8),dimension(:),allocatable   :: find_basisINTERPOL, basisINTERPOL

real(8),dimension(:,:),allocatable :: COORDS, result_COORDS

character(len=72) :: key

! Считываем данные из файла/файлов и делаем расчёт сетки.
call getarg(1, key)

select case(key)
  case('uniform')
    open(12, file='uniform.dat')
      read(12, '(2x, i10)')  numINTERV
      read(12,*)             pointL, pointR

      numPOINTS = numINTERV + 1

      allocate(COORDS(numPOINTS,2))
      write(*,*) 'Reading coordinate Y vector for uniform grid..'
          do i=1,numPOINTS
              read(12,*) COORDS(i,2)
          enddo
    close(12)

  case('chebyshev')
    open(13, file='chebyshev.dat')
      read(13, '(2x, i10)') numINTERV
      read(13,*)            pointL, pointR

      numPOINTS = numINTERV + 1

      allocate(COORDS(numPOINTS,2))

      write(*,*) 'Reading coordinate Y vector for chebyshev grid..'

      do i=1,numPOINTS
        read(13,*) COORDS(i,2)
      enddo
    close(13)
endselect

select case(key)
case('uniform')
  write(*,*) 'Calculation X coordinates of a uniform grid..'
  do i=1, numPOINTS
      COORDS(i,1) = 0.5*(pointR + pointL) +&
      0.5*(pointR - pointL)*(2.d0*(i-1.d0)/numINTERV - 1.d0)
  enddo

case('chebyshev')
  pi = 4.d0*atan(1.d0)
  write(*,*) 'Calculation X coordinates of a chebyshev grid..'
  do i=1,numPOINTS
        COORDS(i,1) = 0.5*(pointR + pointL) + 0.5*(pointR - pointL)*&
        cos(pi*(2*i + 1)/(2.d0*numINTERV + 2))
        write(*,*) COORDS(i, 1)
  enddo
endselect


!--------------------------------------------------------------------
!Вычисление 100N точек, в которых ищем зн-е многочлена и вывод в файл.

numPOINTS100 = 100*numPOINTS
allocate(result_COORDS(numPOINTS100,2))

write(*,*) 'Calculation 100*N intervals of a uniform grid..'
  do i=1, numPOINTS100
    result_COORDS(i,1) = 0.5*(pointR + pointL) +&
    0.5*(pointR - pointL)*(2.d0*(i-1.d0)/(numPOINTS100-1) - 1.d0)
  enddo

!do i = 1, numPOINTS100
!  write(*,*) result_COORDS(i, 1)
!enddo

do j = 1, numPOINTS100
  x = result_COORDS(j, 1)
  result_COORDS(j, 2) = 0
  do k = 1, numPOINTS
    fi = 1
    do i = 1, numPOINTS
      if (i/=k) then
          fi = fi*(x - COORDS(i, 1))/(COORDS(k, 1) - COORDS(i, 1))
      endif
    enddo
    result_COORDS(j, 2) = result_COORDS(j, 2) + COORDS(k, 2)*fi
  enddo
enddo

do i=1, numPOINTS100
  write(*, *) result_COORDS(i, 1), result_COORDS(i, 2)
enddo

end program interpolationLagrange
