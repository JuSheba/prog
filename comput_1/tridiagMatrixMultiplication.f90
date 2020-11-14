program tridiagmatrixMultiplication                                     ! Отличается от первого варианта лишь субрутиной.
implicit none

integer(4) :: n, n_in_b, i
real(8), dimension(:,:), allocatable :: A, B, C

open(11, file='data1.dat')
  read(11, '(2x, i10)') n

  allocate(A(n,3), B(n,3), C(n,5))   ! тут изменение, чтобы не хранить нули

  do i=1,n
    read(11,*) A(i,:)
  enddo
close(11)


open(12, file='data2.dat')
  read(12, '(2x, i10)') n_in_b
  do i=1,n
      read(12,*) B(i,:)
  enddo
close(12)

call multiplication_for_tridiag(A, B, C, n)

open(13, file='result.dat')
  write(13,*) C
close(13)

end program tridiagmatrixMultiplication

! Чтобы не считать произведение для всех элементов, посмотрим, какие
! обнулятся сразу и запишем в них 0. И выведем формулы для оставшихся
! элементов матрицы С. Оказалось, что она имеет 4-диаг. вид, где эл-ты
! на главной, 2 верхних и нижхних диагоналях высчитываются по опред-му
! правилу. Под них не подходял эл-ты С11 и Сnn, их посчитаем отдельно.


subroutine multiplication_for_tridiag(A, B, C, n)
implicit none
integer :: i, n
real(8), intent(in),  dimension(n,3) :: A, B
real(8), intent(out), dimension(n,5) :: C

  do i = 3, n-2
    ! Главная диагональ
    C(i,3) = A(i,1)*B(i-1,3) + A(i,2)*B(i,2) + A(i,3)*B(i+1,1)
    ! Вторая верхняя диагональ
    C(i,2) = A(i,1)*B(i-1,2) + A(i,2)*B(i,1)
    ! Вторая нижняя диагональ
    C(i,4) = A(i,2)*B(i,3) + A(i,3)*B(i+1,2)
    ! Третья верхняя диагональ
    C(i,1) = A(i,1)*B(i-1,1)
    ! Третья нижняя диагональ
    C(i,5) = A(i,3)*B(i+1,3)
  enddo

  C(1,1)   = A(1,1)*B(1,1) + A(1,2)*B(2,1)
  C(1,2)   = A(1,1)*B(1,2) + A(1,2)*B(2,2)
  C(1,3)   = A(1,2)*B(2,3)
  C(n,1)   = A(n,1)*B(n-1,1)
  C(n,2)   = A(n,1)*B(n-1,2) + A(n,2)*B(n,1)
  C(n,3)   = A(n,1)*B(n-1,3) + A(n,2)*B(n,2)

  C(2,1)   = A(2,1)*B(1,1)     + A(2,2)*B(2,1)
  C(2,2)   = A(2,1)*B(1,2)     + A(2,2)*B(2,2)     + A(2,3)*B(3,1)
  C(2,3)   = A(2,2)*B(2,3)     + A(2,3)*B(3,2)
  C(2,4)   = A(2,3)*B(3,3)
  C(n-1,1) = A(n-1,1)*B(n-2,1)
  C(n-1,2) = A(n-1,1)*B(n-2,2) + A(n-1,2)*B(n-1,1)
  C(n-1,3) = A(n-1,1)*B(n-2,3) + A(n-1,2)*B(n-1,2) + A(n-1,3)*B(n,1)
  C(n-1,4) = A(n-1,2)*B(n-1,3) + A(n-1,3)*B(n,2)

end subroutine multiplication_for_tridiag
