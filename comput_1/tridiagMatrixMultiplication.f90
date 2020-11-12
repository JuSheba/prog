program tridiagmatrixMultiplication                                     ! Отличается от первого варианта лишь субрутиной.
implicit none

integer(4) :: n, n_in_b, i
real(8), dimension(:,:), allocatable :: A, B, C

open(11, file='data1.dat')
read(11, '(2x, i10)') n

allocate(A(n,3), B(n,3), C(n,3))   ! тут изменение

do i=1,n
    read(11,*) A(i,1:n)
enddo

close(11)


open(12, file='data2.dat')
read(12, '(2x, i10)') n_in_b
do i=1,n
    read(12,*) B(i,1:n)
enddo
close(12)

call multiplication_for_tridiag(A, B, C, n)

open(13, file='result.dat')
write(*,*) C
close(13)

end program tridiagmatrixMultiplication

! Чтобы не считать произведение для всех элементов, посмотрим, какие
! обнулятся сразу и запишем в них 0. И выведем формулы для оставшихся
! элементов матрицы С. Оказалось, что она имеет 4-диаг. вид, где эл-ты
! на главной, 2 верхних и нижхних диагоналях высчитываются по опред-му
! правилу. Под них не подходял эл-ты С11 и Сnn, их посчитаем отдельно.


subroutine multiplication_for_tridiag(A, B, C, n)
implicit none
real(8), intent(in), dimension(n,n)  :: A, B
real(8), intent(out), dimension(n,n) :: C

! Главная диагональ
do i=2,n-1
     C(i,i) = A(i,i)*B(i,i) + A(i,i-1)*B(i-1,i) + A(i,i+1)*B(i+1,i)
enddo

! Вторая верхняя диагональ
do i=1,n-1
     C(i,i+1) = A(i,i)*B(i,i+1) + A(i,i+1)*B(i+1,i+1)
enddo

! Вторая нижняя диагональ
do i=2,n
      C(i,i-1) = A(i,i)*B(i,i-1) + A(i,i-1)*B(i-1,i-1)
enddo

! Третья верхняя диагональ
do i=1,n-2
      C(i,i+2) = A(i,i+1)*B(i+1,i+2)
enddo

! Третья нижняя диагональ
do i=3,n
      C(i,i-2) = A(i,i-1)*B(i-1,i+2)
enddo

! Первый и последний элементы на гл.диагонали.

      C(i,i-1) = A(i,i)*B(i,i-1) + A(i,i-1)*B(i-1,i-1)
enddo

end subroutine multiplication_for_tridiag
