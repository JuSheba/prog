program problem_Cauchy
  use param_func_mod     ! здесь ф-я и нач.параметры
  use rk_mod             ! метод Рунге-Кутты 4 порядка
  use methods_adams_mod  ! экстр. и интер. методы Адамса, внутри ещё модуль для коэф-ов
 !_____________________________________________________________________________________
  implicit none
  character(2) :: method
  real(8) :: t
  integer :: j
  real(8), dimension(1:size(data)) :: x, x_res
  real(8), dimension(:,:), allocatable :: x_ia, x_ea
  !_______________________________________________________
  ! ЧТЕНИЕ ИЗ ФАЙЛОВ
  !_______________________________________________________
  allocate(x_ea(-e_n+1:0, 1:size(data)),x_ia(-i_n+2:1,1:size(data)))

  call getarg(1,method)

  select case(method)
  case('rk')
  ! интегр. по промежутку до data_t + h
  open(11, file='rk.dat')
      write(11,*) 0d0, data
      t = 0d0    ! З. Коши
      x = data
      do while (t < data_t+h)
          x_res = rk(x, t)
          write(11,*) t, x_res
          t = t+h
          x = x_res
      end do
  close(11)

  case('ea')
    open(12,file='ea.dat')
        write(12,*) 0d0, data
        t = (e_n - 1d0) * h
        x_ea(-e_n+1,:) = data
        do j = -e_n+1, -1
          x_ea(j+1,:) = rk(x_ea(j,:), t+h*j)
          write(12,*) t, x_ea(j+1,:)
        end do

        x = x_ea(0,:)
        do while (t < data_t+h)
            x_res = e_adams(x, t, x_ea)
            do j = -e_n+1, -1
              x_ea(j,:) = x_ea(j+1,:)
            end do
            x_ea(0,:) = x_res
            write(12,*) t, x_res
            t = t + h
            x = x_res
        end do
    close(12)

  case('ia')
    open(13,file='ia.dat')
        write(13,*) 0d0, data
        t = (i_n - 2d0) * h
        x_ia(-i_n+2,:) = data
        do j = -i_n+2, -1
            x_ia(j+1,:) = rk(x_ia(j,:), t+h*j)
            write(13,*) t, x_ia(j+1,:)
        end do

        x = x_ia(0,:)
        do while (t < data_t+h)
            x_res = i_adams(x, t, x_ia)
            do j = -i_n+2, 0
              x_ia(j,:) = x_ia(j+1,:)
            end do
            x_ia(1,:) = x_res
            write(13,*) t, x_res
            t = t + h
            x = x_res
        end do

    close(13)
end select
end program problem_Cauchy
