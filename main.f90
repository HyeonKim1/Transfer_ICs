program main
  use grafic_types
  use grafic_io
  implicit none
  ! 사용자 정의 변수와 상수

  integer :: nx, ny, nz
  real(8) :: InitTime, OmegaLambda, HubbleParam

  type(taille) :: headt
  type(cosmo)  :: headc

  real(sp), allocatable :: ic_array(:)

  integer :: NumPart, TotNumPart ,my_len
  real(8) :: Box, Omega, OmegaBaryon, Hubble
  character(len=20) :: FileBase
  integer :: ThisTask, NTaskWithN

  integer :: ierr, i, ii, jj, kk, extra, start_index, end_index

  character(len=30) :: ic_name(12), message(4)
  real(sp), allocatable :: local_ic_array_x(:) ,local_ic_array_y(:),local_ic_array_z(:)


  integer, parameter :: BUFFER = 10
  type :: header
    integer :: Npart(6)
    integer(8) :: Nall(6)
    real(8) :: Massarr(6)
    real(8) :: Time, Redshift, BoxSize
    integer :: NumFiles
  end type header


  type format
    type(header) :: HEAD
    real,allocatable :: POS(:)
    real,allocatable :: VEL(:)
    integer,allocatable :: ID(:)
    real,allocatable :: U(:)
  end type format  

  type(format) :: my_format

  integer :: blockmaxlen, maxidlen, maxlongidlen
  real(8), allocatable :: block(:), blockid(:)
  integer :: dummy
  character(len=100) :: buf
  integer :: fd=10


    ! MPI 초기화
  call MPI_Init(ierr)

  ! 각 프로세스의 랭크를 얻음
  call MPI_Comm_rank(MPI_COMM_WORLD, ThisTask, ierr)

  ! 총 프로세스 수를 얻음
  call MPI_Comm_size(MPI_COMM_WORLD, NTaskWithN, ierr)


  ic_name(1) = '../ic_posbx'
  ic_name(2) = '../ic_posby'
  ic_name(3) = '../ic_posbz'

  ic_name(4) = '../ic_poscx'
  ic_name(5) = '../ic_poscy'
  ic_name(6) = '../ic_poscz'

  ic_name(7) = '../ic_velbx'
  ic_name(8) = '../ic_velby'
  ic_name(9) = '../ic_velbz'

  ic_name(10) = '../ic_velcx'
  ic_name(11) = '../ic_velcy'
  ic_name(12) = '../ic_velcz'

  message(1) = 'ic_pos_bayron_Transform done'
  message(2) = 'ic_pos_CDM_Transform done'
  message(3) = 'ic_vel_bayron_Transform done'
  message(4) = 'ic_vel_CDM_Transform done'


    ! 헤더 읽기 시도
  call grafic_read_header(trim(ic_name(1)), headt, headc)

  if (ThisTask == 0) then
  ! taille 타입 값들 출력
    print '(A, I6)', "nx =", headt%nx
    print '(A, I6)', "ny =", headt%ny
    print '(A, I6)', "nz =", headt%nz
    print '(A, F8.3)', "dx =", headt%dx
    print '(A, F8.3)', "lx =", headt%lx
    print '(A, F8.3)', "ly =", headt%ly
    print '(A, F8.3)', "lz =", headt%lz

    ! cosmo 타입 값들 출력
    print '(A, F8.5)', "astart =", headc%astart
    print '(A, F8.5)', "omegam =", headc%omegam
    print '(A, F8.5)', "omegav =", headc%omegav
    print '(A, F8.5)', "h0 =", headc%h0

  endif

  nx= headt%nx
  ny= headt%ny
  nz= headt%nz

  Omega=headc%omegam
  InitTime=headc%astart
  OmegaLambda=headc%omegav
  HubbleParam=headc%h0

  OmegaBaryon=0.044
  Hubble=70.
  Box=nx
  TotNumPart=nx*nx*nx
  FileBase='ics'

  write(buf, '(A,".",I0)') trim(FileBase), ThisTask
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  open(unit=fd, file=buf, status='replace', action='write',form='unformatted', iostat=ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
    print *, 'Error. Can''t write in file:', buf
    stop
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  my_len=nz*ny*2*(nx/2+1)
  allocate(ic_array(my_len))

  NumPart = TotNumPart / NTaskWithN
  extra = mod(TotNumPart, NTaskWithN)

  if (ThisTask < extra) then
    NumPart = NumPart + 1
  endif

  allocate(local_ic_array_x(NumPart))
  allocate(local_ic_array_y(NumPart))
  allocate(local_ic_array_z(NumPart))


  my_format%HEAD%Npart = 0
  my_format%HEAD%Nall = 0
  my_format%HEAD%Massarr = 0.0

  my_format%HEAD%Npart(2) = NumPart
  my_format%HEAD%Nall(2) = TotNumPart

  my_format%HEAD%Npart(1) = NumPart
  my_format%HEAD%Nall(1) = TotNumPart
  my_format%HEAD%Massarr(1) = (OmegaBaryon) * 3 * Hubble * Hubble / (8 * pi * 6.67430e-11) * Box**3 / TotNumPart
  my_format%HEAD%Massarr(2) = (Omega - OmegaBaryon) * 3 * Hubble * Hubble / (8 * pi * 6.67430e-11) * Box**3 / TotNumPart

  my_format%HEAD%Time = InitTime
  my_format%HEAD%Redshift = 1.0 / InitTime - 1.0
  my_format%HEAD%NumFiles = NTaskWithN
  my_format%HEAD%BoxSize = Box
  
  dummy = sizeof(my_format%HEAD)
  
  write(fd) dummy
  write(fd) my_format%HEAD
  write(fd) dummy
  

  allocate(my_format%POS(NumPart*3*2))
  allocate(my_format%VEL(NumPart*3*2))
  

  do ii = 1, 4
    do jj = 1, 3
      call grafic_read(ic_array, nz, 0, ny, nx, trim(ic_name(3*(ii-1)+jj))) 

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      if (jj == 1) then
        call MPI_Scatter(ic_array, NumPart, MPI_REAL, local_ic_array_x, NumPart, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      else if (jj == 2) then
        call MPI_Scatter(ic_array, NumPart, MPI_REAL, local_ic_array_y, NumPart, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      else 
        call MPI_Scatter(ic_array, NumPart, MPI_REAL, local_ic_array_z, NumPart, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      endif

    enddo

    start_index = ThisTask * (TotNumPart / NTaskWithN) + min(ThisTask, extra) + 1

    if (ii == 1) then
      do kk = 1, NumPart
        my_format%POS(3 * (kk - 1) + 1) = local_ic_array_x(kk)
        my_format%POS(3 * (kk - 1) + 2) = local_ic_array_y(kk)
        my_format%POS(3 * (kk - 1) + 3) = local_ic_array_z(kk)
      end do
    else if (ii == 2) then
      do kk = 1, NumPart
        my_format%POS(3 * (kk - 1) + 1 + NumPart) = local_ic_array_x(kk)
        my_format%POS(3 * (kk - 1) + 2 + NumPart) = local_ic_array_y(kk)
        my_format%POS(3 * (kk - 1) + 3 + NumPart) = local_ic_array_z(kk)
      end do
      dummy = sizeof(my_format%POS)
      
      write(fd) dummy
      write(fd) my_format%POS
      write(fd) dummy

    else if (ii == 3) then
      do kk = 1, NumPart
        my_format%VEL(3 * (kk - 1) + 1) = local_ic_array_x(kk)
        my_format%VEL(3 * (kk - 1) + 2) = local_ic_array_y(kk)
        my_format%VEL(3 * (kk - 1) + 3) = local_ic_array_z(kk)
      end do
    else 
      do kk = 1, NumPart
        my_format%VEL(3 * (kk - 1) + 1 + NumPart) = local_ic_array_x(kk)
        my_format%VEL(3 * (kk - 1) + 2 + NumPart) = local_ic_array_y(kk)
        my_format%VEL(3 * (kk - 1) + 3 + NumPart) = local_ic_array_z(kk)
      end do
      dummy = sizeof(my_format%VEL)
      
      write(fd) dummy
      write(fd) my_format%VEL
      write(fd) dummy
    endif
    if (ThisTask == 0) print *, message(ii)
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  deallocate(ic_array, local_ic_array_x, local_ic_array_y, local_ic_array_z)

  allocate(my_format%ID(NumPart*2))

  do i = 1, NumPart
    my_format%ID(i) = i + start_index
    my_format%ID(i+NumPart) = i + start_index+ TotNumPart
  end do
  dummy = sizeof(my_format%ID)
  
  write(fd) dummy
  write(fd) my_format%ID
  write(fd) dummy

  if (ThisTask == 0) print *, "ID writing done"
  
  ! 필요 시 제로 온도 쓰기
  allocate(my_format%U(NumPart))


  do i = 1, NumPart
    my_format%U(i) = 0.0
  end do
  dummy = sizeof(my_format%U)
  
  write(fd) dummy
  write(fd) my_format%U
  write(fd) dummy

  close(fd)


  if (ThisTask == 0) print *, "all done"
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  deallocate(my_format%POS, my_format%VEL, my_format%ID, my_format%U)
  
  call MPI_Finalize(ierr)

end program main
