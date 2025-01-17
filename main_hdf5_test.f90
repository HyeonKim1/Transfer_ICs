program main
  use hdf5
  use grafic_types
  use grafic_io
  implicit none
  ! HDF5 관련 변수
  integer(HID_T) :: file_id, group_id, dset_id, space_id, ParticleType0 ,ParticleType1, Coordinates0, Coordinates1, Velocities0, Velocities1
  integer :: error  ! HDF5 error codes
  integer(HSIZE_T), dimension(1) :: dims
  integer(HSIZE_T), dimension(2) :: dims2
  

  ! 사용자 정의 변수와 상수
  logical :: file_exists
  integer :: nx, ny, nz
  real(8) :: InitTime, OmegaLambda, HubbleParam

  type(taille) :: headt
  type(cosmo)  :: headc

  real(sp), allocatable :: ic_array(:)

  integer :: NumPart, TotNumPart ,my_len
  real(8) :: Box, Omega, OmegaBaryon, Hubble
  character(len=20) :: FileBase, hdf5_name
  integer :: ThisTask, NTaskWithN

  integer :: ierr, i, ii, jj, kk, extra, start_index, end_index

  character(len=30) :: ic_name(12), message(4)
  real(sp), allocatable :: local_ic_array_x(:) ,local_ic_array_y(:),local_ic_array_z(:)


  integer, parameter :: BUFFER = 10
  integer :: Npart(6)
  integer(8) :: Nall(6)
  real(8) :: Massarr(6)
  real(8) :: Time, Redshift, BoxSize
  integer :: NumFiles


  real,dimension(:,:),allocatable :: POS
  real,dimension(:,:),allocatable :: VEL
  integer,allocatable :: ID0(:), ID1(:)
  real,allocatable :: U(:)


  integer :: blockmaxlen, maxidlen, maxlongidlen
  real(8), allocatable :: block(:), blockid(:)
  integer :: dummy
  character(len=100) :: buf
  integer :: fd=10
  integer :: type_class
  integer :: ndim, dims11(3), max_dims(3)



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
  hdf5_name='hdf5'

  call H5open_f(error)

  write(buf, '(A,".",I0,".",A)') trim(FileBase), ThisTask , trim(hdf5_name)

! Check if the file exists
  inquire(file=trim(buf), exist=file_exists)

  if (file_exists) then
      call system('rm ' // trim(buf))  ! Delete the existing file
  end if

  call h5fcreate_f(buf, H5F_ACC_TRUNC_F, file_id, error)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  my_len=nz*ny*nx
  allocate(ic_array(my_len))

  NumPart = TotNumPart / NTaskWithN
  extra = mod(TotNumPart, NTaskWithN)

  if (ThisTask < extra) then
    NumPart = NumPart + 1
  endif

  allocate(local_ic_array_x(NumPart))
  allocate(local_ic_array_y(NumPart))
  allocate(local_ic_array_z(NumPart))


  Npart = 0
  Nall = 0
  Massarr = 0.0

  Npart(2) = NumPart
  Nall(2) = TotNumPart

  Npart(1) = NumPart
  Nall(1) = TotNumPart
  Massarr(1) = (OmegaBaryon) * 3 * Hubble * Hubble / (8 * pi * 6.67430e-11) * Box**3 / TotNumPart
  Massarr(2) = (Omega - OmegaBaryon) * 3 * Hubble * Hubble / (8 * pi * 6.67430e-11) * Box**3 / TotNumPart

  Time = InitTime
  Redshift = 1.0 / InitTime - 1.0
  NumFiles = NTaskWithN
  BoxSize = Box
  
  dims(1)=6
  call h5gcreate_f(file_id, "Header",group_id, error)
  call h5screate_simple_f(1, dims, space_id, ierr)
  call h5dcreate_f(group_id, "NumPart_ThisFile", H5T_NATIVE_INTEGER, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, Npart, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(group_id, "NumPart_Total", H5T_NATIVE_INTEGER, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, Nall, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(group_id, "MassTable", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Massarr, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  dims(1)=1
  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(group_id, "Time", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Time, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(group_id, "Redshift", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Redshift, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(group_id, "BoxSize", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, BoxSize, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(group_id, "NumFilesPerSnapshot", H5T_NATIVE_INTEGER, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, NumFiles, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)


  ! 헤더 그룹 닫기
  call h5gclose_f(group_id, error)
  if (error /= 0) then
    print *, "Error Header."
    stop
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  if (ThisTask == 0) print *, "Header writing done"

  dims2(1)=3
  dims2(2)=NumPart

  call h5gcreate_f(file_id, "ParticleType0", ParticleType0 ,error)
  call h5gcreate_f(file_id, "ParticleType1", ParticleType1, error)

  call h5screate_simple_f(2, dims2, space_id, error)
  call h5dcreate_f(ParticleType0, "Coordinates", H5T_NATIVE_REAL, space_id, Coordinates0, error)
  call h5dcreate_f(ParticleType1, "Coordinates", H5T_NATIVE_REAL, space_id, Coordinates1, error)
  call h5dcreate_f(ParticleType0, "Velocities", H5T_NATIVE_REAL, space_id, Velocities0, error)
  call h5dcreate_f(ParticleType1, "Velocities", H5T_NATIVE_REAL, space_id, Velocities1, error)
  

  allocate(POS(3,NumPart))
  allocate(VEL(3,NumPart))


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
        POS(1,kk) = local_ic_array_x(kk)
        POS(2,kk) = local_ic_array_y(kk)
        POS(3,kk) = local_ic_array_z(kk)
      end do    
      
      print *, 1
      call h5dwrite_async(Coordinates0, H5T_NATIVE_REAL, POS, dims2, error)
      print *, 2
      stop
      call h5dclose_f(Coordinates0, error)
      call h5sclose_f(space_id, error)

      if (error /= 0) then
        print *, "Error ParticleType0 Coordinates."
        stop
      end if
    else if (ii == 2) then
      do kk = 1, NumPart
        POS(1,kk) = local_ic_array_x(kk)
        POS(2,kk) = local_ic_array_y(kk)
        POS(3,kk) = local_ic_array_z(kk)
      end do
      
      call h5dwrite_async(Coordinates1, H5T_NATIVE_REAL, POS, dims2, error)
      call h5dclose_f(Coordinates1, error)
      call h5sclose_f(space_id, error)
      if (error /= 0) then
        print *, "Error ParticleType1 Coordinates."
        stop
      end if
    else if (ii == 3) then
      do kk = 1, NumPart
        VEL(1,kk) = local_ic_array_x(kk)
        VEL(2,kk) = local_ic_array_y(kk)
        VEL(3,kk) = local_ic_array_z(kk)
      end do
      
      call h5dwrite_async(Velocities0, H5T_NATIVE_REAL, VEL, dims2, error)
      call h5dclose_f(Velocities0, error)
      call h5sclose_f(space_id, error)
      if (error /= 0) then
        print *, "Error ParticleType0 Velocities."
        stop
      end if
    else 
      do kk = 1, NumPart
        VEL(1,kk) = local_ic_array_x(kk)
        VEL(2,kk) = local_ic_array_y(kk)
        VEL(3,kk) = local_ic_array_z(kk)
      end do
      call h5dwrite_async(Velocities1, H5T_NATIVE_REAL, VEL, dims2, error)
      call h5dclose_f(Velocities1, error)
      call h5sclose_f(space_id, error)

      if (error /= 0) then
        print *, "Error ParticleType1 Velocities."
        stop
      end if
    endif
    if (ThisTask == 0) print *, message(ii)
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  deallocate(ic_array, local_ic_array_x, local_ic_array_y, local_ic_array_z)

  allocate(ID0(NumPart), ID1(NumPart))

  do i = 1, NumPart
    ID0(i) = i + start_index
    ID1(i) = i + start_index+ TotNumPart
  end do

  dims(1)=NumPart
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, ID0, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  if (error /= 0) then
    print *, "Error ParticleType0 ParticleIDs."
    stop
  end if

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(ParticleType1, "ParticleIDs", H5T_NATIVE_REAL, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, ID1, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)
  call h5gclose_f(ParticleType1, error)

  if (error /= 0) then
    print *, "Error ParticleType1 ParticleIDs."
    stop
  end if

  if (ThisTask == 0) print *, "ID writing done"
  
  ! 필요 시 제로 온도 쓰기
  allocate(U(NumPart))

  U = 0.0

  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(ParticleType0, "InternalEnergy", H5T_NATIVE_REAL, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, U, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)
  call h5gclose_f(ParticleType0, error)

  if (error /= 0) then
    print *, "Error InternalEnergy."
    stop
  end if
  

  call h5fclose_f(file_id, error)


  if (ThisTask == 0) print *, "all done"
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  deallocate(POS, VEL, ID0, ID1 ,U)
  
  call MPI_Finalize(ierr)

end program main
