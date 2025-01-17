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
  integer :: file_number=16

  type(taille) :: headt
  type(cosmo)  :: headc


  integer :: NumPart, TotNumPart ,my_len
  real(8) :: Box, Omega, OmegaBaryon, Hubble
  character(len=20) :: FileBase, hdf5_name

  integer :: ierr, i, ii, jj, kk, extra, start_index, end_index
  integer :: n_dim, nz_end, nz_start

  character(len=30) :: ic_name(12), message(4)
  real(sp), allocatable :: ic_array_x(:) ,ic_array_y(:),ic_array_z(:)


  integer, parameter :: BUFFER = 10
  integer :: Npart(2)
  integer(8) :: Nall(2)
  real(8) :: Massarr(2)
  real(8) :: Time, Redshift, BoxSize
  integer :: NumFiles


  real,dimension(:,:),allocatable :: POS
  integer,allocatable :: ID0(:), ID1(:)
  real,allocatable :: U(:)


  integer :: blockmaxlen, maxidlen, maxlongidlen
  real(8), allocatable :: block(:), blockid(:)
  integer :: dummy
  character(len=100) :: buf
  integer :: fd=10
  integer :: type_class
  integer :: ndim, dims11(3), max_dims(3)
  integer(4) :: hdf5_uint32 , hdf5_uint64


  call H5open_f(error)


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



    ! 헤더 읽기 시도
  call grafic_read_header(trim(ic_name(1)), headt, headc)


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
  TotNumPart=nz*ny*nx
  FileBase='ics'
  hdf5_name='hdf5'
  my_len=nz*ny*nx


  if (MOD(TotNumPart, file_number) /= 0) then 
    print *, "The number is not divisible. The total number of particles is ", TotNumPart, "."
    stop
  end if

  NumPart = TotNumPart / file_number

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
  NumFiles = file_number
  BoxSize = Box

  dims2(1)=3
  dims2(2)=NumPart

  start_index=1
  end_index=NumPart
  
  nz_end=nz/file_number
  nz_start=1
  n_dim=nz/file_number

  allocate(ic_array_x(nz*ny*2*(nx/2+1)),ic_array_y(nz*ny*2*(nx/2+1)),ic_array_z(nz*ny*2*(nx/2+1)))
  allocate(POS(3,NumPart))
  allocate(ID0(NumPart), ID1(NumPart))
  allocate(U(NumPart))


  do ii=1, file_number

    write(buf, '(A,".",I0,".",A)') trim(FileBase), ii-1 , trim(hdf5_name)
    print *, buf

    ! Check if the file exists
    inquire(file=trim(buf), exist=file_exists)

    if (file_exists) then
        call system('rm ' // trim(buf))  ! Delete the existing file
    end if

    dims(1)=2

    call h5fcreate_f(buf, H5F_ACC_TRUNC_F, file_id, error)
    call h5gcreate_f(file_id, "Header",group_id, error)

    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(group_id, "NumPart_ThisFile", H5T_NATIVE_INTEGER, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_INTEGER, Npart, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    call h5screate_simple_f(1, dims, space_id, error)
    call h5acreate_f(group_id, "NumPart_Total", H5T_NATIVE_INTEGER, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_INTEGER, Nall, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    call h5screate_simple_f(1, dims, space_id, error)
    call h5acreate_f(group_id, "MassTable", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_DOUBLE, Massarr, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    dims(1)=1
    call h5screate_simple_f(1, dims, space_id, error)
    call h5acreate_f(group_id, "Time", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_DOUBLE, Time, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    call h5screate_simple_f(1, dims, space_id, error)
    call h5acreate_f(group_id, "Redshift", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_DOUBLE, Redshift, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    call h5screate_simple_f(1, dims, space_id, error)
    call h5acreate_f(group_id, "BoxSize", H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_DOUBLE, BoxSize, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    call h5screate_simple_f(1, dims, space_id, error)
    call h5acreate_f(group_id, "NumFilesPerSnapshot", H5T_NATIVE_INTEGER, space_id, dset_id, error)
    call h5awrite_f(dset_id, H5T_NATIVE_INTEGER, NumFiles, dims, error)
    call h5aclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    ! 헤더 그룹 닫기
    call h5gclose_f(group_id, error)
    if (error /= 0) then
      print *, "Error Header."
      stop
    end if

    call h5gcreate_f(file_id, "ParticleType0", ParticleType0 ,error)
    call h5gcreate_f(file_id, "ParticleType1", ParticleType1, error)
    call h5screate_simple_f(2, dims2, space_id, error)

    call h5dcreate_f(ParticleType0, "Coordinates", H5T_NATIVE_REAL, space_id, Coordinates0, error)
    call h5dcreate_f(ParticleType1, "Coordinates", H5T_NATIVE_REAL, space_id, Coordinates1, error)
    call h5dcreate_f(ParticleType0, "Velocities", H5T_NATIVE_REAL, space_id, Velocities0, error)
    call h5dcreate_f(ParticleType1, "Velocities", H5T_NATIVE_REAL, space_id, Velocities1, error)

    print *, "Header writing done"
    
    call grafic_read(ic_array_x, nz, 0, ny, nx, trim(ic_name(1))) 
    call grafic_read(ic_array_y, nz, 0, ny, nx, trim(ic_name(2))) 
    call grafic_read(ic_array_z, nz, 0, ny, nx, trim(ic_name(3))) 
  

    do kk = 1, NumPart
      POS(1,kk) = ic_array_x(kk)
      POS(2,kk) = ic_array_y(kk)
      POS(3,kk) = ic_array_z(kk)
    end do

    call h5dwrite_f(Coordinates0, H5T_NATIVE_REAL, POS, dims2, error)
    call h5dclose_f(Coordinates0, error)
    call h5sclose_f(space_id, error)

    if (error /= 0) then
      print *, "Error ParticleType0 Coordinates."
      stop
    end if

    print *, 'ic_pos_bayron_Transform done'

    call grafic_read(ic_array_x, nz, 0, ny, nx, trim(ic_name(4))) 
    call grafic_read(ic_array_y, nz, 0, ny, nx, trim(ic_name(5))) 
    call grafic_read(ic_array_z, nz, 0, ny, nx, trim(ic_name(6))) 

    do kk = 1, NumPart
      POS(1,kk) = ic_array_x(kk)
      POS(2,kk) = ic_array_y(kk)
      POS(3,kk) = ic_array_z(kk)
    end do
    call h5screate_simple_f(2, dims2, space_id, error)
    call h5dwrite_f(Coordinates1, H5T_NATIVE_REAL, POS, dims2, error)
    call h5dclose_f(Coordinates1, error)
    call h5sclose_f(space_id, error)

    if (error /= 0) then
      print *, "Error ParticleType1 Coordinates."
      stop
    end if

    print *, 'ic_pos_CDM_Transform done'


    call grafic_read(ic_array_x, nz, 0, ny, nx, trim(ic_name(7))) 
    call grafic_read(ic_array_y, nz, 0, ny, nx, trim(ic_name(8))) 
    call grafic_read(ic_array_z, nz, 0, ny, nx, trim(ic_name(9))) 
    
    do kk = 1, NumPart
      POS(1,kk) = ic_array_x(kk)
      POS(2,kk) = ic_array_y(kk)
      POS(3,kk) = ic_array_z(kk)
    end do
    
    call h5screate_simple_f(2, dims2, space_id, error)
    call h5dwrite_f(Velocities0, H5T_NATIVE_REAL, POS, dims2, error)
    call h5dclose_f(Velocities0, error)
    call h5sclose_f(space_id, error)

    if (error /= 0) then
      print *, "Error ParticleType0 Velocities."
      stop
    end if

    print *, 'ic_vel_bayron_Transform done'


    call grafic_read(ic_array_x, nz, 0, ny, nx, trim(ic_name(10))) 
    call grafic_read(ic_array_y, nz, 0, ny, nx, trim(ic_name(11))) 
    call grafic_read(ic_array_z, nz, 0, ny, nx, trim(ic_name(12))) 
   
    do kk = 1, NumPart
      POS(1,kk) = ic_array_x(kk)
      POS(2,kk) = ic_array_y(kk)
      POS(3,kk) = ic_array_z(kk)
    end do

    call h5screate_simple_f(2, dims2, space_id, error)
    call h5dwrite_f(Velocities1, H5T_NATIVE_REAL, POS, dims2, error)
    call h5dclose_f(Velocities1, error)
    call h5sclose_f(space_id, error)

    if (error /= 0) then
      print *, "Error ParticleType1 Velocities."
      stop
    end if

    print *, 'ic_vel_CDM_Transform done'
  

    do i = 1, NumPart
      ID0(i) = i + start_index
      ID1(i) = i + start_index+ TotNumPart
    end do

    dims(1)=NumPart
    call h5screate_simple_f(1, dims, space_id, error)
    call h5dcreate_f(ParticleType0, "ParticleIDs", H5T_NATIVE_REAL, space_id, dset_id, error)
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

    print *, "ID writing done"
  

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

    print *, ii, "file all done"
    start_index=start_index+NumPart
    end_index=end_index+NumPart

  end do 

  deallocate(POS, ic_array_x, ic_array_y, ic_array_z)
  deallocate(ID0, ID1)
  deallocate(U)

  print *, "all done"
  
end program main
