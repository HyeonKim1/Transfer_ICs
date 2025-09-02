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
  integer :: nx, ny, nz , index, index2, n2x
  real(8) :: InitTime, OmegaLambda, HubbleParam
  integer :: file_number=1

  type(taille) :: headt
  type(cosmo)  :: headc


  integer :: NumPart, TotNumPart 
  real(8) :: Box, Omega, OmegaBaryon, dx
  character(len=20) :: FileBase, hdf5_name

  integer :: ierr, ii, jj, kk, extra, nn, i, j,k

  character(len=30) :: message(4)
  real(sp), allocatable :: ic_array_x(:) ,ic_array_y(:),ic_array_z(:)
  real(sp), allocatable :: re_ic_array_x(:) ,re_ic_array_y(:),re_ic_array_z(:)


  integer, parameter :: BUFFER = 10
  integer :: Npart(2)
  integer(8) :: Nall(2)
  real(8) :: Massarr(2)
  real(8) :: Time, Redshift
  integer :: NumFiles


  real,dimension(:,:),allocatable :: POS
  real,dimension(:,:),allocatable :: ini_POS
  real :: shfit_gas, shfit_CDM , grid_size
  real :: pos_x,pos_y,pos_z
  integer,allocatable :: ID0(:), ID1(:)
  real,allocatable :: U(:)

  character(len=100) :: buf
 
  real(8) :: G = 43.0187
  real(8) :: Hubble =100.


  call H5open_f(error)

    ! 헤더 읽기 시도
  call grafic_read_header(trim('../ic_postx'), headt, headc)


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
  dx=headt%dx

  Omega=headc%omegam
  InitTime=headc%astart
  OmegaLambda=headc%omegav
  HubbleParam=headc%h0

  OmegaBaryon=0.044
  TotNumPart=nz*ny*nx
  
  grid_size=dx*(HubbleParam/100)

  Box=grid_size*nx
  print *, "Boxsize =", Box
  print *, "grid_size =", grid_size
  FileBase='ics_tot'
  hdf5_name='hdf5'


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
  Massarr(1) = (OmegaBaryon)  * 3*Hubble*Hubble/(8*PI*G)* Box**3 / TotNumPart
  Massarr(2) = (Omega - OmegaBaryon) * 3*Hubble*Hubble/(8*PI*G)  * Box**3 / TotNumPart


  Time = InitTime
  Redshift = 1.0 / InitTime - 1.0
  NumFiles = file_number

  dims2(1)=3
  dims2(2)=NumPart


  allocate(ic_array_x(nz*ny*2*(nx/2+1)),ic_array_y(nz*ny*2*(nx/2+1)),ic_array_z(nz*ny*2*(nx/2+1)))
  allocate(re_ic_array_x(nz*ny*nz),re_ic_array_y(nz*ny*nz),re_ic_array_z(nz*ny*nz))
  allocate(POS(3,NumPart))
  allocate(ID0(NumPart), ID1(NumPart))
  allocate(U(NumPart))
  allocate(ini_POS(3,TotNumPart))

  n2x = 2*(nx/2+1)

  nn=1

  do kk=1, nz
    do jj=1, ny
      do ii=1, nx
        ini_POS(1,nn) = grid_size*(ii)
        ini_POS(2,nn) = grid_size*(jj)
        ini_POS(3,nn) = grid_size*(kk)
        nn = nn+1
      enddo 
    enddo
  enddo

  shfit_gas= -0.5*grid_size*(Omega-OmegaBaryon)/Omega
  shfit_CDM= 0.5*grid_size*OmegaBaryon/Omega

  print *, "shfit_gas =", shfit_gas
  print *, "shfit_CDM =", shfit_CDM


  write(buf, '(A,".",A)') trim(FileBase), trim(hdf5_name)
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
  call h5acreate_f(group_id, "NumPart_Total", H5T_STD_U64LE, space_id, dset_id, error)
  call h5awrite_f(dset_id, H5T_STD_U64LE, Nall, dims, error)
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
  call h5awrite_f(dset_id, H5T_NATIVE_DOUBLE, Box, dims, error)
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

  call h5gcreate_f(file_id, "PartType0", ParticleType0 ,error)
  call h5gcreate_f(file_id, "PartType1", ParticleType1, error)
  call h5screate_simple_f(2, dims2, space_id, error)

  call h5dcreate_f(ParticleType0, "Coordinates", H5T_NATIVE_REAL, space_id, Coordinates0, error)
  call h5dcreate_f(ParticleType1, "Coordinates", H5T_NATIVE_REAL, space_id, Coordinates1, error)
  call h5dcreate_f(ParticleType0, "Velocities", H5T_NATIVE_REAL, space_id, Velocities0, error)
  call h5dcreate_f(ParticleType1, "Velocities", H5T_NATIVE_REAL, space_id, Velocities1, error)

  print *, "Header writing done"

  
  call grafic_read(ic_array_x, nz, 0, ny, nx, trim('../ic_postx'))
  call grafic_read(ic_array_y, nz, 0, ny, nx, trim('../ic_posty')) 
  call grafic_read(ic_array_z, nz, 0, ny, nx, trim('../ic_postz'))

  ic_array_x=ic_array_x*(HubbleParam/100)
  ic_array_y=ic_array_y*(HubbleParam/100)
  ic_array_z=ic_array_z*(HubbleParam/100)


  index2 = 1
  
  do k=1, nz
    do j=1, ny
      do i=1, nx
        index=((k-1)*ny+j-1)*n2x+i
        pos_x=ini_POS(1,index2)+ic_array_x(index)+shfit_gas
        pos_y=ini_POS(2,index2)+ic_array_y(index)+shfit_gas
        pos_z=ini_POS(3,index2)+ic_array_z(index)+shfit_gas

        POS(1,index2)=MODULO(pos_x,Box)
        POS(2,index2)=MODULO(pos_y,Box)
        POS(3,index2)=MODULO(pos_z,Box)
        index2 = index2 + 1
      enddo
    enddo
  enddo


  call h5dwrite_f(Coordinates0, H5T_NATIVE_REAL, POS, dims2, error)
  call h5dclose_f(Coordinates0, error)
  call h5sclose_f(space_id, error)

  if (error /= 0) then
    print *, "Error ParticleType0 Coordinates."
    stop
  end if


  POS=POS+(grid_size/2)

  call h5screate_simple_f(2, dims2, space_id, error)
  call h5dwrite_f(Coordinates1, H5T_NATIVE_REAL, POS, dims2, error)
  call h5dclose_f(Coordinates1, error)
  call h5sclose_f(space_id, error)

  if (error /= 0) then
    print *, "Error ParticleType1 Coordinates."
    stop
  end if

  print *, 'ic_pos_Transform done'


  call grafic_read(ic_array_x, nz, 0, ny, nx, trim('../ic_veltx'))
  call grafic_read(ic_array_y, nz, 0, ny, nx, trim('../ic_velty')) 
  call grafic_read(ic_array_z, nz, 0, ny, nx, trim('../ic_veltz')) 

  index2 = 1

  do k=1, nz
    do j=1, ny
      do i=1, nx
        index=((k-1)*ny+j-1)*n2x+i
        POS(1,index2) =ic_array_x(index)/sqrt(InitTime)
        POS(2,index2) =ic_array_y(index)/sqrt(InitTime)
        POS(3,index2) =ic_array_z(index)/sqrt(InitTime)
        index2 = index2 + 1
      enddo
    enddo
  enddo

  
  call h5screate_simple_f(2, dims2, space_id, error)
  call h5dwrite_f(Velocities0, H5T_NATIVE_REAL, POS, dims2, error)
  call h5dclose_f(Velocities0, error)
  call h5sclose_f(space_id, error)

  if (error /= 0) then
    print *, "Error ParticleType0 Velocities."
    stop
  end if


  call h5screate_simple_f(2, dims2, space_id, error)
  call h5dwrite_f(Velocities1, H5T_NATIVE_REAL, POS, dims2, error)
  call h5dclose_f(Velocities1, error)
  call h5sclose_f(space_id, error)

  if (error /= 0) then
    print *, "Error ParticleType1 Velocities."
    stop
  end if

  print *, 'ic_vel_Transform done'


  do jj = 1, NumPart
    ID0(jj) = jj 
    ID1(jj) = jj + TotNumPart
  end do


  dims(1)=NumPart
  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(ParticleType0, "ParticleIDs", H5T_NATIVE_INTEGER, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ID0, dims, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(space_id, error)

  if (error /= 0) then
    print *, "Error ParticleType0 ParticleIDs."
    stop
  end if
  call h5screate_simple_f(1, dims, space_id, error)
  call h5dcreate_f(ParticleType1, "ParticleIDs", H5T_NATIVE_INTEGER, space_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ID1, dims, error)
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

  deallocate(POS, ic_array_x, ic_array_y, ic_array_z)
  deallocate(ID0, ID1)
  deallocate(U)
  deallocate(ini_POS)

  print *, "all done"
  
end program main
