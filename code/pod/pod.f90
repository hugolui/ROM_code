!################################################################################################
program spectral_POD
  
  use mod_field, only : logical_primitive,corrFlag, soln_files
  use mod_pod_modes
  implicit none
  
  call read_inputs()
  
  call data_setup

  if (soln_files .eq. 1) then

    call read_grid
    
    call read_soln
  
  end if

  if (soln_files .eq. 2) then

    call read_grid_multi
    
    call read_soln_multi
  
  end if

  call Convert_to_Primitive_Variables(logical_primitive)
  
  call Reynolds_Decomposition
  
  if (corrFlag) then
    
    call compute_cell_area

    call correlation_matrix
    
  endif

  call spatial_modes_POD

  call move_files() 

end program
!#################################################################################  

subroutine read_inputs()
    
    
    use mod_pod_modes
    use mod_field
    use mod_CGNS
    implicit none
    integer, parameter :: fileUnit = 55
    
    open(fileUnit, file='inputs.inp')
    read(fileUnit,'(A)') path_to_soln
    read(fileUnit,*)     soln_files
    read(fileUnit,'(A)') path_to_ROM
    read(fileUnit,'(A)') soln
    read(fileUnit,'(A)') outfile
    read(fileUnit,'(A)') outfile2
    read(fileUnit,*)     corrFlag
    read(fileUnit,*)     svdFlag
    read(fileUnit,*)     idxi, idxf, idxr
    read(fileUnit,*)     idxi_val, idxf_val, idxr_val
    read(fileUnit,*)     logical_primitive
    read(fileUnit,*)     jmaxInput
    read(fileUnit,*)     imin, imax
    read(fileUnit,*)     der

    call remove_comments(outfile, stringLen)
    call remove_comments(outfile2, stringLen)
    call remove_comments(path_to_soln, stringLen)
    call remove_comments(path_to_ROM, stringLen)
    call remove_comments(soln, stringLen)
    
    close(fileUnit)

    if (soln_files .eq. 1) then
      write(CGNS_filename,'(A)') trim(path_to_soln) // trim(soln) // ".cgns"
    end if

    if (soln_files .eq. 2) then
      write(CGNS_filename,'(A,I6.6,A)') trim(path_to_soln) // '/' // trim(soln), idxi, ".cgns"
    end if

end subroutine read_inputs

!#################################################################################  
subroutine remove_comments(string, n)

  implicit none
  integer, intent(in)             :: n
  character(len=n), intent(inout) :: string

  integer                         :: i
  logical                         :: stringCut
  
  stringCut = .false.
  i = 1
  do while(i<n )

    if (string(i:i)== achar(32)) then
      stringCut = .true.
    endif
    
    if (stringCut) string(i:i) = ' '
    i = i + 1
  enddo

end subroutine remove_comments
!#################################################################################  

subroutine data_setup

  use mod_field
  use mod_variable_precision
  use mod_pod_modes  
  implicit none

  call system('mkdir -p ./SVD_files/')
  call system('rm -rf ./SVD_files/V_matrix.dat')  
  call system('rm -rf ./SVD_files/sigma.dat')

  if (VRL .eq. 4) then
    call system('printf "\033[0;31m \n"')
    call system('printf "   WARNING: Working with SINGLE PRECISION... \n"')
    call system('printf "   WARNING: Working with SINGLE PRECISION... \n"')
    call system('printf "   WARNING: Working with SINGLE PRECISION... \n"')
    call system('printf "\033[0m \n"')
  endif  

  write(*,'(A)') '  Grid file: ' // trim(path_to_grid) // achar(10)
  
  write(*,'(A,I0,A,I0,A,I0,A)') '  First file: ', idxi, achar(10) //&
                                '  Last file : ', idxf, achar(10) //&
                                '  Skip      : ', idxr, achar(10)
  
  return
  
end subroutine data_setup
!#################################################################################  

!read grid data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_grid

  use mod_field
  use mod_CGNS
  use cgns
  implicit none
  integer(kind=4) :: m, nxTotal, nyTotal, ijk_min(3), ijk_max(3), nx, ny1, ny2
  
  index_base = 1
  ! Find out if the second zone will be used in calculations
  call cg_open_f(trim(CGNS_filename),CG_MODE_READ,index_grid,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  ! Get dimension of isize
  call cg_index_dim_f(index_grid,1,1,isizeDim, ier)
  allocate(isize(isizeDim,3))
  isize(:,:) = 1

  call cg_zone_read_f(index_grid,1,1,zonename,isize,ier)
  ny1 = isize(2,1)
  call cg_zone_read_f(index_grid,1,2,zonename,isize,ier)
  ny2 = isize(2,1)
  nyTotal = ny1 + ny2
  if ( jmaxInput <= ny1 ) then 
    nZones = 1;
  else
    nZones = 2;
  endif

  ! Set up variables
  allocate(zone(1:nzones))
  allocate(jmax(nzones))
  allocate(kmax(nzones))
  jMax = 0
  nx = imax-imin+1

  write(*,*) ' Reading grid data ...' 
  
  ! Allocate grid variables x and y 
  do m = 1,nzones
    kmax(m) = 1
  
    ! Get dimensions of zone
    call cg_zone_read_f(index_grid,index_base,m,zonename,isize,ier)
    nxTotal = isize(1,1)
    nyTotal = isize(2,1)

    ! Get the appropriate value of jmax
    if (m<=nzones/2) then
      jmax(m) = nyTotal
    else
      jmax(m) = jmaxInput - sum(jmax(:))
    endif
    
    if (nxTotal<nx) stop ' I-range inserted by user not valid'

    zone(m)%nx1 = nx
    zone(m)%nx2 = jmax(m)
    zone(m)%nx3 = 1

    allocate(zone(m)%x(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3))
    allocate(zone(m)%y(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3))

    write(*,'(A,I0,A,I0,A,I0)') '    Zone ',m, '                        - nxTotal = ', nxTotal, ' nyTotal = ', nyTotal
    write(*,'(A,I0,A,I0,A,I0)') '    Domain used in calculations  ' //' : nx      = ', zone(m)%nx1,      ' ny      = ', zone(m)%nx2

  enddo

  ! Read grid
  ijk_min(1) = imin
  ijk_min(2) = 1
  ijk_min(3) = 1
  ijk_max(1) = imax
  ijk_max(3) = 1
  do m = 1, nzones
    ijk_max(2) = jmax(m)
    call cg_coord_read_f(index_grid,index_base,m,'CoordinateX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%x(:,:,:),ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    call cg_coord_read_f(index_grid,index_base,m,'CoordinateY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%y(:,:,:),ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  enddo


  call cg_close_f(index_grid,ier)
  deallocate(isize)
  return
end subroutine
!#################################################################################  

!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_soln

  use mod_field
  use mod_CGNS
  use cgns
  
  implicit none
  integer(kind=4)           :: m, t
  integer(kind=4)           :: iFile
  integer                   :: dataType, dataDimensions, dataLength(12), pos
  real(kind=8), allocatable :: timeGlobal(:)
  character(len=50)         :: tempName
  real(kind=8)              :: lastPrint, percentage
  integer                   :: ijk_min(3), ijk_max(3), iFlow

  write(*,*) achar(10) // '  Reading solution data ...' 

  nsnap = (idxf - idxi)/idxr + 1

  do m = 1,nzones
      
    allocate(zone(m)%r(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%u(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%v(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%w(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%p(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))                

    zone(m)%r(:,:,:,:) = 0.0d0
    zone(m)%u(:,:,:,:) = 0.0d0
    zone(m)%v(:,:,:,:) = 0.0d0
    zone(m)%w(:,:,:,:) = 0.0d0
    zone(m)%p(:,:,:,:) = 0.0d0 
  
  enddo  
  allocate(time(nsnap))
  
  !---- reading the solution files
  write(CGNS_solnname,'(A)') trim(path_to_soln) // trim(soln) // '.cgns'
  call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  index_base = 1

  ! Get index of first FlowSolution inside zones
  call cg_sol_info_f(index_soln,index_base,1,1,solnname,pos,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  read(solnname(14:stringLen),*) iFirstFlowSol
  iFirstFlowSol = iFirstFlowSol - 1

  ! Read all time instants
  call cg_gopath_f(index_soln, '/Base/TimeIterValues', ier)
  if (ier .ne. CG_OK) call cg_error_exit_f 
  call cg_array_info_f(1,tempName,dataType,dataDimensions,dataLength,ier)
  allocate(timeGlobal(dataLength(1)))
  call cg_array_read_f(1,timeGlobal,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  t = 0
  lastPrint = 0d0

  ! Read flow solution on the prompted locations
  ijk_min(1) = imin
  ijk_min(2) = 1
  ijk_min(3) = 1
  ijk_max(1) = imax
  ijk_max(3) = 1
  do iFile = idxi,idxf,idxr

    t = t + 1
    time(t) = timeGlobal(iFile-iFirstFlowSol)

    iFlow = iFile - iFirstFlowSol

    ! Read one zone at a time
    do m = 1, nZones
      ijk_max(2) = jmax(m)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'Density'  ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%r(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%u(:,:,:,t),ier)    
      call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%v(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumZ',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%w(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'Pressure' ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%p(:,:,:,t),ier)   
    enddo

    ! Counter 
    percentage = dble(t) / dble(nsnap) * 100d0
    if (percentage>lastPrint) then
      write(*,'(A,F0.0,A)') '  ',percentage, ' %'
      lastPrint = lastPrint + 10d0
    endif
      
  enddo

  percentage = 100d0
  write(*,'(A,F0.0,A)') '  ',percentage, ' %'

  deallocate(timeGlobal)
  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  
end subroutine

!#################################################################################

!read grid data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_grid_multi

  use mod_field
  use mod_CGNS
  use cgns
  implicit none
  integer(kind=4) :: m, nxTotal, nyTotal, ijk_min(3), ijk_max(3), nx, ny1, ny2
  
  index_base = 1
  ! Find out if the second zone will be used in calculations
  call cg_open_f(trim(CGNS_filename),CG_MODE_READ,index_grid,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  ! Get dimension of isize
  call cg_index_dim_f(index_grid,1,1,isizeDim, ier)
  allocate(isize(isizeDim,3))
  isize(:,:) = 1

  call cg_zone_read_f(index_grid,1,1,zonename,isize,ier)
  ny1 = isize(2,1)
  call cg_zone_read_f(index_grid,1,2,zonename,isize,ier)
  ny2 = isize(2,1)
  nyTotal = ny1 + ny2
  if ( jmaxInput <= ny1 ) then 
    nZones = 1;
  else
    nZones = 2;
  endif

  ! Set up variables
  allocate(zone(1:nzones))
  allocate(jmax(nzones))
  allocate(kmax(nzones))
  jMax = 0
  nx = imax-imin+1

  write(*,*) ' Reading grid data ...' 
  
  ! Allocate grid variables x and y 
  do m = 1,nzones
    kmax(m) = 1
  
    ! Get dimensions of zone
    call cg_zone_read_f(index_grid,index_base,m,zonename,isize,ier)
    nxTotal = isize(1,1)
    nyTotal = isize(2,1)

    ! Get the appropriate value of jmax
    if (m<=nzones/2) then
      jmax(m) = nyTotal
    else
      jmax(m) = jmaxInput - sum(jmax(:))
    endif
    
    if (nxTotal<nx) stop ' I-range inserted by user not valid'

    zone(m)%nx1 = nx
    zone(m)%nx2 = jmax(m)
    zone(m)%nx3 = 1

    allocate(zone(m)%x(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3))
    allocate(zone(m)%y(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3))

    write(*,'(A,I0,A,I0,A,I0)') '    Zone ',m, '                        - nxTotal = ', nxTotal, ' nyTotal = ', nyTotal
    write(*,'(A,I0,A,I0,A,I0)') '    Domain used in calculations  ' //' : nx      = ', zone(m)%nx1,      ' ny      = ', zone(m)%nx2

  enddo

  ! Read grid
  ijk_min(1) = imin
  ijk_min(2) = 1
  ijk_min(3) = 1
  ijk_max(1) = imax
  ijk_max(3) = 1
  do m = 1, nzones
    ijk_max(2) = jmax(m)
    call cg_coord_read_f(index_grid,index_base,m,'CoordinateX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%x(:,:,:),ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    call cg_coord_read_f(index_grid,index_base,m,'CoordinateY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%y(:,:,:),ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  enddo


  call cg_close_f(index_grid,ier)
  deallocate(isize)
  return
end subroutine read_grid_multi
!#################################################################################  

!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_soln_multi

  use mod_field
  use mod_CGNS
  use cgns
  
  implicit none
  integer(kind=4)           :: m, t, idummy
  integer(kind=4)           :: iFile
  integer                   :: dataType, dataDimensions, dataLength(12), pos
  real(kind=8), allocatable :: timeGlobal(:)
  real(kind=8)              :: tmp
  character(len=50)         :: tempName
  real(kind=8)              :: lastPrint, percentage
  integer                   :: ijk_min(3), ijk_max(3), iFlow

  write(*,*) achar(10) // '  Reading solution data ...' 

  nsnap = (idxf - idxi)/idxr + 1
  
  do m = 1,nzones
      
    allocate(zone(m)%r(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%u(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%v(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%w(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))
    allocate(zone(m)%p(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap))                

    zone(m)%r(:,:,:,:) = 0.0d0
    zone(m)%u(:,:,:,:) = 0.0d0
    zone(m)%v(:,:,:,:) = 0.0d0
    zone(m)%w(:,:,:,:) = 0.0d0
    zone(m)%p(:,:,:,:) = 0.0d0 
  
  enddo  
  allocate(time(nsnap))
  
  t = 0
  lastPrint = 0d0
  do idummy = idxi,idxf,idxr

    t = t + 1

    !---- reading the solution files
    write(CGNS_solnname,'(A,I6.6,A)') trim(path_to_soln) // '/' // trim(soln), idummy, '.cgns'
    call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    index_base = 1

    ! Get index of first FlowSolution inside zones
    call cg_sol_info_f(index_soln,index_base,1,1,solnname,pos,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

     ! Read time instant
    call cg_gopath_f(index_soln, '/Base/TimeIterValues', ier)
    if (ier .ne. CG_OK) call cg_error_exit_f 
    call cg_array_read_f(1,tmp,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    time(t) = tmp

    ! Read flow solution on the prompted locations
    ijk_min(1) = imin
    ijk_min(2) = 1
    ijk_min(3) = 1
    ijk_max(1) = imax
    ijk_max(3) = 1

    ! Read one zone at a time
    do m = 1, nZones
      ijk_max(2) = jmax(m)
      call cg_field_read_f(index_soln,index_base,m,1,'Density'  ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%r(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,1,'MomentumX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%u(:,:,:,t),ier)    
      call cg_field_read_f(index_soln,index_base,m,1,'MomentumY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%v(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,1,'MomentumZ',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%w(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,1,'Pressure' ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%p(:,:,:,t),ier)   
    enddo

    ! Counter 
    percentage = dble(t) / dble(nsnap) * 100d0
    if (percentage>lastPrint) then
      write(*,'(A,F0.0,A)') '  ',percentage, ' %'
      lastPrint = lastPrint + 10d0
    endif

    call cg_close_f(index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  end do

  percentage = 100d0
  write(*,'(A,F0.0,A)') '  ',percentage, ' %'

end subroutine read_soln_multi
!#################################################################################  

subroutine Convert_to_Primitive_Variables(flag)

  use mod_field
  implicit none
  logical, intent(in) :: flag
  integer(kind=4) :: t,m  

  if (flag .eqv. .true.) then

    call system('printf "\033[0;35m \n"')
    call system('printf "  REMARK: Working with VELOCITY variables... \n"')

    do m = 1,nzones
       
      do t = 1,nsnap
        zone(m)%u(:,:,:,t) = zone(m)%u(:,:,:,t)/zone(m)%r(:,:,:,t)
        zone(m)%v(:,:,:,t) = zone(m)%v(:,:,:,t)/zone(m)%r(:,:,:,t)
        zone(m)%w(:,:,:,t) = zone(m)%w(:,:,:,t)/zone(m)%r(:,:,:,t)
      enddo

      DEALLOCATE(zone(m)%r)
      
    enddo

  else

    call system('printf "\033[1;33m \n"')
    call system('printf "   REMARK: Working with MOMENTUM variables... \n"')
    call system('printf "\033[0m \n"')

  endif

  return

end subroutine Convert_to_Primitive_Variables
!#################################################################################
subroutine Reynolds_Decomposition

  use mod_field
  implicit none
  integer(kind=4) :: m,t

  write(*,*) achar(10) // '  Computing fluctuation data ...'   

  do m = 1,nzones
  
    zone(m)%r(:,:,:,0) = sum(zone(m)%r(:,:,:,1:nsnap),dim=4)/dble(nsnap)
    zone(m)%u(:,:,:,0) = sum(zone(m)%u(:,:,:,1:nsnap),dim=4)/dble(nsnap)
    zone(m)%v(:,:,:,0) = sum(zone(m)%v(:,:,:,1:nsnap),dim=4)/dble(nsnap)
    zone(m)%w(:,:,:,0) = sum(zone(m)%w(:,:,:,1:nsnap),dim=4)/dble(nsnap)
    zone(m)%p(:,:,:,0) = sum(zone(m)%p(:,:,:,1:nsnap),dim=4)/dble(nsnap)
    
    do t = 1,nsnap
      zone(m)%r(:,:,:,t) = zone(m)%r(:,:,:,t) - zone(m)%r(:,:,:,0)
      zone(m)%u(:,:,:,t) = zone(m)%u(:,:,:,t) - zone(m)%u(:,:,:,0)
      zone(m)%v(:,:,:,t) = zone(m)%v(:,:,:,t) - zone(m)%v(:,:,:,0)
      zone(m)%w(:,:,:,t) = zone(m)%w(:,:,:,t) - zone(m)%w(:,:,:,0)
      zone(m)%p(:,:,:,t) = zone(m)%p(:,:,:,t) - zone(m)%p(:,:,:,0)
    enddo
    
  enddo

  return

end subroutine Reynolds_Decomposition

!###############################################################################################

subroutine Correlation_Matrix

  use mod_pod_modes
  use mod_variable_precision
  use mod_field
  implicit none
 
  integer(kind=4) :: j, k, l, m
  integer(kind=4) :: t1, t2
  integer(kind=4) :: counter, lmax

  real(kind=VRL) :: summ
  real(kind=VRL), allocatable, dimension(:,:,:,:) :: corr_matrix(:,:,:,:)
  real(kind=VRL), allocatable, dimension(:,:,:,:) :: q1
  real(kind=VRL), allocatable, dimension(:,:,:,:) :: q2
  real(kind=VRL), allocatable, dimension(:,:) :: wk

  real(kind=8) percentage, lastPrint
  real(kind=8) start, end

  write(*,*) achar(10) // '  Computing correlation matrix ...'

  call cpu_time(start)

  lmax = 4 ! Number of flow variables
  
  ! Allocate correlation matrix
  allocate(corr_matrix(nsnap,nsnap,1:lmax,0:nzones)) 
  corr_matrix(:,:,:,:) = 0.0_VRL

  ! Write initial time 
  open(998, file='./SVD_files/initial_time.dat')
  write(998,'(9999es21.13)') time(1)
  close(998)

  ! Write the number of snapshots
  open(999,file='./SVD_files/nsnap') 
  write(999,*) nsnap
  close(999)

  counter = 0
  lastPrint = 10d0
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do m = 1,nzones

    ALLOCATE(q1(1:zone(m)%nx1,1:jmax(m),1:kmax(m),lmax)) ! Dummy variable to allocate q(i,j,k,ti)
    ALLOCATE(q2(1:zone(m)%nx1,1:jmax(m),1:kmax(m),lmax)) ! Dummy variable to allocate q(i,j,k,tj)
    ALLOCATE(wk(1:zone(m)%nx1,1:jmax(m))) ! Area

    wk(:,:) = 0.0d0
    

    do l = 1,lmax

      do t1 = 1,nsnap
        counter = counter + 1
        percentage = dble(counter) / dble(nzones*lmax*nsnap) * 100d0
        if (percentage>lastPrint) then
          write(*,'(A,F0.0,A)') '  ',percentage, ' %'
          lastPrint = lastPrint + 10d0
        endif

        wk(:,1:jmax(m)) = zone(m)%area(:,1:jmax(m))
        if (l .eq. 1) q1(:,1:jmax(m),1:kmax(m),l) = real( zone(m)%p(:,1:jmax(m),1:kmax(m),t1), VRL )
        if (l .eq. 2) q1(:,1:jmax(m),1:kmax(m),l) = real( zone(m)%u(:,1:jmax(m),1:kmax(m),t1), VRL )
        if (l .eq. 3) q1(:,1:jmax(m),1:kmax(m),l) = real( zone(m)%v(:,1:jmax(m),1:kmax(m),t1), VRL )
        if (l .eq. 4) q1(:,1:jmax(m),1:kmax(m),l) = real( zone(m)%w(:,1:jmax(m),1:kmax(m),t1), VRL )

        do j = 1,jmax(m)
          do k = 1,kmax(m)
            q1(:,j,k,l) = q1(:,j,k,l)*wk(:,j)
          enddo
        enddo


        do t2 = t1,nsnap

          if (l .eq. 1) q2(:,1:jmax(m),:,l) = real( zone(m)%p(:,1:jmax(m),:,t2),VRL )
          if (l .eq. 2) q2(:,1:jmax(m),:,l) = real( zone(m)%u(:,1:jmax(m),:,t2),VRL )
          if (l .eq. 3) q2(:,1:jmax(m),:,l) = real( zone(m)%v(:,1:jmax(m),:,t2),VRL )
          if (l .eq. 4) q2(:,1:jmax(m),:,l) = real( zone(m)%w(:,1:jmax(m),:,t2),VRL )

          summ = 0.0d0
          do k = 1,kmax(m)

            do j = 1,jmax(m)

              summ = summ + dot_product(q1(:,j,k,l),q2(:,j,k,l))
                              
            enddo ! j
          enddo ! k

          corr_matrix(t1,t2,l,m) = corr_matrix(t1,t2,l,m) + summ

        enddo ! t2

      enddo ! t1
                 
    enddo ! l

    DEALLOCATE(q1)
    DEALLOCATE(q2)
    DEALLOCATE(wk)
    
    corr_matrix(:,:,:,0) = corr_matrix(:,:,:,0) + corr_matrix(:,:,:,m)

  
  enddo ! m
  !end of the correlation matrix construction
  percentage = 100d0
  write(*,'(A,F0.0,A)') '  ',percentage, ' %'
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Exporting the correlation matrix     
  open(1000, file='./SVD_files/corr_matrix_p.dat')
  open(1001, file='./SVD_files/corr_matrix_k.dat')

  do t1 = 1,nsnap ! 

    !mirroring the correlation matrix
    !DIR$ NOVECTOR
    do t2 = 1,t1-1
      corr_matrix(t1,t2,1,0) = corr_matrix(t2,t1,1,0)
      corr_matrix(t1,t2,2,0) = corr_matrix(t2,t1,2,0)
      corr_matrix(t1,t2,3,0) = corr_matrix(t2,t1,3,0)
      corr_matrix(t1,t2,4,0) = corr_matrix(t2,t1,4,0)
    enddo

  enddo
  
  do t1 = 1, nsnap
    !~~~~~~~~~~~~~~~~~~~~~
    ! kinectic energy
    do t2 = 1,nsnap
      corr_matrix(t1,t2,2,0) = corr_matrix(t1,t2,2,0) + corr_matrix(t1,t2,3,0) + corr_matrix(t1,t2,4,0)
    enddo
    
    !~~~~~~~~~~~~~~~~~~~~~

    ! Check for double or single precision    
    if (VRL .eq. DBL) then
      write(1000,'(9999es21.13)') corr_matrix(t1,:,1,0)
      write(1001,'(9999es21.13)') corr_matrix(t1,:,2,0)
    endif
    
    if (VRL .eq. SNG) then
      write(1000,'(9999es14.6)') corr_matrix(t1,:,1,0)
      write(1001,'(9999es14.6)') corr_matrix(t1,:,2,0)
    endif
    !~~~~~~~~~~~~~~~~~~~~~

  enddo

  close(1000)
  close(1001)
  close(1002)

  call cpu_time(end)
  write(*,'(A,F0.1,A)') '  ', end-start, " s" 

  return
  
end subroutine Correlation_Matrix
!################################################################################################
subroutine Spatial_Modes_POD

  use mod_pod_modes
  use mod_CGNS
  use mod_variable_precision
  use mod_field
  implicit none
  integer(kind=4) :: i, j, k, m, n, counter
  real(kind=8) percentage, lastPrint
  character(len=stringLen) :: norm
  
  if (svdFlag) then
    write(*,*) achar(10) // '  Computing SVD ...'
    call system('python3 SVD.py')
  endif

  allocate(pod_zone(1:nzones))

  ! Eigenvalues of the correlation matrix 
  allocate(lambda(nsnap))
  open(10,file='./SVD_files/sigma.dat')
  
  read(10,*) norm
  read(10,*) nmodes
  do i = 1,nmodes
    read(10,*) lambda(i)
  enddo
  close(10)

  ! Eigenvectors of the correlation matrix (Temporal modes) 
  allocate(temporal_modes(nsnap,nmodes))
  open(30,file='./SVD_files/V_matrix.dat')
  read(30,*)
  read(30,*)
  read(30,*)
  read(30,*)
  do i = 1,nsnap
    read(30,*) temporal_modes(i,:)
  enddo
  close(30)

!##############################################################################
! Compute Temporal modes

  write(*,*) achar(10) // '  Computing spatial modes ...'
  
  counter = 0
  lastPrint = 10d0
  do m = 1,nzones

    ALLOCATE(pod_zone(m)%spatial_modes(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nmodes,1:5))
   
    pod_zone(m)%spatial_modes(:,:,:,:,:) = 0.0d0

    do i = 1,zone(m)%nx1
      percentage = dble(counter) / dble(zone(m)%nx1*nzones) * 100d0
      if (percentage>lastPrint) then
        write(*,'(A,F0.0,A)') '  ',percentage, ' %'
        lastPrint = lastPrint + 10d0
      endif
      counter =  counter+1
      do j = 1,jmax(m)
      do k = 1,kmax(m)
      do n = 1,nmodes
        pod_zone(m)%spatial_modes(i,j,k,n,1)  = dot_product(zone(m)%r(i,j,k,1:nsnap),temporal_modes(:,n))/lambda(n)
        pod_zone(m)%spatial_modes(i,j,k,n,2)  = dot_product(zone(m)%u(i,j,k,1:nsnap),temporal_modes(:,n))/lambda(n)
        pod_zone(m)%spatial_modes(i,j,k,n,3)  = dot_product(zone(m)%v(i,j,k,1:nsnap),temporal_modes(:,n))/lambda(n)
        pod_zone(m)%spatial_modes(i,j,k,n,4)  = dot_product(zone(m)%w(i,j,k,1:nsnap),temporal_modes(:,n))/lambda(n)
        pod_zone(m)%spatial_modes(i,j,k,n,5)  = dot_product(zone(m)%p(i,j,k,1:nsnap),temporal_modes(:,n))/lambda(n)
      enddo
      enddo
      enddo
    enddo

  enddo
  percentage = 100d0
  write(*,'(A,F0.0,A)') '  ',percentage, ' %'

  write(*,*) achar(10) // '  Writing cgns file ...'
  
  ! Write mean of variables and POD modes  
  write(*,*) achar(10) // '    Writing mean ...'
  call write_mean_2D_CGNS
  
  ! Create file for sum of POD modes
  
  write(*,*) achar(10) // '    Writing pod modes ...'
  call write_POD_modes_2D_CGNS()
  
  write(*,'(A,I0.0,A)') achar(10) // '    Writing sum of ', nmodes, ' POD_modes'
  call write_sum_POD_modes_2D_CGNS(nmodes)

  return

end subroutine Spatial_Modes_POD
!################################################################################################


!###############################################################################################

subroutine move_files()

  use mod_pod_modes
  use mod_CGNS
  use mod_variable_precision
  use mod_field
  implicit none 
  character(len=stringLen) :: temp
  character(len=stringLen) :: temp1

  write(*,*) achar(10) // '  Move files ...'

  write(temp,"(A)") 'cp -rf SVD_files/V_matrix.dat ' // trim(path_to_ROM) // 'regression/compact_scheme/data/'
  call system (temp)

  open(unit = 60, file = '/' // trim(path_to_ROM) //'regression/compact_scheme/compact_scheme_inputs.dat')
  write(60,'(A)') trim(path_to_ROM)
  write(60,'(I0)') der
  close(60)

  write(temp,"(A)") 'cp -rf inputs.inp ' // trim(path_to_ROM) // 'regression/deep_learning/data/'
  call system (temp)

  write(temp,"(A)") 'cp -rf modfile.f90 ' // trim(path_to_ROM) // 'reconst/'
  call system (temp)

end subroutine

!################################################################################################





