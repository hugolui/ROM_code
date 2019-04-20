PROGRAM reconst

  use mod_pod_modes
  use mod_variable_precision
  use mod_field
  use cgns
  use mod_CGNS
  implicit none

  !!! Read inputs !!!
  call read_inputs()

  !!! Read grid !!!
  call read_grid

  !!! Read mean quantities !!!
  call read_mean

  !!! Read spatial modes !!!
  call read_spatial_modes

  !!! Write the best flow Reconstruction in CGNS file!!!
  if (CGNSFlag) then
    write(*,*) ' Reconstructing the best candidate model ...'
    call reconst_flow_ROM_CGNS(nmodes,nsnap_reconst)
    stop
  end if

  !!! Read validation data !!!!
  if (soln_files .eq. 1) then
    call read_soln
  end if

  if (soln_files .eq. 2) then
    call read_soln_multi
  end if

  ! Reconstruction the flowfield for a determined candidate model
  call reconst_best_models(nmodes,nsnap_reconst,model)

  ! Compute the difference between Validation data and ROM solution
  call diff_FOM_ROM

END PROGRAM reconst

!#################################################################################

subroutine read_inputs()

    use mod_pod_modes
    use mod_field
    use mod_CGNS
    implicit none
    integer, parameter :: fileUnit = 55
    integer            :: POD_scheme
    character(len=4)   :: char

    open(fileUnit, file='inputs_reconst.inp')
    read(fileUnit,'(A)') path_to_soln
    read(fileUnit,*)     soln
    read(fileUnit,*)     soln_files
    read(fileUnit,'(A)') filename
    read(fileUnit,*)     model
    read(fileUnit,*)     idxi, idxf, idxr
    read(fileUnit,*)     idxi_val, idxf_val, idxr_val
    read(fileUnit,*)     jmaxInput
    read(fileUnit,*)     imin, imax
    read(fileUnit,*)     nsnap_reconst
    read(fileUnit,*)     dt
    read(fileUnit,*)     CGNSFlag

    !call remove_comments(outfile, stringLen)
    call remove_comments(path_to_soln, stringLen)
    call remove_comments(soln, stringLen)
    call remove_comments(filename, stringLen)
    call remove_comments(soln, stringLen)

    close(fileUnit)

    write(char,'(i4.4)') model
    write(outfile,'(A,A,A)') 'Candidate model '//char
    write(*,*) achar(10) // trim(outfile)

    if (CGNSFlag) then
      write(char,'(i4.4)') model
      write(outfile,'(A,A,A)') 'rebuiltSol_model_'//char//".cgns"
      write(*,*) achar(10) // trim(outfile)
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
subroutine read_grid

  use mod_pod_modes
  use mod_field
  use mod_CGNS
  use cgns
  implicit none
  integer(kind=4) :: m, nxTotal, nyTotal, ijk_min(3), ijk_max(3), nx, ny1, ny2

  index_base = 1
  call cg_open_f('./data/POD_modes.cgns',CG_MODE_READ,index_grid,ier)
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

end subroutine read_grid

!#################################################################################

!#################################################################################
!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_mean

  use mod_field
  use mod_CGNS
  use cgns

  implicit none
  integer(kind=4)           :: m
  integer                   :: ijk_min(3), ijk_max(3), iFlow

  write(*,*) achar(10) // '  Reading mean quantities ...'

  nsnap = (idxf - idxi)/idxr + 1

  do m = 1,nzones

    allocate(zone(m)%r(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:0))
    allocate(zone(m)%u(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:0))
    allocate(zone(m)%v(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:0))
    allocate(zone(m)%w(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:0))
    allocate(zone(m)%p(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:0))

    zone(m)%r(:,:,:,:) = 0.0d0
    zone(m)%u(:,:,:,:) = 0.0d0
    zone(m)%v(:,:,:,:) = 0.0d0
    zone(m)%w(:,:,:,:) = 0.0d0
    zone(m)%p(:,:,:,:) = 0.0d0

  enddo

  !---- reading the solution files
  call cg_open_f('./data/POD_modes.cgns',CG_MODE_READ,index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  index_base = 1

  ! Read flow solution on the prompted locations
  ijk_min(1) = imin
  ijk_min(2) = 1
  ijk_min(3) = 1
  ijk_max(1) = imax
  ijk_max(3) = 1

  iFlow = 1

  ! Read one zone at a time
  do m = 1, nZones
    ijk_max(2) = jmax(m)
    call cg_field_read_f(index_soln,index_base,m,iFlow,'Density'  ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%r(:,:,:,0),ier)
    call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%u(:,:,:,0),ier)
    call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%v(:,:,:,0),ier)
    call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumZ',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%w(:,:,:,0),ier)
    call cg_field_read_f(index_soln,index_base,m,iFlow,'Pressure' ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%p(:,:,:,0),ier)
  enddo

  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

end subroutine read_mean

!#################################################################################

subroutine read_spatial_modes

  use mod_field
  use mod_CGNS
  use cgns
  use mod_pod_modes
  implicit none
  integer(kind=4)           :: m, i
  integer                   :: ijk_min(3), ijk_max(3), iFlow
  character(len=1) :: norm
  character(len=500) :: char

  write(*,*) achar(10) // '  Reading spatial modes ...'

  !!! Read Eigenvalues !!!
  allocate(lambda(nsnap))
  open(10,file='./data/sigma.dat')
  read(10,*) norm
  read(10,*) nmodes
  do i = 1,nmodes
    read(10,*) lambda(i)
  enddo
  close(10)

  !!! Read temporal modes (Reconstruction) !!!
  allocate(modos_temporais(nsnap_reconst,nmodes))
  write(char,'(A,A)')  '.'//trim(filename)
  open(30,file=trim(char))
  do i = 1,nsnap_reconst
    read(30,*) modos_temporais(i,:)
  enddo
  close(30)

  !---- reading the solution files
  call cg_open_f('./data/POD_modes.cgns',CG_MODE_READ,index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  index_base = 1

  if (allocated(isize) .eqv. .false.) allocate(isize(isizeDim,2))
  isize(:,:) = 1

  allocate(pod_zone(1:nzones))

  ijk_min(1) = imin
  ijk_min(2) = 1
  ijk_min(3) = 1
  ijk_max(1) = imax
  ijk_max(3) = 1

  do m = 1,nzones

    ALLOCATE(pod_zone(m)%spatial_modes(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nmodes,1:5))
    pod_zone(m)%spatial_modes(:,:,:,:,:) = 0.0d0

  end do

  do i=1, nmodes
    write(citer,'(i3.3)') i
    do m = 1,nzones

      isize(1,1) = zone(m)%nx1
      isize(2,1) = zone(m)%nx2
      
      isize(1,2) = isize(1,1) - 1
      isize(2,2) = isize(2,1) - 1  

      ijk_max(2) = jmax(m)
      
      write(zonename,'(a5,i1,a5)') 'Zone_', m, '_0001'

      call cg_goto_f(index_soln, index_base, ier, 'Zone_t', m, 'end')
      
      call cg_goto_f(index_soln, index_base, ier, 'Zone_t', m, "FlowSolution_t", 1, 'end')
      
      call cg_field_read_f(index_soln,index_base,m,1, 'R_' //citer,RealDouble,ijk_min(1:3),ijk_max(1:3),pod_zone(m)%spatial_modes(:,:,:,i,1),ier)
      call cg_field_read_f(index_soln,index_base,m,1, 'Mx_' //citer,RealDouble,ijk_min(1:3),ijk_max(1:3),pod_zone(m)%spatial_modes(:,:,:,i,2),ier)
      call cg_field_read_f(index_soln,index_base,m,1, 'My_' //citer,RealDouble,ijk_min(1:3),ijk_max(1:3),pod_zone(m)%spatial_modes(:,:,:,i,3),ier)
      call cg_field_read_f(index_soln,index_base,m,1, 'Mz_' //citer,RealDouble,ijk_min(1:3),ijk_max(1:3),pod_zone(m)%spatial_modes(:,:,:,i,4),ier)
      call cg_field_read_f(index_soln,index_base,m,1, 'P_' //citer,RealDouble,ijk_min(1:3),ijk_max(1:3),pod_zone(m)%spatial_modes(:,:,:,i,5),ier)

      if (ier .ne. CG_OK) call cg_error_exit_f
      
     enddo
   enddo

  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

end subroutine read_spatial_modes

! !#################################################################################

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

  write(*,*) achar(10) // '  Reading validation data ...' 

  nsnap_val = (idxf_val - idxi_val)/idxr_val + 1

  do m = 1,nzones
      
    allocate(zone(m)%r_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_val))
    allocate(zone(m)%u_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_val))
    allocate(zone(m)%v_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_val))
    allocate(zone(m)%w_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_val))
    allocate(zone(m)%p_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_val))                

    zone(m)%r_val(:,:,:,:) = 0.0d0
    zone(m)%u_val(:,:,:,:) = 0.0d0
    zone(m)%v_val(:,:,:,:) = 0.0d0
    zone(m)%w_val(:,:,:,:) = 0.0d0
    zone(m)%p_val(:,:,:,:) = 0.0d0 
  
  enddo  
  allocate(time(nsnap_val))
  
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
  do iFile = idxi_val,idxf_val,idxr_val

    t = t + 1
    time(t) = timeGlobal(iFile-iFirstFlowSol)

    iFlow = iFile - iFirstFlowSol

    ! Read one zone at a time
    do m = 1, nZones
      ijk_max(2) = jmax(m)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'Density'  ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%r_val(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%u_val(:,:,:,t),ier)    
      call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%v_val(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'MomentumZ',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%w_val(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,iFlow,'Pressure' ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%p_val(:,:,:,t),ier)   
    enddo

    ! Counter 
    percentage = dble(t) / dble(nsnap_val) * 100d0
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

  write(*,*) achar(10) // '  Reading validation data ...' 

  nsnap_val = (idxf_val - idxi_val)/idxr_val + 1
  
  do m = 1,nzones
      
    allocate(zone(m)%r_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nsnap_val))
    allocate(zone(m)%u_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nsnap_val))
    allocate(zone(m)%v_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nsnap_val))
    allocate(zone(m)%w_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nsnap_val))
    allocate(zone(m)%p_val(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nsnap_val))                

    zone(m)%r_val(:,:,:,:) = 0.0d0
    zone(m)%u_val(:,:,:,:) = 0.0d0
    zone(m)%v_val(:,:,:,:) = 0.0d0
    zone(m)%w_val(:,:,:,:) = 0.0d0
    zone(m)%p_val(:,:,:,:) = 0.0d0 
  
  enddo  
  allocate(time(nsnap_val))
  
  t = 0
  lastPrint = 0d0

  do idummy = idxi_val,idxf_val,idxr_val

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
      call cg_field_read_f(index_soln,index_base,m,1,'Density'  ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%r_val(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,1,'MomentumX',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%u_val(:,:,:,t),ier)    
      call cg_field_read_f(index_soln,index_base,m,1,'MomentumY',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%v_val(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,1,'MomentumZ',RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%w_val(:,:,:,t),ier)
      call cg_field_read_f(index_soln,index_base,m,1,'Pressure' ,RealDouble,ijk_min(1:3),ijk_max(1:3),zone(m)%p_val(:,:,:,t),ier)   
    enddo

    ! Counter 
    percentage = dble(t) / dble(nsnap_val) * 100d0
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

subroutine reconst_best_models(nOutModes,nsnap_S,num_model)

  use mod_field
  use mod_pod_modes

  implicit none
  integer(kind=4), intent(in)  :: nsnap_S
  integer(kind=4), intent(in)  :: nOutModes
  integer(kind=4), intent(in)  :: num_model
  real(kind=8), allocatable    :: temp(:,:,:)
  integer(kind=4)              :: m, n, l, nn, ll, i
  real(kind=8), allocatable    :: temporal_modes_model(:,:)
  character(len=500) :: char

  write(*,*) ' Reconstructing the Flowfield ...'

  !!! Read temporal modes (Reconstruction) !!!
  allocate(temporal_modes_model(nsnap_S,nOutModes))
  write(char,'(A,A)')  '.'//trim(filename)
  open(35,file=trim(char))
  do i = 1,nsnap_S
    read(35,*) temporal_modes_model(i,:)
  enddo
  close(35)

  !!! Allocate validation data !!!
  do m = 1,nzones
    
  allocate(zone(m)%r_reconst(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_S))
  allocate(zone(m)%u_reconst(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_S))
  allocate(zone(m)%v_reconst(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_S))
  allocate(zone(m)%w_reconst(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_S))
  allocate(zone(m)%p_reconst(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,0:nsnap_S))                

  zone(m)%r_reconst(:,:,:,:) = 0.0d0
  zone(m)%u_reconst(:,:,:,:) = 0.0d0
  zone(m)%v_reconst(:,:,:,:) = 0.0d0
  zone(m)%w_reconst(:,:,:,:) = 0.0d0
  zone(m)%p_reconst(:,:,:,:) = 0.0d0 
  
  enddo  

  do m = 1,nzones

    allocate(temp(zone(m)%nx1,1:zone(m)%nx2,zone(m)%nx3))

    nn = 0
  ! Save a flow for every timestep
    do n = 1,nsnap_S

      nn = nn + 1
      ! Density
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,1)*temporal_modes_model(nn,l)*lambda(l) + temp
      enddo
      
      zone(m)%r_reconst(:,:,:,nn) = temp + zone(m)%r(:,:,:,0)

      ! X-Momentum
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,2)*temporal_modes_model(nn,l)*lambda(l) + temp
      enddo

      zone(m)%u_reconst(:,:,:,nn) = temp + zone(m)%u(:,:,:,0)

      ! Y-Momentum
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,3)*temporal_modes_model(nn,l)*lambda(l) + temp
      enddo

      zone(m)%v_reconst(:,:,:,nn) = temp + zone(m)%v(:,:,:,0)


      ! Z-Momentum
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,4)*temporal_modes_model(nn,l)*lambda(l) + temp
      enddo

      zone(m)%w_reconst(:,:,:,nn) = temp + zone(m)%w(:,:,:,0)

      ! Pressure
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,5)*temporal_modes_model(nn,l)*lambda(l) + temp
      enddo

      zone(m)%p_reconst(:,:,:,nn) = temp + zone(m)%p(:,:,:,0)

    enddo

    deallocate(temp)

  end do
  

end subroutine reconst_best_models

! !#################################################################################

subroutine diff_FOM_ROM

  use mod_field
  use mod_pod_modes
  
  implicit none
  real(kind=8), allocatable :: diff_p(:,:,:,:)
  real(kind=8) :: diff, L1_error
  integer(kind=4) :: m

  write(*,*) 'Computing ||p_FOM - p_ROM || ...'

  L1_error = 0.0d0
  do m = 1,nzones

    allocate(zone(m)%diff_p(1:zone(m)%nx1,1:zone(m)%nx2,1:zone(m)%nx3,1:nsnap_val))
    zone(m)%diff_p = abs(zone(m)%p_val(:,:,:,:) - zone(m)%p_reconst(:,:,:,nsnap+1:nsnap+nsnap_val))
    L1_error = L1_error + sum(zone(m)%diff_p)

  end do

  open(unit = 37, file = './data/L1_error.dat',position= 'append', action='readwrite')
  write(37,'(9999es21.13,9999es21.13)') L1_error
  close(37)

end subroutine diff_FOM_ROM

! !#################################################################################

subroutine reconst_flow_ROM_CGNS(nOutModes,nsnap_S)

  use cgns
  use mod_CGNS
  use mod_field
  use mod_pod_modes

  implicit none
  integer(kind=4), intent(in)  :: nsnap_S
  integer(kind=4), intent(in)  :: nOutModes
  integer(kind=4)              :: m, n, l, nn, ll
  integer(kind=4), parameter   :: nCharacters = 32
  character(len=nCharacters)   :: flowname
  character(len=nCharacters), allocatable :: flownameVector(:)
  character(len=50)            :: tempString
  real(kind=8), allocatable    :: temp(:,:,:)
  real(kind=8), allocatable    :: time_S(:)
  real(kind=8)                 :: t0
  character(len=1) :: norm

  write(*,*) achar(10) // '  Write CGNS file ...'

  print *, 'Number of snapshots for reconstruction:', nsnap_reconst

  ! Read Initial time
  open(17,file='./data/initial_time.dat')
  read(17,*) t0
  close(17)

  ! Create time array
  allocate(flownameVector(nsnap_S))
  allocate(time_S(nsnap_S))
  time_S = 0.0d0
  do ll = 1,nsnap_S
    time_S(ll) = t0 + dt*(ll-1)
  end do

  if (allocated(isize) .eqv. .false.) allocate(isize(isizeDim,2))

  call cg_open_f(trim(outfile),CG_MODE_WRITE,index_soln,ier)
  ! Create a base
  call cg_base_write_f(index_soln,'Base',2,3, index_base,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  do m = 1,nzones

    isize(1,1)=zone(m)%nx1;      isize(2,1)=zone(m)%nx2;
    isize(1,2)=isize(1,1)-1;    isize(2,2)=isize(2,1)-1;

    write(zonename,'(a5,i1,a5)') 'Zone_', m, '_0001'
    allocate(temp(isize(1,1),isize(2,1),1))

    ! Create zone
    call cg_zone_write_f(index_soln,index_base,trim(zonename),isize,Structured,index_zone,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    nn = 0
    ! Save a flow for every timestep
    do n = 1,nsnap_S

      nn = nn + 1

      write(flowname,'(A,I0.0)') 'FlowSolution_', idxi + (n-1)*idxr  
      flownameVector(nn) = trim(flowname)

      call cg_sol_write_f(index_soln,index_base,m,trim(flowname),Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

      ! R
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,1)*modos_temporais(nn,l)*lambda(l) + temp
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'Density',temp + zone(m)%r(:,:,:,0) ,index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

      ! RU
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,2)*modos_temporais(nn,l)*lambda(l) + temp
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumX',temp + zone(m)%u(:,:,:,0) ,index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

      ! RV
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,3)*modos_temporais(nn,l)*lambda(l)  + temp
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumY',temp + zone(m)%v(:,:,:,0),index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

      ! RW
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,4)*modos_temporais(nn,l)*lambda(l)  + temp
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumZ',temp + zone(m)%w(:,:,:,0),index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

      ! P
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,5)*modos_temporais(nn,l)*lambda(l)  + temp
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'Pressure',temp + zone(m)%p(:,:,:,0),index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    enddo

    call cg_ziter_write_f(index_soln,index_base, m,'ZoneIterativeData',ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    call cg_goto_f(index_soln,index_base,ier,'Zone_t',m,'ZoneIterativeData_t',1,'end')
    if (ier .ne. CG_OK) call cg_error_exit_f

    call cg_array_write_f('FlowSolutionPointers',Character,2,(/ nCharacters, nsnap_S /), flownameVector, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    ! x-coordinate
    call cg_coord_write_f(index_soln,index_base,m, &
                          RealDouble, 'CoordinateX',zone(m)%x,index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    ! y-coordinate
    call cg_coord_write_f(index_soln,index_base,m, &
                          RealDouble,'CoordinateY',zone(m)%y,index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    ! z-coordinate
    temp = 0d0
    call cg_coord_write_f(index_soln,index_base,m, &
                         RealDouble,'CoordinateZ',temp,index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

    deallocate(temp)
  enddo ! nzones

  ! Associate each flow solution to a timestep
  call cg_gopath_f(index_soln, '/Base', ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_biter_write_f(index_soln,index_base,'TimeIterValues',nsnap_S,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_goto_f(index_soln,index_base,ier,'BaseIterativeData_t',1,'end')
  if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_array_write_f('TimeValues',RealDouble,1,nsnap_S,time_S,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_simulation_type_write_f(index_soln,index_base,TimeAccurate,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  ! Create descriptor so that the user may know how many modes were using in the reconstruction
  call cg_goto_f(index_soln,index_base,ier,'end')
  write(tempString,"(I0.0,A)") nmodes, ' used in reconstruction'
  call cg_descriptor_write_f('Additional info',trim(tempString),ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

end subroutine
