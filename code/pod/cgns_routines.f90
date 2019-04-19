!------------------------------------------------------------------------------------------------
subroutine write_POD_modes_2D_CGNS()

  use cgns
  use mod_CGNS
  use mod_field
  use mod_pod_modes
  implicit none
  integer(kind=4) :: m
  integer(kind=4) :: i
  
  call cg_open_f(trim(outfile),CG_MODE_MODIFY,index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  if (allocated(isize) .eqv. .false.) allocate(isize(isizeDim,2))
  isize(:,:) = 1

  do i=1, nmodes
    write(citer,'(i3.3)') i
    do m = 1,nzones
    
      isize(1,1) = zone(m)%nx1
      isize(2,1) = zone(m)%nx2
      
      isize(1,2) = isize(1,1) - 1
      isize(2,2) = isize(2,1) - 1  
      
      write(zonename,'(a5,i1,a5)') 'Zone_', m, '_0001'

      call cg_goto_f(index_soln, index_base, ier, 'Zone_t', m, 'end')
      
      call cg_goto_f(index_soln, index_base, ier, 'Zone_t', m, "FlowSolution_t", 1, 'end')
      
      call cg_field_write_f(index_soln,index_base,m,1, &
                            RealDouble, 'R_' //citer,pod_zone(m)%spatial_modes(:,:,:,i,1),index_field,ier)
      call cg_field_write_f(index_soln,index_base,m,1, &
                            RealDouble, 'Mx_'//citer,pod_zone(m)%spatial_modes(:,:,:,i,2),index_field,ier)
      call cg_field_write_f(index_soln,index_base,m,1, &
                            RealDouble, 'My_'//citer,pod_zone(m)%spatial_modes(:,:,:,i,3),index_field,ier)
      call cg_field_write_f(index_soln,index_base,m,1, &
                            RealDouble, 'Mz_'//citer,pod_zone(m)%spatial_modes(:,:,:,i,4),index_field,ier)
      call cg_field_write_f(index_soln,index_base,m,1, &
                            RealDouble, 'P_' //citer,pod_zone(m)%spatial_modes(:,:,:,i,5),index_field,ier)

      if (ier .ne. CG_OK) call cg_error_exit_f
      
     enddo
   enddo
  
  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f    

end subroutine write_POD_modes_2D_CGNS
!------------------------------------------------------------------------------------------------

subroutine write_mean_2D_CGNS

  use cgns
  use mod_CGNS
  use mod_field
  use mod_pod_modes
  implicit none
  integer(kind=4) :: m
  real(kind=8), allocatable :: temp(:,:,:)

  call cg_open_f(trim(outfile),CG_MODE_WRITE,index_soln,ier)

  ! Create a base
  call cg_base_write_f(index_soln,'Base',2,3, index_base,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  
  if (allocated(isize) .eqv. .false.) allocate(isize(isizeDim,2))
  isize = 1

  do m = 1,nzones
  
    isize(1,1)=zone(m)%nx1;      isize(2,1)=zone(m)%nx2; 
    isize(1,2)=isize(1,1)-1;    isize(2,2)=isize(2,1)-1; 
    
    write(zonename,'(a5,i1,a5)') 'Zone_', m, '_0001'
    allocate(temp(isize(1,1),isize(2,1),1))
    
    ! Create zone
    call cg_zone_write_f(index_soln,index_base,trim(zonename),isize,Structured,index_zone,ier)
    
    call cg_sol_write_f(index_soln,index_base,m,"FlowSolution",Vertex,index_flow,ier)

    call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'Density',zone(m)%r(:,:,:,0),index_field,ier)
    call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumX',zone(m)%u(:,:,:,0),index_field,ier)
    call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumY',zone(m)%v(:,:,:,0),index_field,ier)
    call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumZ',zone(m)%w(:,:,:,0),index_field,ier)
    call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'Pressure',zone(m)%p(:,:,:,0),index_field,ier)  
    call cg_coord_write_f(index_soln,index_base,m, &
                          RealDouble,'CoordinateX',zone(m)%x,index_field,ier)
    call cg_coord_write_f(index_soln,index_base,m, & 
                          RealDouble,'CoordinateY',zone(m)%y,index_field,ier)
    temp = 0d0
    call cg_coord_write_f(index_soln,index_base,m, &
                         RealDouble,'CoordinateZ',temp,index_field,ier)
                      
    if (ier .ne. CG_OK) call cg_error_exit_f
    deallocate(temp)    
   enddo
  
  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f    

end subroutine
!------------------------------------------------------------------------------------------------

subroutine write_sum_POD_modes_2D_CGNS(nOutModes)

  use cgns
  use mod_CGNS
  use mod_field
  use mod_pod_modes
  
  implicit none
  integer(kind=4), intent(in)  :: nOutModes
  integer(kind=4)              :: m, n, l, nn
  integer(kind=4), parameter   :: nCharacters = 32
  character(len=nCharacters)   :: flowname
  character(len=nCharacters)   :: flownameVector(nsnap)
  character(len=50)            :: tempString
  real(kind=8), allocatable    :: temp(:,:,:)
  

  call cg_open_f(trim(outfile2),CG_MODE_WRITE,index_soln,ier)
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
    do n=idxi,idxf, idxr 
      nn = nn + 1
      write(flowname,'(A,I0.0)') 'FlowSolution_', n
      flownameVector(nn) = trim(flowname)
      
      call cg_sol_write_f(index_soln,index_base,m,trim(flowname),Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
      ! R 
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,1)*temporal_modes(nn,l)*lambda(l) + temp  
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'Density',temp + zone(m)%r(:,:,:,0) ,index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
      ! RU 
      temp = 0d0
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,2)*temporal_modes(nn,l)*lambda(l) + temp  
      enddo
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumX',temp + zone(m)%u(:,:,:,0) ,index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
                    
      ! RV   
      temp = 0d0   
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,3)*temporal_modes(nn,l)*lambda(l)  + temp  
      enddo            
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumY',temp + zone(m)%v(:,:,:,0),index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
      ! RW   
      temp = 0d0  
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,4)*temporal_modes(nn,l)*lambda(l)  + temp  
      enddo  
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'MomentumZ',temp + zone(m)%w(:,:,:,0),index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
      ! P   
      temp = 0d0   
      do l=1, nOutModes
        temp = pod_zone(m)%spatial_modes(:,:,:,l,5)*temporal_modes(nn,l)*lambda(l)  + temp  
      enddo  
      call cg_field_write_f(index_soln,index_base,m,index_flow, &
                          RealDouble,'Pressure',temp + zone(m)%p(:,:,:,0),index_field,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    enddo

    
    call cg_ziter_write_f(index_soln,index_base, m,'ZoneIterativeData',ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    
    call cg_goto_f(index_soln,index_base,ier,'Zone_t',m,'ZoneIterativeData_t',1,'end')
    if (ier .ne. CG_OK) call cg_error_exit_f
    
    call cg_array_write_f('FlowSolutionPointers',Character,2,(/ nCharacters, nsnap /), flownameVector, ier)
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
  
  call cg_biter_write_f(index_soln,index_base,'TimeIterValues',nsnap,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f
  
  call cg_goto_f(index_soln,index_base,ier,'BaseIterativeData_t',1,'end')
  if (ier .ne. CG_OK) call cg_error_exit_f 
  
  call cg_array_write_f('TimeValues',RealDouble,1,nsnap,time,ier)
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
