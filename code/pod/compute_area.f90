!#################################################################################
subroutine compute_cell_area
 
  ! the area is for a ficticious cell-vertex element
  ! the "center" is the node from the finite difference scheme
  ! the total area is the sum of portions of the "cell center" elements surrounding

  !this can be easily expanded to 3D by using an extrusion!!
  !   
  !   volume = area*dz

  use mod_field
  implicit none
  integer(kind=4) :: i!,i1,i2
  integer(kind=4) :: j!,j1,j2
  integer(kind=4) :: m
  
  real(kind=8) xaux(4), yaux(4), area

  print*, achar(10) // '  Computing cell area ...'   

  do m = 1,nzones

    allocate(zone(m)%area(1:zone(m)%nx1,1:zone(m)%nx2))
    zone(m)%area(:,:) = 0.0d0

    do j = 1,zone(m)%nx2-1
  
      do i = 1,zone(m)%nx1-1
      
        xaux(1) = zone(m)%x( i , j ,1)
        xaux(2) = zone(m)%x(i+1, j ,1)
        xaux(3) = zone(m)%x(i+1,j+1,1)
        xaux(4) = zone(m)%x( i ,j+1,1)
        
        yaux(1) = zone(m)%y( i , j ,1)
        yaux(2) = zone(m)%y(i+1, j ,1)
        yaux(3) = zone(m)%y(i+1,j+1,1)
        yaux(4) = zone(m)%y( i ,j+1,1)
        
        call compute_quad_area(xaux(1:4),yaux(1:4),area)       
        
        zone(m)%area( i , j ) = zone(m)%area( i , j ) + area*0.25d0
        zone(m)%area(i+1, j ) = zone(m)%area(i+1, j ) + area*0.25d0
        zone(m)%area(i+1,j+1) = zone(m)%area(i+1,j+1) + area*0.25d0
        zone(m)%area( i ,j+1) = zone(m)%area( i ,j+1) + area*0.25d0

              
      enddo
      
    enddo

    zone(m)%area(zone(m)%nx1,:) = zone(m)%area(zone(m)%nx1,:) + zone(m)%area(1,:)
    zone(m)%area( 1 ,:) = zone(m)%area(zone(m)%nx1,:)    

    zone(m)%area(:, 1 ) = 2.0d0*zone(m)%area(:, 1 )
    zone(m)%area(:,zone(m)%nx2) = 2.0d0*zone(m)%area(:,zone(m)%nx2)    
    
  enddo

  ! open(99,file='fort.dat')
  ! write(99,'(A)') 'VARIABLES="X","Y","S"'
  ! do m = 1,2

  !   nx1 = zone(m)%nx
  !   nx2 = zone(m)%ny
  !   nx3 = zone(m)%nz

  !   write(99,'(A,i0,A,i0,A)') 'ZONE T="",I=',nx1,',J=',nx2,',F=POINT'

  !   do j = 1,nx2
  !     do i = 1,nx1
  !       write(99,*) zone(m)%x(i,j,1),zone(m)%y(i,j,1),zone(m)%area(i,j)
  !     enddo
  !   enddo

  ! enddo

  return

end subroutine compute_cell_area
!=================================================================================

!=================================================================================
!
!   Area of a quadrilateral element according to the right-hand rule:
!
!                    A = index 1
!                    B = index 2
!                    C = index 3
!                    D = index 4
!
!              D ____________________ C 
!               /                    |
!              /                     |
!             /                      |
!            /                       |
!           /________________________|
!         A                           B
!       
!
!   https://en.wikipedia.org/wiki/Quadrilateral#Vector_formulas
!
subroutine compute_quad_area(x,y,area)

  implicit none
  real(kind=8), intent(in)  :: x(4), y(4)
  real(kind=8), intent(out) :: area
  real(kind=8) x1,x2,y1,y2

  x1 = x(3) - x(1)
  y1 = y(3) - y(1)
  
  x2 = x(4) - x(2)
  y2 = y(4) - y(2)

  area = 0.5*dabs(x1*y2 - x2*y1)
  
  return

end subroutine
!#################################################################################
