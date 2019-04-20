
!#####################################################################################
module mod_variable_precision

  integer(kind=4), parameter :: SNG = 4
  integer(kind=4), parameter :: DBL = 8
  
  integer(kind=4), parameter :: VRL = 8

end module
!#####################################################################################
module mod_field

  integer(kind=4), parameter :: stringLen = 500
  character(len=stringLen) :: path_to_grid ! Path to grid file
  character(len=stringLen) :: path_to_soln ! Path to solution files
  character(len=stringLen) :: path_to_ROM ! Origin folder
  character(len=stringLen) :: soln ! Name of solution files
  character(len=stringLen) :: outfile ! Name of cgns file with POD modes
  character(len=stringLen) :: outfile2 ! Name of cgns file with reconstructed solution
  character(len=stringLen) :: filename ! Name of the temporal mode file
  
  integer(kind=4)  :: idxi      ! index initial (Training data)
  integer(kind=4)  :: idxf      ! index final (Training data)
  integer(kind=4)  :: idxr      ! index rate (Training data)
  integer(kind=4)  :: idxi_val      ! index initial (Validation data)
  integer(kind=4)  :: idxf_val      ! index final (Validation data)
  integer(kind=4)  :: idxr_val      ! index rate (Validation data)
  integer(kind=4)  :: nsnap     ! Number of snapshots (Training data)
  integer(kind=4)  :: nsnap_val ! Number of snapshots (validation data)
  integer(kind=4)  :: nzones, jmaxInput ! Number of zones and Number of points in j-th direction
  integer(kind=4)  :: der ! Derivative order
  integer(kind=4)  :: soln_files ! CGNS soln in a single file or CGNS soln in multiples files
  integer(kind=4)  :: model      ! Candidate model
  
  
  logical            :: logical_primitive ! Use primitive variables (T/F)
  logical            :: corrFlag, svdFlag ! Calculate correlation matrix (T/F) and Calculate SVD (T/F)
  logical            :: CGNSFlag ! Write CGNS file (T/F)

  real(kind=8) :: dt ! Time step
  real(kind=8), allocatable     :: time(:) ! Time array
  integer(kind=4), allocatable, dimension(:) :: jmax ! Number of points in the j-th direction of each zone
  integer(kind=4), allocatable, dimension(:) :: kmax ! Number of points in the k-th direction of each zone
  integer(kind=4) :: imax, imin ! Range of points in i-th direction
  type pinto

    real(kind=8), allocatable, dimension(:,:,:,:) :: r ! Density (Training data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: u ! X-momentum (Training data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: v ! Y-momentum (Training data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: w ! Z-momentum (Training data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: p ! Pressure (Training data)

    real(kind=8), allocatable, dimension(:,:,:,:) :: r_val ! Density (Validation data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: u_val ! X-momentum (Validation data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: v_val ! Y-momentum (Validation data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: w_val ! Z-momentum (Validation data)
    real(kind=8), allocatable, dimension(:,:,:,:) :: p_val ! Pressure (Validation data)

    real(kind=8), allocatable, dimension(:,:,:,:) :: r_reconst ! Density (Reconstruction)
    real(kind=8), allocatable, dimension(:,:,:,:) :: u_reconst ! X-momentum (Reconstruction)
    real(kind=8), allocatable, dimension(:,:,:,:) :: v_reconst ! Y-momentum (Reconstruction)
    real(kind=8), allocatable, dimension(:,:,:,:) :: w_reconst ! Z-momentum (Reconstruction)
    real(kind=8), allocatable, dimension(:,:,:,:) :: p_reconst ! Pressure (Reconstruction)

    real(kind=8), allocatable, dimension(:,:,:,:) :: diff_p ! Difference between p_val and p_reconst

    real(kind=8), allocatable, dimension(:,:,:) :: x, y, z ! Grid points
    
    real(kind=8), allocatable :: area(:,:) ! Area
    
    integer(kind=4) :: nx, ny, nz 
    integer(kind=4) :: nx1, nx2, nx3 ! Number of points in i-th, j-th and k-th directions
    
  end type pinto
  type(pinto), allocatable :: zone(:) ! Zone array

end module

!#####################################################################################
module mod_CGNS
  use mod_field
  character(len=stringLen) :: CGNS_filename
  character(len=stringLen) :: CGNS_solnname   

  integer(kind=4) :: index_node, index_grid, index_soln, index_field
  integer(kind=4) :: index_base, index_zone, index_flow, index_coord
  
  integer(kind=4) ier
  
  integer(kind=4), allocatable :: isize(:,:)
  integer(kind=4)              :: isizeDim

  character(len=16) :: zonename  
  character(len=stringLen)  :: solnname
  character(len=30)  :: citer
  
  character(len=5)  :: cmode

  integer(kind=4) :: iFirstFlowSol

end module


!#####################################################################################
module mod_pod_modes

  use mod_variable_precision

  !~~~~~ sparsity!
  integer(kind=4), parameter :: istride = 1
  integer(kind=4), parameter :: jstride = 1
  integer(kind=4), parameter :: kstride = 1

  !~~~~~ 
  integer(kind=4) nmodes ! Number of POD modes

  real(kind=VRL), allocatable :: lambda(:) ! Eigenvalues of the correlation matrix
  real(kind=VRL), allocatable :: temporal_modes(:,:) ! Temporal modes (Training data)

  real(kind=8), allocatable , dimension(:,:) :: modos_temporais ! Temporal modes (Reconstruction)
  integer(kind=4) :: nsnap_reconst ! Number of snanpshots (Reconstruction)
  type pod_zones
    real(kind=VRL), allocatable :: spatial_modes(:,:,:,:,:) ! Spatial modes
  end type pod_zones
  type(pod_zones), allocatable :: pod_zone(:) ! Zone array
  
end module
!#####################################################################################
