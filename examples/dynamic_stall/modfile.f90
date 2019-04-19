
!#####################################################################################
module mod_variable_precision

  integer(kind=4), parameter :: SNG = 4
  integer(kind=4), parameter :: DBL = 8
  
  integer(kind=4), parameter :: VRL = 8

end module
!#####################################################################################
module mod_field

  integer(kind=4), parameter :: stringLen = 500
  character(len=stringLen) :: path_to_grid
  character(len=stringLen) :: path_to_soln
  character(len=stringLen) :: soln
  character(len=stringLen) :: outfile
  character(len=stringLen) :: filename
  character(len=5) :: var
  integer(kind=4)  :: soln_files

  integer(kind=4)  :: idxi      ! index initial
  integer(kind=4)  :: idxf      ! index final
  integer(kind=4)  :: idxr      ! index rate
  integer(kind=4)   :: nsnap
  integer(kind=4)   :: nzones, jmaxInput
  integer(kind=4)  :: model      ! model


  logical            :: logical_primitive
  logical            :: removePlunge, corrFlag, svdFlag, additionalCorrelationFlag
  
  
  real(kind=8)     :: yFreq

  real(kind=8), allocatable     :: time(:)
  !real(kind=8), allocatable     :: time_S(:)
  integer(kind=4), allocatable, dimension(:) :: jmax
  integer(kind=4), allocatable, dimension(:) :: kmax
  integer(kind=4) :: imax, imin
  integer(kind=4) :: nx

  type pinto

    real(kind=8), allocatable, dimension(:,:,:,:) :: r
    real(kind=8), allocatable, dimension(:,:,:,:) :: u
    real(kind=8), allocatable, dimension(:,:,:,:) :: v
    real(kind=8), allocatable, dimension(:,:,:,:) :: w
    real(kind=8), allocatable, dimension(:,:,:,:) :: p
        
    real(kind=8), allocatable, dimension(:,:,:) :: x, y, z
    
    real(kind=8), allocatable :: area(:,:)
    
    integer(kind=4) :: nx, ny, nz
    integer(kind=4) :: nx1, nx2, nx3
    
  end type pinto
  type(pinto), allocatable :: zone(:)

end module

module mod_CGNS
  use mod_field
  character(len=stringLen) :: CGNS_filename
  character(len=stringLen) :: CGNS_solnname   

  integer(kind=4) :: index_node, index_grid, index_soln, index_field
  integer(kind=4) :: index_base, index_zone, index_flow, index_coord
  
  integer(kind=4) ier
  
  integer(kind=4), allocatable :: isize(:,:)
  integer(kind=4)              :: isizeDim
  
  integer(kind=4) cell_dim,phys_dim

  character(len=30)  :: basename

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
  integer(kind=4) nmodes

  real(kind=VRL), allocatable :: lambda(:)
  real(kind=VRL), allocatable :: temporal_modes(:,:)
  
  real(kind=8), allocatable , dimension(:,:) :: modos_temporais
  integer(kind=4) :: nsnap_reconst ! Numero de snapshots - reconstrucao

  type pod_zones
    real(kind=VRL), allocatable :: spatial_modes(:,:,:,:,:)
  end type pod_zones
  type(pod_zones), allocatable :: pod_zone(:)
    
  character(len=250) :: path_to_output_matrix

  logical :: Sieber_POD
  logical :: harmonic_POD
  
  real(kind=8) :: dt ! Time step

end module
!#####################################################################################
