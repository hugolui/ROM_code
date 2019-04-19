PROGRAM calc_derivative

	implicit none
	integer(kind=4) :: nt ! Number of snapshots
	integer(kind=4) :: n_modes ! Number of temporal modes
	real(kind=8) :: h ! Time step
	integer(kind=4) :: i,j ! Counters
	real(kind=8), dimension(:,:), allocatable :: X ! Temporal modes
	real(kind=8), dimension(:,:), allocatable :: dX ! Derivative of temporal modes
	character(len=1) :: norm ! Norm
	integer :: der ! Derivative order
	character(len=500) :: path_to_ROM ! Path to working directory
	character(len=500) :: temp

	write(*,*) achar(10) // '  Computing the temporal modes derivatives ...' 

	!!! Read inputs !!!
	open(unit = 60, file = 'compact_scheme_inputs.dat')
	read(60,'(A)') path_to_ROM
  	read(60,*) der
  	close(60)

	!!! Read temporal modes !!!
	open(unit = 12, file = './data/V_matrix.dat')
	read(12,*) norm
	read(12,*) nt
	read(12,*) n_modes
	read(12,*) h
	allocate(X(nt,n_modes)) ! Alocando X
	allocate(dX(nt,n_modes)) ! Alocando dX  
	do i = 1,nt
		read(12,*) X(i,:)
	end do
	close(12)

	!!! Compute the second order derivative !!!
	if (der == 2) then

		do i = 2,nt-1
			dX(i,:) = (X(i+1,:)-X(i-1,:)) / (2.0d0*h)
		end do
		
		dX(1,:) = ((-3.0d0/2.0d0)*X(1,:) + 2*X(2,:) - X(3,:)/2.0d0) / h
		dX(nt,:) = ((3.0d0/2.0d0)*X(nt,:) - 2.0d0*X(nt-1,:) + X(nt-2,:)/2.0d0) / h 

	end if

	!!! Compute the sixth order compact derivative !!!
	if (der == 6) then

		do j = 1,n_modes
			call compact_scheme_6th(nt,h,X(:,j),dX(:,j))
		end do

	end if

	!!! Compute the tenth order compact derivative !!!
	if (der == 10) then
		
		do j = 1, n_modes
			call compact_scheme_10th(nt,h,X(:,j),dX(:,j))
		end do

	end if

	!!! Save the temporal modes derivatives !!!
	open(unit = 13, file = './data/da_'//norm//'.dat') 
	do i = 1,nt
		write(13,'(9999e23.15)') dX(i,:)
	end do
	close(13)

	!!! Move files !!!
	write(temp,"(A)") 'cp -rf data/* ' // trim(path_to_ROM) // 'regression/deep_learning/data/'
  	call system (temp)

	write(*,*) achar(10) // '   Temporal modes derivatives computed ...' 

END PROGRAM calc_derivative

!########################## Sixth order compact scheme ##########################
SUBROUTINE compact_scheme_6th(n,h,vec,diff)

	! This subroutine finds the derivative of a function (vec) using a 6th order compact scheme.
	! It is being considered a non-periodic boundary condition.

	implicit none

	integer, intent(in) :: n
	real(kind=8), intent(in) :: vec(n), diff(n)
	real(kind=8), intent(in) :: h

	integer :: j
	real(kind=8), dimension(n) :: diagA, diagB, diagC
	real(kind=8), dimension(n) :: rhs
	real(kind=8) :: a, b, c
	real(kind=8) :: alpha, beta, alpha1, alpha2


	! Define the matrix coefficients

	! j = 1 and n
	alpha1 = 5.0d0

	! j=2 and n-1
	alpha2 = 2.0d0/11.0d0

	! else
	alpha = 1.0d0/3.0d0
	beta  = 0.0d0
	a = 14.0d0/18.0d0
	b = 1.0d0/36.0d0
	c = 0.0d0

	! Builting tridiagonal vectors (diagonal vectors)

	call diagvec_6th(diagA,diagB,diagC,n,alpha,alpha1,alpha2)

  ! Compute RHS

	! j = 1
	rhs(1) = ((-197.0d0/60.0d0)*vec(1) + (-5.0d0/12.0d0)*vec(2) + 5.0d0*vec(3) + (-5.0d0/3.0d0)*vec(4) + (5.0d0/12.0d0)*vec(5) + (-1.0d0/20.0d0)*vec(6))/h
	j = 2
	rhs(j) = ((-20.0d0/33.0d0)*vec(j-1) + (-35.0d0/132.0d0)*vec(j) + (34.0d0/33.0d0)*vec(j+1) + (-7.0d0/33.0d0)*vec(j+2) + (2.0d0/33.0d0)*vec(j+3) + (-1.0d0/132.0d0)*vec(j+4))/h
 
	! j = n
	rhs(n) = ((197.0d0/60.0d0)*vec(n) + (5.0d0/12.0d0)*vec(n-1) + (-5.0d0)*vec(n-2) + (5.0d0/3.0d0)*vec(n-3) + (-5.0d0/12.0d0)*vec(n-4) + (1.0d0/20.0d0)*vec(n-5))/h
	j=n-1
	rhs(j) = ((20.0d0/33.0d0)*vec(j+1) + (35.0d0/132.0d0)*vec(j) + (-34.0d0/33.0d0)*vec(j-1) + (7.0d0/33.0d0)*vec(j-2) + (-2.0d0/33.0d0)*vec(j-3) + (1.0d0/132.0d0)*vec(j-4))/h  

	do j=3,n-2
		rhs(j) = (vec(j+1) - vec(j-1))*(a/h) + (vec(j+2) - vec(j-2))*(b/h) !+ (vec(j+3) - vec(j-3))*(c/h)
	enddo

	! Solve the linear system to find the derivative
	call tridiagonal(n,diagA,diagB,diagC,rhs,diff) 

	return

END SUBROUTINE compact_scheme_6th
!###################################################################################

!########################## Build tridiagonal matrix ##########################
SUBROUTINE diagvec_6th(A,B,C,n,alpha,alpha1,alpha2)

	implicit none

	integer i ,n
	real(8) A(n), B(n), C(n)
	real(8) alpha, alpha1, alpha2

	B = 1.0d0

	do i=1,n
		A(i) = alpha
		C(i) = alpha
	enddo

	C(1) = alpha1
	C(2) = alpha2

	A(1)  = 0.0d0
	A(2)  = alpha2
	
	A(n) = alpha1
	A(n-1) = alpha2

	C(n)   = 0.0d0
	C(n-1) = alpha2

	return

END SUBROUTINE diagvec_6th
!###################################################################################

!########################## Solve the tridiagonal system ##########################
SUBROUTINE tridiagonal(n,A1,B1,C1,D1,x)

	! This subroutine solve tridiagonal matrix    B[A,D,C] {x} = {b}

	! Ward Cheney and David Kincaid, "Numerical Mathematics and Computing",
	! Thomson Brooks/Cole, 6th ed., chap. 7 - pg 284 (2008).

	! |B(1) C(1)                                        |
	! |A(1) B(2) C(2)                                   |
	! |     A(2) B(3) C(3)                              |
	! |          A(3) B(4) C(4)                         |
	! |                                                 |
	! |                      A(n-3) B(n-2) C(n-2)       |
	! |                             A(n-2) B(n-1) C(n-1)| 
	! |                                    A(n-1) B(N  )|

	Implicit none

	integer, intent(in) :: n
	real(kind=8), intent(in) :: A1(n),B1(n),C1(n),D1(n)

	integer :: i
	real(kind=8) :: A(n),B(n),C(n),D(n),X(n)
		   
  ! Coeficiente c 
  C(1) = C1(1)/B1(1) 
  do i = 2,n-1
    C(i) = C1(i)/(B1(i) - A1(i)*C(i-1))
  end do
      
  ! Coeficiente d 
  D(1) = D1(1)/B1(1)
  do i = 2,n
    D(i) = (D1(i) - A1(i)*D(i-1))/(B1(i) - A1(i)*C(i-1))
  end do

  ! Substituicao regressiva
  x(n) = D(n)
  do i = n-1,1,-1
    x(i) = D(i) - C(i)*x(i+1)
  end do


	return

END SUBROUTINE tridiagonal
!###################################################################################

!########################## Tenth order compact scheme ##########################
SUBROUTINE compact_scheme_10th(n,h,vec,diff)

	! This subroutine finds the derivative of a function (vec) using a 10th order compact scheme.
	! It is being considered a non-periodic boundary condition.

	implicit none

	integer, intent(in) :: n
	real(kind=8), intent(in) :: vec(n), diff(n)
	real(kind=8), intent(in) :: h

	integer :: j
	real(kind=8), dimension(n) :: diagE, diagA, diagD, diagC, diagF
	real(kind=8), dimension(n) :: rhs
	real(kind=8) :: a, b, c, a2, a3, b3, a4, b4, c4
	real(kind=8) :: alpha, beta, alpha1, alpha2, alpha3, beta3, alpha4, beta4


	! Define the matrix coefficients

	! j = 1 and n
	alpha1 = 2.0d0

	! j=2 and n-1
	alpha2 = 0.25d0
	a2 = 3.0d0/4.0d0 

	! j=3 and n-2
	alpha3 = 4.7435d0/10.67175d0
	beta3  = 0.2964375d0/10.67175d0
	a3 = 7.905d0/10.67175d0
	b3 = 1.23515625d0/10.67175d0

	! j=4 and n-3
	alpha4 = 4.63271875d0/9.38146875d0
	beta4  = 0.451390625d0/9.38146875d0
	a4 = 6.66984375d0/9.38146875d0
	b4 = 1.53/9.38146875d0
	c4 = 0.015/9.38146875d0

	! else
	alpha = 0.5d0
	beta  = 0.05d0
	a = 17.0d0/24.0d0
	b = 101.0d0/600.0d0
	c = 0.01d0/6.0d0


	! Builting pentadiagonal vectors (diagonal vectors)

	call diagvec(diagE,diagA,diagD,diagC,diagF,n,beta,beta3,beta4,alpha,alpha1,alpha2,alpha3,alpha4)


	! Compute RHS

	! j = 1
	rhs(1) = (-2.5d0*vec(1) + 2.0d0*vec(2) + 0.5d0*vec(3))/h
	j = 2
	rhs(j) = (vec(j+1) - vec(j-1))*a2/h
	j = 3
	rhs(j) = (vec(j+1) - vec(j-1))*a3/h + (vec(j+2) - vec(j-2))*b3/h
	j = 4
	rhs(j) = (vec(j+1) - vec(j-1))*a4/h + (vec(j+2) - vec(j-2))*b4/h + (vec(j+3) - vec(j-3))*c4/h   

	! j = n
	rhs(n) = (2.5d0*vec(n) - 2.0d0*vec(n-1) - 0.5d0*vec(n-2))/h
	j=n-1
	rhs(j) = (vec(j+1) - vec(j-1))*a2/h  
	j=n-2
	rhs(j) = (vec(j+1) - vec(j-1))*a3/h + (vec(j+2) - vec(j-2))*b3/h
	j=n-3
	rhs(j) = (vec(j+1) - vec(j-1))*a4/h + (vec(j+2) - vec(j-2))*b4/h + (vec(j+3) - vec(j-3))*c4/h

	do j=5,n-4
		rhs(j) = (vec(j+1) - vec(j-1))*a/h + (vec(j+2) - vec(j-2))*b/h + (vec(j+3) - vec(j-3))*c/h
	enddo

	! Solve the linear system to find the derivative

	call PENTA(n,diagE,diagA,diagD,diagC,diagF,rhs,diff) 


	return

END SUBROUTINE compact_scheme_10th
!###################################################################################

!########################## Build the pentadiagonal matrix ##########################
SUBROUTINE diagvec(E1,A2,D3,C4,F5,N,beta,beta3,beta4,alpha,alpha1,alpha2,alpha3,alpha4)

	implicit none

	integer i ,n
	real(8) E1(n), A2(n), D3(n), C4(n), F5(n)
	real(8) alpha, beta, alpha1, alpha2, alpha3, beta3, alpha4, beta4


	d3 = 1.0d0

	do i=1,n
		e1(i) = beta
		a2(i) = alpha
		c4(i) = alpha
		f5(i) = beta
	enddo

	c4(1) = alpha1
	c4(2) = alpha2
	c4(3) = alpha3
	c4(4) = alpha4

	f5(1) = 0.0d0
	f5(2) = 0.0d0
	f5(3) = beta3
	f5(4) = beta4

	a2(1) = alpha2
	a2(2) = alpha3
	a2(3) = alpha4

	e1(1) = beta3
	e1(2) = beta4

	e1(n)   = 0.0d0
	e1(n-1) = 0.0d0
	e1(n-2) = 0.0d0
	e1(n-3) = 0.0d0
	e1(n-4) = beta3
	e1(n-5) = beta4

	a2(n)   = 0.0d0
	a2(n-1) = alpha1
	a2(n-2) = alpha2
	a2(n-3) = alpha3
	a2(n-4) = alpha4

	c4(n)   = 0.0d0
	c4(n-1) = alpha2
	c4(n-2) = alpha3
	c4(n-3) = alpha4

	f5(n)   = 0.0d0
	f5(n-1) = 0.0d0
	f5(n-2) = beta3
	f5(n-3) = beta4


	return

END SUBROUTINE diagvec
!###################################################################################

!########################## Solve the pentadiagonal system ##########################
SUBROUTINE penta(n,E1,A1,D1,C1,F1,b,x)

	! This subroutine solve pentadiagonal matrix    B[E,A,D,C,F] {x} = {b}

	! Ward Cheney and David Kincaid, "Numerical Mathematics and Computing",
	! Thomson Brooks/Cole, 6th ed., chap. 7 - pg 284 (2008).

	! |D(1) C(1) F(1)                                   |
	! |A(1) D(2) C(2) F(2)                              |
	! |E(1) A(2) D(3) C(3) F(3)                         |
	! |     E(2) A(3) D(4) C(4) F(4)                    |
	! |                                                 |
	! |               E(n-4) A(n-3) D(n-2) C(n-2) F(n-2)|
	! |                      E(n-3) A(n-2) D(n-1) C(n-1)| 
	! |                             E(n-1) A(n-1) D(N  )|

	Implicit none

	integer, intent(in) :: n
	real(kind=8), intent(in) :: E1(n),A1(n),D1(n),C1(n),F1(n)

	integer :: i
	real(kind=8) :: E(n),A(n),D(n),C(n),F(n),B(n),X(n), xmult
		    
	E=E1; A=A1; D=D1; C=C1; F=F1

	do i = 2,n-1
		xmult = A(i-1)/D(i-1)
		D(i) = D(i) - xmult*C(i-1)
		C(i) = C(i) - xmult*F(i-1)
		b(i) = b(i) - xmult*b(i-1)
		xmult = E(i-1)/D(i-1)
		A(i) = A(i) - xmult*C(i-1)
		D(i+1) = D(i+1) - xmult*F(i-1)
		B(i+1) = b(i+1) - xmult*b(i-1)
	enddo

	xmult = A(n-1)/D(n-1)
	D(n) = D(n) - xmult*C(n-1)
	x(n) = (b(n) - xmult*b(n-1))/D(n)
	x(n-1) = (b(n-1) - C(n-1)*X(n))/D(n-1)

	do i = n-2,1,-1     
		x(i) = (b(i) - F(i)*x(i+2) - C(i)*x(i+1))/D(i)
	enddo

	return

END SUBROUTINE penta