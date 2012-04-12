module dlg_routines

contains
function dlg_deriv1 ( xPt, xarray, yarray )

	implicit none
	real :: yarray(:), xarray(:), xPt
	real :: dydx(size(yarray))
	integer :: j, nPts, ii1, ii2

	nPts	= size ( yarray )

	!	Calculate derivative

	do j = 1, size ( yarray )

		if ( j > 1 .AND. j < nPts ) dydx(j) = &
			1.0 / (  xarray(j+1) - xarray(j-1) ) * ( yarray(j+1) - yarray(j-1) ) 
		if ( j == 1 ) dydx(j) = 1.0 / ( xarray(j+2) - xarray(j) ) &
			* ( -3.0 * yarray(j) + 4.0 * yarray(j+1) - yarray(j+2) )
		if ( j == nPts ) dydx(j) = 1.0 / ( xarray(j) - xarray(j-2) ) &
			* ( 3.0 * yarray(j) - 4.0 * yarray(j-1) + yarray(j-2) )

	end do

	!	Linear interpolate to desired point

end function dlg_deriv1

subroutine nr_spline ( x, y, yp1, ypn, y2 )
	implicit none
	real, dimension(:), intent(IN) :: x, y
	real, intent(IN) :: yp1, ypn
	real, dimension(:), intent(OUT) :: y2
	integer :: n
	real, dimension(size(x)) :: a, b, c, r
	c(1:n-1)	= x(2:n) - x(1:n-1)
	r(1:n-1)	= 6.0 * ( ( y(2:n) - y(1:n-1) ) / c(1:n-1 ) )
	r(2:n-1)	= r(2:n-1) - r(1:n-2)
	a(2:n-1)	= c(1:n-2)
	b(2:n-1)	= 2.0 * ( c(2:n-1) + a(2:n-1) )
	b(1)	= 1.0
	b(n)	= 1.0
	if ( yp1 > 0.99e30 ) then
		r(1)	= 0.0
		c(1)	= 0.0
	else
		r(1)	= ( 3.0 / ( x(2) - x(1) ) ) * &
			( ( y(2) - y(1) ) / ( x(2) - x(1) ) - yp1 )
		c(1)	= 0.5
	end if
	if ( ypn > 0.99e30 ) then
		r(n)	= 0.0
		a(n)	= 0.0
	else
		r(n)	= (-3.0 / ( x(n) - x(n-1) ) ) * &
			( ( y(n) - y(n-1) ) / ( x(n) - x(n-1) ) - ypn )
	end if
	call nr_tridiag ( a(2:n), b(1:n), c(1:n-1), r(1:n), y2(1:n) )
end subroutine nr_spline

function nr_splint ( xa, ya, y2a, x )
implicit none
real, dimension(:), intent(IN) :: xq, yq, y2a
real, intent(IN) :: x
real :: splint
integer :: khi, klo, n
real :: a, b, h
klo	= max ( min ( locate ( xa, x ), n-1 ), 1 )

end module dlg_routines

