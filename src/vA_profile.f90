module vA_profile
  use constants
  use mhd_grid

  real(kind=DBL), dimension ( nU1, nU3 ) :: epsilon_, vA, magB, rho
  real(kind=DBL), dimension ( nU1 ) :: vA_eq, magB_eq, n0_eq, lVal
  real(kind=DBL) :: k0, tanHm, tTerm
  
  contains
  subroutine make_vA_profile
		implicit none 
		
		integer :: i, k

    lVal  = 1.0d0 / ( cos ( pi / 2.0d0 - theta0 ) ** 2 )
    k0    = 8d15  ! Earth's magnetic moment
    
    u1Loop: do i = 1, nU1
      
      magB_eq(i)  = sqrt ( 4.0d0 - 3.0d0 * cos ( 0.0d0 ) ** 2 ) / &
        cos ( 0.0d0 ) ** 6 * k0 / ( lVal(i) * rE ) ** 3;
      tanHm       = 2.0d0 * ( lVal(i) - 4.5d0 ) ! Plasmapause location
      tTerm       = 55.0d0 / lVal(i) ** ( 2.6d0 ) * ( tanh ( tanHm ) + 1.0d0 ) / 2.0d0;
      vA_eq(i)    = ( 2.09d0 / lVal(i) ** ( 1.1d0 ) + tTerm ) * 990.0d0 - 50.0d0;
      vA_eq(i)    = vA_eq(i) * 1000.0d0;
      n0_eq(i)    = magB_eq(i) ** 2 / ( u0 * rE * vA_eq(i) ** 2 );
    
      u3Loop: do k = 1, nU3
        
       magB(i,k)  = sqrt ( 4.0d0 - 3.0d0 * cos ( pi / 2.0 - theta(i,1,k) ) ** 2 ) / &
                      cos ( pi / 2.0d0 - theta(i,1,k) ) ** 6 * k0 / ( lVal(i) * rE ) ** 3;
        rho(i,k)  = n0_eq(i) * ( lVal(i) / r(i,1,k) ) ** 3 + &
                      8d7 / rE ** 3 * exp ( -( r(i,1,k) - 1.0d0 ) / ( 600d3 / rE ) );
        vA(i,k)   = magB(i,k) / sqrt ( u0 * rE * rho(i,k) ) / rE;
        
        ! Inverted spanken fast vA profile
        !vA(i,k)	= 1.0 / vA(i,k) * 1e-1;
        
        epsilon_(i,k) = e0 * ( 1.0d0 + c ** 2.0d0 / vA(i,k) ** 2 );
        
      end do u3Loop
    end do u1Loop

		!write (*,*) tTerm
		!write (*,*)
		!write (*,*) tanHm
		!write (*,*)
		!write (*,*) r(:,1,:)
		!write (*,*)
		!write (*,*) rE, u0
		 
    end subroutine make_vA_profile

end module vA_profile

