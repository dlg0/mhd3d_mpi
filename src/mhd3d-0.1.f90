program mhd3d
  use constants
  use mhd_grid
  use vA_profile
	use metric
	use spherHarmFns
	use imsl_libraries
	use dlg_timer
	use dislin
	use lin_sol_svd_int
	use linear_operators

  implicit none

	!	Internal variables

	integer :: nBFns, status

	real, allocatable, dimension ( :,: ) :: YBFnArr_h3, &
		hrBFnArr_h3, hThBFnArr_h3, hPhBFnArr_h3, &
		hrBFnArr_e2, hThBFnArr_e2, hPhBFnArr_e2, &
		hrBFnArr_e1, hThBFnArr_e1, hPhBFnArr_e1, &
		hrBFnArr_h3_gnd, hThBFnArr_h3_gnd, hPhBFnArr_h3_gnd, &
		YBFnArr_h3_gnd, YBFnArr_e1, hrBFnArrT, YBFnArr_e2 ! Basis function arrays

	real, dimension ( nU1, nU2 ) :: sigmaP_N, sigmaH_N, sigmaD_N, &
		sigmaP_S, sigmaH_S, sigmaD_S, alpha_N, sigmaZZ_N, &
		sigma11_N, sigma21_N, sigma12_N, sigma22_N ! Ionosphere conductance

	real :: dTArray ( nU1, nU2, nU3 ), &
		derivU1 ( nU1 ), derivU2 ( nU2 ), derivU3 ( nU3 ), &
		innerSum1, innerSum2, innerSum3, totalSum, dT ! dT
	
	integer :: i, j, k, ii, jj, kk
	integer, dimension ( 0:nU2+1 ) :: jjj ! Azimuthally wrapped index array

	real :: u3_ ( nU3 - 1 ), &
		x_ ( nU1, nU2, nU3 - 1 ), &
		y_ ( nU1, nU2, nU3 - 1 ), &
		z_ ( nU1, nU2, nU3 - 1 ), &
		sqrtG_ ( nU1, nU3 - 1 ), &
		g11cov_ ( nU1, nU3 - 1 ), &
		g12cov_ ( nU1, nU3 - 1 ), &
		g13cov_ ( nU1, nU3 - 1 ), &
		g22cov_ ( nU1, nU3 - 1 ), &
		g21cov_ ( nU1, nU3 - 1 ), &
		g23cov_ ( nU1, nU3 - 1 ), &
		g31cov_ ( nU1, nU3 - 1 ), &
		g32cov_ ( nU1, nU3 - 1 ), &
		g33cov_ ( nU1, nU3 - 1 ), &
		g11con_ ( nU1, nU3 - 1 ), &
		g22con_ ( nU1, nU3 - 1 ), &
		g33con_ ( nU1, nU3 - 1 )	! Interpolated metric components

	!	Metric components on the ionosphere

	real :: g33covIono ( nU1, nU2 ), &
					g11conIono ( nU1, nU2 ), &
					g22conIono ( nU1, nU2 ), &
					g33conIono ( nU1, nU2 )

	real, allocatable :: h3conReal (:,:)

	real, allocatable :: alpha (:,:), beta_(:,:), coeffs (:), coeffsOut(:,:)

	real :: u3MinRemoved ( nU3 ), u3MinRemoved_ ( nU3 - 1 )

	type ( timer ) :: startTime, startTimeOLD

	real :: t, driverFreq, driverAmp, u3Variation, u2Variation, &
		u2DiffF ( nU2 ), u2DiffB ( nU2 ), timeTaken, eta

	real :: h1con ( nU1, nU2, nU3 - 1 ), &
		h2con ( nU1, nU2, nU3 - 1 ), &
		h3con ( nU1, nU2, nU3 ), &
		h1cov ( nU1, nU2, nU3 - 1 ), &
		h2cov ( nU1, nU2, nU3 - 1 ), &
		h3cov ( nU1, nU2, nU3 ), &
		e1con ( nU1, nU2, nU3 ), &
		e2con ( nU1, nU2, nU3 ), &
		e3con ( nU1, nU2, nU3 - 1 ), &
		e1cov ( nU1, nU2, nU3 ), &
		e2cov ( nU1, nU2, nU3 ), &
		e3cov ( nU1, nU2, nU3 - 1 ), &
		j2cov ( nU1, nU2, nU3 )

	!	Fitted variables above and within the ionosphere

	real, dimension ( nU1, nU2 ) :: hrFit_h3, &
		hThFit_e1, hThFit_e2, hPhFit_e1, hPhFit_e2, &
		hThFit_h3, hPhFit_h3, YFit_h3

	!	Variables for interpolation at the ionosphere 

	real :: h1cov_h2 ( nU1, nU2 ), h2cov_h1 ( nU1, nU2 ), &
		t1, t2, t3, h1B1, h1B2, h2B1, h2B2, &
		dhPh_e1 ( nU1, nU2 ), dhTh_e1 ( nU1, nU2 ), &
		dhPh_e2 ( nU1, nU2 ), dhTh_e2 ( nU1, nU2 ), &
		JTh_e1 ( nU1, nU2 ), JPh_e1 ( nU1, nU2 ), &
		JTh_e2 ( nU1, nU2 ), JPh_e2 ( nU1, nU2 ), &
		E1covIono ( nU1, nU2 ), E2covIono ( nU1, nU2 )

	!	Plotting variables
	
	integer, parameter ::  nLevsECov=41
	integer :: nTri, nTri_, nTriN, nTriN_, plotSlice
	integer :: i1Ray (3*nU1*nU3+1), i2Ray (3*nU1*nU3+1), &
		i3Ray (3*nU1*nU3+1), i1Ray_ (3*nU1*(nU3-1)+1), i2Ray_ (3*nU1*(nU3-1)+1), &
		i3Ray_ (3*nU1*(nU3-1)+1)
	integer :: i1RayN (3*nU1*nU3+1), i2RayN (3*nU1*nU3+1), &
		i3RayN (3*nU1*nU3+1), i1RayN_ (3*nU1*(nU3-1)+1), i2RayN_ (3*nU1*(nU3-1)+1), &
		i3RayN_ (3*nU1*(nU3-1)+1)
	real :: plotFreq, eCovLevSpacing
	real :: xPlot (nU1*nU3+3), zPlot (nU1*nU3+3), &
		xPlot_ (nU1*(nU3-1)+3), zPlot_ (nU1*(nU3-1)+3), &
		eCovLevs (nLevsECov), eCovColors (nLevsECov)

  call make_grid 
  call make_vA_profile
	call make_metric

	!	Find the number of basis functions and allocate the arrays
	!	Could avoid this if i found an analytical expression for the
	!	number of basis functions, just too lazy

	nBFns	= numberBFns ()

	allocate ( YBFnArr_h3 ( nBFns, nU1 * nU2 ), &
							hrBFnArr_h3 ( nBFns, nU1 * nU2 ), &
							hThBFnArr_h3 ( nBFns, nU1 * nU2 ), &
							hPhBFnArr_h3 ( nBFns, nU1 * nU2 ), &
							hrBFnArr_e2 ( nBFns, nU1 * nU2 ), &
							hThBFnArr_e2 ( nBFns, nU1 * nU2 ), &
							hPhBFnArr_e2 ( nBFns, nU1 * nU2 ), &
							hrBFnArr_e1 ( nBFns, nU1 * nU2 ), &
							hThBFnArr_e1 ( nBFns, nU1 * nU2 ), &
							hPhBFnArr_e1 ( nBFns, nU1 * nU2 ), &
							YBFnArr_e1 ( nBFns, nU1 * nU2 ), &
							YBFnArr_e2 ( nBFns, nU1 * nU2 ), &
							hrBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
							hThBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
							hPhBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
							YBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
							stat = status )

	!	Generate the SH fit basis function arrays

	call setupSHFns ( nU1 * nU2, nBFns,  rI, &
		reshape ( fitTheta0_h3, (/ nU1 * nU2 /) ), reshape ( fitPhi_h3, (/ nU1 * nU2 /) ), &
		hrBFnArr_h3, hThBFnArr_h3, hPhBFnArr_h3, YBFnArr_h3, &
		minCoLat, maxCoLat + dTheta0 / 2.0, dTheta0, dU2 )

!	write ( *,* ) hThBFnArr_h3 ( 1,: )
!	write ( *,* ) hPhBFnArr_h3 ( 1,: )
!	read ( *,* )
	
	call setupSHFns ( nU1 * nU2, nBFns,  rI, &
		reshape ( fitTheta0_e2, (/ nU1 * nU2 /) ), reshape ( fitPhi_e2, (/ nU1 * nU2 /) ), &
		hrBFnArr_e2, hThBFnArr_e2, hPhBFnArr_e2, YBFnArr_e2, &
		minCoLat, maxCoLat + dTheta0 / 2.0, dTheta0, dU2 )

	call setupSHFns ( nU1 * nU2, nBFns,  rI, &
		reshape ( fitTheta0_e1, (/ nU1 * nU2 /) ), reshape ( fitPhi_e1, (/ nU1 * nU2 /) ), &
		hrBFnArr_e1, hThBFnArr_e1, hPhBFnArr_e1, YBFnArr_e1, &
		minCoLat, maxCoLat + dTheta0 / 2.0, dTheta0, dU2 )

	call setupSHFns ( nU1 * nU2, nBFns,  1.0, &
		reshape ( fitTheta0_h3, (/ nU1 * nU2 /) ), reshape ( fitPhi_h3, (/ nU1 * nU2 /) ), &
		hrBFnArr_h3_gnd, hThBFnArr_h3_gnd, hPhBFnArr_h3_gnd, YBFnArr_h3_gnd, &
		minCoLat, maxCoLat + dTheta0 / 2.0, dTheta0, dU2 )

	sigmaP_N	= 1.0 * rE ** 2!
	sigmaH_N	= 2.5 * rE ** 2!
	sigmaD_N	= 1d6 * rE ** 2!

	sigmaP_S	= 1.0 * rE ** 2!
	sigmaH_S	= 2.5 * rE ** 2!
	sigmaD_S	= 1d6 * rE ** 2!

	alpha_N	=	aCos ( -2.0 * cos ( theta(:,:,1) ) / sqrt ( 1.0 + 3.0 * cos ( theta(:,:,1) ) ** 2 ) )! 
	sigmaZZ_N	= sigmaD_N * cos ( alpha_N ) ** 2 + sigmaP_N * sin ( alpha_N ) ** 2!
	
	sigma11_N	= sigmaD_N * sigmaP_N / sigmaZZ_N!
	sigma21_N	= -sigmaD_N * sigmaH_N * cos ( alpha_N ) / sigmaZZ_N!
	sigma12_N	= sigmaD_N * sigmaH_N * cos ( alpha_N ) / sigmaZZ_N!
	sigma22_N	= sigmaP_N + sigmaH_N ** 2 * sin ( alpha_N ) ** 2 / sigmaZZ_N! 

	!	Calculate the maximum stable dT for this grid

	do i = 1, nU1 
		derivU1(i)	= qdder ( 1, real ( i ), (/ ( real ( j ), j = 1, nU1 ) /), u1 )
	end do
	do j = 1, nU2 
		derivU2(j)	= qdder ( 1, real ( j ), (/ ( real ( i ), i = 1, nU2 ) /), u2 )
	end do
	do k = 1, nU3 
		derivU3(k)	= qdder ( 1, real ( k ), (/ ( real ( j ), j = 1, nU3 ) /), u3 )
	end do

	do i = 1, nU1 
		do j = 1, nU2 
			do k = 1, nU3 
				
				innerSum1 = g11con(i,k) / ( derivU1(i) * derivU1(i) ) &
					+ g12con(i,k) / ( derivU1(i) * derivU2(j) ) &
					+ g13con(i,k)	/ ( derivU1(i) * derivU3(k) )!
				innerSum2 = g21con(i,k) / ( derivU2(j) * derivU1(i) ) &
					+ g22con(i,k) / ( derivU2(j) * derivU2(j) ) &
					+ g23con(i,k)	/ ( derivU2(j) * derivU3(k) )!
				innerSum3 = g31con(i,k) / ( derivU3(k) * derivU1(i) ) &
					+ g32con(i,k) / ( derivU3(k) * derivU2(j) ) &
					+ g33con(i,k)	/ ( derivU3(k) * derivU3(k) )!
				
				totalSum	= innerSum1 + innerSum2 + innerSum3!
				dTArray(i,j,k)	= 1.0 / ( vA(i,k) * sqrt ( totalSum ) )!
	
			end do
		end do
	end do

	dT	= minVal ( dTArray ) * 0.9

	!	Create metric components at interpolated grid points

	do j = 1, nU3 - 1 
		u3_(j)	= qdval ( j + 0.5, (/ ( real ( i ), i = 1, nU3 ) /), u3 )
	end do
	
	do i = 1, nU1 
		do j = 1, nU2
			do k = 1, nU3 - 1

				x_(i,j,k)	= qd3vl ( real(i), real(j), k + 0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(jj),jj=1,nU2) /), (/ (real(kk),kk=1,nU3) /), x )
				y_(i,j,k)	= qd3vl ( real(i), real(j), k + 0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(jj),jj=1,nU2) /), (/ (real(kk),kk=1,nU3) /), y )
				z_(i,j,k)	= qd3vl ( real(i), real(j), k + 0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(jj),jj=1,nU2) /), (/ (real(kk),kk=1,nU3) /), z )

			end do
		end do
	end do

	do i = 1, nU1
		do k = 1, nU3 - 1

			sqrtG_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), sqrtG )
			g11cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g11cov )
			g12cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g12cov )
			g13cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g13cov )
			g22cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g22cov )
			g21cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g21cov )
			g23cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g23cov )
			g31cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g31cov )
			g32cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g32cov )
			g33cov_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g33cov )
			
			g11con_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g11con )
			g22con_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g22con )
			g33con_(i,k)	= qd2vl ( real(i), k+0.5, (/ (real(ii),ii=1,nU1) /), (/ (real(kk),kk=1,nU3) /), g33con )
		
		end do
	end do


	!	Create ionosphere copies of the metric components
	do j = 1, nU2 
			
		g33covIono(:,j)	= g33cov(:,1)
		g11conIono(:,j)	= g11con(:,1)
		g22conIono(:,j)	= g22con(:,1)
		g33conIono(:,j)	= g33con(:,1)

	end do

	allocate ( hrBFnArrT ( nU1 * nU2, nBFns ), &
		alpha ( nBFns, nBFns ), &
		coeffs ( nBFns ), &
		coeffsOut ( nBFns, 1 ), &
		beta_ ( nBFns, 1 ), &
		h3conReal ( nU1, nU2 ) )

	hrBFnArrT	= transpose ( hrBFnArr_h3 )
	alpha	= matMul ( hrBFnArr_h3, hrBFnArrT )

	!write ( *,* ) hrBFnArr_h3

	u3MinRemoved	= u3 - minVal ( u3 )
	u3MinRemoved_	= u3_ - minVal ( u3_ )

	!	Create the azimuthally wrapping index array

	jjj(1:nU2)	= (/ (j,j=1,nU2) /)
	jjj(0)	= nU2
	jjj(nU2+1)	= 1

	! Create the azimuthal coord differenece array, 
	!	i.e., u2Diff

	do j = 1, nU2

		u2DiffF(j)	= ( u2(j) - u2( jjj(j+1) ) )
		if ( j == nU2 ) u2DiffF(j) = u2DiffF(j) - 2.0 * pi

		u2DiffB(j)	= ( u2( jjj(j-1) ) - u2( j ) )
		if ( j == 1 ) u2DiffB(j) = u2DiffB(j) - 2.0 * pi

	end do

	!	BEGIN TIME LOOP

	call metaFl ( 'XWIN' )
	call winSiz ( 1200, 400 )
	call page ( 4800, 1600 )
	call disIni () 
	call winMod ( 'NONE' )
	call setVlt ( 'RAIN' )

	plotFreq	= 1.0!100 * dT

	timeLoop: &
	do t = 0.0, 1000.0, dT
		
		if ( mod ( t, plotFreq ) < dT ) then
		
			startTimeOLD	= startTime
			call start_timer ( startTime )
		
		end if

		!	Apply h boundary conditions
	
		!	Driver
	
		driverFreq	= 50e-3!1./ (200.0 * dT )!
		driverAmp	= 10e-9 / u0!	
	
		if ( t <= ( 0.5 / driverFreq ) ) then
	
			!print, 'DRIVING...'!
			
			do k = 1, nU3
				do j = 1, nU2 
			
					u3Variation	= exp ( -u3(k) ** 2 / 0.000001 )!
					u2Variation	= exp	( -( u2(j) / ( 2 * pi ) - 0.5 ) ** 2 / 0.005 )! 
					h3cov(1,j,k)	= driverAmp * sin ( 2.0 * pi * driverFreq ) &
						* u3Variation * u2Variation / sqrt ( g33con(1,k) )! 
	
				end do
			end do
	
		else
			!print, 'NOT DRIVING...'!
			h3cov(1,:,:)	= 0.0
		end if
	
		h1cov(nU1,:,:)	= 0.0!
		e2cov(nU1,:,:)	= 0.0!

	!	Update contravariant e field

	!	e1con except ionospheres
	
		do k = 2, nU3-1
			do j = 1, nU2
				do i = 1, nU1 
				
				!jj	= j + 1
					!if ( jj > nU2 ) jj = jj - nU2
	
					!if ( k > 1 .AND. k < nU3 ) then
					
						!u2Diff	= ( u2(j) - u2( jjj(j+1) ) )!
						!if ( j == nU2 ) u2Diff = u2Diff - 2.0 * pi
	
						e1con(i,j,k)	= e1con(i,j,k) + &
							dT / ( epsilon_(i,k) * sqrtG(i,k) ) * &
							( ( h3cov(i,j,k) - h3cov(i,jjj(j+1),k) ) / u2DiffF(j) &
							- ( h2cov(i,j,k-1) - h2cov(i,j,k) ) / ( u3MinRemoved_(k-1) - u3MinRemoved_(k) ) )!
					
					!end if

				end do
			end do
		end do

		!	e2con except inner boundary and ionospheres

		do k = 2, nU3-1
			do j = 1, nU2
				do i = 1, nU1-1 

					!if ( k > 1 .AND. i < nU1 .AND. k < nU3 ) & 
					e2con(i,j,k)	= e2con(i,j,k) + &
						dT / ( epsilon_(i,k) * sqrtG(i,k) ) * &
						( ( h1cov(i,j,k-1) - h1cov(i,j,k) ) / ( u3MinRemoved_(k-1) - u3MinRemoved_(k) ) &
						- ( h3cov(i,j,k) - h3cov(i+1,j,k) ) / ( u1(i) - u1(i+1) ) )!
	
				end do
			end do
		end do

	!	Calculate covariant e from contravariant e

	!	Do not update the points at the ionosphere as they
	!	have already been filled by the BC
	
	!	e1cov except ionospheres
	
		do k = 2, nU3-1
			do j = 1, nU2
				do i = 1, nU1 

					!if ( k > 1 .AND. k < nU3 ) &
					e1cov(i,j,k)	= 1.0 / g11con(i,k) * e1con(i,j,k)! &
					
					!	Do not update the points at the ionosphere as they
					!	have already been filled by the BC
			
				end do
			end do
		end do
	
		!	e2cov except inner boundary and ionospheres

		do k = 2, nU3-1
			do j = 1, nU2
				do i = 1, nU1-1

					!if ( k > 1 .AND. k < nU3 .AND. i < nU1 ) &
						e2cov(i,j,k)	= g22cov(i,k) * e2con(i,j,k)!
	
				end do
			end do
		end do
		
		!	Update contravariant h field
	
		!	h1con except inner boundary

		do k = 1, nU3-1
			do j = 1, nU2
				do i = 1, nU1-1

					!jj	= j - 1
					!if ( jj < 1 ) jj = jj + nU2

					!if ( k < nU3 .AND. i < nU1 ) then
					
						!u2Diff	= ( u2(jj) - u2(j) )!
						!if ( j == 1 ) u2Diff = u2Diff - 2.0 * pi
	
						h1con(i,j,k)	= h1con(i,j,k) - &
							dT / ( u0 * sqrtG_(i,k) ) * &
							(	-( e2cov(i,j,k) - e2cov(i,j,k+1) ) / ( u3MinRemoved(k) - u3MinRemoved(k+1) ) )!
	
					!end if

				end do
			end do
		end do

		!	h2con all locations

		do k = 1, nU3-1
			do j = 1, nU2
				do i = 1, nU1
		
					!if ( k < nU3 ) &
					h2con(i,j,k)	= h2con(i,j,k) - &
						dT / ( u0 * sqrtG_(i,k) ) * &
						( ( e1cov(i,j,k) - e1cov(i,j,k+1) ) / ( u3MinRemoved(k) - u3MinRemoved(k+1) ) )!&
				
					end do
			end do
		end do
	
		!	h3con except outer boundary

		do k = 1, nU3
			do j = 1, nU2
				do i = 2, nU1
	
					!if ( i > 1 ) then

						!u2Diff	= ( u2(jj) - u2(j) )!
						!if ( j == 1 ) u2Diff = u2Diff - 2.0 * pi
	
						h3con(i,j,k)	= h3con(i,j,k) - &
							dT / ( u0 * sqrtG(i,k) ) * &
							( ( e2cov(i-1,j,k) - e2cov(i,j,k) ) / ( u1(i-1) - u1(i) ) &
							- ( e1cov(i,jjj(j-1),k) - e1cov(i,j,k) ) / u2DiffB(j)	)!
	
					!endif
				
				end do
			end do
		end do

		!	Set the dh3/du1=0 @ outer boundary in the two ionospheres
	
		h3con(1,:,1)	= ( 4.0 * h3con(2,:,1) - h3con(3,:,1) ) / 3.0!
		h3con(1,:,nU3)	= ( 4.0 * h3con(2,:,nU3) - h3con(3,:,nU3) ) / 3.0! 
	
		do k = 1, nU3
			do j = 1, nU2
				do i = 1, nU1

					!	Calculate covariant h from contravariant h
			
					!	h1cov except inner boundary
	
					if ( i < nU1 .AND. k < nU3 ) &
					h1cov(i,j,k)	= g11cov_(i,k) * h1con(i,j,k) &
						+ 0.25 * g13cov(i,k) * &
							( h3con(i,j,k) + h3con(i,j,k+1) + h3con(i+1,j,k) + h3con(i+1,j,k+1) )!
					if ( i == nU1 .AND. k < nU3 ) &
						h1cov(i,j,k)	= 0.0!
	
					!	h2cov everywhere
	
					if ( k < nU3 ) &
					h2cov(i,j,k)	= g22cov_(i,k) * h2con(i,j,k)! 
	
					!	h3cov except outer boundary and ionospheres
	
					if ( i > 1 .AND. k > 1 .AND. k < nU3 ) &
					h3cov(i,j,k)	= g33cov(i,k) * h3con(i,j,k) &
						+ 0.25 * g31cov(i,k) * &
							( h1con(i-1,j,k) + h1con(i-1,j,k-1) + h1con(i,j,k) + h1con(i,j,k-1) )! 
					
					!	h3cov at ionospheres except outer boundary
					
					if ( i > 1 .AND. k == nU3 ) &
					h3cov(i,j,k)	= g33cov(i,k) * h3con(i,j,k) &
						+ 0.5 * g31cov(i,k) * &
							( h1con(i-1,j,k-1) + h1con(i,j,k-1) )! &
					if ( i > 1 .AND. k == 1 ) &
					h3cov(i,j,k)	= g33cov(i,k) * h3con(i,j,k) &
						+ 0.5 * g31cov(i,k) * &
							( h1con(i-1,j,k) + h1con(i,j,k) )! &
				
				end do
			end do
		end do
		
		!	Reconstruct the horizontal b field component below 
		!	the ionosphere using a Spherical Harmonic fit.

	!h3conReal	= h3con(:,:,1) * sqrt ( g33covIono )

	beta_	= matMul ( reshape ( h3con(:,:,1) * sqrt ( g33covIono ), (/ nU1 * nU2, 1 /) ), hrBFnArrT )
	
	!write ( *,* ) h3con(:,:,1)
	call lin_sol_svd ( alpha, beta_, coeffsOut )

	coeffs	= coeffsOut(:,1)
	
	!write ( *,* ) h3con(:,:,1)
	!write ( *,* ) coeffs

	!hrFit_h3	= reshape ( matMul ( coeffs, hrBfnArr_h3 ), (/ nU1, nU2 /) )
	hThFit_e1	= -reshape ( matMul ( coeffs, hThBFnArr_e1 ), (/ nU1, nU2 /) ) 
	hThFit_e2	= -reshape ( matMul ( coeffs, hThBFnArr_e2 ), (/ nU1, nU2 /) )
	hPhFit_e1	= reshape ( matMul ( coeffs, hPhBFnArr_e1 ), (/ nU1, nU2 /) ) 
	hPhFit_e2	= reshape ( matMul ( coeffs, hPhBFnArr_e2 ), (/ nU1, nU2 /) ) 

	!write ( *,* ) alpha
	!write ( *,* ) hrFit_h3

	!read ( *,* )
	
	!hThFit_h3	= -reshape ( hThBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) );
	!hPhFit_h3	= reshape ( hPhBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) ); 
	!YFit_h3	= reshape ( YBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) ); 

!	Calculate the eField boundary values from reconstructed b 
!	and assumed sigma values. First need reconstructed h2 and h1
!	at e1 and e1 locations. Here use a simple averaging ... but 
!	could use the SH reconstruction to get them at those points too.

	!	JUST ABOVE THE IONOSPHERE
	!	Get h1 @ h2 locations


	do i = 1, nU1
		do j = 1, nU2

			jj	= j + 1
			if ( jj > nU2 ) jj = jj - nU2

			if ( i > 1 ) then 

				h1cov_h2(i,j)	= ( h1cov(i,j,1) &
					+ h1cov(i,jj,1) &
					+ h1cov(i-1,j,1) &
					+ h1cov(i-1,jj,1) ) / 4.0!
			
			else 
			
			!	@ i = 0 need to extrapolate h1 to h2 locations.
			!	Use Taylor expansion to get h1 above and below
			!	h2 location, then linear interpolate :)

				t1	= h1cov(1,j,1)!
				t2	= dTheta0 / 2.0 * ( -3.0 * h1cov(1,j,1) + 4.0 * h1cov(2,j,1) &
					- h1cov(3,j,1) ) / ( 2.0 * dTheta0 )!
				t3	= ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h1cov(1,j,1) &
					- 5.0 * h1cov(2,j,1) + 4.0 * h1cov(3,j,1) - h1cov(4,j,1) ) / dTheta0 ** 2!

				h1B1	= t1 + t2 + t3!
	
				t1	= h1cov(1,jj,1)!
				t2	= dTheta0 / 2.0 * ( -3.0 * h1cov(1,jj,1) &
					+ 4.0 * h1cov(2,jj,1) &
					- h1cov(3,jj,1) ) / ( 2.0 * dTheta0 )!
				t3	= ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h1cov(1,jj,1) &
					- 5.0 * h1cov(2,jj,1) &
					+ 4.0 * h1cov(3,jj,1) &
					- h1cov(4,jj,1) ) / dTheta0 ** 2!

				h1B2	= t1 + t2 + t3!
			
				h1cov_h2(1,j)	= ( h1B1 + h1B2 ) / 2.0!
			
			end if

		end do
	end do
	
	!	Get h2 @ h1 locations above the ionosphere

	do i = 1, nU1
		do j = 1, nU2

			jj	= j - 1
			if ( jj < 1 ) jj = jj + nU2

			if ( i < nU1 ) then

				h2cov_h1(i,j)	= ( h2cov(i,j,1) &
					+ h2cov(i,jj,1) &
					+ h2cov(i+1,j,1) &
					+ h2cov(i+1,jj,1) ) / 4.0!
			
			else
			
			!	@ i = nU1 - 1 need to extrapolate h2 to h1 locations.
			!	Use Taylor expansion to get h2 above and below
			!	h1 location, then linear interpolate :)

				t1	= h2cov(i,j,1)!
				t2	= dTheta0 / 2.0 * ( 3.0 * h2cov(i,j,1) - 4.0 * h2cov(i-1,j,1) &
					+ h2cov(i-2,j,1) ) / ( 2.0 * dTheta0 )!
				t3	= ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h2cov(i,j,1) &
					- 5.0 * h2cov(i-1,j,1) + 4.0 * h2cov(i-2,j,1) - h2cov(i-3,j,1) ) / dTheta0 ** 2!

				h2B1	= t1 + t2 + t3!
	
				t1	= h2cov(i,jj,1)!
				t2	= dTheta0 / 2.0 * ( -3.0 * h2cov(i,jj,1) &
					+ 4.0 * h2cov(i-1,jj,1) &
					- h2cov(i-2,jj,1) ) / ( 2.0 * dTheta0 )!
				t3	= ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h2cov(i,jj,1) &
					- 5.0 * h2cov(i-1,jj,1) &
					+ 4.0 * h2cov(i-2,jj,1) &
					- h2cov(i-3,jj,1) ) / dTheta0 ** 2!

				h2B2	= t1 + t2 + t3!
			
				h2cov_h1(i,j)	= ( h2B1 + h2B2 ) / 2.0!
			
			end if

		end do
	end do
	
	dhPh_e1	= (h2cov(:,:,1) * sqrt ( g22conIono ) - hPhFit_e1)!
	dhTh_e1	= (h1cov_h2 * sqrt ( g11conIono ) - hThFit_e1)!

	dhPh_e2	= (h2cov_h1 * sqrt ( g22conIono ) - hPhFit_e2)!
	dhTh_e2	= (h1cov(:,:,1) * sqrt ( g11conIono ) - hThFit_e2)!

	JTh_e1	= -dhPh_e1!
	JPh_e1	= dhTh_e1!

	JTh_e2	= -dhPh_e2!
	JPh_e2	= dhTh_e2!

	E1covIono	= ( dhTh_e1 * sigma21_N + dhPh_e1 * sigma22_N ) / &
				( sigma12_N * sigma21_N - sigma11_N * sigma22_N )!

	E2covIono	= ( dhTh_e2 * sigma11_N + dhPh_e2 * sigma12_N ) / &
				( -sigma12_N * sigma21_N + sigma11_N * sigma22_N )!

	!write ( *,* ) E1covIono
	!write ( *,* ) E2covIono

!	NEED TO SET SIGMAS SLIGHTLY DIFFERENT do GRID OFFSET!

	e1cov(:,:,1)	= E1covIono / sqrt ( g11conIono )! 
	e2cov(:,:,1)	= E2covIono / sqrt ( g22conIono )!
	
		
		!	Plot legendre functions using dislin library

		if ( mod ( t, plotFreq ) < dT ) then  

			plotSlice = nU2 / 2 + 1
			
			xPlot(1:nU1*nU3)	= reshape ( sqrt ( x(:,plotSlice,:)**2 + y(:,plotSlice,:)**2 ), (/ nU1 * nU3 /) )
			zPlot(1:nU1*nU3)	= reshape ( z ( :, plotSlice, : ), (/ nU1 * nU3 /) )
		
			call triang ( xPlot, zPlot, &
					nU1 * nU3, i1Ray, i2Ray, i3Ray, 3 * nU1 * nU3 + 1, nTri )
			call triang ( -xPlot, zPlot, &
					nU1 * nU3, i1RayN, i2RayN, i3RayN, 3 * nU1 * nU3 + 1, nTriN )
			
			xPlot_(1:nU1*(nU3-1))	= reshape ( sqrt ( x_(:,plotSlice,:)**2 + y_(:,plotSlice,:)**2 ), (/ nU1 * (nU3-1) /) )
			zPlot_(1:nU1*(nU3-1))	= reshape ( z_( :, plotSlice, : ), (/ nU1 * (nU3-1) /) )
	
			call triang ( xPlot_, zPlot_, &
					nU1 * (nU3-1), i1Ray_, i2Ray_, i3Ray_, 3 * nU1 * (nU3-1) + 1, nTri_ )
			call triang ( -xPlot_, zPlot_, &
					nU1 * (nU3-1), i1RayN_, i2RayN_, i3RayN_, 3 * nU1 * (nU3-1) + 1, nTriN_ )
					
			call erase ()
			call height ( 18 )
			call pageRa ()
			call pagFll ( 255 ) ! Set page background color
			call frame ( 0 ) ! Set plot frame thickness to 0
			call setGrf ( 'NONE', 'NONE', 'NONE', 'NONE' )

			eCovLevSpacing	= 0.005	! [mVm^-1]
			eCovLevs	= ( (/ ( real(i), i=0,nLevsECov-1 ) /) - real ( int ( nLevsECov / 2 ) ) ) &
				* eCovLevSpacing
			eCovColors	= ( eCovLevs + abs ( minVal ( eCovLevs ) ) ) &
				/ maxVal ( eCovLevs + abs ( minVal ( eCovLevs ) ) ) &
				* 253 + 1.0
		
			call axsPos ( 50, 800 )
			call axsLen ( 1500, 800 )
			call graf ( -9.0, 9.0, -9.0, 1.0, -4.0, 4.0, -4.0, 1.0 )
			call conFll ( xPlot(1:nU1*nU3), zPlot(1:nU1*nU3), &
				reshape ( e1cov (:,plotSlice,:) * sqrt ( g11con ) * rE * 1e3, (/ nU1 * nU3 /) ), &
				nU1 * nU3, i1Ray, i2Ray, i3Ray, nTri, eCovLevs, nLevsECov )
			call conFll ( -xPlot(1:nU1*nU3), zPlot(1:nU1*nU3), &
				reshape ( e1cov (:,1,:) * sqrt ( g11con ) * rE * 1e3, (/ nU1 * nU3 /) ), &
				nU1 * nU3, i1RayN, i2RayN, i3RayN, nTriN, eCovLevs, nLevsECov )
			call endGrf ()

			call axsPos ( 1650, 800 )
			call axsLen ( 1500, 800 )
			call graf ( -9.0, 9.0, -9.0, 1.0, -4.0, 4.0, -4.0, 1.0 )
			call conFll ( xPlot(1:nU1*nU3), zPlot(1:nU1*nU3), &
				reshape ( e2cov (:,plotSlice,:) * sqrt ( g22con ) * rE * 1e3, (/ nU1 * nU3 /) ), &
				nU1 * nU3, i1Ray, i2Ray, i3Ray, nTri, eCovLevs, nLevsECov )
			call conFll ( -xPlot(1:nU1*nU3), zPlot(1:nU1*nU3), &
				reshape ( e2cov (:,1,:) * sqrt ( g22con ) * rE * 1e3, (/ nU1 * nU3 /) ), &
				nU1 * nU3, i1RayN, i2RayN, i3RayN, nTriN, eCovLevs, nLevsECov )
			call endGrf ()

			call axsPos ( 3250, 800 )
			call axsLen ( 1500, 800 )
			call graf ( -9.0, 9.0, -9.0, 1.0, -4.0, 4.0, -4.0, 1.0 )
			call conFll ( xPlot_(1:nU1*(nU3-1)), zPlot_(1:nU1*(nU3-1)), &
				reshape ( e3cov (:,plotSlice,:) * sqrt ( g33con_ ) * rE * 1e3, (/ nU1 * (nU3-1) /) ), &
				nU1 * (nU3-1), i1Ray_, i2Ray_, i3Ray_, nTri_, eCovLevs, nLevsECov )
			call conFll ( -xPlot_(1:nU1*(nU3-1)), zPlot_(1:nU1*(nU3-1)), &
				reshape ( e3cov (:,1,:) * sqrt ( g33con_ ) * rE * 1e3, (/ nU1 * (nU3-1) /) ), &
				nU1 * (nU3-1), i1RayN_, i2RayN_, i3RayN_, nTriN_, eCovLevs, nLevsECov )
			call endGrf ()

			timeTaken	= end_timer ( startTimeOLD )
			eta	= timeTaken / plotFreq * ( 1000.0 - t ) / 3600.0

			if ( eta < 1 ) then 
				write (*,'(a12,f5.1,a6,f4.1,a8)') 'Model time: ', t, &
					' ETA: ', eta * 60.0, &
					' [MINS]'
			else
				write (*,'(a12,f5.1,a6,f4.1,a8)') 'Model time: ', t, &
					' ETA: ', eta, &
					' [HOURS]'
			end if

		end if

	end do timeLoop

	call disFin () 

  read (*,*)

end program mhd3d
