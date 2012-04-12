program mhd3d
    !use mpi
    use constants
    use mhd_grid
    use vA_profile
    use metric
    use dlg_timer
    implicit none
    include "mpif.h" 

    !   Internal variables

    integer, parameter :: SGL_  = selected_real_kind ( p = 6, r = 37 )
    integer, parameter :: DBL_  = selected_real_kind ( p = 13, r = 200 )

    integer :: nBFns, status

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
        g33con_ ( nU1, nU3 - 1 )    ! Interpolated metric components

    real, allocatable :: h3conReal (:,:)

    real :: u3MinRemoved ( nU3 ), u3MinRemoved_ ( nU3 - 1 )

    type ( timer ) :: startTime, startTimeOLD

    real :: driverFreq, driverAmp, u3Variation, u2Variation, &
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

    real :: h, tRun
    integer :: t, nt

    !   MPI variables

    integer :: mpi_iErr, mpi_nP, mpi_pId
    integer :: iStart, iStop, jStart, jStop, kStart, kStop
    integer :: iStart_, iStop_, jStart_, jStop_, kStart_, kStop_
    integer :: mpi_nDims 
    integer :: mpi_dimSize(3)
    integer :: mpi_coords(3)
    logical :: mpi_dimPeriodicity(3)
    integer :: mpi_comm_cart, mpi_type_vec_u1, mpi_type_vec_u2
    integer :: mpi_nU1, mpi_nU2, mpi_nU3
    integer :: mpi_type_2D_U1U3, mpi_type_2D_U1U2, mpi_type_2D_U2U3
    logical :: mpi_reOrder
    integer :: mpi_comm_old
    !   mpi_nP_U1 * mpi_nP_U2 * mpi_nP_U3 == mpi_nP or program will crash
    !   with a "MPI_ERR_ARG: invalid argument of some other kind" error
    !   Also, the grid must be integer multiples of these values in each 
    !   appropriate direction.
    integer :: mpi_nP_U1 = 4, mpi_nP_U2 = 1, mpi_nP_U3 = 1
    integer :: mpi_coords_dest(3), mpi_coords_src(3)
    integer :: mpi_pId_dest, mpi_pId_src
    integer :: mpi_tag_e2cov
    integer, dimension ( MPI_STATUS_SIZE ) :: mpi_status

    real :: progFreq

    !   START PROGRAM

    call mpi_init ( mpi_iErr )
    call mpi_comm_size ( MPI_COMM_WORLD, mpi_nP, mpi_iErr ) 
    call mpi_comm_rank ( MPI_COMM_WORLD, mpi_pId, mpi_iErr )

    mpi_nU1 = nU1 / mpi_nP_U1
    mpi_nU2 = nU2 / mpi_nP_U2
    mpi_nU3 = nU3 / mpi_nP_U3
    mpi_nDims   = 3
    mpi_dimSize(1)  = mpi_nP_U1
    mpi_dimSize(2)  = mpi_nP_U2
    mpi_dimSize(3)  = mpi_nP_U3 
    mpi_dimPeriodicity(1)   = .false.
    mpi_dimPeriodicity(2)   = .true.
    mpi_dimPeriodicity(3)   = .false.
    mpi_reOrder = .true.
    mpi_comm_old    = MPI_COMM_WORLD

    !   Create the mpi topology (cartesian grid in U1,U2,U3)
    !   The cartesian topology in domain decomposition is useful
    !   because using calls to mpi_cart_coords and mpi_cart_rank
    !   we can figure out the rank (process id) of each domain
    !   piece referenced by its coordinates.

    call mpi_cart_create ( mpi_comm_old, mpi_nDims, mpi_dimSize, &
        mpi_dimPeriodicity, mpi_reOrder, mpi_comm_cart, mpi_iErr ) 

    !call mpi_cart_coords ( mpi_comm_cart, mpi_pId, 3, mpi_coords, mpi_iErr )
    !call mpi_cart_rank ( mpi_comm_cart, mpi_coords, mpi_rank1, mpi_iErr )
    !
    !write(*,*) 'coords: ', mpi_pId, mpi_coords, mpi_rank1

    !   Create mpi derived data types for passing boundary plane
    !   data between processes

    !   U1U3 plane
    
    call mpi_type_contiguous ( mpi_nU1, MPI_REAL, mpi_type_vec_u1, &
        mpi_iErr )
    call mpi_type_commit ( mpi_type_vec_u1, mpi_iErr )
    call mpi_type_hVector ( mpi_nU3, 1, ( mpi_nU1 + 1 ) * ( mpi_nU2 + 1 ), &
        mpi_type_vec_u1, mpi_type_2D_U1U3, mpi_iErr ) 
    call mpi_type_commit ( mpi_type_2D_U1U3, mpi_iErr )

    !   U1U2 plane

    call mpi_type_hVector ( mpi_nU2, 1, ( mpi_nU1 + 1 ), &
        mpi_type_vec_u1, mpi_type_2D_U1U2, mpi_iErr )
    call mpi_type_commit ( mpi_type_2D_U1U2, mpi_iErr )

    !   U2U3 plane

    call mpi_type_vector ( mpi_nU2, 1, mpi_nU1, MPI_REAL, &
        mpi_type_vec_u2, mpi_iErr )
    call mpi_type_commit ( mpi_type_vec_u2, mpi_iErr )
    call mpi_type_hVector ( mpi_nU3, 1, mpi_nU1*mpi_nU2, &
        mpi_type_vec_u2, mpi_type_2D_U2U3, mpi_iErr )
    call mpi_type_commit ( mpi_type_2D_U2U3, mpi_iErr )

    write (*,*) 'TaskId: ', mpi_pId, ' of ', mpi_nP
    !read (*,*)

    call make_grid 
    call make_vA_profile
    call make_metric

    !   Calculate the maximum stable dT for this grid

    write (*,*) 'Calculating derivatives of coordinate variables ...'   

    do i = 1, nU1 
        h   = 1.0
        if ( i > 1 .AND. i < nU1 ) derivU1(i) = 1.0 / ( 2.0 * h ) * ( u1(i+1) - u1(i-1) ) 
        if ( i == 1 ) derivU1(i) = 1.0 / ( 2.0 * h ) &
            * ( -3.0 * u1(i) + 4.0 * u1(i+1) - u1(i+2) )
        if ( i == nU1 ) derivU1(i) = 1.0 / ( 2.0 * h ) &
            * ( 3.0 * u1(i) - 4.0 * u1(i-1) + u1(i-2) )
        !derivU1(i) = qdder ( 1, real ( i, kind=DBL_ ), (/ ( real ( j, kind=DBL_ ), j = 1, nU1 ) /), u1 )
    end do

    do j = 1, nU2 
        h   = 1.0
        if ( j > 1 .AND. j < nu2 ) derivU2(j) = 1.0 / ( 2.0 * h ) * ( u2(j+1) - u2(j-1) ) 
        if ( j == 1 ) derivU2(j) = 1.0 / ( 2.0 * h ) &
            * ( -3.0 * u2(j) + 4.0 * u2(j+1) - u2(j+2) )
        if ( j == nU2 ) derivU2(j) = 1.0 / ( 2.0 * h ) &
            * ( 3.0 * u2(j) - 4.0 * u2(j-1) + u2(j-2) )
        !derivU2(j) = qdder ( 1, real ( j, kind=DBL_ ), (/ ( real ( i, kind=DBL_ ), i = 1, nU2 ) /), u2 )
    end do

    do k = 1, nU3 
        h   = 1.0   
        if ( k > 1 .AND. k < nu3 ) derivU3(k) = 1.0 / ( 2.0 * h ) * ( u3(k+1) - u3(k-1) ) 
        if ( k == 1 ) derivU3(k) = 1.0 / ( 2.0 * h ) &
            * ( -3.0 * u3(k) + 4.0 * u3(k+1) - u3(k+2) )
        if ( k == nU3 ) derivU3(k) = 1.0 / ( 2.0 * h ) &
            * ( 3.0 * u3(k) - 4.0 * u3(k-1) + u3(k-2) )
        !derivU3(k) = qdder ( 1, real ( k, kind=DBL_ ), (/ ( real ( j, kind=DBL_ ), j = 1, nU3 ) /), u3 )
    end do

    write (*,*) 'Calculating minimum dT value ... '

    do i = 1, nU1 
        do j = 1, nU2 
            do k = 1, nU3 
                
                innerSum1 = g11con(i,k) / ( derivU1(i) * derivU1(i) ) &
                    + g12con(i,k) / ( derivU1(i) * derivU2(j) ) &
                    + g13con(i,k)   / ( derivU1(i) * derivU3(k) )!
                innerSum2 = g21con(i,k) / ( derivU2(j) * derivU1(i) ) &
                    + g22con(i,k) / ( derivU2(j) * derivU2(j) ) &
                    + g23con(i,k)   / ( derivU2(j) * derivU3(k) )!
                innerSum3 = g31con(i,k) / ( derivU3(k) * derivU1(i) ) &
                    + g32con(i,k) / ( derivU3(k) * derivU2(j) ) &
                    + g33con(i,k)   / ( derivU3(k) * derivU3(k) )!
                
                totalSum    = innerSum1 + innerSum2 + innerSum3!
                dTArray(i,j,k)  = 1.0 / ( vA(i,k) * sqrt ( totalSum ) )!
     
            end do
        end do
    end do

    dT  = minVal ( dTArray ) * 0.9
    tRun    = 1000.0 
    nt  = tRun / dT    

    write (*,*) 'dT: ', dT

    !   Create metric components at interpolated grid points

    write (*,*) 'Calculating metric components at interpolated grid points ... '

    do j = 1, nU3 - 1 
        !write (*,*) qdval ( j + 0.5d0, (/ ( real ( i, kind=DBL_ ), i = 1, nU3 ) /), u3 )
        ! DLG: should probably spline this, and 
        ! those below, but this was quickest to remove IMSL dependence ;-)
        u3_(j)  = ( u3(j) + u3(j+1) ) / 2.0
        !write (*,*) u3_(j)
    end do
    
    do i = 1, nU1 
        do j = 1, nU2
            do k = 1, nU3 - 1

                x_(i,j,k)   = ( x(i,j,k)+x(i,j,k+1) ) / 2.0
                y_(i,j,k)   = ( y(i,j,k)+y(i,j,k+1) ) / 2.0
                z_(i,j,k)   = ( z(i,j,k)+z(i,j,k+1) ) / 2.0

            end do
        end do
    end do
    
    do i = 1, nU1
        do k = 1, nU3 - 1

            sqrtG_(i,k) = ( sqrtG(i,k) + sqrtG(i,k+1) ) / 2.0
            g11cov_(i,k)    = ( g11cov(i,k) + g11cov(i,k+1) ) / 2.0
            g12cov_(i,k)    = ( g12cov(i,k) + g12cov(i,k+1) ) / 2.0
            g13cov_(i,k)    = ( g13cov(i,k) + g13cov(i,k+1) ) / 2.0
            g22cov_(i,k)    = ( g22cov(i,k) + g22cov(i,k+1) ) / 2.0
            g21cov_(i,k)    = ( g21cov(i,k) + g21cov(i,k+1) ) / 2.0
            g31cov_(i,k)    = ( g31cov(i,k) + g31cov(i,k+1) ) / 2.0
            g32cov_(i,k)    = ( g32cov(i,k) + g32cov(i,k+1) ) / 2.0
            g33cov_(i,k)    = ( g33cov(i,k) + g33cov(i,k+1) ) / 2.0

            g11con_(i,k)    = ( g11con(i,k) + g11con(i,k+1) ) / 2.0
            g22con_(i,k)    = ( g22con(i,k) + g22con(i,k+1) ) / 2.0
            g33con_(i,k)    = ( g33con(i,k) + g33con(i,k+1) ) / 2.0
    
        end do
    end do

    !   Create ionosphere copies of the metric components

!   write (*,*) 'Creating ionosphere copies of the metric components ...'
!   
!   do j = 1, nU2 
!           
!       g33covIono(:,j) = g33cov(:,1)
!       g11conIono(:,j) = g11con(:,1)
!       g22conIono(:,j) = g22con(:,1)
!       g33conIono(:,j) = g33con(:,1)
!
!       g33covIonoS(:,j)    = g33cov(:,nU3)
!       g11conIonoS(:,j)    = g11con(:,nU3)
!       g22conIonoS(:,j)    = g22con(:,nU3)
!       g33conIonoS(:,j)    = g33con(:,nU3)
!
!   end do

    !write (*,*) 'Allocating ionosphere fit variables ... ' 

    !allocate ( hrBFnArrT ( nU1 * nU2, nBFns ), &
    !   alpha ( nBFns, nBFns ), &
    !   coeffs ( 1, nBFns ), &
    !   coeffsOut ( nBFns, 1 ), &
    !   coeffs_ ( 0:nBFns-1, 1 ), &
    !   coeffsOutS ( nBFns, 1 ), &
    !   coeffsT ( 1, nBFns ), &
    !   beta_ ( 1, nBFns ), &
    !   beta_S ( nBFns, 1 ), &
    !   h3conReal ( nU1, nU2 ), &
    !   beta_SVSOL ( nBFns, 1 ) )

    !write (*,*) 'Creating svd alpha array ...'

    !hrBFnArrT  = transpose ( hrBFnArr_h3 )
    !alpha  = matMul ( hrBFnArr_h3, transpose ( hrBFnArr_h3 ) )

    !write ( *,* ) hrBFnArr_h3

    u3MinRemoved    = u3 - minVal ( u3 )
    u3MinRemoved_   = u3_ - minVal ( u3_ )

    !   Create the azimuthally wrapping index array

    write (*,*) 'Creating azimuthally wrapped index array ... ' 

    jjj(1:nU2)  = (/ (j,j=1,nU2) /)
    jjj(0)  = nU2
    jjj(nU2+1)  = 1

    ! Create the azimuthal coord differenece array, 
    !   i.e., u2Diff

    do j = 1, nU2

        u2DiffF(j)  = ( u2(j) - u2( jjj(j+1) ) )
        if ( j == nU2 ) u2DiffF(j) = u2DiffF(j) - 2.0 * pi

        u2DiffB(j)  = ( u2( jjj(j-1) ) - u2( j ) )
        if ( j == 1 ) u2DiffB(j) = u2DiffB(j) - 2.0 * pi

    end do

    !   Setup per processor spatial grid chunks
    !   remembere nU3 is odd.
    
    iStart  = mpi_pId * nU1 / mpi_nP + 1
    iStop   = iStart + nU1 / mpi_nP - 1
    
    jStart  = mpi_pId * nU2 / mpi_nP + 1
    jStop   = jStart + nU2 / mpi_nP - 1
    
    kStart  = mpi_pId * (nU3-1) / mpi_nP + 1
    kStop   = kStart + (nU3-1) / mpi_nP - 1
 
    !plotFreq   = 1000d0
    progFreq    = 10d0

    write (*,*) 'STARTING TIME LOOP ... '

    timeLoop: &
    do t = 0, nt
        
        if ( mod ( t*dT, progFreq  ) < dT ) then
        
            startTimeOLD    = startTime
            call start_timer ( startTime )
        
        end if
    
        !   Apply h boundary conditions
    
        !   Driver
    
        driverFreq  = 50d-3!1./ (200.0 * dT )!
        driverAmp   = 10d-9 / u0!   
    
        if ( t*dT <= ( 0.5 / driverFreq ) ) then
    
        !   print, 'DRIVING...'!
            
            do k = nU3/2-1, nU3/2+1
                do j = 1, nU2 
            
                    u3Variation = exp ( -u3(k) ** 2 / 0.000001 )!
                    u2Variation = cos ( 2.0 * u2 ( j ) ) !exp   ( -( u2(j) / ( 2 * pi ) - 0.5 ) ** 2 / 0.005 )! 
                    !h3cov(1,j,k)   = driverAmp * sin ( 2.0 * pi * driverFreq ) &
                    !   * u3Variation * u2Variation / sqrt ( g33con(1,k) )! 
                    h3cov(1,j,k)    = driverAmp * &
                        exp ( - ( t*dT - 10.0 ) ** 2 / ( 2d0 * 10.0 ) ) &
                        * u3Variation * u2Variation / sqrt ( g33con(1,k) )! 

                end do
            end do
    
        else
        !   print, 'NOT DRIVING...'!
            h3cov(1,:,:)    = 0.0
        end if

        h2cov(1,:,:)    = 0.0!
        !h1cov(nU1,:,:) = ( 4.0 * h1cov(nU1-1,:,:) - h1cov(nU1-2,:,:) ) / 3.0! 

        !e2cov(1,:,1)   = ( 4.0 * e2cov(2,:,1) - e2cov(3,:,1) ) / 3.0!
        e2cov(nU1,:,2:nU3-1)    = 0.0!

       !write (*,*) iStart, iStop, jStart, jStop, kStart, kStop, nU1, nU2, nU3

    !   Update contravariant e field

    !   e1con except ionospheres
    
        iStart_ = iStart
        iStop_  = iStop
        jStart_ = jStart
        jStop_  = jStop 
        kStart_ = kStart
        kStop_  = kStop
        
        if ( kStart < 2 ) kStart_ = 2 
        if ( kStop > nU3-1 ) kStop_ = nU3-1 
 
        do k = kStart_, kStop_
            do j = jStart_, jStop_
                do i = iStart_, iStop_
                
                        e1con(i,j,k)    = e1con(i,j,k) + &
                            dT / ( epsilon_(i,k) * sqrtG(i,k) ) * &
                            ( ( h3cov(i,j,k) - h3cov(i,jjj(j+1),k) ) / u2DiffF(j) &
                            - ( h2cov(i,j,k-1) - h2cov(i,j,k) ) / ( u3MinRemoved_(k-1) - u3MinRemoved_(k) ) )!
                    
                end do
            end do
        end do

        !   e2con except inner boundary and ionospheres
        
        iStart_ = iStart
        iStop_  = iStop
        jStart_ = jStart
        jStop_  = jStop 
        kStart_ = kStart
        kStop_  = kStop
    
        if ( iStop > nU1-1 ) iStop_ = nU1-1 
        if ( kStart < 2 ) kStart_ = 2 
        if ( kStop > nU3-1 ) kStop_ = nU3-1 
        
        do k = kStart_, kStop_
            do j = jStart_, jStop_
                do i = iStart_, iStop_

                    e2con(i,j,k)    = e2con(i,j,k) + &
                        dT / ( epsilon_(i,k) * sqrtG(i,k) ) * &
                        ( ( h1cov(i,j,k-1) - h1cov(i,j,k) ) / ( u3MinRemoved_(k-1) - u3MinRemoved_(k) ) &
                        - ( h3cov(i,j,k) - h3cov(i+1,j,k) ) / ( u1(i) - u1(i+1) ) )!

                end do
            end do
        end do

    !   Calculate covariant e from contravariant e

    !   Do not update the points at the ionosphere as they
    !   have already been filled by the BC
    
    !   e1cov except ionospheres
     
        iStart_ = iStart
        iStop_  = iStop
        jStart_ = jStart
        jStop_  = jStop 
        kStart_ = kStart
        kStop_  = kStop
        
        if ( kStart < 2 ) kStart_ = 2 
        if ( kStop > nU3-1 ) kStop_ = nU3-1 
    
        do k = kStart_, kStop_
            do j = jStart_, jStop_
                do i = iStart_, iStop_ 

                    e1cov(i,j,k)    = 1.0 / g11con(i,k) * e1con(i,j,k)! &
                    
                    !   Do not update the points at the ionosphere as they
                    !   have already been filled by the BC
            
                end do
            end do
        end do
    
        !   e2cov except inner boundary and ionospheres

        iStart_ = iStart
        iStop_  = iStop
        jStart_ = jStart
        jStop_  = jStop 
        kStart_ = kStart
        kStop_  = kStop
        
        if ( kStart < 2 ) kStart_ = 2 
        if ( kStop > nU3-1 ) kStop_ = nU3-1 
        if ( iStop > nU1-1 ) iStop_ = nU1-1
 
        do k = kStart_, kStop_ 
            do j = jStart_, jStop_
                do i = iStart_, iStop_ 

                        e2cov(i,j,k)    = g22cov(i,k) * e2con(i,j,k)!
    
                end do
            end do
        end do
       
        !   Pass the covariant e field components to the neighbouring
        !   processor. Also, recieve the e field components from other
        !   neighbour.

        call mpi_barrier ( mpi_comm_cart, mpi_iErr ) 

        call mpi_cart_coords ( mpi_comm_cart, mpi_pId, 3, mpi_coords, mpi_iErr )
        mpi_coords_dest = mpi_coords + (/ 1, 0, 0 /)
        mpi_coords_src  = mpi_coords + (/ -1, 0, 0 /)
        if (mpi_coords_dest(1) <= mpi_nP_U1-1 ) & 
            call mpi_cart_rank ( mpi_comm_cart, mpi_coords_dest, mpi_pId_dest, mpi_iErr )
        if (mpi_coords_src(1) >= 0 ) &
            call mpi_cart_rank ( mpi_comm_cart, mpi_coords_src, mpi_pId_src, mpi_iErr )
       
        if ( mpi_coords(1) == 0 ) then
            e2cov(iStop,:,:) = 4.0
        end if
 
        mpi_tag_e2cov = 1 
        if ( mpi_coords(1) == 0 ) then 
            call mpi_send ( e2cov(iStop,jStart,kStart), 1, mpi_type_2D_U2U3, &
                mpi_pId_dest, mpi_tag_e2cov, mpi_comm_cart, mpi_iErr)
        else if ( mpi_coords(1) == mpi_nP_U1-1 ) then
            call mpi_recv ( e2cov(iStart-1,jStart,kStart), 1, mpi_type_2D_U2U3, &
                mpi_pId_src, mpi_tag_e2cov, mpi_comm_cart, mpi_status, mpi_iErr ) 
        else
            call mpi_sendRecv ( e2cov(iStop,jStart,kStart), 1, mpi_type_2D_U2U3, & 
                mpi_pId_dest, mpi_tag_e2cov, e2cov(iStart-1,jStart,kStart), &
                1, mpi_type_2D_U2U3, mpi_pId_src, mpi_tag_e2cov, mpi_comm_cart, mpi_status, mpi_iErr )
        end if
       
        call mpi_barrier ( mpi_comm_cart, mpi_iErr ) 

        write (*,*)
        if ( mpi_coords(1) == 0 ) then
            write (*,*) mpi_coords
            do k = 1, nU3-1 
                write (*,*) e2cov(iStop,:,k)
            end do 
        end if 
        call mpi_barrier ( mpi_comm_cart, mpi_iErr ) 


        if ( mpi_coords(1) == 1 ) then 
            write (*,*) mpi_coords
            do k = 1, nU3-1
                write (*,*) e2cov(iStart-1,:,k)
            end do 
        end if

 
        !   Update contravariant h field
    
        !   h1con except inner boundary

        do k = 1, nU3-1
            do j = 1, nU2
                do i = 1, nU1-1

                        h1con(i,j,k)    = h1con(i,j,k) - &
                            dT / ( u0 * sqrtG_(i,k) ) * &
                            ( -( e2cov(i,j,k) - e2cov(i,j,k+1) ) / ( u3MinRemoved(k) - u3MinRemoved(k+1) ) )!
    
                end do
            end do
        end do

        !   h2con all locations

        do k = 1, nU3-1
            do j = 1, nU2
                do i = 1, nU1
        
                    h2con(i,j,k)    = h2con(i,j,k) - &
                        dT / ( u0 * sqrtG_(i,k) ) * &
                        ( ( e1cov(i,j,k) - e1cov(i,j,k+1) ) / ( u3MinRemoved(k) - u3MinRemoved(k+1) ) )!&
                
                    end do
            end do
        end do
    
        !   h3con except outer boundary

        do k = 1, nU3
            do j = 1, nU2
                do i = 2, nU1
    
                        h3con(i,j,k)    = h3con(i,j,k) - &
                            dT / ( u0 * sqrtG(i,k) ) * &
                            ( ( e2cov(i-1,j,k) - e2cov(i,j,k) ) / ( u1(i-1) - u1(i) ) &
                            - ( e1cov(i,jjj(j-1),k) - e1cov(i,j,k) ) / u2DiffB(j)   )!
    
                end do
            end do
        end do

    !write (*,*)
        !write (*,*) sqrtG(:,1)
        !read (*,*)

        !   Set the dh3/du1=0 @ outer boundary in the two ionospheres
    
        h3con(1,:,:)    = ( 4.0 * h3con(2,:,:) - h3con(3,:,:) ) / 3.0!
        h3con(1,:,:)    = ( 4.0 * h3con(2,:,:) - h3con(3,:,:) ) / 3.0! 

        e1cov(:,:,1)    = 0.0!E1covIono / sqrt ( g11conIono )! 
        e2cov(:,:,1)    = 0.0!E2covIono / sqrt ( g22conIono )!
    

        if ( mod ( t*dT, progFreq ) < dT ) then  

            timeTaken   = end_timer ( startTimeOLD )
            eta = timeTaken / progFreq * ( tRun - t*dT ) / 3600.0
            
            if ( eta < 1 ) then 
                write (*,'(a12,f5.1,a6,f4.1,a8)') 'Model time: ', t*dT, &
                    ' ETA: ', eta * 60.0, &
                    ' [MINS]'
            else
                write (*,'(a12,f5.1,a6,f4.1,a8)') 'Model time: ', t*dT, &
                    ' ETA: ', eta, &
                    ' [HOURS]'
            end if

        end if

        !read (*,*)

    end do timeLoop

    call mpi_finalize ( mpi_iErr )

end program mhd3d
