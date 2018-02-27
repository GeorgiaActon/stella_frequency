module parallel_streaming

  implicit none

  public :: init_parallel_streaming, finish_parallel_streaming
  public :: advance_parallel_streaming_explicit
  public :: advance_parallel_streaming_implicit
  public :: stream_tridiagonal_solve
  public :: parallel_streaming_initialized
  public :: stream, stream_c, stream_sign
  public :: time_parallel_streaming

  private

  interface center_zed
     module procedure center_zed_segment_real
     module procedure center_zed_segment_complex
     module procedure center_zed_extended
  end interface

  logical :: parallel_streaming_initialized = .false.

  integer, dimension (:), allocatable :: stream_sign
  real, dimension (:,:,:), allocatable :: stream, stream_c
  real, dimension (:,:), allocatable :: stream_tri_a1, stream_tri_a2
  real, dimension (:,:), allocatable :: stream_tri_b1, stream_tri_b2
  real, dimension (:,:), allocatable :: stream_tri_c1, stream_tri_c2
  real, dimension (:,:,:), allocatable :: gradpar_c

  real, dimension (2) :: time_parallel_streaming

contains

  subroutine init_parallel_streaming

    use finite_differences, only: fd3pt
    use stella_time, only: code_dt
    use species, only: spec, nspec
    use vpamu_grids, only: nvgrid, nvpa
    use vpamu_grids, only: vpa
    use zgrid, only: nzgrid, nztot
    use geometry, only: gradpar, nalpha
    use run_parameters, only: stream_implicit
    use run_parameters, only: stream_cell
    use run_parameters, only: include_parallel_streaming

    implicit none

    integer :: iv, is

    if (parallel_streaming_initialized) return
    parallel_streaming_initialized = .true.

    if (.not.allocated(stream)) allocate (stream(-nzgrid:nzgrid,-nvgrid:nvgrid,nspec)) ; stream = 0.
    if (.not.allocated(stream_sign)) allocate (stream_sign(-nvgrid:nvgrid)) ; stream_sign = 0

    ! sign of stream corresponds to appearing on RHS of GK equation
    ! i.e., this is the factor multiplying dg/dz on RHS of equation
    if (include_parallel_streaming) then
       stream = -code_dt*spread(spread(spec%stm,1,nztot),2,nvpa) &
            * spread(spread(vpa,1,nztot)*spread(gradpar(1,:),2,nvpa),3,nspec)
    else
       stream = 0.0
    end if

    ! stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
    ! NB: stream_sign = -1 corresponds to positive advection velocity
    do iv = -nvgrid, nvgrid
       stream_sign(iv) = int(sign(1.0,stream(0,iv,1)))
    end do
    ! vpa = 0 is special case
    stream_sign(0) = 0

    if (stream_implicit) then
       call init_invert_stream_operator
       if (.not.allocated(stream_c)) allocate (stream_c(-nzgrid:nzgrid,-nvgrid:nvgrid,nspec))
       stream_c = stream
       do is = 1, nspec
          do iv = -nvgrid, nvgrid
             call center_zed (iv, stream_c(:,iv,is))
          end do
       end do
       if (.not.allocated(gradpar_c)) allocate (gradpar_c(nalpha,-nzgrid:nzgrid,-1:1))
       gradpar_c = spread(gradpar,3,3)
       ! get gradpar centred in zed for negative vpa (affects upwinding)
       call center_zed(-1,gradpar_c(1,:,-1))
       ! get gradpar centred in zed for positive vpa (affects upwinding)
       call center_zed(1,gradpar_c(1,:,1))
       if (stream_cell) then
          stream = stream_c
       end if
    end if

  end subroutine init_parallel_streaming

  subroutine init_invert_stream_operator

    use zgrid, only: delzed
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: nsegments
    use run_parameters, only: zed_upwind, time_upwind
    use run_parameters, only: stream_cell

    implicit none

    integer :: nz, nseg_max

    nz = maxval(iz_up-iz_low)
    nseg_max = maxval(nsegments)

    if (.not.allocated(stream_tri_a1)) then
       allocate (stream_tri_a1(nz*nseg_max+1,-1:1)) ; stream_tri_a1 = 0.
       allocate (stream_tri_a2(nz*nseg_max+1,-1:1)) ; stream_tri_a2 = 0.
       allocate (stream_tri_b1(nz*nseg_max+1,-1:1)) ; stream_tri_b1 = 1.
       allocate (stream_tri_b2(nz*nseg_max+1,-1:1)) ; stream_tri_b2 = 0.
       allocate (stream_tri_c1(nz*nseg_max+1,-1:1)) ; stream_tri_c1 = 0.
       allocate (stream_tri_c2(nz*nseg_max+1,-1:1)) ; stream_tri_c2 = 0.
    end if

    if (stream_cell) then
       ! corresponds to sign of stream term positive on RHS of equation
       ! i.e., negative parallel advection speed
       ! NB: assumes equal spacing in zed
       stream_tri_b1(:,1) = 0.5*(1.0+zed_upwind)
       stream_tri_b2(:,1) = -1.0/delzed(0)
       stream_tri_c1(:nz*nseg_max,1) = 0.5*(1.0-zed_upwind)
       stream_tri_c2(:nz*nseg_max,1) = 1.0/delzed(0)
       ! corresponds to sign of stream term negative on RHS of equation
       ! NB: assumes equal spacing in zed
       stream_tri_b1(:,-1) = 0.5*(1.0+zed_upwind)
       stream_tri_b2(:,-1) = 1.0/delzed(0)
       stream_tri_a1(2:,-1) = 0.5*(1.0-zed_upwind)
       stream_tri_a2(2:,-1) = -1.0/delzed(0)
    else
       ! corresponds to sign of stream term positive on RHS of equation
       ! i.e., negative parallel advection speed
       ! NB: assumes equal spacing in zed
       stream_tri_a2(2:,1) = -0.5*(1.0-zed_upwind)/delzed(0)
       stream_tri_b2(2:,1) = -zed_upwind/delzed(0)
       stream_tri_c2(2:,1) = (1.0+zed_upwind)*0.5/delzed(0)
       ! must treat boundary carefully
       stream_tri_b2(1,1) = -1.0/delzed(0)
       stream_tri_c2(1,1) = 1.0/delzed(0)
       ! corresponds to sign of stream term negative on RHS of equation
       ! NB: assumes equal spacing in zed
       stream_tri_a2(:nz*nseg_max,-1) = -0.5*(1.0+zed_upwind)/delzed(0)
       stream_tri_b2(:nz*nseg_max,-1) = zed_upwind/delzed(0)
       stream_tri_c2(:nz*nseg_max,-1) = 0.5*(1.0-zed_upwind)/delzed(0)
       ! must treat boundary carefully
       stream_tri_a2(nz*nseg_max+1,-1) = -1.0/delzed(0)
       stream_tri_b2(nz*nseg_max+1,-1) = 1.0/delzed(0)
    end if

    stream_tri_a2 = 0.5*(1.0+time_upwind)*stream_tri_a2
    stream_tri_b2 = 0.5*(1.0+time_upwind)*stream_tri_b2
    stream_tri_c2 = 0.5*(1.0+time_upwind)*stream_tri_c2

  end subroutine init_invert_stream_operator

  subroutine advance_parallel_streaming_explicit (g, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use job_manage, only: time_message
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx, zonal_mode
    use dist_fn_arrays, only: aj0x
    use fields_arrays, only: phi
    use vpamu_grids, only: ztmax, maxwell_mu

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ivmu, iv, imu, is
    complex, dimension (:,:,:,:), allocatable :: g0, g1

    allocate (g0(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate (g1(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! parallel streaming stays in ky,kx,z space with ky,kx,z local
    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')
    ! get dg/dz, with z the parallel coordinate and store in g0
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       g0(:,:,:,ivmu) = aj0x(:,:,:,ivmu)*phi
    end do
    call get_dgdz (g0, g1)
    call get_dgdz (g, g0)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       g0(:,:,:,ivmu) = g0(:,:,:,ivmu) + g1(:,:,:,ivmu)*ztmax(iv,is) &
            * spread(spread(maxwell_mu(1,:,imu),1,naky),2,nakx)
    end do

    ! multiply dg/dz with vpa*(b . grad z) and add to source (RHS of GK equation)
    call add_stream_term (g0, gout)
    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')
    deallocate (g0, g1)

  end subroutine advance_parallel_streaming_explicit

  subroutine get_dgdz (g, dgdz)

    use finite_differences, only: third_order_upwind_zed
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx
    use zgrid, only: nzgrid, delzed
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones
    use kt_grids, only: naky

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dgdz

    integer :: ivmu, iseg, ie, iky, iv
    complex, dimension (2) :: gleft, gright

    ! FLAG -- assuming delta zed is equally spaced below!
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       ! no need to calculate dgdz for vpa=0, as will be multipled by vpa
       if (iv==0) then
          dgdz(:,:,:,ivmu) = 0.
          cycle
       end if
       do iky = 1, naky
          do ie = 1, neigen(iky)
             do iseg = 1, nsegments(ie,iky)
                ! first fill in ghost zones at boundaries in g(z)
                call fill_zed_ghost_zones (iseg, ie, iky, g(:,:,:,ivmu), gleft, gright)
                ! now get dg/dz
                call third_order_upwind_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
                     g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),ivmu), &
                     delzed(0), stream_sign(iv), gleft, gright, &
                     dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),ivmu))
             end do
          end do
       end do
    end do

  end subroutine get_dgdz

  subroutine add_stream_term (g, src)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx!, zonal_mode

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: iv, is, ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       src(:,:,:,ivmu) = src(:,:,:,ivmu) + spread(spread(stream(:,iv,is),1,naky),2,nakx)*g(:,:,:,ivmu)
!       if (zonal_mode(1)) src(1,:,-nzgrid,ivmu) = src(1,:,nzgrid,ivmu)
    end do

  end subroutine add_stream_term

  subroutine advance_parallel_streaming_implicit (g, phi, apar)

    use mp, only: proc0
    use job_manage, only: time_message
    use stella_layouts, only: vmu_lo, iv_idx
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx
    use dist_fn_arrays, only: g1
    use run_parameters, only: stream_matrix_inversion
    use fields, only: advance_fields, fields_updated

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:), intent (in out) :: phi, apar

    integer :: ivmu, iv
    complex, dimension (:,:,:), allocatable :: phi1

    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')

    allocate (phi1(naky,nakx,-nzgrid:nzgrid))

    ! save the incoming g and phi, as they will be needed later
    g1 = g
    phi1 = phi

    ! if iv=0, vpa=0. in this case, dg/dt = 0, so g_out = g_in

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu) ; if (iv==0) cycle

       ! obtain RHS of inhomogeneous GK eqn;
       ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g_{inh}^{n+1} 
       ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n} 
       ! + (1-alph)/2*dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d<phi^{n}>/dz
       call get_gke_rhs (ivmu, g1(:,:,:,ivmu), phi1, phi, g(:,:,:,ivmu), eqn='inhomogeneous')
       
       if (stream_matrix_inversion) then
          ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
          ! g = RHS is input and overwritten by g = g_{inh}^{n+1}
          call invert_parstream (ivmu, g(:,:,:,ivmu))
       else
          call sweep_g_zed (ivmu, g(:,:,:,ivmu))
       end if
    end do

    fields_updated = .false.

    ! we now have g_{inh}^{n+1}
    ! calculate associated fields (phi_{inh}^{n+1})
    call advance_fields (g, phi, apar, dist='gbar')

    ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
    ! phi = phi_{inh}^{n+1} is input and overwritten by phi = phi^{n+1}
    call invert_parstream_response (phi)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu) ; if (iv==0) cycle
    
       ! now have phi^{n+1} for non-negative kx
       ! obtain RHS of GK eqn; 
       ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1} 
       ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n} 
       ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>)
       call get_gke_rhs (ivmu, g1(:,:,:,ivmu), phi1, phi, g(:,:,:,ivmu), eqn='full')
    
       if (stream_matrix_inversion) then
          ! solve (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1} = RHS
          ! g = RHS is input and overwritten by g = g^{n+1}
          call invert_parstream (ivmu, g(:,:,:,ivmu))
       else
          call sweep_g_zed (ivmu, g(:,:,:,ivmu))
       end if
    end do

    deallocate (phi1)

    if (proc0) call time_message(.false.,time_parallel_streaming,' Stream advance')

  end subroutine advance_parallel_streaming_implicit

  subroutine get_gke_rhs (ivmu, gold, phiold, phi, g, eqn)

    use stella_time, only: code_dt
    use zgrid, only: nzgrid
    use species, only: spec
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: naky, nakx
    use kt_grids, only: zonal_mode
    use dist_fn_arrays, only: aj0x
    use vpamu_grids, only: vpa, ztmax, maxwell_mu
    use geometry, only: gradpar, nalpha
    use neoclassical_terms, only: include_neoclassical_terms
    use neoclassical_terms, only: dfneo_dvpa
    use run_parameters, only: stream_cell, time_upwind

    implicit none

    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:), intent (in) :: gold
    complex, dimension (:,:,-nzgrid:), intent (in) :: phiold, phi
    complex, dimension (:,:,-nzgrid:), intent (in out) :: g
    character (*), intent (in) :: eqn

    integer :: iv, imu, is, iz, ikx
    real :: tupwnd1, tupwnd2, fac
    real, dimension (:), allocatable :: vpadf0dE_fac, vpadf0dE_fac_zf
    real, dimension (:,:), allocatable :: gp, gpz
    complex, dimension (:,:,:), allocatable :: dgdz, dphidz

    allocate (vpadf0dE_fac(-nzgrid:nzgrid))
    allocate (vpadf0dE_fac_zf(-nzgrid:nzgrid))
    allocate (gp(nalpha,-nzgrid:nzgrid))
    allocate (gpz(nalpha,-nzgrid:nzgrid))
    allocate (dgdz(naky,nakx,-nzgrid:nzgrid))
    allocate (dphidz(naky,nakx,-nzgrid:nzgrid))

    tupwnd1 = 0.5*(1.0-time_upwind)
    if (eqn=='full') then
       tupwnd2 = 0.5*(1.0+time_upwind)
    else
       tupwnd2 = 0.0
    end if

    ! now have phi^{n+1} for non-negative kx
    ! obtain RHS of GK eqn; 
    ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1} 
    ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n} 
    ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>
    iv = iv_idx(vmu_lo,ivmu)
    imu = imu_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)

    ! obtain dg^{n}/dz and store in dgdz
    ! NB: could eliminate this calculation at the expense of memory
    ! as this was calculated previously
    if (stream_cell) then
       call get_dzed_cell (iv,gold,dgdz)
    else
       call get_dzed (iv,gold,dgdz)
    end if

    ! get <phi> = (1+alph)/2*<phi^{n+1}> + (1-alph)/2*<phi^{n}>
    g = aj0x(:,:,:,ivmu)*(tupwnd1*phiold+tupwnd2*phi)
    ! obtain d<phi>/dz and store in dphidz
    if (stream_cell) then
       call get_dzed_cell (iv,g,dphidz)
    else
       call get_dzed (iv,g,dphidz)
    end if

    fac = code_dt*spec(is)%stm

    ! NB: could do this once at beginning of simulation to speed things up
    ! this is vpa*Z/T*exp(-vpa^2)
    vpadf0dE_fac = vpa(iv)*ztmax(iv,is)*maxwell_mu(1,:,imu)
    ! if including neoclassical correction to equilibrium distribution function
    ! then must also account for -vpa*dF_neo/dvpa
    if (include_neoclassical_terms) then
       do iz = -nzgrid, nzgrid
          vpadf0dE_fac(iz) = vpadf0dE_fac(iz)-0.5*dfneo_dvpa(iz,ivmu)
       end do
    end if

    g = gold
    if (stream_cell) then
       call center_zed (iv,g)
       call center_zed (iv,vpadf0dE_fac)
       if (iv < 0) then
          gp = gradpar_c(:,:,-1)
       else
          gp = gradpar_c(:,:,1)
       end if
       gpz = gp
       vpadf0dE_fac_zf = vpadf0dE_fac
    else
       gp = gradpar
       ! treat zonal flow specially
       if (zonal_mode(1)) then
          do ikx = 1, nakx
             call center_zed (iv,g(1,ikx,:))
          end do
          vpadf0dE_fac_zf = vpadf0dE_fac
          call center_zed (iv,vpadf0dE_fac_zf)
          if (iv < 0) then
             gpz = gradpar_c(:,:,-1)
          else
             gpz = gradpar_c(:,:,1)
          end if
       end if
    end if

    ! construct RHS of GK eqn
    if (zonal_mode(1)) then
       do iz = -nzgrid, nzgrid
          g(1,:,iz) = g(1,:,iz) - fac*gpz(1,iz) &
               * (tupwnd1*vpa(iv)*dgdz(1,:,iz) + vpadf0dE_fac_zf(iz)*dphidz(1,:,iz))
          g(2:,:,iz) = g(2:,:,iz) - fac*gp(1,iz) &
               * (tupwnd1*vpa(iv)*dgdz(2:,:,iz) + vpadf0dE_fac(iz)*dphidz(2:,:,iz))
       end do
    else
       do iz = -nzgrid, nzgrid
          g(:,:,iz) = g(:,:,iz) - fac*gp(1,iz) &
               * (tupwnd1*vpa(iv)*dgdz(:,:,iz) + vpadf0dE_fac(iz)*dphidz(:,:,iz))
       end do
    end if

    deallocate (vpadf0dE_fac, vpadf0dE_fac_zf, gp, gpz)
    deallocate (dgdz, dphidz)

  end subroutine get_gke_rhs

  ! solve (I + (1+alph)/2*dt*vpa . grad)g^{n+1} = RHS
  ! g = RHS is input and overwritten by g = g^{n+1}
  subroutine invert_parstream (ivmu, g)

    use zgrid, only: nzgrid
    use extended_zgrid, only: neigen
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use extended_zgrid, only: map_to_extended_zgrid
    use extended_zgrid, only: map_from_extended_zgrid
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use kt_grids, only: naky
    use kt_grids, only: zonal_mode

    implicit none

    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:), intent (in out) :: g

    integer :: iv, is
    integer :: iky, ie
    integer :: ulim, sgn
    complex, dimension (:), allocatable :: gext

    iv = iv_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)
    sgn = stream_sign(iv)

    do iky = 1, naky
       if (zonal_mode(iky)) then
          call sweep_zed_zonal (iv, is, sgn, g(iky,:,:))
       else
          do ie = 1, neigen(iky)
             allocate (gext(nsegments(ie,iky)*nzed_segment+1))
             ! get g on extended domain in zed
             call map_to_extended_zgrid (ie, iky, g(iky,:,:), gext, ulim)
             ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
             call stream_tridiagonal_solve (iky, ie, iv, is, gext(:ulim))
             ! extract g from extended domain in zed
             call map_from_extended_zgrid (ie, iky, gext, g(iky,:,:))
             deallocate (gext)
          end do
       end if
    end do

  end subroutine invert_parstream

  subroutine stream_tridiagonal_solve (iky, ie, iv, is, g)

    use finite_differences, only: tridag
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment

    implicit none

    integer, intent (in) :: iky, ie, iv, is
    complex, dimension (:), intent (in out) :: g

    integer :: iseg, llim, ulim, n
    integer :: nz, nseg_max, sgn, n_ext
    real, dimension (:), allocatable :: a, b, c

    ! avoid double-counting at boundaries between 2pi segments
    nz = nzed_segment
    nseg_max = nsegments(ie,iky)
    sgn = stream_sign(iv)

    n_ext = nseg_max*nz+1
    allocate (a(n_ext))
    allocate (b(n_ext))
    allocate (c(n_ext))

    iseg = 1
    llim = 1 ; ulim = nz+1
    a(llim:ulim) = stream_tri_a1(llim:ulim,sgn) &
         - stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_a2(llim:ulim,sgn)
    b(llim:ulim) = stream_tri_b1(llim:ulim,sgn) &
         -stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_b2(llim:ulim,sgn)
    c(llim:ulim) = stream_tri_c1(llim:ulim,sgn) &
         -stream(iz_low(iseg):iz_up(iseg),iv,is)*stream_tri_c2(llim:ulim,sgn)

    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          llim = ulim+1
          ulim = llim+nz-1
          a(llim:ulim) = stream_tri_a1(llim:ulim,sgn) &
               - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_a2(llim:ulim,sgn)
          b(llim:ulim) = stream_tri_b1(llim:ulim,sgn) &
               - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_b2(llim:ulim,sgn)
          c(llim:ulim) = stream_tri_c1(llim:ulim,sgn) &
               - stream(iz_low(iseg)+1:iz_up(iseg),iv,is)*stream_tri_c2(llim:ulim,sgn)
       end do
    end if
    n = size(stream_tri_a1,1)
    a(ulim) = stream_tri_a1(n,sgn)-stream(iz_up(nsegments(ie,iky)),iv,is)*stream_tri_a2(n,sgn)
    b(ulim) = stream_tri_b1(n,sgn)-stream(iz_up(nsegments(ie,iky)),iv,is)*stream_tri_b2(n,sgn)
    c(ulim) = 0. ! this line should not be necessary, as c(ulim) should not be accessed by tridag
    call tridag (1, a(:ulim), b(:ulim), c(:ulim), g)

    deallocate (a, b, c)

  end subroutine stream_tridiagonal_solve

  ! g= RHS of gke is input
  ! g = g^{n+1} is output
  subroutine sweep_g_zed (ivmu, g)

    use zgrid, only: nzgrid, delzed
    use extended_zgrid, only: neigen, nsegments, nzed_segment
    use extended_zgrid, only: map_to_extended_zgrid
    use extended_zgrid, only: map_from_extended_zgrid
    use kt_grids, only: naky
    use kt_grids, only: zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use run_parameters, only: zed_upwind, time_upwind

    implicit none

    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:), intent (in out) :: g
    
    integer :: iv, is
    integer :: iky, ie
    integer :: ulim, sgn
    integer :: iz, izext, iz1, iz2
    real :: fac1, fac2
    complex, dimension (:), allocatable :: gext

    iv = iv_idx (vmu_lo,ivmu)
    is = is_idx (vmu_lo,ivmu)
    sgn = stream_sign(iv)
    ! will sweep to right (positive vpa) or left (negative vpa)
    ! and solve for g on the extended z-grid
    do iky = 1, naky
       if (zonal_mode(iky)) then
          call sweep_zed_zonal (iv, is, sgn, g(iky,:,:))
       else
          do ie = 1, neigen(iky)
             allocate (gext(nsegments(ie,iky)*nzed_segment+1))
             ! get g on extended domain in zed
             call map_to_extended_zgrid (ie, iky, g(iky,:,:), gext, ulim)
             if (sgn < 0) then
                iz1 = 1 ; iz2 = ulim
             else
                iz1 = ulim ; iz2 = 1
             end if
             izext = iz1 ; iz = sgn*nzgrid
             fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
             gext(izext) = gext(izext)*2.0/fac1
             do izext = iz1-sgn, iz2, -sgn
                if (iz == -sgn*nzgrid) then
                   iz = sgn*nzgrid-sgn
                else
                   iz = iz - sgn
                end if
                fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
                fac2 = 1.0-zed_upwind-sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
                gext(izext) = (-gext(izext+sgn)*fac2 + 2.0*gext(izext))/fac1
             end do
             ! extract g from extended domain in zed
             call map_from_extended_zgrid (ie, iky, gext, g(iky,:,:))
             deallocate (gext)
          end do
       end if
    end do

  end subroutine sweep_g_zed

  subroutine sweep_zed_zonal (iv, is, sgn, g)

    use zgrid, only: nzgrid, delzed, nztot
    use kt_grids, only: nakx
    use run_parameters, only: zed_upwind, time_upwind

    implicit none

    integer, intent (in) :: iv, is, sgn
    complex, dimension (:,-nzgrid:), intent (in out) :: g

    integer :: iz, iz1, iz2
    real :: fac1, fac2
    complex, dimension (:), allocatable :: gcf
    complex, dimension (:,:), allocatable :: gpi

    allocate (gpi(nakx,-nzgrid:nzgrid))
    allocate (gcf(-nzgrid:nzgrid))
    ! ky=0 is 2pi periodic (no extended zgrid)
    ! decompose into complementary function + particular integral
    ! zero BC for particular integral
    ! unit BC for complementary function (no source)
    if (sgn < 0) then
       iz1 = -nzgrid ; iz2 = nzgrid
    else
       iz1 = nzgrid ; iz2 = -nzgrid
    end if
    gpi(:,iz1) = 0. ; gcf(iz1) = 1.
    do iz = iz1-sgn, iz2, -sgn
       fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       fac2 = 1.0-zed_upwind-sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       gpi(:,iz) = (-gpi(:,iz+sgn)*fac2 + 2.0*g(:,iz))/fac1
       gcf(iz) = -gcf(iz+sgn)*fac2/fac1
    end do
    ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
    g = gpi + (spread(gpi(:,iz2),2,nztot)/(1.-gcf(iz2)))*spread(gcf,1,nakx)
    deallocate (gpi, gcf)

  end subroutine sweep_zed_zonal

  subroutine invert_parstream_response (phi)
 
    use linear_solve, only: lu_back_substitution
    use zgrid, only: nzgrid
    use extended_zgrid, only: neigen
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use extended_zgrid, only: map_to_extended_zgrid
    use extended_zgrid, only: map_from_extended_zgrid
    use extended_zgrid, only: ikxmod
    use kt_grids, only: naky, zonal_mode
    use fields_arrays, only: response_matrix

    implicit none

    complex, dimension (:,:,-nzgrid:), intent (in out) :: phi

    integer :: iky, ie, ulim
    integer :: ikx
    complex, dimension (:), allocatable :: gext
    
    ! need to put the fields into extended zed grid
    do iky = 1, naky
       ! avoid double counting of periodic endpoints for zonal modes
       if (zonal_mode(iky)) then
          do ie = 1, neigen(iky)
             ikx = ikxmod(1,ie,iky)
             call lu_back_substitution (response_matrix(iky)%eigen(ie)%zloc, &
                  response_matrix(iky)%eigen(ie)%idx, phi(iky,ikx,:nzgrid-1))
             phi(iky,ikx,nzgrid) = phi(iky,ikx,-nzgrid)
          end do
       else
          do ie = 1, neigen(iky)
             ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
             allocate (gext(nsegments(ie,iky)*nzed_segment+1))
             call map_to_extended_zgrid (ie, iky, phi(iky,:,:), gext, ulim)
             call lu_back_substitution (response_matrix(iky)%eigen(ie)%zloc, &
                  response_matrix(iky)%eigen(ie)%idx, gext)
             call map_from_extended_zgrid (ie, iky, gext, phi(iky,:,:))
             deallocate (gext)
          end do
       end if
    end do
    
  end subroutine invert_parstream_response

  subroutine get_dzed (iv, g, dgdz)

    use finite_differences, only: fd_variable_upwinding_zed, fd_cell_centres_zed
    use kt_grids, only: naky
    use kt_grids, only: zonal_mode
    use zgrid, only: nzgrid, delzed
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones
    use run_parameters, only: zed_upwind

    implicit none

    integer, intent (in) :: iv
    complex, dimension (:,:,-nzgrid:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:), intent (out) :: dgdz

    integer :: iky, ie, iseg
    complex, dimension (2) :: gleft, gright

    do iky = 1, naky
       do ie = 1, neigen(iky)
          do iseg = 1, nsegments(ie,iky)
             ! first fill in ghost zones at boundaries in g(z)
             call fill_zed_ghost_zones (iseg, ie, iky, g, gleft, gright)
             ! treat zonal flow specially
             if (zonal_mode(iky)) then
                ! get finite difference approximation for dg/dz at cell centres
                ! iv > 0 corresponds to positive vpa, iv < 0 to negative vpa
                call fd_cell_centres_zed (iz_low(iseg), &
                     g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)), &
                     delzed(0), stream_sign(iv), gleft(2), gright(1), &
                     dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)))
             else
                ! get finite difference approximation for dg/dz
                ! with mixture of centered and upwinded scheme
                ! mixture controlled by zed_upwind (0 = centered, 1 = upwind)
                ! iv > 0 corresponds to positive vpa, iv < 0 to negative vpa
                call fd_variable_upwinding_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
                     g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)), &
                     delzed(0), stream_sign(iv), zed_upwind, gleft, gright, &
                     dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)))
             end if
          end do
       end do
    end do

  end subroutine get_dzed

  subroutine get_dzed_cell (iv, g, dgdz)

    use finite_differences, only: fd_cell_centres_zed
    use kt_grids, only: naky
    use zgrid, only: nzgrid, delzed
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones

    implicit none

    integer, intent (in) :: iv
    complex, dimension (:,:,-nzgrid:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:), intent (out) :: dgdz

    integer :: iky, ie, iseg
    complex, dimension (2) :: gleft, gright

    do iky = 1, naky
       do ie = 1, neigen(iky)
          do iseg = 1, nsegments(ie,iky)
             ! first fill in ghost zones at boundaries in g(z)
             call fill_zed_ghost_zones (iseg, ie, iky, g, gleft, gright)
             ! get finite difference approximation for dg/dz at cell centres
             ! iv > 0 corresponds to positive vpa, iv < 0 to negative vpa
             call fd_cell_centres_zed (iz_low(iseg), &
                  g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)), &
                  delzed(0), stream_sign(iv), gleft(2), gright(1), &
                  dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)))
          end do
       end do
    end do

  end subroutine get_dzed_cell

  subroutine center_zed_extended (iv, g)

    use finite_differences, only: cell_centres_zed
    use kt_grids, only: naky, nakx
    use zgrid, only: nzgrid
    use extended_zgrid, only: neigen, nsegments
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: fill_zed_ghost_zones
    use run_parameters, only: zed_upwind

    implicit none

    integer, intent (in) :: iv
    complex, dimension (:,:,-nzgrid:), intent (in out) :: g

    integer :: iky, ie, iseg
    complex, dimension (2) :: gleft, gright
    complex, dimension (:,:), allocatable :: gc

    allocate (gc(nakx,-nzgrid:nzgrid))

    do iky = 1, naky
       do ie = 1, neigen(iky)
          do iseg = 1, nsegments(ie,iky)
             ! first fill in ghost zones at boundaries in g(z)
             call fill_zed_ghost_zones (iseg, ie, iky, g, gleft, gright)
             ! get cell centres values 
             call cell_centres_zed (iz_low(iseg), &
                  g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)), &
                  zed_upwind, stream_sign(iv), gleft(2), gright(1), &
                  gc(ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg)))
          end do
       end do
       g(iky,:,:) = gc
    end do

    deallocate (gc)

  end subroutine center_zed_extended

  subroutine center_zed_segment_real (iv, g)

    use zgrid, only: nzgrid
    use run_parameters, only: zed_upwind

    integer, intent (in) :: iv
    real, dimension (-nzgrid:), intent (in out) :: g

    if (stream_sign(iv) > 0) then
       g(:nzgrid-1) = 0.5*((1.+zed_upwind)*g(:nzgrid-1) + (1.-zed_upwind)*g(-nzgrid+1:))
       g(nzgrid) = g(-nzgrid)
    else
       g(-nzgrid+1:) = 0.5*((1.-zed_upwind)*g(:nzgrid-1) + (1.+zed_upwind)*g(-nzgrid+1:))
       g(-nzgrid) = g(nzgrid)
    end if

  end subroutine center_zed_segment_real

  subroutine center_zed_segment_complex (iv, g)

    use zgrid, only: nzgrid
    use run_parameters, only: zed_upwind

    integer, intent (in) :: iv
    complex, dimension (-nzgrid:), intent (in out) :: g

    if (stream_sign(iv) > 0) then
       g(:nzgrid-1) = 0.5*((1.+zed_upwind)*g(:nzgrid-1) + (1.-zed_upwind)*g(-nzgrid+1:))
       g(nzgrid) = g(-nzgrid)
    else
       g(-nzgrid+1:) = 0.5*((1.-zed_upwind)*g(:nzgrid-1) + (1.+zed_upwind)*g(-nzgrid+1:))
       g(-nzgrid) = g(nzgrid)
    end if

  end subroutine center_zed_segment_complex
  
  subroutine finish_parallel_streaming

    use run_parameters, only: stream_implicit

    implicit none

    if (allocated(stream)) deallocate (stream)
    if (allocated(stream_c)) deallocate (stream_c)
    if (allocated(stream_sign)) deallocate (stream_sign)
    if (allocated(gradpar_c)) deallocate (gradpar_c)

    if (stream_implicit) call finish_invert_stream_operator

    parallel_streaming_initialized = .false.

  end subroutine finish_parallel_streaming

  subroutine finish_invert_stream_operator
    
    implicit none
    
    if (allocated(stream_tri_a1)) then
       deallocate (stream_tri_a1)
       deallocate (stream_tri_a2)
       deallocate (stream_tri_b1)
       deallocate (stream_tri_b2)
       deallocate (stream_tri_c1)
       deallocate (stream_tri_c2)
    end if
    
  end subroutine finish_invert_stream_operator

end module parallel_streaming
