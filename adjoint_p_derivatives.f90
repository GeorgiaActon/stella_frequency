module adjoint_p_derivatives
  
  implicit none

  public :: perturb_p
  public :: allocate_pert
  public :: allocate_unpert
  public :: deallocate_pert
  public :: deallocate_unpert
  public :: calculate_chi_lam
  public :: get_denominator
  private
  
contains 

  subroutine perturb_p! (term_out)
    
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid,ntubes
    use kt_grids, only: naky,nakx
    
    use physics_flags, only: include_parallel_streaming
    use physics_flags, only: include_mirror

    implicit none

    complex, dimension (:,:,:,:,:), allocatable :: g_term
    complex, dimension (:,:,:,:), allocatable :: q_term
!    complex, dimension (:,:), allocatable, intent (out) :: term_out
!    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
!    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: qout

    integer :: iv, imu, is,iz,ivmu, it

    allocate(g_term(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))    
    g_term = 0.
    allocate(q_term(naky,nakx,-nzgrid:nzgrid,ntubes))
    q_term = 0.

!    if (include_mirror) call get_mirror_term (g_term)
!    if (include_parallel_streaming) call get_prll_strm_term (g_term)
    
!    call get_wstar_term (g_term)
!    call get_wdrift_term (g_term)
    call get_omega_term (g_term)

    call get_quasin_term (q_term)

    open(15, file="perturb_p_gnew.txt", status="unknown",action="write",position="append")
    write(15,*) g_term(1,1,1,1,1), q_term(1,1,1,1)
    close(15)


 !   call lagrangian_integrals (g_term, q_term, term_out)

    deallocate(g_term)
    deallocate(q_term)
    
  end subroutine perturb_p

  subroutine get_denominator

    use adjoint_distfn_arrays, only: lam_save, gsave
    use adjoint_field_arrays, only: denominator
    use vpamu_grids, only: integrate_species
    use species, only: spec, nspec
    use kt_grids, only: naky,nakx
    use zgrid, only: nzgrid,ntubes
    use stella_layouts, only: vmu_lo
    use volume_averages, only: fieldline_average

    use mp, only: sum_allreduce
    implicit none

    complex, dimension (:,:,:,:), allocatable :: int_tubes
    real, dimension (:), allocatable :: wgts
    integer :: iz, it 

    allocate(int_tubes(naky,nakx,-nzgrid:nzgrid,ntubes))
    int_tubes = 0.
    allocate(wgts(nspec))
    wgts = 1.

    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          call integrate_species(gsave(:,:,iz,it,:)*lam_save(:,:,iz,it,:),&
               iz,wgts,int_tubes(:,:,iz,it),reduce_in=.false.)
       end do
    end do
    call sum_allreduce(int_tubes)
    
    call fieldline_average(int_tubes, denominator)
    
    deallocate(int_tubes)
    deallocate(wgts)
    
  end subroutine get_denominator

  subroutine lagrangian_integrals (gin, qin, int_term)
    use zgrid, only: nzgrid,ntubes
    use kt_grids, only: naky,nakx
    use species, only: spec, nspec
    use vpamu_grids, only: integrate_species

    use adjoint_distfn_arrays, only: pert_l
    use adjoint_field_arrays, only: pert_q, pert_term
    use adjoint_distfn_arrays, only: lam_save
    use adjoint_field_arrays, only: chi_save

    use stella_layouts, only: vmu_lo
    use volume_averages, only: fieldline_average

    use mp, only: sum_allreduce
    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gin
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: qin
    complex, dimension (:,:), allocatable, intent (out) :: int_term
    
    complex, dimension (:,:), allocatable :: qout, gout
    complex, dimension (:,:,:,:), allocatable :: int_tubes
    real, dimension (:), allocatable :: wgts
    integer :: it, iz,ivmu, iv,imu,is

    allocate(int_term(naky,nakx))
    int_term = 0.
    allocate(gout(naky,nakx))
    gout = 0.
    allocate(int_tubes(naky,nakx,-nzgrid:nzgrid,ntubes))
    int_tubes = 0.
    allocate(qout(naky,nakx))
    qout = 0.
    allocate(wgts(nspec))
    wgts = 1.

    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          call integrate_species(gin(:,:,iz,it,:)*lam_save(:,:,iz,it,:),&
               iz,wgts,int_tubes(:,:,iz,it),reduce_in=.false.)
       end do
    end do

    call sum_allreduce(int_tubes)
    
    call fieldline_average(int_tubes, gout)
    call fieldline_average(qin*chi_save, qout)

    int_term = gout + qout

    deallocate(gout)
    deallocate(int_tubes)
    deallocate(qout)
    deallocate(wgts)
    
  end subroutine lagrangian_integrals

  subroutine calculate_chi_lam
    use zgrid, only: nzgrid,ntubes
    use kt_grids, only: naky,nakx
    use vpamu_grids, only: nvpa, nmu
    use stella_layouts, only: kxkyz_lo
    use stella_geometry, only: gradpar, jacob
    use species, only: spec
    use redistribute, only: scatter, gather
    use dist_redistribute, only: kxkyz2vmu
    use adjoint_distfn_arrays, only: lam_save
    use adjoint_field_arrays, only: chi_save
    use stella_geometry, only: bmag, bmag_psi0
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx
    use stella_time, only: code_time
    use adjoint_field_arrays, only: omega
    
    implicit none

    complex, dimension (:,:,:), allocatable :: lamv, lamvmu
    
    integer :: ia, iz, ivmu, iv, imu

    allocate(lamv(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc))
    allocate(lamvmu(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_proc))
    
    call scatter (kxkyz2vmu, lam_save, lamv)
    
    do iv = 1, nvpa
       lamvmu(iv,:,:) = lamv(nvpa-iv+1,:,:)
    end do
    ia = 1.
    call gather (kxkyz2vmu, lamvmu, lam_save)
    do iz = -nzgrid,nzgrid
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          lam_save(:,:,iz,ia,ivmu) = lam_save(:,:,iz,ia,ivmu)*minval(bmag_psi0)/(gradpar(iz)*jacob(ia,iz)*&
               bmag(ia,iz))!*(1-omega(:,:)*code_time)*exp(omega(:,:)*code_time))
       end do
       chi_save(:,:,iz,ia) = chi_save(:,:,iz,ia)*minval(bmag_psi0)/(gradpar(iz)*jacob(ia,iz)&
            *bmag(ia,iz))!*(1-omega(:,:)*code_time)*exp(omega(:,:)*code_time))
    end do

    deallocate(lamv)
    deallocate(lamvmu)

  end subroutine calculate_chi_lam

  !! Mirror Term
  subroutine get_mirror_term (gout)
    
    use vpamu_grids, only: nmu, mu
    use zgrid, only: nzgrid,ntubes
    use species, only: spec,nspec
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx
    
    use stella_geometry, only: bmag, gradpar, dbdzed
    use kt_grids, only: naky,nakx

    use physics_flags, only: full_flux_surface
    use mirror_terms, only: get_dgdvpa_explicit, get_dgdvpa_global
    use dist_fn_arrays, only: gvmu
    use adjoint_distfn_arrays, only: gsave
    
    use redistribute, only: gather, scatter
    use dist_redistribute, only: kxkyz2vmu
    
    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    
    real, dimension (:,:,:,:), allocatable :: mirror
    complex, dimension (:,:,:,:,:), allocatable :: dgdvpa
    integer :: iz, imu, ia
    
    allocate(mirror(ntubes,-nzgrid:nzgrid,nmu,nspec))
    allocate(dgdvpa(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) 
    
    ia = 1
    
    call scatter (kxkyz2vmu, gsave, gvmu)
    if (full_flux_surface) then
       call get_dgdvpa_global (gvmu)
    else
      call get_dgdvpa_explicit (gvmu)
   end if
   call gather (kxkyz2vmu, gvmu, dgdvpa)

   mirror = 0.
    do iz = -nzgrid, nzgrid
       do imu =	1, nmu
          mirror(ia,iz,imu,:) = - spec%stm_psi0*gradpar(iz)*mu(imu)*dbdzed(ia,iz)
       end do
    end do

    call add_mirror_term (dgdvpa, mirror, gout)
    
    deallocate(mirror)
    deallocate(dgdvpa)
  end subroutine get_mirror_term

  subroutine add_mirror_term (g, mirror, mirr_out)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, is_idx
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension(:,-nzgrid:,:,:), intent (in) :: mirror
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: mirr_out
    
    integer :: imu, is, ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       mirr_out(:,:,:,:,ivmu) = mirr_out(:,:,:,:,ivmu) + &
            spread(spread(spread(mirror(1,:,imu,is),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_mirror_term

  !! Parallel streaming term !!
  subroutine get_prll_strm_term (gout)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use vpamu_grids, only: nvpa, nmu, vpa
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid, ntubes, nztot
    use species, only: spec, nspec
    
    use stella_geometry, only: gradpar
    use vpamu_grids, only:  maxwell_vpa, maxwell_mu, maxwell_fac

    use gyro_averages, only: gyro_average
    use parallel_streaming, only: get_dzed
    
    use adjoint_distfn_arrays, only: gsave
    use adjoint_field_arrays, only: phi_save
    
    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    
    real, dimension(:,:,:), allocatable :: stream
    complex, dimension (:,:,:,:), allocatable :: par    
    complex, dimension (:,:,:,:), allocatable :: gyrophi
    complex, dimension (:,:,:,:), allocatable :: dphidz, dgdz    
    integer :: is, imu, iv, ia, ivmu

    allocate(stream(-nzgrid:nzgrid,nvpa,nspec))

    allocate(gyrophi(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate(dphidz(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate(dgdz(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate(par(naky,nakx,-nzgrid:nzgrid,ntubes))
    
    ia = 1

    stream = 0.
    stream = spread(spread(spec%stm_psi0,1,nztot),2,nvpa) &
         *spread(spread(vpa,1,nztot)*spread(gradpar,2,nvpa),3,nspec)
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc       
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       gyrophi = 0.
       call gyro_average (phi_save, ivmu, gyrophi)

       dphidz = 0.
       dgdz = 0.
       
       call get_dzed (iv, gyrophi, dphidz)
       call get_dzed (iv, gsave(:,:,:,:,ivmu), dgdz)
       
       par = 0.
       par(:,:,:,ia) = dgdz(:,:,:,ia) + &
            spec(is)%zt*maxwell_vpa(iv,is)*spread(spread(maxwell_mu(ia,:,imu,is),1,naky),2,nakx)&
            *maxwell_fac(is)*dphidz(:,:,:,ia)

       call add_stream_term (par, stream, ivmu, gout(:,:,:,:,ivmu))
    end do
    
    deallocate(stream)
    deallocate(gyrophi)
    deallocate(dphidz)
    deallocate(dgdz)
    deallocate(par)
    
  end subroutine get_prll_strm_term
  
  subroutine add_stream_term (g, stream, ivmu, stream_out)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx                                                                                                     

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    real, dimension (:,:,:), intent (in) :: stream
    integer, intent (in) :: ivmu

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: stream_out

    integer :: iv, is

    iv = iv_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)
    
    stream_out(:,:,:,:) = stream_out(:,:,:,:) + &
         spread(spread(spread(stream(:,iv,is),1,naky),2,nakx),4,ntubes)*g(:,:,:,:)

  end subroutine add_stream_term


  !! Diamagnetic drift term !!
  subroutine get_wstar_term (gout)

    use zgrid, only: nzgrid, ntubes, nztot
    use kt_grids, only: naky, nakx
    use stella_layouts, only: vmu_lo

    use dist_fn_arrays, only: wstar
    use fields_arrays, only: apar
    use adjoint_field_arrays, only: phi_save
    use stella_time, only: code_dt

    use fields, only: get_dchidy

    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use kt_grids, only: aky
    use constants, only: zi
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use species, only: spec
    use adjoint_field_arrays, only: omega_g
    use adjoint_distfn_arrays, only: gsave

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    complex, dimension (:,:,:,:,:), allocatable :: dchidy
    real, dimension (:,:,:), allocatable :: star
    complex, dimension (:,:,:,:,:), allocatable :: check
    
    integer :: ivmu,is,imu,iv
    allocate(dchidy(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(check(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(star(ntubes,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    
    check = 0.
    dchidy = 0.
    star = -wstar/code_dt

    call get_dchidy (phi_save, apar, dchidy)
    apar = 0.
    call add_wstar_term (dchidy, star, gout)

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       check(:,:,:,:,ivmu) = spread(spread(omega_g, 3, nztot),4, ntubes)*gsave(:,:,:,:,ivmu) + &
            gout(1,1,1,1,ivmu)
    end do

    deallocate(star)
    deallocate(check)
    deallocate(dchidy)
  end subroutine get_wstar_term

  subroutine add_wstar_term (g, wstar, wstar_out)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none
    
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension(:,-nzgrid:,vmu_lo%llim_proc:), intent(in) :: wstar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: wstar_out
    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
       wstar_out(:,:,:,:,ivmu) = wstar_out(:,:,:,:,ivmu) + &
            spread(spread(spread(wstar(1,:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_wstar_term
    
!!!! Fluid Drifts !!!!
  
  subroutine get_wdrift_term (gout)

    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use adjoint_distfn_arrays, only: gsave
    use adjoint_field_arrays, only: phi_save
    
    use gyro_averages, only: gyro_average
    
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use stella_time, only: code_dt
    use species, only: nspec, spec
    
    use time_advance, only: get_dgdy, get_dgdx
    use constants, only: zi

    implicit none
    
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    
    complex, dimension(:,:,:,:,:), allocatable :: gx, gy
    complex, dimension(:,:,:,:), allocatable :: dphidx, dphidy

    integer :: ivmu
    
    allocate(dphidx(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate(dphidy(naky,nakx,-nzgrid:nzgrid,ntubes))
    
    allocate(gx(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate(gy(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    
    wdriftx_g = -wdriftx_g/code_dt
    wdrifty_g = -wdrifty_g/code_dt
    wdriftx_phi = -wdriftx_phi/code_dt
    wdrifty_phi = -wdrifty_phi/code_dt

    !! x-derivatives
    gx = 0.
    call get_dgdx (gsave, gx)
    !! y-derivatives
    gy = 0.
    call get_dgdy (gsave, gy)
    
    call add_drift (gx, wdriftx_g(1,:,:), gout)
    call add_drift (gy, wdrifty_g(1,:,:), gout)

    !! phi x-derivatives
    call get_dgdx (phi_save, dphidx)
    !! phi y-derivatives
    call get_dgdy (phi_save, dphidy)

    gx = 0.
    gy = 0.
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
       call gyro_average (dphidx, ivmu, gx(:,:,:,:,ivmu))
       call gyro_average (dphidy, ivmu, gy(:,:,:,:,ivmu))
    end do

    call add_drift (gx, wdriftx_phi(1,:,:), gout)
    call add_drift (gy, wdrifty_phi(1,:,:), gout)

    deallocate(dphidx, dphidy)
    deallocate(gx, gy)

  end subroutine get_wdrift_term


  !! Add component of drift to drift terms !!
  subroutine add_drift (g, wdrift, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) + &
            spread(spread(spread(wdrift(:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_drift


  subroutine get_omega_term (gout)
    use adjoint_distfn_arrays, only: gsave
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes, nztot
    use adjoint_field_arrays, only: omega_g

    use gyro_averages, only: gyro_average
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use adjoint_field_arrays, only: phi_save
    use adjoint_distfn_arrays, only: g_omega
    use stella_time, only: code_dt
    use species, only: spec
    use kt_grids, only: naky, nakx
    use dist_fn_arrays, only: gnew
    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    integer :: ivmu

    complex, dimension (:,:,:,:,:), allocatable :: gyro
    integer :: iv,imu,is, ia
    
    allocate(gyro(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gyro = 0.
    ia = 1 
    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       call gyro_average (phi_save, ivmu, gyro(:,:,:,:,ivmu))
       gout(:,:,:,:,ivmu) = gout(:,:,:,:,ivmu) &
            + spread(spread(omega_g, 3, nztot),4, ntubes)* &
            (gsave(:,:,:,:,ivmu) - spec(is)%zt*maxwell_fac(is)*maxwell_vpa(iv,is) & 
            *spread(spread(spread(maxwell_mu(ia,:,imu,is),1,naky),2,nakx),4,ntubes)&
            * gyro(:,:,:,:,ivmu))
    end do

    ivmu = 1 
    iv = iv_idx(vmu_lo,ivmu)
    imu = imu_idx(vmu_lo,ivmu)
    is = is_idx(vmu_lo,ivmu)
    open(10, file="omega_g2.txt", status="unknown",action="write",position="append")
    write(10,*) omega_g(1,1)*gsave(1,1,1,1,1), (gnew(1,1,1,1,1)-g_omega(1,1,1,1,1))/code_dt, &
         spec(is)%zt*omega_g(1,1)*maxwell_fac(is)*maxwell_vpa(iv,is)*maxwell_mu(1,1,imu,is)*&
         gyro(1,1,1,1,1),&
         (gnew(1,1,1,1,1)-g_omega(1,1,1,1,1))/code_dt - omega_g(1,1)*gsave(1,1,1,1,1)
    close(10)
    
    deallocate(gyro)         

  end subroutine get_omega_term
  
  !! Quasineutrality term!!

  subroutine get_quasin_term (qout)
    
!    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
!    use kt_grids, only: naky, nakx

    use adjoint_distfn_arrays, only: gsave
    use adjoint_field_arrays, only: phi_save
    use fields_arrays, only: gamtot, phi, apar
!    use species, only: spec
    !use gyro_averages, only: gyro_average
    !use vpamu_grids, only: integrate_species
    
    use mp, only: sum_allreduce
    use fields, only: advance_fields, fields_updated

    implicit none
    
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: qout
!    complex, dimension (:,:,:), allocatable :: gyro_g
    integer :: ivmu, iz, it
    
    qout = 0.
    ! allocate (gyro_g(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! gyro_g = 0.
    ! do it = 1, ntubes
    !    do iz = -nzgrid, nzgrid
    !       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_alloc
    !          call gyro_average (gsave(:,:,iz,it,ivmu),iz, ivmu, gyro_g(:,:,ivmu))
    !       end do
    !       call integrate_species (gyro_g,iz,spec%z*spec%dens_psi0, qout(:,:,iz,it),reduce_in=.false.)
    !    end do
    ! end do
    ! call sum_allreduce(qout)
    
    ! deallocate(gyro_g)
    
    fields_updated = .false.
    call advance_fields (gsave, phi, apar, dist='gbar')

    qout = (phi - phi_save)*spread(gamtot,4,ntubes)
    
!    qout = qout - spread(gamtot,4,ntubes)*phi_save

!    deallocate (gyro_g)
  end subroutine get_quasin_term

  subroutine integrate_zed (f, intf)

    use zgrid, only: nzgrid, delzed, ntubes
    
    implicit none
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: f
    complex, dimension (:,:), intent (out) :: intf

    integer :: iz

    intf = 0.
    do iz = -nzgrid, nzgrid
       intf(:,:) = intf(:,:) + sum(f(:,:,iz,:)*delzed(iz), dim = 3)
    end do
    
    intf = intf/real(ntubes)
!    intf = 0.5*delzed(0)*(f(:,:,-nzgrid) + f(:,:,nzgrid) + intf)

  end subroutine integrate_zed

  !! Allocate and Deallocate !!
  subroutine allocate_unpert (np)
    
    use stella_layouts, only: vmu_lo
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid, ntubes

 !   use adjoint_field_arrays, only: unpert_q, denominator
  !  use adjoint_distfn_arrays, only: unpert_l

    use adjoint_field_arrays, only: derivative, denominator
    use adjoint_field_arrays, only: unpert_term
    implicit none

    integer, intent (in) :: np

    ! if(.not. allocated(unpert_q)) then
    !    allocate(unpert_q(naky,nakx,-nzgrid:nzgrid,ntubes))
    !    unpert_q = 0.
    ! end if
    ! if(.not. allocated(unpert_l)) then
    !    allocate(unpert_l(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    !    unpert_l = 0.
    ! end if
    if(.not. allocated(unpert_term)) then
       allocate(unpert_term(naky,nakx))
       unpert_term = 0.
    end if
    
    if(.not. allocated(denominator)) then
       allocate(denominator(naky,nakx))
       denominator = 0.
    end if
    
    if(.not. allocated(derivative)) then
       allocate(derivative(naky,nakx,np))
       derivative = 0.
    end if
 !   call allocate_chi_lam
    
  end subroutine allocate_unpert

  subroutine allocate_pert
    use stella_layouts, only: vmu_lo
    use kt_grids, only: nakx, naky
    use zgrid, only: nzgrid, ntubes
    
!    use adjoint_distfn_arrays, only: pert_l
    use adjoint_field_arrays, only: pert_term

    implicit none
    
    ! if(.not. allocated(pert_q)) then
    !    allocate(pert_q(naky,nakx,-nzgrid:nzgrid,ntubes))
    !    pert_q = 0.
    ! end if
    ! if(.not. allocated(pert_l)) then
    !    allocate(pert_l(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    !    pert_l =	0.
    ! end if
    if(.not. allocated(pert_term)) then
       allocate(pert_term(naky,nakx))
       pert_term = 0.
    end if

  end subroutine allocate_pert
  
  ! subroutine allocate_chi_lam
  !   use stella_layouts, only: vmu_lo
  !   use kt_grids, only: nakx, naky
  !   use zgrid, only: nzgrid, ntubes
  !   use adjoint_distfn_arrays, only: lam_save
  !   use adjoint_field_arrays, only: chi_save

  !   implicit none

  !   if(.not. allocated(lam_save)) then
  !      allocate(lam_save(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
  !      lam_save = 0.
  !   end if
  !   if(.not. allocated(chi_save)) then
  !      allocate(chi_save(naky,nakx,-nzgrid:nzgrid,ntubes))
  !      chi_save = 0.
  !   end if

  ! end subroutine allocate_chi_lam
  
  subroutine deallocate_pert

!    use adjoint_distfn_arrays, only: pert_l
!    use adjoint_field_arrays, only: pert_q, pert_term, chi_save
    use adjoint_field_arrays, only: pert_term
    implicit none

!    deallocate(pert_q)
!    deallocate(pert_l)
    deallocate(pert_term)
    
  end subroutine deallocate_pert
  
  subroutine deallocate_unpert

!    use adjoint_field_arrays, only: unpert_q
!    use adjoint_distfn_arrays, only: unpert_l
    use adjoint_field_arrays, only: unpert_term
    use adjoint_field_arrays, only: derivative, denominator
    
    implicit none

!    deallocate(unpert_q)
!    deallocate(unpert_l)
    deallocate(unpert_term)
    deallocate(denominator)
    deallocate(derivative)

!    call deallocate_chi_lam
    
  end subroutine deallocate_unpert

  ! subroutine deallocate_chi_lam

  !   use adjoint_distfn_arrays, only: lam_save
  !   use adjoint_field_arrays, only: chi_save

  !   implicit none

  !   deallocate(lam_save)
  !   deallocate(chi_save)
    
  ! end subroutine deallocate_chi_lam
end module adjoint_p_derivatives
