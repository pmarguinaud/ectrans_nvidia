! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program transform_test

#ifdef USE_GPU
#define PINNED_TAG , pinned
#else
#define PINNED_TAG
#endif

!
! Spectral transform test
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
! 1) One "surface" field is always transformed:
!      zspsc2(1,1:nspec2) <-> zgmvs(1:nproma,1:1,1:ngbplk)
!
! 2) A Multiple "3d" fields are transformed and can be disabled with "--nfld 0"
!
!      zspsc3a(1:nlev,1:nspec2,1:nfld) <-> zgp3a(1:nproma,1:nlev,1:nfld,1:ngpblk)
!
! 3) Optionally a "3d" vorticity/divergence field is transformed to uv (wind) and
!   can be enabled with "--vordiv"
!
!      zspvor(1:nlev,1:nspec2) / zspdiv(1:nlev,1:nspec2) <-> zgpuv(1:nproma,1:nlev,1:2,1:ngpblk)
!
! 4) Optionally scalar derivatives can be computed for the fields described in 1) and 2)
!    This must be enabled with "--scders"
!
! 5) Optionally uv East-West derivate can be computed from vorticity/divergence.
!    This must be enabled with "--vordiv --uvders"
!
!
! Authors : George Mozdzynski
!           Willem Deconinck
!           Ioan Hadade
!           Sam Hatfield
!

use parkind1, only: jpim, jprb, jprd
use oml_mod ,only : oml_max_threads
use mpl_module
use yomgstats, only: jpmaxstat, gstats_lstats => lstats
use yomhook, only : dr_hook_init

implicit none

! Number of points in top/bottom latitudes
integer(kind=jpim), parameter :: min_octa_points = 20

integer(kind=jpim) :: istack, getstackusage
real(kind=jprb), dimension(1) :: zmaxerr(5), zerr(5)
real(kind=jprb) :: zmaxerrg

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output

! Default parameters
integer(kind=jpim) :: nsmax   = 79  ! Spectral truncation
integer(kind=jpim) :: iters   = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of scalar fields 
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels

integer(kind=jpim) :: nflevg
integer(kind=jpim) :: ndgl ! Number of latitudes
integer(kind=jpim) :: nspec2
integer(kind=jpim) :: ngptot
integer(kind=jpim) :: ngptotg
integer(kind=jpim) :: ifld
integer(kind=jpim) :: jroc
integer(kind=jpim) :: jb
integer(kind=jpim) :: nspec2g
integer(kind=jpim) :: i
integer(kind=jpim) :: ja
integer(kind=jpim) :: ib
integer(kind=jpim) :: jprtrv

integer(kind=jpim), allocatable :: nloen(:), nprcids(:)
integer(kind=jpim) :: myproc, jj
integer :: jstep

real(kind=jprd) :: ztinit, ztloop, timef, ztstepmax, ztstepmin, ztstepavg, ztstepmed
real(kind=jprd) :: ztstepmax1, ztstepmin1, ztstepavg1, ztstepmed1
real(kind=jprd) :: ztstepmax2, ztstepmin2, ztstepavg2, ztstepmed2
real(kind=jprd), allocatable :: ztstep(:), ztstep1(:), ztstep2(:)

real(kind=jprb), allocatable :: znormsp(:), znormsp1(:), znormdiv(:), znormdiv1(:)
real(kind=jprb), allocatable :: znormvor(:), znormvor1(:), znormt(:), znormt1(:)
real(kind=jprd) :: zaveave(0:jpmaxstat)

! Grid-point space data structures
real(kind=jprb), allocatable, target PINNED_TAG :: zgmv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), allocatable, target PINNED_TAG :: zgmvs  (:,:,:)   ! Single level fields at t and t-dt
real(kind=jprb), pointer :: zgp3a (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgpuv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgp2 (:,:,:) ! Single level fields at t and t-dt

! Spectral space data structures
real(kind=jprb), allocatable, target PINNED_TAG :: sp3d(:,:,:)
real(kind=jprb), pointer :: zspvor(:,:) => null()
real(kind=jprb), pointer :: zspdiv(:,:) => null()
real(kind=jprb), pointer :: zspsc3a(:,:,:) => null()
real(kind=jprb), allocatable PINNED_TAG :: zspsc2(:,:)

logical :: lstack = .false. ! Output stack info
logical :: luserpnm = .false.
logical :: lkeeprpnm = .false.
logical :: luseflt = .false. ! Use fast legendre transforms
logical :: ltrace_stats = .false.
logical :: lstats_omp = .false.
logical :: lstats_comms = .false.
logical :: lstats_mpl = .false.
logical :: lstats = .true. ! gstats statistics
logical :: lbarrier_stats = .false.
logical :: lbarrier_stats2 = .false.
logical :: ldetailed_stats = .false.
logical :: lstats_alloc = .false.
logical :: lsyncstats = .false.
logical :: lstatscpu = .false.
logical :: lstats_mem = .false.
logical :: lxml_stats = .false.
logical :: lfftw = .true. ! Use FFTW for Fourier transforms
logical :: lvordiv = .false.
logical :: lscders = .false.
logical :: luvders = .false.
logical :: lprint_norms = .false. ! Calculate and print spectral norms
logical :: lmeminfo = .false. ! Show information from FIAT routine ec_meminfo at the end

integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1

! The multiplier of the machine epsilon used as a tolerance for correctness checking
! ncheck = 0 (the default) means that correctness checking is disabled
integer(kind=jpim) :: ncheck = 0

logical :: lmpoff = .false. ! Message passing switch

! Verbosity level (0 or 1)
integer :: verbosity = 0

real(kind=jprb) :: zra = 6371229._jprb

integer(kind=jpim) :: nmax_resol = 37 ! Max number of resolutions
integer(kind=jpim) :: npromatr = 0 ! nproma for trans lib
integer(kind=jpim) :: ncombflen = 1800000 ! Size of comm buffer

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: nspecresmin = 80 ! Minimum spectral resolution, for controlling nprtrw
integer(kind=jpim) :: mysetv
integer(kind=jpim) :: mysetw
integer(kind=jpim) :: mp_type = 2 ! Message passing type
integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size

integer(kind=jpim), allocatable :: numll(:), ivset(:)
integer(kind=jpim) :: ivsetsc(1)

integer(kind=jpim) :: nflevl

! sumpini
integer(kind=jpim) :: isqr
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag


integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev

integer(kind=jpim) :: ndimgmv  = 0 ! Third dim. of gmv "(nproma,nflevg,ndimgmv,ngpblks)"
integer(kind=jpim) :: ndimgmvs = 0 ! Second dim. gmvs "(nproma,ndimgmvs,ngpblks)"

integer(kind=jpim) :: jbegin_uv = 0
integer(kind=jpim) :: jend_uv   = 0
integer(kind=jpim) :: jbegin_sc = 0
integer(kind=jpim) :: jend_sc   = 0
integer(kind=jpim) :: jbegin_scder_NS = 0
integer(kind=jpim) :: jend_scder_NS = 0
integer(kind=jpim) :: jbegin_scder_EW = 0
integer(kind=jpim) :: jend_scder_EW = 0
integer(kind=jpim) :: jbegin_uder_EW = 0
integer(kind=jpim) :: jend_uder_EW = 0
integer(kind=jpim) :: jbegin_vder_EW = 0
integer(kind=jpim) :: jend_vder_EW = 0

logical :: ldump_values = .false.

integer, external :: ec_mpirank
logical :: luse_mpi = .true.

character(len=16) :: cgrid = ''

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"
#include "ec_meminfo.intfb.h"

!===================================================================================================
stop

! Setup
nflevg = nlev

!===================================================================================================

if (luse_mpi) then
  call mpl_init(ldinfo=(verbosity>=1))
  nproc  = mpl_nproc()
  myproc = mpl_myrank()
else
  nproc = 1
  myproc = 1
  mpl_comm = -1
endif
nthread = oml_max_threads()

call dr_hook_init()

!===================================================================================================

if( lstats ) call gstats(0,0)
ztinit = timef()

! only output to stdout on pe 1
if (nproc > 1) then
  if (myproc /= 1) then
    open(unit=nout, file='/dev/null')
  endif
endif

if (ldetailed_stats) then
  lstats_omp    = .true.
  lstats_comms  = .true.
  lstats_mpl    = .true.
  lstatscpu     = .true.
  nprnt_stats   = nproc
!  lstats_mem   = .true.
!  lstats_alloc = .true.
endif


!===================================================================================================

allocate(nprcids(nproc))
do jj = 1, nproc
  nprcids(jj) = jj
enddo

if (nproc <= 1) then
  lmpoff = .true.
endif

! Compute nprgpns and nprgpew
! This version selects most square-like distribution
! These will change if leq_regions=.true.
if (nproc == 0) nproc = 1
isqr = int(sqrt(real(nproc,jprb)))
do ja = isqr, nproc
  ib = nproc/ja
  if (ja*ib == nproc) then
    nprgpns = max(ja,ib)
    nprgpew = min(ja,ib)
    exit
  endif
enddo

! From sumpini, although this should be specified in namelist
if (nspecresmin == 0) nspecresmin = nproc

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv > 0 .or. nprtrw > 0) then
  if (nprtrv == 0) nprtrv = nproc/nprtrw
  if (nprtrw == 0) nprtrw = nproc/nprtrv
  if (nprtrw*nprtrv /= nproc) call abor1('transform_test:nprtrw*nprtrv /= nproc')
  if (nprtrw > nspecresmin) call abor1('transform_test:nprtrw > nspecresmin')
else
  do jprtrv = 4, nproc
    nprtrv = jprtrv
    nprtrw = nproc/nprtrv
    if (nprtrv*nprtrw /= nproc) cycle
    if (nprtrv > nprtrw) exit
    if (nprtrw > nspecresmin) cycle
    if (nprtrw <= nspecresmin/(2*oml_max_threads())) exit
  enddo
  ! Go for approx square partition for backup
  if (nprtrv*nprtrw /= nproc .or. nprtrw > nspecresmin .or. nprtrv > nprtrw) then
    isqr = int(sqrt(real(nproc,jprb)))
    do ja = isqr, nproc
      ib = nproc/ja
      if (ja*ib == nproc) then
        nprtrw = max(ja, ib)
        nprtrv = min(ja, ib)
        if (nprtrw > nspecresmin ) then
          call abor1('transform_test:nprtrw (approx square value) > nspecresmin')
        endif
        exit
      endif
    enddo
  endif
endif

! Create communicators for mpi groups
if (.not.lmpoff) then
  call mpl_groups_create(nprtrw, nprtrv)
endif

if (lmpoff) then
  mysetw = (myproc - 1)/nprtrv + 1
  mysetv = mod(myproc - 1, nprtrv) + 1
else
  call mpl_cart_coords(myproc, mysetw, mysetv)

  ! Just checking for now...
  iprtrv = mod(myproc - 1, nprtrv) + 1
  iprtrw = (myproc - 1)/nprtrv + 1
  if (iprtrv /= mysetv .or. iprtrw /= mysetw) then
    call abor1('transform_test:inconsistency when computing mysetw and mysetv')
  endif
endif

if (.not. lmpoff) then
  call mpl_buffer_method(kmp_type=mp_type, kmbx_size=mbx_size, kprocids=nprcids, ldinfo=(verbosity>=1))
endif

! Determine number of local levels for fourier and legendre calculations
! based on the values of nflevg and nprtrv
allocate(numll(nprtrv+1))

! Calculate remainder
iprused = min(nflevg+1, nprtrv)
ilevpp = nflevg/nprtrv
irest = nflevg -ilevpp*nprtrv
do jroc = 1, nprtrv
  if (jroc <= irest) then
    numll(jroc) = ilevpp+1
  else
    numll(jroc) = ilevpp
  endif
enddo
numll(iprused+1:nprtrv+1) = 0

nflevl = numll(mysetv)

ivsetsc(1) = iprused
ifld = 0

!===================================================================================================
! Setup gstats
!===================================================================================================

if (lstats) then
  call gstats_setup(nproc, myproc, nprcids,                                            &
    & lstats, lstatscpu, lsyncstats, ldetailed_stats, lbarrier_stats, lbarrier_stats2, &
    & lstats_omp, lstats_comms, lstats_mem, nstats_mem, lstats_alloc,                  &
    & ltrace_stats, ntrace_stats, nprnt_stats, lxml_stats)
  call gstats_psut

  ! TODO
  ! Assign labels to GSTATS regions
  !call gstats_labels
endif

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call gstats(1, 0)
call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1),                &
  &               kmax_resol=nmax_resol, kpromatr=npromatr, kprgpns=nprgpns, kprgpew=nprgpew, &
  &               kprtrw=nprtrw, kcombflen=ncombflen, ldsync_trans=lsync_trans,               &
  &               ldeq_regions=leq_regions, prad=zra, ldalloperm=.true., ldmpoff=.not.luse_mpi)
call gstats(1, 1)

call gstats(2, 0)
call setup_trans(ksmax=nsmax, kdgl=ndgl, kloen=nloen, ldsplit=.true.,          &
  &                 ldusefftw=lfftw, lduserpnm=luserpnm, ldkeeprpnm=lkeeprpnm, &
  &                 lduseflt=luseflt)
call gstats(2, 1)

call trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbosity >= 0) then
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'(" ")')
  write(nout,'("nsmax     ",i0)') nsmax
  write(nout,'("grid      ",a)') trim(cgrid)
  write(nout,'("ndgl      ",i0)') ndgl
  write(nout,'("nproc     ",i0)') nproc
  write(nout,'("nthread   ",i0)') nthread
  write(nout,'("nprgpns   ",i0)') nprgpns
  write(nout,'("nprgpew   ",i0)') nprgpew
  write(nout,'("nprtrw    ",i0)') nprtrw
  write(nout,'("nprtrv    ",i0)') nprtrv
  write(nout,'("ngptot    ",i0)') ngptot
  write(nout,'("ngptotg   ",i0)') ngptotg
  write(nout,'("nfld      ",i0)') nfld
  write(nout,'("nlev      ",i0)') nlev
  write(nout,'("nproma    ",i0)') nproma
  write(nout,'("ngpblks   ",i0)') ngpblks
  write(nout,'("nspec2    ",i0)') nspec2
  write(nout,'("nspec2g   ",i0)') nspec2g
  write(nout,'("luseflt   ",l)') luseflt
  write(nout,'("lvordiv   ",l)') lvordiv
  write(nout,'("lscders   ",l)') lscders
  write(nout,'("luvders   ",l)') luvders
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end if

!===================================================================================================
! Allocate and Initialize spectral arrays
!===================================================================================================

! Allocate spectral arrays
! Try to mimick IFS layout as much as possible
nullify(zspvor)
nullify(zspdiv)
nullify(zspsc3a)
allocate(sp3d(nflevl,nspec2,2+nfld))
allocate(zspsc2(1,nspec2))

! Point convenience variables to storage variable sp3d
zspvor  => sp3d(:,:,1)
zspdiv  => sp3d(:,:,2)
zspsc3a => sp3d(:,:,3:3+(nfld-1))

!===================================================================================================
! Allocate gridpoint arrays
!===================================================================================================

allocate(ivset(nflevg))

! Compute spectral distribution
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

! Allocate grid-point arrays
if (lvordiv) then
  jbegin_uv = 1
  jend_uv = 2
endif
if (luvders) then
  jbegin_uder_EW  = jend_uv + 1
  jend_uder_EW    = jbegin_uder_EW + 1
  jbegin_vder_EW  = jend_uder_EW + 1
  jend_vder_EW    = jbegin_vder_EW + 1
else
  jbegin_uder_EW = jend_uv
  jend_uder_EW   = jend_uv
  jbegin_vder_EW = jend_uv
  jend_vder_EW   = jend_uv
endif

jbegin_sc = jbegin_vder_EW + 1
jend_sc   = jbegin_vder_EW + nfld

if (lscders) then
  ndimgmvs = 3
  jbegin_scder_NS = jend_sc + 1
  jend_scder_NS   = jend_sc + nfld
  jbegin_scder_EW = jend_scder_NS + 1
  jend_scder_EW   = jend_scder_NS + nfld
else
  ndimgmvs = 1
  jbegin_scder_NS = jend_sc
  jend_scder_NS   = jend_sc
  jbegin_scder_EW = jend_sc
  jend_scder_EW   = jend_sc
endif

ndimgmv = jend_scder_EW

allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
allocate(zgmvs(nproma,ndimgmvs,ngpblks))

zgpuv => zgmv(:,:,1:jend_vder_EW,:)
zgp3a => zgmv(:,:,jbegin_sc:jend_scder_EW,:)
zgp2  => zgmvs(:,:,:)

!===================================================================================================
! Allocate norm arrays
!===================================================================================================

if (lprint_norms .or. ncheck > 0) then
  allocate(znormsp(1))
  allocate(znormsp1(1))
  allocate(znormvor(nflevg))
  allocate(znormvor1(nflevg))
  allocate(znormdiv(nflevg))
  allocate(znormdiv1(nflevg))
  allocate(znormt(nflevg))
  allocate(znormt1(nflevg))

  call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor1, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv1, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt1,   kvset=ivset(1:nflevg))
  call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp1,  kvset=ivsetsc)

  if (verbosity >= 1) then
    do ifld = 1, nflevg
      write(nout,'("norm zspvor( ",i4,",:)   = ",f20.15)') ifld, znormvor1(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv1(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormt1(ifld)
    enddo
    do ifld = 1, 1
      write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15)') ifld, znormsp1(ifld)
    enddo
  endif
endif

!===================================================================================================
! Setup timers
!===================================================================================================

ztinit = (timef() - ztinit)/1000.0_jprd

if (verbosity >= 0) then
  write(nout,'(" ")')
  write(nout,'(a,i6,a,f9.2,a)') "transform_test initialisation, on",nproc,&
                                & " tasks, took",ztinit," sec"
  write(nout,'(" ")')
endif

if (iters <= 0) call abor1('transform_test:iters <= 0')

allocate(ztstep(iters))
allocate(ztstep1(iters))
allocate(ztstep2(iters))

ztstepavg  = 0._jprd
ztstepmax  = 0._jprd
ztstepmin  = 9999999999999999._jprd
ztstepavg1 = 0._jprd
ztstepmax1 = 0._jprd
ztstepmin1 = 9999999999999999._jprd
ztstepavg2 = 0._jprd
ztstepmax2 = 0._jprd
ztstepmin2 = 9999999999999999._jprd

write(nout,'(a)') '======= Start of spectral transforms  ======='
write(nout,'(" ")')

ztloop = timef()

!===================================================================================================
! Do spectral transform loop
!===================================================================================================

gstats_lstats = .false.
do jstep = 1, iters+2
  if (jstep == 3) gstats_lstats = .true.

  call gstats(3,0)
  ztstep(jstep) = timef()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  ztstep1(jstep) = timef()
  call gstats(4,0)
    call inv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  & ! scalar derivatives
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgp3a=zgp3a)

  call gstats(4,1)

enddo


contains

end program transform_test

!===================================================================================================
