! PROGRAM TO ESTIMATE THE bmim ring center -- anion COM and anion atoms -- OPENMP parallel
! ALSO PRINT THE SDF OF THE SELECTED COORDINATION ENVIRONMENT

! NEEDS A REDESIGN: ARRANGE ACCORDING TO SOLVATION SHELL
! 1ST SOLVATION SHELL   : CL ,TFSI -- I. COM ; II. ATOMS  III. CHARGE || AND THE CORRESPONDING NORMALIZED DATA
! CONDITIONAL 1ST SHELL : CL, TFSI COM 
! CONDITIONAL 2ND SHELL : HA, HB, RING-CENTER

! FIRST STARTED ON 02 AUG 2021
! V2    STARTED ON 30 SEP 2021
! V3    STARTED ON 30 SEP 2021
! LAST MODIFIED ON 05 OCT 2021

! THE TOTAL CHARGE PRINTED SEEMED TO BE NORMALIZED TWICE! CORRECTED IN THIS VERSION

! THE CONDITIONAL MATRIX CONVENTION (TFSI,CL) COORD NUMs

! DEFINITION OF COORDINATE AXES:
! 1. RING CENTER OF GEOMETRY: CENTER
! 2. RING-CENTER--CA: Y-AXIS
! 3. RING-CENTER--NA(BUTYL SIDE): TEMP X-AXIS
! 4. X_TEMP .CROSS. Y : Z-AXIS
! 5. Y .CROSS. X : X-AXIS

! ------------------ CUBE FILE FORMAT --------------------------
! Introduction
! The cube file describes volumetric data as well as atom positions, it originates from the Gaussian software package. 
! The file consists of a header which includes the atom information and the size as well as orientation of the volumetric data. 
! This is followed by the volumetric data, one scalar per voxel element. All aspects of the file are text (human readable), 
! originally the numerical values were 5 wide for integers that started each header line (after the first) and 
! floating point values were formatted 12.6, that is, 12 characters wide with 6 decimal places.

! Header
! The first two lines of the header are comments, they are generally ignored by parsing packages or used as two default labels.
! The third line has the number of atoms included in the file followed by the position of the origin of the volumetric data.
! The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector.
! Note this means the volume need not be aligned with the coordinate axis, indeed it also means it may be sheared,
! although most volumetric packages won't support that. The length of each vector is the length of the side of the voxel thus allowing non cubic volumes.
! --> If the sign of the number of voxels in a dimension is positive then the units are Bohr, if negative then Angstroms.

! The last section in the header is one line for each atom consisting of 5 numbers, 
! the first is the atom number, second (?), the last three are the x,y,z coordinates of the atom center.

! Volumetric data
! The volumetric data is straightforward, one floating point number for each volumetric element. 
! The original Gaussian format arranged the values in the format shown below in the example, most parsing programs can read any white space separated format. 
! Traditionally the grid is arranged with the x axis as the outer loop and the z axis as the inner loop


MODULE vector_ops
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    INTERFACE OPERATOR (.cross.)
        PROCEDURE cross_prod
    END INTERFACE OPERATOR (.cross.)
    CONTAINS
    FUNCTION cross_prod(r1, r2)
        REAL(KIND=REAL32), INTENT(IN) :: r1(3), r2(3)
        REAL(KIND=REAL32)             :: cross_prod(3)
        cross_prod(1) = r1(2)*r2(3) - r1(3)*r2(2)
        cross_prod(2) = r1(3)*r2(1) - r1(1)*r2(3)
        cross_prod(3) = r1(1)*r2(2) - r1(2)*r2(1)
    END FUNCTION cross_prod
    SUBROUTINE v_normalize(r1)
        REAL(KIND=REAL32), INTENT(INOUT) :: r1(3)
        r1 = r1/SQRT(DOT_PRODUCT(r1,r1))
    END SUBROUTINE v_normalize
END MODULE vector_ops

PROGRAM sdf
    USE vector_ops
    USE gmxfort_trajectory
    USE ISO_FORTRAN_ENV 
    IMPLICIT NONE

    ! CONSTANTS
    INTEGER, PARAMETER :: PRINTEVERY    = 100
    REAL(KIND=REAL64), PARAMETER :: pi = 4.D0*ATAN(1.D0)
    REAL(KIND=REAL64), PARAMETER :: au_conv = 0.1889726125D+1 ! convert angstrom to atomic units

    ! TRAJECTORY RELATED
    TYPE(Trajectory)                 :: trj, test
    REAL(KIND=REAL32)                :: mybox(3,3), mytime
    INTEGER                          :: mystep, myatoms
    REAL(KIND=REAL32)                :: volume = 0.0
    INTEGER                          :: proc_frames=4, skip_frames=10
    
    TYPE :: coordtype
        REAL(KIND=REAL32), ALLOCATABLE :: r(:,:,:)
    END TYPE coordtype
    TYPE(coordtype), ALLOCATABLE       :: coord(:)

    ! I/O FILES RELATED
    INTEGER            :: inpunit, sdf_unit
    CHARACTER(LEN=500) :: inpfile_name, suffix_name, outfilename
    CHARACTER(LEN=500) :: xtcfile_name = '/home/nikhil/projects/&
                          8_extend_BIL/pre-may2021/welton-structures/mb-v1/bmim-cl-ntf2/&
                          0.50/3-nvt/nvt_mol.xtc'

    ! SDF RELATED VARIABLES 
    INTEGER, ALLOCATABLE           :: tfsi_com_hist(:,:,:), cl_com_hist(:,:,:)
    
    REAL(KIND=REAL32), ALLOCATABLE :: tfsi_com_dens(:,:,:), tfsi_com_norm_dens(:,:,:) 
    REAL(KIND=REAL32), ALLOCATABLE :: cl_com_dens(:,:,:), cl_com_norm_dens(:,:,:) 
    REAL(KIND=REAL32), ALLOCATABLE :: tfsi_charge_dens(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: cl_charge_dens(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: print_dens(:,:,:)

    ! CONDITIONAL DENSITY REALATED VARIABLES
    INTEGER, ALLOCATABLE           :: cond_tfsi_com_hist(:,:,:,:,:), cond_cl_com_hist(:,:,:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: cond_tfsi_com_dens(:,:,:,:,:), cond_cl_com_dens(:,:,:,:,:) 
    INTEGER, ALLOCATABLE           :: cond_bmim_cor_hist(:,:,:,:,:), cond_bmim_ha_hist(:,:,:,:,:), cond_bmim_hb_hist(:,:,:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: cond_bmim_cor_dens(:,:,:,:,:), cond_bmim_ha_dens(:,:,:,:,:), cond_bmim_hb_dens(:,:,:,:,:) 
    INTEGER, ALLOCATABLE           :: cond_hist_count(:,:)
    REAL(KIND=REAL32)              :: rcut_cl, rcut_tfsi, rcut_cl_sq, rcut_tfsi_sq
    INTEGER                        :: cond_maxcl_coordnum, cond_maxtfsi_coordnum, cond_num_cl_coordnum, cond_num_tfsi_coordnum

    INTEGER                        :: temp_cl_coordnum, temp_tfsi_coordnum, temp_cl_coordid(20), temp_tfsi_coordid(20)
    REAL(KIND=REAL32)              :: cond_norm_const

    ! TFSI ATOMIC SDF
    INTEGER                        :: num_N_atoms = 1, num_O_atoms = 4, num_F_atoms = 6
    INTEGER                        :: index_N_atoms(1), index_O_atoms(4), index_F_atoms(6)
    INTEGER, ALLOCATABLE           :: tfsi_N_hist(:,:,:),tfsi_O_hist(:,:,:),tfsi_F_hist(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: tfsi_N_dens(:,:,:),tfsi_O_dens(:,:,:),tfsi_F_dens(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: tfsi_N_norm_dens(:,:,:),tfsi_O_norm_dens(:,:,:),tfsi_F_norm_dens(:,:,:)

    ! INTER-BMI SDF
    INTEGER                        :: num_HA_atoms = 1, num_HB_atoms = 2
    INTEGER                        :: index_HA_atoms(1), index_HB_atoms(2)
    INTEGER                        :: tot_atoms_ha, tot_atoms_hb
    INTEGER, ALLOCATABLE           :: bmim_cor_hist(:,:,:), bmim_ha_hist(:,:,:), bmim_hb_hist(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: bmim_cor_dens(:,:,:), bmim_ha_dens(:,:,:), bmim_hb_dens(:,:,:)
    REAL(KIND=REAL32)              :: avg_ha_dens, avg_hb_dens


    REAL(KIND=REAL32), ALLOCATABLE :: vx_bmim_ring(:,:,:),vy_bmim_ring(:,:,:),vz_bmim_ring(:,:,:) 
    REAL(KIND=REAL32), ALLOCATABLE :: vx_ring_temp(:,:), vy_ring_temp(:,:), vz_ring_temp(:,:)
    
    REAL(KIND=REAL32), ALLOCATABLE :: box_size(:,:)
    REAL(KIND=REAL32)              :: rcut_x, rcut_y, rcut_z
    REAL(KIND=REAL32)              :: rbin_size=0.01, box_temp(3)
    INTEGER                        :: num_bin_x, num_bin_y, num_bin_z
    INTEGER                        :: max_bin_x, max_bin_y, max_bin_z
    LOGICAL                        :: select_rcut
    REAL(KIND=REAL32)              :: sphere_cut_off

    REAL(KIND=REAL32)              :: norm_const
    REAL(KIND=REAL32)              :: vox_vol
    REAL(KIND=REAL32)              :: delrv(3), delsq
    REAL(KIND=REAL32)              :: xproj, yproj, zproj
    INTEGER                        :: xbin, ybin, zbin
    REAL(KIND=REAL32)              :: vtemp(3)

    ! BMIM-CL-TFSI SYSTEM RELATED
    REAL(KIND=REAL32), ALLOCATABLE :: coord_bmim(:,:,:), coord_tfsi(:,:,:), coord_cl(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: coord_ha(:,:,:), coord_hb(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: coord_bmim_temp(:,:), coord_tfsi_temp(:,:), coord_cl_temp(:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: cor_bmim(:,:,:), com_tfsi(:,:,:), com_cl(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: cor_bmim_temp(:,:), com_tfsi_temp(:,:), com_cl_temp(:,:)
    REAL(KIND=REAL32)              :: cube_coord(3,25)
    INTEGER                        :: bmim_atomic_num(25)
    
    INTEGER           :: num_res_bmim, num_res_tfsi, num_res_cl 
    INTEGER           :: num_atoms_bmim=25, num_atoms_tfsi=15, num_atoms_cl=1
    INTEGER           :: tot_atoms_bmim, tot_atoms_tfsi, tot_atoms_cl
    REAL(KIND=REAL32) :: avg_bmim_dens, avg_cl_dens, avg_tfsi_dens
    REAL(KIND=REAL32) :: avg_tfsi_N_dens, avg_tfsi_O_dens, avg_tfsi_F_dens
    
    REAL(KIND=REAL32) :: mass_bmim(25), mass_cl(1), mass_tfsi(15)
    REAL(KIND=REAL32) :: sum_mass_bmim, sum_mass_cl, sum_mass_tfsi
    REAL(KIND=REAL32) :: one_by_mass_tfsi, one_by_mass_bmim, one_by_mass_cl
    
    REAL(KIND=REAL32) :: charge_bmim(25), charge_cl(1), charge_tfsi(15)
    REAL(KIND=REAL32) :: sum_charge_bmim, sum_charge_tfsi, sum_charge_cl

    ! REF AXIS DEFINITION
    INTEGER           :: num_ring_atoms = 5, ring_atom_indices(5), index_CR, index_NA_butyl
    REAL(KIND=REAL32) :: ring_atoms(3,5) = 0.0, ring_center(3)

    ! CLOCK RELATED VARIABLES
    INTEGER(KIND=8)   :: clock_start, clock_end, clock_rate, clock_max

    ! BOOK-KEEPING VARS
    INTEGER           :: i, j, k, l, ios=0, itemp=0, itemp1=0, itemp2=0
    REAL(KIND=REAL32) :: rvtemp(3)=0.0, rtemp, dum1, dum2

    CALL SYSTEM_CLOCK(count_rate=clock_rate)

    !! ENTER THE MASS AND THE CHARGE INFORMATION !!
    bmim_atomic_num(1)  = 6
    bmim_atomic_num(2)  = 1
    bmim_atomic_num(3)  = 6
    bmim_atomic_num(4)  = 1
    bmim_atomic_num(5)  = 7
    bmim_atomic_num(6)  = 7
    bmim_atomic_num(7)  = 6
    bmim_atomic_num(8)  = 1
    bmim_atomic_num(9)  = 6
    bmim_atomic_num(10) = 1
    bmim_atomic_num(11) = 1
    bmim_atomic_num(12) = 6
    bmim_atomic_num(13) = 1
    bmim_atomic_num(14) = 1
    bmim_atomic_num(15) = 6
    bmim_atomic_num(16) = 1
    bmim_atomic_num(17) = 1
    bmim_atomic_num(18) = 6
    bmim_atomic_num(19) = 1
    bmim_atomic_num(20) = 1
    bmim_atomic_num(21) = 1
    bmim_atomic_num(22) = 6
    bmim_atomic_num(23) = 1
    bmim_atomic_num(24) = 1
    bmim_atomic_num(25) = 1

    mass_bmim(1)  = 12.011000
    mass_bmim(2)  =  1.008000
    mass_bmim(3)  = 12.011000
    mass_bmim(4)  =  1.008000
    mass_bmim(5)  = 14.007000
    mass_bmim(6)  = 14.007000
    mass_bmim(7)  = 12.011000
    mass_bmim(8)  =  1.008000
    mass_bmim(9)  = 12.011000
    mass_bmim(10) =  1.008000
    mass_bmim(11) =  1.008000
    mass_bmim(12) = 12.011000
    mass_bmim(13) =  1.008000
    mass_bmim(14) =  1.008000
    mass_bmim(15) = 12.011000
    mass_bmim(16) =  1.008000
    mass_bmim(17) =  1.008000
    mass_bmim(18) = 12.011000
    mass_bmim(19) =  1.008000
    mass_bmim(20) =  1.008000
    mass_bmim(21) =  1.008000
    mass_bmim(22) = 12.011000
    mass_bmim(23) =  1.008000
    mass_bmim(24) =  1.008000
    mass_bmim(25) =  1.008000
    sum_mass_bmim = SUM(mass_bmim)
    one_by_mass_bmim = 1.0/sum_mass_bmim

    mass_cl(1) = 35.453000
    sum_mass_cl = SUM(mass_cl)
    one_by_mass_cl = 1.0/sum_mass_cl

    mass_tfsi(1)  = 14.007000
    mass_tfsi(2)  = 32.065000
    mass_tfsi(3)  = 15.999000
    mass_tfsi(4)  = 15.999000
    mass_tfsi(5)  = 12.011000
    mass_tfsi(6)  = 18.998000
    mass_tfsi(7)  = 18.998000
    mass_tfsi(8)  = 18.998000
    mass_tfsi(9)  = 32.065000
    mass_tfsi(10) = 15.999000
    mass_tfsi(11) = 15.999000
    mass_tfsi(12) = 12.011000
    mass_tfsi(13) = 18.998000
    mass_tfsi(14) = 18.998000
    mass_tfsi(15) = 18.998000
    sum_mass_tfsi = SUM(mass_tfsi)
    one_by_mass_tfsi = 1.0/sum_mass_tfsi

    charge_bmim(1)  = -0.125000
    charge_bmim(2)  =  0.155000
    charge_bmim(3)  = -0.125000
    charge_bmim(4)  =  0.155000
    charge_bmim(5)  =  0.137500
    charge_bmim(6)  =  0.137500
    charge_bmim(7)  = -0.007500
    charge_bmim(8)  =  0.162500
    charge_bmim(9)  = -0.250000
    charge_bmim(10) =  0.120000
    charge_bmim(11) =  0.120000
    charge_bmim(12) = -0.076000
    charge_bmim(13) =  0.098000
    charge_bmim(14) =  0.098000
    charge_bmim(15) = -0.196000
    charge_bmim(16) =  0.098000
    charge_bmim(17) =  0.098000
    charge_bmim(18) = -0.294000
    charge_bmim(19) =  0.098000
    charge_bmim(20) =  0.098000
    charge_bmim(21) =  0.098000
    charge_bmim(22) = -0.250000
    charge_bmim(23) =  0.120000
    charge_bmim(24) =  0.120000
    charge_bmim(25) =  0.120000
    sum_charge_bmim = SUM(charge_bmim)

    charge_cl(1) = -0.6400
    sum_charge_cl = SUM(charge_cl)

    charge_tfsi(1)  =  -0.740000
    charge_tfsi(2)  =   1.090000
    charge_tfsi(3)  =  -0.545000
    charge_tfsi(4)  =  -0.545000
    charge_tfsi(5)  =   0.445000
    charge_tfsi(6)  =  -0.155000
    charge_tfsi(7)  =  -0.155000
    charge_tfsi(8)  =  -0.155000
    charge_tfsi(9)  =   1.090000
    charge_tfsi(10) =  -0.545000
    charge_tfsi(11) =  -0.545000
    charge_tfsi(12) =   0.445000
    charge_tfsi(13) =  -0.155000
    charge_tfsi(14) =  -0.155000
    charge_tfsi(15) =  -0.155000
    sum_charge_tfsi = SUM(charge_tfsi)

!********** READ FROM INPUT FILE **********

    CALL GET_COMMAND_ARGUMENT(1, inpfile_name)
    OPEN(NEWUNIT=inpunit, FILE=inpfile_name, STATUS='OLD', ACTION ='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING INPUT FILE", TRIM(inpfile_name)
    END IF
    PRINT '(A,A)', TRIM(inpfile_name), " OPENED SUCCESSFULLY"

    READ(inpunit,*)
    READ(inpunit,'(A)') xtcfile_name
    WRITE(*,'(A,A)')          "INPUT) XTCFILE                : ", TRIM(xtcfile_name)
    
    READ(inpunit,*)
    READ(inpunit,'(A)') suffix_name
    WRITE(*,'(A,A)')          "INPUT) SUFFIX NAME            : ", TRIM(suffix_name)
    
    READ(inpunit,*)
    READ(inpunit,*) rbin_size
    WRITE(*,'(A,F6.3)')       "INPUT) SDF BIN SIZE (nm)      : ", rbin_size
    
    READ(inpunit,*)
    READ(inpunit,*) rcut_x, rcut_y, rcut_z
    WRITE(*,'(A,3(F6.3,2X))') "INPUT) SDF RCUT  (nm)         : ", rcut_x, rcut_y, rcut_z

    READ(inpunit,*)
    READ(inpunit,*) ring_atom_indices(:)
    WRITE(*,'(A,5(I0,2X))')   "INPUT) RING INDICES           : ", ring_atom_indices(:)

    READ(inpunit,*)
    READ(inpunit,*) index_CR, index_NA_butyl
    WRITE(*,'(A,2(I0,2X))')   "INPUT) INDICES CR,NA          : ", index_CR, index_NA_butyl

    READ(inpunit,*)
    READ(inpunit,*) index_HA_atoms , index_HB_atoms(:)
    WRITE(*,'(A,3(I0,2X))')   "INPUT) INDICES HA,HB          : ", index_HA_atoms , index_HB_atoms(:)
    
    READ(inpunit,*)
    READ(inpunit,*) index_N_atoms
    WRITE(*,'(A,(I0,2X))')    "INPUT) INDICES N (TFSI)       : ", index_N_atoms
    
    READ(inpunit,*)
    READ(inpunit,*) index_O_atoms(:)
    WRITE(*,'(A,4(I0,2X))')   "INPUT) INDICES O (TFSI)       : ", index_O_atoms(:)
    
    READ(inpunit,*)
    READ(inpunit,*) index_F_atoms(:)
    WRITE(*,'(A,6(I0,2X))')   "INPUT) INDICES F (TFSI)       : ", index_F_atoms(:)

    READ(inpunit,*)
    READ(inpunit,*) num_res_bmim, num_res_cl, num_res_tfsi
    WRITE(*,'(A,3(I0,2X))')   "INPUT) BMIM, CL, TFSI         : ", num_res_bmim, num_res_cl, num_res_tfsi

    READ(inpunit,*)
    READ(inpunit,*) rcut_cl, rcut_tfsi
    WRITE(*,'(A,2(F6.3,2X))') "INPUT) RCUT (CL, TFSI) (nm)   : ", rcut_cl, rcut_tfsi
    rcut_cl_sq = rcut_cl**2 ; rcut_tfsi_sq = rcut_tfsi**2

    READ(inpunit,*)
    READ(inpunit,*) cond_maxtfsi_coordnum,  cond_maxcl_coordnum
    WRITE(*,'(A,2(I0,2X))')   "INPUT) CONDITIONAL MAXCOORDNO : ",  cond_maxtfsi_coordnum,  cond_maxcl_coordnum
    cond_num_tfsi_coordnum = cond_maxtfsi_coordnum + 1 ; cond_num_cl_coordnum = cond_maxcl_coordnum + 1 
    
    READ(inpunit,*)
    READ(inpunit,*) skip_frames
    WRITE(*,'(A,I0)')         "INPUT) SKIP FRAMES            : ", skip_frames
    
    READ(inpunit,*)
    READ(inpunit,*) proc_frames
    WRITE(*,'(A,I0)')         "INPUT) PROCESS FRAMES         : ", proc_frames

    
    IF(num_res_bmim .NE. num_res_cl+num_res_tfsi) THEN
          WRITE(*,*) "ERROR) NUMBER OF BMIM RESIDUES DOESNOT MATCH THE NUMBER OF CL+TFSI RESIDUES..."
          STOP
    END IF

    tot_atoms_bmim = num_res_bmim*num_atoms_bmim
    tot_atoms_cl = num_res_cl*num_atoms_cl
    tot_atoms_tfsi = num_res_tfsi*num_atoms_tfsi
    tot_atoms_ha = num_res_bmim*num_HA_atoms
    tot_atoms_hb = num_res_bmim*num_HB_atoms
    ALLOCATE(box_size(3,proc_frames))

    WRITE(*,*) "************* SYSTEM DETAILS ************* "
    WRITE(*,'(A,I0)')   "num cations  : ", num_res_bmim
    WRITE(*,'(A,I0)')   "num cl       : ", num_res_cl
    WRITE(*,'(A,I0)')   "num tfsi     : ", num_res_tfsi
    WRITE(*,'(A,F8.4)') "Cation Charge: ", sum_charge_bmim
    WRITE(*,'(A,F8.4)') "Cl  Charge   : ", sum_charge_cl
    WRITE(*,'(A,F8.4)') "tfsi  Charge : ", sum_charge_tfsi
    WRITE(*,'(A)') "Cation charges are for the 0.50 system. They are not used in the code. Just for info."
    WRITE(*,*) "******************************************** " 

!********** INITIALIZE MORE VARIABLES **********
    CALL test%open(TRIM(xtcfile_name))

    ios = test%read_next()

    myatoms = test%natoms()
    mystep  = test%step(1)
    mytime  = test%time(1)

    mybox = test%box(1)
    box_size(1,1) = mybox(1,1)
    box_size(2,1) = mybox(2,2)
    box_size(3,1) = mybox(3,3)
    CALL test%close()

    ! WRITE SDF SETTINGS TO OUTPUT   
    WRITE(*,*) "************* XTC FILE DETAILS ************* " 
    WRITE(*,*) "READ IN ", ios, "TEST FRAMES"
    WRITE(*,*) "NUMBER OF ATOMS: ", myatoms
    WRITE(*,*) "SIMU STEP: ", mystep
    WRITE(*,*) "TIME STEP: ", mytime

    WRITE(*,'(A,3(F12.6,2X))') "BOX SIZE: ", box_size(:,1)
    
    volume = box_size(1,1)*box_size(2,1)*box_size(3,1)
    WRITE(*,'(A,F12.4)') "VOLUME  : ", volume

    max_bin_x = INT(rcut_x/(rbin_size))
    max_bin_y = INT(rcut_y/(rbin_size))
    max_bin_z = INT(rcut_z/(rbin_size))
    num_bin_x = 2*max_bin_x
    num_bin_y = 2*max_bin_y
    num_bin_z = 2*max_bin_z
    vox_vol = rbin_size**3
    sphere_cut_off = (rcut_x**2+rcut_y**2+rcut_z**2)
    WRITE(*,*) "******************************************** " 

    WRITE(*,*) "************* SDF GRID DETAILS ************* " 
    WRITE(*,'(A,3(F12.6,2X))') "RCUT(nm)     : ", rcut_x, rcut_y, rcut_z
    WRITE(*,'(A,F12.6)') "RBIN(nm)     : ", rbin_size
    WRITE(*,'(A,F12.6)') "Vox Vol(nm^3): ", vox_vol
    WRITE(*,'(A,3(I0,3X))')    "MAX BINS     : ", max_bin_x, max_bin_y, max_bin_z 
    WRITE(*,'(A,3(I0,3X))')    "NUM BINS     : ", num_bin_x+2, num_bin_y+2, num_bin_z+2
    WRITE(*,'(A,3(I0,3X))')    "NUM BINS OUT : ", num_bin_x, num_bin_y, num_bin_z
    WRITE(*,'(A,F12.6)')       "SPHERE CUTOFF: ", SQRT(sphere_cut_off)
    WRITE(*,'(A,F12.6)')       "CL CUTOFF  nm: ", rcut_cl
    WRITE(*,'(A,F12.6)')       "TFSI CUTOFFnm: ", rcut_tfsi
    IF (num_bin_x < 4) THEN
        WRITE(*,*) "TOO COARSE A GRID FOR SDF.. STOPPING"
        STOP
    END IF
    WRITE(*,*) "******************************************** " 

!********** SKIP FRAMES **********
    CALL trj%open(TRIM(xtcfile_name))

    CALL SYSTEM_CLOCK(count=clock_start)
    IF(skip_frames > 0) THEN
        PRINT *, "***** SKIPPING FRAMES ***** : ", skip_frames
        ios = trj%read_next(skip_frames)
    END IF
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
!********** END SKIP FRAMES **********

!********** PROCESS FRAMES **********
    ALLOCATE(cor_bmim(3,num_res_bmim,proc_frames))
    ALLOCATE(com_cl(3,num_res_cl,proc_frames))
    ALLOCATE(com_tfsi(3,num_res_tfsi,proc_frames))
    
    ALLOCATE(coord_cl(3,tot_atoms_cl,proc_frames))
    ALLOCATE(coord_tfsi(3,tot_atoms_tfsi,proc_frames))

    ALLOCATE(coord_ha(3,tot_atoms_ha,proc_frames))
    ALLOCATE(coord_hb(3,tot_atoms_hb,proc_frames))
    
    ALLOCATE(vx_bmim_ring(3,num_res_bmim,proc_frames))
    ALLOCATE(vy_bmim_ring(3,num_res_bmim,proc_frames))
    ALLOCATE(vz_bmim_ring(3,num_res_bmim,proc_frames))

    cor_bmim(:,:,:)     = 0.0 ; com_cl(:,:,:)       = 0.0  ; com_tfsi(:,:,:) = 0.0 
    coord_cl(:,:,:)     = 0.0 ; coord_tfsi(:,:,:)   = 0.0
    coord_ha(:,:,:)     = 0.0 ; coord_hb(:,:,:)     = 0.0 

    vx_bmim_ring(:,:,:) = 0.0 ; vy_bmim_ring(:,:,:) = 0.0  ; vz_bmim_ring(:,:,:) = 0.0
    cube_coord(:,:) = 0.0 ; volume = 0.0

    PRINT *, "***** READING FRAMES ***** : ", proc_frames
    CALL SYSTEM_CLOCK(count=clock_start)
    DO i = 1, proc_frames
        ios = trj%read_next()
        IF (MOD(i,PRINTEVERY)==0) THEN
          WRITE(*,'(A,I0,A,I0,A,F10.2)') "PROCESSING FRAME : ", i, " ;; MD FRAME : ",&
                                       trj%step(1), " ;; MD TIME : ", trj%time(1)
        END IF

        mybox = trj%box(1)
        box_size(1,i) = mybox(1,1)
        box_size(2,i) = mybox(2,2)
        box_size(3,i) = mybox(3,3)
        volume = volume + box_size(1,i)*box_size(2,i)*box_size(3,i)

        !!! FOR THE BMIM-X1-X2 SYSTEMS, COORIDNATES ARE LISTED IN THAT ORDER !!
        ! STEP 1: STORE THE ANION COORDS (FOR CHARGE DENSITY CALCULATION)
        ! DO j = 1, tot_atoms_bmim
        !     coord_bmim(:,j,i) = trj%x(1,j) ! NO NEED TO STORE ALL BMIM COORDS
        ! END DO
        itemp = tot_atoms_bmim
        DO j = 1, tot_atoms_cl
            coord_cl(:,j,i) = trj%x(1,itemp+j)
        END DO
        itemp = itemp + tot_atoms_cl
        DO j = 1, tot_atoms_tfsi
            coord_tfsi(:,j,i) = trj%x(1,itemp+j)
        END DO

        ! STEP 2: STORE THE COM/COR COORDINATES
        ! LOOP OVER ALL BMIM - CENTER OF RING CALCULATION
        itemp = 0 ; itemp1 = 0 ; itemp2 = 0
        DO j = 1, num_res_bmim
            ! STORE HA/HB coordinates
            DO k = 1, num_HA_atoms
                itemp1 = itemp1 + 1
                coord_ha(:,itemp1,i) = trj%x(1,itemp+index_HA_atoms(k))
            END DO
            DO k = 1, num_HB_atoms
                itemp2 = itemp2 + 1
                coord_hb(:,itemp2,i) = trj%x(1,itemp+index_HB_atoms(k))
            END DO

            ! COMPUTE THE COM OF BMIM
            rvtemp(:) = 0.0
            DO k = 1, num_ring_atoms
                rvtemp(:) = rvtemp(:) + trj%x(1,itemp+ring_atom_indices(k))
            END DO
            cor_bmim(:,j,i) = rvtemp(:)*0.2D0
            
            !********** IDENTIFY AXES **********
            ! Y-axis: along the center-CA-HA
            vy_bmim_ring(:,j,i) = trj%x(1,itemp+index_CR) - cor_bmim(:,j,i)
        
            ! Z-axis is the cross product of two non-collinear vectors on the ring
            vtemp(:) = trj%x(1,itemp+index_NA_butyl) - cor_bmim(:,j,i)
            vz_bmim_ring(:,j,i) = vtemp(:) .cross. vy_bmim_ring(:,j,i)

            ! X-axis is the cross product of Y cross Z
            vx_bmim_ring(:,j,i) = vy_bmim_ring(:,j,i) .cross. vz_bmim_ring(:,j,i)

            CALL v_normalize(vx_bmim_ring(:,j,i))
            CALL v_normalize(vy_bmim_ring(:,j,i))
            CALL v_normalize(vz_bmim_ring(:,j,i))

            itemp = itemp + num_atoms_bmim
        END DO

        ! LOOP OVER ALL CL - COM CALCULATION
        DO j = 1, num_res_cl
            rvtemp(:) = 0.0
            DO k = 1, num_atoms_cl
                rvtemp(:) = rvtemp(:) + trj%x(1,itemp+k)*mass_cl(k)
            END DO
            com_cl(:,j,i) = rvtemp(:)
            itemp = itemp + num_atoms_cl 
        END DO

        ! LOOP OVER ALL TFSI - COM CALCULATION
        DO j = 1, num_res_tfsi
            rvtemp(:) = 0.0
            DO k = 1, num_atoms_tfsi
                rvtemp(:) = rvtemp(:) + trj%x(1,itemp+k)*mass_tfsi(k)
            END DO
            com_tfsi(:,j,i) = rvtemp(:)
            itemp = itemp + num_atoms_tfsi 
        END DO
        
        ! STEP 3: STORE THE FIRST FRAME BMIM STRUCTURE 
        IF (i ==1) THEN
            itemp = 0
            DO j = 1, num_atoms_bmim
                rvtemp(:) = trj%x(1,itemp+j) - cor_bmim(:,1,i)
                cube_coord(1,j) = DOT_PRODUCT(vx_bmim_ring(:,1,i),rvtemp(:))
                cube_coord(2,j) = DOT_PRODUCT(vy_bmim_ring(:,1,i),rvtemp(:))
                cube_coord(3,j) = DOT_PRODUCT(vz_bmim_ring(:,1,i),rvtemp(:))
            END DO
        END IF
    END DO
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)

    ! deallocate raw trajectory -- need this to 
    CALL trj%close()
    DEALLOCATE(trj%frameArray)

    volume = volume/REAL(proc_frames)
    com_cl(:,:,:)   = com_cl(:,:,:)*one_by_mass_cl
    com_tfsi(:,:,:) = com_tfsi(:,:,:)*one_by_mass_tfsi
    avg_bmim_dens   = REAL(num_res_bmim)/volume
    avg_cl_dens     = REAL(num_res_cl)/volume
    avg_tfsi_dens   = REAL(num_res_tfsi)/volume
    avg_tfsi_N_dens = REAL(num_res_tfsi*num_N_atoms)/volume
    avg_tfsi_O_dens = REAL(num_res_tfsi*num_O_atoms)/volume
    avg_tfsi_F_dens = REAL(num_res_tfsi*num_F_atoms)/volume
    avg_ha_dens     = REAL(tot_atoms_ha)/volume
    avg_hb_dens     = REAL(tot_atoms_hb)/volume

    WRITE(*,*) "************* AVERAGE DENSITY DETAILS ************* " 
    WRITE(*,'(A,F20.8)') "AVERAGE BOX VOLUME (nm^3)               : ", volume
    WRITE(*,'(A,F20.8)') "AVERAGE BMIM PARTICLE DENSITY (nm^-3)   : ", avg_bmim_dens
    WRITE(*,'(A,F20.8)') "AVERAGE CL PARTICLE DENSITY (nm^-3)     : ", avg_cl_dens
    WRITE(*,'(A,F20.8)') "AVERAGE TFSI PARTICLE DENSITY (nm^-3)   : ", avg_tfsi_dens
    WRITE(*,'(A,F20.8)') "AVERAGE TFSI-N PARTICLE DENSITY (nm^-3) : ", avg_tfsi_N_dens
    WRITE(*,'(A,F20.8)') "AVERAGE TFSI-O PARTICLE DENSITY (nm^-3) : ", avg_tfsi_O_dens
    WRITE(*,'(A,F20.8)') "AVERAGE TFSI-F PARTICLE DENSITY (nm^-3) : ", avg_tfsi_F_dens
    WRITE(*,'(A,F20.8)') "AVERAGE HA     PARTICLE DENSITY (nm^-3) : ", avg_ha_dens
    WRITE(*,'(A,F20.8)') "AVERAGE HB     PARTICLE DENSITY (nm^-3) : ", avg_hb_dens
    WRITE(*,*) "******************************************** " 

!********** END PROCESS FRAMES **********

!********** POPULATE GRID **********
   ALLOCATE(cl_com_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   cl_com_hist(:,:,:) = 0
   ALLOCATE(tfsi_com_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   tfsi_com_hist(:,:,:) = 0
   ALLOCATE(tfsi_N_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   tfsi_N_hist(:,:,:) = 0
   ALLOCATE(tfsi_O_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   tfsi_O_hist(:,:,:) = 0
   ALLOCATE(tfsi_F_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   tfsi_F_hist(:,:,:) = 0

   ALLOCATE(cl_charge_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   cl_charge_dens(:,:,:) = 0.0
   ALLOCATE(tfsi_charge_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1))
   tfsi_charge_dens(:,:,:) = 0.0

   ALLOCATE(bmim_cor_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   bmim_cor_hist(:,:,:) = 0
   ALLOCATE(bmim_ha_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   bmim_ha_hist(:,:,:) = 0
   ALLOCATE(bmim_hb_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
   bmim_hb_hist(:,:,:) = 0

   ALLOCATE(cond_cl_com_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
            0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
   cond_cl_com_hist(:,:,:,:,:) = 0
   ALLOCATE(cond_tfsi_com_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
            0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
   cond_tfsi_com_hist(:,:,:,:,:) = 0
   ALLOCATE(cond_bmim_cor_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
            0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
   cond_bmim_cor_hist(:,:,:,:,:) = 0
   ALLOCATE(cond_bmim_ha_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
            0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
   cond_bmim_ha_hist(:,:,:,:,:) = 0
   ALLOCATE(cond_bmim_hb_hist(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
            0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
   cond_bmim_hb_hist(:,:,:,:,:) = 0
   ALLOCATE(cond_hist_count(0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
   cond_hist_count(:,:) = 0

   ALLOCATE(coord_cl_temp(3,tot_atoms_cl))
   ALLOCATE(coord_tfsi_temp(3,tot_atoms_tfsi))
   ALLOCATE(cor_bmim_temp(3,num_res_bmim))
   ALLOCATE(com_cl_temp(3,num_res_cl))
   ALLOCATE(com_tfsi_temp(3,num_res_tfsi))
   ALLOCATE(vx_ring_temp(3,num_res_bmim))
   ALLOCATE(vy_ring_temp(3,num_res_bmim))
   ALLOCATE(vz_ring_temp(3,num_res_bmim))
   
    PRINT *, "***** POPULATING SDF ***** : ", proc_frames
    CALL SYSTEM_CLOCK(count=clock_start)
!$OMP PARALLEL  DO SCHEDULE(DYNAMIC) DEFAULT(NONE) &
!$OMP& PRIVATE(i,j,k,l,itemp,itemp2,rtemp)&
!$OMP& PRIVATE(xbin,ybin,zbin,vx_ring_temp,vy_ring_temp,vz_ring_temp,xproj,yproj,zproj) &
!$OMP& PRIVATE(cor_bmim_temp,coord_cl_temp,coord_tfsi_temp,delrv,box_temp) &
!$OMP& PRIVATE(temp_cl_coordnum,temp_tfsi_coordnum,temp_cl_coordid,temp_tfsi_coordid) &
!$OMP& SHARED(rcut_x,rcut_y,rcut_z,rcut_cl_sq,rcut_tfsi_sq,rbin_size,proc_frames,sphere_cut_off) &
!$OMP& SHARED(num_atoms_cl,num_atoms_tfsi,num_res_bmim,num_res_cl,num_res_tfsi) &
!$OMP& SHARED(vx_bmim_ring,vy_bmim_ring,vz_bmim_ring,cor_bmim) &
!$OMP& SHARED(com_cl,com_tfsi,coord_cl,coord_tfsi,coord_ha,coord_hb) &
!$OMP& SHARED(charge_cl,charge_tfsi,box_size) &
!$OMP& SHARED(index_N_atoms,index_O_atoms,index_F_atoms,num_HA_atoms,num_HB_atoms) &
!$OMP& SHARED(cond_maxcl_coordnum,cond_maxtfsi_coordnum) &
!$OMP& REDUCTION(+:cl_com_hist,tfsi_com_hist,cl_charge_dens,tfsi_charge_dens)&
!$OMP& REDUCTION(+:tfsi_N_hist,tfsi_O_hist,tfsi_F_hist)&
!$OMP& REDUCTION(+:bmim_ha_hist,bmim_hb_hist,bmim_cor_hist)&
!$OMP& REDUCTION(+:cond_cl_com_hist,cond_tfsi_com_hist,cond_bmim_cor_hist)&
!$OMP& REDUCTION(+:cond_bmim_ha_hist,cond_bmim_hb_hist,cond_hist_count)
    DO i = 1, proc_frames
    	IF (MOD(i,PRINTEVERY)==0) THEN
          WRITE(*,'(A,I0)') "PROCESSING FRAME : ", i
        END IF

        cor_bmim_temp(:,:)  = cor_bmim(:,:,i)
        vx_ring_temp(:,:) = vx_bmim_ring(:,:,i)
        vy_ring_temp(:,:) = vy_bmim_ring(:,:,i)
        vz_ring_temp(:,:) = vz_bmim_ring(:,:,i)
        box_temp(:)       = box_size(:,i)

        coord_cl_temp(:,:)   = coord_cl(:,:,i)
        coord_tfsi_temp(:,:) = coord_tfsi(:,:,i)

        ! LOOP OVER ALL BMIM
        DO l = 1, num_res_bmim
            ! COORDSHELL DETAILS OF THIS BMIM
            temp_cl_coordnum = 0 ; temp_tfsi_coordnum = 0
            temp_cl_coordid(:) = 0 ; temp_tfsi_coordid(:) = 0
            
            ! LOOP OVER ALL CL
            itemp =0
        	DO j = 1, num_res_cl
                ! CHARGE HISTOGRAM
                itemp = num_atoms_cl*(j-1)
            	DO k = 1, num_atoms_cl
                    delrv(:) = coord_cl_temp(:,itemp+k) - cor_bmim_temp(:,l)
                    delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                    rtemp = DOT_PRODUCT(delrv,delrv)
                    IF(rtemp .LT. sphere_cut_off) THEN
                        xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                        yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                        zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                        IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                            xbin  = FLOOR(xproj/rbin_size) + 1
                            ybin  = FLOOR(yproj/rbin_size) + 1
                            zbin  = FLOOR(zproj/rbin_size) + 1
                            cl_charge_dens(zbin,ybin,xbin) = cl_charge_dens(zbin,ybin,xbin) + charge_cl(k)
                        END IF
                    END IF
            	END DO

                ! COM HISTOGRAM
                delrv(:) = com_cl(:,j,i) - cor_bmim_temp(:,l)
                delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                rtemp = DOT_PRODUCT(delrv,delrv)
                IF(rtemp .LT. sphere_cut_off) THEN
                    xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                    yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                    zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                    IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                        xbin  = FLOOR(xproj/rbin_size) + 1
                        ybin  = FLOOR(yproj/rbin_size) + 1
                        zbin  = FLOOR(zproj/rbin_size) + 1
                        cl_com_hist(zbin,ybin,xbin) = cl_com_hist(zbin,ybin,xbin) + 1
                    END IF
                END IF

                ! CHECK THE SOLVATION SHELL
                IF (rtemp .LT. rcut_cl_sq) THEN 
                    temp_cl_coordnum = temp_cl_coordnum + 1
                    temp_cl_coordid(temp_cl_coordnum) = j
                END IF
            END DO
            
            ! LOOP OVER ALL TFSI
            itemp =0
        	DO j = 1, num_res_tfsi
                ! CHARGE HISTOGRAM
                itemp = num_atoms_tfsi*(j-1)
            	DO k = 1, num_atoms_tfsi
                    delrv(:) = coord_tfsi_temp(:,itemp+k) - cor_bmim_temp(:,l)
                    delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                    rtemp = DOT_PRODUCT(delrv,delrv)
                    IF(rtemp .LT. sphere_cut_off) THEN
                        xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                        yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                        zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                        IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                            xbin  = FLOOR(xproj/rbin_size) + 1
                            ybin  = FLOOR(yproj/rbin_size) + 1
                            zbin  = FLOOR(zproj/rbin_size) + 1
                            tfsi_charge_dens(zbin,ybin,xbin) = tfsi_charge_dens(zbin,ybin,xbin) + charge_tfsi(k)
                            
                            IF(k == index_N_atoms(1)) THEN
                                tfsi_N_hist(zbin,ybin,xbin) = tfsi_N_hist(zbin,ybin,xbin) + 1
                            ELSE IF((k==index_O_atoms(1)).OR.(k==index_O_atoms(2)).OR.&
                                    (k==index_O_atoms(3)).OR.(k==index_O_atoms(4))) THEN
                                tfsi_O_hist(zbin,ybin,xbin) = tfsi_O_hist(zbin,ybin,xbin) + 1
                            ELSE IF((k==index_F_atoms(1)).OR.(k==index_F_atoms(2)).OR.&
                                    (k==index_F_atoms(3)).OR.(k==index_F_atoms(4)).OR.&
                                    (k==index_F_atoms(5)).OR.(k==index_F_atoms(6))) THEN
                                tfsi_F_hist(zbin,ybin,xbin) = tfsi_F_hist(zbin,ybin,xbin) + 1
                            END IF
                        END IF
                    END IF
            	END DO

                ! COM HISTOGRAM
                delrv(:) = com_tfsi(:,j,i) - cor_bmim_temp(:,l)
                delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                rtemp = DOT_PRODUCT(delrv,delrv)
                IF(rtemp .LT. sphere_cut_off) THEN
                    xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                    yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                    zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                    IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                        xbin  = FLOOR(xproj/rbin_size) + 1
                        ybin  = FLOOR(yproj/rbin_size) + 1
                        zbin  = FLOOR(zproj/rbin_size) + 1
                        tfsi_com_hist(zbin,ybin,xbin) = tfsi_com_hist(zbin,ybin,xbin) + 1
                    END IF
                END IF

                ! CHECK THE SOLVATION SHELL
                IF (rtemp .LT. rcut_tfsi_sq) THEN 
                    temp_tfsi_coordnum = temp_tfsi_coordnum + 1
                    temp_tfsi_coordid(temp_tfsi_coordnum) = j
                END IF
            END DO

            ! CONDITIONAL SDFs and BMIM-CAT SDFs
            IF(temp_cl_coordnum .LE. cond_maxcl_coordnum .AND. temp_tfsi_coordnum .LE. cond_maxtfsi_coordnum) THEN
                ! BMIM-CL SDF
                DO j = 1, temp_cl_coordnum
                    itemp2 = temp_cl_coordid(j)
                    delrv(:) = com_cl(:,itemp2,i) - cor_bmim_temp(:,l)
                    delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                    xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                    yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                    zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                    xbin  = FLOOR(xproj/rbin_size) + 1
                    ybin  = FLOOR(yproj/rbin_size) + 1
                    zbin  = FLOOR(zproj/rbin_size) + 1
                    cond_cl_com_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) = &
                                    cond_cl_com_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) + 1
                END DO

                ! BMIM-TFSI SDF
                DO j = 1, temp_tfsi_coordnum
                    itemp2 = temp_tfsi_coordid(j)
                    delrv(:) = com_tfsi(:,itemp2,i) - cor_bmim_temp(:,l)
                    delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                    xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                    yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                    zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                    xbin  = FLOOR(xproj/rbin_size) + 1
                    ybin  = FLOOR(yproj/rbin_size) + 1
                    zbin  = FLOOR(zproj/rbin_size) + 1
                    cond_tfsi_com_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) = &
                                    cond_tfsi_com_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) + 1
                END DO
                

                ! LOOP OVER ALL OTHER BMI (CANNOT DO N(N-1)/2)
                DO j = 1, num_res_bmim
                    IF (j == l) CYCLE     ! AVOID SELF COUNTING

                    !COR HISTOGRAM
                    delrv(:) = cor_bmim_temp(:,j) - cor_bmim_temp(:,l)
                    delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                    rtemp = DOT_PRODUCT(delrv,delrv)
                    IF(rtemp .LT. sphere_cut_off) THEN
                        xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                        yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                        zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                        IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                            xbin  = FLOOR(xproj/rbin_size) + 1
                            ybin  = FLOOR(yproj/rbin_size) + 1
                            zbin  = FLOOR(zproj/rbin_size) + 1
                            cond_bmim_cor_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) = &
                                    cond_bmim_cor_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) + 1
                            bmim_cor_hist(zbin,ybin,xbin) = bmim_cor_hist(zbin,ybin,xbin) + 1
                        END IF
                    END IF

                    ! HA HISTOGRAM
                    DO k = 1, num_HA_atoms
                        delrv(:) = coord_ha(:,(j-1)*num_HA_atoms+k,i) - cor_bmim_temp(:,l)
                        delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                        rtemp = DOT_PRODUCT(delrv,delrv)
                        IF(rtemp .LT. sphere_cut_off) THEN
                            xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                            yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                            zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                            IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                                xbin  = FLOOR(xproj/rbin_size) + 1
                                ybin  = FLOOR(yproj/rbin_size) + 1
                                zbin  = FLOOR(zproj/rbin_size) + 1
                                cond_bmim_ha_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) = &
                                    cond_bmim_ha_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) + 1
                                bmim_ha_hist(zbin,ybin,xbin) = bmim_ha_hist(zbin,ybin,xbin) + 1
                            END IF
                        END IF
                    END DO

                    ! HB HISTOGRAM
                    DO k = 1, num_HB_atoms
                        delrv(:) = coord_hb(:,(j-1)*num_HB_atoms+k,i) - cor_bmim_temp(:,l)
                        delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                        rtemp = DOT_PRODUCT(delrv,delrv)
                        IF(rtemp .LT. sphere_cut_off) THEN
                            xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                            yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                            zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                            IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                                xbin  = FLOOR(xproj/rbin_size) + 1
                                ybin  = FLOOR(yproj/rbin_size) + 1
                                zbin  = FLOOR(zproj/rbin_size) + 1
                                cond_bmim_hb_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) = &
                                    cond_bmim_hb_hist(zbin,ybin,xbin,temp_tfsi_coordnum,temp_cl_coordnum) + 1
                                bmim_hb_hist(zbin,ybin,xbin) = bmim_hb_hist(zbin,ybin,xbin) + 1
                            END IF
                        END IF
                    END DO
                END DO

                ! COUNT THE NUM OF CATIONS SATISFYING THE CONDITION TO NORAMLIZE LATER
                cond_hist_count(temp_tfsi_coordnum,temp_cl_coordnum) = cond_hist_count(temp_tfsi_coordnum,temp_cl_coordnum) + 1 
            ELSE
                ! LOOP OVER ALL OTHER BMI (CANNOT DO N(N-1)/2)
                DO j = 1, num_res_bmim
                    IF (j == l) CYCLE     ! AVOID SELF COUNTING

                    !COR HISTOGRAM
                    delrv(:) = cor_bmim_temp(:,j) - cor_bmim_temp(:,l)
                    delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                    rtemp = DOT_PRODUCT(delrv,delrv)
                    IF(rtemp .LT. sphere_cut_off) THEN
                        xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                        yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                        zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                        IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                            xbin  = FLOOR(xproj/rbin_size) + 1
                            ybin  = FLOOR(yproj/rbin_size) + 1
                            zbin  = FLOOR(zproj/rbin_size) + 1
                            bmim_cor_hist(zbin,ybin,xbin) = bmim_cor_hist(zbin,ybin,xbin) + 1
                        END IF
                    END IF

                    ! HA HISTOGRAM
                    DO k = 1, num_HA_atoms
                        delrv(:) = coord_ha(:,(j-1)*num_HA_atoms+k,i) - cor_bmim_temp(:,l)
                        delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                        rtemp = DOT_PRODUCT(delrv,delrv)
                        IF(rtemp .LT. sphere_cut_off) THEN
                            xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                            yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                            zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                            IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                                xbin  = FLOOR(xproj/rbin_size) + 1
                                ybin  = FLOOR(yproj/rbin_size) + 1
                                zbin  = FLOOR(zproj/rbin_size) + 1
                                bmim_ha_hist(zbin,ybin,xbin) = bmim_ha_hist(zbin,ybin,xbin) + 1
                            END IF
                        END IF
                    END DO

                    ! HB HISTOGRAM
                    DO k = 1, num_HB_atoms
                        delrv(:) = coord_hb(:,(j-1)*num_HB_atoms+k,i) - cor_bmim_temp(:,l)
                        delrv(:) = delrv(:) - box_temp(:)*NINT(delrv(:)/box_temp(:))
                        rtemp = DOT_PRODUCT(delrv,delrv)
                        IF(rtemp .LT. sphere_cut_off) THEN
                            xproj = DOT_PRODUCT(vx_ring_temp(:,l),delrv(:))
                            yproj = DOT_PRODUCT(vy_ring_temp(:,l),delrv(:))
                            zproj = DOT_PRODUCT(vz_ring_temp(:,l),delrv(:))
                            IF((ABS(xproj) < rcut_x) .AND. (ABS(yproj) < rcut_y) .AND. (ABS(zproj) < rcut_z)) THEN
                                xbin  = FLOOR(xproj/rbin_size) + 1
                                ybin  = FLOOR(yproj/rbin_size) + 1
                                zbin  = FLOOR(zproj/rbin_size) + 1
                                bmim_hb_hist(zbin,ybin,xbin) = bmim_hb_hist(zbin,ybin,xbin) + 1
                            END IF
                        END IF
                    END DO
                END DO
            END IF 

        END DO
    END DO
!$OMP END PARALLEL DO
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)

    DEALLOCATE(cor_bmim,com_cl,com_tfsi)
    DEALLOCATE(coord_cl,coord_tfsi)
    DEALLOCATE(coord_cl_temp,coord_tfsi_temp)

    ALLOCATE(print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x))
    print_dens = 0.0

    ! ALL THE CALCULATIONS ARE DONE IN NM AND THEN FINALLY SCALED TO ANGSTROM HERE -- through vox_vol
    ! UPDATE ON 02 AUG 2021 -- NM^-3 SEEM TO BE A NATURAL UNITS FOR THIS PROBLEM. CHANGING TO THEM
    ALLOCATE(cl_com_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_com_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_N_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_O_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_F_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_N_norm_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_O_norm_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(tfsi_F_norm_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(bmim_cor_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(bmim_ha_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    ALLOCATE(bmim_hb_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1)) 
    
    norm_const = 1.D0/(REAL(proc_frames*num_res_bmim)*vox_vol)
    cl_com_dens = REAL(cl_com_hist)*norm_const
    tfsi_com_dens = REAL(tfsi_com_hist)*norm_const
    tfsi_N_dens = REAL(tfsi_N_hist)*norm_const
    tfsi_O_dens = REAL(tfsi_O_hist)*norm_const
    tfsi_F_dens = REAL(tfsi_F_hist)*norm_const

    cl_charge_dens = cl_charge_dens*norm_const
    tfsi_charge_dens = tfsi_charge_dens*norm_const

    bmim_cor_dens = REAL(bmim_cor_hist)*norm_const
    bmim_ha_dens = REAL(bmim_ha_hist)*norm_const
    bmim_hb_dens = REAL(bmim_hb_hist)*norm_const

    ALLOCATE(cond_cl_com_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
             0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
    ALLOCATE(cond_tfsi_com_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
             0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
    ALLOCATE(cond_bmim_cor_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
             0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum)) 
    ALLOCATE(cond_bmim_ha_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
             0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum))
    ALLOCATE(cond_bmim_hb_dens(-1*max_bin_z:max_bin_z+1, -1*max_bin_y:max_bin_y+1, -1*max_bin_x:max_bin_x+1, &
             0:cond_maxtfsi_coordnum, 0:cond_maxcl_coordnum))

    DO i = 0, cond_maxtfsi_coordnum
        DO j = 0, cond_maxcl_coordnum
            WRITE(*,*) 
            WRITE(*,'(A,F16.4)') "PERCENT CONDITIONAL BMI      (NORMALIZED by frames and cations): ", & 
                                100.D0*REAL(cond_hist_count(i,j))/REAL(proc_frames*num_res_bmim)
            WRITE(*,'(A)') "-------------------------------------------------------------------------------"
            IF (cond_hist_count(i,j) .GT. 0) THEN
                cond_norm_const = 1.D0/(REAL(cond_hist_count(i,j))*vox_vol)
                cond_cl_com_dens(:,:,:,i,j) = REAL(cond_cl_com_hist(:,:,:,i,j))*cond_norm_const
                cond_tfsi_com_dens(:,:,:,i,j) = REAL(cond_tfsi_com_hist(:,:,:,i,j))*cond_norm_const
                cond_bmim_cor_dens(:,:,:,i,j) = REAL(cond_bmim_cor_hist(:,:,:,i,j))*cond_norm_const
                cond_bmim_ha_dens(:,:,:,i,j) = REAL(cond_bmim_ha_hist(:,:,:,i,j))*cond_norm_const
                cond_bmim_hb_dens(:,:,:,i,j) = REAL(cond_bmim_hb_hist(:,:,:,i,j))*cond_norm_const

                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL CL HISTOGRAM      ( ",i," , ",j," ): ", & 
                                SUM(cond_cl_com_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL TFSI HISTOGRAM    ( ",i," , ",j," ): ", &
                                SUM(cond_tfsi_com_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL BMI-COR HISTOGRAM ( ",i," , ",j," ): ", & 
                                SUM(cond_bmim_cor_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL BMI-HA HISTOGRAM  ( ",i," , ",j," ): ", & 
                                SUM(cond_bmim_ha_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL BMI-HB HISTOGRAM  ( ",i," , ",j," ): ", & 
                                SUM(cond_bmim_hb_dens(:,:,:,i,j))*vox_vol
            ELSE
                cond_cl_com_dens(:,:,:,i,j) = 0.D0
                cond_tfsi_com_dens(:,:,:,i,j) = 0.D0
                cond_bmim_cor_dens(:,:,:,i,j) = 0.D0
                cond_bmim_ha_dens(:,:,:,i,j) = 0.D0
                cond_bmim_hb_dens(:,:,:,i,j) = 0.D0

                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL CL HISTOGRAM      ( ",i," , ",j," ): ", & 
                                SUM(cond_cl_com_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL TFSI HISTOGRAM    ( ",i," , ",j," ): ", &
                                SUM(cond_tfsi_com_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL BMI-COR HISTOGRAM ( ",i," , ",j," ): ", & 
                                SUM(cond_bmim_cor_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL BMI-HA HISTOGRAM  ( ",i," , ",j," ): ", & 
                                SUM(cond_bmim_ha_dens(:,:,:,i,j))*vox_vol
                WRITE(*,'(A,I2,A,I2,A,F16.4)') "TOTAL CONDITIONAL BMI-HB HISTOGRAM  ( ",i," , ",j," ): ", & 
                                SUM(cond_bmim_hb_dens(:,:,:,i,j))*vox_vol
            END IF
            WRITE(*,'(A)') "-------------------------------------------------------------------------------"
            WRITE(*,*)
        END DO
    END DO

    WRITE(*,'(A,F16.4)') "TOTAL CL HISTOGRAM (NORMALIZED by frames and cations)         : ", & 
                                SUM(REAL(cl_com_hist))/REAL(proc_frames*num_res_bmim)
    WRITE(*,'(A,F16.4)') "TOTAL TFSI HISTOGRAM (NORMALIZED by frames and cations)       : ", & 
                                SUM(REAL(tfsi_com_hist))/REAL(proc_frames*num_res_bmim)                                
    WRITE(*,'(A,F16.4)') "TOTAL TFSI-N HISTOGRAM (NORMALIZED by frames and cations)     : ", & 
                                SUM(REAL(tfsi_N_hist))/REAL(proc_frames*num_res_bmim)                                
    WRITE(*,'(A,F16.4)') "TOTAL TFSI-O HISTOGRAM (NORMALIZED by frames and cations)     : ", & 
                                SUM(REAL(tfsi_O_hist))/REAL(proc_frames*num_res_bmim)                                
    WRITE(*,'(A,F16.4)') "TOTAL TFSI-F HISTOGRAM (NORMALIZED by frames and cations)     : ", & 
                                SUM(REAL(tfsi_F_hist))/REAL(proc_frames*num_res_bmim)                                                                                                                                

    WRITE(*,'(A,F16.8)') "TOTAL CL CHARGE HISTOGRAM (NORMALIZED by frames and cations)  : ", & 
                                SUM(cl_charge_dens)*vox_vol
    WRITE(*,'(A,F16.8)') "TOTAL TFSI CHARGE HISTOGRAM (NORMALIZED by frames and cations): ", & 
                                SUM(tfsi_charge_dens)*vox_vol


    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 1 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-ring-cl-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 1"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-CL COM"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        cl_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-ring-cl-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 1n"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-CL COM NORMALIZED"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1n"
    print_dens = print_dens/avg_cl_dens
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 2 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-ring-tfsi-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 2"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI COM"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 2"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        tfsi_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-ring-tfsi-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 2n"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI COM NORMALIZED"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 2n"
    print_dens = print_dens/avg_tfsi_dens
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 3 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!

    OPEN(NEWUNIT=sdf_unit, FILE="charge-sdf-ring-cl-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 3"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-CL-CHARGE"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 3"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        cl_charge_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 4 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!

    OPEN(NEWUNIT=sdf_unit, FILE="charge-sdf-ring-tfsi-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 4"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-CHARGE"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 4"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        tfsi_charge_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 5 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-ring-tfsi-N-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 5"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-N-ATOMS"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 5"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        tfsi_N_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-ring-tfsi-N-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 5n"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-N-ATOMS NORMALIZED"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 5n"
    print_dens = print_dens/avg_tfsi_N_dens
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 6 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-ring-tfsi-O-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 6"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-O-ATOMS"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 6"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        tfsi_O_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-ring-tfsi-O-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 6n"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-O-ATOMS NORMALIZED"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 6n"
    print_dens = print_dens/avg_tfsi_O_dens
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE 7 ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-ring-tfsi-F-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 7"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-F-ATOMS"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 7"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        tfsi_F_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-ring-tfsi-F-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE 7n"
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-TFSI-F-ATOMS NORMALIZED"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 7n"
    print_dens = print_dens/avg_tfsi_F_dens
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE ********************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-bmim-cor-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE "
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-RING COG"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        bmim_cor_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !*********************************************************
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-bmim-ha-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE "
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-HA"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        bmim_ha_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !*********************************************************
    OPEN(NEWUNIT=sdf_unit, FILE="sdf-bmim-hb-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE "
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-HB"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        bmim_hb_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!****************************** PRINTING TO CUBE FILE ********************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-bmim-cor-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE "
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-RING COG"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    norm_const = volume/(REAL(num_res_bmim-1))
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        bmim_cor_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*norm_const
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !*********************************************************
    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-bmim-ha-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE "
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-HA"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    norm_const = volume/(REAL(tot_atoms_ha-1))
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        bmim_ha_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*norm_const
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !*********************************************************
    OPEN(NEWUNIT=sdf_unit, FILE="norm-sdf-bmim-hb-"//TRIM(suffix_name)//".cube")
    PRINT *, "PRINTING HEADER OF CUBE FILE "
    WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: RING-HB"
    WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    dum1 = 0.D0
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                -max_bin_y*rbin_size*10.0*au_conv, &
                                                -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
    WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
    DO i = 1, num_atoms_bmim
        WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
    END DO

    PRINT *, "PRINTING BODY OF CUBE FILE 1"
    norm_const = volume/(REAL(tot_atoms_hb-2))
    print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
        bmim_hb_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*norm_const
    DO i = -1*max_bin_x+1, max_bin_x
        DO j = -1*max_bin_y+1, max_bin_y
            WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
        END DO
    END DO
    CLOSE(sdf_unit)

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!*************************** PRINTING TO COND CUBE FILE ******************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!
    DO k = 0, cond_maxtfsi_coordnum
        DO l = 0, cond_maxcl_coordnum

            WRITE(outfilename,'(A,I0,A,I0,A,A,A)') "cond-sdf-cl-com-", k, "-", l, "-", TRIM(suffix_name), ".cube"
            OPEN(NEWUNIT=sdf_unit, FILE=TRIM(outfilename))
            PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
            WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL SDF"
            WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
            dum1 = 0.D0
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                     -max_bin_y*rbin_size*10.0*au_conv, &
                                                     -max_bin_z*rbin_size*10.0*au_conv              ! nm to Ang conversion here 
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1    !  nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1    ! nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv    !  nm to Ang conversion here
            
            DO i = 1, num_atoms_bmim
                WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
            END DO

            PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE"
            print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
                 cond_cl_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x,k,l)
            DO i = -1*max_bin_x+1, max_bin_x
                DO j = -1*max_bin_y+1, max_bin_y
                    WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
                END DO
            END DO
            CLOSE(sdf_unit)

            !*********************************************************
            WRITE(outfilename,'(A,I0,A,I0,A,A,A)') "cond-sdf-tfsi-com-", k, "-", l, "-", TRIM(suffix_name), ".cube"
            OPEN(NEWUNIT=sdf_unit, FILE=TRIM(outfilename))
            PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
            WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL SDF"
            WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
            dum1 = 0.D0
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                     -max_bin_y*rbin_size*10.0*au_conv, &
                                                     -max_bin_z*rbin_size*10.0*au_conv              ! nm to Ang conversion here 
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1    !  nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1    ! nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv    !  nm to Ang conversion here
            
            DO i = 1, num_atoms_bmim
                WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
            END DO

            PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE"
            print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
                 cond_tfsi_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x,k,l)
            DO i = -1*max_bin_x+1, max_bin_x
                DO j = -1*max_bin_y+1, max_bin_y
                    WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
                END DO
            END DO
            CLOSE(sdf_unit)

            !*********************************************************
            WRITE(outfilename,'(A,I0,A,I0,A,A,A)') "cond-sdf-bmim-cor-", k, "-", l, "-", TRIM(suffix_name), ".cube"
            OPEN(NEWUNIT=sdf_unit, FILE=TRIM(outfilename))
            PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
            WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL SDF"
            WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
            dum1 = 0.D0
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                     -max_bin_y*rbin_size*10.0*au_conv, &
                                                     -max_bin_z*rbin_size*10.0*au_conv              ! nm to Ang conversion here 
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1    !  nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1    ! nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv    !  nm to Ang conversion here
            
            DO i = 1, num_atoms_bmim
                WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
            END DO

            PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE"
            print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
                 cond_bmim_cor_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x,k,l)
            DO i = -1*max_bin_x+1, max_bin_x
                DO j = -1*max_bin_y+1, max_bin_y
                    WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
                END DO
            END DO
            CLOSE(sdf_unit)

            !*********************************************************
            WRITE(outfilename,'(A,I0,A,I0,A,A,A)') "cond-sdf-bmim-ha-", k, "-", l, "-", TRIM(suffix_name), ".cube"
            OPEN(NEWUNIT=sdf_unit, FILE=TRIM(outfilename))
            PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
            WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL SDF"
            WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
            dum1 = 0.D0
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                     -max_bin_y*rbin_size*10.0*au_conv, &
                                                     -max_bin_z*rbin_size*10.0*au_conv              ! nm to Ang conversion here 
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1    !  nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1    ! nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv    !  nm to Ang conversion here
            
            DO i = 1, num_atoms_bmim
                WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
            END DO

            PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE"
            print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
                 cond_bmim_ha_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x,k,l)
            DO i = -1*max_bin_x+1, max_bin_x
                DO j = -1*max_bin_y+1, max_bin_y
                    WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
                END DO
            END DO
            CLOSE(sdf_unit)

            !*********************************************************
            WRITE(outfilename,'(A,I0,A,I0,A,A,A)') "cond-sdf-bmim-hb-", k, "-", l, "-", TRIM(suffix_name), ".cube"
            OPEN(NEWUNIT=sdf_unit, FILE=TRIM(outfilename))
            PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
            WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL SDF"
            WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
            dum1 = 0.D0
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
                                                     -max_bin_y*rbin_size*10.0*au_conv, &
                                                     -max_bin_z*rbin_size*10.0*au_conv              ! nm to Ang conversion here 
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1    !  nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1    ! nm to Ang conversion here
            WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv    !  nm to Ang conversion here
            
            DO i = 1, num_atoms_bmim
                WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
            END DO

            PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE"
            print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
                 cond_bmim_hb_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x,k,l)
            DO i = -1*max_bin_x+1, max_bin_x
                DO j = -1*max_bin_y+1, max_bin_y
                    WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
                END DO
            END DO
            CLOSE(sdf_unit)

        END DO
    END DO

    !**********************************************************************************************************************************!
    !!!!!!!!!!!!!!!!!!!!!!!!********************************** END PRINTING  ************************************!!!!!!!!!!!!!!!!!!!!!!!
    !**********************************************************************************************************************************!

    WRITE(*,*) "$$$$$$$$$$$$$$$$$$$ END OF CODE $$$$$$$$$$$$$$$$$$$$$"
    WRITE(*,*) "                   ALL THE BEST :)                   "

    CONTAINS
    SUBROUTINE get_DDMonYY(date)
        CHARACTER(len=11), INTENT(out) :: date

        CHARACTER(len=2) :: dd
        CHARACTER(len=3) :: mons(12)
        CHARACTER(len=4) :: yyyy
        INTEGER :: values(8)

        mons = ['Jan','Feb','Mar','Apr','May','Jun',&
          'Jul','Aug','Sep','Oct','Nov','Dec']

        CALL DATE_AND_TIME(VALUES=values)

        WRITE(  dd,'(i2)') values(3)
        WRITE(yyyy,'(i4)') values(1)

        date = dd//' '//mons(values(2))//' '//yyyy
    END SUBROUTINE get_DDMonYY

END PROGRAM sdf

!     
!     OPEN(NEWUNIT=sdf_unit, FILE="cond-sdf-ring-cl-"//TRIM(suffix_name)//".cube")
!     PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
!     WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL RING-CL COM"
!     WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
!     dum1 = 0.D0
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
!                                                 -max_bin_y*rbin_size*10.0*au_conv, &
!                                                 -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
!     DO i = 1, num_atoms_bmim
!         WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
!     END DO

!     PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE 1"
!     print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
!         cond_cl_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
!     DO i = -1*max_bin_x+1, max_bin_x
!         DO j = -1*max_bin_y+1, max_bin_y
!             WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
!         END DO
!     END DO
!     CLOSE(sdf_unit)

!     !*********************************************************

!     OPEN(NEWUNIT=sdf_unit, FILE="cond-norm1-sdf-ring-cl-"//TRIM(suffix_name)//".cube")
!     PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
!     WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL NORM1 RING-CL COM"
!     WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
!     dum1 = 0.D0
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
!                                                 -max_bin_y*rbin_size*10.0*au_conv, &
!                                                 -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
!     DO i = 1, num_atoms_bmim
!         WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
!     END DO

!     PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE 1"
!     cond_norm_const = 1.D0/(REAL(cond_maxcl_coordnum))
!     print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
!         cond_cl_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*cond_norm_const
!     DO i = -1*max_bin_x+1, max_bin_x
!         DO j = -1*max_bin_y+1, max_bin_y
!             WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
!         END DO
!     END DO
!     CLOSE(sdf_unit)

!     !*********************************************************

!     OPEN(NEWUNIT=sdf_unit, FILE="cond-norm2-sdf-ring-cl-"//TRIM(suffix_name)//".cube")
!     PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
!     WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL NORM2 RING-CL COM"
!     WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
!     dum1 = 0.D0
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
!                                                 -max_bin_y*rbin_size*10.0*au_conv, &
!                                                 -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
!     DO i = 1, num_atoms_bmim
!         WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
!     END DO

!     PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE 1"
!     cond_norm_const = ((4.0/3.0)*pi*rcut_cl**3)/REAL(cond_maxcl_coordnum)
!     print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
!         cond_cl_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*cond_norm_const
!     DO i = -1*max_bin_x+1, max_bin_x
!         DO j = -1*max_bin_y+1, max_bin_y
!             WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
!         END DO
!     END DO
!     CLOSE(sdf_unit)

!     !*********************************************************!*********************************************************

!     OPEN(NEWUNIT=sdf_unit, FILE="cond-sdf-ring-tfsi-"//TRIM(suffix_name)//".cube")
!     PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
!     WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL NORM1 RING-TFSI COM"
!     WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
!     dum1 = 0.D0
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
!                                                 -max_bin_y*rbin_size*10.0*au_conv, &
!                                                 -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
!     DO i = 1, num_atoms_bmim
!         WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
!     END DO

!     PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE 1"
!     print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
!         cond_tfsi_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)
!     DO i = -1*max_bin_x+1, max_bin_x
!         DO j = -1*max_bin_y+1, max_bin_y
!             WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
!         END DO
!     END DO
!     CLOSE(sdf_unit)

!     !*********************************************************

!     OPEN(NEWUNIT=sdf_unit, FILE="cond-norm1-sdf-ring-tfsi-"//TRIM(suffix_name)//".cube")
!     PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
!     WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL RING-TFSI COM"
!     WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
!     dum1 = 0.D0
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
!                                                 -max_bin_y*rbin_size*10.0*au_conv, &
!                                                 -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
!     DO i = 1, num_atoms_bmim
!         WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
!     END DO

!     PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE 1"
!     cond_norm_const = 1.D0/(REAL(cond_maxtfsi_coordnum))
!     print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
!         cond_tfsi_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*cond_norm_const
!     DO i = -1*max_bin_x+1, max_bin_x
!         DO j = -1*max_bin_y+1, max_bin_y
!             WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
!         END DO
!     END DO
!     CLOSE(sdf_unit)

!     !*********************************************************

!     OPEN(NEWUNIT=sdf_unit, FILE="cond-norm2-sdf-ring-tfsi-"//TRIM(suffix_name)//".cube")
!     PRINT *, "PRINTING HEADER OF CONDITIONAL CUBE FILE "
!     WRITE(sdf_unit,*) "MD GENERATED CUBE FILE: CONDITIONAL NORM2 RING-TFSI COM"
!     WRITE(sdf_unit,*) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
!     dum1 = 0.D0
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') num_atoms_bmim, -max_bin_x*rbin_size*10.0*au_conv, &
!                                                 -max_bin_y*rbin_size*10.0*au_conv, &
!                                                 -max_bin_z*rbin_size*10.0*au_conv         ! nm to Ang conversion here 
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_x, rbin_size*10.0*au_conv, dum1, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_y, dum1, rbin_size*10.0*au_conv, dum1 ! nm to Ang conversion here
!     WRITE(sdf_unit,'(I0, 2X, 3(F10.6,2X))') 2*max_bin_z, dum1, dum1, rbin_size*10.0*au_conv ! nm to Ang conversion here
!     DO i = 1, num_atoms_bmim
!         WRITE(sdf_unit,'(I3, 4(F10.6, 2X))') bmim_atomic_num(i), dum1, cube_coord(:,i)*10.0*au_conv ! nm to Ang conversion here
!     END DO

!     PRINT *, "PRINTING BODY OF CONDITIONAL CUBE FILE 1"
!     cond_norm_const = ((4.0/3.0)*pi*rcut_tfsi**3)/REAL(cond_maxtfsi_coordnum)
!     print_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x) = &
!         cond_tfsi_com_dens(-1*max_bin_z+1:max_bin_z, -1*max_bin_y+1:max_bin_y, -1*max_bin_x+1:max_bin_x)*cond_norm_const
!     DO i = -1*max_bin_x+1, max_bin_x
!         DO j = -1*max_bin_y+1, max_bin_y
!             WRITE(sdf_unit,'(6(E13.6, 2X))') print_dens(:, j, i)
!         END DO
!     END DO
!     CLOSE(sdf_unit)