! PROGRAM TO ESTIMATE THE CONDUCTIVITY OF DEME-TFSI SYSTEM

! FIRST STARTED ON 06 NOV 2020
! LAST MODIFIED ON 06 NOV 2020

! The prefactor PF(V,T) ;; where: cond(t) = PF(V,T)*MSD_coll in {ps S/m}

!                    e^2                              1
! The PF(V,T)  = ---------- = 3.098 747 964 * 10^5  ----- {S/m} {ps/nm^2}
!                 6 Kb V T                           V T
! V in nm^3
! T in K

! As of 04 Nov 2020 the net ion charge of DEME-TFSI system is 0.8024

PROGRAM iconduct
    USE ISO_FORTRAN_ENV
    USE gmxfort_trajectory 
    IMPLICIT NONE

    REAL(KIND=REAL32), PARAMETER :: pi = 4.D0*ATAN(1.D0)
    REAL(KIND=REAL32), PARAMETER :: prefactor = 309874.7964
    CHARACTER(LEN=100), PARAMETER :: last_edited = '06 NOV 2020'
    CHARACTER(LEN=100), PARAMETER :: source_file = 'ic-mj-dsp-v1.F08'
    CHARACTER(LEN=200), PARAMETER :: source_dir = '/home/nikhil/projects/11_codes/&
                                                   conductivity/einstein/deme-tfsi/'
    INTEGER, PARAMETER            :: PRINTEVERY = 10

! CLOCK RELATED VARIABLES
    INTEGER(KIND=8)   :: clock_start, clock_end, clock_rate, clock_max

! DUMP FILE RELATED VARIABLE
    TYPE(Trajectory)                 :: trj, test
    REAL(KIND=REAL32)                :: mybox(3,3), box_size(3), mytime
    INTEGER                          :: mystep, myatoms
    REAL(KIND=REAL32), ALLOCATABLE   :: xyz(:,:)
    REAL(KIND=REAL32)                :: time_per_frame
    REAL(KIND=REAL32)                :: volume = 0.0

! DEME-TFSI RELATED VARIABLES
    INTEGER :: num_res_deme, num_res_tfsi
    INTEGER :: num_atoms_deme=30, num_atoms_tfsi=15
    INTEGER :: tot_atoms_deme, tot_atoms_tfsi
    REAL(KIND=REAL32)    :: mass_deme(30), mass_tfsi(15)
    REAL(KIND=REAL32)    :: sum_mass_deme, sum_mass_tfsi, one_by_mass_deme, one_by_mass_tfsi

! MSD RELATED VARIABLES
    REAL(KIND=REAL32), ALLOCATABLE, DIMENSION(:,:)   :: mj_cat, mj_ani
    REAL(KIND=REAL32), ALLOCATABLE, DIMENSION(:)     :: msdmj_cc, msdmj_aa
    REAL(KIND=REAL32), ALLOCATABLE, DIMENSION(:)     :: msdmj_ca
    REAL(KIND=REAL32), ALLOCATABLE, DIMENSION(:)     :: msdmj

! MSD ALGO RELATED VARIABLES
    REAL(KIND=REAL32)                                :: delmj_cat(3),mj_cat_t0(3)
    REAL(KIND=REAL32)                                :: delmj_ani(3),mj_ani_t0(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)               :: msd_count
    REAL(KIND=REAL32), ALLOCATABLE, DIMENSION(:)     :: one_by_msd_count
    INTEGER                                          :: num_restarts


! CONDUCTIVITY RELATED VARIABLES
    REAL(KIND=REAL32)                                :: ion_charge, temperature
    REAL(KIND=REAL32)                                :: mult_factor

! INTERNAL VARIABLES
    REAL(KIND=REAL32)                  :: delr(3)
    INTEGER                            :: mycount=0
! TIME RELATED VARIABLES
    INTEGER                            :: t0, tf, delt
    INTEGER                            :: tot_frames, skip_frames, skip_t0_frames, proc_frames
! BOOK KEEPING VARIABLES
    INTEGER            :: i=0, j=0, k=0, ios=0, itemp=0
    REAL(KIND=REAL32)  :: rvtemp(3)=0.0
    CHARACTER(LEN=500) :: inpfile_name, xtcfile_name, outfile_name
    CHARACTER(LEN=100) :: suffix_name
    CHARACTER(LEN=11) :: mydate

    CALL SYSTEM_CLOCK(count_rate=clock_rate)

    mass_deme(1)  = 14.010000
    mass_deme(2)  = 12.010000
    mass_deme(3)  =  1.008000
    mass_deme(4)  =  1.008000
    mass_deme(5)  = 12.010000
    mass_deme(6)  =  1.008000
    mass_deme(7)  =  1.008000
    mass_deme(8)  =  1.008000
    mass_deme(9)  = 12.010000
    mass_deme(10) =  1.008000
    mass_deme(11) =  1.008000
    mass_deme(12) = 12.010000
    mass_deme(13) =  1.008000
    mass_deme(14) =  1.008000
    mass_deme(15) =  1.008000
    mass_deme(16) = 12.010000
    mass_deme(17) =  1.008000
    mass_deme(18) =  1.008000
    mass_deme(19) =  1.008000
    mass_deme(20) = 12.010000
    mass_deme(21) =  1.008000
    mass_deme(22) =  1.008000
    mass_deme(23) = 12.010000
    mass_deme(24) =  1.008000
    mass_deme(25) =  1.008000
    mass_deme(26) = 16.000000
    mass_deme(27) = 12.010000
    mass_deme(28) =  1.008000
    mass_deme(29) =  1.008000
    mass_deme(30) =  1.008000
    sum_mass_deme = SUM(mass_deme)
    one_by_mass_deme = 1.0/sum_mass_deme

    mass_tfsi(1)  = 12.010000
    mass_tfsi(2)  = 32.060000
    mass_tfsi(3)  = 14.010000
    mass_tfsi(4)  = 32.060000
    mass_tfsi(5)  = 12.010000
    mass_tfsi(6)  = 16.000000
    mass_tfsi(7)  = 16.000000
    mass_tfsi(8)  = 16.000000
    mass_tfsi(9)  = 16.000000
    mass_tfsi(10) = 19.000000
    mass_tfsi(11) = 19.000000
    mass_tfsi(12) = 19.000000
    mass_tfsi(13) = 19.000000
    mass_tfsi(14) = 19.000000
    mass_tfsi(15) = 19.000000
    sum_mass_tfsi = SUM(mass_tfsi)
    one_by_mass_tfsi = 1.0/sum_mass_tfsi


!********** READ FROM INPUT FILE **********

    CALL GET_COMMAND_ARGUMENT(1, inpfile_name)
    OPEN(UNIT=10, FILE=inpfile_name, STATUS='OLD', ACTION ='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING INPUT FILE", TRIM(inpfile_name)
    END IF
    PRINT '(A,A)', TRIM(inpfile_name), " OPENED SUCCESSFULLY"

    READ(10,*)
    READ(10,'(A)') xtcfile_name
    READ(10,*)
    READ(10,'(A)') suffix_name
    READ(10,*)
    READ(10,*) num_res_deme, num_res_tfsi
    READ(10,*)
    READ(10,*) ion_charge
    READ(10,*)
    READ(10,*) tot_frames
    READ(10,*)
    READ(10,*) skip_frames
    READ(10,*)
    READ(10,*) skip_t0_frames
    READ(10,*)
    READ(10,*) time_per_frame
    READ(10,*)
    READ(10,*) temperature

    proc_frames = tot_frames - skip_frames
    num_restarts = INT((proc_frames-2)/skip_t0_frames) + 1

    IF(num_res_deme .NE. num_res_tfsi) THEN
          WRITE(*,*) "ERROR) NUMBER OF DEME RESIDUES DOESNOT MATCH THE NUMBER OF TFSI RESIDUES..."
          STOP
    END IF

    tot_atoms_deme = num_res_deme*num_atoms_deme
    tot_atoms_tfsi = num_res_tfsi*num_atoms_tfsi

!********** INITIALIZE MORE VARIABLES **********
    CALL test%open(TRIM(xtcfile_name))

    ios = test%read_next()

    myatoms = test%natoms()
    mystep  = test%step(1)
    mytime  = test%time(1)

    mybox = test%box(1)
    box_size(1) = mybox(1,1)
    box_size(2) = mybox(2,2)
    box_size(3) = mybox(3,3)

    WRITE(*,*) "READ IN ", ios, "TEST FRAMES"
    WRITE(*,*) "NUMBER OF ATOMS: ", myatoms
    WRITE(*,*) "SIMU STEP: ", mystep
    WRITE(*,*) "TIME STEP: ", mytime

    CALL test%close()

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

    ALLOCATE(mj_cat(3,proc_frames))
    ALLOCATE(mj_ani(3,proc_frames))

    mj_cat(:,:) = 0.0 ; mj_ani(:,:) = 0.0
    volume = 0.0

    PRINT *, "***** READING FRAMES ***** : ", proc_frames
    CALL SYSTEM_CLOCK(count=clock_start)
    DO i = 1, proc_frames
        ios = trj%read_next()
        IF (MOD(i,PRINTEVERY)==0) THEN
          WRITE(*,'(A,I0,A,I0,A,F8.2)') "PROCESSING FRAME : ", i, " ;; MD FRAME : ",&
                                       trj%step(1), " ;; MD TIME : ", trj%time(1)
        END IF

        mybox = trj%box(1)
        box_size(1) = mybox(1,1)
        box_size(2) = mybox(2,2)
        box_size(3) = mybox(3,3)
        volume = volume + box_size(1)*box_size(2)*box_size(3)
        
        !********** IDENTIFY COM and ESTIMATE MJ **********
        ! LOOP OVER ALL DEME
        itemp = 0 ; 
        DO j = 1, num_res_deme
            rvtemp(:) = 0.0
            DO k = 1, num_atoms_deme
                rvtemp(:) = rvtemp(:) + trj%x(1,itemp+k)*mass_deme(k)
            END DO
            mj_cat(:,i) = mj_cat(:,i) + ion_charge*rvtemp(:)*one_by_mass_deme 
            itemp = itemp + num_atoms_deme
        END DO

        ! LOOP OVER ALL TFSI
        DO j = 1, num_res_tfsi
            rvtemp(:) = 0.0
            DO k = 1, num_atoms_tfsi
                rvtemp(:) = rvtemp(:) + trj%x(1,itemp+k)*mass_tfsi(k)
            END DO
            mj_ani(:,i) = mj_ani(:,i) - ion_charge*rvtemp(:)*one_by_mass_tfsi
            itemp = itemp + num_atoms_tfsi 
        END DO
    END DO
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)

    ! deallocate raw trajectory
    CALL trj%close()
    DEALLOCATE(trj%frameArray)

    volume = volume/REAL(proc_frames)
    WRITE(*,'(A,F12.3)') " Volume (nm^3) : ", volume

!********** END PROCESS FRAMES **********

!********** COLLECTIVE MSD CALCULATION **********

    ALLOCATE(msdmj_cc(proc_frames-1))
    ALLOCATE(msdmj_aa(proc_frames-1))
    ALLOCATE(msdmj_ca(proc_frames-1))
    ALLOCATE(msd_count(proc_frames-1))

    msdmj_cc = 0.0 ; msdmj_aa = 0.0 ; msdmj_ca = 0.0 ;    
    msd_count = 0
    
    CALL SYSTEM_CLOCK(count=clock_start)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) DEFAULT(NONE) PRIVATE(t0,tf,delt,j,k) &
!$OMP& PRIVATE(delmj_cat,delmj_ani,mj_cat_t0,mj_ani_t0) &
!$OMP& SHARED(time_per_frame,proc_frames,mj_cat,mj_ani) &
!$OMP& REDUCTION(+:msdmj_cc,msdmj_aa,msdmj_ca,msd_count)
    DO t0 = 1, proc_frames-1, skip_t0_frames
        IF (MOD(t0,PRINTEVERY)==0) THEN
            PRINT '(A,F13.2,A)', "***** ESTIMATING COLL MSD AT T0 : ", REAL(t0)*time_per_frame, " ps"
        END IF

        mj_cat_t0(:) = mj_cat(:,t0)
        mj_ani_t0(:) = mj_ani(:,t0)
        
        DO tf = t0+1, proc_frames
            ! delt goes from 1 to proc_frames-1 when t0=1 and
            !                1 to 1             when t0=proc_frames-1
            delt = tf-t0
            msd_count(delt) = msd_count(delt) + 1

            ! calculate delta Mj
            delmj_cat(:) = mj_cat(:,tf) - mj_cat_t0(:)
            delmj_ani(:) = mj_ani(:,tf) - mj_ani_t0(:)

            ! calculate cation - cation
            msdmj_cc(delt) = msdmj_cc(delt) + DOT_PRODUCT(delmj_cat,delmj_cat)

            ! calculate anion - anion
            msdmj_aa(delt) = msdmj_aa(delt) + DOT_PRODUCT(delmj_ani,delmj_ani)

            ! calculate cation - anion
            msdmj_ca(delt) = msdmj_ca(delt) + DOT_PRODUCT(delmj_cat,delmj_ani)

        END DO
    END DO
!$OMP END PARALLEL DO
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)

    DEALLOCATE(mj_cat)
    DEALLOCATE(mj_ani)

    ! These non diagonal terms are symmetric and hence calculated explicitly only half the times
    msdmj_ca(:) = 2.0*msdmj_ca(:)

    ALLOCATE(one_by_msd_count(proc_frames-1))
    one_by_msd_count(:) = 1.0/REAL(msd_count(:))
    msdmj_cc(:) = msdmj_cc(:)*one_by_msd_count(:)
    msdmj_aa(:) = msdmj_aa(:)*one_by_msd_count(:)
    msdmj_ca(:) = msdmj_ca(:)*one_by_msd_count(:)
!********** COLLECTIVE MSD CALCULATION **********


!********** UNIT CONVERSION **********

    mult_factor = prefactor/(volume*temperature)

    msdmj_cc(:) = msdmj_cc(:)*mult_factor
    msdmj_aa(:) = msdmj_aa(:)*mult_factor
    msdmj_ca(:) = msdmj_ca(:)*mult_factor

    ALLOCATE(msdmj(proc_frames-1))
    msdmj(:) = msdmj_cc(:) + msdmj_aa(:) + msdmj_ca(:)

!********** END UNIT CONVERSION **********



!********** PRINT TO FILES **********    

    PRINT *, "***** PRINTING TO FILES *****"

    CALL get_DDMonYY(mydate)

    outfile_name = "ic-cat-cat-"//TRIM(suffix_name)//".dat"
    OPEN(UNIT=20, FILE=outfile_name)
    ! cation - cation
    WRITE(20,'(A,A)')  "# This file was created on   : ", mydate
    WRITE(20,'(A,A)')  "# Created by                 : ", TRIM(source_file)
    WRITE(20,'(A,A)')  "# Source at                  : ", TRIM(source_dir)
    WRITE(20,'(A,A)')  "# Last edited                : ", TRIM(last_edited)
    WRITE(20,'(A,A)')  "# Dump file used             : ", TRIM(xtcfile_name)
    WRITE(20,'(A,I0)') "# Number of frames processed : ", proc_frames
    WRITE(20,'(A,I0)') "# Number of frames total     : ", tot_frames
    WRITE(20,'(A,I0)') "# Number of frames skipped   : ", skip_frames
    WRITE(20,'(A)')    "# Units: x-axis: {ps} ;; y-axis: {ps S/m} "
    WRITE(20,'(A)')    "# time ;; cat-cat-total "
    WRITE(20,'(F13.3,2X,3(E12.6,2X))')  0.0, 0.0
    DO t0 = 1, proc_frames -1
        WRITE(20,'(F13.3,2X,E12.6,2X)')  REAL(t0)*time_per_frame, msdmj_cc(t0)
    END DO
    CLOSE(20)

    ! anion - anion
    outfile_name = "ic-ani-ani-"//TRIM(suffix_name)//".dat"
    OPEN(UNIT=20, FILE=outfile_name)
    WRITE(20,'(A,A)')  "# This file was created on   : ", mydate
    WRITE(20,'(A,A)')  "# Created by                 : ", TRIM(source_file)
    WRITE(20,'(A,A)')  "# Source at                  : ", TRIM(source_dir)
    WRITE(20,'(A,A)')  "# Last edited                : ", TRIM(last_edited)
    WRITE(20,'(A,A)')  "# Dump file used             : ", TRIM(xtcfile_name)
    WRITE(20,'(A,I0)') "# Number of frames processed : ", proc_frames
    WRITE(20,'(A,I0)') "# Number of frames total     : ", tot_frames
    WRITE(20,'(A,I0)') "# Number of frames skipped   : ", skip_frames
    WRITE(20,'(A)')    "# Units: x-axis: {ps} ;; y-axis: {ps S/m} "
    WRITE(20,'(A)')    "# time ;; ani-ani-total "
    WRITE(20,'(F13.3,2X,3(E12.6,2X))')  0.0, 0.0
    DO t0 = 1, proc_frames -1

        WRITE(20,'(F13.3,2X,E12.6,2X)')  REAL(t0)*time_per_frame, msdmj_aa(t0)
    END DO
    CLOSE(20)

    ! cation - anion
    outfile_name = "ic-cat-ani-"//TRIM(suffix_name)//".dat"
    OPEN(UNIT=20, FILE=outfile_name)
    WRITE(20,'(A,A)')  "# This file was created on   : ", mydate
    WRITE(20,'(A,A)')  "# Created by                 : ", TRIM(source_file)
    WRITE(20,'(A,A)')  "# Source at                  : ", TRIM(source_dir)
    WRITE(20,'(A,A)')  "# Last edited                : ", TRIM(last_edited)
    WRITE(20,'(A,A)')  "# Dump file used             : ", TRIM(xtcfile_name)
    WRITE(20,'(A,I0)') "# Number of frames processed : ", proc_frames
    WRITE(20,'(A,I0)') "# Number of frames total     : ", tot_frames
    WRITE(20,'(A,I0)') "# Number of frames skipped   : ", skip_frames
    WRITE(20,'(A)')    "# Units: x-axis: {ps} ;; y-axis: {ps S/m} "
    WRITE(20,'(A)')    "# time ;; cat-ani-total "
    WRITE(20,'(F13.3,2X,E12.6)')  0.0, 0.0
    DO t0 = 1, proc_frames -1
        WRITE(20,'(F13.3,2X,E12.6)')  REAL(t0)*time_per_frame, msdmj_ca(t0)
    END DO
    CLOSE(20)

    ! total
    outfile_name = "ic-total-"//TRIM(suffix_name)//".dat"
    OPEN(UNIT=20, FILE=outfile_name)
    WRITE(20,'(A,A)')  "# This file was created on   : ", mydate
    WRITE(20,'(A,A)')  "# Created by                 : ", TRIM(source_file)
    WRITE(20,'(A,A)')  "# Source at                  : ", TRIM(source_dir)
    WRITE(20,'(A,A)')  "# Last edited                : ", TRIM(last_edited)
    WRITE(20,'(A,A)')  "# Dump file used             : ", TRIM(xtcfile_name)
    WRITE(20,'(A,I0)') "# Number of frames processed : ", proc_frames
    WRITE(20,'(A,I0)') "# Number of frames total     : ", tot_frames
    WRITE(20,'(A,I0)') "# Number of frames skipped   : ", skip_frames
    WRITE(20,'(A)')    "# Units: x-axis: {ps} ;; y-axis: {ps S/m} "
    WRITE(20,'(A)')    "# time ;; total conductivity(t) ;; cat-cat ;; ani-ani ;; cat-ani "
    WRITE(20,'(F13.3,2X,4(E12.6,2X))')  0.0, 0.0, 0.0, 0.0, 0.0
    DO t0 = 1, proc_frames -1
        WRITE(20,'(F13.3,2X,4(E12.6,2X))')  REAL(t0)*time_per_frame, msdmj(t0), msdmj_cc(t0), msdmj_aa(t0), msdmj_ca(t0)
    END DO
    CLOSE(20)

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

END PROGRAM iconduct