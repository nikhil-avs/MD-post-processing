!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM TO CALCULATE THE ALL GOFR (INTER AND AGNOSTIC)

! FIRST STARTED ON 21 JAN 2020
! LAST MODIFIED ON 13 FEB 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! ----------------------------- INDEX FILE FORMAT ---------------------------------------------!!
! ASSUMES THAT ALL ATOMS ARE ACCOUNTED FOR IN THE INDEX FILE PROVIDED
! INDEX GROUP NAME SHOULD HAVE THE FOLLOWING FOUR FIELDS ==
! 1. PARTICLE NAME         (EX: CAT_HC)
! 2. ELEMENT NAME          (EX: H)
! 3. SPECIES NUMBER        (EX: 1)
! 4. PARTICLE MULTIPLICITY (EX: 3)
! EX: [ CAT_HC = H  = 1 = 3 ]
!! ---------------------------------------------------------------------------------------------!!


! CAREFUL WITH UNITS -- GROMACS PREFERS nm, BUT I FEEL angstrom is intuitive
! DO ALL THE PROCESSING IN nm and CONVERT JUST BEFORE PRINTING

! UNITS --
! DISTANCE INPUT (XTC,INP) : NM
! DISTANCE OUTPUT          : ANGSTROM

! PERFORMANCE TESTS --
! 1. INTEGER HISTOGRAM VS REAL HISTOGRAM DOESN'T MAKE MUCH DIFFERENCE
! 2. USING LOCAL VARIABLES SAVES A LOT OF TIME

! PERFORMANCE OPTIONS --
    ! THERE ARE TWO OPTIONS --
    ! 1. READ IN ALL THE TRAJECTORY INTO MEMORY AND THEN ANALYSE -> NEEDS LARGE RAM
    ! 2. READ FRAME BY FRAME AND ANALYSE ON THE FLY -> MIGHT BE SLOW AS READ/WRITE IS A BOTTLE NECK

MODULE myglobal
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    INTEGER, PARAMETER :: maxflags      = 4
    INTEGER, PARAMETER :: maxcharlength = 200
    INTEGER, PARAMETER :: maxgrpname    = 200
    INTEGER, PARAMETER :: maxgrps       = 200
    INTEGER, PARAMETER :: maxlinechars  = 1000
    INTEGER, PARAMETER :: maxindices    = 100000
    INTEGER, PARAMETER :: PRINTEVERY    = 10
    REAL(KIND=REAL64), PARAMETER :: pi = 4.D0*ATAN(1.D0)
    CHARACTER(LEN=20), PARAMETER :: gofr_inter_dir='data-inter/'
    CHARACTER(LEN=20), PARAMETER :: gofr_agn_dir='data-agnostic/'
END MODULE myglobal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! A MODULE TO READ THE INDICES OF A GROUP FROM AN INDEX FILE
 ! LAST EDITED 16 DEC 2019

    !!!!! THE API !!!!!
    
    ! INTEGER :: numgrps
    ! TYPE(ndxgrp), ALLOCATABLE :: mygroups(:)
        ! CHARACTER(LEN=maxgrpname) :: gname
        ! INTEGER :: gnumatoms
        ! INTEGER, ALLOCATABLE :: gindices(:)
    
    ! CALL ndxopenfile("index.ndx")
    ! CALL ndxgetnumgrps(numgrps)
    ! ALLOCATE(mygroups(numgrps))
    ! CALL ndxgetindices(mygroups)
    ! CALL ndxclosefile()
    
    !!!!! THE API !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE readndx
    USE myglobal
    IMPLICIT NONE
    
    TYPE, PUBLIC :: ndxgrp
        CHARACTER(LEN=maxgrpname) :: gname
        INTEGER :: gnumatoms
        INTEGER, ALLOCATABLE :: gindices(:)
    END TYPE ndxgrp

    INTEGER :: grploc(maxgrps)
    INTEGER :: numgrps
    INTEGER :: FUNIT

    PRIVATE :: grploc, numgrps, FUNIT, skiplines
    PUBLIC  :: ndxopenfile, ndxgetnumgrps, ndxgetindices, ndxclosefile

    CONTAINS

    SUBROUTINE ndxopenfile(fname)
        CHARACTER(LEN=*), INTENT(IN) :: fname
        INTEGER :: ios
        OPEN(NEWUNIT=FUNIT,FILE=fname,IOSTAT=ios,STATUS='OLD',ACTION='READ')   ! --> ADD FILE CHECK DETAILS
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR OPENING THE FILE..", TRIM(fname), ".. ERROR CODE: ", ios
            STOP
        ENDIF
    END SUBROUTINE ndxopenfile

    ! This should be used to detect the number of groups and also the corresponding positions
    SUBROUTINE ndxgetnumgrps(mynumgrps)
        INTEGER, INTENT(OUT) :: mynumgrps
        CHARACTER(LEN=maxlinechars) :: thisline
        INTEGER :: i, ios, temp1, lineno

        ios = 0 ; numgrps = 0 ; lineno = 0 ; grploc = -1

        DO
            READ(FUNIT,'(A)',IOSTAT=ios) thisline
            IF(IS_IOSTAT_END(ios)) THEN
                EXIT
            END IF

            lineno = lineno + 1
            thisline = ADJUSTL(thisline)
            
            IF(thisline(1:1)=='[') THEN
                temp1 = LEN_TRIM(thisline)
                IF(thisline(temp1:temp1) == ']') THEN
                    numgrps = numgrps + 1
                    grploc(numgrps) = lineno
                    ! WRITE(*,*) "found a group : [", thisline(2:temp1-1), "] at location : ", grploc(numgrps)
                ELSE
                    WRITE(*,*) "found an improper group heading "
                END IF
            END IF
        END DO
        REWIND(FUNIT)
        grploc(numgrps+1) = lineno + 1
        mynumgrps = numgrps
    END SUBROUTINE ndxgetnumgrps

    ! This should populate the ndxgrp data structure  
    SUBROUTINE ndxgetindices(mygroups)
        TYPE(ndxgrp), INTENT(INOUT) :: mygroups(:)
        CHARACTER(LEN=maxlinechars) :: thisline
        INTEGER :: i, j, ios, lineno, temp1
        INTEGER :: prevloc, intcount
        INTEGER :: largeintarr(maxlinechars), indices(maxindices)

        ios = 0 ; lineno = 0 ; indices = -1 ; largeintarr = -1

        REWIND(FUNIT)
        DO i = 1, numgrps
            CALL skiplines(FUNIT,grploc(i)-lineno-1)
            READ(FUNIT,'(A)',IOSTAT=ios) thisline
            lineno = grploc(i)
            thisline = ADJUSTL(thisline)
            temp1 = LEN_TRIM(thisline)
            mygroups(i)%gname = thisline(2:temp1-1)

            ! WRITE(*,*) "Reading indices of group ", TRIM(mygroups(i)%gname)

            prevloc = 0
            DO WHILE(lineno .LT. grploc(i+1)-1)
                READ(FUNIT,'(A)',IOSTAT=ios) thisline
                READ(thisline,*,IOSTAT=ios) largeintarr
                lineno = lineno + 1
                IF(.NOT. IS_IOSTAT_END(ios)) THEN
                    WRITE(*,*) "Illegal line detected in index file : ", lineno
                    STOP
                END IF

                DO j = 1, 1000
                    IF(largeintarr(j) == -1) THEN
                        intcount = j-1
                        EXIT
                    END IF
                END DO

                IF(intcount == 0) EXIT     !EXIT CASE 2

                indices(prevloc+1:prevloc+intcount) = largeintarr(1:intcount)
                prevloc = prevloc + intcount
                largeintarr(1:intcount) = -1
            END DO
            mygroups(i)%gnumatoms = prevloc
            ALLOCATE(mygroups(i)%gindices(prevloc))
            mygroups(i)%gindices(1:prevloc) = indices(1:prevloc)
            indices(1:prevloc) = -1
            ! WRITE(*,*) "finished with group ", i, "on line ", lineno, " with atoms ", prevloc
        END DO
    END SUBROUTINE ndxgetindices

    SUBROUTINE skiplines(sunit, nskip)
        INTEGER, INTENT(IN) :: sunit, nskip
        INTEGER :: i
        DO i = 1, nskip
            READ(sunit,*)
        END DO
    END SUBROUTINE skiplines

    SUBROUTINE ndxclosefile()
        CLOSE(FUNIT)   ! --> ADD FILE CHECK DETAILS
    END SUBROUTINE ndxclosefile
END MODULE readndx

PROGRAM gofrall
    USE myglobal
    USE readndx
    USE gmxfort_trajectory 
    IMPLICIT NONE

    !!!!---------------------------------------- ALPHA GROUP RELATED VARIABLES ------------------------------------!!!!
    ! Assign a molecule number to an atom 
    !(For ex: there are two HBs in cation, we also have to know which two HBs belong to cation1, cation2 etc..)
    ! Now we know the species and the molecule number in that species to which a certain atom belongs to 
    TYPE :: moltype
        INTEGER, ALLOCATABLE :: molnum(:)
    END TYPE moltype
    INTEGER                                :: numgrps, totgrpatoms
    TYPE(ndxgrp), ALLOCATABLE              :: alphagrps(:)
    TYPE(moltype), ALLOCATABLE             :: alpha_moltype(:)
    CHARACTER(LEN=2), ALLOCATABLE          :: alphagrp_elname(:)
    INTEGER, ALLOCATABLE                   :: alphagrp_species(:), alphagrp_mult(:)
    CHARACTER(LEN=500)                     :: ndxfile_name='alphas-bmim-cl-ntf2-0p50.ndx'

    ! INDEX GROUP NAME RELATED
    INTEGER :: alphagrp_namesplit(3), count_split
 
    !!!!--------------------------------------------- I/O RELATED VARIABLES ---------------------------------------!!!!
    TYPE(Trajectory)   :: trj
    INTEGER            :: proc_frames=4, skip_frames=10
    CHARACTER(LEN=500)   :: inpfile_name, suffix_name
    CHARACTER(LEN=500)   :: xtcfile_name = '/home/nikhil/projects/8_extend_BIL/welton-structures/&
                                          &mb-v1/bmim-cl-ntf2/0.50/3-nvt/nvt.xtc'
    INTEGER            :: inpunit

    ! ASYNCHRONOUS I/O    
    INTEGER, ALLOCATABLE            :: gofr_unit(:,:)
    CHARACTER(LEN=500), ALLOCATABLE :: gofr_fname(:,:)


    !!!!----------------------------------------- TRAJECTORY RELATED VARIABLES ------------------------------------!!!!
    TYPE :: coordtype
        REAL(KIND=REAL32), ALLOCATABLE :: r(:,:,:)
    END TYPE coordtype
    TYPE(coordtype), ALLOCATABLE :: coord(:)
    REAL(KIND=REAL32)            :: box_size(3), mybox(3,3), volume

    !!!!--------------------------------------------- GOFR RELATED VARIABLES --------------------------------------!!!!
    REAL(KIND=REAL32), ALLOCATABLE :: gofr_inter(:,:,:), gofr_agn(:,:,:), cn_inter(:,:,:), cn_agn(:,:,:), dist(:)
    REAL(KIND=REAL32) :: rcut, rbin_size=0.01
    INTEGER           :: num_rbin, max_rbin
    LOGICAL           :: select_rcut

    REAL(KIND=REAL32), ALLOCATABLE :: tempra(:,:), temprb(:,:)
    REAL(KIND=REAL32) :: rcut_sq, rlow, rhigh, dvol
    REAL(KIND=REAL32) :: delrv(3), delr, tempr(3)
    INTEGER           :: rbin, max_compare_rbin
    LOGICAL           :: select_inter

    !!!!---------------------------------------- CLOCK RELATED VARIABLES ------------------------------------------!!!!
    INTEGER(KIND=8)   :: clock_start, clock_end, clock_rate, clock_max

    !!!!----------------------------------------- BOOK KEEPING VARIABLES ------------------------------------------!!!!
    INTEGER :: i, j, k, l, alpha, beta
    INTEGER :: ios


    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    select_rcut = .FALSE.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!***************************************** READ IN INDICES OF ALPHA GROUPS ******************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*)
    CALL GET_COMMAND_ARGUMENT(1, inpfile_name)
    OPEN(NEWUNIT=inpunit, FILE=inpfile_name, STATUS='OLD', ACTION ='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        WRITE(*,*) "ERROR IN OPENENING INPUT FILE : ", TRIM(inpfile_name)
        STOP
    END IF
    WRITE(*,'(A,A,A)') "INPUT) ", TRIM(inpfile_name), " OPENED SUCCESSFULLY"

    READ(inpunit,*)
    READ(inpunit,'(A)') xtcfile_name
    WRITE(*,'(A,A)')    "INPUT) XTCFILE        : ", TRIM(xtcfile_name)
    READ(inpunit,*)
    READ(inpunit,'(A)') ndxfile_name
    WRITE(*,'(A,A)')    "INPUT) NDXFILE        : ", TRIM(ndxfile_name)
    READ(inpunit,*)
    READ(inpunit,'(A)') suffix_name
    WRITE(*,'(A,A)')    "INPUT) SUFFIX NAME    : ", TRIM(suffix_name)
    READ(inpunit,*)
    READ(inpunit,*) rbin_size
    WRITE(*,'(A,F6.3)') "INPUT) GOFR BIN SIZE  : ", rbin_size
    READ(inpunit,*)
    READ(inpunit,*) rcut
    WRITE(*,'(A,F6.3)') "INPUT) GOFR RCUT      : ", rcut
    READ(inpunit,*)
    READ(inpunit,*) skip_frames
    WRITE(*,'(A,I0)')   "INPUT) SKIP FRAMES    : ", skip_frames
    READ(inpunit,*)
    READ(inpunit,*) proc_frames
    WRITE(*,'(A,I0)')   "INPUT) PROCESS FRAMES : ", proc_frames
    IF(rcut <= rbin_size) THEN
        WRITE(*,'(A,F6.3)') "INPUT) WILL SET rcut TO 0.5*BOX_SIZE "
        select_rcut = .TRUE.
    END IF
    WRITE(*,*)

    !!!!********************************************************************************************************************!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!***************************************** READ IN INDICES OF ALPHA GROUPS ******************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL ndxopenfile(ndxfile_name)
    CALL ndxgetnumgrps(numgrps)
    ALLOCATE(alphagrps(numgrps),alphagrp_elname(numgrps),alphagrp_species(numgrps),alphagrp_mult(numgrps))
    ALLOCATE(alpha_moltype(numgrps))
    CALL ndxgetindices(alphagrps)
    CALL ndxclosefile()

    totgrpatoms = 0 ; ios = 0
    WRITE(*,'(A,I0)')   "----> TOTAL GROUPS : ", numgrps
    DO i = 1, numgrps
        ! GET ELEMENT NAME, MULTIPLICITY ETC FROM INDEX GROUP NAME
        count_split = 0
        DO j = 1, maxgrpname
            IF(alphagrps(i)%gname(j:j) == '=') THEN
                count_split = count_split + 1
                alphagrp_namesplit(count_split) = j
                IF(count_split == 3) EXIT
            END IF
        END DO
        IF (count_split .NE. 3) THEN
            WRITE(*,'(A,A)') TRIM(alphagrps(i)%gname), " -- NDX GROUP DOESN'T PROPER NAMING.. STOPPING"
            STOP
        END IF

        alphagrp_elname(i) = TRIM(ADJUSTL(alphagrps(i)%gname(alphagrp_namesplit(1)+1:alphagrp_namesplit(2)-1)))
        READ(alphagrps(i)%gname(alphagrp_namesplit(2)+1:alphagrp_namesplit(3)-1),*) alphagrp_species(i)
        READ(alphagrps(i)%gname(alphagrp_namesplit(3)+1:),*) alphagrp_mult(i)
        alphagrps(i)%gname(alphagrp_namesplit(1):) = ' '

        ! GET MOLECULE NUMBER FROM MULTIPLICITY
        IF( MOD(alphagrps(i)%gnumatoms,alphagrp_mult(i)) .NE. 0) THEN
            WRITE(*,*) "ERROR) MULTIPLICITY AND NUMATOMS DON'T MATCH FOR THE INDEX GROUP : ", &
                alphagrp_mult(i), alphagrps(i)%gnumatoms, TRIM(alphagrps(i)%gname)
            STOP
        END IF
        ALLOCATE(alpha_moltype(i)%molnum(alphagrps(i)%gnumatoms))
        DO j = 1, alphagrps(i)%gnumatoms
            alpha_moltype(i)%molnum(j) = INT(REAL(j-1)/REAL(alphagrp_mult(i))) + 1
        END DO

        ! WRITE NDX FILE INFORMATION TO OUTPUT
        WRITE(*,'(A)') "****************************************"
        WRITE(*,'(I3,A,A)') i, " GROUP         : ", TRIM(alphagrps(i)%gname)
        WRITE(*,'(A,A)')    "--> ELEMENT       : ", alphagrp_elname(i)
        WRITE(*,'(A,I0)')   "--> SPECIES       : ", alphagrp_species(i)
        WRITE(*,'(A,I0)')   "--> MULTIPLICITY  : ", alphagrp_mult(i)
        WRITE(*,'(A,I0)')   "--> NUM-MOLS      : ", alpha_moltype(i)%molnum(alphagrps(i)%gnumatoms)
        WRITE(*,'(A,I0)')   "--> NUMATOMS      : ", alphagrps(i)%gnumatoms
        WRITE(*,'(A,I0)')   "--> ATOM1    : ", alphagrps(i)%gindices(1)
        WRITE(*,'(A,I0)')   "--> ATOMLAST : ", alphagrps(i)%gindices(alphagrps(i)%gnumatoms)
        WRITE(*,'(A)') "****************************************"
        totgrpatoms = totgrpatoms + alphagrps(i)%gnumatoms
    END DO
    WRITE(*,'(A,I0)')   "----> TOTAL ATOMS : ", totgrpatoms
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!********************************************** READ IN TRAJECTORY **************************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CALL trj%open(TRIM(xtcfile_name))

    WRITE(*,*) "***** OPENED TRAJECTORY TO READ FRAMES *****"
    WRITE(*,*) "TOTAL NUMBER OF FRAMES PRESENT IN TRAJECTORY: ", trj%nframes
    WRITE(*,*) "TOTAL NUMBER OF ATOMS IN A FRAME            : ", trj%natoms()

    IF(totgrpatoms .NE. trj%natoms()) THEN
        WRITE(*,*) "WARNING! THE NUMBER OF ATOMS IN INDEX FILE AND TRAJECTORY DO NOT MATCH"
    END IF

    IF(trj%nframes .LT. skip_frames+proc_frames) THEN
        WRITE(*,*) "TOTAL NUMBER OF FRAMES IN THE TRAJECTORY LESS THAN PROC_FRAMES+SKIP_FRAMES "
        STOP
    END IF

    WRITE(*,*) 

    
    !********** SKIP FRAMES **********
    IF(skip_frames > 0) THEN
        CALL SYSTEM_CLOCK(count=clock_start)
        WRITE(*,*) "***** SKIPPING FRAMES ***** : ", skip_frames
        
        ios = trj%read_next(skip_frames)
        mybox = trj%box(1)
        box_size(1) = mybox(1,1)
        box_size(2) = mybox(2,2)
        box_size(3) = mybox(3,3)
        
        CALL SYSTEM_CLOCK(count=clock_end)
        WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
        WRITE(*,*)
    END IF

    !********** PROCESS FRAMES **********
    WRITE(*,*) "***** READING FRAMES ***** : ", proc_frames
    CALL SYSTEM_CLOCK(count=clock_start)

    ALLOCATE(coord(numgrps))
    DO i = 1, numgrps
        ALLOCATE(coord(i)%r(3,alphagrps(i)%gnumatoms,proc_frames))
    END DO
    
    DO i = 1, proc_frames
        ios = trj%read_next()
        IF (MOD(i,PRINTEVERY)==0) THEN
          WRITE(*,'(A,I0,A,I0,A,F8.2)') "PROCESSING FRAME : ", i, " ;; MD FRAME : ",&
                                       trj%step(1), " ;; MD TIME : ", trj%time(1)
        END IF

        DO j = 1, numgrps
            DO k = 1, alphagrps(j)%gnumatoms
                coord(j)%r(:,k,i) = trj%x(1,alphagrps(j)%gindices(k))
            END DO
        END DO
    END DO

    ! GET BOX SIZE IF NO FRAMES WERE SKIPPED
    IF(skip_frames .LE. 0) THEN
        mybox = trj%box(1)
        box_size(1) = mybox(1,1)
        box_size(2) = mybox(2,2)
        box_size(3) = mybox(3,3)
    END IF
    
    CALL trj%close()
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!************************************************ CALCULATE G(R) ****************************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) "***** CALCULATING G(R) ***** : "
    CALL SYSTEM_CLOCK(count=clock_start)

    ! WRITE GOFR SETTINGS TO OUTPUT    
    WRITE(*,'(A,3(F12.6,2X))') "BOX SIZE: ", box_size(:)
    
    volume = box_size(1)*box_size(2)*box_size(3)
    WRITE(*,'(A,F12.4)') "VOLUME  : ", volume

    IF(select_rcut) THEN    
        rcut = SUM(box_size(:))/6.0
    END IF
    rcut = rcut - rbin_size
    rcut_sq = rcut*rcut 
    WRITE(*,'(A,F12.6)') "RCUT(nm): ", rcut
    WRITE(*,'(A,F12.6)') "RBIN(nm): ", rbin_size
    WRITE(*,'(A,F12.6)') "RCUT SQ : ", rcut_sq
    
    max_rbin = INT(rcut/rbin_size) + 1
    num_rbin = max_rbin
    WRITE(*,'(A,I0)') "NUM BINS: ", num_rbin
    IF (num_rbin < 4) THEN
        WRITE(*,*) "TOO SPARCE A GRID FOR G(R).. STOPPING"
        STOP
    END IF
    
    ALLOCATE(gofr_agn(num_rbin,numgrps,numgrps)) ; gofr_agn(:,:,:) = 0.0
    ALLOCATE(gofr_inter(num_rbin,numgrps,numgrps)) ; gofr_inter(:,:,:) = 0.0
    select_inter = .FALSE.

    ! MUCH OF THE SHARED DATA IS ALSO MADE PRIVATE TO AVOID RACE CONDITIONS
!$OMP PARALLEL DO SCHEDULE(STATIC, 1) DEFAULT(NONE) PRIVATE(i,alpha,beta,j,k) &
!$OMP& PRIVATE(tempra,temprb,tempr,delrv,delr,rbin) &
!$OMP& FIRSTPRIVATE(proc_frames,numgrps,rcut_sq,rbin_size,box_size) &
!$OMP& FIRSTPRIVATE(alphagrps,alphagrp_species,alpha_moltype,select_inter) &
!$OMP& SHARED(coord) &
!$OMP& REDUCTION(+:gofr_agn,gofr_inter)
    DO i = 1, proc_frames
        IF (MOD(i,PRINTEVERY)==0) THEN
            WRITE(*,'(A,I0)') "PROCESSING FRAME : ", i
        END IF

        ! g_ab(r)
        DO alpha = 1, numgrps-1

            ALLOCATE(tempra(3,alphagrps(alpha)%gnumatoms))
            tempra(:,:) =  coord(alpha)%r(:,:,i)
            
            DO beta = alpha+1, numgrps
                ALLOCATE(temprb(3,alphagrps(beta)%gnumatoms))
                temprb(:,:) = coord(beta)%r(:,:,i)
                IF(alphagrp_species(alpha) == alphagrp_species(beta)) select_inter = .TRUE.

                DO j = 1, alphagrps(alpha)%gnumatoms
                    tempr(:) = tempra(:,j)
                    DO k = 1, alphagrps(beta)%gnumatoms
                        delrv(:) = temprb(:,k) - tempr(:)
                        delrv(:) = delrv(:) - box_size(:)*NINT(delrv(:)/box_size(:))
                        delr = DOT_PRODUCT(delrv,delrv)
                        IF(delr < rcut_sq) THEN
                            delr = SQRT(delr)
                            rbin = INT(delr/rbin_size) + 1
                            gofr_agn(rbin,alpha,beta) = gofr_agn(rbin,alpha,beta) + 1.0
                            IF(select_inter) THEN
                                IF(alpha_moltype(alpha)%molnum(j) == alpha_moltype(beta)%molnum(k))THEN
                                    CYCLE
                                END IF
                            END IF
                            gofr_inter(rbin,alpha,beta) = gofr_inter(rbin,alpha,beta) + 1.0
                        END IF
                    END DO
                END DO
                DEALLOCATE(temprb)
                select_inter = .FALSE.
            END DO
            DEALLOCATE(tempra)
        END DO

        ! g_aa(r)
        DO alpha = 1, numgrps
            ALLOCATE(tempra(3,alphagrps(alpha)%gnumatoms))
            tempra(:,:) =  coord(alpha)%r(:,:,i)
            DO j = 1, alphagrps(alpha)%gnumatoms-1
                tempr(:) = tempra(:,j)
                DO k = j+1, alphagrps(alpha)%gnumatoms
                    delrv(:) = tempra(:,k) - tempr(:)
                    delrv(:) = delrv(:) - box_size(:)*NINT(delrv(:)/box_size(:))
                    delr = DOT_PRODUCT(delrv,delrv)
                    IF(delr < rcut_sq) THEN
                        delr = SQRT(delr)
                        rbin = INT(delr/rbin_size) + 1
                        gofr_agn(rbin,alpha,alpha) = gofr_agn(rbin,alpha,alpha) + 2.0
                        IF(alpha_moltype(alpha)%molnum(j) == alpha_moltype(alpha)%molnum(k))THEN
                            CYCLE
                        END IF
                        gofr_inter(rbin,alpha,alpha) = gofr_inter(rbin,alpha,alpha) + 2.0
                    END IF
                END DO
            END DO
            DEALLOCATE(tempra)
        END DO
    END DO
!$OMP END PARALLEL DO

    DEALLOCATE(coord)

    !! CONVERT GOFR TO NUMBER DENSITY || CALCULATE COORDINATION NUMBER || NORMALIZE GOFR !!
    ALLOCATE(cn_agn(num_rbin,numgrps,numgrps)) ; cn_agn(:,:,:) = 0.0        
    ALLOCATE(cn_inter(num_rbin,numgrps,numgrps)) ; cn_inter(:,:,:) = 0.0     
    select_inter = .FALSE.

    DO alpha = 1, numgrps-1
        DO beta = alpha+1, numgrps
            IF(alphagrp_species(alpha) == alphagrp_species(beta)) select_inter = .TRUE.

            gofr_agn(:,beta,alpha)   = gofr_agn(:,alpha,beta)
            gofr_inter(:,beta,alpha) = gofr_inter(:,alpha,beta)

            ! ENSEMBLE AVERAGING (AVERAGE OVER PARTICLE ORIGINS AND TIME ORIGINS)
            gofr_agn(:,alpha,beta) = gofr_agn(:,alpha,beta)/(REAL(proc_frames)*REAL(alphagrps(alpha)%gnumatoms))
            gofr_agn(:,beta,alpha) = gofr_agn(:,beta,alpha)/(REAL(proc_frames)*REAL(alphagrps(beta)%gnumatoms))

            gofr_inter(:,alpha,beta) = gofr_inter(:,alpha,beta)/(REAL(proc_frames)*REAL(alphagrps(alpha)%gnumatoms))
            gofr_inter(:,beta,alpha) = gofr_inter(:,beta,alpha)/(REAL(proc_frames)*REAL(alphagrps(beta)%gnumatoms))

            ! CALCULATE THE COORDINATION NUMBER
            cn_agn(1,alpha,beta) = gofr_agn(1,alpha,beta) ;  cn_agn(1,beta,alpha) = gofr_agn(1,beta,alpha)
            cn_inter(1,alpha,beta) = gofr_inter(1,alpha,beta) ; cn_inter(1,beta,alpha) = gofr_inter(1,beta,alpha)
            DO i = 2, num_rbin
                cn_agn(i,alpha,beta) = cn_agn(i-1,alpha,beta) + gofr_agn(i,alpha,beta)
                cn_agn(i,beta,alpha) = cn_agn(i-1,beta,alpha) + gofr_agn(i,beta,alpha)

                cn_inter(i,alpha,beta) = cn_inter(i-1,alpha,beta) + gofr_inter(i,alpha,beta)
                cn_inter(i,beta,alpha) = cn_inter(i-1,beta,alpha) + gofr_inter(i,beta,alpha)
            END DO

            DO i = 1, num_rbin
                rlow = REAL(i-1)*rbin_size
                rhigh = rlow + rbin_size
                dvol = (rhigh**3 - rlow**3)
                gofr_agn(i,alpha,beta) = (0.75*volume)*gofr_agn(i,alpha,beta)/(pi*dvol*alphagrps(beta)%gnumatoms)
                gofr_agn(i,beta,alpha) = (0.75*volume)*gofr_agn(i,beta,alpha)/(pi*dvol*alphagrps(alpha)%gnumatoms)

                IF(select_inter) THEN
                    gofr_inter(i,alpha,beta) = (0.75*volume)*gofr_inter(i,alpha,beta)/ &
                                               (pi*dvol*(alphagrps(beta)%gnumatoms-alphagrp_mult(beta)))
                    gofr_inter(i,beta,alpha) = (0.75*volume)*gofr_inter(i,beta,alpha)/ &
                                               (pi*dvol*(alphagrps(alpha)%gnumatoms-alphagrp_mult(alpha)))
                ELSE
                    gofr_inter(i,alpha,beta) = (0.75*volume)*gofr_inter(i,alpha,beta)/(pi*dvol*alphagrps(beta)%gnumatoms)
                    gofr_inter(i,beta,alpha) = (0.75*volume)*gofr_inter(i,beta,alpha)/(pi*dvol*alphagrps(alpha)%gnumatoms)
                END IF
            END DO
            select_inter = .FALSE.
        END DO
    END DO
    
    !! CONVERT GOFR TO NUMBER DENSITY || CALCULATE COORDINATION NUMBER || NORMALIZE GOFR !!
    DO alpha = 1, numgrps
        gofr_agn(:,alpha,alpha) = gofr_agn(:,alpha,alpha)/(REAL(proc_frames)*REAL(alphagrps(alpha)%gnumatoms))
        gofr_inter(:,alpha,alpha) = gofr_inter(:,alpha,alpha)/(REAL(proc_frames)*REAL(alphagrps(alpha)%gnumatoms))

        cn_agn(1,alpha,alpha) = gofr_agn(1,alpha,alpha) 
        cn_inter(1,alpha,alpha) = gofr_inter(1,alpha,alpha)
        DO i = 2, num_rbin
            cn_agn(i,alpha,alpha) = cn_agn(i-1,alpha,alpha) + gofr_agn(i,alpha,alpha)
            cn_inter(i,alpha,alpha) = cn_inter(i-1,alpha,alpha) + gofr_inter(i,alpha,alpha)
        END DO

        DO i = 1, num_rbin
            rlow = REAL(i-1)*rbin_size
            rhigh = rlow + rbin_size
            dvol = (rhigh**3 - rlow**3)
            ! NUMBER OF TARGET ATOMS ARE N-1 (BECAUSE WE DON'T COUNT SELF CORRELATION)
            gofr_agn(i,alpha,alpha) = (0.75*volume)*gofr_agn(i,alpha,alpha)/(pi*dvol*(alphagrps(alpha)%gnumatoms-1))
            ! NUMBER OF TARGET ATOMS ARE N-MULT (BECAUSE WE DON'T COUNT ATOMS FROM THE SAME MOLECULE)
            gofr_inter(i,alpha,alpha) = (0.75*volume)*gofr_inter(i,alpha,alpha)/ &
                                        (pi*dvol*(alphagrps(alpha)%gnumatoms-alphagrp_mult(alpha)))
        END DO
    END DO

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    ALLOCATE(dist(num_rbin))
    DO i = 1, num_rbin
        dist(i) = REAL(i-1)*rbin_size + 0.50*rbin_size
    END DO

    !!---------------------->>> WRITE THE GOFR AND CN DATA TO FILES <<<----------------------!!
    WRITE(*,'(A)') " **** WRITING GOFR TO FILES SYNCHRONOUSLY **** "
    CALL SYSTEM_CLOCK(count=clock_start)

    ALLOCATE(gofr_unit(numgrps,numgrps),gofr_fname(numgrps,numgrps))
    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            gofr_fname(alpha,beta) = TRIM(gofr_agn_dir)//'gofr-agn-'//TRIM(ADJUSTL(alphagrps(alpha)%gname))//&
                                     '-'//TRIM(ADJUSTL(alphagrps(beta)%gname))//'-'//TRIM(suffix_name)//'.dat'
            OPEN(NEWUNIT=gofr_unit(alpha,beta), FILE=gofr_fname(alpha,beta))
            DO i = 1, num_rbin
                ! WRITE TO FILES WITH NM TO ANGSTROM CONVERSION
                WRITE(gofr_unit(alpha,beta),'(E15.7,2X,E15.7,2X,E15.7)') &
                                 dist(i)*10.0, gofr_agn(i,alpha,beta), cn_agn(i,alpha,beta)
            END DO
            CLOSE(gofr_unit(alpha,beta))
        END DO
    END DO

    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            gofr_fname(alpha,beta) = TRIM(gofr_inter_dir)//'gofr-inter-'//TRIM(ADJUSTL(alphagrps(alpha)%gname))//&
                                     '-'//TRIM(ADJUSTL(alphagrps(beta)%gname))//'-'//TRIM(suffix_name)//'.dat'
            OPEN(NEWUNIT=gofr_unit(alpha,beta), FILE=gofr_fname(alpha,beta))
            DO i = 1, num_rbin
                ! WRITE TO FILES WITH NM TO ANGSTROM CONVERSION
                WRITE(gofr_unit(alpha,beta),'(E15.7,2X,E15.7,2X,E15.7)') &
                                 dist(i)*10.0, gofr_inter(i,alpha,beta), cn_inter(i,alpha,beta)
            END DO
        END DO
    END DO

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    !!---------------------->>> WRITE SUMMARY TO LOG FILE <<<----------------------!!
    WRITE(*,'(A)') " **** WRITING SUMMARY **** "
    CALL SYSTEM_CLOCK(count=clock_start)

    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            gofr_fname(alpha,beta) = 'gofr-*-'//TRIM(ADJUSTL(alphagrps(alpha)%gname))//&
                                     '-'//TRIM(ADJUSTL(alphagrps(beta)%gname))//'.dat'
            WRITE(*,'(A)') "****************************************"
            WRITE(*,'(A,A)')    " FILE : ", TRIM(gofr_fname(alpha,beta))
            WRITE(*,'(A,F9.2)') "--> MAX COORD NUM AGN   : ", MAXVAL(cn_agn(:,alpha,beta))
            WRITE(*,'(A,F9.2)') "--> MAX COORD NUM INTER : ", MAXVAL(cn_inter(:,alpha,beta))
            WRITE(*,'(A,F9.6)') "--> LAST GOFR VAL AGN   : ", gofr_agn(num_rbin-1,alpha,beta)
            WRITE(*,'(A,F9.6)') "--> LAST GOFR VAL INTER : ", gofr_inter(num_rbin-1,alpha,beta)
            IF(ABS(gofr_agn(num_rbin-1,alpha,beta)-1.0) .GT. 0.05) THEN
                WRITE(*,*) "WARNING! LAST BINS OF THIS GOFR_AGN DEVIATE SIGNIFICATNLY FROM 1"
            END IF
            IF(ABS(gofr_inter(num_rbin-1,alpha,beta)-1.0) .GT. 0.05) THEN
                WRITE(*,*) "WARNING! LAST BINS OF THIS GOFR_INTER DEVIATE SIGNIFICATNLY FROM 1"
            END IF
            WRITE(*,'(A)') "****************************************"
        END DO
    END DO

    DEALLOCATE(gofr_unit,gofr_fname)

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    WRITE(*,'(A)') " **** ALL THE BEST :) **** "

END PROGRAM gofrall