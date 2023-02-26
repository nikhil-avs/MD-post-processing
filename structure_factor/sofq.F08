!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM TO CALCULATE THE S(Q) [MARGULIS JCP 2011 FORMULA] 

! FIRST STARTED ON 08 FEB 2020
! LAST MODIFIED ON 14 FEB 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ASSUMES THAT ALL ATOMS ARE ACCOUNTED FOR IN THE INDEX FILE PROVIDED

!! ----------------------------- INDEX FILE FORMAT ---------------------------------------------!!
! ASSUMES THAT ALL ATOMS ARE ACCOUNTED FOR IN THE INDEX FILE PROVIDED
! INDEX GROUP NAME SHOULD HAVE THE FOLLOWING FOUR FIELDS ==
! 1. PARTICLE NAME         (EX: CAT_HC)
! 2. ELEMENT NAME          (EX: H)
! 3. SPECIES NUMBER        (EX: 1) --> used to recognize which molecule they belong to
! 4. PARTICLE MULTIPLICITY (EX: 3)
! EX: [ CAT_HC = H  = 1 = 3 ]
!! ---------------------------------------------------------------------------------------------!!

! UNITS --
! DISTANCE INPUT (GOFR,INP): ANGSTROM
! DISTANCE OUTPUT          : ANGSTROM
! Q (RECIPROCAL VECTOR)    : 1/ANGSTROM (DEFINED AS 2*PI/R ## NOT 1/R ##)

MODULE myglobal
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    INTEGER, PARAMETER :: maxflags      = 4
    INTEGER, PARAMETER :: maxcharlength = 200
    INTEGER, PARAMETER :: maxgrpname    = 200
    INTEGER, PARAMETER :: maxgrps       = 200
    INTEGER, PARAMETER :: maxlinechars  = 1000
    INTEGER, PARAMETER :: maxindices    = 100000
    INTEGER, PARAMETER :: PRINTEVERY    = 2
    REAL(KIND=REAL64), PARAMETER :: pi = 4.D0*ATAN(1.D0)
    REAL(KIND=REAL64), PARAMETER :: one_by_2pi = 0.5D0/pi
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

PROGRAM sofq
    USE myglobal
    USE readndx
    IMPLICIT NONE

    !!!!---------------------------------------- ALPHA GROUP RELATED VARIABLES ------------------------------------!!!!
    TYPE :: moltype
        INTEGER, ALLOCATABLE :: molnum(:)
    END TYPE moltype
    INTEGER                                :: numgrps, totgrpatoms
    TYPE(ndxgrp), ALLOCATABLE              :: alphagrps(:)
    TYPE(moltype), ALLOCATABLE             :: alpha_moltype(:)
    CHARACTER(LEN=2), ALLOCATABLE          :: alphagrp_elname(:)
    INTEGER, ALLOCATABLE                   :: alphagrp_species(:), alphagrp_mult(:)
    CHARACTER(LEN=500)                     :: ndxfile_name

    TYPE :: grptype
        INTEGER, ALLOCATABLE :: grpnum(:)
    END TYPE grptype
    INTEGER                                :: numspecies
    INTEGER, ALLOCATABLE                   :: num_grp_species(:) !--> number of groups in each species
    TYPE(grptype), ALLOCATABLE             :: grp_of_species(:) !--> store the groups that belong a given species
    CHARACTER(LEN=maxgrpname), ALLOCATABLE :: species_name(:)

    ! INDEX GROUP NAME RELATED
    INTEGER :: alphagrp_namesplit(3), count_split
 
    !!!!--------------------------------------------- I/O RELATED VARIABLES ---------------------------------------!!!!
    CHARACTER(LEN=500) :: inpfile_name
    CHARACTER(LEN=100) :: suffix
    INTEGER            :: inpunit
    CHARACTER(LEN=500) :: sabq_dir='sabq', sabs_dir='sabs', ff_dir='ffs', qpart_dir='sqpart', spart_dir='sspart'


    !!!!----------------------------------------- TRAJECTORY RELATED VARIABLES ------------------------------------!!!!
    REAL(KIND=REAL32)               :: box_size, volume

    !!!!--------------------------------------------- GOFR RELATED VARIABLES --------------------------------------!!!!
    REAL(KIND=REAL32), ALLOCATABLE :: gofr(:,:,:), cn(:,:,:), dist(:)
    REAL(KIND=REAL32)  :: rcut, rbin_size, rmax
    INTEGER            :: num_rbin, tot_rbin, skip_rbin
    CHARACTER(LEN=500) :: gofr_dir
    CHARACTER(LEN=100) :: gofr_prefix, gofr_suffix
    CHARACTER(LEN=500) :: gofr_fname
    INTEGER            :: gofr_unit

    REAL(KIND=REAL32), ALLOCATABLE :: tempra(:,:), temprb(:,:)
    INTEGER           :: rbin

    !!!!------------------------------------------ S_ab(q) RELATED VARIABLES --------------------------------------!!!!
    REAL(KIND=REAL32), ALLOCATABLE :: iq(:), sq(:), sq_ab(:,:,:), iq_lorch(:), sq_lorch(:), sq_ab_lorch(:,:,:), q(:)
    REAL(KIND=REAL32), ALLOCATABLE :: sq_part(:,:,:), iq_part(:,:,:), sq_part_lorch(:,:,:), iq_part_lorch(:,:,:)
    REAL(KIND=REAL32), ALLOCATABLE :: sq_cat_cat(:), sq_an_an(:), sq_cat_an(:), iq_cat_cat(:), iq_an_an(:), iq_cat_an(:)
    REAL(KIND=REAL32), ALLOCATABLE :: sq_cat_cat_lorch(:), sq_an_an_lorch(:), sq_cat_an_lorch(:)
    REAL(KIND=REAL32), ALLOCATABLE :: iq_cat_cat_lorch(:), iq_an_an_lorch(:), iq_cat_an_lorch(:)
    REAL(KIND=REAL32), ALLOCATABLE :: iq_verify(:)

    REAL(KIND=REAL32)              :: qmin, qmin_box, qmin_gofrdata, qmax, delq=0.01
    REAL(KIND=REAL32)              :: twopi_by_l
    REAL(KIND=REAL32)              :: rtemp1, rtemp2
    INTEGER                        :: num_qbin

    INTEGER                        :: sabq_unit
    CHARACTER(LEN=500)             :: sabq_fname
    INTEGER                        :: sq_unit

    !!!!-------------------------------------- FORM FACTOR RELATED VARIABLES --------------------------------------!!!!
    REAL(KIND=REAL32), ALLOCATABLE :: ff_rawdata(:,:,:), ff_intdata(:,:), ff_avg(:), grpconc(:), ff_sq_avg(:)
    INTEGER                        :: num_ff = 62 ! NUMBER OF STANDARD F(Q) VALUES
    CHARACTER(LEN=2), ALLOCATABLE  :: ff_elname(:)
    CHARACTER(LEN=500)             :: ff_path = '/home/nikhil/projects/form_factors/'
    CHARACTER(LEN=500)             :: ff_fname
    INTEGER                        :: funit

    !!!!---------------------------------------- CLOCK RELATED VARIABLES ------------------------------------------!!!!
    INTEGER(KIND=8)   :: clock_start, clock_end, clock_rate, clock_max

    !!!!----------------------------------------- BOOK KEEPING VARIABLES ------------------------------------------!!!!
    INTEGER :: i, j, k, l, alpha, beta, itemp1, itemp2
    INTEGER :: ios


    CALL SYSTEM_CLOCK(count_rate=clock_rate)

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
    READ(inpunit,'(A)') gofr_dir
    WRITE(*,'(A,A)') "INPUT) GOFR DIRECTORY: ", TRIM(gofr_dir)
    
    READ(inpunit,*)
    READ(inpunit,'(A)') ndxfile_name
    WRITE(*,'(A,A)') "INPUT) NDXFILE: ", TRIM(ndxfile_name)

    READ(inpunit,*)
    READ(inpunit,*) box_size
    WRITE(*,'(A,F9.4)') "INPUT) BOX SIZE (ANG): ", box_size

    READ(inpunit,*)
    READ(inpunit,*) tot_rbin
    WRITE(*,'(A,I0)') "INPUT) NUMBER OF DIST CENTERS: ", tot_rbin 

    READ(inpunit,*)
    READ(inpunit,*) skip_rbin
    WRITE(*,'(A,I0)') "INPUT) NUMBER OF SKIP DIST CENTERS: ", skip_rbin 

    !! READ IN NUMBER OF SPECIES !!
    READ(inpunit,*)
    READ(inpunit,*) numspecies
    WRITE(*,'(A,I0)') "INPUT) NUMBER OF SPECIES: ", numspecies

    !! READ IN THE NAMES OF THE SPECIES !!
    ALLOCATE(species_name(numspecies))
    READ(inpunit,*)
    DO i = 1, numspecies
        READ(inpunit,*) species_name(i)
        WRITE(*,*) "INPUT) NAME OF SPECIES: ", i, " -- ", TRIM(ADJUSTL(species_name(i)))
    END DO

    !! READ IN THE NUMBER OF GROUPS IN EACH SPECIES
    ALLOCATE(num_grp_species(numspecies))
    READ(inpunit,*)
    READ(inpunit,*) num_grp_species(:)
    WRITE(*,*) "INPUT) NUMBER OF GROUPS IN EACH SPECIES: ", num_grp_species(:)

    !! READ IN THE GROUPS IN EACH SPECIES THEMSELVES !!
    ALLOCATE(grp_of_species(numspecies))
    DO i = 1, numspecies
        ALLOCATE(grp_of_species(i)%grpnum(num_grp_species(i)))
    END DO

    READ(inpunit,*)
    DO i = 1, numspecies
        READ(inpunit,*) grp_of_species(i)%grpnum(:) 
        WRITE(*,*) "INPUT) GROUPS IN SPECIES: ", i, " -- ", grp_of_species(i)%grpnum(:)
    END DO

    READ(inpunit,*)
    READ(inpunit,*) rbin_size
    WRITE(*,'(A,F9.4)') "INPUT) DISTANCE BIN SIZE (ANG): ", rbin_size

    READ(inpunit,*)
    READ(inpunit,*) delq
    WRITE(*,'(A,F6.3)') "INPUT) DELTA Q: ", delq
    
    READ(inpunit,*)
    READ(inpunit,*) qmin, qmax
    WRITE(*,'(A,F6.3,A,F6.3)') "INPUT) QMIN, QMAX: ", qmin," , ", qmax

    READ(inpunit,*)
    READ(inpunit,'(A)') ff_path
    WRITE(*,'(A,A)') "INPUT) FORM FACTOR PATH: ", TRIM(ff_path)
    
    READ(inpunit,*)
    READ(inpunit,*) num_ff
    WRITE(*,'(A,I0)') "INPUT) NUMBER OF FORM FACTOR DATA: ", num_ff

    READ(inpunit,*)
    READ(inpunit,'(A)') gofr_prefix
    WRITE(*,'(A,A)') "INPUT) GOFR PREFIX: ", TRIM(gofr_prefix)

    READ(inpunit,*)
    READ(inpunit,'(A)') gofr_suffix
    WRITE(*,'(A,A)') "INPUT) GOFR SUFFIX: ", TRIM(gofr_suffix)

    READ(inpunit,*)
    READ(inpunit,'(A)') suffix
    WRITE(*,'(A,A)') "INPUT) SUFFIX: ", TRIM(suffix)

    sabq_dir = TRIM(sabq_dir) // "-" // TRIM(ADJUSTL(suffix)) // "/"
    sabs_dir = TRIM(sabs_dir) // "-" // TRIM(ADJUSTL(suffix)) // "/"
    ff_dir = TRIM(ff_dir) // "-" // TRIM(ADJUSTL(suffix)) // "/"
    qpart_dir = TRIM(qpart_dir) // "-" // TRIM(ADJUSTL(suffix)) // "/"
    spart_dir = TRIM(spart_dir) // "-" // TRIM(ADJUSTL(suffix)) // "/"
    WRITE(*,'(A,A)') "OUTPUT) PARTIAL Sab(Q) INTO DIRECTORY: ", TRIM(sabq_dir)
    WRITE(*,'(A,A)') "OUTPUT) PARTIAL Sab(S) INTO DIRECTORY: ", TRIM(sabs_dir)
    WRITE(*,*)

    CLOSE(inpunit)
    !!!!********************************************************************************************************************!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!***************************************** READ IN INDICES OF ALPHA GROUPS ******************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL ndxopenfile(ndxfile_name)
    CALL ndxgetnumgrps(numgrps)
    ALLOCATE(alphagrps(numgrps),alphagrp_elname(numgrps),alphagrp_species(numgrps),alphagrp_mult(numgrps))
    ALLOCATE(alpha_moltype(numgrps))
    ALLOCATE(ff_elname(numgrps))
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

        !! ELEMENT NAME !!
        alphagrp_elname(i) = TRIM(ADJUSTL(alphagrps(i)%gname(alphagrp_namesplit(1)+1:alphagrp_namesplit(2)-1)))
        !! SPECIES NUMBER (SPECIES TO WHICH AN INDEX GROUP BELONGS TO !!
        READ(alphagrps(i)%gname(alphagrp_namesplit(2)+1:alphagrp_namesplit(3)-1),*) alphagrp_species(i)
        !! MULTIPLICITY (NUMBER OF ATOMS OF THE SAME TYPE IN A MOLECULE/RESIDUE) !!
        READ(alphagrps(i)%gname(alphagrp_namesplit(3)+1:),*) alphagrp_mult(i)
        alphagrps(i)%gname(alphagrp_namesplit(1):) = ' '

        !! CHECK THE TOTAL NUMBER OF ATOMS OF A GROUP IS DIVISIBLE BY MULTIPLICITY !!
        IF( MOD(alphagrps(i)%gnumatoms,alphagrp_mult(i)) .NE. 0) THEN
            WRITE(*,*) "ERROR) MULTIPLICITY AND NUMATOMS DON'T MATCH FOR THE INDEX GROUP : ", &
                alphagrp_mult(i), alphagrps(i)%gnumatoms, TRIM(alphagrps(i)%gname)
            STOP
        END IF

        !! ASSIGN A MOLECULE NUMBER TO AN ATOM WITHIN A SPECIES !!
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
    ff_elname(:) = alphagrp_elname(:)
    WRITE(*,'(A,I0)')   "----> TOTAL ATOMS : ", totgrpatoms
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!*********************************************** READING IN G(R) ****************************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) "***** READING IN G(R) ***** : "
    CALL SYSTEM_CLOCK(count=clock_start)

100 FORMAT(A,F9.4)
    !! SETUP R AND Q GRIDS !!    
    num_rbin = tot_rbin - skip_rbin
    volume = box_size**3
    rcut = box_size*0.5 
    rmax = REAL(2*num_rbin-1)*0.5*rbin_size
    twopi_by_l = 2.0*pi/box_size
    
    qmin_box = 4.0*pi/box_size 
    qmin_gofrdata = 2.0*pi/rmax
    IF(qmin .LT. delq) THEN
        WRITE(*,100) "QMIN WAS TOO LOW.. ", qmin
        qmin = qmin_box
        WRITE(*,100) "QMIN IS NOW SET TO QMIN_CUT.. ", qmin
    END IF
    num_qbin = INT((qmax-qmin)/delq) + 1


    WRITE(*,*)
    WRITE(*,100) "BOX SIZE                 : ", box_size
    WRITE(*,'(A,F12.2)') "VOLUME                   : ", volume
    WRITE(*,100) "MAX RCUT (BOX/2)         : ", rcut
    WRITE(*,*)
    
    WRITE(*,100) "RBIN SIZE                : ", rbin_size
    WRITE(*,'(A,I0)') "NUMBER OF DISTANCE BINS  : ", num_rbin
    WRITE(*,100) "RMAX (FROM GOFR DATA)    : ", rmax
    WRITE(*,*)

    WRITE(*,100) "QBIN SIZE                : ", delq
    WRITE(*,100) "QMIN                     : ", qmin
    WRITE(*,100) "QMAX                     : ", qmax
    WRITE(*,100) "QMIN_CUT (FROM BOX)      : ", qmin_box
    WRITE(*,100) "QMIN_CUT (FROM GOFR DATA): ", qmin_gofrdata
    WRITE(*,100) "2*PI/L                   : ", twopi_by_l
    WRITE(*,'(A,I0)') "NUMBER OF Q BINS         : ", num_qbin
    WRITE(*,*)

    ALLOCATE(dist(num_rbin))
    ALLOCATE(q(num_qbin))
    DO i = 1, num_qbin
        q(i) = qmin + REAL(i-1)*delq + 0.50*delq
    END DO

    ALLOCATE(gofr(num_rbin,numgrps,numgrps)) ; gofr(:,:,:) = 0.0

    !! READ ONE FILE TO GET DISTANCE POINTS AND CHECK FILE LENGTH !!
    gofr_fname = TRIM(ADJUSTL(gofr_dir))//'/'//&
                 TRIM(ADJUSTL(gofr_prefix))//'-'//&
                 TRIM(ADJUSTL(alphagrps(1)%gname))//'-'//&
                 TRIM(ADJUSTL(alphagrps(1)%gname))//'-'//&
                 TRIM(ADJUSTL(gofr_suffix))//'.dat'

    OPEN(NEWUNIT=gofr_unit, FILE=gofr_fname)
    DO i = 1, num_rbin
        READ(gofr_unit,*,IOSTAT=ios) dist(i), rtemp1, rtemp2
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR) UNABLE TO COMPLETE READING THE GOFR DATA FILE : ", TRIM(gofr_fname)
            STOP
        END IF
    END DO
    CLOSE(gofr_unit)

    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            gofr_fname = TRIM(ADJUSTL(gofr_dir))//'/'//&
                         TRIM(ADJUSTL(gofr_prefix))//'-'// &
                         TRIM(ADJUSTL(alphagrps(alpha)%gname))//'-'//&
                         TRIM(ADJUSTL(alphagrps(beta)%gname))//'-'//&
                         TRIM(ADJUSTL(gofr_suffix))//'.dat'

            OPEN(NEWUNIT=gofr_unit, FILE=gofr_fname)
            DO i = 1, num_rbin
                READ(gofr_unit,*, IOSTAT=ios) rtemp1, gofr(i,alpha,beta), rtemp2
                IF(ios .NE. 0) THEN
                    WRITE(*,'(A,A)') "ERROR) UNABLE TO COMPLETE READING THE GOFR DATA FILE : ", &
                               TRIM(gofr_fname)
                    STOP
                END IF
            END DO
            CLOSE(gofr_unit)
        END DO
    END DO
   
    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*) "***** WE HAVE G(R) NOW! ***** : "
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!*********************************************** CALCULATE S_AB(q) **************************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) "***** CALCULATING S_AB(q) ***** : "
    
    CALL SYSTEM_CLOCK(count=clock_start)

    ALLOCATE(sq_ab(num_qbin,numgrps,numgrps)) ; sq_ab(:,:,:) = 0.0
    ALLOCATE(sq_ab_lorch(num_qbin,numgrps,numgrps)) ; sq_ab_lorch(:,:,:) = 0.0
    ! check if making sq_ab into a reduction clause helps speed up -- i'm worried about clashes 

!$OMP PARALLEL DO SCHEDULE(STATIC, 1) DEFAULT(NONE) PRIVATE(i,alpha,beta,rtemp1,rtemp2) &
!$OMP& FIRSTPRIVATE(q,dist,num_rbin,twopi_by_l) &
!$OMP& SHARED(numgrps,num_qbin,alphagrps,sq_ab,sq_ab_lorch,gofr)
    DO alpha = 1, numgrps-1
        DO beta = alpha+1, numgrps
            DO i = 1, num_qbin
                ! INTEGRATE USING TRAPEZOID RULE
                 rtemp1 = q(i)

                rtemp2 = SUM(dist(2:num_rbin-1)*(gofr(2:num_rbin-1,alpha,beta)-1.0)*sin(rtemp1*dist(2:num_rbin-1)))
                rtemp2 = rtemp2 + 0.50*(dist(1)*(gofr(1,alpha,beta)-1.0)*sin(rtemp1*dist(1)))
                sq_ab(i,alpha,beta) = rtemp2 + &
                        & 0.50*(dist(num_rbin)*(gofr(num_rbin,alpha,beta)-1.0)*sin(rtemp1*dist(num_rbin)))
                 
                rtemp2 = SUM(dist(2:num_rbin-1)*(gofr(2:num_rbin-1,beta,alpha)-1.0)*sin(rtemp1*dist(2:num_rbin-1)))
                rtemp2 = rtemp2 + 0.50*(dist(1)*(gofr(1,beta,alpha)-1.0)*sin(rtemp1*dist(1)))
                sq_ab(i,beta,alpha) = rtemp2 + &
                        & 0.50*(dist(num_rbin)*(gofr(num_rbin,beta,alpha)-1.0)*sin(rtemp1*dist(num_rbin)))

                rtemp2 = SUM((gofr(2:num_rbin-1,alpha,beta)-1.0)*sin(rtemp1*dist(2:num_rbin-1)) &
                     *sin(twopi_by_l*dist(2:num_rbin-1)))
                rtemp2 = rtemp2 + 0.50*((gofr(1,alpha,beta)-1.0)*sin(rtemp1*dist(1))*sin(twopi_by_l*dist(1)))
                sq_ab_lorch(i,alpha,beta) = rtemp2 + 0.50*(dist(num_rbin)* &
                    (gofr(num_rbin,alpha,beta)-1.0)*sin(rtemp1*dist(num_rbin))*sin(twopi_by_l*dist(num_rbin)))

                rtemp2 = SUM((gofr(2:num_rbin-1,beta,alpha)-1.0)*sin(rtemp1*dist(2:num_rbin-1)) &
                     *sin(twopi_by_l*dist(2:num_rbin-1)))
                rtemp2 = rtemp2 + 0.50*((gofr(1,beta,alpha)-1.0)*sin(rtemp1*dist(1))*sin(twopi_by_l*dist(1)))
                sq_ab_lorch(i,beta,alpha) = rtemp2 + 0.50*(dist(num_rbin)* &
                    (gofr(num_rbin,beta,alpha)-1.0)*sin(rtemp1*dist(num_rbin))*sin(twopi_by_l*dist(num_rbin)))
            END DO
            ! DIVIDE BY Q(:)
            sq_ab(:,alpha,beta) = sq_ab(:,alpha,beta)/q(:)
            sq_ab(:,beta,alpha) = sq_ab(:,beta,alpha)/q(:)
            sq_ab_lorch(:,alpha,beta) = sq_ab_lorch(:,alpha,beta)/q(:)
            sq_ab_lorch(:,beta,alpha) = sq_ab_lorch(:,beta,alpha)/q(:)
        END DO
    END DO
!$OMP END PARALLEL DO

    DO alpha = 1, numgrps
        DO i = 1, num_qbin
            rtemp1 = q(i)
            rtemp2 = SUM(dist(2:num_rbin-1)*(gofr(2:num_rbin-1,alpha,alpha)-1.0)*sin(rtemp1*dist(2:num_rbin-1)))
            rtemp2 = rtemp2 + 0.50*(dist(1)*(gofr(1,alpha,alpha)-1.0)*sin(rtemp1*dist(1)))
            sq_ab(i,alpha,alpha) = rtemp2 + &
                    & 0.50*(dist(num_rbin)*(gofr(num_rbin,alpha,alpha)-1.0)*sin(rtemp1*dist(num_rbin)))
            
            rtemp2 = SUM((gofr(2:num_rbin-1,alpha,alpha)-1.0)*sin(rtemp1*dist(2:num_rbin-1)) &
                     *sin(twopi_by_l*dist(2:num_rbin-1)))
            rtemp2 = rtemp2 + 0.50*((gofr(1,alpha,alpha)-1.0)*sin(rtemp1*dist(1))*sin(twopi_by_l*dist(1)))
            sq_ab_lorch(i,alpha,alpha) = rtemp2 + 0.50*(dist(num_rbin)* &
                    (gofr(num_rbin,alpha,alpha)-1.0)*sin(rtemp1*dist(num_rbin))*sin(twopi_by_l*dist(num_rbin)))
        END DO
        ! DIVIDE BY Q(:)
        sq_ab(:,alpha,alpha) = sq_ab(:,alpha,alpha)/q(:)
        sq_ab_lorch(:,alpha,alpha) = sq_ab_lorch(:,alpha,alpha)/q(:)
    END DO

    ! MULTIPLY THE CONSTANT : 4*PI*DELTA-R*NUMBER-DENSITY
    rtemp1 = 4.0*pi*rbin_size*totgrpatoms/volume
    sq_ab(:,:,:) = sq_ab(:,:,:)*rtemp1

    ! 4*PI*delta-r*number-density/vol ==> divided by 2*PI/BOXSIZE (redundant, but increases clarity and security)
    rtemp1 = 4.0*pi*rbin_size*totgrpatoms/volume
    rtemp1 = rtemp1*box_size/(2.0*pi)
    sq_ab_lorch(:,:,:) = sq_ab_lorch(:,:,:)*rtemp1

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*) "***** WE HAVE S_ab NOW! ***** : "
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    !!---------------------->>> WRITE THE PARTIAL STRUCTURE FACTOR DATA TO FILES <<<----------------------!!
    WRITE(*,'(A)') " **** WRITING SAB(Q) TO FILES SYNCHRONOUSLY **** "
    CALL SYSTEM_CLOCK(count=clock_start)

    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            sabq_fname = TRIM(sabq_dir)//'sab-'//&
                         TRIM(ADJUSTL(alphagrps(alpha)%gname))//'-'//&
                         TRIM(ADJUSTL(alphagrps(beta)%gname))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO i = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(i), sq_ab(i,alpha,beta)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(sabq_dir)//'sab-lorch-'//&
                         TRIM(ADJUSTL(alphagrps(alpha)%gname))//'-'//&
                         TRIM(ADJUSTL(alphagrps(beta)%gname))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO i = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(i), sq_ab_lorch(i,alpha,beta)
            END DO
            CLOSE(sabq_unit)
        END DO
    END DO

    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            sabq_fname = TRIM(sabs_dir)//'sab-'//&
                         TRIM(ADJUSTL(alphagrps(alpha)%gname))//'-'//&
                         TRIM(ADJUSTL(alphagrps(beta)%gname))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO i = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(i)*one_by_2pi, sq_ab(i,alpha,beta)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(sabs_dir)//'sab-lorch-'//&
                         TRIM(ADJUSTL(alphagrps(alpha)%gname))//'-'//&
                         TRIM(ADJUSTL(alphagrps(beta)%gname))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO i = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(i)*one_by_2pi, sq_ab_lorch(i,alpha,beta)
            END DO
            CLOSE(sabq_unit)
        END DO
    END DO

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!******************************************* PREPARE THE FORM FACTORS ***********************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) "***** PREPARE THE FORM FACTORS ***** "
    CALL SYSTEM_CLOCK(count=clock_start)

120 FORMAT(A,2(F9.4,2X))
    !! READ IN THE RAW FORM FACTOR DATA !!
    WRITE(*,*) "--> READ IN THE RAW FORM FACTOR DATA "
    ALLOCATE(ff_rawdata(2,num_ff,numgrps))
    DO i = 1, numgrps
        ff_fname = TRIM(ff_path) // 'temp-ff-' // TRIM(ff_elname(i)) // '.dat'
        OPEN(NEWUNIT=funit, FILE=ff_fname, STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF(ios .NE. 0) THEN
            WRITE(*,*) "ERROR IN OPENING .. ", TRIM(ff_fname)
            STOP
        END IF
        ! DISCARD THE COMMENT IN FORM FACTOR FILE
        READ(funit,*)
        DO j = 1, num_ff
            READ(funit,*) ff_rawdata(1,j,i), ff_rawdata(2,j,i)
        END DO
        CLOSE(funit)
    END DO
    ! form factor data is given w.r.t to sin(theta)/lambda ; multiply 4pi to convert it to q
    ff_rawdata(1,:,:) = 4.0*pi*ff_rawdata(1,:,:)
    WRITE(*,120) "-> RAW FF DATA (Q) -- alpha 1 -- : ", MINVAL(ff_rawdata(1,:,1)), MAXVAL(ff_rawdata(1,:,1))
    WRITE(*,120) "-> RAW FF DATA (F) -- alpha 1 -- : ", MINVAL(ff_rawdata(2,:,1)), MAXVAL(ff_rawdata(2,:,1))

    !! INTERPOLATE THE FORM FACTOR DATA !!
    WRITE(*,*) "-> INTERPOLATING THE FORM FACTOR DATA... "
    ALLOCATE(ff_intdata(num_qbin,numgrps))
    DO alpha = 1, numgrps
        DO i = 1, num_qbin
            rtemp1 = q(i)
            IF((rtemp1 .LT. ff_rawdata(1,1,alpha)) .OR. (rtemp1 .GT. ff_rawdata(1,num_ff,alpha))) THEN
                WRITE(*,'(A,F6.2,A)') " THE Q VALUE : ", rtemp1, " IS OUT OF BOUNDS.. CANNOT EXTRAPOLATE YET"
                STOP
            END IF

            DO j = 2, num_ff
                IF(rtemp1 .LT. ff_rawdata(1,j,alpha)) THEN
                    ff_intdata(i,alpha) = ff_rawdata(2,j-1,alpha) + (rtemp1-ff_rawdata(1,j-1,alpha))* &
                    (ff_rawdata(2,j,alpha)-ff_rawdata(2,j-1,alpha))/(ff_rawdata(1,j,alpha)-ff_rawdata(1,j-1,alpha))
                    EXIT
                END IF
            END DO
        END DO
    END DO
    WRITE(*,120) "->               Q -- alpha * -- : ", MINVAL(q(:)), MAXVAL(q(:))
    WRITE(*,120) "-> INT FF DATA (F) -- alpha 1 -- : ", MINVAL(ff_intdata(:,1)), MAXVAL(ff_intdata(:,1))
    WRITE(*,120) "-> INT FF DATA (F) -- alpha L -- : ", MINVAL(ff_intdata(:,numgrps)), MAXVAL(ff_intdata(:,numgrps))

    ALLOCATE(grpconc(numgrps))
    DO alpha = 1, numgrps
        grpconc(alpha) = REAL(alphagrps(alpha)%gnumatoms)/REAL(totgrpatoms)
        WRITE(*,'(A,I0,A,A,A,I0,2X,I0,F9.4)') "-> GROUP ", alpha, " - " ,TRIM(alphagrps(alpha)%gname) ," Conc : ", &
                                          alphagrps(alpha)%gnumatoms, totgrpatoms, grpconc(alpha) 
    END DO

    !! CALCULATE THE AVERAGE FORM FACTOR !!
    ALLOCATE(ff_avg(num_qbin)) ; ff_avg(:) = 0.0
    ALLOCATE(ff_sq_avg(num_qbin)) ; ff_sq_avg(:) = 0.0
    DO i = 1, num_qbin
        DO alpha = 1, numgrps
            ff_avg(i) = ff_avg(i) + grpconc(alpha)*ff_intdata(i,alpha)
            ff_sq_avg(i) = ff_sq_avg(i) + grpconc(alpha)*ff_intdata(i,alpha)**2
        END DO
    END DO
    WRITE(*,120) "-> AVG FF ---------------------- : ", MINVAL(ff_avg(:)), MAXVAL(ff_avg(:))

    !! WRITE THE INTERPOLATED FORM FACTORS TO FILES !!
    DO i = 1, numgrps
        sabq_fname = TRIM(ADJUSTL(ff_dir))//"ff-"// &
                     TRIM(ADJUSTL(alphagrps(i)%gname))//"-"// &
                     TRIM(ADJUSTL(suffix)) // ".dat"
        OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
        DO j = 1, num_qbin
            WRITE(sq_unit,'(E14.6,2X,E14.6)') q(j), ff_intdata(j,i)
        END DO
        CLOSE(sq_unit)
    END DO

    !! WRITE THE AVERAGE FORM FACTORS TO FILES !!
    sabq_fname = "ff-avg-"//TRIM(ADJUSTL(suffix))//".dat"
    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO j = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(j), ff_avg(j)
    END DO
    CLOSE(sq_unit)

    !! WRITE THE AVERAGE OF SQUARED FORM FACTORS TO FILES !!
    sabq_fname = "ff-sq-avg-"//TRIM(ADJUSTL(suffix))//".dat"
    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO j = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(j), ff_sq_avg(j)
    END DO
    CLOSE(sq_unit)

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*) "***** FORM FACTORS ARE NOW PREPARED! ***** "
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!************************************************ CALCULATE S(q) ****************************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) "***** FINALLY CALCULATING S(q) ***** : "
    CALL SYSTEM_CLOCK(count=clock_start)

    ALLOCATE(sq(num_qbin)) ; sq(:) = 0.0
    ALLOCATE(sq_lorch(num_qbin)) ; sq_lorch(:) = 0.0
    DO alpha = 1, numgrps
        DO beta = 1, numgrps
            rtemp1 = grpconc(alpha)*grpconc(beta)
            DO i = 1, num_qbin
                ! S(Q) = {{SUM-A-B}} CA*CB*FA(Q)*FB(Q)*SAB(Q)/(FAVG(Q))**2
                sq(i) = sq(i) + rtemp1*ff_intdata(i,alpha)*ff_intdata(i,beta)*sq_ab(i,alpha,beta)/(ff_avg(i)**2)
                sq_lorch(i) = sq_lorch(i) + rtemp1*ff_intdata(i,alpha)* &
                              ff_intdata(i,beta)*sq_ab_lorch(i,alpha,beta)/(ff_avg(i)**2)
            END DO
        END DO
    END DO

    ALLOCATE(iq(num_qbin)) ; iq(:) = 0.0
    ALLOCATE(iq_lorch(num_qbin)) ; iq_lorch(:) = 0.0
    !! I(Q) = <F**2(Q)> + S(q)*<FF(Q)>**2
    iq(:) = ff_sq_avg(:) + sq(:)*ff_avg(:)**2
    iq_lorch(:) = ff_sq_avg(:) + sq_lorch(:)*ff_avg(:)**2

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)
    !!!!********************************************************************************************************************!!!!

    sabq_fname = "sofq-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i), sq(i)
    END DO
    CLOSE(sq_unit)

    sabq_fname = "sofs-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i)*one_by_2pi, sq(i)
    END DO
    CLOSE(sq_unit)

    sabq_fname = "sofq-lorch-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i), sq_lorch(i)
    END DO
    CLOSE(sq_unit)

    sabq_fname = "sofs-lorch-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i)*one_by_2pi, sq_lorch(i)
    END DO
    CLOSE(sq_unit)

    !!!!********************************************************************************************************************!!!!

    sabq_fname = "iofq-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i), iq(i)
    END DO
    CLOSE(sq_unit)

    sabq_fname = "iofs-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i)*one_by_2pi, iq(i)
    END DO
    CLOSE(sq_unit)

    sabq_fname = "iofq-lorch-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i), iq_lorch(i)
    END DO
    CLOSE(sq_unit)

    sabq_fname = "iofs-lorch-" // TRIM(ADJUSTL(suffix)) // ".dat"

    OPEN(NEWUNIT=sq_unit, FILE=sabq_fname)
    DO i = 1, num_qbin
        WRITE(sq_unit,'(E14.6,2X,E14.6)') q(i)*one_by_2pi, iq_lorch(i)
    END DO
    CLOSE(sq_unit)

    !!!!********************************************************************************************************************!!!!

!!<----------    HARD CODING STARTS HERE -------- RECTIFY IT LATER !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!************************************************ CALCULATE PARTIAL S(q) *********************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) "***** CALCULATING PARTIAL S(q) ***** : "
    CALL SYSTEM_CLOCK(count=clock_start)

    !  INTEGER                                :: numspecies
    ! INTEGER, ALLOCATABLE                   :: num_grp_species(:) !--> number of groups in each species
    ! TYPE(grptype), ALLOCATABLE             :: grp_of_species(:) !--> store the groups that belong a given species

    ALLOCATE(sq_part(num_qbin,numspecies,numspecies)) ; sq_part(:,:,:) = 0.0
    ALLOCATE(sq_part_lorch(num_qbin,numspecies,numspecies)) ; sq_part_lorch(:,:,:) = 0.0
    DO i = 1, numspecies
        DO j = 1, numspecies
            DO itemp1 = 1, num_grp_species(i)
                DO itemp2 = 1, num_grp_species(j)
                    ! itemp1 and itemp2 are just counters; alpha and beta are the group numbers
                    alpha = grp_of_species(i)%grpnum(itemp1)
                    beta = grp_of_species(j)%grpnum(itemp2)
                    rtemp1 = grpconc(alpha)*grpconc(beta)
                    DO k = 1, num_qbin
                    ! S(Q) = {{SUM-A-B}} CA*CB*FA(Q)*FB(Q)*SAB(Q)/(FAVG(Q))**2
                        sq_part(k,i,j) = sq_part(k,i,j) + rtemp1*ff_intdata(k,alpha)*ff_intdata(k,beta)* &
                                         sq_ab(k,alpha,beta)/(ff_avg(k)**2)
                        sq_part_lorch(k,i,j) = sq_part_lorch(k,i,j) + rtemp1*ff_intdata(k,alpha)* &
                              ff_intdata(k,beta)*sq_ab_lorch(k,alpha,beta)/(ff_avg(k)**2)
                    END DO
                END DO
            END DO
        END DO
    END DO

    ALLOCATE(iq_part(num_qbin,numspecies,numspecies)) ; iq_part(:,:,:) = 0.0
    ALLOCATE(iq_part_lorch(num_qbin,numspecies,numspecies)) ; iq_part_lorch(:,:,:) = 0.0
    ALLOCATE(iq_verify(num_qbin)) ; iq_verify(:) = 0.0
    DO i = 1, numspecies
        DO j = 1, numspecies
            !! I(Q) = <F**2(Q)> + S(q)*<FF(Q)>**2
            iq_part(:,i,j) = sq_part(:,i,j)*ff_avg(:)**2
            iq_part_lorch(:,i,j) = sq_part_lorch(:,i,j)*ff_avg(:)**2
            iq_verify(:) = iq_verify(:) + iq_part(:,i,j)
        END DO
    END DO
    iq_verify(:) = iq_verify(:) + ff_sq_avg(:)


    ALLOCATE(sq_cat_cat(num_qbin), sq_an_an(num_qbin), sq_cat_an(num_qbin))
    ALLOCATE(sq_cat_cat_lorch(num_qbin), sq_an_an_lorch(num_qbin), sq_cat_an_lorch(num_qbin))

    sq_cat_cat(:) = 0.0 ;  sq_an_an(:) = 0.0 ;  sq_cat_an(:) = 0.0
    sq_cat_cat_lorch(:) = 0.0 ;  sq_an_an_lorch(:) = 0.0 ;  sq_cat_an_lorch(:) = 0.0

    sq_cat_cat(:) = sq_part(:,1,1)
    sq_cat_cat_lorch(:) = sq_part_lorch(:,1,1)
    IF (numspecies == 2) THEN
        sq_an_an(:) = sq_part(:,2,2)
        sq_an_an_lorch(:) = sq_part_lorch(:,2,2)
    ELSE 
        sq_an_an(:) = sq_part(:,2,2) + sq_part(:,3,3)
        sq_an_an(:) = sq_an_an(:) + sq_part(:,2,3) + sq_part(:,3,2)

        sq_an_an_lorch(:) = sq_part_lorch(:,2,2) + sq_part_lorch(:,3,3)
        sq_an_an_lorch(:) = sq_an_an_lorch(:) + sq_part_lorch(:,2,3) + sq_part_lorch(:,3,2)
    END IF

    IF (numspecies == 2) THEN
        sq_cat_an(:) = sq_part(:,1,2) + sq_part(:,2,1)
        sq_cat_an_lorch(:) = sq_part_lorch(:,1,2) + sq_part_lorch(:,2,1)
    ELSE 
        sq_cat_an(:) = sq_part(:,1,2) + sq_part(:,2,1)
        sq_cat_an(:) = sq_cat_an(:) + sq_part(:,1,3) + sq_part(:,3,1)

        sq_cat_an_lorch(:) = sq_part_lorch(:,1,2) + sq_part_lorch(:,2,1)
        sq_cat_an_lorch(:) = sq_cat_an_lorch(:) + sq_part_lorch(:,1,3) + sq_part_lorch(:,3,1)
    END IF

    ALLOCATE(iq_cat_cat(num_qbin), iq_an_an(num_qbin), iq_cat_an(num_qbin))
    ALLOCATE(iq_cat_cat_lorch(num_qbin), iq_an_an_lorch(num_qbin), iq_cat_an_lorch(num_qbin))

    iq_cat_cat(:) = sq_cat_cat(:)*ff_avg(:)**2
    iq_an_an(:) = sq_an_an(:)*ff_avg(:)**2
    iq_cat_an(:) = sq_cat_an(:)*ff_avg(:)**2

    iq_cat_cat_lorch(:) = sq_cat_cat_lorch(:)*ff_avg(:)**2
    iq_an_an_lorch(:) = sq_an_an_lorch(:)*ff_avg(:)**2
    iq_cat_an_lorch(:) = sq_cat_an_lorch(:)*ff_avg(:)**2

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!************************************************ WRITE PARTIAL S(q) ************************************************!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!---------------------->>> WRITE THE PARTIAL STRUCTURE FACTOR DATA TO FILES <<<----------------------!!
    WRITE(*,'(A)') " **** WRITING SPART(Q) TO FILES SYNCHRONOUSLY **** "
    CALL SYSTEM_CLOCK(count=clock_start)

    DO i = 1, numspecies
        DO j = 1, numspecies
            sabq_fname = TRIM(qpart_dir)//'sqpart-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_part(k,i,j)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(qpart_dir)//'sqpart-lorch-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_part_lorch(k,i,j)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(qpart_dir)//'iqpart-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_part(k,i,j)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(qpart_dir)//'iqpart-lorch-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_part_lorch(k,i,j)
            END DO
            CLOSE(sabq_unit)
        END DO
    END DO

    DO i = 1, numspecies
        DO j = 1, numspecies
            sabq_fname = TRIM(spart_dir)//'sspart-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_part(k,i,j)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(spart_dir)//'sspart-lorch-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_part_lorch(k,i,j)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(spart_dir)//'ispart-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_part(k,i,j)
            END DO
            CLOSE(sabq_unit)

            sabq_fname = TRIM(spart_dir)//'ispart-lorch-'//&
                         TRIM(ADJUSTL(species_name(i)))//'-'//&
                         TRIM(ADJUSTL(species_name(j)))//'.dat'

            OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
            DO k = 1, num_qbin
                WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_part_lorch(k,i,j)
            END DO
            CLOSE(sabq_unit)
        END DO
    END DO

    ! PRINT OTHER USEFUL FILES TO HELP DEBUGING !
    sabq_fname = 'iofq-nosqff-'// TRIM(ADJUSTL(suffix)) //'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq(k) - ff_sq_avg(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = 'iofs-nosqff-'// TRIM(ADJUSTL(suffix)) //'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq(k) - ff_sq_avg(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = 'iofq-verify-'// TRIM(ADJUSTL(suffix)) //'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_verify(k)
    END DO
    CLOSE(sabq_unit)

    
    !!---------------------->>> WRITE THE ION-ION PARTIAL STRUCTURE FACTOR DATA TO FILES <<<----------------------!!
    !!---------------------->>> WRITE THE ION-ION PARTIAL STRUCTURE FACTOR DATA TO FILES <<<----------------------!!

        ! REAL(KIND=REAL32), ALLOCATABLE :: sq_cat_cat(:), sq_an_an(:), sq_cat_an(:), iq_cat_cat(:), iq_an_an(:), iq_cat_an(:)
    ! REAL(KIND=REAL32), ALLOCATABLE :: sq_cat_cat_lorch(:), sq_an_an_lorch(:), sq_cat_an_lorch(:)
    ! REAL(KIND=REAL32), ALLOCATABLE :: iq_cat_cat_lorch(:), iq_an_an_lorch(:), iq_cat_an_lorch(:)

    sabq_fname = TRIM(qpart_dir)//'sq-cat-cat-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_cat_cat(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'sq-an-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_an_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'sq-cat-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_cat_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'sq-cat-cat-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_cat_cat_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'sq-an-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_an_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'sq-cat-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), sq_cat_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    !!!! !!!! !!!! 

    sabq_fname = TRIM(qpart_dir)//'iq-cat-cat-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_cat_cat(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'iq-an-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_an_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'iq-cat-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_cat_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'iq-cat-cat-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_cat_cat_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'iq-an-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_an_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(qpart_dir)//'iq-cat-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k), iq_cat_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    !!!!!! !!!!! !!!!!

    sabq_fname = TRIM(spart_dir)//'ss-cat-cat-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_cat_cat(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'ss-an-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_an_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'ss-cat-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_cat_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'ss-cat-cat-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_cat_cat_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'ss-an-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_an_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'ss-cat-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, sq_cat_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    !!!! !!!! !!!!

    sabq_fname = TRIM(spart_dir)//'is-cat-cat-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_cat_cat(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'is-an-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_an_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'is-cat-an-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_cat_an(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'is-cat-cat-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_cat_cat_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'is-an-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_an_an_lorch(k)
    END DO
    CLOSE(sabq_unit)

    sabq_fname = TRIM(spart_dir)//'is-cat-an-lorch-'//TRIM(ADJUSTL(suffix))//'.dat'
    OPEN(NEWUNIT=sabq_unit, FILE=sabq_fname)
    DO k = 1, num_qbin
        WRITE(sabq_unit,'(E14.6,2X,E14.6)') q(k)*one_by_2pi, iq_cat_an_lorch(k)
    END DO
    CLOSE(sabq_unit)


    !!---------------------->>> WRITE THE ION-ION PARTIAL STRUCTURE FACTOR DATA TO FILES <<<----------------------!!
    !!---------------------->>> WRITE THE ION-ION PARTIAL STRUCTURE FACTOR DATA TO FILES <<<----------------------!!

    CALL SYSTEM_CLOCK(count=clock_end)
    WRITE(*,'(A,F9.2)') " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    WRITE(*,*)

    WRITE(*,*) "***** ALL THE BEST :) ***** "

END PROGRAM sofq
