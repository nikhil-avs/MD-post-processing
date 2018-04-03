! FIRST VERSION OF INTERMITTENT HBOND TIME CORRELATION FUNCTION -- BMIM-CL-BF4 system

! Last modified : 21 Feb 2018 by Nikhil 

! Varible description :
! Some information about the topology is hard-coded
! c_ring stores the centroid position of the ring
! vy_ring vector --> centroid to CR atom
! vz_ring vector --> is the plane normal
! 

MODULE read_dump
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    ! num_frames and num_lines are not used now
    INTEGER :: num_atoms, num_frames, num_lines, curr_num_atoms, curr_frame_no, curr_line_no, curr_timestep
    INTEGER, ALLOCATABLE, DIMENSION(:) :: id, mol, atype
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: r
    REAL(KIND=REAL64) :: box_size(3), box_dim(6)
    CHARACTER(LEN=100) :: dum1, dum2

    PUBLIC :: open_file, read_next, close_file
    PUBLIC :: get_data, get_boxinfo

    PRIVATE :: num_atoms, num_frames, num_lines, curr_num_atoms
    PRIVATE :: curr_line_no, curr_frame_no, curr_timestep
    PRIVATE :: id, mol, atype, r, dum1, dum2
    PRIVATE :: box_size, box_dim
    
    CONTAINS
    SUBROUTINE open_file(file_name, file_handle, file_type)
        CHARACTER(LEN=200), INTENT(IN) :: file_name
        INTEGER, INTENT(IN) :: file_handle, file_type
        INTEGER :: ios

        OPEN(UNIT=file_handle, FILE=file_name, STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios .NE. 0) THEN
            PRINT *, "ERROR IN OPENENING", TRIM(file_name)
        END IF
        PRINT '(A,A,I0)', TRIM(file_name), " OPENED SUCCESSFULLY WITH HANDLE : ", file_handle
 
        ! ALLOCATE AND INITIALIZE VARIABLES HERE
        num_atoms =0 ; num_frames = 0 ; num_lines = 0
        curr_num_atoms = 0 ; curr_frame_no = 0 ; curr_line_no = 0
        curr_timestep = 0 
        
        READ(file_handle,'(A)') dum1
        READ(file_handle,*) curr_timestep
        READ(file_handle,'(A)') dum1
        READ(file_handle,*) num_atoms
        READ(file_handle,'(A)') dum1
        READ(file_handle,*) box_dim(1), box_dim(2)
        READ(file_handle,*) box_dim(3), box_dim(4)
        READ(file_handle,*) box_dim(5), box_dim(6)
        READ(file_handle,'(A)') dum1
        REWIND(file_handle)
        
        ALLOCATE(id(num_atoms), mol(num_atoms), atype(num_atoms))
        ALLOCATE(r(3,num_atoms))

        curr_num_atoms = num_atoms
        box_size(1) = box_dim(2) - box_dim(1)
        box_size(2) = box_dim(4) - box_dim(3)
        box_size(3) = box_dim(6) - box_dim(5)
       
        ! PERFORM CHECKS HERE
    END SUBROUTINE open_file

    SUBROUTINE read_next(file_handle)
        INTEGER, INTENT(IN) :: file_handle
        INTEGER :: ios, i
        
        curr_frame_no = curr_frame_no + 1

        READ(file_handle,'(A)') dum1
        READ(file_handle,*) curr_timestep
        READ(file_handle,'(A)') dum1
        READ(file_handle,*) curr_num_atoms
        READ(file_handle,'(A)') dum1
        READ(file_handle,*) box_dim(1), box_dim(2)
        READ(file_handle,*) box_dim(3), box_dim(4)
        READ(file_handle,*) box_dim(5), box_dim(6)
        READ(file_handle,'(A)') dum1
        curr_line_no = curr_line_no + 9

        IF(curr_num_atoms .NE. num_atoms) THEN
            PRINT *, "WARNING : NUMBER OF ATOMS CHANGING WITH TIME!"
        END IF

        DO i = 1, curr_num_atoms
            READ(UNIT=file_handle,FMT=*,IOSTAT=ios) id(i), mol(i), atype(i), r(:,i)
            curr_line_no = curr_line_no + 1
        END DO
    END SUBROUTINE read_next

    SUBROUTINE get_data(my_id, my_mol, my_type, my_r, n)
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(INOUT) :: my_id(n), my_mol(n), my_type(n)
        REAL(KIND=REAL64), INTENT(INOUT) :: my_r(3,n)
        INTEGER :: i

        IF(n .NE. num_atoms) THEN
            PRINT *, "MISMATCH BETWEEN NUMBER OF ATOMS"
        END IF

        DO i = 1, num_atoms
            my_id(i) = i
            my_mol(id(i)) = mol(i)
            my_type(id(i)) = atype(i)
            my_r(:,id(i)) = r(:,i)
        END DO
    END SUBROUTINE get_data

    SUBROUTINE get_boxinfo(my_num_atoms, my_curr_frame_no, my_curr_line_no, my_curr_timestep, my_box_dim)
        INTEGER, INTENT(OUT) :: my_num_atoms, my_curr_frame_no, my_curr_line_no, my_curr_timestep
        REAL(KIND=REAL64) :: my_box_dim(6)

        my_num_atoms = curr_num_atoms
        my_curr_frame_no = curr_frame_no
        my_curr_line_no = curr_line_no
        my_curr_timestep = curr_timestep
        my_box_dim(1:6) = box_dim(1:6)
    END SUBROUTINE get_boxinfo

    SUBROUTINE close_file(file_handle)
        INTEGER, INTENT(IN) :: file_handle
        CLOSE(UNIT=file_handle)
        DEALLOCATE(r, id, mol, atype)
        PRINT *, "FILE CLOSED SUCCESSFULLY"
    END SUBROUTINE close_file
END MODULE read_dump

MODULE vector_ops
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    INTERFACE OPERATOR (.cross.)
        PROCEDURE cross_prod
    END INTERFACE OPERATOR (.cross.)
    CONTAINS
    FUNCTION cross_prod(r1, r2)
        REAL(KIND=REAL64), INTENT(IN) :: r1(3), r2(3)
        REAL(KIND=REAL64)             :: cross_prod(3)
        cross_prod(1) = r1(2)*r2(3) - r1(3)*r2(2)
        cross_prod(2) = r1(3)*r2(1) - r1(1)*r2(3)
        cross_prod(3) = r1(1)*r2(2) - r1(2)*r2(1)
    END FUNCTION cross_prod
    SUBROUTINE v_normalize(r1)
        REAL(KIND=REAL64), INTENT(INOUT) :: r1(3)
        r1 = r1/SQRT(DOT_PRODUCT(r1,r1))
    END SUBROUTINE v_normalize
END MODULE vector_ops

PROGRAM sdf
    USE read_dump
    USE vector_ops
    IMPLICIT NONE

    REAL(KIND=REAL64), PARAMETER :: pi = 4.D0*ATAN(1.D0)
    REAL(KIND=REAL64), PARAMETER :: au_conv = 0.1889726125D+1
    INTEGER, PARAMETER           :: file_hno = 666
    INTEGER, PARAMETER           :: dump_freq = 5   ! 5ps per frame

! DUMP FILE RELATED VARIABLE
    INTEGER                                      :: atom_no, frame_no, line_no, timestep
    INTEGER, ALLOCATABLE, DIMENSION(:)           :: id, mol, at_type
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: r
    REAL(KIND=REAL64)                            :: box_size(3), box_dim(6), volume
! IL RELATED VARIABLES
    INTEGER :: no_il1, no_il2, atoms_cat, atoms_an1, atoms_an2
    INTEGER :: id_ring(5), id_ha, id_hb(2), id_tail, id_an1, id_an2B, id_an2F(4)
    INTEGER :: no_ring, no_ha, no_hb, no_tail, no_an1, no_an2B, no_an2F
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: c_ring, c_an1, c_an2B, c_an2F
! TCF RELATED VARIABLES
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:)         :: sij_an1, sij_an2F
    INTEGER, ALLOCATABLE, DIMENSION(:)             :: denom_an1, denom_an2F
    INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: numer_an1, numer_an2F
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: dnumer_an1, dnumer_an2F
    REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:)   :: ci_an1, ci_an2F
    INTEGER                                        :: mycount=0
! H-BOND RELATED VARIABLES
    REAL(KIND=REAL64)                  :: rcut_an1, rcut_an2B, rcut_an2F
    REAL(KIND=REAL64)                  :: acut_an1, acut_an2B, acut_an2F
    INTEGER                            :: hb_state
    REAL(KIND=REAL64)                  :: rd(3), rh(3), ra(3)
! TIME RELATED VARIABLES
    INTEGER                            :: frame_zero, t0, tf, delt
    INTEGER, ALLOCATABLE, DIMENSION(:) :: run_time
    INTEGER                            :: tot_frames, skip_frames, proc_frames
! BOOK KEEPING VARIABLES
    
    INTEGER            :: i=0, j=0, k=0, l=0, m=0, ios=0, temp1=0, temp2=0
    REAL(KIND=REAL64)  :: dum1=0.D0, dum2=0.D0
    CHARACTER(LEN=200) :: inpfile_name, dumpfile_name

!********** READ FROM INPUT FILE **********
    CALL GET_COMMAND_ARGUMENT(1, inpfile_name)
    OPEN(UNIT=10, FILE=inpfile_name, STATUS='OLD', ACTION ='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING INPUT FILE", TRIM(inpfile_name)
    END IF
    PRINT '(A,A,I0)', TRIM(inpfile_name), " OPENED SUCCESSFULLY"

    READ(10,*)
    READ(10,'(A)') dumpfile_name
    READ(10,*)
    READ(10,*) rcut_an1
    READ(10,*)
    READ(10,*) rcut_an2B
    READ(10,*)
    READ(10,*) rcut_an2F
    READ(10,*)
    READ(10,*) acut_an1
    READ(10,*)
    READ(10,*) tot_frames
    READ(10,*)
    READ(10,*) skip_frames
    READ(10,*)
    READ(10,*) no_il1
    READ(10,*)
    READ(10,*) no_il2
    READ(10,*)
    READ(10,*) atoms_cat
    READ(10,*)
    READ(10,*) atoms_an1
    READ(10,*)
    READ(10,*) atoms_an2
    READ(10,*)
    READ(10,*) id_ring(:)
    READ(10,*)
    READ(10,*) id_ha
    READ(10,*)
    READ(10,*) id_hb(:)
    READ(10,*)
    READ(10,*) id_tail
    READ(10,*)
    READ(10,*) id_an1
    READ(10,*)
    READ(10,*) id_an2B
    READ(10,*)
    READ(10,*) id_an2F(:)

    acut_an2F = acut_an1
    acut_an2B = acut_an1

    proc_frames = tot_frames - skip_frames

    PRINT *, "ANGLE CRITERIA : ", acut_an1

    OPEN(UNIT=50,FILE='ci-an1.out', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING OUTPUT FILE."
    END IF
    PRINT *, "ci-an1 OUTPUT FILE OPENED SUCCESSFULLY"

    OPEN(UNIT=60,FILE='ci-an2F.out', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING OUTPUT FILE."
    END IF
    PRINT *, "ci-an2F OUTPUT FILE OPENED SUCCESSFULLY"

    CALL open_file(dumpfile_name, file_hno, 1)
    CALL get_boxinfo(atom_no, frame_no, line_no, timestep, box_dim)

!********** INITIALIZE MORE VARIABLES **********

    no_ring = no_il1+no_il2
    no_ha   = no_il1+no_il2
    no_hb   = 2*(no_il1+no_il2)
    no_tail = no_il1+no_il2
    no_an1  = no_il1
    no_an2B = no_il2
    no_an2F = no_il2*4  
        

!********** SKIP FRAMES **********
    PRINT *, "***** SKIPPING FRAMES ***** : ", skip_frames
    DO i = 1, skip_frames
        CALL read_next(file_hno)
        PRINT *, "skipping FRAME : ", i
    END DO
!********** END SKIP FRAMES **********

    volume = 0.D0 
    box_size(1) = box_dim(2) - box_dim(1)
    box_size(2) = box_dim(4) - box_dim(3)
    box_size(3) = box_dim(6) - box_dim(5)

    ALLOCATE(id(atom_no), mol(atom_no), at_type(atom_no), r(3,atom_no))
    ALLOCATE(c_ring(3,no_ring), c_an1(3,no_an1), c_an2B(3,no_an2B), c_an2F(3,no_an2F))

    ALLOCATE(sij_an1(no_an1, no_ha, proc_frames), sij_an2F(no_an2F, no_ha, proc_frames))
    sij_an1(:,:,:)  = 0
    sij_an2F(:,:,:) = 0
    ALLOCATE(run_time(proc_frames))

!********** PROCESS FRAMES **********
    PRINT *, "***** READING FRAMES ***** : ", proc_frames
    DO i = 1, proc_frames
        CALL read_next(file_hno)
        CALL get_boxinfo(atom_no, frame_no, line_no, timestep, box_dim)
        CALL get_data(id, mol, at_type, r, atom_no)
        IF(i == 1) frame_zero = frame_no
        PRINT *, "READING FRAME : ", frame_no
        box_size(1) = box_dim(2) - box_dim(1)
        box_size(2) = box_dim(4) - box_dim(3)
        box_size(3) = box_dim(6) - box_dim(5)
        volume = volume + (box_size(1)*box_size(2)*box_size(3))

!********** IDENTIFY CENTERS **********
        DO l = 1, no_il1
            c_ring(:,l) = 0.D0
            temp1 = (l-1)*(atoms_cat+atoms_an1)
            DO m = 1, 5
                c_ring(:,l) = c_ring(:,l) + r(:,temp1+id_ring(m))
            END DO
            c_ring(:,l) = c_ring(:,l)/5.D0
            c_an1(:,l) = r(:,temp1+atoms_cat+id_an1)
        END DO

        DO l= no_il1+1, no_il1+no_il2
            c_ring(:,l) = 0.D0
            temp1 = no_il1*(atoms_cat+atoms_an1) + (l-1-no_il1)*(atoms_cat+atoms_an2)
            DO m = 1, 5
                c_ring(:,l) = c_ring(:,l) + r(:,temp1+id_ring(m))
            END DO
            c_ring(:,l) = c_ring(:,l)/5.D0
            c_an2B(:,l-no_il1) = r(:,temp1+atoms_cat+id_an2B)
            DO m = 1, 4
                c_an2F(:,(l-no_il1-1)*4+m) = r(:,temp1+atoms_cat+id_an2F(m))
            END DO
        END DO

        DO j = 1, no_ring
            
            IF (j <= no_il1) THEN
                temp1 = (atoms_cat+atoms_an1)*(j-1) 
            ELSE
                temp1 = (atoms_cat+atoms_an1)*no_il1 + (atoms_cat+atoms_an2)*(j-1-no_il1)
            END IF
            
            rh(:) = r(:,temp1+id_ha)
            rd(:) = r(:,temp1+id_ring(2))

            DO k = 1, no_an1
                ra(:) = c_an1(:,k)
                CALL id_hbond(rd, rh, ra, rcut_an1, acut_an1, hb_state)
                IF (hb_state == 1) THEN
                    sij_an1(k,j,i) = 1
                END IF
            END DO

            DO k = 1, no_an2F
                ra(:) = c_an2F(:,k)
                CALL id_hbond(rd, rh, ra, rcut_an2F, acut_an2F, hb_state)
                IF (hb_state == 1) THEN
                    sij_an2F(k,j,i) = 1
                END IF
            END DO
        END DO

        run_time(i) = dump_freq*(frame_no-frame_zero)

    END DO
    DEALLOCATE(id, mol, at_type, r)
    CALL close_file(file_hno)

    ALLOCATE(denom_an1(proc_frames-1), denom_an2F(proc_frames-1))
    ALLOCATE(dnumer_an1(proc_frames-1,proc_frames-1), dnumer_an2F(proc_frames-1,proc_frames-1))
    
    denom_an1(:)     = 0
    denom_an2F(:)    = 0
    dnumer_an1(:,:)  = 0.D0
    dnumer_an2F(:,:) = 0.D0

    DO t0 = 1, proc_frames-1
        PRINT '(A,I6,A)', "***** ESTIMATING TCF AT T0 : ", run_time(t0), " ps"
        DO j = 1, no_ring
            DO k = 1, no_an1
                denom_an1(t0) = denom_an1(t0) + sij_an1(k,j,t0)**2
            END DO

            DO k = 1, no_an2F
                denom_an2F(t0) = denom_an2F(t0) + sij_an2F(k,j,t0)**2
            END DO
        END DO

        IF(t0 == 1) PRINT *, "DENOMINATOR : ", denom_an1(t0)

        DO tf = t0+1, proc_frames
            ! delt goes from 1 to proc_frames-1 when t0=1 and
            !                1 to 1             when t0=proc_frames-1
            delt = tf-t0
            DO j = 1, no_ring
                DO k = 1, no_an1
                    dnumer_an1(delt,t0) = dnumer_an1(delt,t0) + REAL(sij_an1(k,j,t0)*sij_an1(k,j,tf))
                END DO
                DO k = 1, no_an2F
                    dnumer_an2F(delt,t0) = dnumer_an2F(delt,t0) + REAL(sij_an2F(k,j,t0)*sij_an2F(k,j,tf))
                END DO
            END DO
        END DO
        IF(t0 == 1) PRINT *, "NUMERATOR AT DELT 5ps : ", dnumer_an1(1,t0)
        dnumer_an1(:,t0)  = REAL(dnumer_an1(:,t0),KIND=REAL64)/REAL(denom_an1(t0),KIND=REAL64)  
        dnumer_an2F(:,t0) = REAL(dnumer_an2F(:,t0),KIND=REAL64)/REAL(denom_an2F(t0),KIND=REAL64)  
    END DO
    DEALLOCATE(sij_an1, sij_an2F)
    DEALLOCATE(denom_an1, denom_an2F)

    ALLOCATE(ci_an1(proc_frames-1), ci_an2F(proc_frames-1))
    ci_an1(:)  = 0.D0
    ci_an2F(:) = 0.D0

    PRINT *, "***** AVERAGING TCF *****"

    DO delt = 1, proc_frames-1
        mycount = 0
        DO t0 = 1, proc_frames-delt
            mycount = mycount + 1
            ci_an1(delt)  = ci_an1(delt) + dnumer_an1(delt,t0)
            ci_an2F(delt) = ci_an2F(delt) + dnumer_an2F(delt,t0)
        END DO
        ci_an1(delt)  = ci_an1(delt)/REAL(mycount,KIND=REAL64)
        ci_an2F(delt) = ci_an2F(delt)/REAL(mycount,KIND=REAL64)
    END DO
    DEALLOCATE(dnumer_an1, dnumer_an2F)

    PRINT *, "***** PRINTING TCF TO FILE *****"

    12 FORMAT(I6, 2X, E15.8)
    WRITE(50,12) run_time(1), 1.D0
    WRITE(60,12) run_time(1), 1.D0
    DO delt = 1, proc_frames-1
        WRITE(50,12) run_time(delt+1), ci_an1(delt)
        WRITE(60,12) run_time(delt+1), ci_an2F(delt)
    END DO
    DEALLOCATE(ci_an1, ci_an2F, run_time)

    ! CALL close_file(file_hno)
    CLOSE(10) ; CLOSE(50) ; CLOSE(60)
    PRINT *, "END OF CODE"

    CONTAINS
    SUBROUTINE id_hbond (myrd, myrh, myra, myrcut, myacut, iout)
        REAL(KIND=REAL64), INTENT(IN) :: myrd(3), myrh(3), myra(3)
        REAL(KIND=REAL64), INTENT(IN) :: myrcut, myacut
        INTEGER, INTENT(INOUT)        :: iout
        REAL(KIND=REAL64)             :: v_ch(3), v_da(3), v_ha(3), myr
        REAL(KIND=REAL64)             :: mod1, mod2, dotp

        v_ha(:) = myra(:) - myrh(:)
        ! minimum image distance
        v_ha(:) = v_ha(:) - box_size(:)*NINT(v_ha(:)/box_size(:))
        myr = SQRT(DOT_PRODUCT(v_ha,v_ha))
        iout = 0
        IF(myr < myrcut) THEN
            v_ch(:) = myrh(:) - myrd(:)
            v_ch(:) = v_ch(:) - box_size(:)*NINT(v_ch(:)/box_size(:))
            v_da(:) = myra(:) - myrd(:)
            v_da(:) = v_da(:) - box_size(:)*NINT(v_da(:)/box_size(:))
            mod1 = SQRT(DOT_PRODUCT(v_ch,v_ch))
            mod2 = SQRT(DOT_PRODUCT(v_da,v_da))
            dotp = DOT_PRODUCT(v_ch,v_da)
            ! angle < angle_cut
            ! COS(angle) > COS(angle_cut)
            ! A.B > |A|.|B|.COS(angle_cut) 
            IF(dotp > mod1*mod2*COS(myacut*pi/180.D0)) THEN
                iout = 1
            END IF
        END IF
    END SUBROUTINE id_hbond
END PROGRAM sdf