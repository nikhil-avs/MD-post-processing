! ACCELERATED VERSION OF THE VISCOSITY GREEN-KUBO 9(6) component CODE 

! PROGRAM started on  : 05 Nov 19
! LAST updated        : 18 Jul 2020
! LAST updated        : 26 May 2022
! LAST updated        : 05 Jun 2022 -- corrected printing error (DELT was not multiplied before)

!!!!!!!!! GROMACS  
! GROMACS OUTPUT: TIMESTEP ;; PXX ;; PXY ;; PXZ ;; PYX ;; PYY ;; PYZ ;; PZX ;; PZY ;; PZZ -- in bar
!                              1      2      3      4      5      6      7      8      9
! DIAGONAL COMPS      : 1 , 5 , 9
! UNIQUE NON-DIAGONAL : 2 , 3 , 6
!                       |   |   | 
!                     : 4 , 7 , 8
! THIS SCALING CONSTANT converts to cP UNITS, ASSUMING THAT VOLUME = ANG^3 ;; delT = fs ;; PRESSURE = bar

! THIS SCALING CONSTANT is (vol * 10^-9)/(1.380649 * T)
! Multiplying SCALING CONSTANT just multiplies V/(KbT)
! The Boltzmann constant is defined to be exactly 1.380649 * 10^-23 J/K
!!!!!!!!! GROMACS 

!!!!!!!!! LAMMPS 

! LAMMPS OUTPUT: TIMESTEP ;; PXY ;; PXZ ;; PYZ ;; PXX ;; PYY ;; PZZ -- in atmospheres
! UNIQUE NON-DIAGONAL : 1 , 2 , 3
! DIAGONAL COMPS      : 4 , 5 , 6
! THIS SCALING CONSTANT converts to cP UNITS, ASSUMING THAT VOLUME (vol) = ANG^3 ;; delT = fs ;; PRESSURE = atm
! Multiplying SCALING CONSTANT just multiplies V/(KbT)

! THIS SCALING CONSTANT is (vol/T) 10^-9 (1.01325^2/1.380649)

! The Boltzmann constant is defined to be exactly 1.380649 * 10^-23 J/K

!!!!!!!!! LAMMPS 


! fort.17 prints the 6 independent components of the pressure tensor 
! fort.19 prints the 6 independent components of the symmetric traceless pressure tensor

PROGRAM viscosity_gk
    IMPLICIT NONE
    
    ! ALL THE CONSTANTS HERE
    CHARACTER(LEN=10), PARAMETER  :: last_edit = '26 may 22'
    CHARACTER(LEN=200), PARAMETER :: source_dir = '/home/nikhil/projects/11_codes/&
                                     &viscosity/gpu-version/gk-9comp/may2022/'
    CHARACTER(LEN=100), PARAMETER :: source_name = 'gk-gpu-water.F08'
    INTEGER, PARAMETER            :: dp = 8
    REAL(KIND=8), PARAMETER       :: kb = 1.380649
    REAL(KIND=8), PARAMETER       :: atm_to_pascal = 1.01325
    REAL(KIND=8), PARAMETER       :: bar_to_pascal = 1.00

    ! UNITS RELATED
    REAL(KIND=8)               :: beta, temperature, scaling_const
    REAL(KIND=8)               :: box_x, box_y, box_z, box_volume
    INTEGER                    :: pchoice

    ! INPUT FILE RELATED
    INTEGER                    :: skip_vals, proc_vals, delt=1
    INTEGER                    :: num_timesteps, num_proc_timesteps
    REAL(KIND=8)               :: start_md_time, end_md_time, proc_frac
    
    ! CORE DATA-STRUCTURES
    REAL(KIND=4), ALLOCATABLE  :: p_tens(:,:), ptcf(:,:)
    REAL(KIND=8), ALLOCATABLE  :: run_integral(:,:)
    REAL(KIND=8), ALLOCATABLE  :: visc_diag(:), visc_nondiag(:), visc_3comp(:), visc_9comp(:)
    REAL(KIND=8)               :: pvalues(6), ptrace
    
    ! TIMING RELATED
    INTEGER(KIND=8)            :: clock_start, clock_end, clock_rate, clock_max
    REAL(KIND=8)               :: start_time, end_time
    
    CHARACTER(LEN=200)         :: inpfile_name, pres_tens_file, suffix_name
    INTEGER                    :: i, ios=0, t, t0, tf, delta_t
    REAL(KIND=8)               :: rdum1, rdum2, rdum3
    CHARACTER(LEN=1)           :: pstatus(0:10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** READ FROM INPUT FILE **********
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    
    CALL GET_COMMAND_ARGUMENT(1, inpfile_name)
    OPEN(UNIT=10, FILE=inpfile_name, STATUS='OLD', ACTION ='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING INPUT FILE", TRIM(inpfile_name)
        STOP
    END IF
    PRINT *, TRIM(inpfile_name), " OPENED SUCCESSFULLY"

    READ(10,*)
    READ(10,'(A)') pres_tens_file
    READ(10,*) 
    READ(10,'(A)') suffix_name
    ! box volume should be in ang^3
    READ(10,*) 
    READ(10,*) box_volume
    READ(10,*) 
    READ(10,*) temperature
    ! DELT should be in fs -- it is the difference between two successive ptens dumps
    READ(10,*) 
    READ(10,*) delt
    READ(10,*) 
    READ(10,*) skip_vals
    READ(10,*) 
    READ(10,*) proc_vals
    READ(10,*) 
    READ(10,*) proc_frac

    num_timesteps = proc_vals-1

    ! THIS SCALING CONSTANT IS IN cP UNITS, ASSUMING THAT VOLUME = ANG^3 ;; delT = fs ;; PRESSURE = bar
    scaling_const = 1e-9_dp* (box_volume/temperature) * (atm_to_pascal**2 / kb)

    PRINT '(A,E12.5)', "BOX VOLUME       : ", box_volume
    PRINT '(A,E12.5)', "SCALING CONSTANT : ", scaling_const

    WRITE(*,*) " -------- INPUT FILE STARTS -------" 
    WRITE(*,'(A)') "1. PRESSURE TENSOR FILE: "
    WRITE(*,'(A)') TRIM(pres_tens_file)
    WRITE(*,'(A)') "2. SUFFIX NAME         : " 
    WRITE(*,'(A)') TRIM(suffix_name)
    WRITE(*,'(A)') "3. VOLUME              : "
    WRITE(*,*) box_volume
    WRITE(*,'(A)') "4. TEMPERATURE         : "
    WRITE(*,*) temperature
    WRITE(*,'(A)') "5. DELT (fs UNITS)     : "
    WRITE(*,*) delt
    WRITE(*,'(A)') "6. SKIP VALS           : "
    WRITE(*,*) skip_vals
    WRITE(*,'(A)') "7. PROCESS VALS        : "
    WRITE(*,*) proc_vals
    WRITE(*,*) " -------- INPUT FILE ENDS   -------"


    OPEN(UNIT=20, FILE=pres_tens_file, STATUS='OLD', ACTION ='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING PRESSURE TENSOR FILE", TRIM(pres_tens_file)
        STOP
    END IF
    PRINT *, TRIM(pres_tens_file), " OPENED SUCCESSFULLY"

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** skipping values **********
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
10  FORMAT(A1,A,11A,A,F6.2)
    pstatus(:) = ' '
    PRINT *, "SKIPPING VALUES... "
    DO i = 1, skip_vals
        READ(20,*)
        IF(MOD(i,10000) == 0) THEN
            pstatus(i*10/skip_vals) = '='
            WRITE(*,10,ADVANCE='NO') ACHAR(13), "---> STATUS : [", pstatus(:), "]", REAL(i)*100.0/REAL(skip_vals)
        END IF
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** Reading in the values **********
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pstatus(:) = ' '
    ALLOCATE(p_tens(6,0:num_timesteps))

    PRINT *, " "
    PRINT *, "READING VALUES... "
    READ(20,*,IOSTAT=ios) rdum1, pvalues(:)
    start_md_time = rdum1
    
    ptrace = (pvalues(4) + pvalues(5) + pvalues(6))/3.0_dp
    p_tens(1,0) = pvalues(4) - ptrace
    p_tens(2,0) = pvalues(5) - ptrace
    p_tens(3,0) = pvalues(6) - ptrace
    p_tens(4,0) = pvalues(1)
    p_tens(5,0) = pvalues(2)
    p_tens(6,0) = pvalues(3)


    CALL SYSTEM_CLOCK(count=clock_start)
    DO i = 1, num_timesteps

        READ(20,*,IOSTAT=ios) rdum1, pvalues(:)

        ptrace = (pvalues(4) + pvalues(5) + pvalues(6))/3.0_dp
        p_tens(1,i) = pvalues(4) - ptrace
        p_tens(2,i) = pvalues(5) - ptrace
        p_tens(3,i) = pvalues(6) - ptrace
        p_tens(4,i) = pvalues(1)
        p_tens(5,i) = pvalues(2)
        p_tens(6,i) = pvalues(3)

        ! To handle broken files
        IF (IS_IOSTAT_END(ios)) THEN
            PRINT *, "WARNING) EOF encountered. TERMINATING at time : ", rdum1, " proc_timestep : ", rdum1-start_md_time
            end_md_time = rdum1 ; num_timesteps = i-1
            EXIT
        END IF
        ! Printing in an elegant way
        IF(MOD(i,10000) == 0) THEN
            pstatus(i*10/proc_vals) = '='
            WRITE(*,10,ADVANCE='NO') ACHAR(13), "---> STATUS : [", pstatus(:), "]", REAL(i)*100.0/REAL(proc_vals)
        END IF
    END DO
    end_md_time = rdum1 
    CALL SYSTEM_CLOCK(count=clock_end)
    PRINT '(A,F9.2)', " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF(proc_frac > 1.0 .OR. proc_frac < 0.0 ) THEN
        proc_frac = 1.0_dp
    END IF
    num_proc_timesteps = INT(REAL(num_timesteps,KIND=8)*proc_frac)
    
    PRINT *, " "
    PRINT '(A,F12.4)', "START TIME   : ", start_md_time
    PRINT '(A,F12.4)', "END TIME     : ", end_md_time
    PRINT '(A,I0)',    "MAX DELTA Ts : ", num_timesteps
    PRINT '(A,I0)',    "NUM DELTA Ts : ", num_proc_timesteps

    PRINT '(A)', "CALCULATING THE TCF... "
    ALLOCATE(ptcf(6,0:num_proc_timesteps)) ; ptcf(:,:) = 0.0
    pstatus(:) = ' '
    CALL SYSTEM_CLOCK(count=clock_start)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** Calculating the TCF **********<-----------	MOST TIME CONSUMING PART ----------- OPTIMIZE -------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$acc parallel loop copyin(p_tens(1:6,0:num_timesteps)) copy(ptcf(1:6,0:num_proc_timesteps))
    DO delta_t = 0, num_proc_timesteps
        DO t0 = 0, num_timesteps-delta_t
            tf = t0 + delta_t
            DO i = 1, 6
                ptcf(i,delta_t) = ptcf(i,delta_t) + p_tens(i,t0)*p_tens(i,tf)
            END DO
        END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** Calculating the TCF **********<-----------	MOST TIME CONSUMING PART ----------- OPTIMIZE -------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CALL SYSTEM_CLOCK(count=clock_end)
    PRINT '(A,F9.2)', " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    PRINT *, " "
    DEALLOCATE(p_tens)

!$acc kernels
    DO delta_t = 0, num_proc_timesteps
        ptcf(:,delta_t) = ptcf(:,delta_t)/REAL(num_timesteps-delta_t+1,KIND=8)
    END DO
!$acc end kernels

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** printing to file **********    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PRINT '(A)', "PRINTING TO TCF FILES... "

    OPEN(UNIT=18, FILE="tot-corr-"//TRIM(suffix_name)//".xvg", IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING OUTPUT FILE", "tot-corr-"//TRIM(suffix_name)//".xvg"
        STOP
    END IF
    PRINT *, "tot-corr-"//TRIM(suffix_name)//".xvg", " OPENED SUCCESSFULLY"
    CALL WriteXVGHeader("Total Pressure Correlation Function","time (ps)", "tot-corr (atm^2)", 18)

    pstatus(:) = ' '
    rdum1 = 1.0/SUM(ptcf(1:3,0))
    rdum2 = 1.0/SUM(ptcf(4:6,0))
    rdum3 = REAL(delt,KIND=8)*0.001_dp

    CALL SYSTEM_CLOCK(count=clock_start)
    DO i = 0, num_proc_timesteps

        WRITE(18,'(F13.3,4X,E16.9)') REAL(i,KIND=8)*rdum3, SUM(ptcf(1:3,i)) + 2.d0*SUM(ptcf(4:6,i))

        IF(MOD(i+1,10000) == 0) THEN
            pstatus(i*10/num_proc_timesteps) = '='
            WRITE(*,10,ADVANCE='NO') ACHAR(13), "---> STATUS : [", pstatus(:), "]", REAL(i)*100.0/REAL(num_proc_timesteps)
        END IF
    END DO
    CALL SYSTEM_CLOCK(count=clock_end)
    PRINT '(A,F9.2)', " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    PRINT *, " "

    DO i = 30, 45
        CLOSE(i)
    END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** Integrating the TCF **********
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Trapezoidal rule 
    PRINT *, "EVALUATING RUNNING INTEGRAL... "
    ALLOCATE(run_integral(6,0:num_proc_timesteps)) 
    run_integral(:,:) = 0.d0
    
    ! This uses Trapezoidal rule for the numerical integration
    DO t = 1, num_proc_timesteps
        run_integral(:,t) = run_integral(:,t-1) + 0.50_dp*(ptcf(:,t)+ptcf(:,t-1))
    END DO
    run_integral(:,:) = run_integral(:,:)*scaling_const*REAL(delt,KIND=8)

    ALLOCATE(visc_diag(0:num_proc_timesteps),visc_nondiag(0:num_proc_timesteps)) 
    ALLOCATE(visc_3comp(0:num_proc_timesteps),visc_9comp(0:num_proc_timesteps)) 
    
    visc_diag(:) = (run_integral(1,:) + run_integral(2,:) + run_integral(3,:))/3.0_dp
    visc_nondiag(:) = (run_integral(4,:) + run_integral(5,:) + run_integral(6,:))/3.0_dp
    
    visc_3comp(:) = visc_nondiag(:)
    visc_9comp(:) = 0.10_dp * (3.0_dp*visc_diag(:) + 6.0_dp*visc_nondiag(:))

    PRINT *, " "

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !********** printing to file **********    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PRINT *, "PRINTING TO viscosity FILES... "

    ! 3-component viscosity = (nxx + nyy + nzz)/3
    OPEN(UNIT=38, FILE="tot-visc-3-comp-"//TRIM(suffix_name)//".xvg", IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING OUTPUT FILE", "tot-visc-3-comp-"//TRIM(suffix_name)//".xvg"
        STOP
    END IF
    PRINT *, "tot-visc-3-comp-"//TRIM(suffix_name)//".xvg", " OPENED SUCCESSFULLY"
    CALL WriteXVGHeader("VISCOSITY 3 components","time (fs)", "visc-3-comp (cP)", 38)

    ! 9-component viscosity = 1/10 * {(nxx + nyy + nzz) + 2*(nxy + nxz + nyz)}
    OPEN(UNIT=39, FILE="tot-visc-9-comp-"//TRIM(suffix_name)//".xvg", IOSTAT=ios)
    IF (ios .NE. 0) THEN
        PRINT *, "ERROR IN OPENENING OUTPUT FILE", "tot-visc-9-comp-"//TRIM(suffix_name)//".xvg"
        STOP
    END IF
    PRINT *, "tot-visc-9-comp-"//TRIM(suffix_name)//".xvg", " OPENED SUCCESSFULLY"
    CALL WriteXVGHeader("VISCOSITY 9 components","time (fs)", "visc-9-comp (cP)", 39)

    rdum3 = REAL(delt,KIND=8)*1.0_dp
    pstatus(:) = ' '
    CALL SYSTEM_CLOCK(count=clock_start)
    DO i = 0, num_proc_timesteps
        WRITE(38,'(F13.4,4X,E13.6)') REAL(i)*rdum3, visc_3comp(i)
        WRITE(39,'(F13.4,4X,E13.6)') REAL(i)*rdum3, visc_9comp(i)
       
        IF(MOD(i,10000) == 0) THEN
            pstatus(i*10/num_proc_timesteps) = '='
            WRITE(*,10,ADVANCE='NO') ACHAR(13), "---> STATUS : [", pstatus(:), "]", REAL(i)*100.0/REAL(num_proc_timesteps)
        END IF
    END DO

    CALL SYSTEM_CLOCK(count=clock_end)
    PRINT '(A,F9.2)', " -- TIME TAKEN (S) : ", REAL(clock_end-clock_start)/REAL(clock_rate)
    PRINT *, " "

    DEALLOCATE(run_integral)
    DEALLOCATE(visc_diag)
    DEALLOCATE(visc_nondiag)
    DEALLOCATE(visc_3comp)
    DEALLOCATE(visc_9comp)

    PRINT *, " "
    PRINT *, "$$$$$$$$$ END OF CODE $$$$$$$$"
    PRINT *, "$$$$$$$$$ ALL THE BEST :) $$$$$$$$"

    CONTAINS
    SUBROUTINE WriteXVGHeader(title, xaxis, yaxis, filehandle)
        CHARACTER(LEN=*), INTENT(IN) :: title
        CHARACTER(LEN=*), INTENT(IN) :: xaxis
        CHARACTER(LEN=*), INTENT(IN) :: yaxis
        INTEGER, INTENT(IN) :: filehandle

        CHARACTER(LEN=11) :: mydate
        CALL get_DDMonYY(mydate)

        WRITE(filehandle,'(A,A)') "# This file was created on ", mydate
        WRITE(filehandle,'(A,A)') "# Created by  : ", TRIM(source_name)
        WRITE(filehandle,'(A,A)') "# Source at   : ", TRIM(source_dir)
        WRITE(filehandle,'(A,A)') "# Last edited : ", TRIM(last_edit)

        WRITE(filehandle,'(A,A,A)') '@ title "', TRIM(title), '"' 
        WRITE(filehandle,'(A,A,A)') '@ xaxis label "', TRIM(xaxis), '"' 
        WRITE(filehandle,'(A,A,A)') '@ yaxis label "', TRIM(yaxis), '"' 
        ! sets the type of the set
        WRITE(filehandle,'(A)') '@ TYPE xy'
        ! sets the view port
        ! WRITE(filehandle,'(A)') '@ view 0.175, 0.1, 0.9, 0.9'
        ! switches on the legen display
        WRITE(filehandle,'(A)') '@ legend on'
        ! legend location type -- world or view -- world coords or viewport coords
        WRITE(filehandle,'(A)') '@ legend loctype view'
        ! legend location 
        WRITE(filehandle,'(A)') '@ legend 0.5,0.8'
        ! legen line length
        WRITE(filehandle,'(A)') '@ legend length 5'
    END SUBROUTINE WriteXVGHeader

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

END PROGRAM viscosity_gk
