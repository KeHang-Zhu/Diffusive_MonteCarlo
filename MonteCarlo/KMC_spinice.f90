
  !*******************************************************************
  !

  !  codes for Diffsive Monte Carlo.

  !  Author    : Kehang Zhu
  !  Date      : February 27th, 2021
  !  Last Edit : April     4th, 2022
  !*******************************************************************


  
  
  
  !====================== START VARIABLE MODULE =====================
  MODULE my_vrbls
    IMPLICIT NONE

    !-- common parameters and variables ------------------------------
    ! THIS IS ALMOST PROJECT-INDEPENDENT 
    real(8)         ,parameter  :: tm32   = 1.d0/(2.d0**32)
    real(8)         ,parameter  :: eps    = 1.d-14      ! very small number
    logical                     :: prtflg               ! flag for write2file
    integer         ,parameter  :: Mxint  = 2147483647  ! maximum integer
    integer         ,parameter  :: Mnint  =-2147483647  ! minimum integer

    integer(8)                  :: Nsamp                ! # samples in unit 'NBlck'
    integer(8)                  :: Totsamp
    integer(8)                  :: Ntoss                ! # samples to be thrown away
    integer(8)                  :: itoss
    
    integer(8)                  :: i1=0,i2=0,i3=0,kw=0
    integer(8)                  :: N_measure
    integer(8)                  :: N_print
    integer(8)                  :: N_write

    character(10)               :: cnf_file  = 'cnf.txt'
    character(10)               :: stat_file = 'stat.txt'
    character(10)               :: print_file= 'print.txt'
    character(10)               :: field1_file= 'field1.txt'
    character(20)               :: field2_file= 'field2.txt'
    character(20)               :: magnet1_file= 'magnet1.txt'
    character(20)               :: magnet2_file= 'magnet2.txt'
    character(20)               :: monopole_file= 'monopole.txt'
    character(20)               :: density1_file= 'density_mono1.txt'
    character(20)               :: density2_file= 'density_anti1.txt'
    character(20)               :: density3_file= 'density_mono2.txt'
    character(20)               :: density4_file= 'density_anti2.txt'

    !-- observe ----------------------------------------------------------
    ! THIS IS ALMOST PROJECT-INDEPENDENT
    integer         ,parameter  :: Max_block = 4000
    integer         ,parameter  :: N_block = 1000
    real            ,parameter  :: Max_cor=0.1
    logical                     :: flg_cor
    type obs                                            ! define obs type
        character(8)            :: nam                  ! the name of observables
        real(8)                 :: vnr                  ! current value (integer)
        real(8)                 :: val                  ! the accumulated value of observables
        real(8)                 :: cor                  ! correlation
        real(8)                 :: erb                  ! errobar
        real(8)                 :: blk(Max_block)       ! blocks for calculating the errorbars
        real(8)                 :: a                    !
        real(8)                 :: b                    ! obser=a*val+b
    end type

    integer         ,parameter  :: NObs_b=19             ! #basic observables
    integer         ,parameter  :: NObs_c=0             ! #composite observables
    integer         ,parameter  :: NObs=NObs_b+NObs_c
    type(Obs)                   :: Obser(NObs)
    integer(8)                  :: Npb                  ! # samples per block
    integer(8)                  :: Nbk                  ! the {Nbk}th block is the current one
    integer(8)                  :: Np                   ! the pointer

! collect the time sequence
    integer(8)                    :: sequce                   ! the pointer
    double precision ,allocatable :: timesequce(:,:)
    integer(8)                    :: Max_time

 
    !-----------------------------------------------------------------

    !-- parameters and variables -------------------------------------
    !! THIS IS PROJECT-DEPENDENT 
    real(8)         ,parameter  :: Pi=3.1415926536
    integer(1)      ,parameter  :: D = 3                    ! dimensionality
    integer(8)                  :: dim, Nsite
    integer(8)                  :: Lx, Ly, Lz, Vol          ! the system size
    integer         ,allocatable:: Ns(:)

    integer(1)      ,allocatable:: sp(:, :)
    
    real(8)                     :: J_coup
    real(8)                     :: Dipo
    real(8)                     :: h, magnet_x
    real(8)                     :: beta
    real(8)                     :: E0, v0

    double precision         ,allocatable :: p_flip(:), prob_list(:)
    double precision                      :: Accu_prob
    integer(8)                            :: flipspin
    double precision                      :: current_t
    integer(8)                            :: num_sensor


    !-----------------------------------------------------------------

    !-- Lattice and State --------------------------------------------
    !! THIS IS PROJECT-DEPENDENT 


    integer(1)                  :: dir                 ! # neighboring vertices
    integer(1)      ,parameter  :: sublatnum = 4
    integer(8)      ,allocatable:: nnsite1(:,:)          ! nearest neighbor site
    integer(8)      ,allocatable:: ms(:,:)          ! nearest neighbor site
    integer(8)      ,allocatable:: vector(:,:)
    integer(8)      ,allocatable:: latt(:)
    integer(8)      ,allocatable:: trahe(:)
    integer         ,allocatable:: backdir1(:, :)           ! the backward direction
    integer         ,allocatable:: back(:)
    integer         ,allocatable:: sublatt(:)
    integer         ,allocatable:: lattno(:, :)
    integer         ,allocatable:: neigh_Num(:)
    double precision,allocatable:: locs(:, :)
    integer         ,allocatable:: Pat(:)
    double precision,allocatable:: d_spin(:, :)            !distance matrix of the spins
    double precision,allocatable:: d_sensor(:, :)       !distance matrix of spin and sensor
    double precision,allocatable:: sensor_loc(:, :)
    INTEGER, ALLOCATABLE :: inv_sp(:), inv_trahe(:, :)


    double precision            :: dx, dy, dz
    integer(8)      ,allocatable:: break(:), break_inv(:)
    integer(8)                  :: mono1, mono2
    double precision,allocatable:: cent(:, :), cent_inv(:, :)
    double precision,allocatable:: density(:, :)

    logical                     :: pyrochl     ! whether pyro or not
    logical                     :: thermal     ! whether collect data
    
    real(8)                     :: p, dE
    double precision            :: enS_tot, enC_tot
    double precision            :: Mz
    double precision,allocatable:: enS(:), enC(:)
    double precision,allocatable:: Bz_sensor(:, :)


    real(8)                     :: wormstep
    integer(8)                  :: step, accept
    logical                     :: accept0, apply



    !-----------------------------------------------------------------


    !-- Random-number generator---------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    integer, parameter           :: mult=32781
    integer, parameter           :: mod2=2796203, mul2=125
    integer, parameter           :: len1=9689,    ifd1=471
    integer, parameter           :: len2=127,     ifd2=30
    integer, dimension(1:len1)   :: inxt1
    integer, dimension(1:len2)   :: inxt2
    integer, dimension(1:len1)   :: ir1
    integer, dimension(1:len2)   :: ir2
    integer                      :: ipnt1, ipnf1
    integer                      :: ipnt2, ipnf2
    integer, parameter           :: mxrn = 10000
    integer, dimension(1:mxrn)   :: irn(mxrn)

    integer                      :: Seed                   ! random-number seed
    integer                      :: nrannr                 ! random-number counter
    !-----------------------------------------------------------------

    !-- time-checking variables --------------------------------------
    ! THIS IS PROJECT-INDEPENDENT 
    character( 8)           :: date
    character(10)           :: time
    character(5 )           :: zone
    integer, dimension(8)   :: tval
    double precision        :: t_prev, t_curr, t_elap
    integer                 :: h_prev, h_curr
    double precision        :: t_init, t_simu, t_meas, t_toss
    !-----------------------------------------------------------------
  END MODULE my_vrbls
  !======================== END VARIABLE MODULE =====================
  

  !======================== START MAIN PROGRAM ======================
  
  PROGRAM Ising
    use my_vrbls
    implicit none
    integer :: isamp, isamp1,r_cnf_stat,i_each,N_each,iblck, time_i

    read *, r_cnf_stat
    read *, pyrochl
    read *, Lx
    read *, Ly
    read *, Lz
    read *, J_coup
    read *, Dipo
    read *, h
    read *, E0
    read *, v0
    read *, beta
    read *, Ntoss
    read *, Nsamp
    read *, N_print
    read *, N_write
    read *, N_each
    read *, N_measure
    read *, Seed
    read *, Max_time

    print *, 'r_cnf_stat, pyrochl , Lx, Ly, Lz, J, Dipo, h, E0, v0, beta, Ntoss, Nsamp, N_print, N_write, N_each, N_measure, Seed'
!    read  *,  r_cnf_stat, Lx, Ntoss, Nsamp, N_print, N_write, N_each, N_measure, t, h, beta, Seed
    print *,  r_cnf_stat, pyrochl , Lx, Ly, Lz, J_coup, Dipo, h, E0, v0, beta, Ntoss, Nsamp, N_print, &
            N_write, N_each, N_measure, Seed

    !--- Initialization ----------------------------------------------    
    call set_time_elapse

    call initialize
    if(r_cnf_stat==0) then
        print*,"create a new simulation:"
        call init_stat
        call system("rm "//trim(print_file))
        call system("rm "//trim(field1_file))
        call system("rm "//trim(field2_file))
        call system("rm "//trim(magnet1_file))
        call system("rm "//trim(magnet2_file))
        call system("rm "//trim(monopole_file))
        call system("rm "//trim(density1_file))
        call system("rm "//trim(density2_file))
        call system("rm "//trim(density3_file))
        call system("rm "//trim(density4_file))
    else if(r_cnf_stat == 1) then
        print*,"start a simulation on an old configuration:"
        call read_cnf()
        call init_stat
        call system("rm "//trim(print_file))
    else if(r_cnf_stat == 2) then
        print*,"resume an old simulation:"
        call read_cnf()
        call init_stat
        call read_stat()
        call printing(6)
    else
        stop "r_cnf_stat should be 0, 1 or 2"
    end if

    call time_elapse
    t_init = t_elap
    write(6,50) t_init
    50 format(/'        set up time:',f16.7,2x,'s')



do time_i = 1, 100
    !--- Thermialization ---------------------------------------------
    ! to make the sample fully thermalized and reach equilibrium
!do isamp1 = 1, Nsamp
    call set_Obs
    if(r_cnf_stat /= 2) then
        prtflg = .False.
        t_toss=0.d0
        beta = beta ;
        print*,'Initial energy', enC_tot, enS_tot
        thermal = .false.
        apply = .false.
print*,'num',num_sensor
print*, 'Bz', Bz_sensor(: , 3)
        ! thermailize without the magnetic field
        magnet_x = 0
        do isamp = 1, Ntoss
            do i_each = 1, N_each
                call markov
                call measure
            enddo
            i2 = i2 + 1
            if(i2 >= N_print) then
                beta = beta ;
                call time_elapse;   t_toss = t_toss + t_elap
                itoss = itoss + N_print

                call printing(6)
                open(1, file = trim(print_file), access = "append")
                call printing(1)
                close(1)
                call check
                i2=0
            end if
            i3 = i3 + 1
            if(i3 >= N_write) then
                print*, "write cnf. and stat."
                call write_cnf()
                call write_stat()
                print*, "write done"
                i3=0
            end if
        enddo

        call time_elapse
        t_toss = t_toss + t_elap
        write(6,51) t_toss
        51 format(/'thermalization time:',f16.7,2x,'s')
        t_simu   = 0.d0;   t_meas = 0.d0
        wormstep = 0; step = 0
        !Call check()
        print*, 'Totol energy', enC_tot, enS_tot
    else
        print*,"no need to thermalizate."
    end if

    !--- Simulation --------------------------------------------------

        prtflg=.True.
        apply = .true.

        print*, 'start simulation---------------------------------------------'
        print*, 'Apply magnetic field along [001] direction','B_x = ', h
        magnet_x = h
        enS_tot = get_enS_tot()
        enS = get_enS()
        thermal = .true.
        print*, 'Totol energy', enC_tot, enS_tot
        print*, 'Bz', Bz_sensor(:, 3)
        do i_each = 1, Max_time ![0,N_each] is the specific relaxtion path
           ! if (i_each > Max_time * Ly* Lz)then
           !     thermal = .false.
           ! endif
            call markov
            call measure
        end do
        print*, 'Energy after relaxation', enC_tot, enS_tot
        print*, 'initialization under magnetic field done.'

    ! record one specific relaxation path
       call writejjfile1()
       sequce = 1
       timesequce = 0
        
        apply = .false.
        call check()
        print*, '---------------------------------------------'
        print*, 'withdraw the magnetic field'
       ! thermal = .true.
        current_t = 0
        magnet_x = 0
        enS_tot = get_ens_tot()
        enS = get_ens()
        print*, 'Totol energy', enC_tot, enS_tot
        do i_each = 1, Max_time !N_each is the specific relaxtion path
            call markov
            call measure
        end do
        print*, 'relaxation done.'
        print*, 'Energy after relaxation', enC_tot, enS_tot


        ! record one specific relaxation path
        call writejjfile2()

        print*, 'write relaxation path done.'
seed = seed + 1
enddo
        stop
        
       ! i2 = i2 + 1
       ! if(i2 >= N_print) then
            !call cal_ComObs()
            !flg_cor = cal_cor_dev()
        !    call time_elapse;    t_simu = t_simu+t_elap
            
        !    call printing(6)
        !    open(1,file = trim(print_file), access = "append")
        !    call printing(1)
       !     close(1)
        !    call check
         !   i2 = 0
       ! end if
        !i3=i3+1
        !if(i3>=N_write) then
           ! print*,"write cnf. and stat."
          !  call write_cnf()
          !  call write_stat()
          !  print*,"write done"
        !   i3=0
        !end if

!! collet multiple data
   ! Seed = Seed + 1 ! change the ramdon number seed
    !end do
    stop

CONTAINS

    !--- PROJECT-DEPENDENT ---------------------------------------------
  !============== START INITIALIZATION ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE initialize
    implicit none

    !-- order should be changed --------------------------------------

    call allo_aux
    call def_latt
    call set_RNG
    call init_cnf
    call set_Obs
    !call def_prob
    call tst_and_prt

    return
  END SUBROUTINE initialize
 
  !==============Test and Print ======================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE tst_and_prt
    implicit none

    !-- Test and Print -----------------------------------------------

    write(6,*)
    write(6,*) "Spin ice relaxation model"
    write(6,40) Lx, Ly, Lz
    40 format('on ',i4,1x,"*",i4,1x,"*",i4,1x,'pyrochlore lattice')

    write(6,41) J_coup, Dipo, h, beta
    41 format(',where J_coup=',f12.8,' Dipo=',f12.8,' h=',f12.8,' beta=',f12.8)

    write(6,42) E0, v0
    42 format(',where E0=',f12.8,' v0=',f12.8)

    write(6,43) Nsamp, N_each
    43 format(' Will simulate      ',i10,'*',i10,'steps ')

    write(6,44) Ntoss, N_each
    44 format(' Throw away         ',i10,'*',i10,'steps')

    return
  END SUBROUTINE tst_and_prt


  !==============initialize configuration ======================================
  !! THIS IS PROJECT_DEPENDENT
  SUBROUTINE init_cnf
    implicit none
    integer(8) :: i,Vc
    !allocate(sp(Vol, 3))
    do i = 1, Vol
        if (sublatt(i) == 1) then
            sp(i,1) = 1; sp(i,2) = 1; sp(i,3) = 1;
        elseif (sublatt(i) == 2) then
            sp(i,1) = -1; sp(i,2) = -1; sp(i,3) = 1;
        elseif (sublatt(i) == 3) then
            sp(i,1) = 1; sp(i,2) = -1; sp(i,3) = -1;
        elseif (sublatt(i) == 4) then
            sp(i,1) = -1; sp(i,2) = 1; sp(i,3) = -1;
        endif
        !create the 2in-2out configuration
        if (sublatt(i) == 2) then
            sp(i, :) = - sp(i, :)
        elseif (sublatt(i) == 4) then
            sp(i, :) = - sp(i, :)
        endif
    enddo
    !define the direction of the spin moment
    Pat(1) = 1
    Pat(2) = -1
    Pat(3) = 1
    Pat(4) = -1


    !intial the time
    current_t = 0
    accept = 0

    return
  END SUBROUTINE


!===============allocate auxiliary arrays ======================================
 !! THIS IS PROJECT_DEPENDENT
 SUBROUTINE allo_aux
   implicit none

       dim = D
       allocate(Ns(dim))
       dir = 2 * dim
        Nsite=1

        Ns(1) = Lx; Ns(2) = Ly; Ns(3) = Lz
        Nsite = Lx * Ly * Lz

       if (pyrochl) then
           Nsite = Nsite * sublatnum
       endif
       Vol = Nsite
             print*, 'directions', dir
             print*, 'number of sites', Nsite

       num_sensor = 5  ! the sensors are placed by Ly sites
       allocate(vector(Nsite,dim))
       allocate(backdir1(Nsite, 7), nnsite1(Nsite,7))
       allocate(back(dim))
       allocate(latt(Nsite),trahe(Nsite), neigh_Num(Vol))
       allocate(sublatt(Nsite),lattno(Nsite,3))
       allocate(locs(Nsite, 3))
       allocate(sensor_loc(num_sensor, 3))
       allocate(d_spin(Vol, Vol))
       allocate(d_sensor(Vol, num_sensor))
       allocate(Pat(sublatnum))
       allocate(timesequce(NObs, Max_time))
       allocate(break(Vol/sublatnum), break_inv(Vol/sublatnum))
       allocate(cent(Vol/sublatnum, 3), cent_inv(Vol/sublatnum, 3))
       allocate(density(Lx, 2))
       ALLOCATE(inv_sp(Nsite), inv_trahe(4, Lx*Ly*Lz))
       allocate(sp(Vol, 3))
       allocate(enS(Vol), enC(Vol))
       allocate(ms(dir, Vol))

       allocate(p_flip(Vol))
       allocate(prob_list(Vol))
       
       allocate(Bz_sensor(num_sensor, 3))

        dx = Lx/2;
        dy = Ly/2;
        dz = Lz/2;
       !Max_block = N_each
        timesequce = 0
       print*, 'parameters done'

       
   return
 END SUBROUTINE



  
  !==============define flipping probability =========================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE def_prob()
    implicit none
    integer :: i,j
        do i = 1, Vol
            dE = - (2 * enS(i) + 2 * enC(i)) / 2.d0 + E0  !Arrhenius law
            p_flip(i) = dsqrt(v0) * dexp( - beta * dE )
            if (p_flip(i) > 1e10) then
                print*, p_flip(i), i, dE,  2 * enS(i) , 2 * enC(i)
                stop 'too large'
            endif
        end do
    return
  END SUBROUTINE def_prob


!! calculate the accumulation of probability in KMC
  SUBROUTINE accumulate_prob()
    implicit none
    integer :: i,j

        Accu_prob = 0;
        prob_list = 0
        ! call the def_prob subroutine for N times
        do i = 1, Vol
            Accu_prob = Accu_prob + p_flip(i)
            prob_list(i) = Accu_prob
        enddo

    return
  END SUBROUTINE accumulate_prob

!!using Binary search for the probavility in KMC
  SUBROUTINE search_prob(val)
    implicit none
    integer :: range, start, finish, mid
    integer :: i, n
    double precision, intent(in) :: val
    double precision :: x(Vol)
        
        x = prob_list
        start =  1
        finish = Vol
        range = finish - start
        mid = (start + finish) /2

        do while( (x(mid + 1) - val) * (x(mid) - val) > 0 .and. range >  0)
          if (val > x(mid)) then
            start = mid
          else if (val < x(mid)) then
            finish = mid
          else if (val > x(mid) .and. val < x(mid+1)) then
            exit
          end if
          range = finish - start
          mid = (start + finish)/2
        end do

        if( (x(mid + 1) - val) * (x(mid) - val) > 0 .and. mid /= 1) then
          print *, val, x(mid + 1), x(mid), start, finish, mid
            do i = 1, Vol
                print*, prob_list(i)
            enddo
          stop 'NOT FOUND'
        elseif (val < x(mid) )then
            flipspin = 1
        else
           flipspin = mid + 1
        end if

    return
  END SUBROUTINE search_prob


!! record the evolution time in KMC
  SUBROUTINE record_time()
    implicit none
    double precision :: delta_t

        delta_t = 1.d0 / Accu_prob /dsqrt(v0)  *log ( 1.d0/rn() )
        current_t = current_t + delta_t

    return
  END SUBROUTINE record_time
 
  !==============definite Lattice ====================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE def_latt
    implicit none
    integer(8)  :: i,j
    integer :: ir,ix,iy,iz,kf,kb,kcf,kcb,irr
    integer :: j1,n12,n1,nn(dim),ii,Ny,nn1, temporary_data, Vc,nsp, counter(Vol/sublatnum)
    integer :: neigh_spin, back_spin, site0, mm, n0
    double precision :: site1(Nsite,10), deltay, deltaz, deltax
    !if((Lx/2)*2/=Lx) then
    !    write(6,*) 'L be even?';   stop
   ! end if


       IF (pyrochl) then  ! use pyrochlore lattice
            open(1,file='backdir1.dat');  ! read the data of site
            DO j=1, Nsite;
                read(1,*) backdir1(j,1:7)
            END DO
            close(1)

            open(7,file='nnsite1.dat');  ! read the data of nniste
            DO j=1, Nsite
                read(7,*) nnsite1(j,1:7);
            ENDDO
            close(7)

            open(1,file='site.dat');  ! read the data of site
                READ(1, *) temporary_data

            DO j = 1, Nsite;
                read(1,*) site1(j,1:10);
               ! sublatt(j,1)=site1(j,1) ;sublatt(j,2)=site1(j,2) ;sublatt(j,3)=site1(j,3)!site num
    !site num
                latt(j) = site1(j,1) ;
    !belong to trahetra
                trahe(j) = site1(j,2)
    !sublattice num
                sublatt(j) = site1(j,3)
    !in plane lattice site number
                !lattno(1,j)=site1(j,5);lattno(2,j)=site1(j,6) ;lattno(3,j)=site1(j,7)
                lattno(j,1) = site1(j,5);lattno(j,2)=site1(j,6) ;lattno(j,3)=site1(j,7)
    !absolute location
                locs(j,1) = site1(j,8);locs(j,2)=site1(j,9);locs(j,3)=site1(j,10)

             !!reshape the lattice into cubic lattice
                if ((locs(j, 1) - Lx) > 1e-3 .or. locs(j, 1) == Lx) then
                    locs(j, 1) = locs(j, 1) - int( locs(j, 1) / Lx) * Lx
                endif
                if ((locs(j, 2) - Ly) > 1e-3 .or. locs(j, 2) == Ly) then
                    locs(j, 2) = locs(j, 2) - int( locs(j, 2) / Ly) * Ly
                endif
                if ((locs(j, 3) - Lz) > 1e-3) then
                    locs(j, 3) = locs(j, 3) - int( locs(j, 3) / Lz) * Lz
                endif
            END DO
            close(1)


            Do ii = 1, dir
                do j = 1, Nsite
                    ms(ii, j) = nnsite1(j, ii + 1)
                end do
            ENDDO
            
            do n1 = 1, Nsite  !check the consistancy
                DO kf=1,dir;kb=backdir1(n1,kf+1);if (ms(kb,ms(kf,n1))/=n1) then
                    print*,n1
                    stop 'nnsite arrange wrong'
                endif;ENDDO
            enddo

      !check the consistancy
            DO n1 = 1, Nsite
               DO kf = 1, dir
                  kb = backdir1(n1, kf + 1)
                  IF (ms(kb, ms(kf, n1)) /= n1) THEN
                     PRINT *, n1
                     STOP 'nnsite arrange wrong'
                  END IF
               END DO
            END DO

     !define the inverse tetrahedra
     !inv_trahe documents the spin component of each inverse tetrahedron
     !inv_sp documents the spin belonging
            DO ii = 1, Lx*Ly*Lz
               n0 = 4 * (ii - 1) + 1
               site0 = trahe(n0)
               inv_trahe(1, ii) = n0
               inv_sp(n0) = ii
               DO kf = 1, dir

                  mm = ms(kf, n0)
                  IF (trahe(mm) /= site0) THEN
                     IF (sublatt(mm) == 2) THEN
                        inv_trahe(2, ii) = mm
                        inv_sp(mm) = ii
                     ELSE IF (sublatt(mm) == 3) THEN
                        inv_trahe(3, ii) = mm
                        inv_sp(mm) = ii
                     ELSE IF (sublatt(mm) == 4) THEN
                        inv_trahe(4, ii) = mm
                        inv_sp(mm) = ii
                     END IF
                  END IF
               END DO
            END DO

        ! check the consistency
            counter = 0
            do ii = 1, Vol
                nsp = inv_sp(ii)
                counter(nsp) = counter(nsp)+1
            enddo

            do Vc = 1, Vol/sublatnum
                
                if (counter(Vc) /= 4)then
                    print*, counter(Vc), Vc
                    stop 'i wrong'
                endif
            enddo

            ! assign the number of the neighbors (initially 6)
            do ii = 1, Vol
                neigh_Num(ii) = dir
            enddo
            !! get rid of the periodic BC at x direction
             ! also change the number of neighbor
            do ii = 1, Vol
                if (locs(ii,1) == 0) then
                    do n1 = 1, dir
                        neigh_spin = ms(n1, ii)
                        if (locs(neigh_spin, 1) > Lx/2) then
                            back_spin = backdir1(ii, n1+1) ! the spin number connected by the bond
                            ms(back_spin, neigh_spin) = 0 !cut the bond
                            ms(n1, ii) = 0 ! 0 means no neighbors
                            neigh_Num(ii) = neigh_Num(ii) - 1
                        endif
                    enddo
                endif
                if (locs(ii,1) == Lx) then
                        do n1 = 1, dir
                            neigh_spin = ms(n1, ii)
                            if (locs(neigh_spin, 1) < Lx/2 ) then
                                back_spin = backdir1(ii, n1+1) ! the spin number connected by the bond
                                ms(back_spin, neigh_spin) = 0 !cut the bond
                                ms(n1, ii) = 0 ! 0 means no neighbors
                                neigh_Num(ii) = neigh_Num(ii) - 1
                            endif
                        enddo
                endif
                if (locs(ii,3) == 0) then
                        do n1 = 1, dir
                            neigh_spin = ms(n1, ii)
                            if (neigh_spin /= 0 )then
                                if (locs(neigh_spin, 3) > Lz/2) then
                                    back_spin = backdir1(ii, n1+1) ! the spin number connected by the bond
                                    ms(back_spin, neigh_spin) = 0 !cut the bond
                                    ms(n1, ii) = 0 ! 0 means no neighbors
                                    neigh_Num(ii) = neigh_Num(ii) - 1
                            endif;endif
                        enddo
                endif
                if (locs(ii,3) == Lz) then
                        do n1 = 1, dir
                            neigh_spin = ms(n1, ii)
                            if (neigh_spin /= 0 )then
                                if (locs(neigh_spin, 3) < Lz/2 ) then
                                    back_spin = backdir1(ii, n1+1) ! the spin number connected by the bond
                                    ms(back_spin, neigh_spin) = 0 !cut the bond
                                    ms(n1, ii) = 0 ! 0 means no neighbors
                                    neigh_Num(ii) = neigh_Num(ii) - 1
                            endif;endif
                        enddo
            
                endif
            enddo
do ii=1,Vol
    if (sum(ms(:,ii))==0) then
        print*, ii, locs(ii,:)
        stop 'neighbor wrong!'
    endif
enddo

            !! take care of the boundary effect
                ! todo
        ELSE
             ! simple cubic lattice
             !directions
                     DO i=1,dim; back(i)=i+dim; back(i+dim)=i; ENDDO
             ! lattice connections
             vector=0
             IF(dim>1) n12=Ns(1)*Ns(2)
                             IF(dim==1) then;  n1=0
             DO ix=1,Ns(1); n1=n1+1
               vector(1,n1)=ix-1
               nn(1)=n1+1; IF(ix==Ns(1)) nn(1)=n1-Ns(1)+1;
               DO kf=1,dim;kb=back(kf);ms(kf,n1)=nn(kf);ms(kb,nn(kf))=n1;ENDDO
             ENDDO
                          ELSE IF(dim==2) then;  n1=0
             DO iy=1,Ns(2); DO ix=1,Ns(1); n1=n1+1
               vector(1,n1)=ix-1; vector(2,n1)=iy-1;
               nn(1)=n1+1;     IF(ix==Ns(1))   nn(1)=n1-Ns(1)+1
               nn(2)=n1+Ns(1); IF(iy==Ns(2))   nn(2)=n1-n12+Ns(1)
               DO kf=1,dim;kb=back(kf); ms(kf,n1)=nn(kf); ms(kb,nn(kf))=n1;ENDDO
               ENDDO; ENDDO
                        ELSE IF(dim==3) then;  n1=0
             DO iz=1,Ns(3); DO iy=1,Ns(2); DO ix=1,Ns(1); n1=n1+1;
               vector(1,n1)=ix-1; vector(2,n1)=iy-1;   vector(3,n1)=iz-1;
               nn(1)=n1+1;     IF(ix==Ns(1))   nn(1)=n1-Ns(1)+1
               nn(2)=n1+Ns(1); IF(iy==Ns(2))   nn(2)=n1-n12+Ns(1)
               nn(3)=n1+n12;   IF(iz==Ns(3))   nn(3)=n1-Nsite+n12
             DO kf=1,dim;kb=back(kf); ms(kf,n1)=nn(kf); ms(kb,nn(kf))=n1;ENDDO
               ENDDO; ENDDO; ENDDO
                        ELSE; stop ' I did not plan for four dimensions'
                        ENDIF
         ENDIF

        !! distance matrix is defined here
        d_spin = 0
        do i = 1, Vol
            do j = 1, Vol
            deltay = locs(i,2) - locs(j,2)
            deltaz = locs(i,3) - locs(j,3)
            ! take care of tye Periodic BC
            if (deltay > Ly/2) then
                deltay = deltay - Ly
            else if (deltay < -Ly/2) then
                deltay = deltay + Ly
            endif
           ! if (deltaz > Lz/2) then
           !     deltaz = deltaz - Lz
           ! else if (deltaz < -Lz/2) then
           !     deltaz = deltaz + Lz
           ! endif
            deltax = locs(i,1) - locs(j,1)
            d_spin(i,j) = dsqrt( deltax * deltax + deltay * deltay + deltaz * deltaz  )

        enddo;enddo

        !! define the sensor location
        do i = 1, num_sensor
            sensor_loc(i, 3) = - 30 ! the sensor is 10nm below the sample
            sensor_loc(i, 2) = Ly / 2 ! the sensor is located at the center of the sample from above
            sensor_loc(i, 1) = (Lx / num_sensor) * (i - 1 + 0.50)!Ly * (i - 1) + Ly / 2.d0
        enddo
print*,sensor_loc(: , 1)

        !! distance between the spin and the sensor
        do i = 1, Vol
            do j = 1, num_sensor
            deltay = locs(i,2) - sensor_loc(j,2)
            deltaz = locs(i,3) - sensor_loc(j,3)
            ! take care of tye Periodic BC
            if (deltay > Ly/2) then
                deltay = deltay - Ly
            else if (deltay < -Ly/2) then
                deltay = deltay + Ly
            endif
        !    if (deltaz > Lz/2) then
        !        deltaz = deltaz - Lz
        !   else if (deltaz < -Lz/2) then
         !       deltaz = deltaz + Lz
        !    endif
            deltax = locs(i,1) - sensor_loc(j,1)
            d_sensor(i, j) = dsqrt( deltax * deltax + deltay * deltay + deltaz * deltaz )

        enddo;enddo

        print*, 'lattice done'
    return
  END SUBROUTINE def_latt
  !===================================================================

  !==========initialize state ======================
  !! THIS IS PROJECT-DEPENDENT

    ! get the magnetization at z direction
    ! this quantity need to be divided by a factor of 1/sqrt(3)
  FUNCTION get_Mz()
    implicit none
    integer(8) :: get_Mz
    integer(8) :: Vc
    get_Mz = 0
    do Vc = 1, Vol
        get_Mz = get_Mz + sp(Vc, 3)
    end do
    return
  END FUNCTION

    ! get the Ising interaction
  FUNCTION get_ens_tot()
    implicit none
    double precision :: get_ens_tot
    integer(8) :: Vc
    integer(1) :: ii

    get_ens_tot = 0
    do Vc = 1, Vol
        get_ens_tot = get_ens_tot - magnet_x * sp(Vc, 1) !take care of the magnetic field at X direction

        do ii = 1, dir
            if (ms(ii, Vc) /= 0) then
                
                get_ens_tot = get_ens_tot + J_coup * ( sp(Vc, 1) * sp(ms(ii, Vc),1) + sp(Vc, 2) &
                    * sp(ms(ii, Vc),2) + sp(Vc, 3) * sp(ms(ii, Vc),3)) /2 / 3.d0!divided by 2 because we add up twice
            endif
        end do
    end do
  END FUNCTION

      ! get the Ising interaction
    FUNCTION get_ens()
      implicit none
      double precision :: get_ens(Vol)
      integer(8) :: Vc
      integer(1) :: ii

      get_ens = 0
      do Vc = 1, Vol
          get_ens(Vc) = get_ens(Vc) - magnet_x * sp(Vc,1) !take care of the magnetic field at X direction

          do ii = 1, dir
              if (ms(ii, Vc) /= 0) then
                  get_ens(Vc) = get_ens(Vc) + J_coup * ( sp(Vc, 1) * sp(ms(ii, Vc),1) + sp(Vc, 2) &
                      * sp(ms(ii, Vc), 2) + sp(Vc, 3) * sp(ms(ii, Vc), 3)) /3.d0
              endif
          end do
      end do
    END FUNCTION
    
    ! get the Coulomb interaction O(N**2/2)
  FUNCTION get_enc_tot()
    implicit none
    double precision :: get_enc_tot
    integer(8) :: Vc, ii, j
    double precision :: amp, amp1, amp2
    real(8) :: r(3)

    get_enc_tot = 0
    amp1 = 0
    do Vc = 1, Vol
        do ii = Vc + 1, Vol

            do j = 1, 3
                r(j) = locs(Vc, j) - locs(ii, j)
            enddo

            if (r(2) > Ly / 2.d0) then
                r(2) = r(2) - Ly
            elseif(r(2) < -Ly / 2.d0) then
                r(2) = r(2) + Ly
            endif

            if (d_spin(Vc, ii) == 0) then
                print*, 'distance wrong',Vc, ii
                stop
            endif

            amp1 = ( sp(Vc, 1) * r(1) + sp(Vc, 2) * r(2) + sp(Vc, 3) * r(3) ) &
                * ( sp(ii, 1) * r(1) + sp(ii, 2) * r(2) + sp(ii, 3) * r(3) )
            amp1 = amp1 / d_spin(Vc, ii)**5 /3.d0  / dsqrt(2.d0)**3
            amp2 = sp(Vc, 1) * sp(ii, 1) + sp(Vc, 2) * sp(ii, 2) + sp(Vc, 3) * sp(ii, 3)
            amp2 = amp2 / d_spin(Vc, ii)**3 /3.d0  / dsqrt(2.d0)**3
            amp = amp2 - 3 * amp1
 
            get_enc_tot = get_enc_tot + Dipo * amp
        end do
    end do
  END FUNCTION

    ! get the Coulomb interaction O(N**2/2)
     FUNCTION get_enc()
       implicit none
       double precision :: energy2(Vol, Vol), get_enc(Vol)
       integer(8) :: Vc, ii, j
       double precision :: amp, amp1, amp2
       real(8) :: r(3)

       get_enc = 0
       energy2 = 0
       do Vc = 1, Vol
           do ii = Vc + 1, Vol

               do j = 1, 3
                   r(j) = locs(Vc, j) - locs(ii, j)
               enddo

                if (r(2) > Ly / 2.d0) then
                    r(2) = r(2) - Ly
                elseif(r(2) < -Ly / 2.d0) then
                    r(2) = r(2) + Ly
                endif

               amp1 = ( sp(Vc,1) * r(1) + sp(Vc,2) * r(2) + sp(Vc,3) * r(3) ) &
                * ( sp(ii,1) * r(1) + sp(ii,2) * r(2) + sp(ii,3) * r(3) )
               amp1 = amp1 / d_spin(Vc, ii)**5 /3.d0 / dsqrt(2.d0)**3
               amp2 = sp(Vc,1) * sp(ii,1) + sp(Vc,2) * sp(ii,2) + sp(Vc,3) * sp(ii,3)
               amp2 = amp2 / d_spin(Vc, ii)**3 / 3.d0 / dsqrt(2.d0)**3
               amp = amp2 - 3 * amp1
   ! if (Vc == 6418 .and. ii == 6429 )then
   !     print*, r(1), r(2), r(3)
   !     print*,'time to gg ....', locs(Vc, 2) , locs(ii, 2)
   !     print*, amp1, amp2, amp
   !     print*,'gte ythe energy',Dipo * amp
   ! endif
               energy2(Vc, ii) = Dipo * amp
               energy2(ii, Vc) = energy2(Vc, ii)
           end do
       end do
!print*, 'yejoo', energy2(6418,6429)
        do Vc = 1, Vol
            get_enc(Vc) = sum (energy2(Vc, :))
        enddo
     END FUNCTION

    ! get the magnetic field outside the material
    FUNCTION get_Bz()
      implicit none
      double precision :: get_Bz(num_sensor, 3)
      integer(8) :: Vc, j, ii, spin_count(5)
      double precision :: amp, Bz, Bx, By
      double precision :: r(3)

      get_Bz = 0
      spin_count=0

    do Vc = 1, Vol
        do ii = 1, num_sensor
            do j = 1, 3
                r(j) = locs(Vc, j) - sensor_loc(ii, j)
            enddo
            
            if (r(2) > Ly / 2.d0) then
                r(2) = r(2) -Ly
            elseif(r(2) < -Ly / 2.d0) then
                r(2) = r(2) + Ly
            endif

            !if (locs(Vc, 1) >= (ii-1) * Ly .and.  locs(Vc, 1) <= ii * Ly) then
                spin_count(ii) = spin_count(ii) + 1

                amp = sp(Vc,1) * r(1) + sp(Vc,2) * r(2) + sp(Vc,3) * r(3)
                Bz = 3 * r(3) * amp / dsqrt(3.d0) / d_sensor(Vc, ii) **5 - sp(Vc, 3)/ dsqrt(3.d0) / d_sensor(Vc, ii) **3
                get_Bz(ii, 3) = get_Bz(ii, 3) + Bz !/dsqrt(3.d0)
                
                Bx = 3 * r(1) * amp / dsqrt(3.d0) / d_sensor(Vc, ii) **5 - sp(Vc, 1)/ dsqrt(3.d0) / d_sensor(Vc, ii) **3
                get_Bz(ii, 1) = get_Bz(ii, 1) + Bx !/dsqrt(3.d0)
    
                By = 3 * r(2) * amp / dsqrt(3.d0) / d_sensor(Vc, ii) **5 - sp(Vc, 2)/ dsqrt(3.d0) / d_sensor(Vc, ii) **3
                get_Bz(ii, 2) = get_Bz(ii, 2) + By !/dsqrt(3.d0)

            !ENDIF

      enddo;enddo
!print*, 'spin_count',spin_count
      return
    END FUNCTION


        ! get the monopole distribution in each tetrahegen
        FUNCTION get_mono1()
        implicit none
        integer(8) :: get_mono1, break0(Vol/sublatnum)
        integer(8) :: Vc, jj, ii
        integer(8) :: sub, nsp

        get_mono1 = 0
        break0 = 0
        
        
        do Vc = 1, Vol/sublatnum
            
            jj = (Vc - 1)* sublatnum
            do ii = 1, sublatnum
                sub = sublatt(jj + ii)
                break0(Vc) = break0(Vc) + sp(jj + ii, 1) / Pat(sub)
            enddo
        enddo

        do Vc = 1, Vol/sublatnum
            if (break0(Vc) /= 0)then
                get_mono1 = get_mono1 + 1
            endif
        end do
        END FUNCTION

    ! get the monopole distribution at each inverse tetrahegen
        FUNCTION get_mono2()
        implicit none
        integer(8) :: get_mono2, break_inv0(Vol/sublatnum)
        integer(8) :: counter(Vol/sublatnum)
        integer(8) :: Vc1, Vc
        integer(8) :: sub, nsp

        get_mono2 = 0
        break_inv0 = 0
        counter = 0
        do Vc1 = 1, Vol
            
            sub = sublatt(Vc1)
            nsp = inv_sp(Vc1)

            break_inv0(nsp) = break_inv0(nsp) - sp(Vc1, 1) / Pat(sub)
            counter(nsp) = counter(nsp) + 1
        enddo


        do Vc = 1, Vol/sublatnum
            if (break_inv0(Vc) /= 0)then
                get_mono2 = get_mono2 + 1
            endif
        end do

! check the consistency
        do Vc = 1, Vol/sublatnum
            
            if (counter(Vc) /= 4)then
                print*, counter, Vc
                stop 'inverse lattice wrong'
            endif
        enddo

        END FUNCTION

    ! get the monopole density
        FUNCTION get_density()
        implicit none
        double precision :: get_density(Lx, 2)
        integer(8) :: Vc, jj, ii
        integer(8) :: sub, nsp
        integer(8) :: fl_int1, fl_int2

        get_density = 0

            do jj = 1, Vol/sublatnum
                fl_int1 = cent(jj, 1) + 1 !+ 0.00001;
                fl_int2 = cent_inv(jj, 1) + 1  !+ 0.0001;
                if (break(jj) == 2)then
                    
                    get_density(fl_int1, 1) = get_density(fl_int1, 1) + 1
                    
                 elseif (break(jj) == -2)then
                    
                    get_density(fl_int1, 2) = get_density(fl_int1, 2) + 1

                 endif

                if (break_inv(jj) == 2)then
                   
                    get_density(fl_int2, 1) = get_density(fl_int2, 1) + 1
                    
                 elseif (break_inv(jj) == -2)then
                    
                    get_density(fl_int2, 2) = get_density(fl_int2, 2) + 1
                    
                 endif
            enddo


        END FUNCTION


      SUBROUTINE init_stat
        implicit none
        integer(8) :: Vc,Vcn, sub
        integer(8) :: ii, jj, ll, nsp, density_index(2), B2P(35,2), spin_num(4)

        Mz = get_Mz()

        enC_tot = get_enc_tot()

        enS_tot = get_ens_tot()
        enS = get_ens()
        enC = get_enc()
        Bz_sensor = get_Bz()
        mono1 = get_mono1()
        mono2 = get_mono2()
        
    ! initialize the tetrehedra structure
        break = 0
        break_inv = 0
        cent = 0
        cent_inv = 0
        do Vc = 1, Vol/sublatnum
            
            jj = (Vc - 1)* sublatnum
            do ii = 1, 4
                sub = sublatt(jj + ii)
                break(Vc) = break(Vc) + sp(jj + ii, 1) / Pat(sub)
                do ll = 1, 3
                    cent(Vc, ll) = cent(Vc, ll) + locs(jj + ii, ll)
                enddo
            enddo
        enddo
        cent = cent/4.d0 ;

   !     do Vc = 1, Vol
   !
   !         sub = sublatt(Vc)
   !         nsp = inv_sp(Vc)
   !         break_inv(nsp) = break_inv(nsp) - sp(Vc, 1) / Pat(sub)
   !         do ll = 1, 3
   !             cent_inv(nsp, ll) = cent_inv(nsp, ll) + locs(Vc, ll)
   !         enddo
   !         cent_inv(nsp, :) = cent_inv(nsp, :) / 4.d0;
   !     enddo

!correct the cross boundary tetrahegons
        do nsp = 1, Vol/sublatnum
            
            spin_num = inv_trahe(:, nsp)
            do Vc = 1, sublatnum
                if (sublatt(spin_num(Vc)) == 3) then
                    do ll = 1, 3
                        cent_inv(nsp, ll) = locs(spin_num(Vc), ll)
                    enddo
                endif
            enddo
        enddo
        
! check the anti-monopoles are evenly distributed in the lattice
B2P = 0
do Vc = 1, Vol/sublatnum
        density_index(1) = cent(Vc, 1) + 1;
        density_index(2) = cent_inv(Vc, 1) + 1;
        B2P(density_index(1), 1) = B2P(density_index(1), 1) + 1
        B2P(density_index(2), 2) = B2P(density_index(2), 2) + 1
enddo
!print*,'cent_inv~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!print*, B2P(:,2)



    !density distribution of the monopole and antimonopole
        density = get_density()
        kw=0
        wormstep=0.d0;  step=0


        return
      END SUBROUTINE init_stat
  
  !==================== END INITIALIZATION =====================

  !==================== Update module during KMC =====================
    ! This is project dependent

    !! update the spin state
        SUBROUTINE up_spin(targ)
          implicit none
          integer :: i, j, jj, ll, loc_
          integer(8) ::  targ, broken, broken_inv, nsp, density_index(2)
              
            do i = 1, 3
                sp(targ , i) = - sp(targ , i)
            enddo
            

            I = sublatt( targ )
            ll = mod(targ, 4)
            if (ll == 0) then
                jj = targ /4
            else
                jj = targ /4 + 1
            endif

            density_index(1) = cent(jj, 1) + 1;
!print*, 'density_index(1)', density_index(1)
            broken = break(jj)
            break(jj) = break(jj) + 2 * sp(targ ,1) / Pat(I)
            if (break(jj) /= 0 .and. broken == 0 ) then
                mono1 = mono1 + 1
                
            elseif (break(jj) == 0 .and. broken /=0 )then
                mono1 = mono1 - 1
                
            endif

            nsp = inv_sp(targ)
            density_index(2) = cent_inv(nsp, 1) + 1 ;
            broken_inv = break_inv(nsp)

            break_inv(nsp) = break_inv(nsp) - 2 * sp(targ ,1) / Pat(I)
!print*, 'density_index(2)', density_index(2), broken_inv, break_inv(nsp )
            if (break_inv(nsp) /= 0 .and. broken_inv == 0 ) then
                mono2 = mono2 + 1
            elseif (break_inv(nsp) == 0 .and. broken_inv /=0 )then
                mono2 = mono2 - 1
            endif

            if (break(jj) == 2 .and. broken /= 2 )then
                density(density_index(1), 1) = density(density_index(1), 1) + 1
            elseif(break(jj) /= 2 .and. broken == 2)then
                density(density_index(1), 1) = density(density_index(1), 1) - 1
            endif

            if (break_inv(nsp) == 2 .and. broken_inv /= 2 )then
                density(density_index(2), 1) = density(density_index(2), 1) + 1
            elseif(break_inv(nsp) /= 2 .and. broken_inv == 2)then
                density(density_index(2), 1) = density(density_index(2), 1) - 1
            endif

            if (break(jj) == -2 .and. broken /= -2 )then
                density(density_index(1), 2) = density(density_index(1), 2) + 1
            elseif(break(jj) /= -2 .and. broken == -2)then
                density(density_index(1), 2) = density(density_index(1), 2) - 1
            endif

            if (break_inv(nsp) == -2 .and. broken_inv /= -2 )then
                density(density_index(2), 2) = density(density_index(2), 2) + 1
            elseif(break_inv(nsp) /= -2 .and. broken_inv == -2)then
                density(density_index(2), 2) = density(density_index(2), 2) - 1
            endif

          return
        END SUBROUTINE up_spin


    !! calculate the energy change in NN term
      SUBROUTINE up_energy_nn(targ)
        implicit none
        integer :: i, j, spin
        integer(8) :: targ
        double precision :: amp

            enS(targ) = - enS(targ)
            enS_tot = enS_tot + 2 * enS(targ)
            do i = 1, dir
                if (ms(i, targ) /= 0) then

                    spin = ms(i, targ)
!print*,'now update',ms(i, targ), enS(spin)
                    amp = sp(spin,1) * sp(targ,1) + sp(spin,2) * sp(targ,2) + sp(spin,3) * sp(targ,3)
                    amp = + J_coup * amp / 3.d0 !/3.d0
                    enS(spin) = enS(spin) + 2 * amp
!print*,'after update', enS(spin)
                    !enS_tot = enS_tot + amp
                endif
                
            enddo
            
        return
      END SUBROUTINE up_energy_nn

    !! calculate the energy change in Dipolar term O(N)
      SUBROUTINE up_energy_coulomb(targ)
        implicit none
        integer :: i, j
        integer(8) :: targ
        double precision :: amp1, amp2, amp
        double precision :: r(3)

        do i = 1, Vol
            if (i /= targ)then
                do j = 1, 3
                    r(j) = locs(i,j) - locs(targ,j)
                enddo

                if (r(2) > Ly / 2.d0) then
                    r(2) = r(2) - Ly
                elseif(r(2) < -Ly / 2.d0) then
                    r(2) = r(2) + Ly
                endif
!if (i ==6429)then
!    print*, enc(i), locs(i,2) , locs(targ,2)
!endif
                amp1 = ( sp(i,1) * r(1) + sp(i,2) * r(2) + sp(i,3) * r(3) ) &
                * ( sp(targ,1) * r(1) + sp(targ,2) * r(2) + sp(targ,3) * r(3) )
                amp1 = amp1 / d_spin(i,targ)**5 /3.d0 / dsqrt(2.d0)**3
                amp2 = sp(i,1) * sp(targ,1) + sp(i,2) * sp(targ,2) + sp(i,3) * sp(targ,3)
                amp2 = amp2 / d_spin(i,targ)**3 /3.d0 / dsqrt(2.d0)**3
                amp = amp2 - 3 * amp1
                enC(i) = enC(i) + 2 * amp * Dipo
                enC_tot = enC_tot + amp * Dipo
            else
                enC(i) = - enC(i)
                    enC_tot = enC_tot + enC(i)
            endif
!if (i ==6429)then
 !   print*, r(1), r(2), r(3)
 !   print*, amp1, amp2, amp
!    print*, 'energy updated', amp * Dipo
!endif
        enddo

        return
      END SUBROUTINE up_energy_coulomb

    !! update the magnetization at z direction
        SUBROUTINE up_Mz(targ)
          implicit none
          integer :: ii, j
          integer(8) :: targ
        
            Mz = Mz + 2 * sp(targ, 3)

          return
        END SUBROUTINE up_Mz

    !! update the magnetic field outside
        SUBROUTINE up_Bz(targ)
          implicit none
          integer :: ii, j
          integer(8) :: targ
          double precision :: amp, ene1, ene2, ene3
          double precision :: r(3)

            do ii = 1, num_sensor
                do j = 1, 3
                      r(j) = locs(targ, j) - sensor_loc(ii, j)
                enddo

                if (r(2) > Ly / 2.d0) then
                    r(2) = r(2) - Ly
                elseif(r(2) < -Ly / 2.d0) then
                    r(2) = r(2) + Ly
                endif

                !if (locs(targ, 1) >= (ii-1) * Ly .and.  locs(targ, 1) <= ii * Ly) then
                      amp = sp(targ,1) * r(1) + sp(targ,2) * r(2) + sp(targ,3) * r(3)
                      ene3 = 3 * r(3) * amp /dsqrt(3.d0) / d_sensor(targ, ii) **5 - sp(targ,3)/ dsqrt(3.d0) / d_sensor(targ, ii)**3
                      Bz_sensor(ii, 3) = Bz_sensor(ii, 3) + 2 * ene3 !/ dsqrt(3.d0)

                    ene1 = 3 * r(1) * amp /dsqrt(3.d0) / d_sensor(targ, ii) **5 - sp(targ,1)/ dsqrt(3.d0) / d_sensor(targ, ii)**3
                    Bz_sensor(ii, 1) = Bz_sensor(ii, 1) + 2 * ene1 !/ dsqrt(3.d0)

                    ene2 = 3 * r(2) * amp /dsqrt(3.d0) / d_sensor(targ, ii) **5 - sp(targ,2)/ dsqrt(3.d0) / d_sensor(targ, ii)**3
                    Bz_sensor(ii, 2) = Bz_sensor(ii, 2) + 2 * ene2 !/ dsqrt(3.d0)

                !ENDIF

            enddo

          return
        END SUBROUTINE up_Bz
  
  !====================== START MARKOV =========================
  
  
  !! THIS IS PROJECT-DEPENDENT
  SUBROUTINE markov
    implicit none
    integer :: istep
    do istep= 1, N_measure

        !call check

        if (prtflg) then
            call KMC()
        else
            call worm()
            !if (accept0) then
            !    call check()
            !endif
        endif
    enddo

    return
  END SUBROUTINE markov


  !==============KMC ==========================================
  !! THIS IS PROJECT-DEPENDENT

  SUBROUTINE KMC
    implicit none
    real(8) :: p
    integer(8) :: I
    integer :: s

    call def_prob()
    call accumulate_prob()
    call search_prob( rn()*Accu_prob )
    call record_time()

    !update the quantities
    call up_spin( flipspin )
    call up_energy_nn( flipspin )
    call up_energy_coulomb( flipspin )
    call up_Mz( flipspin )
    call up_Bz( flipspin )

    !print*, 'one step KMC', Accu_prob, flipspin
   ! if (Accu_prob >1e6) then
    !    stop 'Prob wrong!'
    !endif

    return
  END SUBROUTINE

  SUBROUTINE worm
    implicit none
    real(8) :: p
    integer(8) ::  I

    I = dint(Vol * rn()) + 1

    dE = - (2 * enS(i) + 2 * enC(i)) !- 2 * magnet_x * sp(I, 1)
    p = dexp( - beta * dE ) !detailed balance
!if (prtflg)then
!print*, dE, ens(i), enc(i),p,rn()
!print*, mono
!endif
    if(p > 1 .or. rn() < p) then
        !print*,'accepted', I
        call up_spin( I )
        call up_energy_nn( I )
        call up_energy_coulomb( I )
        call up_Mz( I )
        call up_Bz( I )
        accept0 = .true.
        accept = accept + 1
    else
        accept0 = .false.
    end if

!if (accept0) then
!print*,'spin number', I, sp(I, :)
!print*, 'location' , locs(I, :), cent(I/4 + 1, 1), cent(I/4, 1)
!print*, 'compose',inv_trahe(:, 1)
!endif
    return
  END SUBROUTINE


  subroutine check()
    implicit none
    integer(8) :: Mz0
    double precision :: Es0, Ec0, Bz0(num_sensor, 3)
    double precision :: each_Es0(Vol), each_Ec0(Vol)
    integer :: ii, jj
    integer(8) :: mono10, mono20, density0(Lx,2)

    print*, "checking ..."

    Es0 = get_ens_tot()
    Ec0 = get_enc_tot()
    Mz0 = get_Mz()
    Bz0 = get_Bz()
    each_Ec0 = get_enc()
    each_Es0 = get_ens()
    mono10 = get_mono1()
    mono20 = get_mono2()
    density0 = get_density()


    !check each state of Spin-spin energy
        do ii = 1, Vol
            if (dabs(each_Es0(ii) - Ens(ii)) > 1.d-5) then
                print*, 'right', each_Es0(ii)
                print*, 'current', enS(ii)
                print*, 'spin number', ii
                stop "wrong ens!"
            end if
        end do

    !check each state of coulomb energy
        do ii = 1, Vol
            if (dabs(each_Ec0(ii) - enc(ii)) > 1.d-5) then
                print*, each_Ec0(ii), enc(ii)
                print*, 'spin number', ii
                stop "wrong enC!!"
            end if
        enddo

    do jj =1, 2
        do ii = 1, Lx
        if (dabs(density0(ii, jj) - density(ii, jj)) > 1.d-5) then
            print*, density0(ii, jj), density(ii, jj)
            print*, 'tetra number', ii, jj
            stop "wrong density!!"
        end if
    enddo; enddo

    if(abs(mono10 - mono1) > 1.d-5) then
        print*,mono1, mono10
        stop "wrong monopole1!"
    end if
    if(abs(mono20 - mono2) > 1.d-5) then
        print*,mono2, mono20
        stop "wrong monopole2!"
    end if

    if(dabs(enS_tot - Es0) > 1.d-5) then
        print*,enS_tot, Es0
        stop "wrong Spin-spin!"
    end if
    if(dabs(enC_tot - Ec0) > 1.d-5) then
        print*,enC_tot, Ec0
        stop "wrong Coulomb!"
    end if
    if(Mz0 /= Mz) then
        print*, Mz, Mz0
        stop "wrong Mz!"
    end if
    do ii = 1, num_sensor
        do jj =1 , 3
            if (dabs(Bz0(ii, jj) - Bz_sensor(ii, jj)) > 1.d-5) then
                print*, Bz0(ii, jj), Bz_sensor(ii, jj), jj
                stop "wrong Bz!"
            end if
        enddo
    end do

    print*, "check done"

    return
  end subroutine
    
  !====================== END MARKOV =========================


  !================== START MEASUREMENT ======================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE measure
    implicit none
    integer :: ii


!    i1=i1+1
!    if(i1>=kw) then
!        call sigma()
!        i1=0
!    end if



	!call check()


    Obser(1)%vnr = mono1
!if(abs(mono) > 100) then
!    print*,mono
!endif
    Obser(2)%vnr = mono2
    Obser(3)%vnr = enS_tot
    Obser(4)%vnr = current_t
    do ii = 1, num_sensor
        Obser(4 + ii)%vnr = Bz_sensor(ii, 3)
        Obser(4 + num_sensor + ii)%vnr = Bz_sensor(ii, 1)
        Obser(4 + 2 * num_sensor + ii)%vnr = Bz_sensor(ii, 2)
    enddo


    !if(prtflg .and. thermal )then
    if( prtflg )then
        call coll_data()
        sequce = sequce + 1
        step = step + 1
        if ( step >= 400) then
            if (apply) then
            !call shot_monopole()
                call shot_density1()
            else
                call shot_density2()
            endif
            
            step = 0
     !   if (sequce == Max_time)then
     !       call writejjfile()
      !      sequce = 1
     !       timesequce = 0
        endif
    end if;! endif

  !  if(prtflg) then
   !     call get_equil()
   ! end if




!    wormstep=wormstep+step
!    step=0

    return
  END SUBROUTINE measure
  



  SUBROUTINE sigma()
    implicit none
    integer :: i,k,n,it

    return
  END SUBROUTINE

  !================== END MEASUREMENT =======================
  
  !============== START WRITE TO FILES =====================================
  !! THIS IS PROJECT-DEPENDENT 
  SUBROUTINE writejjfile1
    implicit none
    integer     :: j
    real(8)     :: xx

      !  open(10,file = 'magnet1.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
      !!  write(10, *) 'path'
       ! do j = 1, Max_time
      !      write(10, *) j, timesequce(4,j), timesequce(1,j)
      ! enddo
       ! close(10)

        open(10,file = 'field1.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
        write(10, *) 'path'
        do j = 1, Max_time
            write(10, *) j, timesequce(4,j), timesequce(5 : 4 + 3 * num_sensor,j)
        enddo
        close(10)
    print*, 'write the path done------------------------------'

    return
  END SUBROUTINE writejjfile1

    SUBROUTINE writejjfile2
      implicit none
      integer     :: j
      real(8)     :: xx

          open(10,file = 'field2.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
          write(10, *) 'path'
          do j = 1, Max_time
              write(10, *) j, timesequce(4,j), timesequce(5 : 4 + 3 * num_sensor,j)
          enddo
          close(10)
      print*, 'write the path done------------------------------'

      return
    END SUBROUTINE writejjfile2


SUBROUTINE shot_monopole
  implicit none
  integer     :: j
  real(8)     :: xx

      open(10,file = 'monopole.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
      write(10, *) 'shot'

      write(10, *) current_t
      do j = 1, Vol/sublatnum
        write(10, *) cent(j,1), cent(j,2), cent(j,3), break(j)
      enddo
      do j = 1, Vol/sublatnum
        write(10, *) cent_inv(j,1), cent_inv(j,2), cent_inv(j,3), break_inv(j)
      enddo

      close(10)

  !print*, 'shot the monopole distribution done------------------------------'

  return
END SUBROUTINE shot_monopole

    SUBROUTINE shot_density1
      implicit none
      integer     :: j
      real(8)     :: xx

          open(10,file = 'density_mono1.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
          write(10, *) 'mono',current_t
          do j = 1, Lx
            write(10, *) (j-1)+0.5, density(j,1)*1.d0/Vol
          enddo

          close(10)
            open(10,file = 'density_anti1.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
            write(10, *) 'anti', current_t
            do j = 1, Lx
              write(10, *) (j-1)+0.5, density(j,2)*1.d0/Vol
            enddo

            close(10)
            

      !print*, 'shot the density done------------------------------'

      return
    END SUBROUTINE shot_density1

SUBROUTINE shot_density2
  implicit none
  integer     :: j
  real(8)     :: xx

      open(10,file = 'density_mono2.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
      write(10, *) 'mono',current_t
      do j = 1, Lx
        write(10, *) (j-1)+0.5, density(j,1)*1.d0/Vol
      enddo

      close(10)
        open(10,file = 'density_anti2.txt', access = "append")!,form='UNFORMATTED',access='sequential',status="replace")
        write(10, *) 'anti', current_t
        do j = 1, Lx
          write(10, *) (j-1)+0.5, density(j,2)*1.d0/Vol
        enddo

        close(10)
        

  !print*, 'shot the density done------------------------------'

  return
END SUBROUTINE shot_density2

  SUBROUTINE printing(id)
    implicit none
    integer :: id,i
    real(8) :: temp
    write(id,*) "========================================"
    write(id,*) "Now	:"
    write(id,*) "Z=",Npb*(Nbk-1)+Np
    do i=1,NObs_b
        write(id,*) i,trim(Obser(i)%nam)//'	',Obser(i)%vnr*Obser(i)%a+Obser(i)%b
    end do
    if(prtflg) then
        write(id,*) "Average:"
        write(id,*) "ZZ=",Npb*(Nbk-1)
        do i=1,NObs
            write(id,*) i,'<'//trim(Obser(i)%nam)//'>	',Obser(i)%val*Obser(i)%a+Obser(i)%b,Obser(i)%erb*Obser(i)%a,Obser(i)%cor
        end do
        write(id,*) "flg_cor=",flg_cor,"Nbk=",Nbk!"Npb=",Npb,"Nbk=",Nbk,"Np=",Np
        write(id,*) "wormstep=",wormstep
        write(id,*) "simulation time :",t_simu
    else
        write(id,*) "therm step:",itoss * N_each
        write(id,*) "therm time:",t_toss
        write(id,*) "accepted step:",accept
    end if
!    write(id,*) "measure time    :",t_meas
    write(id,*) "========================================"
    return
  END SUBROUTINE


  SUBROUTINE write_cnf()
    implicit none
    integer(8) :: i,j
    open(10,file=trim(cnf_file), form='UNFORMATTED',access = 'sequential',status="replace")
    write(10) D, Lx, Ly, Lz
    write(10) sp(1:Vol, :)
    close(10)

!    print*,rtal

    return
  END SUBROUTINE

  SUBROUTINE read_cnf()
    implicit none
    integer(1) :: D1
    integer(8) :: Lx1, ly1, Lz1, i, j

    open(1,file=trim(cnf_file),form='UNFORMATTED',access = 'sequential',status="old")
    read(1) D1, Lx1, Ly1, Lz1
    if(D1 /= D .or. Lx1 /= Lx .or. Ly1 /= Ly .or. Lz1 /= Lz ) stop "the cnf is not suitable for the mismatched Ly."
    read(1) sp(1:Vol, :)
    close(1)
    return
  END SUBROUTINE

  SUBROUTINE write_stat()
    implicit none
    open(10,file=trim(stat_file),form='UNFORMATTED',access='sequential',status="replace")
    write(10) D, Lx, Ly, Lz, J_coup, Dipo, h, beta, N_measure, Seed
    write(10) i1,i2,prtflg,t_simu,t_meas,step,wormstep
    write(10) Npb,Nbk,Np
    write(10) Obser(1:NObs)
    write(10) ir1,ir2,ipnt1,ipnf1,ipnt2,ipnf2,irn,nrannr
    close(10)
    return
  END SUBROUTINE

  SUBROUTINE read_stat()
    implicit none
    integer(1) :: D1
    integer(8) :: Lx1, Ly1, Lz1, N_measure1
    integer :: Seed1
    real(8) :: J1, Dipo1, h1, beta1
    open(1,file=trim(stat_file),form='UNFORMATTED',access='sequential',status="old")
    read(1) D1, Lx1, Ly1, Lz1, J1, Dipo1, h1, beta1, N_measure1, Seed1
    if(D1/=D) stop "the stat is not suitable for the mismatched Dimension."
    if(Lx1/=Lx) stop "the stat is not suitable for the mismatched Lx."
    if(N_measure1/=N_measure)  stop "the stat is not suitable for the mismatched N_measure."
    if(Seed1/=Seed) stop "the stat is not suitable for the mismatched Seed."
    read(1) i1,i2,prtflg,t_simu,t_meas,step,wormstep
    read(1) Npb,Nbk,Np
    read(1) Obser(1:NObs)

    read(1) ir1,ir2,ipnt1,ipnf1,ipnt2,ipnf2,irn,nrannr
    close(1)
    return
  END SUBROUTINE

  !====================== END WRITE TO FILES =============================================
  
  
  !================= Obs ========================================

  SUBROUTINE set_Obs
    implicit none
    integer :: i, ii
    call init_Obs(1,"mono1    ",1.d0/Vol,0.d0)
    call init_Obs(2,"mono2    ",1.d0/Vol,0.d0)
    call init_Obs(3,"energy2 ",1.d0/Vol,0.d0)
    call init_Obs(4,"time    ",1.d0,0.d0)
    do i = 1, num_sensor
        call init_Obs(4 + i,"B_z     ",1.d0,0.d0)
    enddo
    do ii = 1, num_sensor
        call init_Obs(4 + 5 + ii,"B_x     ",1.d0,0.d0)
    enddo
    do ii = 1, num_sensor
        call init_Obs(4 + 10 + ii,"B_y     ",1.d0,0.d0)
    enddo



    Npb=1
    Nbk=1
    Np=0
    sequce = 1
    return
  END SUBROUTINE
  
  SUBROUTINE init_Obs(i,nam0,a0,b0)
    implicit none
    integer :: i
    character(8) :: nam0
    real(8) :: a0,b0
    Obser(i)%nam=trim(nam0)
    Obser(i)%a=a0
    Obser(i)%b=b0
    Obser(i)%vnr=0.d0
    Obser(i)%val=0.d0
    Obser(i)%cor=1.d100
    Obser(i)%erb=1.d100
    Obser(i)%blk=0.d0
    return
  END SUBROUTINE
  
  SUBROUTINE coll_data()
    implicit none
    integer :: i,j

    do i = 1, NObs
        !Obser(i)%blk(Nbk)=Obser(i)%blk(Nbk)+Obser(i)%vnr
        timesequce(i, sequce) = Obser(i)%vnr
    end do
    return
  END SUBROUTINE

    ! get the equilibirum state expactation value and variance

    SUBROUTINE get_equil()
      implicit none
      integer :: i,j

      if(Np==Npb) then
          if(flg_cor) then
!print*,'here'
              if(Nbk==Max_block) then

                  call merge_blk()

              end if
          else
              if(Nbk>=N_block .and. Nbk/2*2==Nbk) then
         !         call cal_ComObs()
                  if(cal_cor_dev()) then
                      flg_cor=.True.
                  else
                      call merge_blk()
                  end if
              end if
          end if
          Nbk=Nbk+1
          Np=0
      end if
      Np=Np+1
      do i=1, NObs_b
          Obser(i)%blk(Nbk)=Obser(i)%blk(Nbk)+Obser(i)%vnr
      end do
      return
    END SUBROUTINE
  
  SUBROUTINE merge_blk()
    implicit none
    integer :: i,j
    Nbk=Nbk/2
    Npb=2*Npb
    do i=1,NObs_b
        do j=1,Nbk
            Obser(i)%blk(j)=Obser(i)%blk(2*j-1)+Obser(i)%blk(2*j)
        end do
        Obser(i)%blk(Nbk+1:2*Nbk)=0.d0
    end do
    return
  END SUBROUTINE
  
  FUNCTION cal_cor_dev()
    implicit none
    logical :: cal_cor_dev
    integer :: i,iblck,Nbk0
    real :: cor,dev,devp,devn
    cal_cor_dev=.True.
    if(Np/=Npb) then
        Nbk0=Nbk-1
    else
        Nbk0=Nbk
    end if
    do i=1,NObs
        Obser(i)%val=Sum(Obser(i)%blk(1:Nbk0))/Nbk0
        cor=0.d0;   dev=0.d0;   devp=0.d0
        do iblck=1,Nbk0
            devn=Obser(i)%blk(iblck)-Obser(i)%val
            dev=dev+devn*devn
            cor=cor+devn*devp
            devp=devn
        end do
        Obser(i)%cor=cor/dev
        Obser(i)%val=Obser(i)%val/Npb
        Obser(i)%erb=dsqrt(dev/Nbk0/(Nbk0-1.d0))/Npb
        if(dabs(Obser(i)%cor)>Max_cor) then
            cal_cor_dev=.False.
        end if
    end do
    return
  END FUNCTION
  
  SUBROUTINE cal_ComObs()
    implicit none
    integer :: iObs(NObs),N
    real(8) ,external :: spec_heat,binder,suscept   !define variable name
    iObs(1)=1;  iObs(2)=2;  N=2                    !select what kind of basic variable is needed
    call cal_ComObs_func(NObs_b+1,spec_heat,iObs,N)
    iObs(1)=4;  iObs(2)=5;  N=2
    call cal_ComObs_func(NObs_b+2,binder,iObs,N)
    iObs(1)=3;  iObs(2)=4;  N=2
    call cal_ComObs_func(NObs_b+3,suscept,iObs,N)
    return
  END SUBROUTINE

  SUBROUTINE cal_ComObs_func(indx,func,iObs,N)
        implicit none
        integer :: indx,N,iblck,i
        integer :: iObs(N)
        real(8) :: a(N)
        real(8) ,external :: func
        do iblck=1,Nbk
            do i=1,N
                a(i)=Obser(iObs(i))%blk(iblck)/Npb
            end do
            Obser(indx)%blk(iblck)=func(a,N)*Npb
        end do
        return
  END SUBROUTINE


  !===============Shift register random number generator =============
  !  very long period sequential version
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_RNG
    implicit none
    integer :: i_r,k_r,k1_r
    integer :: iseed

    nrannr = mxrn
    iseed  = iabs(Seed)+1
    k_r    = 3**18+2*iseed
    k1_r   = 1313131*iseed
    k1_r   = k1_r-(k1_r/mod2)*mod2

    do i_r = 1, len1
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir1(i_r) = k_r+k1_r*8193
    enddo

    do i_r = 1, len2
      k_r  = k_r *mult
      k1_r = k1_r*mul2
      k1_r = k1_r-(k1_r/mod2)*mod2
      ir2(i_r) = k_r+k1_r*4099
    enddo

    do i_r = 1, len1
      inxt1(i_r) = i_r+1
    enddo
    inxt1(len1) = 1
    ipnt1 = 1
    ipnf1 = ifd1+1

    do i_r = 1, len2
      inxt2(i_r) = i_r+1
    enddo
    inxt2(len2) = 1
    ipnt2 = 1
    ipnf2 = ifd2 + 1
    return
  END SUBROUTINE set_RNG 
  !===================================================================

 !===============Calculate next random number =======================
  !! THIS IS ALMOST PROJECT-INDEPENDENT 
  double precision function rn()
  !integer function rn()
    implicit none
    integer   :: i_r, l_r, k_r
    nrannr = nrannr +1
    if(nrannr>=mxrn) then
      nrannr = 1
      do i_r= 1, mxrn
        l_r = ieor(ir1(ipnt1),ir1(ipnf1))
        k_r = ieor(ir2(ipnt2),ir2(ipnf2))
        irn(i_r) = ieor(k_r,l_r)
        ir1(ipnt1)=l_r
        ipnt1 = inxt1(ipnt1)
        ipnf1 = inxt1(ipnf1)
        ir2(ipnt2) = k_r
        ipnt2 = inxt2(ipnt2)
        ipnf2 = inxt2(ipnf2)
      enddo
    endif 
    !rn = irn(nrannr)
    rn = irn(nrannr)*tm32+0.5d0
  end function rn
  !===================================================================


  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_time_elapse
    implicit none
    !-- read and calculate time (in seconds) -------------------------
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
    h_curr = tval(5)
    t_prev = t_curr
    h_prev = h_curr
    return
  END SUBROUTINE set_time_elapse
    

  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE time_elapse
    implicit none
    
    !-- read and calculate time (in seconds) -------------------------
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
    h_curr = tval(5)

    t_elap = t_curr-t_prev
    if(h_curr<h_prev) t_elap = t_elap+24*3600.d0
    t_prev = t_curr
    h_prev = h_curr 
    return
  END SUBROUTINE time_elapse
  !===================================================================

END PROGRAM

function spec_heat(a,N)
    implicit none
    integer :: N
    real(8) :: spec_heat,a(N)
    spec_heat=a(2)-a(1)*a(1)
    return
end function

function binder(a,N)
    implicit none
    integer :: N
    real(8) :: binder,a(N)
    binder=a(2)/a(1)/a(1)
    return
end function

function suscept(a,N)
    implicit none
    integer :: N
    real(8) :: suscept,a(N)
    suscept=a(2)
    return
end function

