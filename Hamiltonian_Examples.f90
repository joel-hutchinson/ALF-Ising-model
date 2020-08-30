      Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables

      
      Implicit none

!>    Public variables. Have to be set by user 
      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot
!>    Variables for updating scheme
      Logical              :: Propose_S0
      

!>    Private variables 
      Type (Lattice),       private :: Latt 
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Logical,              private :: One_D_Ising = .false.  !!!!!!!!!!!!!! New to 1D Ising model !!!!!!!!!
      Logical,              private :: PiFlux = .false.  !!!!!!!!!!!!!! Added Jan 9, 2019 !!!!!!!!!
      Integer,              private :: N_coord, Norb
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell



!>    Private Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)),        private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer, allocatable, private :: L_bond(:), L_bond_inv(:), Ising_nnlist(:,:)

      
    contains 


      Subroutine Ham_Set

          Implicit none

#ifdef MPI
          include 'mpif.h'
#endif   

          integer :: ierr


          !!  Namelist for  various models. 
          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, Dtau, Beta

          NAMELIST /VAR_Hubbard_Ising/  ham_T, ham_chem, ham_U, Dtau, Beta,  Ham_h, Ham_J, Ham_xi

          
#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          

#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
 
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          Propose_S0 = .false.
          !!!!!!!!!! New to 1D !!!!!!!!!
          If ( model == "Hubbard_SU2_Ising_1D") then
             One_D_Ising = .true.
             If (L2 .NE. 1) then
                Write(6,*) "1D Ising requires L2=1"
                Stop
             Endif
             model = "Hubbard_SU2_Ising"
          endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!Added Jan 9, 2019!!!!!!!!!!!!!
          If ( model == "Hubbard_SU2_PiFlux") then
            PiFlux = .true.
            model = "Hubbard_SU2"
          endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ADJUST WHEN ISING SPINS ARE ADDED

          If ( Model == "Hubbard_Mz") then
             N_FL = 2
             N_SUN = 1
          elseif  ( Model == "Hubbard_SU2" ) then
             N_FL = 1
             N_SUN = 2
          elseif  ( Model == "Hubbard_SU2_Ising" ) then
             N_FL = 1
             N_SUN = 2
             If ( Lattice_type == "Honeycomb" ) then 
                Write(6,*) "Hubbard_SU2_Ising is only implemented for a square lattice"
                Stop
             Endif
          else
             Write(6,*) "Model not yet implemented!"
             Stop
          endif
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             If (Model == "Hubbard_SU2" .or. Model == "Hubbard_Mz" )  Read(5,NML=VAR_Hubbard)
             If (Model == "Hubbard_SU2_Ising"                      )  Read(5,NML=VAR_Hubbard_Ising)
             CLOSE(5)
             
#ifdef MPI
          endif
          If (Model == "Hubbard_SU2" .or. Model == "Hubbard_Mz" ) then 
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          elseif ( Model == "Hubbard_SU2_Ising" ) then
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          Endif
#endif
          Call Ham_hop
          Ltrot = nint(beta/dtau)

          If  ( Model == "Hubbard_SU2_Ising" )  Call Setup_Ising_action
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             !!!!!!!!!! New to 1D !!!!!!!!!!
             If (One_D_Ising ) then
                Write(50,*) 'Model is       : ', 'Hubbard_SU2_Ising_1D'
             Else
                Write(50,*) 'Model is      : ', Model
             Endif
             !!!!!!!!!!!!!!!!!!!!!!!!!!
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) 'L1, L2        : ', L1, L2
             Write(50,*) '# of orbitals : ', Ndim
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_chem      : ', Ham_chem
             If (Model == "Hubbard_SU2" .or. Model == "Hubbard_Mz" ) then 
                Write(50,*) 'Ham_U         : ', Ham_U
             elseif  ( Model == "Hubbard_SU2_Ising" ) then
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_xi        : ', Ham_xi
                Write(50,*) 'Ham_J         : ', Ham_J
                Write(50,*) 'Ham_h         : ', Ham_h
             Endif
#if defined(STAB1) 
             Write(50,*) 'STAB1 is defined '
#endif
#if defined(QRREF) 
             Write(50,*) 'QRREF is defined '
#endif
             close(50)
#ifdef MPI
          endif
#endif
          call Ham_V

        end Subroutine Ham_Set
!=============================================================================
        Subroutine Ham_Latt
          Implicit none

!!!!!!!   In this routine one sets up the lattice.           
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no

          If ( Lattice_type =="Square" ) then
             Norb = 1
             N_coord   = 2
             One_dimensional = .false.
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             If ( L1 == 1 .or. L2 == 1 ) then 
                One_dimensional = .true.
                N_coord   = 1
                If (L1 == 1 ) then 
                   Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                   Stop
                endif
             endif
          elseif (Lattice_type=="Honeycomb" ) then
             Norb    = 2
             N_coord = 3
             a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
             a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0

             !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             L1_p    =  dble(L1) * a1_p
             L2_p    =  dble(L2) * a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif

          ! This is for the orbital structure.
          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                ! For the Honeycomb lattice no = 1,2 corresponds to the A,and B sublattices.
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !Setup the hopping
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, J1, I2, n, Ncheck,nc, nc1, no, Iy

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If ( Lattice_type =="Square" ) then
                   If (One_dimensional ) then 
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I, 1, 0)
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   else
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I,1,0)
                         I2 = Latt%nnlist(I,0,1)
                         If (PiFlux) then !!!!! Added Jan 9,2019 !!!!!!!!
                            Iy = Latt%list(I,2)
                            !Write(6,*) I, " ", Latt%list(I,1)
                            Op_T(nc,n)%O(I,I1) = cmplx(-((-1.0d0)**Iy)*Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I1,I) = cmplx(-((-1.0d0)**Iy)*Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                         else
                            Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                        endif
                      ENDDO
                   endif
                elseif ( Lattice_type=="Honeycomb" ) then
                   DO I = 1, Latt%N
                      do no = 1,Norb
                         I1 = Invlist(I,no)
                         Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      enddo
                      I1 = Invlist(I,1)
                      J1 = I1
                      Do nc1 = 1,N_coord
                         select case (nc1)
                         case (1)
                            J1 = invlist(I,2) 
                         case (2)
                            J1 = invlist(Latt%nnlist(I,1,-1),2) 
                         case (3)
                            J1 = invlist(Latt%nnlist(I,0,-1),2) 
                         case default
                            Write(6,*) ' Error in  Ham_Hop '  
                            Stop
                         end select
                         Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Enddo
                   Enddo
                endif

                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6 .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g = 0.d0
                else
                   Op_T(nc,n)%g = -Dtau
                endif
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n)) 
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 

          !Setup the interaction

          
          Integer :: nf, I, I1, I2,  nc, nc1,  J
          Real (Kind=Kind(0.d0)) :: X
          
          If (Model == "Hubbard_SU2")  then
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Ndim,N_FL))
             do nf = 1,N_FL
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                Do i = 1,Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Elseif (Model == "Hubbard_Mz")  then
             Allocate(Op_V(Ndim,N_FL))
             do nf = 1,N_FL
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                X = 1.d0
                if (nf == 2) X = -1.d0
                Do i = 1,Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Elseif  ( Model == "Hubbard_SU2_Ising" ) then
             If ( Abs(Ham_U) >  1.D-6) then
                Allocate(Op_V(3*Ndim,N_FL))
                do nf = 1,N_FL
                   do i  =  1, N_coord*Ndim
                      call Op_make(Op_V(i,nf),2)
                   enddo
                   do i  = N_coord*Ndim +1 ,  N_coord*Ndim + Ndim ! For Hubbatd
                      Call Op_make(Op_V(i,nf),1)
                   enddo
                enddo
             else
                Allocate(Op_V(N_coord*Ndim,N_FL)) !!!!!!!!! Changed for 1D (2->N_coord) !!!!!!!!!!!
                do nf = 1,N_FL
                   do i  =  1, N_coord*Ndim
                      call Op_make(Op_V(i,nf),2)
                   enddo
                enddo
             endif
             Do nc = 1,Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                !!!!!!!!!!!!!!!New to 1D!!!!!!!!!!!!!!
                I1 = L_bond_inv(nc)
                I2 = I1
                I2 = Latt%nnlist(I1,1,0)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Op_V(nc,1)%P(1) = I1
                Op_V(nc,1)%P(2) = I2
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%g = cmplx(-dtau*Ham_xi,0.D0,kind(0.D0))
                Op_V(nc,1)%alpha = cmplx(0d0,0.d0, kind(0.D0)) 
                Op_V(nc,1)%type =1
                Call Op_set( Op_V(nc,1) )
             Enddo

             If ( Abs(Ham_U) >  1.D-6) then
                Do i = 1,Ndim
                   nc1 = N_coord*Ndim + i
                   Op_V(nc1,1)%P(1)   = i
                   Op_V(nc1,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc1,1)%g      = sqrt(cmplx(-dtau*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                   Op_V(nc1,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.d0))
                   Op_V(nc1,1)%type   = 2
                   Call Op_set( Op_V(nc1,1) )
                Enddo
             Endif
             
          Endif
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)
          Implicit none
          Integer, Intent(IN) :: n,nt
          Integer :: nt1,I
          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then
             !!!!!!!!!!!! New to 1D !!!!!!!!!!!!!
             If (One_D_Ising) then
                do i = 1,2
                   S0 = S0*DW_Ising_space(nsigma(n,nt)*nsigma(Ising_nnlist(n,i),nt))
                enddo
             Else
                do i = 1,4
                   S0 = S0*DW_Ising_space(nsigma(n,nt)*nsigma(Ising_nnlist(n,i),nt))
                enddo
             Endif
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             nt1 = nt +1
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             nt1 = nt - 1
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif
        end function S0

!===================================================================================
        Subroutine Setup_Ising_action

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
        !!!!!!!!!!!!!!!!!!New to 1D !!!!!!!!!!!!!!!!!!!!!!
          Allocate (L_Bond(Latt%N),  L_bond_inv(Latt%N*N_coord) )
          nc = 0
          do nth = 1,2*N_coord
             Do n1= 1, L1/2
                Do n2 = 1,L2
                   nc = nc + 1
                   I1 = 1
                   If (nth == 1 ) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                   elseif (nth == 2) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p  + latt%a1_p
                      I1 = Inv_R(X_p,Latt)
                   endif
                   L_bond(I1) = nc
                   L_bond_inv(nc) = I1
                   ! The bond is given by  I1, I1 + a_(n_orientation).
                Enddo
             Enddo
          Enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Setup the nearest neigbour lists for the Ising spins.
          !!!!!!!!New to 1D!!!!!!!!!!!!!
          If (One_D_Ising) then
            allocate(Ising_nnlist(Latt%N,2))
            do I  = 1,Latt%N
                n  = L_bond(I)
                n1 = L_bond(Latt%nnlist(I, 1, 0))
                n2 = L_bond(Latt%nnlist(I, -1, 0))
                Ising_nnlist(n,1) = n1
                Ising_nnlist(n,2) = n2
            enddo
          Else
            Print*, "Cannot do 2D Ising with code as is"
           ! allocate(Ising_nnlist(2*Latt%N,4))
            !do I  = 1,Latt%N
             !   n  = L_bond(I,1)
              !  n1 = L_bond(Latt%nnlist(I, 1, 0),2)
               ! n2 = L_bond(Latt%nnlist(I, 0, 0),2)
                !n3 = L_bond(Latt%nnlist(I, 0,-1),2)
                !n4 = L_bond(Latt%nnlist(I, 1,-1),2)
                !Ising_nnlist(n,1) = n1
                !Ising_nnlist(n,2) = n2
                !Ising_nnlist(n,3) = n3
                !Ising_nnlist(n,4) = n4
                !n  = L_bond(I,2)
                !n1 = L_bond(Latt%nnlist(I, 0, 1),1)
                !n2 = L_bond(Latt%nnlist(I,-1, 1),1)
                !n3 = L_bond(Latt%nnlist(I,-1, 0),1)
                !n4 = L_bond(Latt%nnlist(I, 0, 0),1)
                !Ising_nnlist(n,1) = n1
                !Ising_nnlist(n,2) = n2
                !Ising_nnlist(n,3) = n3
                !Ising_nnlist(n,4) = n4
            !enddo
          Endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!! New to 1D !!!!!!!!!
          If(Abs(Ham_h)<0.001) then
            DW_Ising_tau  ( 1) = 1.d0  !tanh(Dtau*Ham_h)
            DW_Ising_tau  (-1) = 1.d0  !1.D0/DW_Ising_tau(1)
            DW_Ising_Space( 1) = exp(-2.d0*beta*Ham_J)
            DW_Ising_Space(-1) = exp( 2.d0*beta*Ham_J)
          Else
            DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
            DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
            DW_Ising_Space( 1) = exp(-2.d0*dtau*Ham_J)
            DW_Ising_Space(-1) = exp( 2.d0*dtau*Ham_J)
          Endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        End Subroutine Setup_Ising_action
!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename

          ! Scalar observables
          If  ( model == "Hubbard_SU2_Ising" ) Then
             Allocate ( Obs_scal(8) )       !!!!!!!!!!!!!!!New to 1D (7->8)!!!!!!!!!!!
             Do I = 1,Size(Obs_scal,1)
                select case (I)
                case (1)
                   N = 1;   Filename ="Kin"
                case (2)
                   N = 1;   Filename ="Pot"
                case (3)
                   N = 1;   Filename ="Part"
                case (4)
                   N = 1;   Filename ="Ener"
                case (5)
                   N = 1;   Filename ="X"
                case (6)
                   N = 1;   Filename ="M2"
                case (7)
                   N = 2;   Filename ="M4"
                !!!!!!!!!!!!New to 1D!!!!!!!!!!!!!!!
                case (8)
                   N = 1;   Filename ="Ener_Ising"
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Call Obser_Vec_make(Obs_scal(I),N,Filename)
             enddo
          Else
             Allocate ( Obs_scal(4) )
             Do I = 1,Size(Obs_scal,1)
                select case (I)
                case (1)
                   N = 1;   Filename ="Kin"
                case (2)
                   N = 1;   Filename ="Pot"
                case (3)
                   N = 1;   Filename ="Part"
                case (4)
                   N = 1;   Filename ="Ener"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Call Obser_Vec_make(Obs_scal(I),N,Filename)
             enddo
          endif

          ! Equal time correlators
          If  ( model == "Hubbard_SU2_Ising" ) Then
             Allocate ( Obs_eq(5) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Ns = Latt%N;  No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N;  No = Norb;  Filename ="Den"
                case (5)
                   Ns = Latt%N;  No = N_coord;  Filename ="ZZ"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = 1
                Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
             enddo
          Else
             Allocate ( Obs_eq(4) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Ns = Latt%N;  No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N;  No = Norb;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = 1
                Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
             enddo
          Endif
          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(4) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N; No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N; No = Norb;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
        end Subroutine Alloc_obs

!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J, no, ntau1, N, nt, ns
          Real    (Kind=Kind(0.d0)) :: X_ave, X, X2_ave !!!!!!!New to 1D, added X2_ave!!!!!!!!!

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          

          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables. 
          
          Obs_scal(1)%N         =  Obs_scal(1)%N + 1
          Obs_scal(1)%Ave_sign  =  Obs_scal(1)%Ave_sign + Real(ZS,kind(0.d0))
          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do J = 1,Ndim
                Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
             ENddo
          Enddo
          Zkin = Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          Obs_scal(2)%N         =  Obs_scal(2)%N + 1
          Obs_scal(2)%Ave_sign  =  Obs_scal(2)%Ave_sign + Real(ZS,kind(0.d0))
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
             Do I = 1,Ndim
                ZPot = ZPot + Grc(i,i,1) * Grc(i,i, 1)
             Enddo
             Zpot = Zpot*ham_U
          elseif ( Model == "Hubbard_Mz")  then
             Do I = 1,Ndim
                ZPot = ZPot + Grc(i,i,1) * Grc(i,i,2)
             Enddo
             Zpot = Zpot*ham_U
          Endif
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Obs_scal(3)%N         =  Obs_scal(3)%N + 1
          Obs_scal(3)%Ave_sign  =  Obs_scal(3)%Ave_sign + Real(ZS,kind(0.d0))
          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS


          Obs_scal(4)%N         =  Obs_scal(4)%N + 1
          Obs_scal(4)%Ave_sign  =  Obs_scal(4)%Ave_sign + Real(ZS,kind(0.d0))
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          If ( Model == "Hubbard_SU2_Ising" ) then
             Obs_scal(5)%N         =  Obs_scal(5)%N + 1
             Obs_scal(5)%Ave_sign  =  Obs_scal(5)%Ave_sign + Real(ZS,kind(0.d0))
             X_ave = 0.d0
             ntau1 = ntau + 1
             If (ntau == Ltrot)  ntau1 = 1
            do nt = 1,Ltrot
                Do I = 1,Latt%N
                !!!!!!!!!!!!!!New to 1D!!!!!!!!!!!!!!! This is the spin-flip observable X
                    do no = 1,1
                    X_ave = X_ave + DW_Ising_tau( nsigma(L_bond(I),ntau)*nsigma(L_bond(I),ntau1) )
                    Enddo
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Enddo
            enddo
             Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + cmplx(X_ave,0.d0,kind(0.d0))*ZP*ZS
             if (Ntau == Ltrot/2) then
                Obs_scal(6)%N         =  Obs_scal(6)%N + 1
                Obs_scal(6)%Ave_sign  =  Obs_scal(6)%Ave_sign + Real(ZS,kind(0.d0))
                Obs_scal(7)%N         =  Obs_scal(7)%N + 1
                Obs_scal(7)%Ave_sign  =  Obs_scal(7)%Ave_sign + Real(ZS,kind(0.d0))
                ns = 1
                if (Ham_J < 0) ns = -1
                !!!!!!!!!!!!!New to 1D, this is for h=0!!!!!!!!!!!!!! This is the magnetization
                if (Abs(Ham_h)<0.001) then
                    X_ave=0.d0
                    X2_ave=0.d0
                    do nt = 1,Ltrot
                        N = 0
                        do I  = 1,Latt%N
                            N = N  +  (ns**I)*nsigma(L_bond(I),nt) !* ns*nsigma(L_bond(I,2),nt)
                        Enddo
                        X_ave = X_ave+(Real(N,kind(0.d0))/Real(Latt%N,kind(0.d0)))**2
                        X2_ave = X2_ave+(Real(N,kind(0.d0))/Real(Latt%N,kind(0.d0)))**4
                    Enddo
                    Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + cmplx(X_ave/Real(Ltrot, kind(0.d0)),0.d0,kind(0.d0))*ZP*ZS
                    Obs_scal(7)%Obs_vec(1) = Obs_scal(7)%Obs_vec(1) + cmplx(X2_ave/Real(Ltrot, kind(0.d0)),0.d0,kind(0.d0))*ZP*ZS
                else
                    N = 0
                    do nt = 1,Ltrot
                        do I  = 1,Latt%N
                            N = N  +  (ns**I)*nsigma(L_bond(I),nt) !!!!New to 1d!!!!!!!! This is the magnetization
                        Enddo
                        X_ave = Real(N,kind(0.d0))/(Real( Latt%N * Ltrot, Kind(0.d0) ))
                    Enddo
                    Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + cmplx(X_ave**2,0.d0,kind(0.d0))*ZP*ZS
                    Obs_scal(7)%Obs_vec(1) = Obs_scal(7)%Obs_vec(1) + cmplx(X_ave**4,0.d0,kind(0.d0))*ZP*ZS
                endif
                !!!!!!!!!!!!!!!New to 1D!!!!!!!!!!!!!!!!Ising energy
                X_ave=0.0
                do nt = 1,Ltrot
                    Do n = 1, Latt%N
                        if (ham_h<0.001) then
                            X_ave = X_Ave -ham_J*nsigma(n,nt)*nsigma(Ising_nnlist(n,1),nt)
                        else
                            X_ave = X_Ave - ham_J*nsigma(n,nt)*nsigma(Ising_nnlist(n,1),nt)&
                            & - ham_h*DW_Ising_tau( nsigma(L_bond(I),ntau)*nsigma(L_bond(I),ntau1) )
                        endif
                    enddo
                enddo
                Obs_scal(8)%N         =  Obs_scal(8)%N + 1
                Obs_scal(8)%Ave_sign  =  Obs_scal(8)%Ave_sign + Real(ZS,kind(0.d0))
                Obs_scal(8)%Obs_vec(1)  =    Obs_scal(8)%Obs_vec(1) + X_ave/(Real( Latt%N * Ltrot, Kind(0.d0) ))*ZP*ZS
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             Endif
          endif
          
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          If ( Model == "Hubbard_SU2_Ising" ) then
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   do no_I = 1, N_coord
                      Do no_J = 1,N_coord
                        !!!!!!!!!!!!!!!!New to 1D!!!!!!!!!!!!!!!
                         X = Real(nsigma(L_bond(I),ntau)*nsigma(L_bond(J),ntau), Kind(0.d0))
                         Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J)  +  &
                              &               cmplx(X,0.d0,Kind(0.d0)) * ZP*ZS
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         
                      Enddo
                   Enddo
                Enddo
             Enddo
          endif
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then
             ! SU(N) symmetry is present.
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) *  ZP*ZS 
                   ! SpinZ
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! SpinXY
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                        &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                        &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                        &                                   ) * Z* ZP*ZS
                ENDDO
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
             ENDDO
          elseif (Model == "Hubbard_Mz" ) Then
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &                          ( GRC(I1,J1,1) + GRC(I1,J1,2)) *  ZP*ZS 
                   ! SpinZ
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))    ) * ZP*ZS
                   ! SpinXY
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                     & (   GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)    ) * ZP*ZS
                   !Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &   (GRC(I1,I1,2) + GRC(I1,I1,1))*(GRC(J1,J1,2) + GRC(J1,J1,1))     ) * ZP*ZS
                enddo
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  (GRC(I1,I1,1) + GRC(I1,I1,2)) * ZP*ZS
             enddo
          Endif
                

        end Subroutine Obser
!=====================================================
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  Z * GT0(I1,J1,1) * ZP* ZS
                   
                   ! SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                        &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                        &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                Enddo
                Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
                     &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
             Enddo
          Elseif ( Model == "Hubbard_Mz"  ) then 
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   !Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &   +   ( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS

                   !SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                       & +  ( &
                       &    (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                       &  - (G0T(J1,I1,1) * GT0(I1,J1,1)  +  G0T(J1,I1,2) * GT0(I1,J1,2) )    )*ZP*ZS

                   !SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &  - &
                        &   (G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1))*ZP*ZS
                   !Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  (                                        &  
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                        & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS     

                enddo
             
                Obs_tau(4)%Obs_Latt0(no_I) =  Obs_tau(4)%Obs_Latt0(no_I) + &
                     &       (cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
             Enddo
          Endif
          
        end Subroutine OBSERT

!==========================================================        
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau
          
          !Local 
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call  Print_bin_Vec(Obs_scal(I))
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call  Print_bin_Latt(Obs_eq(I),Latt,dtau)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I),Latt,dtau)
             enddo
          endif

        end Subroutine Pr_obs
!===================================================================================           
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          ! Local 
          Integer :: I

          Do I = 1,Size(Obs_scal,1)
             Call Obser_vec_Init(Obs_scal(I))
          Enddo

          Do I = 1,Size(Obs_eq,1)
             Call Obser_Latt_Init(Obs_eq(I))
          Enddo

          If (Ltau == 1) then
             Do I = 1,Size(Obs_tau,1)
                Call Obser_Latt_Init(Obs_tau(I))
             Enddo
          Endif

        end Subroutine Init_obs

      end Module Hamiltonian
