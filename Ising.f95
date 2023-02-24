
Program Ising

  external CalcE
  external EnergyNN
  external mag_av
  external sweep

  integer NMax,NN                  !/* Maximal lattice size, # neighbors*/
  real*8 JJ, HH                      !/* Coupling constant, magnetic field */
  parameter (NMax=50,NN=4)
  parameter (JJ=1.d0,HH=0.d0)


  integer mag(Nmax,Nmax)            !/* 2D Ising Lattice */
  integer i, j, k, m, n             !/* Loop counters */
  integer s,d                       !/* Lattice spin variables */  
  real*8 Energy                     !/* Total lattice energy */
  integer Inn(NN)                    !/* Nearest neighbor array I */
  parameter (Inn=(/1,-1,0,0/))
  integer Jnn(NN)                       !/* Nearest neighbor array J */
  parameter (Jnn=(/0,0,1,-1/)) 
  integer Inew, Jnew                !/* Nearest neighbot indices */ 
  real*8 Etemp, deltaE, E             !/* Temp energy variables for MC moves */ 
  integer :: accept=0                   !/* Number of accepted moves */   
  integer :: move=0                      !/* Number of moves total */   
  real*8 mag_sum, mag_sum2, mag_sum4, mag_av1, binder

  real*8 T                          !/* temperature (in units of J/k_B) */
  real*8 beta
  real*8 r
  integer sweeps                    !/* number of measurement sweeps */
  integer warm                      !/* number of warm-up sweeps */
  integer L                             !/* lattice dimension */
  integer :: seed
  CALL SYSTEM_CLOCK(seed)


 ! /***************************
 !  * Initialization          *
 !  ***************************/

!  write (*,*) "Temperatur?"
!  read(*,*) T
  T=5.
  beta=1.d0/T

!  write (*,*) "Gittergroesse (LxL) ?"
!  read(*,*) L
  L=5

!  write (*,*) "# sweeps?"
!  read(*,*) sweeps
  sweeps=500000
  
!  write (*,*) "# warm up  sweeps?"
!  read(*,*) warm
  warm=300000		


!  Energy = CalcE(mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
!  call  outputmag(mag,L,Nmax)
!  write (*,*) 'Average energy',Energy

 
! /***************************
! /* warum up sweeps */
! /***************************
 open(unit=1, file='Magnetisierung5.dat')
 do n=1,60,1
  T=n*0.1
  beta=1./T
  call srand(seed)

!  call random_init(mag,L,Nmax)  
  call alternate_init(mag,L,Nmax)
  call  outputmag(mag,L,Nmax)
  do m=1,warm,1
	
    call sweep(mag,L,Jnn,Inn,NN,JJ,HH,Nmax,beta,accept)
  
  end do

  Energy = CalcE(mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
 

! /***************************
! /* end warum up sweeps */ 
! /***************************
 

! /***************************
! /*  sweeps */
! /***************************

  accept=0
  move=0
  mag_sum=0.
  mag_sum2=0.
  mag_sum4=0.
  binder=0.

  do m=1,sweeps,1
	
    call sweep(mag,L,Jnn,Inn,NN,JJ,HH,Nmax,beta,accept)
    call mag_av(mag,L,Nmax,mag_av1)
    mag_sum=mag_sum+mag_av1
    mag_sum2=mag_sum2+mag_av1**2
    mag_sum4=mag_sum4+mag_av1**4
!    write(*,*) mag_av1
    move=move+1
  end do
  write(*,*) 'mittlere Magnetisierung=',mag_sum/sweeps
  
  write(*,*)accept, 'von', move, 'gewechselt'
  Energy = CalcE(mag,L,Jnn,Inn,NN,JJ,HH,Nmax)


! /***************************
! /*  end sweeps */ 
! /***************************
  mag_sum=mag_sum/sweeps
  mag_sum2=mag_sum2/sweeps
  mag_sum4=mag_sum4/sweeps

  binder=1.-(mag_sum4)/(3.*mag_sum2**2)
  call  outputmag(mag,L,Nmax)
  write (*,*) 'Average energy',Energy
  write (1,*) T, mag_sum, mag_sum2, mag_sum4, binder
 end do
 close(1)
end Program Ising



subroutine random_init(mag,L,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax), L
 
  integer i,j
  
  do i=1,L
     do j=1,L
       	if (rand().gt.0.5d0) then
	   mag(i,j) = 1
        else
           mag(i,j) = -1
        endif
     enddo
  enddo

  return
end subroutine random_init

subroutine alternate_init(mag,L,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax), L
 
  integer i,j
!Anfang abwechselnd Spin up und down
  do i=1,L,1
    do j=1,L,1
      if (mod(i,2)==0) then
        if ((mod(j,2)==0)) then
           mag(i,j)=1
        else
	   mag(i,j)=-1
        end if
      else
        if ((mod(j,2)==0)) then
           mag(i,j)=-1
        else
	   mag(i,j)=1
        end if
       end if
     end do
  end do
  return
end subroutine alternate_init

subroutine outputmag(mag,L,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax), L
 
  integer i,j
  write (*,*) 'Spin snapshot'

  do i=1,L
     do j=1,L
        if(mag(i,j).eq. 1) then
	   write(*,'(a)',ADVANCE='NO') '-> '
        else
           write(*,'(a)',ADVANCE='NO') '<- '
        endif
      
     enddo
     !     Newline to complete
     write (*,*) 
     write (*,*)
  enddo
  return 
end subroutine outputmag

subroutine sweep(mag,L,Jnn,Inn,NN,JJ,HH,Nmax,beta,accept)
  integer Nmax
  integer mag(Nmax,Nmax),Jnn(NN),Inn(NN),L,NN, accept
  
  real*8 JJ, HH, E, Etemp, deltaE, r, beta
  integer i,j,k,Inew,Jnew
  integer seed
  

    i=int(rand()*L+1)			!Zufallszahl zwischen (0 und L) +1 damit alle gleichwahrscheinlich sind, danach auf Integer konvertieren
    j=int(rand()*L+1)
  
    E=EnergyNN(i,j,mag,L,Jnn,Inn,NN,JJ,HH,Nmax)

    mag(i,j)=mag(i,j)*(-1)		!Spin flippen

    Etemp=EnergyNN(i,j,mag,L,Jnn,Inn,NN,JJ,HH,Nmax)

    deltaE=Etemp-E
    r=exp(-deltaE*beta)

    if(r>1) then			!min(r,1)
      r=1.
    end if

    if(rand() < r) then
      mag(i,j)=mag(i,j)
      accept=accept+1
    else
      mag(i,j)=mag(i,j)*(-1)
    end if

  return 
end subroutine sweep


real function CalcE(mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax),Jnn(NN),Inn(NN),L,NN
  
  real*8 JJ, HH
  integer i,j,k,Inew,Jnew
  real*8 Energy

  !/* Determine the initial energy */
  Energy = 0.0
  
  do i=1,L
     do j=1,L
        
	!/* Loop over nearest neighbors */
	do k=1, NN  
            Inew = i + Inn(k)       
            Jnew = j + Jnn(k)
    
	    !/* Check periodic boundary conditions */
	    if (Inew .le. 0) then
	      Inew = L
            else 
               if(Inew .gt. L) then 
                  Inew = 1
               endif
            endif    
            if (Jnew .le. 0) then
	      Jnew = L
            else 
               if(Jnew .gt. L) then 
                  Jnew = 1
               endif
            endif
	    
	    !/* Update the energy */
	    Energy = Energy-JJ * mag(i,j) * mag(Inew,Jnew)
        enddo
	!/*Calculate the contribution from the field H */
	Energy = Energy - 2.d0*HH*mag(i,j)
     enddo
   enddo

   !/* Account for double counting */
   Energy = Energy/2.d0

   CalcE=Energy
   return
 end function CalcE

! Energie der nächsten Nachbarn
real function EnergyNN(i,j,mag,L,Jnn,Inn,NN,JJ,HH,Nmax)
  integer Nmax
  integer mag(Nmax,Nmax),Jnn(NN),Inn(NN),L,NN
  
  real*8 JJ, HH
  integer i,j,k,Inew,Jnew
  real*8 Energy

  !/* Determine the initial energy */
  Energy = 0.0
  
	!/* Loop over nearest neighbors */
	do k=1, NN  
            Inew = i + Inn(k)       
            Jnew = j + Jnn(k)
    
	    !/* Check periodic boundary conditions */
	    if (Inew .le. 0) then
	      Inew = L
            else 
               if(Inew .gt. L) then 
                  Inew = 1
               endif
            endif    
            if (Jnew .le. 0) then
	      Jnew = L
            else 
               if(Jnew .gt. L) then 
                  Jnew = 1
               endif
            endif
	    
	    !/* Update the energy */
	    Energy = Energy-JJ * mag(i,j) * mag(Inew,Jnew)
        enddo
	!/*Calculate the contribution from the field H */
	Energy = Energy - HH*mag(i,j)

   EnergyNN=Energy
   return
 end function EnergyNN

subroutine mag_av(mag,L,Nmax, mag_av1)
!/* computes average magnetization */

  integer i,j,Nmax,L
  real*8 x
  integer mag(Nmax,Nmax)
  real*8 mag_av1
  x=0.0

  do i=1,L,1
    do j=1,L,1
      x=x+1.0*mag(i,j)

    end do
  end do
  x=x/L**2

  mag_av1=abs(x)
return 
end subroutine mag_av
