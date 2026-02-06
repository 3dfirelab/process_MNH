module lambda2

    implicit none 
    integer, parameter :: dp = kind(1.d0)  
    Real(kind=dp),dimension(:,:,:),allocatable :: values

    contains

    subroutine clear()
    deallocate(values)
    end subroutine clear 
    
    subroutine compute(nyi,nzi,nxi,xi,yi,zi,u,v,w)
!
!     F. Laporte, Fev. 1999 / CERFACS
!______________________________________________________________________
!
!     Variables Declaration
!______________________________________________________________________
!
      implicit none
!
!***  In
      Integer(kind=dp),intent(in)    :: nyi,nzi,nxi
!***  velocity field with boundaries
      Real(kind=dp),dimension(nxi,nyi,nzi),intent(in) :: u,v,w    
      Real(kind=dp),dimension(nxi,nyi,nzi),intent(in)     :: xi
      Real(kind=dp),dimension(nxi,nyi,nzi),intent(in)     :: zi
      Real(kind=dp),dimension(nxi,nyi,nzi),intent(in)     :: yi
!***  Out
      !Real,dimension(nxi,nyi,nzi),intent(out) :: lambdaM

      ! local variable
      integer*4 :: n1max,n2max,n3max,&
             numax
      real*4 big,small,half
      logical nsmb,fast
!
!      parameter (n1max=132,n2max=117,n3max=560,numax=5,
!     *           fast=.false.,nsmb=.flase.)
      parameter(numax=5)
      parameter(big=1.e+12,small=1.e-30,half=0.5)      
!     
      integer*4 :: i,j,k,ii,jj,kk, &
           nblock,nx,ny,nz,nu,   &
           umesh,ures,uout,ucent,&
           INFO,N                  
!
     Real*4  :: dummy
     real*4,dimension(:,:,:),allocatable  :: x, &
          y, &
          z
     !real*4,dimension(:),allocatable  :: xv1,xv2,&
     !     yv1,yv2
     real*4,dimension(:,:,:),allocatable :: axialvort,&
          vlambda2 
     real*4,dimension(:,:,:,:) ,allocatable :: worka
     real*4 :: mach,alpha,re,time,&
          irho,acoeff,tmpv,tmpx,tmpy, &
          aminl2,rminl2
      
!
      real,dimension(3,3)  ::  gradient,matrice,A
      real,dimension(3)    :: E
      real,dimension(6)    :: WORK(6)
      real*4                 :: dx,dy,dz
!
      real     :: MACHEP
      logical  :: FIRST

!**** initi ronan
      FIRST = .TRUE.

!
!
! allocate array
     allocate(x(nyi,nzi,nxi))
     allocate(y(nyi,nzi,nxi))
     allocate(z(nyi,nzi,nxi))
     !allocate(xv1(nxi))
     !allocate(xv2(nxi))
     !allocate(yv1(nxi))
     !allocate(yv2(nxi))
     allocate(axialvort(nyi,nzi,nxi))
     allocate(vlambda2(nyi,nzi,nxi))
     allocate(worka(nyi,nzi,nxi,numax))
    
     !out 
     allocate(values(nxi,nyi,nzi))

!     Define some constants
      N = 3
!
      n1max = nyi
      n2max = nzi
      n3max = nxi    
      fast = .false.
!______________________________________________________________________
!
!     Program beginning
!______________________________________________________________________
!
      umesh = 17
      ures  = 18
      uout  = 19
      ucent = 20
!     
!      write(*,*) ' '
!      write(*,*) '---------------- LAMBDA2 ------------------'
!      write(*,*) '-------------------------------------------'
!      write(*,*) '---------- extraction de Lambda2 ----------'
!      write(*,*) '-------- a partir des champs MesoNH -------'
!      write(*,*) '-------------------------------------------'
!
!
!______________________________________________________________________
!
!     Read Input - we turn the axis from mesoNH to NTMIX grid (x to z)
!                  then we can use the subroutine of Laporte 
!______________________________________________________________________
!
      
!***  mesh
      do i=1,nyi
          do j=1,nzi
              do k=1,nxi
                 x(i,j,k)=yi(k,i,j)
                 y(i,j,k)=zi(k,i,j)
                 z(i,j,k)=xi(k,i,j)
              enddo
          enddo
      enddo
      
!***  we just use nu=2,3,4 of worka for u,v,w
!***  nevertheless worka is initiate to zero                                
      do nu=1,5
         do k=1,nxi
            do j=1,nzi
               do i=1,nyi
                  worka(i,j,k,nu) = 0.
               enddo
            enddo
         enddo
      enddo
!
      do k=1,nxi
         do j=1,nzi
            do i=1,nyi
               worka(i,j,k,2) = v(k,i,j)
            enddo
         enddo
      enddo
      do k=1,nxi
         do j=1,nzi
            do i=1,nyi
               worka(i,j,k,3) = w(k,i,j)
            enddo
         enddo
      enddo      
      do k=1,nxi
         do j=1,nzi
            do i=1,nyi
               worka(i,j,k,4) = u(k,i,j)
            enddo
         enddo
      enddo


!
!______________________________________________________________________
!
!     Calculate Lambda2
!______________________________________________________________________
!
!      
      do k=2,nxi-1
         do j=2,nzi-1
            do i=2,nyi-1
               
               dx = x(i+1,j+1,k+1) - x(i-1,j-1,k-1)
               dy = y(i+1,j+1,k+1) - y(i-1,j-1,k-1)
               dz = z(i+1,j+1,k+1) - z(i-1,j-1,k-1)
               
               gradient(1,1) = & 
                   (worka(i+1,j,k,2)-worka(i-1,j,k,2))/(dx)
               gradient(1,2) = & 
                   (worka(i,j+1,k,2)-worka(i,j-1,k,2))/(dy)
               gradient(1,3) =  &
                   (worka(i,j,k+1,2)-worka(i,j,k-1,2))/(dz)
               gradient(2,1) = &
                   (worka(i+1,j,k,3)-worka(i-1,j,k,3))/(dx)
               gradient(2,2) =  &
                    (worka(i,j+1,k,3)-worka(i,j-1,k,3))/(dy)
               gradient(2,3) = &
                    (worka(i,j,k+1,3)-worka(i,j,k-1,3))/(dz)
               gradient(3,1) = &
                    (worka(i+1,j,k,4)-worka(i-1,j,k,4))/(dx)
               gradient(3,2) = &
                    (worka(i,j+1,k,4)-worka(i,j-1,k,4))/(dy)
               gradient(3,3) = &
                    (worka(i,j,k+1,4)-worka(i,j,k-1,4))/(dz)
!     

               do ii=1,3
                  do jj=ii,3
                     matrice(ii,jj)=0
                     do kk=1,3
                        matrice(ii,jj) =  &
                            matrice(ii,jj) &
                            +(0.5)*(         &
                                   gradient(ii,kk)*gradient(kk,jj)     &
                                 + gradient(kk,ii)*gradient(jj,kk)     &
                                   )
                     end do
                  end do
               end do
!
               do jj=1,3
                  do ii=1,3
                     A(ii,jj) = matrice(ii,jj)
                  end do 
               end do

!     
               E(1) = A(1,1)
               INFO = 0
!     
               do  jj=1,N
                  do  ii=1,jj
                     A(jj,ii) = A(ii,jj)
                  enddo
               enddo
!     
               CALL TRED1(N,A,E,WORK(1:N),WORK(N+1:2*N))
               CALL TQLRAT(N,E,WORK(N+1:2*N),INFO,FIRST,MACHEP)
               vlambda2(i,j,k)=E(2)
!     
            enddo
!     
         enddo
      enddo
!
!______________________________________________________________________
!
!     compute lambda2 in MesoNH grid
!______________________________________________________________________
      
      do k=1,nxi
         do j=1,nzi
            do i=1,nyi
               values(k,i,j) = vlambda2(i,j,k)
            enddo
         enddo
      enddo

!     
!______________________________________________________________________
!
!     End
!______________________________________________________________________
!
!      write(*,*) ' '
!      write(*,*) 'Tout est bien... (all brights)'
!      write(*,*) 'and thank you for all the door'
!      write(*,*) '-------------------------------------------'
!      write(*,*) ' '


      end subroutine compute


!#################
      !real FUNCTION PYTHAG (A, B)
      SUBROUTINE PYTHAG (A, B, res)

!***PURPOSE  Compute the complex square root of a complex number without
!            destructive overflow or underflow.

      real(kind=dp),intent(in)  :: A,B
      real(kind=dp),intent(out)  :: res
      
      real             :: P,Q,R,S,T

      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
      10 CONTINUE
         if(p.gt.1.e10) then
         r=0.
         else
         R = (Q/P)**2
         end if
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
      20 res = P
    END SUBROUTINE PYTHAG

!   20 PYTHAG = P
      !RETURN
      !END FUNCTION PYTHAG
!#################
      
     ! real FUNCTION R1MACH (I)
    SUBROUTINE R1MACH(I,res)
!***PURPOSE  Return floating point machine dependent constants.
      integer,intent(in)  :: I
      real(kind=dp),intent(out)  :: res


      INTEGER   :: SMALL(2)
      INTEGER   :: LARGE(2)
      INTEGER   :: RIGHT(2)
      INTEGER   :: DIVER(2)
      INTEGER   :: LOG10(2)

      real   :: RMACH(5)
      SAVE RMACH

      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))

!     MACHINE CONSTANTS FOR THE SUN
!
      DATA RMACH(1) / Z'00800000' /
      DATA RMACH(2) / Z'7F7FFFFF' /
      DATA RMACH(3) / Z'33800000' /
      DATA RMACH(4) / Z'34000000' /
      DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /

!      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'R1MACH',
!     +   'I OUT OF BOUNDS', 1, 2)

      res = RMACH(I)

     ! END FUNCTION R1MACH
    END SUBROUTINE R1MACH

!#################
      SUBROUTINE TQLRAT (N, D, E2, IERR,FIRST,MACHEP)

!***PURPOSE  Compute the eigenvalues of symmetric tridiagonal matrix
!            using a rational variant of the QL method.

        implicit none

        Integer,intent(in) :: N
        Real(kind=dp),dimension(N),intent(inout) :: D,E2
        Integer,intent(inout):: IERR
        Logical,intent(inout)    :: FIRST
        real(kind=dp),intent(inout)       :: MACHEP

        Integer*4    :: I,J,L,M,II,L1,MML
        Real*4       :: B,C,F,G,H,S !,MACHEP
        Real :: P,R
        !Real*4       :: PYTHAG

        !Real       :: R1MACH

        
      IF (FIRST) THEN
         !MACHEP = R1MACH(4)
         CALL R1MACH(4, MACHEP )
      ENDIF
      FIRST = .FALSE.
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
!
      F = 0.0E0
      B = 0.0E0
      E2(N) = 0.0E0
!
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         IF (B .GT. H) GO TO 105
        B = H
         C = B * B
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * S)
         CALL PYTHAG(P,1.0D0,R)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
!
         DO 140 I = L1, N
  140    D(I) = D(I) - H
!
         F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0E0) G = B
         H = G
         S = 0.0E0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0E0) G = B
            H = G * P / R
  200    CONTINUE
!
         E2(L) = S * G
         D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0E0) GO TO 210

         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0E0) GO TO 130
  210    P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE

  250    I = 1
  270    D(I) = P
  290 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001  CONTINUE

      END SUBROUTINE TQLRAT 

!#################
      SUBROUTINE TRED1 (N, A, D, E, E2)

        implicit none

!      INTEGER I,J,K,L,N,II,NM,JP1
!      real*4 A(3,3),D(3),E(3),E2(3)
!      real*4 F,G,H,SCALE

        Integer,intent(in)     ::  N 
        Real(kind=dp),dimension(N,N),intent(inout)    :: A
        Real(kind=dp),dimension(N),intent(inout)      :: D,E,E2

        Integer*4                :: I,J,K,L,II,NM,JP1
        Real*4                   :: F,G,H,SCALE

        NM=N

!
!***FIRST EXECUTABLE STATEMENT  TRED1

      DO 100 I = 1, N
  100 D(I) = A(I,I)

!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
         
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K))
!
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 290
!
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
!
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0E0
!
         DO 240 J = 1, L
            G = 0.0E0
!     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
!
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
!     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
!
         H = F / (H + H)
!     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
!
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
!
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
!
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
!
      
      END SUBROUTINE TRED1

    

end module lambda2

