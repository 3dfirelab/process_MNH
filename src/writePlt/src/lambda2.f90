
      subroutine lambda2(nyi,nzi,nxi,xi,yi,zi,u,v,w,lambdaM)
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
      Integer,intent(in)    :: nyi,nzi,nxi
!***  velocity field with boundaries
      Real,dimension(nxi,nyi,nzi),intent(in) :: u,v,w    
      Real,dimension(nxi,nyi,nzi),intent(in)     :: xi
      Real,dimension(nxi,nyi,nzi),intent(in)     :: zi
      Real,dimension(nxi,nyi,nzi),intent(in)     :: yi
!***  Out
      Real,dimension(nxi,nyi,nzi),intent(out) :: lambdaM

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
      real*4,dimension(3,3)  ::  gradient,matrice,A
      real*4,dimension(3)    :: E
      real*4,dimension(6)    :: WORK(6)
      real*4                 :: dx,dy,dz
!
      real*4     :: MACHEP
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
               CALL TRED1(N,A,E,WORK(1),WORK(N+1))
               CALL TQLRAT(N,E,WORK(N+1),INFO,FIRST,MACHEP)
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
               lambdaM(k,i,j) = vlambda2(i,j,k)
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


      end subroutine


