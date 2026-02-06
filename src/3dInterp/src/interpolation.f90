subroutine interpolation(ngx,ngy,ngz,xg,yg,zg, &
                         ni,nj,nk,x,y,z,         &
                         nvar,q,QQ)
  use global_var

  implicit none

  integer,intent(in) :: ngx,ngy,ngz
  integer,intent(in) :: ni,nj,nk
  integer,intent(in) :: nvar

  real*8,dimension(ngx),intent(in) :: xg
  real*8,dimension(ngy),intent(in) :: yg
  real*8,dimension(ngz),intent(in) :: zg

  real*8,dimension(ni,nj,nk,nvar),intent(in) :: q  

  real*8,dimension(ni),intent(in) :: x
  real*8,dimension(nj),intent(in) :: y
  real*8,dimension(ni,nj,nk),intent(in) :: z

  real*8,dimension(ngx,ngy,ngz,nvar),intent(out) :: QQ  

  
  integer :: i,j,k,jj,kk,ij,ik
  integer :: ig,ige,jg,jge
  integer :: ivar
  real,dimension(2) :: xintp,yintp,fintp,zintp
  real :: xp,yp,zp
  real :: fxp,fxp1,fxp2
  
  real,dimension(:,:,:,:),allocatable :: QQi

!f2py intent(in) ngx,ngy,ngz,q,xg,yg,zg
!f2py intent(in) ni,nj,nk,x,y,z,nvar
!f2py intent(out) QQ

    flag_merde = .false.

  allocate(QQi(ngx,ngy,nk,nvar+1))
  if(flag_merde) write(*,*) 'alloc QQi'

!  write(*,*) xg(1:10)

!  write(*,*) 'in inter 1'
!  call min_max_3d(ni,nj,nk,q(:,:,:,1))

  !---------------
  ! interpolation over xy plan
  !---------------
  if(flag_merde)then
      write(*,*) 'merde interpolation'
  endif
    

  ig=1
  ige = ngx

  jg=1
  jge = ngy

  if(flag_merde) write(*,*) xg(ig),xg(ige), x(1),x(ni)
  if(flag_merde) write(*,*) ig,ige,ngx
  if(flag_merde) write(*,*) yg(jg),yg(jge), y(1),y(nj)
  if(flag_merde) write(*,*) jg,jge,ngy

  do  i = ig, ige

     xp = xg(i)

     ij = 1
     do while( x(ij) .le. xp )
        ij = ij+1
        if(ij .gt. 1.e5) then 
           write(*,*) 'problem in interpolation over x'
           stop
        endif
     enddo
     ij = ij - 1

     xintp(1) = x(ij)
     xintp(2) = x(ij+1)

!     write(*,*) 'xp',xintp(1),xp, xintp(2),ij

     do  j = jg, jge

        yp = yg(j)
!        write(*,*) yp, jg,j,jge

        jj = 1
        do while(  y(jj) .le. yp)
!           write(*,*) jj
           if(jj .gt. nj) then 
              write(*,*) 'problem in interpolation over y'
              write(*,*) jj,yp
              write(*,*) (y(kk),kk=1,nj)
              stop
           endif
           jj = jj + 1
        enddo
        jj = jj - 1

!        write(*,*) 'merde',jj


        yintp(1) = y(jj)
        yintp(2) = y(jj+1)

!        write(*,*) 'yp',yintp(1),yp, yintp(2),jj


        do k=1,nk

           do ivar=1,nvar+1

              if(xp .eq. xintp(1) .and. yp .eq. yintp(1)) then

                  if (ivar == nvar+1) then 
                      fxp = z(ij,jj,k)
                  else
                      fxp = q(ij,jj,k,ivar)
                  endif
              elseif ( yp .eq. yintp(1) .and. xp .ne. xintp(1) ) then

                 !                       write(*,*) yintp(1),yp, yintp(2),jj
                 !                       write(*,*) xintp(1),xp, xintp(2),ij                                  

                 if (ivar == nvar+1) then                  
                    fintp(1)= z(ij ,jj,k)
                    fintp(2)= z(ij+1,jj,k)
                 else
                    fintp(1)= q(ij  ,jj,k,ivar) 
                    fintp(2)= q(ij+1,jj,k,ivar)
                 endif

                 fxp = (fintp(2)-fintp(1))/(xintp(2)-xintp(1)) * (xp-xintp(1)) + fintp(1)

                 !                       write(*,*) fintp(1),fintp(2),fxp

              elseif ( xp .eq. xintp(1) .and.  yp .ne. yintp(1) ) then

                  if (ivar == nvar+1) then                  
                    fintp(1)= z(ij,jj  ,k)
                    fintp(2)= z(ij,jj+1,k)
                 else
                    fintp(1)= q(ij,jj  ,k,ivar) 
                    fintp(2)= q(ij,jj+1,k,ivar) 
                 endif
                 fxp = (fintp(2)-fintp(1))/(yintp(2)-yintp(1)) * (yp-yintp(1)) + fintp(1)

              elseif( xp .ne. xintp(1) .and.  yp .ne. yintp(1) ) then

                  if (ivar == nvar+1) then                  
                     fintp(1)= z(ij,  jj,k)
                     fintp(2)= z(ij+1,jj,k)
                 else
                     fintp(1)= q(ij  ,jj,k,ivar) 
                     fintp(2)= q(ij+1,jj,k,ivar) 
                 endif

                 fxp1 = (fintp(2)-fintp(1))/(xintp(2)-xintp(1)) * (xp-xintp(1)) + fintp(1)

                 if (ivar == nvar+1) then                  
                     fintp(1)= z(ij,  jj+1,k)
                     fintp(2)= z(ij+1,jj+1,k)
                 else
                     fintp(1)= q(ij  ,jj+1,k,ivar) 
                     fintp(2)= q(ij+1,jj+1,k,ivar) 
                 endif

                 fxp2 = (fintp(2)-fintp(1))/(xintp(2)-xintp(1)) * (xp-xintp(1)) + fintp(1)                       
                 fxp = (fxp2-fxp1)/(yintp(2)-yintp(1)) * (yp-yintp(1)) + fxp1

              else

                 write(*,*) 'merde interpolation'
                 stop

              endif

              QQi(i,j,k,ivar) = fxp

!              write(*,*) ige,jge,nk,nvar
!              write(*,*) 'm',i,j,k,ivar
           enddo
        enddo

     enddo
  enddo

  if(flag_merde) write(*,*) 'interolation over z'

  !---------------
  ! interpolation over z
  !---------------

!  write(*,*) 'zg',(zg(kk),kk=1,ngz)
!  write(*,*) 'z',(z(kk),kk=1,nk)
 if (ngz .ne. 1) then 

  do ivar=1,nvar

    !print*, ivar
    do  i = ig, ige
       do  j = jg, jge
         do k = 1,ngz
       
              zp = zg(k)
              if( (zp .gt. QQi(i,j,nk,nvar+1)) .or. (zp .lt. QQi(i,j,1,nvar+1)) )then
                QQ(i,j,k,ivar) = -999
                cycle
              endif

              ik = 1
              do while( QQi(i,j,ik,nvar+1) .le. zp )
                  if( ik .ge. nk+1 .and. zp .gt. QQi(i,j,ik,nvar+1) ) then 
                     write(*,*) 'problem in interpolation over z'
                     write(*,*) ik,QQi(i,j,ik,nvar+1),zp
                     stop
                  endif
                  ik = ik+1
              enddo
              ik = ik - 1
              
              zintp(1) = QQi(i,j,ik,nvar+1)
              zintp(2) = QQi(i,j,ik+1,nvar+1)

              if(zp .eq. zintp(1) ) then

                 fxp = QQi(i,j,ik,ivar)

              else

                 fintp(1)= QQi(i,j,ik,  ivar) 
                 fintp(2)= QQi(i,j,ik+1,ivar) 

                 fxp = (fintp(2)-fintp(1))/(zintp(2)-zintp(1)) * (zp-zintp(1)) + fintp(1)

              endif

              QQ(i,j,k,ivar) = fxp

           enddo
        enddo

     enddo
     
     !print*, 'ivar done'

  enddo
  else  ! 2-D
      do ivar=1,nvar
        QQ(ig:ige,jg:jge,1,ivar) = QQi(ig:ige,jg:jge,1,ivar)
      enddo
  endif

  deallocate(QQi)

!  write(*,*) 'in inter 2'
! call min_max_3d(ngx,ngy,ngz,QQ(:,:,:,1))


end subroutine interpolation
