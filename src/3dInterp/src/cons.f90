module global_var

  integer :: LU_rfM,LU_cv,LU_maxminTf,Lu_mass,Lu_frp_f,Lu_err,Lu_frp_compa,Lu_simlPic
!  integer,dimension(:),allocatable :: LU_frp_a

  real,dimension(:),allocatable :: phi_pt,theta_pt
  
  real :: MassFtmdt,MassFini

  logical :: flag_firstP,flag_merde,flag_part_Glo
  logical,dimension(:),allocatable :: flag_part_G

  real :: xc,yc,zc

  real :: pi, piov2, piov8

  integer :: timeStr2D

  integer :: lu_open

  real :: pict_xb,pict_xe,pict_yb,pict_ye ! dimension of the picture in meter

end module global_var
