module teclib

    type section
        integer :: ni
        integer :: nj
        integer :: nvar
        real,dimension(:,:),allocatable :: x
        real,dimension(:,:),allocatable :: y
        real,dimension(:,:),allocatable :: z
        real,dimension(:,:,:),allocatable :: value
    end type section
    
    type(section),dimension(:),allocatable :: sections

    real,dimension(:,:,:),allocatable :: envelop

    contains

    !######################################################
    subroutine init_sections(ns)
        
        integer,intent(in) :: ns
        
        allocate(sections(ns))
    
    end subroutine init_sections

    !######################################################
    subroutine init_section(id_section,ni,nj,nvar,x,y,z,value)

        integer,intent(in) :: id_section
        integer,intent(in) :: ni,nj,nvar
        real(8),dimension(ni,nj),intent(in) :: x,y,z
        real(8),dimension(ni,nj,nvar),intent(in) :: value

        integer :: ivar

        sections(id_section)%ni   = ni
        sections(id_section)%nj   = nj
        sections(id_section)%nvar = nvar

        allocate(sections(id_section)%x(ni,nj))
        allocate(sections(id_section)%y(ni,nj))
        allocate(sections(id_section)%z(ni,nj))
        allocate(sections(id_section)%value(ni,nj,nvar))

        sections(id_section)%x = x
        sections(id_section)%y = y
        sections(id_section)%z = z

        do ivar=1,nvar
            sections(id_section)%value(:,:,ivar) = value(:,:,ivar)
        enddo

    end subroutine init_section


    !######################################################
    subroutine clear_sections()
        deallocate(sections)
    end subroutine clear_sections


    !######################################################
    subroutine save_sections(simulname,numberletter,      &
                            length_outputDir,outputDir,   &
                            time,                         &
                            length_list_var,list_var,     &
                            numberSection                 )

        use global_var, only: timeStr_id

        implicit none 

        include 'tecio.inc'

        ! in var 
        integer,intent(in)           :: numberletter
        character(len=numberletter),intent(in) :: simulname
        
        integer,intent(in)           :: length_outputDir
        character(len=length_outputDir),intent(in) :: outputDir

        real(8),intent(in) :: time
        
        integer,intent(in) :: length_list_var
        character(len=length_list_var),intent(in) :: list_var

        integer,intent(in) :: numberSection

        ! loc var
        CHARACTER*1 NULCHAR
        Real*4,dimension(:,:,:),allocatable ::  XX, YY, ZZ
        Real*8    :: SolTime
        Integer*4 :: VIsDouble, FileType, FileFormat, DIsDouble
        Integer*4 :: ZoneType,StrandID,ParentZn,IsBlock
        Integer*4 :: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
        POINTER   (NullPtr,Null)
        Integer*4 :: Null(*)
        Integer*4 :: IMax, JMax, Kmax
        character*1 NULLCHR
        Integer*4   Debug,III
        integer :: igrass, ivar, i, j, k, kk, i_section,nvar
        character(len=2) :: numero
        integer :: tempint

        real*4,dimension(:),allocatable :: temprealA
        real*4,dimension(:,:,:),allocatable :: temprealA3d
        
        character(len=200) :: pltFile
        integer :: lpltFile, len_ZoneName
        character(len=34) :: ZoneName

        ! set cte
        !----------------------------------
        Debug     = 0
        VIsDouble = 0
        FileType  = 0
        FileFormat= 0
        DIsDouble = 0
        NULCHAR   = CHAR(0)
        NullPtr   = 0
        NULLCHR = NULCHAR

        ! define the name of the plt file
        !---------------------------------
        if ( time > 0 ) then 
            lpltFile = numberletter+length_outputDir+12+3+5+2
            tempint = int((time-int(time)) * 1.e2)
            write(pltfile,'(a,a,a,i6.6,a,i2.2,a)') outputDir(1:length_outputDir),simulname(1:numberletter),'_section_',&
                                                   int(time),'_',tempint,'.plt'
        else 
            lpltFile = numberletter+length_outputDir+12+3+5+3
            tempint = int((time-int(time)) * 1.e2)
            write(pltfile,'(a,a,a,i7.6,a,i2.2,a)') outputDir(1:length_outputDir),simulname(1:numberletter),'_section_',&
                                                   int(time),'_',tempint,'.plt'
        endif 

        I = TecIni111(simulname(1:numberletter)//NULLCHR, &
            'x y z '//list_var(1:length_list_var)//NULLCHR, &
            pltfile(1:lpltFile)//NULLCHR, &
            "."//NULLCHR, &
            FileType, &
            Debug, &
            VIsDouble)
       
        !loop over section:
        do i_section = 1, numberSection

            IMax = sections(i_section)%ni
            JMax = sections(i_section)%nj
            KMax = 1
            nvar = sections(i_section)%nvar

            allocate(XX(IMax,JMax,KMax))
            allocate(YY(IMax,JMax,KMax))
            allocate(ZZ(IMax,JMax,KMax))    

            do k =1,KMax
               do j =1,JMax
                  do i =1,IMax
                     XX(i,j,k) = sections(i_section)%x(i,j)
                     YY(i,j,k) = sections(i_section)%y(i,j)
                     ZZ(i,j,k) = sections(i_section)%z(i,j)
                  enddo
               enddo
            enddo

            ! write 3D
            !----------------
            ! define name of the zone

            tempint = int((time-int(time)) * 1.e1)
            if ( time > 0 ) then
                len_ZoneName = 22
                write(ZoneName,'(a,i6.6,a,i2.2,a,i2.2)')  't=', int(time),'.',tempint,'_section_',i_section
            else
                len_ZoneName = 23
                write(ZoneName,'(a,i7.6,a,i2.2,a,i2.2)')  't=', int(time),'.',tempint,'_section_',i_section
            endif
            SolTime =  time
            StrandID = i_section !timeStr_id
            ParentZn = 0
            I = TECZNE111(ZoneName(1:len_ZoneName)//NULCHAR, &
                 0,    & ! ZONETYPE
                 IMax, &
                 JMax, &
                 KMax, &
                 0,    &
                 0, &
                 0, &
                 SolTime, &
                 StrandID, &
                 ParentZn, &
                 1, &     ! ISBLOCK
                 0, &     ! NumFaceConnections
                 0, &     ! FaceNeighborMode
                 0, &     ! TotalNumFaceNodes
                 0, &     ! NumConnectedBoundaryFaces
                 0, &     ! TotalNumBoundaryConnections
                 Null, &  ! PassiveVarList
                 Null, &  ! ValueLocation
                 Null, &  ! ShareVarFromZone
                 0)       ! ShareConnectivityFromZone)

            III = IMax*JMax*KMax
            I   = TECDAT111(III,XX,DIsDouble)
            I   = TECDAT111(III,YY,DIsDouble)
            I   = TECDAT111(III,ZZ,DIsDouble)
    
            do ivar=1,nvar
              allocate(temprealA3d(IMax,JMax,KMax))
              temprealA3d(:,:,1) = sections(i_section)%value(:,:,ivar)
              I   = TECDAT111(III,temprealA3d,DIsDouble)
              deallocate(temprealA3d)
            enddo
            
            deallocate(XX)
            deallocate(YY)
            deallocate(ZZ)    

        enddo
        I = TecEnd111()


    end subroutine save_sections


    !######################################################
    subroutine save_flow(simulname,numberletter,flag_compute_lambda, &
                            length_outputDir,outputDir,                 &
                            time,                                       &
                            ngx,ngy,ngz,nvar,xg,yg,zg,                  &
                            QQ,                                    &
                            length_list_var_plt,list_var_plt            )

        use global_var, only: timeStr_id

        implicit none 

        include 'tecio.inc'

        ! in var 
        integer,intent(in)           :: numberletter
        character(len=numberletter),intent(in) :: simulname
        
        integer,intent(in)           :: length_outputDir
        character(len=length_outputDir),intent(in) :: outputDir

        logical,intent(in) :: flag_compute_lambda

        real(8),intent(in) :: time
        
        integer,intent(in) :: ngx,ngy,ngz,nvar
        real(8),dimension(ngx,ngy,ngz),intent(in) :: xg
        real(8),dimension(ngx,ngy,ngz),intent(in) :: yg
        real(8),dimension(ngx,ngy,ngz),intent(in) :: zg
        !character(len=30),dimension(nvar),intent(in)  :: vara

        real(8),dimension(ngx,ngy,ngz,nvar),intent(in) :: QQ

        integer,intent(in) :: length_list_var_plt
        character(len=length_list_var_plt),intent(in) :: list_var_plt


        ! loc var
        CHARACTER*1 NULCHAR
        Real*4,dimension(:,:,:),allocatable ::  XX, YY, ZZ
        Real*8    :: SolTime
        Integer*4 :: VIsDouble, FileType, FileFormat,  DIsDouble
        Integer*4 :: ZoneType,StrandID,ParentZn,IsBlock
        Integer*4 :: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
        POINTER   (NullPtr,Null)
        Integer*4 :: Null(*)
        Integer*4 :: IMax, JMax, Kmax
        character*1 NULLCHR
        Integer*4   Debug,III
        integer :: igrass, ivar, i, j, k, kk
        character(len=2) :: numero
        integer :: nvarPLT, tempint


        real*4,dimension(:),allocatable :: temprealA
        real*4,dimension(:,:,:),allocatable :: temprealA3d
        
        character(len=200) :: pltFile
        integer :: lpltFile, len_ZoneName
        character(len=18) :: ZoneName
        character(len=16) :: ZoneNamePart
        character(len=12) :: ZoneNameHem
        character(len=50) :: ZoneName_plan
        real,dimension(:,:,:),allocatable :: lambda

        !print*, vara

        ! compute lambda
        !----------------------------------
        if (flag_compute_lambda .eqv. .TRUE.) then 
            allocate(lambda(ngx,ngy,ngz))        
            lambda(:,:,:) = 0.e0
            call lambda2(ngy,ngz,ngx,xg,yg,zg,QQ(:,:,:,1),QQ(:,:,:,2),QQ(:,:,:,3),lambda)
        endif

        ! open plt file
        !----------------------------------
        Debug     = 0
        VIsDouble = 0
        FileType  = 0
        FileFormat= 0
        DIsDouble = 0
        NULCHAR   = CHAR(0)
        NullPtr   = 0
        NULLCHR = NULCHAR

        ! define the name of the plt file
        !---------------------------------
        if ( time > 0 ) then 
            lpltFile = numberletter+length_outputDir+12+3+2
            tempint = int((time-int(time)) * 1.e2)
            write(pltfile,'(a,a,a,i6.6,a,i2.2,a)') outputDir(1:length_outputDir),simulname(1:numberletter),'_3d_',&
                                                   int(time),'_',tempint,'.plt'
        else 
            lpltFile = numberletter+length_outputDir+12+3+3
            tempint = int((time-int(time)) * 1.e2)
            write(pltfile,'(a,a,a,i7.6,a,i2.2,a)') outputDir(1:length_outputDir),simulname(1:numberletter),'_3d_',&
                                                   int(time),'_',tempint,'.plt'
        endif 

        if (flag_compute_lambda .eqv. .TRUE.) then 
            I = TecIni111(simulname(1:numberletter)//NULLCHR, &
                'x y z '//list_var_plt(1:length_list_var_plt)//' lambda'//NULLCHR, &
                pltfile(1:lpltFile)//NULLCHR, &
                "."//NULLCHR, &
                FileType, &
                Debug, &
                VIsDouble)
            nvarPLT = nvar 
        else
            I = TecIni111(simulname(1:numberletter)//NULLCHR, &
                'x y z '//list_var_plt(1:length_list_var_plt)//NULLCHR, &
                pltfile(1:lpltFile)//NULLCHR, &
                "."//NULLCHR, &
                FileType, &
                Debug, &
                VIsDouble)
            nvarPLT = nvar
        endif
        IMax = ngx
        JMax = ngy
        KMax = ngz

        allocate(XX(IMax,JMax,KMax))
        allocate(YY(IMax,JMax,KMax))
        allocate(ZZ(IMax,JMax,KMax))    

        do k =1,KMax
           do j =1,JMax
              do i =1,IMax
                 XX(i,j,k) = xg(i,j,k)
                 YY(i,j,k) = yg(i,j,k)
                 ZZ(i,j,k) = zg(i,j,k)
              enddo
           enddo
        enddo

        ! write Flow data
        !----------------
        ! define name of the zone

        tempint = int((time-int(time)) * 1.e1)
        if (time > 0) then 
            len_ZoneName = 16
            write(ZoneName,'(a,i6.6,a,i2.2,a)')  't=', int(time),'.',tempint,'_flow'
        else
            len_ZoneName = 17
            write(ZoneName,'(a,i7.6,a,i2.2,a)')  't=', int(time),'.',tempint,'_flow'
        endif
        !timeStr_id = timeStr_id + 1

        SolTime =  time
        StrandID = timeStr_id
        ParentZn = 0
        I = TECZNE111(ZoneName(1:16)//NULCHAR, &
             0,    & ! ZONETYPE
             IMax, &
             JMax, &
             KMax, &
             0,    &
             0, &
             0, &
             SolTime, &
             StrandID, &
             ParentZn, &
             1, &     ! ISBLOCK
             0, &     ! NumFaceConnections
             0, &     ! FaceNeighborMode
             0, &     ! TotalNumFaceNodes
             0, &     ! NumConnectedBoundaryFaces
             0, &     ! TotalNumBoundaryConnections
             Null, &  ! PassiveVarList
             Null, &  ! ValueLocation
             Null, &  ! ShareVarFromZone
             0)       ! ShareConnectivityFromZone)

        III = IMax*JMax*KMax
        I   = TECDAT111(III,XX,DIsDouble)
        I   = TECDAT111(III,YY,DIsDouble)
        I   = TECDAT111(III,ZZ,DIsDouble)

        do ivar=1,nvarPLT
          allocate(temprealA3d(IMax,JMax,KMax))
          temprealA3d(:,:,:) = QQ(:,:,:,ivar)
          I   = TECDAT111(III,temprealA3d,DIsDouble)
          deallocate(temprealA3d)
        enddo

        if (flag_compute_lambda .eqv. .TRUE.) then 
          allocate(temprealA3d(IMax,JMax,KMax))
          temprealA3d(:,:,:) = lambda
          I   = TECDAT111(III,temprealA3d,DIsDouble)
          deallocate(temprealA3d)
        endif

        I = TecEnd111()

        deallocate(XX)
        deallocate(YY)
        deallocate(ZZ)    

    end subroutine save_flow


    !######################################################
    subroutine save_surface(simulname,numberletter,                     &
                            length_outputDir,outputDir,                 &
                            time,                                       &
                            ngx,ngy,nvar,xg,yg,orography,               &
                            QQ2d,                                       &
                            length_list_var_plt,list_var_plt            )

        use global_var, only: timeStr_id

        implicit none 

        include 'tecio.inc'

        ! in var 
        integer,intent(in)           :: numberletter
        character(len=numberletter),intent(in) :: simulname
        
        integer,intent(in)           :: length_outputDir
        character(len=length_outputDir),intent(in) :: outputDir

        real(8),intent(in) :: time
        
        integer,intent(in) :: ngx,ngy,nvar
        real(8),dimension(ngx,ngy),intent(in) :: xg
        real(8),dimension(ngx,ngy),intent(in) :: yg
        real(8),dimension(ngx,ngy),intent(in) :: orography
        !character(len=30),dimension(nvar),intent(in)  :: vara

        real(8),dimension(ngx,ngy,nvar),intent(in) :: QQ2d

        integer,intent(in) :: length_list_var_plt
        character(len=length_list_var_plt),intent(in) :: list_var_plt


        ! loc var
        CHARACTER*1 NULCHAR
        Real*4,dimension(:,:,:),allocatable ::  XX, YY, ZZ
        Real*8    :: SolTime
        Integer*4 :: VIsDouble, FileType, FileFormat, DIsDouble
        Integer*4 :: ZoneType,StrandID,ParentZn,IsBlock
        Integer*4 :: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
        POINTER   (NullPtr,Null)
        Integer*4 :: Null(*)
        Integer*4 :: IMax, JMax, Kmax
        character*1 NULLCHR
        Integer*4   Debug,III
        integer :: igrass, ivar, i, j, k, kk
        character(len=2) :: numero
        integer :: nvarPLT, tempint, ZoneName_len


        real*4,dimension(:),allocatable :: temprealA
        real*4,dimension(:,:),allocatable :: temprealA2d
        
        character(len=200) :: pltFile
        integer :: lpltFile
        character(len=18) :: ZoneName
        character(len=16) :: ZoneNamePart
        character(len=12) :: ZoneNameHem
        character(len=50) :: ZoneName_plan
        real,dimension(:,:,:),allocatable :: lambda


        ! open plt file
        !----------------------------------
        Debug     = 0
        VIsDouble = 0
        FileType  = 0
        FileFormat= 0
        DIsDouble = 0
        NULCHAR   = CHAR(0)
        NullPtr   = 0
        NULLCHR = NULCHAR

        ! define the name of the plt file
        !---------------------------------
        if (time > 0) then 
            lpltFile = numberletter+length_outputDir+12+3+2
            tempint = int((time-int(time)) * 1.e2)
            write(pltfile,'(a,a,a,i6.6,a,i2.2,a)') outputDir(1:length_outputDir),simulname(1:numberletter),'_2d_',&
                                                   int(time),'_',tempint,'.plt'
        else
            lpltFile = numberletter+length_outputDir+12+3+3
            tempint = int((time-int(time)) * 1.e2)
            write(pltfile,'(a,a,a,i7.6,a,i2.2,a)') outputDir(1:length_outputDir),simulname(1:numberletter),'_2d_',&
                                                   int(time),'_',tempint,'.plt'
        endif

        I = TecIni111(simulname(1:numberletter)//NULLCHR, &
            'x y z '//list_var_plt(1:length_list_var_plt)//NULLCHR, &
            pltfile(1:lpltFile)//NULLCHR, &
            "."//NULLCHR, &
            FileType, &
            Debug, &
            VIsDouble)
        nvarPLT = nvar
        
        IMax = ngx
        JMax = ngy
        KMax = 1

        allocate(XX(IMax,JMax,KMax))
        allocate(YY(IMax,JMax,KMax))
        allocate(ZZ(IMax,JMax,KMax))    

        do k =1,KMax
           do j =1,JMax
              do i =1,IMax
                 XX(i,j,k) = xg(i,j)
                 YY(i,j,k) = yg(i,j)
                 ZZ(i,j,k) = orography(i,j) 
              enddo
           enddo
        enddo

        ! write Surface data
        !----------------
        ! define name of the zone

        tempint = int((time-int(time)) * 1.e1)
        if (int(time) > 0 ) then 
            ZoneName_len = 16
            write(ZoneName,'(a,i6.6,a,i2.2,a)')  't=', int(time),'.',tempint,'_grd'
        else 
            ZoneName_len = 17
            write(ZoneName,'(a,i7.6,a,i2.2,a)')  't=', int(time),'.',tempint,'_grd'
        endif 

        !timeStr_id = timeStr_id + 1

        SolTime =  time
        StrandID = timeStr_id
        ParentZn = 0
        I = TECZNE111(ZoneName(1:ZoneName_len)//NULCHAR, &
             0,    & ! ZONETYPE
             IMax, &
             JMax, &
             KMax, &
             0,    &
             0, &
             0, &
             SolTime, &
             StrandID, &
             ParentZn, &
             1, &     ! ISBLOCK
             0, &     ! NumFaceConnections
             0, &     ! FaceNeighborMode
             0, &     ! TotalNumFaceNodes
             0, &     ! NumConnectedBoundaryFaces
             0, &     ! TotalNumBoundaryConnections
             Null, &  ! PassiveVarList
             Null, &  ! ValueLocation
             Null, &  ! ShareVarFromZone
             0)       ! ShareConnectivityFromZone)

        III = IMax*JMax*KMax
        I   = TECDAT111(III,XX,DIsDouble)
        I   = TECDAT111(III,YY,DIsDouble)
        I   = TECDAT111(III,ZZ,DIsDouble)

        do ivar=1,nvarPLT
          allocate(temprealA2d(IMax,JMax))
          temprealA2d(:,:) = QQ2d(:,:,ivar)
          I   = TECDAT111(III,temprealA2d,DIsDouble)
          deallocate(temprealA2d)
        enddo

        I = TecEnd111()

        deallocate(XX)
        deallocate(YY)
        deallocate(ZZ)    

    end subroutine save_surface


end module teclib

