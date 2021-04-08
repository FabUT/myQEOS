module mod_ensight

use mod_data, only : iout, PR

implicit none
!!!----- Objets liés au maillage
type nod
   real(PR) :: x=0_PR
   real(PR) :: y=0_PR
   real(PR) :: z=0_PR
   integer :: id
   integer :: Nel=0
   logical :: tag
end type nod

type elem
   integer :: n(1:8) !!! liste des noeuds
   integer :: Nnod=0 !! 1: point, 2: line ...
   integer :: part, coul
   integer :: typ
   integer :: id
end type elem

!!!----- Objet maillage :

type mesh

   !!!---dimension du maillage 1:1D;2:2D;3:3D;4:2Daxi
   integer :: Ndim=3
   integer :: Npart !!! nombre de partitions

   !!!---Nombres de noeuds, d'éléments, d'edges et de faces 
   integer :: Nnod =0, Nel =0, Nedge =0, Nfac =0 ! local 
   integer :: Nnodg =0, Nelg =0, Nedgeg =0, Nfacg =0 ! global  
   integer :: Nnodt=0, Nelt=0 !!! total
   integer :: Nnodf,  Nelf=0  !!! fantôme
   integer :: Np=0 , Nl=0 , Ntri=0 , Nquad=0 , Ntetra=0, Nhexa=0, Npenta=0, Npyra=0
   integer :: Npg=0 , Nlg=0 , Ntrig=0 , Nquadg=0 , Ntetrag=0, Nhexag=0, Npentag=0, Npyrag=0
   integer :: Npf=0 , Nlf=0 , Ntrif=0 , Nquadf=0 , Ntetraf=0, Nhexaf=0, Npentaf=0, Npyraf=0
   integer :: Npt=0 , Nlt=0 , Ntrit=0 , Nquadt=0 , Ntetrat=0, Nhexat=0, Npentat=0, Npyrat=0
  
   !!!---elements geo:
   type(elem), allocatable :: el(:)
  
   !!!---connecteurs:
   type(nod), allocatable :: nods(:)
  
   !!!---grandeurs geometriques
   real(PR), allocatable :: xyzel(:,:)
 
   !!!---scalaires
   real(PR), allocatable :: P(:), T(:), E(:), Cv(:), Cp(:), c(:), nu(:), g(:), w(:), Z(:)
   real(PR), allocatable :: sig(:), lambda(:), Seebeck(:)
 
end type mesh


type output 
   character(len=50) :: folder, folder_loc
   character(len=50) :: fcase='CHR'
   integer :: Nscal, Nvec
   character(len=20) :: scal(1:20)
   character(len=20) :: vec(1:20)
end type output

integer :: it=1

contains

!!!================================================================================================
!!!
!!!                                     SORTIES ENSIGHT
!!!
!!!================================================================================================


!!!===================================================================================
!!!                         SORTIES SEQUENTIELLE ENSIGHT AUX ELEMENTS
!!!===================================================================================

subroutine OUTPUT_Ensight(N_T,N_nu,nu_sat,T,nu,P,E,Cv,Cp,c,g,w,Z,sig,lambda,Seebeck,imode)

    implicit none
    integer, intent(in) :: N_T, N_nu
    real(PR), intent(in) :: nu_sat(1:2,1:N_T)
    real(PR), intent(in) :: T(1:N_T)
    real(PR), intent(in) :: nu(1:N_nu)
    real(PR), intent(in) :: P(1:N_nu,1:N_T)
    real(PR), intent(in) :: E(1:N_nu,1:N_T)
    real(PR), intent(in) :: Cv(1:N_nu,1:N_T)
    real(PR), intent(in) :: Cp(1:N_nu,1:N_T)
    real(PR), intent(in) :: c(1:N_nu,1:N_T)
    real(PR), intent(in) :: g(1:N_nu,1:N_T)
    real(PR), intent(in) :: w(1:N_nu,1:N_T)
    real(PR), intent(in) :: Z(1:N_nu,1:N_T)
    real(PR), intent(in) :: sig(1:N_nu,1:N_T)
    real(PR), intent(in) :: lambda(1:N_nu,1:N_T)
    real(PR), intent(in) :: Seebeck(1:N_nu,1:N_T)

    integer, intent(in) :: imode
 
    real(PR) :: t1,t2
    integer :: i,j,ie,iel
    character(len=100) :: dump, repfile, sorfolder
    integer :: Nel
    real(PR), allocatable :: rtemp(:), vtemp(:,:)

    type(mesh) :: M
    type(output) :: Sor

    write(iout,*) ''
    write(iout,*) '  sortie aux élements'

    call cpu_time(t1)
 
    call create_mesh_from_tab(N_T,N_nu,nu_sat,T,nu,P,E,Cv,Cp,c,g,w,Z,sig,lambda,Seebeck,M,imode)


    if(imode.eq.1)then
      Sor%folder='out_ensight_P'
    else
      Sor%folder='out_ensight_s'
    endif
    Sor%fcase='CHR'
    Sor%Nvec=0
    Sor%Nscal=13
    Sor%scal(1)='nu'
    Sor%scal(2)='T'
    Sor%scal(3)='P'
    Sor%scal(4)='E'
    Sor%scal(5)='Cv'
    Sor%scal(6)='Cp'
    Sor%scal(7)='c'
    Sor%scal(8)='g'
    Sor%scal(9)='w'
    Sor%scal(10)='Z'
    Sor%scal(11)='sig'
    Sor%scal(12)='lambda'
    Sor%scal(13)='Seebeck'



    if(it.eq.1)then 
        call WriteEnsightMeshASCII(M,trim(adjustl(Sor%folder)),trim(adjustl(Sor%fcase)))
    endif

    call EnsightCase(Sor,it)


     if(Sor%Nscal.gt.0)then

      !!!---scalaires
      do i=1,Sor%Nscal

        if(trim(adjustl(Sor%scal(i))).eq.'T')then
            call WriteScal(M%T(1:Nel),M,trim(adjustl(Sor%folder)),'T',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'nu')then
            call WriteScal(M%nu(1:Nel),M,trim(adjustl(Sor%folder)),'nu',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'P')then
            call WriteScal(M%P(1:Nel),M,trim(adjustl(Sor%folder)),'P',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'E')then
            call WriteScal(M%E(1:Nel),M,trim(adjustl(Sor%folder)),'E',it)
        endif

        if(trim(adjustl(Sor%scal(i))).eq.'Cv')then
            call WriteScal(M%Cv(1:Nel),M,trim(adjustl(Sor%folder)),'Cv',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'Cp')then
            call WriteScal(M%Cp(1:Nel),M,trim(adjustl(Sor%folder)),'Cp',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'c')then
            call WriteScal(M%c(1:Nel),M,trim(adjustl(Sor%folder)),'c',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'g')then
            call WriteScal(M%g(1:Nel),M,trim(adjustl(Sor%folder)),'g',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'w')then
            call WriteScal(M%w(1:Nel),M,trim(adjustl(Sor%folder)),'w',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'Z')then
            call WriteScal(M%Z(1:Nel),M,trim(adjustl(Sor%folder)),'Z',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'sig')then
            call WriteScal(M%sig(1:Nel),M,trim(adjustl(Sor%folder)),'sig',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'lambda')then
            call WriteScal(M%lambda(1:Nel),M,trim(adjustl(Sor%folder)),'lambda',it)
        endif
        if(trim(adjustl(Sor%scal(i))).eq.'Seebeck')then
            call WriteScal(M%Seebeck(1:Nel),M,trim(adjustl(Sor%folder)),'Seebeck',it)
        endif





      enddo

     endif

   deallocate(M%nods)
   deallocate(M%el)
   deallocate(M%T)
   deallocate(M%nu)
   deallocate(M%P)
   deallocate(M%E)
   deallocate(M%Cv)
   deallocate(M%Cp)
   deallocate(M%c)
   deallocate(M%g)
   deallocate(M%w)
   deallocate(M%Z)
   deallocate(M%sig)
   deallocate(M%lambda)
   deallocate(M%Seebeck)

   call cpu_time(t2) 
   write(iout,*) '  sortie aux éléments fini:', t2-t1, 's'
   write(iout,*) '' 

end subroutine OUTPUT_Ensight

subroutine WriteEnsightMeshASCII(M,folder,fcase)

     implicit none
     type(mesh), intent(in) :: M
     character(len=*),intent(in) :: folder,fcase
     integer :: i, iens, ke, ierr
     character(len=100) :: f
     integer :: Nnod,Nel,Np,Nl,Ntri,Nquad,Ntetra,Nhexa,Npenta,Npyra
     !..............................................................................
      iens = 1
      
      f=trim(adjustl(folder))//'/'//trim(adjustl(fcase))//'.geo'
     
     ! write(iout,*) 'opening file: ', trim(adjustl(f))
     
      open (unit=iens,file=trim(adjustl(f)),form='formatted',&
      status='replace',action='write',iostat=ierr)
     
     !     Header
                            
         write(iout,'(5x,A)') '> Exporting mesh'
         write(iout,'(10x,A)') trim(adjustl(f))


         write (iens,'(A)') 'description line 1'
         write (iens,'(A)') 'description line 2'
         write (iens,'(A)') 'node id assign'
         write (iens,'(A)') 'element id assign'
         write (iens,'(A)') 'part'
         write (iens,'(I10)') 1
         write (iens,'(A)') 'description line'
         write (iens,'(A)') 'coordinates'
     
         Nnod  =M%Nnod 
         Nel   =M%Nel  
         Np    =M%Np   
         Nl    =M%Nl   
         Ntri  =M%Ntri 
         Nquad =M%Nquad 
         Ntetra=M%Ntetra
         Nhexa =M%Nhexa 
         Npenta=M%Npenta
         Npyra =M%Npyra 
     
     
         write (iens,'(I10)') Nnod
     !    do i=1,Nnod
     !        write (iens,'(I10)') i 
     !    enddo
         do i=1,Nnod
             write (iens,'(E12.5)') M%nods(i)%x  
         enddo
         do i=1,Nnod
             write (iens,'(E12.5)') M%nods(i)%y 
         enddo
         do i=1,Nnod
             write (iens,'(E12.5)') M%nods(i)%z
         enddo
     
     
     !!!-------Points
         if(Np.gt.0)then
         write (iens,'(A)') 'point'
         write (iens,'(I10)') Np
          do i=1,Nel
            if(M%el(i)%typ.eq.15)then
             write (iens,'(I10)') M%el(i)%n(1)
            endif
          enddo
         endif
     
     !!!-------lignes
         if(Nl.gt.0)then
         write (iens,'(A)') 'bar2'
         write (iens,'(I10)') Nl
          do i=1,Nel
            if(M%el(i)%typ.eq.1)then
             write (iens,'(2I10)') M%el(i)%n(1) , M%el(i)%n(2)
            endif
          enddo
         endif
     
     !!!-------triangles
         if(Ntri.gt.0)then
         write (iens,'(A)') 'tria3'
         write (iens,'(I10)') Ntri
          do i=1,Nel
            if(M%el(i)%typ.eq.2)then
             write (iens,'(3I10)') M%el(i)%n(1) , M%el(i)%n(2), M%el(i)%n(3)
            endif
          enddo
         endif
     
     !!!-------quads
         if(Nquad.gt.0)then
         write (iens,'(A)') 'quad4'
         write (iens,'(I10)') Nquad
          do i=1,Nel
            if(M%el(i)%typ.eq.3)then
             write (iens,'(4I10)') M%el(i)%n(1) , M%el(i)%n(2), M%el(i)%n(3), M%el(i)%n(4)
            endif
          enddo
         endif
     
     !!!-------tetras
         if(Ntetra.gt.0)then
         write (iens,'(A)') 'tetra4'
         write (iens,'(I10)') Ntetra
          do i=1,Nel
            if(M%el(i)%typ.eq.4)then
             write (iens,'(4I10)') M%el(i)%n(1) , M%el(i)%n(2), M%el(i)%n(3), M%el(i)%n(4)
            endif
          enddo
         endif
     
     !!!-------hexas
         if(Nhexa.gt.0)then
         write (iens,'(A)') 'hexa8'
         write (iens,'(I10)') Nhexa
          do i=1,Nel
            if(M%el(i)%typ.eq.5)then
             write (iens,'(8I10)') M%el(i)%n(1),M%el(i)%n(2), M%el(i)%n(3), M%el(i)%n(4),&
                                   M%el(i)%n(5),M%el(i)%n(6), M%el(i)%n(7), M%el(i)%n(8)
            endif
          enddo
         endif
     
     !!!-------pentas
         if(Npenta.gt.0)then
         write (iens,'(A)') 'penta6'
         write (iens,'(I10)') Npenta
          do i=1,Nel
            if(M%el(i)%typ.eq.6)then
             write (iens,'(6I10)') M%el(i)%n(1),M%el(i)%n(2), M%el(i)%n(3), M%el(i)%n(4),&
                                   M%el(i)%n(5),M%el(i)%n(6)
            endif
          enddo
         endif
     
     !!!-------pyras
         if(Npyra.gt.0)then
         write (iens,'(A)') 'pyra5'
         write (iens,'(I10)') Npyra
          do i=1,Nel
            if(M%el(i)%typ.eq.7)then
             write (iens,'(5I10)') M%el(i)%n(1),M%el(i)%n(2), M%el(i)%n(3), M%el(i)%n(4),&
                                   M%el(i)%n(5)
            endif
          enddo
         endif
     
        close(iens)
     
end subroutine WriteEnsightMeshASCII

subroutine EnsightCase(Sor,it)
     
     implicit none
     
     !  EnsightCase helps to write a Ensight's case file
     !  radical...: word used to build filenames
     !  nf........: number of result files (time steps)
     !      character(*), intent(in) ::  radical
     
     type(output), intent(in) :: Sor
     integer, intent(in) :: it
     integer :: icase, i
     character(len=100) :: arch1, arch2
     
     write(iout,'(5x,A)') '> Exporting casefile'

     icase = 35
     open(icase,file=trim(adjustl(Sor%folder))//'/'//trim(adjustl(Sor%fcase))//'.case',status='replace')
     
     !    Writing Introdutory lines
     write(icase,10) trim(adjustl(Sor%fcase))//'.geo'
     
     if(Sor%Nscal.gt.0)then
      !!!---scalaires
      do i=1,Sor%Nscal
        arch2 = trim(adjustl(Sor%scal(i)))//"****"//".scl"    !!!!'0001.scl'
        write(icase,16) trim(adjustl(Sor%scal(i))),'       ', trim(adjustl(arch2))
      enddo
     endif

     if(Sor%Nvec.gt.0)then
     !!!---vecteurs
      do i=1,Sor%Nvec
        arch2 = trim(adjustl(Sor%vec(i)))//"****"//".vec"
        write(icase,26) trim(adjustl(Sor%vec(i))),'       ', trim(adjustl(arch2))
      enddo
     endif
 
     !!!---TIME sequence: 
     write(icase,*) ''
     write(icase,'(A)') 'TIME'
     write(icase,'(A)') 'time set:              1'
     write(icase,'(A17,I7)') 'number of steps: ', it
     write(icase,'(A)') 'filename start number: 1'
     write(icase,'(A)') 'filename increment:    1'
     write(icase,'(A)') 'time values:'
     do i=1,it
     write(icase,'(A6,I7)') '            ', i
     enddo
     
     close(icase)
     
     !!!------formats :
     10 format(&
       'FORMAT'        ,/ ,&
       'type:	ensight gold',//,&
       'GEOMETRY'      ,/ ,&
       'model:	'       ,A20,//,&
       'VARIABLE'      )
     
     !15 format('scalar per node:    ',A20,A7,A)
     !25 format('vector per node:    ',A20,A7,A)
     16 format('scalar per element: ',A20,A7,A)
     26 format('vector per element: ',A20,A7,A)
     
     
     
end subroutine EnsightCase

subroutine WriteScal(p,M,folder,radical,it)

     implicit none

     type(mesh), intent(in) :: M
     integer, intent(in) :: it
     real(PR), dimension(1:M%Nelt), intent(in) :: p 
     character(*), intent(in) :: folder,radical
     character(len=100) :: arch1, arch2
     character(len=10) :: dump
     integer :: iens, i
     integer :: Nnod,Nel,Np,Nl,Ntri,Nquad,Ntetra,Nhexa,Npenta,Npyra

     Nnod  = M%Nnod   
     Nel   = M%Nel    
     Np    = M%Np     
     Nl    = M%Nl     
     Ntri  = M%Ntri   
     Nquad = M%Nquad  
     Ntetra= M%Ntetra 
     Nhexa = M%Nhexa  
     Npenta= M%Npenta 
     Npyra = M%Npyra  

     iens = 21

     write(iout,'(5x,A,A)') '> Exporting scalar ', trim(adjustl(radical))

     ! SCALAR FILE (PRESSURE)
     write(dump,'(I0.4)') it
     arch2 = trim(adjustl(folder))//'/'//trim(radical)//trim(adjustl(dump))//'.scl'
     write(iout,'(5x,A)') arch2
     open (unit=iens,file=arch2, form='formatted')
     write(iens,'(A)') 'description line 1'
     write(iens,'(A)') 'part'
     write(iens,'(I10)') 1
     
     !!!-------Points
     if(Np.gt.0)then
       write (iens,'(A)') 'point'
       do i=1,Nel
        if(M%el(i)%typ.eq.15)then
         write (iens,'(E12.5)') real(p(i))
        endif
       enddo
     endif

     !!!-------lignes
     if(Nl.gt.0)then
       write (iens,'(A)') 'bar2'
       do i=1,Nel
         if(M%el(i)%typ.eq.1)then
          write (iens,'(E12.5)') real(p(i))
         endif
       enddo
     endif

     !!!-------triangles
     if(Ntri.gt.0)then
       write (iens,'(A)') 'tria3'
       do i=1,Nel
         if(M%el(i)%typ.eq.2)then
          write (iens,'(E12.5)') real(p(i))
         endif
       enddo
     endif

     !!!-------quads
     if(Nquad.gt.0)then
       write (iens,'(A)') 'quad4'
       do i=1,Nel
         if(M%el(i)%typ.eq.3)then
          write (iens,'(E12.5)') real(p(i))
         endif
       enddo
     endif
 
     !!!-------tetra4
     if(Ntetra.gt.0)then
       write (iens,'(A)') 'tetra4'
       do i=1,Nel
        if(M%el(i)%typ.eq.4)then
         write (iens,'(E12.5)') real(p(i))
        endif
       enddo
     endif
 
     !!!-------hexa8
     if(Nhexa.gt.0)then
       write (iens,'(A)') 'hexa8'
       do i=1,Nel
        if(M%el(i)%typ.eq.5)then
         write (iens,'(E12.5)') real(p(i))
        endif
       enddo
     endif
 
     !!!-------penta
     if(Npenta.gt.0)then
       write (iens,'(A)') 'penta6'
       do i=1,Nel
        if(M%el(i)%typ.eq.6)then
         write (iens,'(E12.5)') real(p(i))
        endif
       enddo
     endif
 
     !!!-------pyra
     if(Npyra.gt.0)then
       write (iens,'(A)') 'pyramid5'
       do i=1,Nel
        if(M%el(i)%typ.eq.7)then
         write (iens,'(E12.5)') real(p(i))
        endif
       enddo
     endif
 
     close(iens)

end subroutine WriteScal

subroutine WriteVect(p,M,folder,radical,it)

      implicit none

      type(mesh), intent(in) :: M
      integer, intent(in) :: it
      real(PR), dimension(1:3,1:M%Nelt), intent(in) :: p 
      character(*), intent(in) :: folder,radical
        
      character(len=50) :: arch1, arch2
      character(len=10) :: dump
      integer :: iens, i,k
      integer :: Nnod,Nel,Np,Nl,Ntri,Nquad,Ntetra,Nhexa,Npenta,Npyra
!..............................................................................

      Nnod  = M%Nnodt   !max(M%Nnod,M%Nnodt)
      Nel   = M%Nelt    !max(M%Nel,M%Nelt)
      Np    = M%Npt     !max(M%Np,M%Npt)
      Nl    = M%Nlt     !max(M%Nl,M%Nlt)
      Ntri  = M%Ntrit   !max(M%Ntri,M%Ntrit)
      Nquad = M%Nquadt  !max(M%Nquad,M%Nquadt)
      Ntetra= M%Ntetrat !max(M%Ntetra,M%Ntetrat)
      Nhexa = M%Nhexat  !max(M%Nhexa,M%Nhexat)
      Npenta= M%Npentat !max(M%Npenta,M%Npentat)
      Npyra = M%Npyrat  !max(M%Npyra,M%Npyrat)


      iens = 21

      write(iout,'(5x,A,A)') '> Exporting vector ', trim(adjustl(radical))


      ! SCALAR FILE (PRESSURE)
      write(dump,'(I0.4)') it
      arch2 = trim(adjustl(folder))//'/'//trim(radical)//trim(adjustl(dump))//'.vec'
      write(iout,'(10x,A)') arch2
      open (unit=iens,file=arch2, form='formatted')

      write(iens,'(A)') 'description line 1'
      write(iens,'(A)') 'part'
      write(iens,'(I10)') 1

      if(Np.gt.0)then
      write(iens,'(A)') 'point'
      do k=1,3
        do i=1,Nel
          if(M%el(i)%typ.eq.15) write(iens,'(E12.5)') real(p(k,i))
        enddo
      enddo
      endif

      if(Nl.gt.0)then
      write(iens,'(A)') 'bar2'
      do k=1,3
        do i=1,Nel
          if(M%el(i)%typ.eq.1) write(iens,'(E12.5)') real(p(k,i))
        enddo
      enddo
      endif

      if(Ntri.gt.0)then
      write(iens,'(A)') 'tria3'
      do k=1,3
         do i=1,Nel
           if(M%el(i)%typ.eq.2) write(iens,'(E12.5)') real(p(k,i))
         enddo
      enddo
      endif

      if(Nquad.gt.0)then
      write(iens,'(A)') 'quad4'
      do k=1,3
         do i=1,Nel
           if(M%el(i)%typ.eq.3) write(iens,'(E12.5)') real(p(k,i))
         enddo
      enddo
      endif

     !!!-------tetra4
     if(Ntetra.gt.0)then
     write (iens,'(A)') 'tetra4'
     do k=1,3
        do i=1,Nel
         if(M%el(i)%typ.eq.4) write (iens,'(E12.5)') real(p(k,i))
        enddo
     enddo
     endif

     !!!-------hexa8
     if(Nhexa.gt.0)then
     write (iens,'(A)') 'hexa8'
     do k=1,3
        do i=1,Nel
         if(M%el(i)%typ.eq.5) write (iens,'(E12.5)') real(p(k,i))
        enddo
     enddo
     endif

     !!!-------penta
     if(Npenta.gt.0)then
     write (iens,'(A)') 'penta6'
     do k=1,3
        do i=1,Nel
         if(M%el(i)%typ.eq.6) write (iens,'(E12.5)') real(p(k,i))
        enddo
     enddo
     endif

     !!!-------pyra
     if(Npyra.gt.0)then
     write (iens,'(A)') 'pyramid5'
     do k=1,3
        do i=1,Nel
         if(M%el(i)%typ.eq.7) write (iens,'(E12.5)') real(p(k,i))
        enddo
     enddo
     endif


     !..............................................................................

     close(iens)


endsubroutine WriteVect


!!!===================================================================================
!!!                         REPRISES ENSIGHT AUX ELEMENTS
!!!===================================================================================

subroutine OpenScal(p,Np,Nl,Nt,Nq,folder,radical,it)

      implicit none

      integer, intent(in) :: Np,Nl,Nt,Nq,it
      real(PR), dimension(Np+Nl+Nt+Nq), intent(out) :: p 
      character(*), intent(in) :: folder,radical
      character(len=100) :: arch1, arch2
      character(len=10) :: dump
      integer :: iscl1, i
!..............................................................................


      iscl1 = 21

      write(iout,'(/,A   )') '(ENSIGHTIN) Importing Ensight result files'

      ! SCALAR FILE (PRESSURE)
      write(dump,'(I0.4)') it
      arch2 = trim(adjustl(folder))//'/'//trim(radical)//trim(adjustl(dump))//'.scl'
      write(iout,'(5x,A)') arch2
      open (unit=iscl1,file=arch2, form='formatted',status='old',action='read')
      read(iscl1,*) dump !write(iscl1,'(A)') 'description line 1'
      read(iscl1,*) dump !write(iscl1,'(A)') 'part'
      read(iscl1,*) dump !write(iscl1,'(I10)') 1
      
      !!!-------Points
      if(Np.gt.0)then

      read(iscl1,*) dump !!!!!write (iscl1,'(A)') 'point'
      do i=1,Np
       read(iscl1,'(E12.5)') p(i)
      enddo
      endif

      !!!-------lignes
      if(Nl.gt.0)then
      read(iscl1,*) dump !write (iscl1,'(A)') 'bar2'
      do i=Np+1,Np+Nl
        read(iscl1,'(E12.5)') p(i)
      enddo
      endif

      !!!-------triangles
      if(Nt.gt.0)then
      read(iscl1,*) dump !write (iscl1,'(A)') 'tria3'
      do i=Np+Nl+1,Np+Nl+Nt
        read(iscl1,'(E12.5)') p(i)
      enddo
      endif

     !!!-------quads
     if(Nq.gt.0)then
     read(iscl1,*) dump  !write (iscl1,'(A)') 'quad4'
     do i=Np+Nl+Nt+1,Np+Nl+Nt+Nq
       read(iscl1,'(E12.5)') p(i)
     enddo
     endif



      !write(iscl1,'(A)') 'coordinates'
      !do i=1,nnos
      !  write(iscl1,'(E12.5)') p(i)
      !enddo
     
      close(iscl1)

      !..............................................................................

endsubroutine OpenScal

subroutine OpenVect (vx,vy,vz,Np,Nl,Nt,Nq,folder,radical,it)

      implicit none

      integer, intent(in) :: Np,Nl,Nt,Nq, it
      real(PR), dimension(Np+Nl+Nt+Nq), intent(out) :: vx,vy,vz 
      character(*), intent(in) :: folder,radical
        
      character(len=50) :: arch1, arch2
      character(len=10) :: dump
      integer :: iscl1, i
!..............................................................................

      iscl1 = 21

      write(iout,'(/,A   )') '(ENSIGHTOUT) Importing Ensight result files'

      ! SCALAR FILE (PRESSURE)
      write(dump,'(I0.4)') it
      arch2 = trim(adjustl(folder))//'/'//trim(radical)//trim(adjustl(dump))//'.vec'
      write(iout,'(5x,A)') arch2
      open (unit=iscl1,file=arch2, form='formatted',status='old',action='read')

      read(iscl1,*) dump ! write(iscl1,'(A)') 'description line 1'
      read(iscl1,*) dump ! write(iscl1,'(A)') 'part'
      read(iscl1,*) dump ! write(iscl1,'(I10)') 1

      if(Np.gt.0)then
      read(iscl1,*) dump !write(iscl1,'(A)') 'point'
      do i=1,Np
        read(iscl1,'(E12.5)') vx(i)
      enddo
      do i=1,Np
        read(iscl1,'(E12.5)') vy(i)
      enddo
      do i=1,Np
        read(iscl1,'(E12.5)') vz(i)
      enddo
      endif

      if(Nl.gt.0)then
      read(iscl1,*) dump !write(iscl1,'(A)') 'bar2'
      do i=Np+1,Nl+Np
        read(iscl1,'(E12.5)') vx(i)
      enddo
      do i=Np+1,Nl+Np
        read(iscl1,'(E12.5)') vy(i)
      enddo
      do i=Np+1,Nl+Np
        read(iscl1,'(E12.5)') vz(i)
      enddo
      endif

      if(Nt.gt.0)then
      read(iscl1,*) dump !write(iscl1,'(A)') 'tria3'
      do i=Np+Nl+1,Nl+Np+Nt
        read(iscl1,'(E12.5)') vx(i)
      enddo
      do i=Np+Nl+1,Nl+Np+Nt
        read(iscl1,'(E12.5)') vy(i)
      enddo
      do i=Np+Nl+1,Nl+Np+Nt
        read(iscl1,'(E12.5)') vz(i)
      enddo
      endif

      if(Nq.gt.0)then
      read(iscl1,*) dump !write(iscl1,'(A)') 'quad4'
      do i=Np+Nl+Nt+1,Nl+Np+Nt+Nq
        read(iscl1,'(E12.5)') vx(i)
      enddo
      do i=Np+Nl+Nt+1,Nl+Np+Nt+Nq
        read(iscl1,'(E12.5)') vy(i)
      enddo
      do i=Np+Nl+Nt+1,Nl+Np+Nt+Nq
        read(iscl1,'(E12.5)') vz(i)
      enddo
      endif


!      write(iscl1,'(A)') 'coordinates'
!      do i=1,nnos
!        write(iscl1,'(E12.5)') vx(i)
!      enddo
!      do i=1,nnos
!        write(iscl1,'(E12.5)') vy(i)
!      enddo
!      do i=1,nnos
!        write(iscl1,'(E12.5)') vz(i)
!      enddo
!
      close(iscl1)

      !..............................................................................

endsubroutine OpenVect



!!!===================== CONVERTIR UNE TABLE EN MAILLAGE !========================

subroutine create_mesh_from_tab(N_T,N_nu,nu_sat,T,nu,P,E,Cv,Cp,c,g,w,Z,sig,lambda,Seebeck,M,imode)

   implicit none
   integer, intent(in) :: N_T, N_nu
   real(PR), intent(in) :: nu_sat(1:2,1:N_T)
   real(PR), intent(in) :: T(1:N_T)
   real(PR), intent(in) :: nu(1:N_nu)
   real(PR), intent(in) :: P(1:N_nu,1:N_T)
   real(PR), intent(in) :: E(1:N_nu,1:N_T)
   real(PR), intent(in) :: Cv(1:N_nu,1:N_T)
   real(PR), intent(in) :: Cp(1:N_nu,1:N_T)
   real(PR), intent(in) :: c(1:N_nu,1:N_T)
   real(PR), intent(in) :: g(1:N_nu,1:N_T)
   real(PR), intent(in) :: w(1:N_nu,1:N_T)
   real(PR), intent(in) :: Z(1:N_nu,1:N_T)
   real(PR), intent(in) :: sig(1:N_nu,1:N_T)
   real(PR), intent(in) :: lambda(1:N_nu,1:N_T)
   real(PR), intent(in) :: Seebeck(1:N_nu,1:N_T)
   integer, intent(in) :: imode
   type(mesh), intent(out) :: M
   real(PR) :: numin, numax, Tmin,Tmax,Pmin,Pmax,Sigmin,Sigmax
   integer :: inod, iel, it, inu, iSE,iSW,iNE,iNW

   M%Nnod  = 0.0_PR 
   M%Nel   = 0.0_PR
   M%Np    = 0.0_PR
   M%Nl    = 0.0_PR
   M%Ntri  = 0.0_PR
   M%Nquad = 0.0_PR
   M%Ntetra= 0.0_PR
   M%Nhexa = 0.0_PR
   M%Npenta= 0.0_PR
   M%Npyra = 0.0_PR
 


   M%Nnod=(N_T)*(N_nu)

   M%Nel=(N_T-1)*(N_nu-1)
   M%Nquad=M%Nel

   allocate(M%nods(1:M%Nnod))
   allocate(M%el(1:M%Nel))
   allocate(M%T(1:M%Nel))
   allocate(M%nu(1:M%Nel))
   allocate(M%P(1:M%Nel))
   allocate(M%E(1:M%Nel))
   allocate(M%Cv(1:M%Nel))
   allocate(M%Cp(1:M%Nel))
   allocate(M%c(1:M%Nel))
   allocate(M%g(1:M%Nel))
   allocate(M%w(1:M%Nel))
   allocate(M%Z(1:M%Nel))
   allocate(M%sig(1:M%Nel))
   allocate(M%lambda(1:M%Nel))
   allocate(M%Seebeck(1:M%Nel))

   iel=0

   numin=log(minval(nu(:)))
   numax=log(maxval(nu(:)))
   Tmin=minval(T(:))
   Tmax=maxval(T(:))
   !Pmin=log(minval(P(:,:)))
   !Pmax=log(maxval(P(:,:)))  
   Pmin=minval(P(:,:))
   Pmax=1.0e3_PR !!! Al: 0.5e3_PR
   Sigmin=log10(1.0_PR)  !!!minval(P(:,:))
   Sigmax=log10(maxval(sig(:,:)))


   do it=1,N_T
    
     do inu=1,N_nu

         inod=(it-1)*N_nu+inu

         M%nods(inod)%id=inod
         M%nods(inod)%x=(log(nu(inu))-numin)/(numax-numin)
         M%nods(inod)%y=(T(it)-Tmin)/(Tmax-Tmin)
         !M%nods(inod)%z=(log(P(inu,it))-Pmin)/(Pmax-Pmin)
         if(imode.eq.1)then
            M%nods(inod)%z=(min(P(inu,it),Pmax)-Pmin)/(Pmax-Pmin)
         else
            M%nods(inod)%z=(min(max(log10(sig(inu,it)),sigmin),sigmax)-sigmin)/(sigmax-sigmin)
         endif

         if(it.gt.1.and.inu.gt.1)then
 
            iel=iel+1
        
            iSW=(it-2)*N_nu+(inu-1)
            iSE=(it-2)*N_nu+inu
            iNW=(it-1)*N_nu+(inu-1)
            iNE=(it-1)*N_nu+inu

            if(nu(inu-1).lt.nu_sat(1,it-1).and.&
               nu(inu).gt.nu_sat(1,it-1))then
               iSE=iSW
            endif
            if(nu(inu-1).lt.nu_sat(1,it).and.&
                   nu(inu).gt.nu_sat(1,it))then
               iNW=iNE
            endif
            if(nu(inu-1).lt.nu_sat(2,it-1).and.&
               nu(inu).gt.nu_sat(2,it-1))then
               iSW=iSE
            endif
            if(nu(inu-1).lt.nu_sat(2,it).and.&
               nu(inu).gt.nu_sat(2,it))then
               iNE=iNW
            endif

            M%el(iel)%id=iel  
            M%el(iel)%Nnod=4
            M%el(iel)%typ=3
            M%el(iel)%n(1)=iSW
            M%el(iel)%n(2)=iSE
            M%el(iel)%n(3)=iNE
            M%el(iel)%n(4)=iNW
           
            M%T(iel)=0.5_PR*(T(it-1)+T(it))
            M%nu(iel)=0.5_PR*(nu(inu-1)+nu(inu))
            M%P(iel)=0.25_PR*(P(inu-1,it-1)+P(inu,it-1)+P(inu-1,it)+P(inu,it))
            M%Cv(iel)=0.25_PR*(Cv(inu-1,it-1)+Cv(inu,it-1)+Cv(inu-1,it)+Cv(inu,it))
            M%Cp(iel)=0.25_PR*(Cp(inu-1,it-1)+Cp(inu,it-1)+Cp(inu-1,it)+Cp(inu,it))
            M%c(iel)=0.25_PR*(c(inu-1,it-1)+c(inu,it-1)+c(inu-1,it)+c(inu,it))
            M%E(iel)=0.25_PR*(E(inu-1,it-1)+E(inu,it-1)+E(inu-1,it)+E(inu,it))
            M%g(iel)=0.25_PR*(g(inu-1,it-1)+g(inu,it-1)+g(inu-1,it)+g(inu,it))
            M%w(iel)=0.25_PR*(w(inu-1,it-1)+w(inu,it-1)+w(inu-1,it)+w(inu,it))
            M%Z(iel)=0.25_PR*(Z(inu-1,it-1)+Z(inu,it-1)+Z(inu-1,it)+Z(inu,it))
            M%sig(iel)=0.25_PR*(sig(inu-1,it-1)+sig(inu,it-1)+sig(inu-1,it)+sig(inu,it))
            M%lambda(iel)=0.25_PR*(lambda(inu-1,it-1)+lambda(inu,it-1)+lambda(inu-1,it)+lambda(inu,it))
            M%Seebeck(iel)=0.25_PR*(Seebeck(inu-1,it-1)+Seebeck(inu,it-1)+Seebeck(inu-1,it)+Seebeck(inu,it))

         endif

     enddo

   enddo

end subroutine create_mesh_from_tab


end module mod_ensight


