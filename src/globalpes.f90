module utility
    implicit none
    public :: calc_Jacobi_from_dist, calc_weights, switchfac
    contains

    subroutine calc_Jacobi_from_dist(rcor,jacCoord,mass)
         !-------------------------------------------------------------------------------------------------------
         !   ****
         !     rcor(1) = 1-2 distance in bohr (diatom system) [r12]
         !     rcor(2) = 2-3 distance in bohr [r23]
         !     rcor(3) = 3-1 distance in bohr [r13]
         !       
         !  done by Otto
         !  jacobi coordinates is 
         !     "Jacobi coordinates: ", "    r    ", "    R    ", "    alpha    "
         !     "(surface 1)       : ", jacCoord(1,1), jacCoord(1,2), jacCoord(1,3)
         !     "(surface 2)       : ", jacCoord(2,1), jacCoord(2,2), jacCoord(2,3)
         !     "(surface 3)       : ", jacCoord(3,1), jacCoord(3,2), jacCoord(3,3)
         !
         !   mass(1) correspond to C
         !-------------------------------------------------------------------------------------------------------
         !-------------------------------------------------------------------------------------------------------
         implicit none
         real*8, dimension(3), intent(in) :: mass, rcor
         real*8, dimension(3,3), intent(out) :: jacCoord
         integer :: i
         real*8 :: pis,factor,tem
         if((rcor(1)+rcor(2)+0.00001d0).lt.rcor(3).or.(rcor(1)+rcor(3)+0.00001d0).lt.rcor(2).or.  &
           (rcor(3)+rcor(2)+0.00001d0).lt.rcor(1))then
           print*,'no se cumple la desigualdad triangular rcor(1),rcor(2),rcor(3)',rcor(1),rcor(2),rcor(3)
           stop
         endif
         
         !show  singularity ocurrence r12=0 or r13=0 or r23=0 with result NaN in global energy
         if(rcor(1)==0d0) then
           print*,"singularity r12"
         end if
         if(rcor(2)==0d0) then
           print*,"singularity r13"
         end if
         if(rcor(3)==0d0) then
           print*,"singularity r23"
         end if

         
         pis = acos(-1.d0)
         factor=180.d0/pis
         !    
         !   for pes 1**************************************
         jacCoord(1,1) = rcor(1)
         
         jacCoord(1,2) = sqrt( mass(1)*(rcor(3))**2/(mass(1)+mass(2)) + mass(2)*(rcor(2))**2/(mass(1)+mass(2)) &
                      - mass(1)*mass(2)*rcor(1)*rcor(1)/((mass(1)+mass(2))**2))
         
         tem = (mass(1)+mass(2))*jacCoord(1,2)/(2.d0*mass(1)*rcor(1)) &
                         + mass(1)*rcor(1)/(2.d0*(mass(1)+mass(2))*jacCoord(1,2)) &
                         - (mass(1)+mass(2))*(rcor(2))**2/( 2.d0*mass(1)*rcor(1)*jacCoord(1,2) )
         
         if(tem.gt.1.d0.and.tem.lt.1.0000001d0)then
            tem=1.d0   
         end if
         if(tem.lt.(-1.d0).and.tem.gt.(-1.0000001d0))then
            tem=-1.d0   
         end if
         
         jacCoord(1,3) =  factor*acos(tem)                                      
         
         !   print*,'SEP1', tem
         
         !   for pes 2**********************************
         jacCoord(2,1) = rcor(3)
         
         jacCoord(2,2) = sqrt( mass(1)*(rcor(1))**2/(mass(1)+mass(3)) + mass(3)*(rcor(2))**2/(mass(1)+mass(3)) &
                      - mass(1)*mass(3)*rcor(3)*rcor(3)/((mass(1)+mass(3))**2))
         
         tem = (mass(1)+mass(3))*jacCoord(2,2)/(2.d0*mass(1)*rcor(3)) &
                         + mass(1)*rcor(3)/(2.d0*(mass(1)+mass(3))*jacCoord(2,2)) &
                         - (mass(1)+mass(3))*(rcor(2))**2/( 2.d0*mass(1)*rcor(3)*jacCoord(2,2) )
         
         if(tem.gt.1.d0.and.tem.lt.1.0000001d0)then
            tem=1.d0   
         end if
         if(tem.lt.(-1.d0).and.tem.gt.(-1.0000001d0))then
            tem=-1.d0   
         end if
         
         jacCoord(2,3) =  factor*acos(tem) 
         !   print*,'SEP2',tem
         !   for pes 3*******************************************
         jacCoord(3,1) = rcor(2)
         
         jacCoord(3,2) = sqrt( mass(2)*(rcor(1))**2/(mass(3)+mass(2)) + mass(3)*(rcor(3))**2/(mass(3)+mass(2)) &
                      - mass(3)*mass(2)*rcor(2)*rcor(2)/((mass(3)+mass(2))**2))
         
         tem = (mass(3)+mass(2))*jacCoord(3,2)/(2.d0*mass(3)*rcor(2)) &
                         + mass(3)*rcor(2)/(2.d0*(mass(3)+mass(2))*jacCoord(3,2)) &
                         - (mass(3)+mass(2))*(rcor(1))**2/( 2.d0*mass(3)*rcor(2)*jacCoord(3,2) )
         
         if(tem.gt.1.d0.and.tem.lt.1.0000001d0)then
            tem=1.d0   
         end if
         if(tem.lt.(-1.d0).and.tem.gt.(-1.0000001d0))then
            tem=-1.d0   
         end if
         
         jacCoord(3,3) =  factor*acos( tem ) 
         
         RETURN
    end subroutine calc_Jacobi_from_dist

    subroutine calc_weights(w,r)
        implicit none
        real*8, intent(in)  :: r(3)
        real*8, intent(out) :: w(3)
        real*8              :: w0(3),sumw0
        integer             :: i

        real*8, dimension(3), save :: w_save = (/1.d0,0.d0,0.d0/) 
        do i=1,3
            call switchfac(r(i),w0(i))
        end do
        if((w0(1) == 0.d0).and.(w0(2) == 0.d0).and.(w0(3) == 0.d0)) then
            w = w_save
            sumw0 = 4.0d-306
        else
            sumw0 = sum(w0)
             w = w0/sumw0
        end if
    
        !store old weights
        w_save = w

    end subroutine calc_weights

    subroutine switchfac(x,f,dfdx)
        !This is the subroutine to calculate the switching factors
        implicit none

        !input
        real*8, intent(in) :: x

        !output
        real*8, intent(out) :: f
        real*8, intent(out), optional :: dfdx
        
        !local
        !####CHANGE CUTOFF FOR SWITCHING FUNCTION HERE####
        real*8, parameter :: dx = 0.8d0 !determines the "cutoff" of the switching 
        real*8, parameter :: a  = 4.0d0 !determines plateau of switching function 
        
        !HERE THE SWITCHING FUNCTION IS SPECIFIED. 
        
        f= dexp(-(x/dx)**a)
        
        if (present(dfdx))  &
             dfdx= -(a/x)*f*(x/dx)**a
         

    end subroutine switchfac


end module
!----------------end modules------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BEGIN OF THE MAIN PROGRAM!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program globalpes
    use myRKHS
    use utility

    implicit none
    integer                     :: nalpha, nr, nvdw, ntot, pesidx

    integer                     :: i,j,k,l,rc,indexi
    integer                     :: ios, pl1_n, pl2_n, pl3_n
    real*8                      :: alpha,r,vdw,a
    real*8                      :: disv(3)! (r12,r23,r13)
    real*8                      :: n(2), E(3)
    real*8                      :: compl,f,EE,minn
    real*8                      :: mass(3)
    real*8,allocatable,save     :: dataArray1(:,:), dataArray2(:,:), dataArray3(:,:)
    real*8,allocatable          :: r_grid1(:),vdw_grid1(:),alpha_grid1(:), ene_grid1(:)
    real*8,allocatable          :: r_grid2(:),vdw_grid2(:),alpha_grid2(:), ene_grid2(:)
    real*8,allocatable          :: r_grid3(:),vdw_grid3(:),alpha_grid3(:), ene_grid3(:)
    real*8,allocatable          :: asymparray(:) 
    real*8,allocatable,save     :: Ql1_coef(:),Qr1_coef(:)
    real*8,allocatable,save     :: Ql2_coef(:),Qr2_coef(:)
    real*8,allocatable,save     :: Ql3_coef(:),Qr3_coef(:)

    real*8,allocatable          :: gridvalues(:,:)
    real*8,allocatable          :: currvalues(:)
    integer,allocatable         :: dimcount(:)
    character(len=1)            :: dummy
    character(len=30)           :: fname

    ! physical constants, converting factors
    real*8, parameter :: hartree2kcalpmol =   627.510d0    !kcal/mol to hartree
    real*8, parameter :: hartree2kjoulepmol =   2625.500d0    !kcal/mol to hartree
    real*8, parameter :: hartree2eV       =   27.211396132d0
    real*8            :: eV2Kcalpmol     = hartree2kcalpmol/hartree2eV
    real*8            :: eV2Kjoulepmol     = hartree2kjoulepmol/hartree2eV


!determine the total number of gridpoints
open(20,file="data/data1.txt", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/data1.txt was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, GRID HAVE NOT BEEN CALCULATED."
        stop
end if
i = 0
do while(.true.)
    read(20,*,iostat = ios) dummy
    if(ios /= 0) exit
    i = i + 1
end do
rewind(20)
ntot = i
write(*,*) "#",ntot
!allocation of grid...define the dimension inluding the E value
allocate(gridvalues(ntot,4),currvalues(4),dimcount(4))
dimcount=1
read(20,*,iostat=ios) gridvalues(1,:)

!read the gridvalues
do i = 2,ntot
    read(20,*,iostat=ios) currvalues(:)
    if(ios /= 0) write(*,*) "File data/data1.txt could not be read properly."
    do j = 1,3
        if(currvalues(j) > gridvalues(dimcount(j),j)) then
            dimcount(j) = dimcount(j) + 1
            gridvalues(dimcount(j),j) = currvalues(j)
        end if
    end do
end do    
close(20)
write(*,*) "#",dimcount

nalpha=dimcount(1)
nr=dimcount(2)
nvdw=dimcount(3)


!determine the total number of wpoints to determinate npl
open(20,file="data/w_tot1.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_tot1.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
i = 0
do while(.true.)
    read(20,*,iostat = ios) dummy
    if(ios /= 0) exit
    i = i + 1
end do
close(20)
pl1_n = i/(nr*nvdw)-1
write(*,*) "#Pl1_n=",pl1_n+1


allocate(dataArray1(ntot,4), dataArray2(ntot,4), dataArray3(ntot,4))
allocate(r_grid1(nr),vdw_grid1(nvdw),alpha_grid1(nalpha), ene_grid1(ntot))
allocate(r_grid2(nr),vdw_grid2(nvdw),alpha_grid2(nalpha), ene_grid2(ntot))
allocate(r_grid3(nr),vdw_grid3(nvdw),alpha_grid3(nalpha), ene_grid3(ntot))
allocate(asymparray(nr))
allocate(Ql1_coef(ntot),Qr1_coef(nr))
allocate(Ql2_coef(ntot),Qr2_coef(nr))
allocate(Ql3_coef(ntot),Qr3_coef(nr))


!------------------ read datas1----------
open(20,file="data/data1.txt", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/data1.txt was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, GRID HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,ntot
  read(20,*) dataarray1(i,:)
end do
close(20)

open(20,file="data/w_tot1.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_tot1.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,(pl1_n+1)*nr*nvdw
  read(20,*) Ql1_coef(i)
end do
close(20)

open(20,file="data/w_r1.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_r1.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,nr
  read(20,*) Qr1_coef(i)
end do
close(20)

!read VdW
VdW_grid1=dataArray1(1:nvdw,3)
!read r
do i=1,nR
    indexi = 1+(i-1)*nVdW
    r_grid1(i)=dataArray1(indexi,2)
end do
!read aplha
do i=1,nAlpha
    indexi = 1+(i-1)*nVdW*nR
    alpha_grid1(i)=dataArray1(indexi,1)
end do


!determine the total number of wpoints to determinate npl2
open(20,file="data/w_tot2.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_tot2.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
i = 0
do while(.true.)
    read(20,*,iostat = ios) dummy
    if(ios /= 0) exit
    i = i + 1
end do
close(20)
pl2_n = i/(nr*nvdw)-1
write(*,*) "#Pl2_n=",pl2_n+1


!------------------ read datas2----------

open(20,file="data/data2.txt", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/data2.txt was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, GRID HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,ntot
  read(20,*) dataarray2(i,:)
end do
close(20)

open(20,file="data/w_tot2.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_tot2.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,(pl2_n+1)*nr*nvdw
  read(20,*) Ql2_coef(i)
end do
close(20)

open(20,file="data/w_r2.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_r2.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,nr
  read(20,*) Qr2_coef(i)
end do
close(20)

!read VdW
VdW_grid2=dataArray2(1:nvdw,3)
!read r
do i=1,nR
    indexi = 1+(i-1)*nVdW

    r_grid2(i)=dataArray2(indexi,2)
end do
!read aplha
do i=1,nAlpha
    indexi = 1+(i-1)*nVdW*nR
    alpha_grid2(i)=dataArray2(indexi,1)
end do



!determine the total number of wpoints to determinate npl
open(20,file="data/w_tot3.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_tot3.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
i = 0
do while(.true.)
    read(20,*,iostat = ios) dummy
    if(ios /= 0) exit
    i = i + 1
end do
close(20)
pl3_n = i/(nr*nvdw)-1
write(*,*) "#Pl3_n (only even)=",pl3_n+1

!------------------ read datas3----------

open(20,file="data/data3.txt", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/data2.txt was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, GRID HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,ntot
  read(20,*) dataarray3(i,:)
end do
close(20)

open(20,file="data/w_tot3.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_tot2.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,(pl3_n+1)*nr*nvdw
  read(20,*) Ql3_coef(i)
end do
close(20)

open(20,file="data/w_r3.dat", status='old',action='read',iostat=ios)
if(IOS /= 0) then
        print*, "WARNING: data/w_r3.dat was not found! Please check the filename."
        print*, "PROGRAM TERMINATED, WEIGHTS HAVE NOT BEEN CALCULATED."
        stop
end if
do i=1,nr
  read(20,*) Qr3_coef(i)
end do
close(20)

!read VdW
VdW_grid3=dataArray3(1:nvdw,3)
!read r
do i=1,nR
    indexi = 1+(i-1)*nVdW

    r_grid3(i)=dataArray3(indexi,2)
end do
!read aplha
do i=1,nAlpha
    indexi = 1+(i-1)*nVdW*nR
    alpha_grid3(i)=dataArray3(indexi,1)
end do

!-------------------end read data--------------



!change on the mas order define the initial jacobi coordinates
mass(1) = 12.0096d0  !this corresponds to atom 1 (C)
mass(2) = 14.00643d0 !this corresponds to atom 2 (N)
mass(3) = 14.00643d0 !this corresponds to atom 3 (N)


!TESTING
do i=0,1000
  disv(1)=1.17d0 !1-2 DISTANCE
  disv(3)=.7d0+0.01d0*i !2-3 DISTANCE
  disv(2)=sqrt(disv(1)**2+disv(3)**2) !1-3 DISTANCE
  a=Etot(mass,disv)
  write(*,"(F16.8,F24.16)") disv(3),a
end do


!----------------functions--------------------------
contains
    real*8 recursive function Epot(jacv,surf) result(Ep)
        implicit none
        real*8, intent(in)  :: jacv(3)
        real*8              :: r,alpha,vdw! 
        integer,intent(in)  :: surf
        integer             :: j,l,pl_n, bnd=0
        real*8              :: y
        real*8              :: a,b,delta,e0,e1,aux
        real*8              :: dataa(ntot,4)
        real*8,allocatable  :: ql(:)

        delta=0.05
        Ep=0d0
        e0=0d0
        e1=0d0
        a=0d0
        b=0d0
        alpha=jacv(3)
        r=jacv(1)
        vdw=jacv(2)
        y=dcosd(alpha)

        if(surf==1) then
          pl_n=pl1_n
          allocate(ql((pl_n+1)*nr*nvdw))
          ql=Ql1_coef
          dataa=dataArray1
        else if(surf==2) then
          pl_n=pl2_n
          allocate(ql((pl_n+1)*nr*nvdw))
          ql=Ql2_coef
          dataa=dataArray2          
        else if(surf==3) then
          pl_n=pl3_n
          allocate(ql((pl_n+1)*nr*nvdw))
          ql=Ql3_coef
          dataa=dataArray3
        else 
          print *,"error Epot"
        end if
        !1/r^2 extrapolation in repulsive region
        if(1.eq.1) then
            if(r<dataa(1,2) .or. vdW<dataa(1,3)) then
                if(r<dataa(1,2)) then
                  e0=Epot((/dataa(1,2),vdw,alpha/),surf)
                  e1=Epot((/dataa(1,2)+delta,vdw,alpha/),surf)
                  if(e0<e1) then
                    aux=e0
                    e0=e1
                    e1=aux
                  end if
                  a=(e0-e1)/( dataa(1,2)**(-2) - (dataa(1,2)+delta)**(-2))
                  b=e1-a*(dataa(1,2)+delta)**(-2)
                  Ep=a*(r**(-2))+b
                else
                  e0=Epot((/r,dataa(1,3),alpha/),surf)
                  e1=Epot((/r,dataa(1,3)+delta,alpha/),surf)
                  if(e0<e1) then
                    aux=e0
                    e0=e1
                    e1=aux
                  end if
                  a=(e0-e1)/( dataa(1,3)**(-2) - (dataa(1,3)+delta)**(-2) )
                  b=e1-a*(dataa(1,3)+delta)**(-2)
                  Ep=a*(vdw**(-2))+b
                end if
                bnd=1
            else
            !lagendre-kernel interpolation in remain regions
                if(surf==3) then
                  do l=0, pl_n
                      do j=1, nr*nvdw
                          Ep =Ep + ql(j+l*nR*nVdW)*q_ker(r,dataa(j,2))*q_ker(vdw,dataa(j,3))*pl(y,l*2)
                      end do
                  end do
                else
                  do l=0, pl_n
                      do j=1, nr*nvdw
                          Ep =Ep + ql(j+l*nR*nVdW)*q_ker(r,dataa(j,2))*q_ker(vdw,dataa(j,3))*pl(y,l)
                      end do
                  end do
                end if
            end if
        end if
        

        if(bnd==0) then
            Ep=Ep+asymp(r,surf)
        else
          bnd=0
        end if
        
    end function Epot


    real*8 function Etot(mass,disv)
        implicit none
        real*8, intent(in)   :: mass(3),disv(3)! r12,r23,r13
        integer             :: i
        real*8              :: w(3)
        real*8              :: jacm(3,3)! matrix with jaccoord generated from disv
        !jacm: r, R, alpha
        Etot=0d0
        
        call calc_Jacobi_from_dist(disv,jacm,mass)

        call calc_weights(w,jacm(:,1))
        
        do i=1,3
          E(i)=Epot(jacm(i,:),i)
        end do
        
        Etot=w(1)*E(1)+w(2)*E(2)+w(3)*E(3)
        
    end function Etot

    real*8 function asymp(r,surf)
        implicit none
        real*8, intent(in)  :: r
        integer,intent(in)  :: surf
        integer             :: i
        real*8             :: qr(nr),dataa(ntot,4)
        asymp = 0d0
        
        if(surf==1) then
          qr=Qr1_coef
          dataa=dataArray1
        else if(surf==2) then
          qr=Qr2_coef
          dataa=dataArray2
        else if(surf==3) then
          qr=Qr3_coef
          dataa=dataArray3
        else
          print *,"error:asymp"
          stop
        end if

        do i=1,nr
          asymp = asymp + qr(i)*q_ker(r,dataa(1+(i-1)*nvdw,2))
        end do
    end function asymp 

    ! *********************************************************************
    function pl(x,n)
    ! Legendre Polynomial
    ! download from
    ! http://ww2.odu.edu/~agodunov/computing/programs/f90/Legendre.f90
    ! INPUT:
    !      n: the number of Legendre Polynomial terms
    !      x: the coordinate (-1<x<1)
    ! OUTPUT:
    !     pl: Pl_0(x),Pl_1(x)... Pl_n-1(x)
    !         p(0)= Pl_0(x), p(1)=Pl_1(x)...
    ! calculates Legendre polynomials Pn(x)
    ! using the recurrence relation
    ! if n > 100 the function retuns 0.0

    double precision pl
    double precision x
    double precision pln(0:n)
    integer n, k

    pln(0) = 1.0d0
    pln(1) = x

    if (n <= 1) then
      pl = pln(n)
    else
      do k=1,n-1
        pln(k+1) = ((2.0d0*k+1.0d0)*x*pln(k) - float(k)*pln(k-1))/(float(k+1))
      enddo
     pl = pln(n)
    endif
    return
    END
!
! *********************************************************************


end program globalpes
