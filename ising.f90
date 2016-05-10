
module configuration
save
integer(4)Nx,Ny,nh,which_lattice,nn
integer(4), allocatable :: ivic(:,:),ivict(:,:),ixpiy(:) ! nearest neighbors
integer(4), allocatable :: spin(:) ! spin (-1 or 1 ) or (boson 0 or 1 respectively) configuration
real(8) temperature,beta,Jh
end module configuration

module measurement
save
integer(4)lm
real(8), allocatable :: Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),mea(:)
real(8)E
end module measurement


program ising
use configuration
implicit none
integer(4) i,j,k,N,iseedd,nther,nrun,msteps
real(8) :: drand1

write(6,*)'Nx,Ny,which,temperature,Jh'
read(5,*)Nx,Ny,which_lattice,temperature,Jh
write(6,*)'iseedd,nther,nrun,msteps'
read(5,*)iseedd,nther,nrun,msteps

nh=Nx*Ny
beta=1.0d0/temperature

if(which_lattice==1)then
 nn=4
 allocate(ivic(nh,nn)) !square
 allocate(ixpiy(nh))
 call square()  
elseif(which_lattice==2)then
 nn=6
 allocate(ivic(nh,nn),ivict(nh,4)) !triangular
 call triangular()  
end if

!initializations
call rand_init(iseedd)
call initial()
!files to write measurements and configurations
open(15,file='measurement.dat',status='unknown',form='formatted')
open(20,file='configuration.dat',status='unknown',form='formatted')


!# thermalization
do i = 1,nther
   call runbin(msteps)
end do

write(6,*)"Thermalization ready..."

!# run
do i = 1,nrun
   call runbin(msteps)
   call accmeasure(i)
end do


close(15)
close(20)

end program ising

subroutine square()
use configuration
implicit none
integer(4)i,j,k,ii
ii=1

do j=1,Ny
   do i=1,Nx
   
     !#right  (1)
       if (i==Nx)then
          ivic(ii,1)=ii-(Nx-1)
       else
          ivic(ii,1)=ii+1
       end if
      ! # up (2)
       if (j==Ny)then
          ivic(ii,2)=ii-Ny*(Nx-1)
       else
          ivic(ii,2)=ii+Nx
       end if
       !#left (3)
       if (i==1)then
          ivic(ii,3)=ii+(Nx-1)
       else
        ivic(ii,3)=ii-1
       end if
       !#down (4)
       if (j==1)then
          ivic(ii,4)=ii+(Ny-1)*(Nx)
       else
          ivic(ii,4)=ii-Nx
       end if


       ixpiy(ii)=i+j
  
       ii=ii+1      
   

   end do  
end do

!do i=1,nh
!write(6,*)i,ivic(i,:)
!end do

end subroutine square

subroutine triangular()
use configuration
implicit none
integer(4)i,j,k,ii
ii=1

do j=1,Ny
   do i=1,Nx
   
     !#right  (1)
       if (i==Nx)then
          ivict(ii,1)=ii-(Nx-1)
       else
          ivict(ii,1)=ii+1
       end if
      ! # up (2)
       if (j==Ny)then
          ivict(ii,2)=ii-Ny*(Nx-1)
       else
          ivict(ii,2)=ii+Nx
       end if
       !#left (3)
       if (i==1)then
          ivict(ii,3)=ii+(Nx-1)
       else
          ivict(ii,3)=ii-1
       end if
       !#down (4)
       if (j==1)then
          ivict(ii,4)=ii+(Ny-1)*(Nx)
       else
          ivict(ii,4)=ii-Nx
       end if
       ii=ii+1       

   end do  
end do

ivic=0

ii=1
do j=1,Ny
    do i=1,Nx

        ivic(ii,1)=ivict(ii,1)
        ivic(ii,3)=ivict(ii,2)
        ivic(ii,4)=ivict(ii,3)
        ivic(ii,6)=ivict(ii,4)
        ivic(ii,2)=ivict(ivict(ii,1),2)
        ivic(ii,5)=ivict(ivict(ii,3),4)
        ii=ii+1
    end do
end do

!do i=1,nh
!write(6,*)i,ivic(i,:)
!end do

end subroutine triangular


!------------------------
!initial configuration
!-----------------------
subroutine initial()

use configuration
use measurement
implicit none
integer :: i,j 
real(8) :: drand1,bu

allocate(spin(nh))

do i=1,nh
 ! random initialization
  spin(i)=2*int(2.*drand1())-1
end do

beta=1/temperature

!#initial energy
E=0.0
do i=1,nh
    bu=0.0
    do j=1,nn/2
        bu=bu+spin(ivic(i,j))
    end do
    E=E-bu*Jh*spin(i)
  
end do

!#iinitialize arryas to accumulate measurements
lm=2
allocate(mea(lm),Vm(lm),Vm2(lm),Vmtemp(lm),Vm2temp(lm))
Vm=0.0d0
Vm2=0.0d0
mea=0.0d0

!write(6,*)'initial energy', E/nh
end subroutine initial


subroutine runbin(msteps)
use configuration
use measurement
implicit none
integer(4)i,j,k,kk,msteps
real(8) DE,Eest,r,drand1,mag,magt

mea=0.0
Eest=0.0
mag=0.0
do i=1,msteps
    do j=1,nh
        k=int(drand1()*nh+1) 
        DE=0
        do kk=1,nn
            DE=DE+spin(ivic(k,kk))
        end do
        DE=2*DE*Jh*spin(k)
        
        !#metropolis
        if (DE<0)then
           spin(k)=-spin(k)
           E=E+DE
        else
           r=drand1()   
           if (r<exp(-beta*DE))then
             spin(k)=-spin(k)
             E=E+DE
           end if
        end if
    end do

    if(which_lattice==1.and.Jh>0.0)then
      mag=mag+abs(sum(dble(spin)))/(Nx*Ny)
    else if(which_lattice==1.and.Jh<0.0)then
      magt=0.0d0 
      do k=1,nh
       magt=magt+dble(spin(k))*(-1.0d0)**ixpiy(k)
       !write(6,*)'checking',spin(k),ixpiy(k),(-1.0d0)**ixpiy(k),dble(spin(k))*(-1.0d0)**ixpiy(k)
      end do
      mag=mag+abs(magt)/(Nx*Ny)
    else if(which_lattice==2)then
        mag=mag+abs(sum(dble(spin)))/(Nx*Ny)
    end if 

    Eest=Eest+E/(Nx*Ny)
end do
mea(1)=mag/msteps
mea(2)=Eest/msteps

!write(6,*)"mea",mea

end subroutine runbin


subroutine accmeasure(i)
use configuration
use measurement
implicit none
integer(4)i,j,k
real(8)magt

Vm=Vm+mea
Vm2=Vm2+mea**2




!write(15,*)mea(1),mea(2)

if(which_lattice==1.and.Jh>0.0)then
      magt=abs(sum(dble(spin)))/(Nx*Ny)
else if(which_lattice==1.and.Jh<0.0)then
      magt=0.0d0
      do k=1,nh
       magt=magt+dble(spin(k))*(-1.0d0)**ixpiy(k)
       !write(6,*)'checking',spin(k),ixpiy(k),(-1.0d0)**ixpiy(k),dble(spin(k))*(-1.0d0)**ixpiy(k)
      end do
      magt=abs(magt)/(Nx*Ny)
else if(which_lattice==2)then
        magt=abs(sum(dble(spin)))/(Nx*Ny)
end if


write(15,*)E/(Nx*Ny),magt



!#Current estimates
Vmtemp=Vm/i
Vm2temp=Vm2/i
Vm2temp=sqrt(abs(Vm2temp-Vmtemp**2  )/dble(i))

open(10,file='results.txt',status='replace')
write(10,*)' Number of bins completed : ',i
write(10,*)' ========================================='
write(10,10)'  M         : ',Vmtemp(1),Vm2temp(1)
write(10,10)'  E         : ',Vmtemp(2),Vm2temp(2)  
write(10,*)' ========================================='
10 format(1x,a,2f14.8)
close(10)


write(20,*)spin
!write(6,*)sum(dble(spin)/nh)

end subroutine accmeasure
