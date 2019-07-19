!subroutine preri(eri,nao,dm,singleatom,singleitem,num1,twoatom,twoitem,num2,threeatom,threeitem,num3,fouratom,fouritem,num4,listA)
!subroutine preri(eri,nao,dm,singleatom,singleitem,num1,twoatom,twoitem,num2,listA)
subroutine preri(atm,atml,bas,basl,env,envl,nao,nshls,dm,singleatom,singleitem,num1,atom_ej, atom_ek)
implicit none 
!integer findloc
integer,external :: Ainclude,body2,CINTcgto_spheric
!real*8,external :: cint2e_sph
integer i,j,k,l,tem,a,b,c,d,flag,di,dj,dk,dl
integer i0,i1,n
integer m,II,JJ,KK,LL 
integer singleitem,num1,num
integer nao
integer singleatom(singleitem,num1)
integer single(singleitem*num1)
integer basis(nao)
integer slice
!real*8 eri(1,nao,nao,nao)
real*8 dm(nao,nao)
real*8 e_j, e_k
!real*8 vjk(nao,nao)
real*8 sumE
!real*8 listA(nao,nao,nao,nao)
real*8 atom_ej(singleitem+1)
real*8 atom_ek(singleitem+1)
real*8 tem1,tem2,t1,t2
integer time1,time2
integer atml,basl,envl
integer atm(6,atml),bas(8,basl)
!integer atm1(6,atml),bas1(8,basl)
real*8 env(envl)
real*8 buf(:,:,:,:)
real*8 aoshell(:,:,:,:)
real*8 energy(singleitem+1)
real*8 schw(:,:)
integer ncf_sh(:),flcf_sh(:,:),cf2sh(:)
allocatable buf
allocatable ncf_sh
allocatable flcf_sh
allocatable cf2sh
allocatable aoshell
allocatable schw
integer shls(4)
integer nshls,nbas
integer ifi,ils,jfi,jls,kfi,kls,lfi,lls,ibas,kbas,jbas,lbas
integer natm
integer ni,nj
real*8 thresh
real*8 hf_ene
logical kl,lj
integer int2e
real*8 ddm
character*2 whatjk
!atm1 = reshape((/1, 20,  1, 23,  0,  0,1, 24,  1, 27,  0,  0/),(/6,atml/))
!bas1 = reshape((/0,  0,  3,  1,  0, 28, 31,  0,1,  0,  3,  1,  0, 28, 31,0/),(/8,basl/))

! -----------------------------------------------------------------------
!  f2py -m frame_small5 -c frame_small5.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
!  f2py -m frame_small5 -c frame_small5.f90 -L/home/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint
!  f2py -m frame_small5 -c frame_small5.f90 -L/home/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
!  f2py -m frame_small6 -c frame_small6.f90 -L/home/wsr/pyscf/lib/ -lcint --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
!  f2py -m jkeda -c jkeda.f90 -L/home/wsr/pyscf/lib/ -lcint --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
! -----------------------------------------------------------------------

!f2py intent(in) :: atm,atml,bas,basl,env,envl,nshls
!f2py intent(in) :: nao,dm,singleatom
!f2py intent(in) :: singleitem,num1
!f2py intent(in) :: whatjk

!f2py intent(out) :: atom_ej, atom_ek
!f2py depend(nao) :: dm
!f2py depend(basl) :: bas
!f2py depend(envl) :: env
!f2py depend(atml) :: atm


!vjk = 0.0d0
!listA = 0.0d0
call system_clock(time1)
t1 = time()
!e_coul = 0.0d0
e_j = 0.0d0
e_k = 0.0d0
atom_ej = 0.0d0
atom_ek = 0.0d0
energy = 0.0d0

flag = 1
do i=1,singleitem
   do j=1,num1
      single(flag) = singleatom(i,j)
      flag = flag + 1
   enddo
enddo

!write(*,*) "atm",atm

num = singleitem * num1
do i=1,nao
  tem = Ainclude(i,num,single)
  if (mod(tem,num1)==0) then
      basis(i) = tem/num1
  else 
      basis(i) = tem/num1 + 1 
  endif
enddo

sumE =0.0d0
hf_ene = 0.0d0

nbas = nao
natm = singleitem

!------------------------------------------------------------ 
! obtain the map between contracted functions and contracted shells
!------------------------------------------------------------

allocate(ncf_sh(nshls))
allocate(flcf_sh(2,nshls))
allocate(schw(nshls,nshls))

!call schwarz(schw,nshl)
!---------------------------------------------------

do i=1,nshls
   ncf_sh(i) = cintcgto_spheric(i-1,bas)
enddo

allocate(cf2sh(nbas))
i0 = 1 
do i = 1, nshls
   n = ncf_sh(i)
   cf2sh(i0:i0+n-1) = i
   flcf_sh(1,i) = i0
   flcf_sh(2,i) = i0 + n -1 
   i0 = i0 + n
enddo
!------------------------------------------------------------------
do i=1,nshls
    ni=ncf_sh(i)
    shls(1)=i-1
    shls(3)=i-1
    do j=1,i
        nj=ncf_sh(j)
        shls(2)=j-1
        shls(4)=j-1
        allocate(buf(ni,nj,ni,nj))
        call cint2e_sph(buf, shls, atm, 0, bas, 0, env, 0_8)
        schw(i,j)=maxval(buf)
        if (i/=j) schw(j,i)=schw(i,j)
        deallocate(buf)
    enddo
enddo
thresh = 1.0D-10
int2e = 0
!do i=1,nao
!   do j = 1,nao
!       ddm = dm(i,j)-dm(j,i)
!       if (ddm > 1.0d-10) then
!           write(*,*) i,j,dm(i,j),dm(j,i)
!       endif
!   enddo
!enddo
!------------------------------------------------------------------
!open(unit=16,file='c12h26-ccpvtz-debug.out')
!nshls = 62

!$OMP PARALLEL DO schedule(guided) &
!$omp default(private) &
!$omp shared(nshls,atm,bas,env,flcf_sh,nbas, &
!$omp dm,basis,nao,natm,schw,thresh,ncf_sh) reduction(+:atom_ej, atom_ek,int2e)
do i=1,nshls
  shls(1) = i - 1 
  ifi = flcf_sh(1,i)
  ils = flcf_sh(2,i)
  di = ncf_sh(i) 
  do l=1,nshls 
    shls(3) = l - 1
    lfi = flcf_sh(1,l)
    lls = flcf_sh(2,l)
    dl = ncf_sh(l)
    allocate(aoshell(ifi:ils,nbas,lfi:lls,nbas))
    aoshell = 0.0d0
    do j=1,nshls
      shls(2) = j - 1 
      jfi = flcf_sh(1,j)
      jls = flcf_sh(2,j)
      dj = ncf_sh(j)
      do k=1,nshls
        shls(4) = k - 1 
        kfi = flcf_sh(1,k)
        kls = flcf_sh(2,k)
        dk = ncf_sh(k)
      !write(*,*) di,dj,dk,dl
        if (dsqrt(schw(i,j))*dsqrt(schw(k,l))>thresh) then
            allocate(buf(di,dj,dl,dk))
            call cint2e_sph(buf, shls, atm, 0, bas, 0, env, 0_8)
            aoshell(:,jfi:jls,:,kfi:kls) = buf(:,:,:,:)
      !write(*,*) "buf\n",buf
            deallocate(buf)
        endif
      enddo
    enddo
    ! vj
    
     do lbas = lfi,lls
       d = basis(lbas)
       do k = lbas,nao
         c = basis(k)
         kl = (k/=lbas)
         tem1 = dm(lbas,k)
         do ibas = ifi,ils
           !if (ibas > kbas) then
           !  exit
           !endif
           a = basis(ibas)
           do j = ibas,nao
               !if (j > l) then
               !  exit
               !endif 
               b = basis(j)
               e_j = 0.5d0*tem1*dm(j,ibas)*aoshell(ibas,j,lbas,k)
               !write(*,*) ibas,j,kbas,l,e_j
               !if (dabs(e_j) > 1D-12) then
                   if (kl) then
                   !    !write(*,*)
                       e_j = e_j*2.0d0
                   endif
                   if (ibas/=j) then
                   !    !write(*,*) "ij",ibas,j
                       e_j = e_j*2.0d0
                   endif
                   !if ((ibas/=kbas) .or. (j/=l)) then
                   !    e_j = e_j*2.0d0
                   !endif
               int2e = int2e+1     
               !write(*,*) ibas,j,kbas,l,e_j
               !if (dabs(e_j) > 1D-12) then
               if (dabs(e_j) > 1D-12) then
               !hf_ene = hf_ene + e_coul
               !  tem = body2(a,b,c,d,natm,m,II,JJ,KK,LL)
                   II = natm+1
                   JJ = natm+1
                   KK = natm+1
                   LL = natm+1
                   m = 0
                   if (a/=b) then
                       II = a
                       JJ = b
                       m = m+2
                   else
                       II = a
                       m = m+1
                   endif
                
                   if (c /= II) then
                       if (c /= JJ) then
                           KK = c
                            m = m+1
                       endif
                   endif
                   if (d /= II) then
                       if (d /= JJ) then
                           if (d /= KK) then
                               LL = d
                               m = m+1
                           endif
                       endif
                   endif
                   e_j = e_j/m 
                   atom_ej(II) = atom_ej(II) + e_j
                   atom_ej(JJ) = atom_ej(JJ) + e_j
                   atom_ej(KK) = atom_ej(KK) + e_j
                   atom_ej(LL) = atom_ej(LL) + e_j
               !else
               !    write(*,*) ibas,j,kbas,l,e_j
               endif
          enddo
        enddo
      enddo
    
    enddo
    !stop
    !if (whatjk=='jk') then
    ! vk
    do lbas = lfi,lls
       d = basis(lbas)
       do j = 1,nao
          b = basis(j)
          lj = (lbas/=j)
          do ibas = ifi,ils
            a = basis(ibas)
            tem1 = dm(ibas,lbas)
            do k = 1,nao
               c = basis(k)
               e_k = 0.5d0*tem1*dm(j,k)*(-0.5d0)*aoshell(ibas,k,lbas,j)
               !if (dabs(e_k) > 1D-12) then
                   !if (lj) then
                   !    e_k = e_k*2.0d0
                   !endif
                   !if (ibas/=k) then
                   !    e_k = e_k*2.0d0
                   !endif
               if (dabs(e_k) > 1D-12) then
               !hf_ene = hf_ene + e_coul
               !  tem = body2(a,b,c,d,natm,m,II,JJ,KK,LL)
                   II = natm+1
                   JJ = natm+1
                   KK = natm+1
                   LL = natm+1
                   m = 0
                   if (a/=b) then
                       II = a
                       JJ = b
                       m = m+2
                   else
                       II = a
                       m = m+1
                   endif

                   if (c /= II) then
                       if (c /= JJ) then
                           KK = c
                            m = m+1
                       endif
                   endif
                   if (d /= II) then
                       if (d /= JJ) then
                           if (d /= KK) then
                               LL = d
                               m = m+1
                           endif
                       endif
                   endif                 
                   e_k = e_k/m
                   write(*,*) ibas,k,lbas,j, m, e_k
                   atom_ek(II) = atom_ek(II) + e_k
                   atom_ek(JJ) = atom_ek(JJ) + e_k
                   atom_ek(KK) = atom_ek(KK) + e_k
                   atom_ek(LL) = atom_ek(LL) + e_k
               endif
          enddo
        enddo
      enddo
    enddo
    !endif

!  !energy(:) = atom_energy(:)
!  !energy = energy + atom_energy
  deallocate(aoshell)
  enddo
enddo
!$omp end parallel do

write(*,*) "int2e", int2e


call system_clock(time2)
t2 = time()
write(*,*) "time",t2-t1,time2-time1

return 
end subroutine preri

function Ainclude(k,n,array)
  implicit none
  integer :: i
  integer :: k,n
  integer :: Ainclude
  integer :: array(n)

  do i = 1,n
     if (array(i) == k) then
         Ainclude =i
         return 
     endif
  enddo

  if (i==n+1) Ainclude =0

  return
end function

function body(a,b,c,d,m,II,JJ,KK,LL)
    implicit none
    integer :: a,b,c,d,m,II,JJ,KK,LL,n
    integer :: body
    integer :: atom(4),btom(4)
    integer :: i,j,temp
 
    if ((a==b) .and. (b==c) .and. (c==d)) then                                                                                      
            m = 1 
    else if  (((a==b) .and. (b==c)) .or. &
              ((a==b) .and. (b==d)) .or. &
              ((a==c) .and. (c==d)) .or. &
              ((b==c) .and. (c==d)) .or. &
              ((a==b) .and. (c==d)) .or. &
              ((a==c) .and. (b==d)) .or. &
              ((a==d) .and. (b==c))) then
            m = 2
    else if ((a==b) .or. (b==c) .or. (c==d) .or. (a==c) .or. (a==d) .or.(b==d)) then
            m = 3
    else 
            m = 4
    endif
    II = 0
    JJ = 0
    KK = 0 
    LL = 0
    atom(1) = a 
    atom(2) = b
    atom(3) = c
    atom(4) = d
    do i=1,3
        do j=i+1,4
            if (atom(i) .ge. atom(j)) then
                temp = atom(i)
                atom(i) = atom(j)
                atom(j) = temp
            endif
        enddo
    enddo
    btom = 0
    btom(1) = atom(1)
    n = 0
    do i=1,4
       do j=i+1,4
          if (atom(i)==atom(j)) then
             n = n + 1
             continue
          else
              btom(i+1)= atom(i+1+n)
          endif 
       enddo
    enddo

    if(m==1) then 
        II = btom(1)
    else if (m==2) then
        II = btom(1)
        JJ = btom(2)
    else if (m==3) then 
        II = btom(1)
        JJ = btom(2)
        KK = btom(3)
    else 
        II = btom(1)
        JJ = btom(2)
        KK = btom(3)
        LL = btom(4)
    endif
    body = 0
return 
end function body

function body2(a,b,c,d,natm,m,II,JJ,KK,LL)
   implicit none
   integer :: a,b,c,d,natm,m,II,JJ,KK,LL
   integer :: body2
   II = natm+1
   JJ = natm+1
   KK = natm+1
   LL = natm+1
   m = 0
   if (a/=b) then
       II = a
       JJ = b
       m = m+2
   else
       II = a
       m = m+1
   endif

   if (c /= II) then
       if (c /= JJ) then
           KK = c
            m = m+1
       endif
   endif

    if (d /= II) then
        if (d /= JJ) then
            if (d/= KK) then
                LL = d
                m = m+1
            endif
        endif
    endif
    body2 = 0
    !PRINT *,a,b,c,d,'->', II, JJ ,KK ,LL,m,natm

return
end function body2

