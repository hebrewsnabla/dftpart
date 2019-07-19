!subroutine preri(eri,nao,dm,singleatom,singleitem,num1,twoatom,twoitem,num2,threeatom,threeitem,num3,fouratom,fouritem,num4,listA)
!subroutine preri(eri,nao,dm,singleatom,singleitem,num1,twoatom,twoitem,num2,listA)
subroutine preri(atm,atml,bas,basl,env,envl,nao,nshls,dm,singleatom,singleitem,num1,atom_energy)
implicit none 
!integer findloc
integer,external :: Ainclude,body2,CINTcgto_spheric
!real*8,external :: cint2e_sph
integer i,j,k,l,tem,a,b,c,d,flag,di,dj,dk,dl
integer i0,i1,n
integer m,II,JJ,KK,LL 
integer mij,mijk,mijkl
integer singleitem,num1,num
integer nao
integer singleatom(singleitem,num1)
integer single(singleitem*num1)
integer basis(nao)
integer slice
!real*8 eri(1,nao,nao,nao)
real*8 dm(nao,nao)
real*8 e_coul
!real*8 vjk(nao,nao)
real*8 sumE
!real*8 listA(nao,nao,nao,nao)
real*8 atom_energy(singleitem+1)
real*8 tem1,tem2,t1,t2
integer time1,time2,time3,time4
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
logical keql,ieqj
!atm1 = reshape((/1, 20,  1, 23,  0,  0,1, 24,  1, 27,  0,  0/),(/6,atml/))
!bas1 = reshape((/0,  0,  3,  1,  0, 28, 31,  0,1,  0,  3,  1,  0, 28, 31,0/),(/8,basl/))

! -----------------------------------------------------------------------
!  f2py -m frame_small5 -c frame_small5.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
!  f2py -m frame_small5 -c frame_small5.f90 -L/home/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint
!  f2py -m frame_small6 -c frame_small6.f90 -L/home/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint --fcompiler=gfortran --f90flags='-fopenmp' -lgomp

! -----------------------------------------------------------------------

!f2py intent(in) :: atm,atml,bas,basl,env,envl,nshls
!f2py intent(in) :: nao,dm,singleatomi
!f2py intent(in) :: singleitem,num1

!f2py intent(out) :: atom_energy
!f2py depend(nao) :: dm
!f2py depend(basl) :: bas
!f2py depend(envl) :: env
!f2py depend(atml) :: atm

!f2py depend(singleitem) :: atom_energy

!vjk = 0.0d0
!listA = 0.0d0
call system_clock(time1)
t1 = time()
e_coul = 0.0d0
atom_energy = 0.0d0
energy = 0.0d0

flag = 1
do i=1,singleitem
   do j=1,num1
      single(flag) = singleatom(i,j)
      flag = flag + 1
   enddo
enddo

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
!------------------------------------------------------------------
!open(unit=16,file='c12h26-ccpvtz-debug.out')
!nshls = 62

!$OMP PARALLEL DO schedule(guided) &
!$omp default(private) &
!$omp shared(nshls,atm,bas,env,flcf_sh,nbas, &
!$omp dm,basis,nao,natm,schw,thresh,ncf_sh) reduction(+:atom_energy)
do i=1,nshls
  shls(4) = i - 1 
  ifi = flcf_sh(1,i)
  ils = flcf_sh(2,i)
  di = ncf_sh(i) 
  do k=1,nshls 
    shls(2) = k - 1
    kfi = flcf_sh(1,k)
    kls = flcf_sh(2,k)
    dk = ncf_sh(k)
    allocate(aoshell(nbas,kfi:kls,nbas,ifi:ils))
    call system_clock(time1)
    aoshell = 0.0d0
    do l=1,nshls
      shls(3) = l - 1 
      lfi = flcf_sh(1,l)
      lls = flcf_sh(2,l)
      dl = ncf_sh(l)
      do j=1,nshls
        shls(1) = j - 1 
        jfi = flcf_sh(1,j)
        jls = flcf_sh(2,j)
        dj = ncf_sh(j)
      !write(*,*) di,dj,dk,dl
        if (dsqrt(schw(j,k))*dsqrt(schw(i,l))>thresh) then
            allocate(buf(dj,dk,dl,di))
            call cint2e_sph(buf, shls, atm, 0, bas, 0, env, 0_8)
            aoshell(jfi:jls,:,lfi:lls,:) = buf(:,:,:,:)
            !write(*,*) "buf\n",buf
            deallocate(buf)
        endif
      enddo
    enddo
    call system_clock(time2) 
    write(*,*) "Calc eri time=",time2-time1
    do ibas = ifi,ils
       a = basis(ibas)
       do j =ibas,nao 
          b = basis(j)
          ieqj = (ibas/=j)
          tem1 = dm(j,ibas)
       !   II = natm + 1 
       !   JJ = natm + 1 
       !   mij = 0
       !   if (a/=b) then 
       !      II = a
       !      JJ = b 
       !      mij = mij + 2
       !   else 
       !      II = a 
       !      mij = mij + 1
       !   endif
          do kbas = kfi,kls
             c = basis(kbas)
       !      mijk = mij
       !      KK = natm + 1 
       !      if (c/=II) then 
       !          if (c/=JJ) then 
       !              KK = c
       !              mijk = mijk + 1 
       !          endif
       !      endif
             do l=kbas,nao
                d = basis(l)
       !         mijkl = mijk
       !         LL = natm + 1
                e_coul = 0.5d0*tem1*dm(l,kbas)*aoshell(l,kbas,j,ibas)
                !write(*,*) "a,b,c,d,m",a,b,c,d,e_coul
                if (dabs(e_coul) > 1D-12) then
                  if (kbas/=l) then
                    e_coul = e_coul*2
                  endif 
                  if (ieqj) then 
                    e_coul = e_coul*2 
                  endif

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

             !hf_ene = hf_ene + e_coul
             !  tem = body2(a,b,c,d,natm,m,II,JJ,KK,LL)
                  if (d /= II) then
                      if (d /= JJ) then
                          if (d/= KK) then
                              LL = d
                              m = m + 1
                          endif
                      endif
                  endif
                  e_coul = e_coul/m

                  atom_energy(II) = atom_energy(II) + e_coul
                  atom_energy(JJ) = atom_energy(JJ) + e_coul
                  atom_energy(KK) = atom_energy(KK) + e_coul
                  atom_energy(LL) = atom_energy(LL) + e_coul
                  !write(*,*) "a,b,c,d,m",a,b,c,d,m,II,JJ,KK,LL,e_coul
                  !write(*,*) "atom_energy(1)",atom_energy(1)
                endif
          enddo
        enddo
      enddo
    enddo
    
    call system_clock(time3) 
    write(*,*) "ej time=",time3 - time2

    do ibas = ifi,ils
      a = basis(ibas)
      do l=1,nao
         b = basis(l)
         do kbas = kfi,kls
            c = basis(kbas)
            tem1 = dm(kbas,l)
            do j =1,nao 
               d = basis(j)
               e_coul = 0.5d0*tem1*dm(j,ibas)*((-0.5d0)*aoshell(j,kbas,l,ibas))
               if (dabs(e_coul) > 1D-12) then
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
                           if (d/= KK) then
                               LL = d
                               m = m+1
                           endif
                       endif
                   endif

                   e_coul = e_coul/m 
                   atom_energy(II) = atom_energy(II) + e_coul
                   atom_energy(JJ) = atom_energy(JJ) + e_coul
                   atom_energy(KK) = atom_energy(KK) + e_coul
                   atom_energy(LL) = atom_energy(LL) + e_coul
               endif
          enddo
        enddo
      enddo
    enddo
    
    call system_clock(time4) 
    write(*,*) "ek time=",time4 - time3
  !energy(:) = atom_energy(:)
  !energy = energy + atom_energy
  deallocate(aoshell)
  enddo
enddo
!$omp end parallel do

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

