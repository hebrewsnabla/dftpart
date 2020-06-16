subroutine preri(atm,atml,bas,basl,env,envl,cart, nbas,nshls,dm, bas2atm, bas2frg, natm, nfrg, &
atom_energy, & !energyj_one, energyk_one, energyj_two, energyk_two,energy_three,energy_four
e1, e2, e3, e4)
implicit none 
!integer findloc
integer,external :: Ainclude,body3,body4,CINTcgto_spheric, CINTcgto_cart
!real*8,external :: cint2e_sph
integer i,j,k,l,tem,a,b,c,d,flag,di,dj,dk,dl
integer af,bf,cf,df                     ! frags
integer i0,i1,n,p
integer m,II,JJ,KK,LL 
!integer mij,mijk,mijkl
integer natm,nfrg,num
logical cart                            ! cartesian basis or not
!integer nao
integer nshls, nbas
integer bas2atm(nbas)
integer bas2frg(nbas)
integer slice
!real*8 eri(1,nao,nao,nao)
real*8 dm(nbas,nbas)
real*8 e_coul,e_coulA,e_coulB,e_coulC
real*8 e_coul_j,e_coul_k
!real*8 vjk(nao,nao)
!real*8 sumE
real*8 atom_energy(natm+1)
real*8 tem1,tem2
integer time1,time2,time3,time4,t1,t2,t3
integer atml,basl,envl
integer atm(6,atml),bas(8,basl)
!integer atm1(6,atml),bas1(8,basl)
real*8 env(envl)
real*8 buf(:,:,:,:)
real*8 aoshell(:,:,:,:)
!real*8 energy(natm+1)
real*8 schw(:,:)
! --------------------------------------------
!real*8 energyj_one(natm)                                  
!real*8 energyk_one(natm)                                  
!real*8 energyj_two(natm,natm)
!real*8 energyk_two(natm,natm)
!real*8 energy_three(natm,natm,natm)          
!real*8 energy_three          
!real*8 energy_four(natm,natm,natm,natm)
!real*8 energy_four
real*8 e1(nfrg+2)
real*8 e2(nfrg+2,nfrg+2)
real*8 e3(nfrg+2,nfrg+2,nfrg+2)
real*8 e4(nfrg+2,nfrg+2,nfrg+2,nfrg+2)
! --------------------------------------------
integer ncf_sh(:),flcf_sh(:,:),cf2sh(:)
allocatable buf
allocatable ncf_sh
allocatable flcf_sh
allocatable cf2sh
allocatable aoshell
allocatable schw
integer shls(4)
integer ifi,ils,jfi,jls,kfi,kls,lfi,lls,ibas,kbas,jbas,lbas
integer ni,nj
real*8 thresh
real*8 hf_ene
logical ieql,ieqj,keql

! -----------------------------------------------------------------------
!  f2py -m frame_small5 -c frame_small5.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
!  f2py -m frame_small5 -c frame_small5.f90 -L/home/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint
!  f2py -m frame_small6 -c frame_small6.f90 -L/home/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
!  f2py -m new_eda -c new_eda.f90  -L/home/liwei01/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib/ -lcint --fcompiler=intelem --compiler=intelem -liomp5
!  f2py -m new_eda -c new_eda.f90  -L/home/liwei01/zcheng/pyscf/pyscf/lib/ -lcint --fcompiler=intelem --compiler=intelem -liomp5
! -----------------------------------------------------------------------

!f2py intent(in) :: atm,atml,bas,basl,env,envl,cart, nshls
!f2py intent(in) :: nbas,dm,bas2atm,bas2frg,natm,nfrg
!f2py intent(out) :: atom_energy, e1,e2,e3,e4
!energyj_one,energyk_one,energyj_two,energyk_two,energy_three,energy_four

!f2py depend(nbas) :: dm, bas2atm, bas2frg
!f2py depend(basl) :: bas
!f2py depend(envl) :: env
!f2py depend(atml) :: atm
!f2py depend(natm) :: atom_energy
!f2py depend(nfrg) :: e1,e2,e3,e4
!,energyj_one,energyk_one,energyj_two,energyk_two,energy_three,energy_four

!vjk = 0.0d0
!listA = 0.0d0
call system_clock(t1)
e_coul = 0.0d0
atom_energy = 0.0d0
! --------------------------------------
!energyj_one = 0.0d0
!energyk_one = 0.0d0
!energyj_two = 0.0d0
!energyk_two = 0.0d0
!energy_three = 0.0d0
!energy_four = 0.0d0
e1 = 0.0d0
e2 = 0.0d0
e3 = 0.0d0
e4 = 0.0d0
! --------------------------------------

!energy = 0.0d0

!sumE =0.0d0
!hf_ene = 0.0d0

!nbas = nao
!natm = natm
!write(*,*) 'nbas',nbas
!do b=1,nbas
!    write(*,*) b,bas2frg(b)
!enddo


!------------------------------------------------------------ 
! obtain the map between contracted functions and contracted shells
!------------------------------------------------------------

allocate(ncf_sh(nshls))
allocate(flcf_sh(2,nshls))
allocate(schw(nshls,nshls))

!call schwarz(schw,nshl)
!---------------------------------------------------

do i=1,nshls
    if (cart) then
        ncf_sh(i) = cintcgto_cart(i-1,bas)
    else
        ncf_sh(i) = cintcgto_spheric(i-1,bas)
    endif
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
        if (cart) then
            call cint2e_cart(buf, shls, atm, 0, bas, 0, env, 0_8)
        else
            call cint2e_sph(buf, shls, atm, 0, bas, 0, env, 0_8)
        endif
        schw(i,j)=maxval(buf)
        if (i/=j) schw(j,i)=schw(i,j)
        deallocate(buf)
    enddo
enddo
thresh = 1.0D-10

call system_clock(t1)

!------------------------------------------------------------------
!$OMP PARALLEL DO schedule(guided) &
!$omp default(private) &
!$omp shared(nshls,atm,bas,env,flcf_sh, &
!$omp dm,bas2atm,bas2frg,nbas,natm,schw,thresh,ncf_sh) &
!$omp reduction(+:atom_energy, e1,e2,e3,e4)

! reduction(+:atom_energy,energyj_one, energyk_one, energyj_two, energyk_two, energy_three, energy_four)
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
    do l=i,nshls
      shls(3) = l - 1 
      lfi = flcf_sh(1,l)
      lls = flcf_sh(2,l)
      dl = ncf_sh(l)
      do j=k,nshls
        shls(1) = j - 1 
        jfi = flcf_sh(1,j)
        jls = flcf_sh(2,j)
        dj = ncf_sh(j)
      !write(*,*) di,dj,dk,dl
        if (dsqrt(schw(j,k))*dsqrt(schw(i,l))>thresh) then
            allocate(buf(dj,dk,dl,di))
            if (cart) then
                call cint2e_cart(buf, shls, atm, 0, bas, 0, env, 0_8)
            else
                call cint2e_sph(buf, shls, atm, 0, bas, 0, env, 0_8)
            endif
            aoshell(jfi:jls,:,lfi:lls,:) = buf(:,:,:,:)
            !write(*,*) "buf\n",buf
            deallocate(buf)
        endif
      enddo
    enddo
    !call system_clock(time2) 
    !write(*,*) "Calc eri time=",time2-time1
    do ibas = ifi,ils
       a = bas2atm(ibas)
       af = bas2frg(ibas)
       do j =ibas,nbas
          b = bas2atm(j)
          bf = bas2frg(j)
          !write(*,*) j, bf
          ieqj = (ibas/=j)
          tem1 = dm(j,ibas)
          do kbas = kfi,kls
             c = bas2atm(kbas)
             cf = bas2frg(kbas)
             do l=kbas,nbas
                d = bas2atm(l)
                df = bas2frg(l)
                keql = (kbas/=l)
                e_coul_j = 0.5d0*tem1*dm(l,kbas)*aoshell(l,kbas,j,ibas)
                e_coulA =  0.5d0*dm(kbas,j)*dm(l,ibas)*((-0.5d0)*aoshell(l,kbas,j,ibas))
                e_coulB =  0.5d0*dm(kbas,ibas)*dm(l,j)*((-0.5d0)*aoshell(l,kbas,j,ibas))
                e_coulC =  e_coulA + e_coulB 
                if (ibas/=j) then
                    if (kbas/=l) then 
                        e_coul_k = e_coulC*2 
                        e_coul_j = e_coul_j*4
                    else 
                        e_coul_k = e_coulC 
                        e_coul_j = e_coul_j*2
                    endif
                else
                    if (kbas/=l) then 
                        e_coul_k = e_coulC 
                        e_coul_j = e_coul_j*2
                    else 
                        e_coul_k = e_coulA
                        e_coul_j = e_coul_j
                    endif
                endif

                e_coul = e_coul_j + e_coul_k
                !if (bf==4) then
                !    write(*,*) af,bf,cf,df, e_coul
                !endif
                
                !if (dabs(e_coul) > 1D-12) then
                if (.true.) then
                  
                  !tem = body3(a,b,c,d,natm,m,II,JJ,KK,LL)

                  atom_energy(a) = atom_energy(a) + e_coul/4.0
                  atom_energy(b) = atom_energy(b) + e_coul/4.0
                  atom_energy(c) = atom_energy(c) + e_coul/4.0
                  atom_energy(d) = atom_energy(d) + e_coul/4.0

                  !write(*,*) "1111"  
                  !if ( m ==1 ) then
                  !    energyj_one(II) = energyj_one(II) + e_coul_j 
                  !    energyk_one(II) = energyk_one(II) + e_coul_k 
                  !else if ( m ==2 ) then 
                  !    energyj_two(II,JJ) = energyj_two(II,JJ) + e_coul_j 
                  !    energyk_two(II,JJ) = energyk_two(II,JJ) + e_coul_k 
                  !    energyj_two(JJ,II) = energyj_two(JJ,II) + e_coul_j 
                  !    energyk_two(JJ,II) = energyk_two(JJ,II) + e_coul_k 
                  !else if ( m ==3 ) then 
                  !    !energy_three(II,JJ,KK) = energy_three(II,JJ,KK) + e_coul 
                  !    energy_three = energy_three + e_coul 
                  !else 
                  !    !energy_four(II,JJ,KK,LL) = energy_four(II,JJ,KK,LL) + e_coul 
                  !    energy_four = energy_four + e_coul 
                  !endif 
                tem = body4(af,bf,cf,df,m,II,JJ,KK,LL)
                !if (bas2frg(1)==2) then
                !    stop
                !endif
                               
                if (m==1) then
                    e1(II) = e1(II) + e_coul
                else if (m==2) then
                    if ((II==3) .and. (JJ==4)) then
                        !if (dabs(aoshell(l,kbas,j,ibas)) > 1d-15) then
                        !write(*,*) af,bf,cf,df, e_coul, aoshell(l,kbas,j,ibas)
                        !endif
                    endif
                    e2(II,JJ) = e2(II,JJ) + e_coul
                else if (m==3) then
                    e3(II,JJ,KK) = e3(II,JJ,KK) + e_coul
                else
                    e4(II,JJ,KK,LL) = e4(II,JJ,KK,LL) + e_coul
                endif
                  

                endif
          enddo
        enddo
      enddo
    enddo
  
    !call system_clock(time3) 
    !write(*,*) "ej-ek time=",time3 - time2

  deallocate(aoshell)
  enddo
enddo
!$omp end parallel do
    
call system_clock(t2)
write(*,*) "Run-time",t2-t1
!write(*,*) "enrgy_four",energy_four
!write(*,*) "size",size(energy_four)

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

function body3(a,b,c,d,natm,m,II,JJ,KK,LL)
   implicit none
   integer :: a,b,c,d,natm,m,II,JJ,KK,LL
   integer :: body3
   II = natm+1
   JJ = natm+1
   KK = natm+1
   LL = natm+1
   m = 0
   if (a/=b) then
       II = a
       JJ = b
       m = m+2
       if ((c /= II) .and. (c /= JJ)) then
           KK = c
           m = m+1
           if ((d /= II) .and. (d /= JJ) .and. (d/=KK)) then
               LL = d
               m = m+1
           endif
       else    
           if ((d /= II) .and. (d /= JJ) ) then
               KK = d
               m = m+1
           endif
       endif
   else
       II = a
       m = m+1
       if (c /= II) then
           JJ = c
           m = m+1
           if ((d /= II) .and. (d /= JJ) ) then
               KK = d
               m = m+1
           endif
       else
           if (d /= II)  then
               JJ = d
               m = m+1
           endif
       endif
   endif


    body3 = 0
    !PRINT *,a,b,c,d,'->', II, JJ ,KK ,LL,m,natm

return
end function body3


function body4(a,b,c,d,m,II,JJ,KK,LL)
   implicit none
   integer :: a,b,c,d,m,II,JJ,KK,LL
   integer :: tmp
   integer :: body4
   II = 0
   JJ = 0
   KK = 0
   LL = 0
   m = 0
   if (a/=b) then
       II = a
       JJ = b
       m = m+2
       if ((c /= II) .and. (c /= JJ)) then
           KK = c
           m = m+1
           if ((d /= II) .and. (d /= JJ) .and. (d/=KK)) then
               LL = d
               m = m+1
           endif
       else    
           if ((d /= II) .and. (d /= JJ) ) then
               KK = d
               m = m+1
           endif
       endif
   else
       II = a
       m = m+1
       if (c /= II) then
           JJ = c
           m = m+1
           if ((d /= II) .and. (d /= JJ) ) then
               KK = d
               m = m+1
           endif
       else
           if (d /= II)  then
               JJ = d
               m = m+1
           endif
       endif
   endif
    ! sorting to make II < JJ < KK < LL
    if (m>1) then
        if (JJ < II) then
            tmp = JJ
            JJ = II
            II = tmp
        endif
    endif
    if (m>2) then
        if (KK < JJ) then
            if (KK < II) then
                tmp = KK
                KK = JJ
                JJ = II
                II = tmp
            else
                tmp = KK
                KK = JJ
                JJ = tmp
            endif
        endif
    endif
    if (m==4) then
        if (LL < KK) then
            if (LL < JJ) then
                if (LL < II) then
                    tmp = LL
                    LL = KK
                    KK = JJ
                    JJ = II
                    II = tmp
                else
                    tmp = LL
                    LL = KK
                    KK = JJ
                    JJ = tmp
                endif
            else
                tmp = LL
                LL = KK
                KK = tmp
            endif
        endif
    endif               

    body4 = 0
    !PRINT *,a,b,c,d,'->', II, JJ ,KK ,LL,m,natm

return
end function body4

