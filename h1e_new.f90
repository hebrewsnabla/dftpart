subroutine h1e(dm,bas2atm,hcore,natm,nao,energy_h1e)
implicit none

integer i,j,l
integer a,b
integer time1,time2
integer natm,nao
integer bas2atm(nao)
real*8 dm(nao,nao)
real*8 hcore(natm,nao,nao)
real*8 energy_h1e(natm)
real*8 tem

!f2py intent(in) :: dm,hcore,bas2atm
!f2py intent(in) :: natm,nao
!f2py intent(out) :: energy_h1e

! -----------------------------------------------------------------------
! f2py -m h1e_new -c h1e_new.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
! -----------------------------------------------------------------------

energy_h1e = 0.0d0
call system_clock(time1)

do i=1,nao 
  a = bas2atm(i)+1
  do j=1,nao 
    b = bas2atm(j)+1
    do l=1,natm
      tem = dm(i,j) * hcore(l,j,i)
      energy_h1e(l) = energy_h1e(l) + tem*1/2.0d0
      energy_h1e(a) = energy_h1e(a) + tem*1/4.0d0
      energy_h1e(b) = energy_h1e(b) + tem*1/4.0d0
    enddo
  enddo
enddo
call system_clock(time2)
!write(*,*) "Time_h1e=",time2-time1

return
end subroutine h1e

subroutine bgh1e(dm,bas2atm,hcore,natm,nao,energy_h1e)
implicit none

integer i,j,l
integer a,b
integer time1,time2
integer natm,nao
integer bas2atm(nao)
real*8 dm(nao,nao)
real*8 hcore(nao,nao)
real*8 energy_h1e(natm)
real*8 tem

!f2py intent(in) :: dm,hcore,bas2atm
!f2py intent(in) :: natm,nao
!f2py intent(out) :: energy_h1e

! -----------------------------------------------------------------------
! f2py -m h1e_new -c h1e_new.f90 --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
! -----------------------------------------------------------------------

energy_h1e = 0.0d0
call system_clock(time1)

do i=1,nao 
  a = bas2atm(i)+1
  do j=1,nao 
    b = bas2atm(j)+1
      tem = dm(i,j) * hcore(j,i)
      energy_h1e(a) = energy_h1e(a) + tem*1/4.0d0
      energy_h1e(b) = energy_h1e(b) + tem*1/4.0d0
  enddo
enddo
call system_clock(time2)
!write(*,*) "Time_h1e=",time2-time1

return
end subroutine bgh1e

subroutine h1e_inter(dm,bas2atm,bas2frg,atm2frg,hcore,natm,nao,nfrag,energy_h1e, e1_1, e1_2, e1_3)
implicit none

integer i,j,l
integer a,b,aa,bb,ll
integer time1,time2
integer natm,nao,nfrag
integer bas2atm(nao)
integer bas2frg(nao)
integer atm2frg(natm)
real*8 dm(nao,nao)
real*8 hcore(natm,nao,nao)
real*8 energy_h1e(natm)
real*8 e1_1(nfrag)
real*8 e1_2(nfrag,nfrag)
real*8 e1_3(nfrag,nfrag,nfrag)
!real*8 e1_3(natm,natm,natm)
real*8 tem
real*8 H1E_THRESH

call system_clock(time1)
!f2py intent(in) :: dm,hcore,bas2atm,bas2frg,atm2frg
!f2py intent(in) :: natm,nao,nfrag
!f2py intent(out) :: energy_h1e,e1_1, e1_2, e1_3

!f2py depend(natm) :: atm2frg, energy_h1e, hcore
!f2py depend(nao) :: dm, hcore, bas2atm,bas2frg
!f2py depend(nfrag) :: e1_1, e1_2, e1_3

energy_h1e = 0.0d0
e1_1=0.0d0 
e1_2=0.0d0
e1_3=0.0d0
H1E_THRESH = 1.0d-12

write(*,*) hcore(1,1,2), hcore(1,2,3)
  

!$OMP PARALLEL DO schedule(guided) &
!$omp default(private) &
!$omp shared(dm,bas2atm,bas2frg,atm2frg,hcore,nao,natm,H1E_THRESH) &
!$omp reduction(+:energy_h1e,e1_1,e1_2, e1_3)
do i=1,nao 
  a = bas2atm(i)+1
  aa = bas2frg(i)
  do j=1,nao 
    b = bas2atm(j)+1
    bb = bas2frg(j)
    do l=1,natm
      ll = atm2frg(l)
      !write(*,*) aa,bb,ll
      tem = dm(i,j) * hcore(l,j,i)
      !if (dabs(tem) > H1E_THRESH) then 
      if (.true.) then
          energy_h1e(l) = energy_h1e(l) + tem*1/2.0d0
          energy_h1e(a) = energy_h1e(a) + tem*1/4.0d0
          energy_h1e(b) = energy_h1e(b) + tem*1/4.0d0
      
          if (aa==bb) then
            if (aa==ll) then
              e1_1(aa) = e1_1(aa) + tem
            else
              e1_2(aa,ll) = e1_2(aa,ll) + tem
              e1_2(ll,aa) = e1_2(ll,aa) + tem
            endif
          else 
            if ((aa==ll) .or. (bb==ll)) then
              e1_2(aa,bb) = e1_2(aa,bb) + tem
              e1_2(bb,aa) = e1_2(bb,aa) + tem
            else
              !write(*,*) aa,bb,ll,tem
              e1_3(aa,bb,ll) = e1_3(aa,bb,ll) + tem
              e1_3(aa,ll,bb) = e1_3(aa,ll,bb) + tem
              e1_3(bb,aa,ll) = e1_3(bb,aa,ll) + tem
              e1_3(bb,ll,aa) = e1_3(bb,ll,aa) + tem
              e1_3(ll,aa,bb) = e1_3(ll,aa,bb) + tem
              e1_3(ll,bb,aa) = e1_3(ll,bb,aa) + tem
            endif
          endif
      
          !if (a==b) then
          !  if (a==l) then
          !    e1_1(a) = e1_1(a) + tem
          !  else
          !    e1_2(a,l) = e1_2(a,l) + tem
          !    e1_2(l,a) = e1_2(l,a) + tem
          !  endif
          !else 
          !  if ((a==l) .or. (b==l)) then
          !    e1_2(a,b) = e1_2(a,b) + tem
          !    e1_2(b,a) = e1_2(b,a) + tem
          !  else
          !    e1_3 = e1_3 + tem
          !  endif
          !endif
        endif
    enddo
  enddo
enddo
!$omp end parallel do
call system_clock(time2)
!write(*,*) "Time_h1e=",time2-time1

return
end subroutine h1e_inter



