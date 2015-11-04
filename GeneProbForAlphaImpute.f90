!GeneProbForAlphaPhase
!Description:	GeneProb code that is run for multiple SNPs, by SNP
!Date:			30 Mar 2011
!Author:		Brian Kinghorn (this code put together by MA Cleveland)
!Usage::
!Notes:			Input genotype file: IDENT, SIRE, DAM, MID1,...,MID(n)
!				Assumes Geno format is 0,1,2
!				Unknown genotypes should be "9", but can be "3"
!
!				No header
!
!				Need paramter file=GeneProbSpec.txt
!					nAnis
!					nSnp
!					inputFile
!					startSnp
!					endSnp
!
!References:	This is same as GeneProb_sub.f90, but requires a starting and stopping SNP
!Updates:

module Global_GP
	implicit none

	character(len=1000) :: inputFile,outputFile
	integer :: nSnp,nAnis,startSnp,endSnp

end module

module common_GP
  !implicit double precision (a-h,o-z)
  implicit none
  INTEGER mm,nn, mxeq,mxaneq
  INTEGER, ALLOCATABLE :: prog(:),MATE(:),NEXT(:),IFIRST(:)
  REAL(KIND=8), ALLOCATABLE :: POST(:,:)
end module

module length_of_ID
      INTEGER, PARAMETER :: lengan=16  ! needs to be big enough to handle the internally generated dummy IDs, eg DUM00001 => 8 digits.  Line~814
end module

module commonbits
 use length_of_ID
 CHARACTER*(lengan), allocatable:: id(:),sire(:),dam(:)
 INTEGER, allocatable:: seqid(:),seqsire(:),seqdam(:),passedorder(:),phenhold(:)
 INTEGER :: nobs,Imprinting,PauseAtEnd,nfreq_max,phenotypes,maxfs,maxmates,nfamilies
 REAL (KIND=8)         :: pprior, qprior, g(25,0:2),StopCrit
 integer :: phenotype(25)

   REAL (KIND=8)         :: pprior_hold 
   REAL (KIND=8)         :: qprior_hold 
   REAL (KIND=8)         :: StopCrit_hold 

   INTEGER         :: phenotypes_hold
   INTEGER :: Imprinting_hold 
   INTEGER :: PauseAtEnd_hold 
   INTEGER :: nobs_hold 
   INTEGER :: nfreq_max_hold 
end module

module GPinput
	implicit none
	integer (KIND=1),allocatable,dimension(:,:) :: InputGenos,tmpInputGenos
	real(kind=4),allocatable,dimension(:,:) :: Probs00,Probs01,Probs10,Probs11,GPI
	real(kind=4),allocatable,dimension(:) :: OutputMaf
end module

program GeneProb_sub
	use Global_GP
	use GPinput

	integer :: i
	character(len=7) :: cm !use for formatting output - allows for up to 1 million SNPs

	call ReadParms

	open(UNIT=19,FILE=trim(inputFile),STATUS="old") !INPUT FILE
	open (UNIT=111,FILE=trim(outputFile),STATUS="unknown") !OUTPUT FILE
  open (UNIT=222,FILE='GPI.txt',STATUS="unknown")

	allocate(Probs00(nAnis,endSnp-startSnp+1))
 	allocate(Probs01(nAnis,endSnp-startSnp+1))
 	allocate(Probs10(nAnis,endSnp-startSnp+1))
	allocate(Probs11(nAnis,endSnp-startSnp+1))
 	allocate(GPI(nAnis,endSnp-startSnp+1))
 	allocate(OutputMaf(endSnp-startSnp+1))
 	
	call preprocessGeneprob
	do i=1,endSnp-startSnp+1 !calculate geneprobs for all individuals
		call geneprob(i)
 	end do

 	write(cm,'(I7)') nSnp !for formatting
 	cm = adjustl(cm)
 	do i=1,nAnis
!		write(111,'(I7,'//cm//'(1x,F8.4))') i, (Probs00(i,:))
!		write(111,'(I7,'//cm//'(1x,F8.4))') i, (Probs01(i,:))
!		write(111,'(I7,'//cm//'(1x,F8.4))') i, (Probs10(i,:))
!		write(111,'(I7,'//cm//'(1x,F8.4))') i, (Probs11(i,:))
		write(111,'(i16,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4)') i, (Probs00(i,:))
		write(111,'(i16,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4)') i, (Probs01(i,:))
		write(111,'(i16,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4)') i, (Probs10(i,:))
		write(111,'(i16,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4,20000f8.4)') i, (Probs11(i,:))

    ! write(111,'(I7,'//cm//'(1x,F8.4))') i, (GPI(i,:))
		write(222,'(i16,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4,20000f9.4)') i, (GPI(i,:))
 	end do

 	deallocate(Probs00)
 	deallocate(Probs01)
 	deallocate(Probs10)
	deallocate(Probs11)
 	deallocate(GPI)
 	

 	open(unit=1002,file="GpDone.txt",status="unknown")
 	write(1002,*) "Geneprobs done"

 	open(unit=2002,file="MinorAlleleFrequency.txt",status="unknown")
	do i=1,((endSnp-startSnp)+1)
 		write(2002,*) OutputMaf(i)
 	enddo	
	deallocate(OutputMaf)

end program

subroutine ReadParms
	use Global_GP

	character :: tmp

	open(unit=1001,file="GeneProbSpec.txt",status="old")
	read(1001,*) tmp,nAnis
	read(1001,*) tmp,nSnp
	read(1001,*) tmp,inputFile
	read(1001,*) tmp,outputFile
	read(1001,*) tmp,startSnp
	read(1001,*) tmp,endSnp
	close(1001)

end subroutine ReadParms

subroutine preprocessGeneprob
	use Global_GP
	use common_GP
	use commonbits
	use GPinput
	implicit none

	integer :: i,j,k,dum

	pprior_hold = 0.5
	qprior_hold = 1 - pprior
	phenotypes_hold = 3
	Imprinting_hold = 1
	PauseAtEnd_hold = 0
	nfreq_max_hold = 50
	StopCrit_hold = 0.0001

	nobs_hold = nAnis
	nobs=nobs_hold


	ALLOCATE(id(0:nobs),sire(nobs),dam(nobs),seqid(nobs),&
         seqsire(nobs),seqdam(nobs),passedorder(nobs))

	phenotype(1) = 0
	phenotype(2) = 1
	phenotype(3) = 2
	!g(1,0) = 1.0
	!g(1,1) = 0.0
	!g(1,2) = 0.0
	!g(2,0) = 0.0
	!g(2,1) = 1.0
	!g(2,2) = 0.0
	!g(3,0) = 0.0
	!g(3,1) = 0.0
	!g(3,2) = 1.0

	g(1,0) = 0.99
	g(1,1) = 0.005
	g(1,2) = 0.005
	g(2,0) = 0.005
	g(2,1) = 0.99
	g(2,2) = 0.005
	g(3,0) = 0.005
	g(3,1) = 0.005
	g(3,2) = 0.99

	allocate(tmpInputGenos(nAnis,nSnp))
	allocate(InputGenos(nAnis,endSnp-StartSnp+1))
	!read(19,*)
	do i=1,nobs
		read (19,*) seqid(i), seqsire(i), seqdam(i), tmpInputGenos(i,:)
		passedorder(i) = i
	end do

	do i=1,nobs
		do j=1,nSnp
			if(tmpInputGenos(i,j)==3) tmpInputGenos(i,j)=9
		end do
	end do

	do i=1,nobs
		k=1
		do j=startSnp,endSnp
			InputGenos(i,k)=tmpInputGenos(i,j)
			k=k+1
		end do
	end do

	call GetMaxFS!!(maxfs, maxmates, nfamilies)

end subroutine preprocessGeneprob

subroutine geneprob(currentSnp)
	  use Global_GP
	  use common_GP
      USE commonbits
      use GPinput
      implicit none

      integer, intent(in)	:: currentSnp
      integer               :: i, j, k, l, i2, i3, iticks1, iticks2
      integer               :: ia, is, idd, ifreq_iterate, maxint, maxiter, itersused, kl, kc, kd, kj, nfams, last
      integer               :: nf, im, ns, mf, iaa, ii, ms, m, n, maxvalspost, ierrors, iflag, nwritten
      integer               :: f,ff,nonzed(3,3),ntype(3,3,3)
      integer               :: maxRegpoints,LeastPositive, LeastNegative, HoldInt, LimitAnimals, LimitNumber

      real (kind=8)         :: tsum, prod, ProbFit, SumFreq, IMPratio, p12, p21, LeastPositiveValue, LeastNegativeValue
      REAL (KIND=8)         :: spost(3),dpost(3),fpost(3),tpost(3),temp(3),sum1(3),sum2(3),sum3(3),pt(3,3,3)
      REAL (KIND=8)         :: phethw,phomhw,probindex  ! this is needed for info - or compile with dble.  Don't know why!
      REAL (KIND=8)         :: s0,s1,s2, areg, breg, meanX, meanY, sumY, sumXY, sumX, sumX2

      INTEGER, ALLOCATABLE  :: phen(:)
      INTEGER, ALLOCATABLE  :: nmem(:),isib(:,:),damGP(:),p1(:),p2(:)  ! note this is a different p1,p2 to sequence's

      REAL (KIND=8), allocatable :: ant(:,:),term(:,:),phom(:),phet(:),pnor(:),work(:,:,:),freq(:,:)
      REAL (KIND=8), allocatable :: pHold(:), pResult(:), pDev(:)


      CHARACTER(256)        :: cl_arg, infile, outfile, HoldStr

phenotypes=phenotypes_hold
pprior=pprior_hold
qprior=qprior_hold
nfreq_max=nfreq_max_hold
Imprinting=Imprinting_hold
PauseAtEnd=PauseAtEnd_hold
StopCrit=StopCrit_hold


! ----------------------------------------------------------------
!  P-MATRIX: PROB. OF OFFSPRING GENOTYPE GIVEN GENOTYPE OF PARENTS
! ----------------------------------------------------------------

LimitAnimals = 0 ! 1 to invoke Limit
LimitNumber = 250

      pt=0.  !  =log(1)!  the log(zero) elements should not be required ar nonzed and ntype below control addressing
      pt(1,2,1)=log(.5)
      pt(2,1,1)=log(.5)
      pt(2,2,1)=log(.25)
      pt(1,2,2)=log(.5)
      pt(2,1,2)=log(.5)
      pt(2,2,2)=log(.5)
      pt(2,3,2)=log(.5)
      pt(3,2,2)=log(.5)
      pt(2,2,3)=log(.25)
      pt(2,3,3)=log(.5)
      pt(3,2,3)=log(.5)
      nonzed(1,1)=1
      ntype(1,1,1)=1
      nonzed(1,2)=2
      ntype(1,2,1)=1
      ntype(1,2,2)=2
      nonzed(1,3)=1
      ntype(1,3,1)=2
      nonzed(2,1)=2
      ntype(2,1,1)=1
      ntype(2,1,2)=2
      nonzed(2,2)=3
      ntype(2,2,1)=1
      ntype(2,2,2)=2
      ntype(2,2,3)=3
      nonzed(2,3)=2
      ntype(2,3,1)=2
      ntype(2,3,2)=3
      nonzed(3,1)=1
      ntype(3,1,1)=2
      nonzed(3,2)=2
      ntype(3,2,1)=2
      ntype(3,2,2)=3
      nonzed(3,3)=1
      ntype(3,3,1)=3

!!PRINT*, ''
!!PRINT*, 'Program GENEPROB Version 2.7'
!!PRINT*, 'Written by Brian Kinghorn and Richard Kerr'
!!PRINT*, 'Copyright Authors/University of New England'
!!PRINT*, 'Please do not provide a copy of this program to other parties.'
!!
!!
!!if (LimitAnimals==1) then
!!  PRINT*, ''
!!  PRINT*, '******************************************************************************'
!!  PRINT*, '     This is a demonstration version. Commercial use is prohibited.'
!!  PRINT*, '  It is limited to results on ', LimitNumber, ' animals, including unlisted parents.'
!!  PRINT*, '******************************************************************************'
!!  PRINT*, ''
!!endif
!!
!!PRINT*, 'This version expects ID fields of maximum length ',lengan
!!
!!***** call getcl(cl_arg) !only supported by Lahey
!!	 call getarg(2,cl_arg) !PGI command - actually shoudn't need this, there won't be any commands
!!    cl_arg = TRIM(cl_arg)
!!
!!     infile=''
!!     outfile=''
!!     HoldStr=''
!!     j=1
!!     k=1
!!     DO i=1, LEN_TRIM(cl_arg)
!!        if (cl_arg(i:i)==',' .or. cl_arg(i:i)==';') then
!!          infile=TRIM(HoldStr)
!!          j=j+1
!!          IF(j>2) STOP ' Too many command line arguments - Stopping now.'
!!          HoldStr=''
!!          k=1
!!         else
!!          if ( .NOT. (k==1 .and. cl_arg(i:i)==' ') ) then    !ignore leading spaces
!!            HoldStr(k:k)=cl_arg(i:i)
!!            k=k+1
!!          end if
!!         end if
!!     ENDDO
!!     outfile=TRIM(HoldStr)
!!
!!     PRINT'(1x,a10,a60)', ' Infile: ', infile
!!     PRINT'(1x,a10,a60)', 'Outfile: ', outfile


!!    open(3,file=infile,status='unknown')

!!call system_clock(iticks1,i2,i3)
!!   read(3,*)
!!    read(3,*)
!!    read(3,*) pprior,phenotypes, Imprinting, PauseAtEnd

!!    nobs=0
!!101 READ(3,'(a13)',END=102) HoldStr
!!    nobs=nobs+1
!!    GOTO 101
!!
!!102  nobs=nobs-(4+phenotypes)  ! to account for header lines

!!write(*,*) ' Number of records: ',nobs

!!ALLOCATE(id(0:nobs),sire(nobs),dam(nobs),seqid(nobs),&
!!         seqsire(nobs),seqdam(nobs),phenhold(0:nobs)) !,passedorder(nobs)
!!
!!REWIND(3)
!!
!!       read(3,*)
!!       read(3,*)
!!       read(3,*) pprior,phenotypes
!!       read(3,*)
!!
!!        qprior=1-pprior
!!
!!       do i = 1, phenotypes
!!        read(3,*) phenotype(i), g(i, 0), g(i, 1), g(i, 2)
!!!       write(*,*) phenotype(i), g(i, 0), g(i, 1), g(i, 2)
!!       enddo
!!
!!       read(3,*) nfreq_max
!!       read(3,*) StopCrit
!!
!!       if (StopCrit < 0.000001) StopCrit = 0.000001  ! The regression get wobbly below this because of low X variance
!!
!!       read(3,*)
!!
!!
!!! 1    read(3,*,end=9) ia,is,idd,tag,phen(ia)
!!do i=1,nobs
!!  read(3,*) id(i),sire(i),dam(i),phenhold(i)
!!enddo
!!CLOSE(3)

!!****************************************************************
!!MC: set all parameters that are normally in the input file
!!MC: need to put this somewhere else
!!MC: need to populate these arrays from the data already stored
allocate(phenhold(0:nobs))
do i=1,nobs
		phenhold(i) = InputGenos(i,currentSnp)
end do


!!************************************************************************
!!call pvseq(1)
!!call GetMaxFS(maxfs, maxmates, nfamilies)
!!Print*, 'Maxfs:     ', maxfs
!!Print*, 'Maxmates:  ', maxmates
!!Print*, 'nFamilies: ', nFamilies


!maxmates = maxmates!*20     ! why ??  Fixed in 2.7
!maxfs = maxfs

mxaneq= nobs   ! was multiplied up.  For old sequencer I think.

mm=0
nn=mxaneq

ALLOCATE(MATE(0:2*MXANEQ),NEXT(2*MXANEQ),IFIRST(0:MXANEQ),POST(3,0:2*MXANEQ),phen(mxaneq), &
         prog(0:2*MXANEQ),&
         p1(mxaneq),p2(mxaneq),ant(3,0:mxaneq),phom(0:mxaneq),phet(0:mxaneq), &
         freq(3,0:mxaneq),pnor(0:mxaneq))

phen=9 ! covers unlisted parents
mate=0
next=0
ifirst=0
post=0.
prog=0.
phom=0.
phet=0.
pnor=0.

! nmem needs 2*maxfs.  If future problem look for additional need related to maxmates etc.
HoldInt = MAX(2*maxfs,maxmates)
! work dimension is the max value of nmem itself


ALLOCATE(nmem(HoldInt),isib(maxmates,2*maxfs),damGP(maxmates),work(3,2*maxfs,2*maxfs),term(3,HoldInt))
HoldInt=0



!ALLOCATE(nmem(maxmates),isib(maxmates,maxfs),damGP(maxmates),work(3,maxfs,maxfs),term(3,maxmates))
!maxmates used, as nFamilies is done within sire.
!ALLOCATE(nmem(nfamilies),isib(nfamilies,maxfs),damGP(nfamilies),work(3,maxfs,maxfs),term(3,nfamilies))

isib=0

phenhold(0)=9  ! unknown

do ia=1,nobs
  ! if (seqid(ia).ne.ia) then
  !  STOP 'blowout'
  ! end if
   is=seqsire(ia)
   idd=seqdam(ia)
   p1(ia)=is
   p2(ia)=idd
   call LNKLST(is,idd,ia,1)
   call LNKLST(idd,is,ia,0)
end do

do ia=1,nobs
   phen(ia)=phenhold(passedorder(ia))
   iflag=0
   if (phen(ia).eq.9) THEN
    iflag=1
    freq(1,ia) =log(1.)
    freq(2,ia) =log(1.)
    freq(3,ia) =log(1.)
   endif
    do i = 1, phenotypes
     IF (phen(ia).eq.phenotype(i)) THEN
      iflag=1
      IF(g(i,0).lt..000000001)then
        freq(1,ia) =-9999
        else
        freq(1,ia) =log(g(i, 0))
      endif
      IF(g(i,1).lt..000000001)then
        freq(2,ia) =-9999
        else
        freq(2,ia) =log(g(i, 1))
      endif
      IF(g(i,2).lt..000000001)then
        freq(3,ia) =-9999
        else
        freq(3,ia) =log(g(i, 2))
      endif
     endif
    enddo
    if(iflag==0) then
       print*, "Unregistered phenotype ", phen(ia), " for individual ", id(ia)
       stop "Aborting run"
    endif
end do

deallocate (phenhold)

j=0
Sumfreq=0.
do ia=1,nobs
    do i = 1, phenotypes
     IF (phen(ia).eq.phenotype(i)) THEN
        j=j+1
        SumFreq = SumFreq + g(i, 1) + 2*g(i, 2)
     endif
    enddo
end do

!!PRINT'(a32,i6)', ' Number of phenotyped animals: ', j

if (phenotypes==3) then
  if(g(1,0)>.99999 .and. g(2,1)>.99999 .and.g(3,2)>.99999) then   ! Identity - eg DNA test
    if(j>0)then
!!      PRINT'(a32,f12.8)', ' Raw observed frequency: ', SumFreq/(2.*float(j))
    else
!!      PRINT'(a32,a36)', ' Raw observed frequency: ', ' No genotypes/incidence observed'
    endif
  endif
end if

maxRegpoints=5
!!PRINT'(a32,i6)', ' Max Regression points: ', maxRegpoints
!!PRINT*, ' Segregation analysis ...'
!!IF(nfreq_max>0) PRINT'(2a8,5a12)', 'iter','Shells','Prior used', 'Result', 'Difference', 'Next prior'


ifreq_iterate = -1
if (nfreq_max==1) nfreq_max=2
ALLOCATE (pHold(0:nfreq_max), pResult(0:nfreq_max), pDev(0:nfreq_max))

if (nfreq_max==0)then
  pHold(0)=pprior  ! one hit only
else
  pHold(0)=0.001
  pHold(1)=0.999
endif

pDev(0)=0.
LeastPositiveValue =  999.
LeastNegativeValue = -999.
LeastPositive = 0
LeastNegative = 0


do WHILE (ifreq_iterate < nfreq_max)

      ifreq_iterate = ifreq_iterate + 1

      pprior = pHold(ifreq_iterate)

      qprior = 1-pprior

!   initialise
      post=0.
      phet=0.
      do i=1,nobs
         ant(1,i)=log(qprior*qprior)
         ant(2,i)=log(2.0*pprior*qprior)
         ant(3,i)=log(pprior*pprior)
      enddo


! ----------------------------------------------
! Prob(Gi) = ( Ai f(Yi|Gi) PROD - mates Pi ) / L
!     where L = SUM-Gi Ai f(Yi|Gi) PROD - mates Pi
!    Ai is the joint probability of phenotypes of members anterior to i
!       genotype Gi for i
!    PROD over mates Pi is the conditional probability of phenotypes of
!     posterior to i, given i has genotype Gi
! ----------------------------------------------
! BUILD A LINKLIST. EACH PARENT ANIMAL HAS A ROW, ALONG THE COLUMNS ARE
! OF THEIR MATES AND PROGENY. MATES CAN BE REPEATED AT SUCCESSIVE NODES
! REFLECT FULL SIB FAMILIES. NOTE, THE FIRST NODE FOR A PARTICULAR MATE
! CONTAIN THE POST(i,j) term, (the jth mate of the ith animal)
! To initialise the iterative peeling up and peeling down cycles, the an
! term for founder animals is set equal to the HW probs. Post. terms for
! animals and ant. terms for non founders are set to 1 (reflecting no
! information)
      maxint=9999999
       maxiter=7
       itersUsed=maxiter
 888  format(3(f9.2,1x,f9.2,3x))
      do kl=1,maxiter
!!      IF(nfreq_max==0) print*,'  Iteration number ',kl
      !print*,'  Iteration number ',kl
! ----------------------------------------------------------------
! PEEL DOWN, IE CONDENSE INFO ON MUM AND DAD ONTO PROGENY
! START WITH OLDEST ANIMAL IN THE LIST. WHEN DESCENDING WE CALCULATE
! ANT() TERMS FOR ALL PROGENY. PICK OUT PROGENY FROM THE ROWS OF
! THE PARENT WITH ON AVERAGE THE MOST MATES OR PROGENY, USUALLY THE SIRE
! IGNORE THE ROWS OF THE OTHER PARENT.
! ----------------------------------------------------------------
         do is=1,nobs
            kc=ifirst(is)
            if(kc.eq.0.or.kc.gt.nobs) goto 50
!     load up the mates of this sire and collect posterior terms
            idd=mate(kc)
            do i=1,3
               spost(i)=0.0
            enddo
            nfams=0
            last=maxint
            do while (kc.ne.0)
               ia=prog(kc)
               if(idd.ne.last) then
                  nfams=nfams+1
! if(nfams>HoldInt) then  ! nfams goes to maxmates here
!  print*, nfams
!  HoldInt=nfams
! endif
 				  nmem(nfams)=1
                  damGP(nfams)=idd
                  isib(nfams,1)=ia
                  term(1,nfams)=post(1,kc)
                  term(2,nfams)=post(2,kc)
                  term(3,nfams)=post(3,kc)
                  spost(1)=spost(1)+term(1,nfams)
                  spost(2)=spost(2)+term(2,nfams)
                  spost(3)=spost(3)+term(3,nfams)
               else
                  nmem(nfams)=nmem(nfams)+1

!if(nfams>maxmates) then
! print*, 'nfams',nfams, maxmates
!endif

!if(nmem(nfams)>2*maxfs-5) then  ! nmem(nfams) goes to 2*maxfs
! print*, 'nmem', nfams, nmem(nfams),maxfs
!endif

                  isib(nfams,nmem(nfams))=ia
               endif
               last=idd
               kc=next(kc)
               idd=mate(kc)
            enddo

            do nf=1,nfams     ! GO THROUGH THROUGH THE MATES "idd" OF "is
               idd=damGP(nf)
!     collect posterior term for damGP "idd" through all its mates
               do i=1,3
                  dpost(i)=0.
               enddo
               kc=ifirst(idd)
               im=mate(kc)
               last=maxint
               do while (kc.ne.0)
                  if(im.ne.is.and.im.ne.last) then
                     dpost(1)=dpost(1)+post(1,kc)
                     dpost(2)=dpost(2)+post(2,kc)
                     dpost(3)=dpost(3)+post(3,kc)
                  endif
                  last=im
                  kc=next(kc)
                  im=mate(kc)
               enddo
!     correct posterior prob of "is" for "idd"
               tpost(1)=spost(1)-term(1,nf)
               tpost(2)=spost(2)-term(2,nf)
               tpost(3)=spost(3)-term(3,nf)
!    for this damGP "idd", and for each of her progeny "ia" to "is"
!    mark out the full sibs "iaa" to "ia" all the k mates of "iaa"
!    and store the term prod-k post(ia,iaa) in a work vector
               do ns=1,nmem(nf)
                  ia=isib(nf,ns)
                  do mf=1,nmem(nf)
                     iaa=isib(nf,mf)
                     if(iaa.ne.ia) then
                        do ii=1,3
                           fpost(ii)=0.
                        enddo
                        kc=ifirst(iaa) ! collect mates of "iaa"
                        ms=mate(kc)
                        last=maxint
                        do while (kc.ne.0)
                           if(ms.ne.last) then
                              fpost(1)=fpost(1)+post(1,kc)
                              fpost(2)=fpost(2)+post(2,kc)
                              fpost(3)=fpost(3)+post(3,kc)
                           endif
                           last=ms
                           kc=next(kc)
                           ms=mate(kc)
                        enddo
                        do i=1,3
                           work(i,ns,mf)=fpost(i)
                        enddo
                     endif
                  enddo
               enddo
!  NOW WE ARE READY To CALCULATE THE BLOODY THING
               do ns=1,nmem(nf)
                  ia=isib(nf,ns)
!	Here we will get the anterior probability for the progeny in question.  For appendix equations' m, f, s and i:
!	m is here is   - the sire of i as in the outer loop
!	f is here imum - the dam of i
!	s is here iaa  - the sibs of i
!	i is here ia   - the progeny animal
                  do i=1,3
                     mm=0
                     do m=1,3
                        if(i.eq.1.and.m.eq.3)goto 10
                        if(i.eq.3.and.m.eq.1)goto 10
                        mm=mm+1
                        ff=0
                        do f=1,3
                           if(i.eq.1.and.f.eq.3)goto 20
                           if(i.eq.2.and.(m+f.lt.3.or.m+f.gt.5))goto 20
                           if(i.eq.3.and.f.eq.1)goto 20
                           ff=ff+1
                           prod=0.
                           do mf=1,nmem(nf) ! go through this animal' si
                              iaa=isib(nf,mf)
                              if(iaa.ne.ia) then
                                 do j=1,nonzed(m,f)
                                    k=ntype(m,f,j)
                                    sum3(j)=pt(m,f,k)+freq(k,iaa)+work(k,ns,mf)
                                 enddo
                                 call LOGADD(sum3,nonzed(m,f))
                                 prod=prod+sum3(1)
                              endif
                           enddo ! line 4 calculated
                           sum2(ff)=ant(f,idd)+freq(f,idd)+dpost(f)+pt(m,f,i)+prod
 20                        continue
                        enddo
                        call LOGADD(sum2,ff)
                        sum1(mm)=ant(m,is)+freq(m,is)+tpost(m)+sum2(1)
 10                     continue
                     enddo
                     call LOGADD(sum1,mm)
                     ant(i,ia)=sum1(1)
                     temp(i)=sum1(1)
                  enddo
                  call LOGADD(temp,3)
                  do i=1,3
                     ant(i,ia)=ant(i,ia)-temp(1)
!                     print*, ia,i,ant(i,ia)
                  enddo
               enddo
            enddo
 50         continue
         enddo

! ----------------------------------------------------------------
! PEEL UP, IE CONDENSE INFO ON MATE AND PROGENY ONTO INDIVIDUAL
! START WITH YOUNGEST IN THE LIST
! ----------------------------------------------------------------
         call flippt

         do i=nobs,1,-1
            kc=ifirst(i)
 100        if(kc.eq.0) goto 200
            im=mate(kc)
!     collect posterior probability for "im" through all mates "is" of "im"
            kd=ifirst(im)
            is=mate(kd)
            do ii=1,3
               spost(ii)=0.
            enddo
           last=maxint
            do while (kd.ne.0)
               if(is.ne.last.and.is.ne.i) then
                  do ii=1,3
                     spost(ii)=spost(ii)+post(ii,kd)
                  enddo
               endif
               last=is
               kd=next(kd)
              is=mate(kd)
            enddo
!     pick up all the offspring of "i" and "im", say "idd" and
!     go through through the l mates "ms" of "idd", and calculate
!     prod-l post(idd,ms), store in work vector
            is=im
            kd=kc
            k=0
            do while (is.eq.im.and.kd.ne.0)
               k=k+1
               idd=prog(kd)

!if (k>HoldInt) then  ! k goes to 2*maxfs here
!  print*, "k",k
!  HoldInt=k
!endif
               nmem(k)=idd
               kj=ifirst(idd)
               ms=mate(kj)
               do ii=1,3
                  dpost(ii)=0.
               enddo
               last=maxint
               do while (kj.ne.0)        !Go through all mates of "idd"
                  if(ms.ne.last) then
                     do ii=1,3
                        dpost(ii)=dpost(ii)+post(ii,kj)
                     enddo
                  endif
                  last=ms
                  kj=next(kj)
                  ms=mate(kj)
               enddo
               do ii=1,3
                  term(ii,k)=dpost(ii)
               enddo
               kd=next(kd)
               is=mate(kd)
            enddo
!     NOW we are ready to calculate the posterior prob for match "i" and "im"
            do ii=1,3                   !for 3 g'type of "i"
               do j=1,3                 !for 3 g'types of "im"
                  prod=0.
                  do ns=1,k             ! through k offspring
                     idd=nmem(ns)
                     do m=1,nonzed(ii,j)
                        l=ntype(ii,j,m)
                        sum2(m)=pt(ii,j,l)+freq(l,idd)+term(l,ns)
                     enddo
                     call LOGADD(sum2,nonzed(ii,j))
                     prod=prod+sum2(1)
                  enddo
                  sum1(j)=ant(j,im)+freq(j,im)+spost(j)+prod
               enddo
               call LOGADD(sum1,3)
               temp(ii)=sum1(1)
               post(ii,kc)=sum1(1)
            enddo
            call LOGADD(temp,3)
            do ii=1,3
               post(ii,kc)=post(ii,kc)-temp(1)
!               print*,i,im,ii,post(ii,kc)
            enddo
            kc=kd
            goto 100
 200        continue
         enddo
! ------------------------------------------------------------------
!     NOW CALCULATE GENOTYPE PROBABILITIES
         sum1(1)=0.
         do i=1,nobs
            do ii=1,3
               spost(ii)=0.
            enddo
            kc=ifirst(i)
            if(kc.eq.0) goto 250
            im=mate(kc)
            last=maxint
            do while (kc.ne.0)
               if(im.ne.last) then
                  do ii=1,3
                     spost(ii)=spost(ii)+post(ii,kc)
                  enddo
               endif
               last=im
               kc=next(kc)
               im=mate(kc)
            enddo
 250        tsum=0.
            do ii=1,3
               spost(ii)=ant(ii,i)+freq(ii,i)+spost(ii)  ! just store it all in spost for convenience
            enddo
            maxvalspost=MAXVAL(spost)
            do ii=1,3
               spost(ii)=spost(ii)-maxvalspost + 20   ! bring them all up towards zero, then some (e^20 = 485 million = OK)
            enddo
            do ii=1,3
               if (spost(ii).LT.-100.) then  !e^-100  =~ 10^-44
                temp(ii)=0.
               else
                temp(ii)=exp(spost(ii))
               end if
               tsum=tsum+temp(ii)
            enddo
            pnor(i)=temp(1)/tsum
            prod=temp(2)/tsum
            sum1(1)=sum1(1)+abs(prod-phet(i))
            phet(i)=prod
            phom(i)=temp(3)/tsum
         enddo
         sum1(1)=sum1(1)/nobs
!         write(*,'(f13.7)') sum1(1)
         call flippt
         if(sum1(1).le.StopCrit) then
          itersUsed=kl
          goto 300   ! eg 0.0000001
         endif
      enddo
 300  continue

s0=0.
s1=0.
s2=0.
n=0
do i=1,nobs

  if (ABS(1.-pnor(i)-phet(i)-phom(i)) > .0000001) then
   PRINT*, 'Error: ',ABS(1.-pnor(i)-phet(i)-phom(i)),pnor(i),phet(i),phom(i)
  end if

  if(p1(i).eq.0 .and. p2(i).eq.0)then
     s0=s0+pnor(i)
     s1=s1+phet(i)
     s2=s2+phom(i)
     n=n+1
  endif
enddo


phethw=2*pprior*qprior ! do these here before reset pprior so GPI comes from freq used to get probabilities
phomhw=pprior*pprior

pResult(ifreq_iterate) = (s1+2*s2)/(2*(s0+s1+s2))
   pDev(ifreq_iterate) = pResult(ifreq_iterate) - pHold(ifreq_iterate)

if(pDev(ifreq_iterate) > 0.0 .and. pDev(ifreq_iterate) < LeastPositiveValue ) LeastPositive = ifreq_iterate
if(pDev(ifreq_iterate) < 0.0 .and. pDev(ifreq_iterate) > LeastNegativeValue ) LeastNegative = ifreq_iterate

!print'(i4,f12.8,a8)', 0,pDev(0),'zero2'
!print'(i4,f12.8)', ifreq_iterate,pDev(ifreq_iterate)

!Print'(i6,3f16.8)', ifreq_iterate, pHold(ifreq_iterate), pResult(ifreq_iterate), pDev(ifreq_iterate)

    if (ifreq_iterate==0) then
! do nothing here
elseif (ifreq_iterate==1) then
! 2-point regression
    breg = (pDev(1) - pDev(0)) / (pHold(1) - pHold(0))
    areg = pDev(1) - breg*pHold(1)
    pHold(ifreq_iterate+1)= -1*areg/breg
!    print'(4f12.8)', pDev(1) , pDev(0),  pHold(1) , pHold(0)
!    print'(4f12.8)', (pDev(1) - pDev(0)) / (pHold(1) - pHold(0))
!    print'(2i4,3f12.8)', 2, ifreq_iterate+1, areg, breg, pHold(ifreq_iterate+1)
	if (pDev(0)<0.  .AND. pDev(1)>0. ) then
           nfreq_max=ifreq_iterate+1 ! Blowing out already - get sensible midpoint
!!           PRINT*, 'Diverging solution - taking 2-point regression'
    endif
elseif (ifreq_iterate>1) then
! j-point regression
  j=MIN(maxRegpoints,ifreq_iterate+1)  ! 3 point seems fastest.  More points could be more robust.
        meanX=0
        meanY=0
        sumY =0
        sumXY=0
        sumX =0
        sumX2=0
        do i = ifreq_iterate-(j-1), ifreq_iterate
                meanX = meanX + pHold(i)
                meanY = meanY + pDev(i)
                sumY  = sumY  + pDev(i)
                sumX  = sumX  + pHold(i)
                sumX2 = sumX2 + pHold(i)*pHold(i)
                sumXY = sumXY + pHold(i)*pDev(i)
        end do
        breg  = (sumXY - sumX*sumY/j) / (sumX2 - sumX*sumX/j)
        meanX = meanX/j
        meanY = meanY/j
        areg  = meanY - breg*meanX
        pHold(ifreq_iterate+1)= -1*areg/breg
        !print'(2i4,3f12.8)', j, ifreq_iterate+1, areg, breg, pHold(ifreq_iterate+1)

	if(pHold(ifreq_iterate+1) < pHold(LeastPositive) .or. pHold(ifreq_iterate+1) > pHold(LeastNegative)) then ! outside the best bound, probably curlivilear.
	  pHold(ifreq_iterate+1) = pHold(LeastPositive) + (pHold(LeastNegative)-pHold(LeastPositive))&
	  *   pDev(LeastPositive)/(pDev(LeastPositive)-pDev(LeastNegative))
	  !print*,' '
	  !print'(2i3,ff12.8)', LeastPositive , LeastNegative,    pDev(LeastPositive)/(pDev(LeastPositive)-pDev(LeastNegative)), pHold(ifreq_iterate+1)
	endif

end if

if (ifreq_iterate>=1) then
   !IF(pHold(ifreq_iterate+1) > 0.9999) pHold(ifreq_iterate+1) = 0.9999
   !IF(pHold(ifreq_iterate+1) < 0.0001) pHold(ifreq_iterate+1) = 0.0001
   if(pHold(ifreq_iterate+1) > 0.9999 .OR. pHold(ifreq_iterate+1) < 0.0001) then
     pHold(ifreq_iterate+1) = pHold(ifreq_iterate) ! take last value
     nfreq_max=ifreq_iterate+1 ! and stop after getting probs from that (redundant - previous line makes a StopCrit stop)
!!     PRINT*, 'Unstable projection - taking current prior ...'
   endif
endif


!!IF(nfreq_max>0) PRINT'(2i8,5f12.8)', ifreq_iterate, itersused, pHold(ifreq_iterate), pResult(ifreq_iterate), pDev(ifreq_iterate), pHold(ifreq_iterate+1)

!print'(2f20.15)', pHold(ifreq_iterate)-pHold(ifreq_iterate+1) , StopCrit

if(nfreq_max>0) then
 IF ( ABS(pHold(ifreq_iterate)-pHold(ifreq_iterate+1)) < StopCrit ) ifreq_iterate=nfreq_max
endif

ierrors=1
IF ( ifreq_iterate==nfreq_max) then
!!       open(7,file=outfile,status='unknown')
!!       	if(Imprinting==0) Write(7,*)'ID SireID DamID Tag p(11) p(Het) p(22) Phen Index'
!!     		if(Imprinting==1) Write(7,*)'ID SireID DamID Tag p(11) p(12) p(21) p(22) Phen Index'
!!       	if(Imprinting==2) Write(7,*)'ID SireID DamID Tag p(11) p(het) p(21)-p(12) p(22) Phen Index'
       pnor(0)=1.- phethw - phomhw
       phet(0)=phethw
       phom(0)=phomhw

	if (LimitAnimals==1 .and. nobs>LimitNumber) then
	  PRINT*, 'Results for only', LimitNumber, ' animals will be written.'
	  nwritten=LimitNumber
	else
	  nwritten=nobs
	endif


       do i=1,nwritten
			call info(phet(i), phom(i), phethw, phomhw, probindex)

			iflag=0
			if (ABS(1. - pnor(i)-phet(i)-phom(i)) > .0001 ) iflag=iflag+10

			do j=1, phenotypes
				if (phen(i)==phenotype(j)) then
					if (g(j,0)<0.000001 .and. pnor(i)>0.01 ) iflag=iflag+1
					if (g(j,1)<0.000001 .and. phet(i)>0.01 ) iflag=iflag+1
					if (g(j,2)<0.000001 .and. phom(i)>0.01 ) iflag=iflag+1
				endif
			enddo

			!  Only for identity penetrance matrix ....
			!			if (phen(i)==0 .and. ABS(1. - pnor(i))  > .0001 ) iflag=iflag+1
			!			if (phen(i)==1 .and. ABS(1. - phet(i))  > .0001 ) iflag=iflag+1
			!			if (phen(i)==2 .and. ABS(1. - phom(i))  > .0001 ) iflag=iflag+1
			! for 0, 1 condition only here  1 (affected) = phenotype 2  ....
			!ProbFit = ABS( phen(i) - ( g(2, 0)*pnor(i) + g(2, 1)*phet(i) + g(2, 2)*phom(i) ))
			!ProbFit = 1 - Probfit ! Fit on 0 to 1 scale


			if (ierrors*iflag > 0 ) then
				open(8,file='geneprob_err.txt',status='unknown')
				Write(8,*)'ID SireID DamID Tag p(11) p(Het) p(22) Phen Index Error'
!!				PRINT*,' Errors found.  See Geneprob_err.ped  *********************** Errors '
				ierrors = -1
			end if
			if (iflag/=0) then
			 write(8,56) i,p1(i),p2(i),'"',id(i),'"',pnor(i),phet(i),phom(i),phen(i),probindex,'Error ', iflag
			end if

       	if(Imprinting>0) then
				if(phet(i)<0.0000001) then
!!						write(7,551) i,p1(i),p2(i),'"',id(i),'"',pnor(i),phet(i),phet(i),phom(i),phen(i),probindex
						Probs00(i,currentSnp) = pnor(i)
						Probs01(i,currentSnp) = phet(i)
						Probs10(i,currentSnp) = phet(i)
						Probs11(i,currentSnp) = phom(i)
						GPI(i,currentSnp) = probindex
				else
					p12= (pnor(seqsire(i))+0.5*phet(seqsire(i))) * (phom( seqdam(i))+0.5*phet( seqdam(i)))  ! extra safe due to the above
					p21= (pnor( seqdam(i))+0.5*phet( seqdam(i))) * (phom(seqsire(i))+0.5*phet(seqsire(i)))
					if(p12+p21>0.00000001) then
						IMPratio= p12 / (p12+p21)
					else
						IMPratio= 0.5 ! neutral but should not be invked anyway
					endif


!					IMPratio=                       (pnor(seqsire(i))+0.5*phet(seqsire(i))) * (phom( seqdam(i))+0.5*phet( seqdam(i)))
!					IMPratio= IMPratio / (IMPratio+ (pnor( seqdam(i))+0.5*phet( seqdam(i))) * (phom(seqsire(i))+0.5*phet(seqsire(i))))
					if(Imprinting==2)then
!!						write(7,552) i,p1(i),p2(i),'"',id(i),'"',pnor(i),   phet(i), (1.-2.*IMPratio)*phet(i),  phom(i),phen(i),probindex
						Probs00(i,currentSnp) = pnor(i)
						Probs01(i,currentSnp) = phet(i)
						Probs10(i,currentSnp) = (1.-2.*IMPratio)*phet(i)
						Probs11(i,currentSnp) = phom(i)
						GPI(i,currentSnp) = probindex
					else
!!						write(7,551) i,p1(i),p2(i),'"',id(i),'"',pnor(i),IMPratio*phet(i),(1.-IMPratio)*phet(i),phom(i),phen(i),probindex
						Probs00(i,currentSnp) = pnor(i)
						Probs01(i,currentSnp) = IMPratio*phet(i)
						Probs10(i,currentSnp) = (1.-IMPratio)*phet(i)
						Probs11(i,currentSnp) = phom(i)
						GPI(i,currentSnp) = probindex

					endif
	       	endif
       	else
!!				   write(7,55 ) i,p1(i),p2(i),'"',id(i),'"',pnor(i),             phet(i),                  phom(i),phen(i),probindex
						Probs00(i,currentSnp) = pnor(i)
						Probs01(i,currentSnp) = phet(i)
						Probs10(i,currentSnp) = phet(i)
						Probs11(i,currentSnp) = phom(i)
						GPI(i,currentSnp) = probindex

			endif

      enddo

      !!*****************************************************
      !!temp write geneprobs to file for testing
      !!open (UNIT=102,FILE="GenProbs.txt",STATUS="unknown")
      !!do i=1,nwritten
		!!write(102,*) seqid(i), Probs0(i,:)
		!!write(102,*) seqid(i), Probs1(i,:)
		!!write(102,*) seqid(i), Probs2(i,:)
      !!end do
      !!close(102)
      !!*******************************************************

 55   format(3i8,1x,a1,a16,a1,3f9.6,i2,f6.1)  ! make id field long enough!
 551   format(3i8,1x,a1,a16,a1,4f9.6,i2,f6.1)  ! make id field long enough!
 552   format(3i8,1x,a1,a16,a1,2f9.6,f9.5,f9.6,i2,f6.1)  ! make id field long enough!
 56   format(3i8,1x,a1,a16,a1,3f9.6,i2,f6.1,a15,i3)  ! make id field long enough!
      CLOSE(7)
!!      print*,' '
!      PRINT'(a43,f14.7)', 'Frequency estimated from base individuals:', pprior  ! The value used for the set of probabilities written.
	OutputMaf(currentSnp)=pprior
      call system_clock(iticks2,i2,i3)
      !IF(iticks2.lt.iticks1)iticks2=iticks2+i3
!!      print'(a43,f14.7)','Number of seconds for this run:',(float(iticks2-iticks1))/i2
!!      print*,'Results are in: ', outfile
endif


ENDDO ! ifreq_iterate

close (7)
if (ierrors==-1) close (8)


!if (PauseAtEnd==1) Pause
print *,'done with geneprob', currentSnp  !!*******************************************************!!

!deallocate(phenhold)
deALLOCATE(MATE,NEXT,IFIRST,POST,phen,prog,&
         p1,p2,ant,phom,phet, &
         freq,pnor)
!deALLOCATE(nmem,isib,damGP,work,term)
!print *,'Start'
deallocate(nmem)
!print *,'**1'
deallocate(isib)
!print *,'***2'
deallocate(damGP)
!print *,'***3'
deallocate(work)
!print *, '***4'
deallocate(term)
!print *,'***5'

deALLOCATE (pHold, pResult, pDev)

end subroutine geneprob


SUBROUTINE LNKLST(I,J,NP,IFLAG)
      use common_GP

!      integer prog(0:2*MXANEQ)

      IF (I.le.0.or.J.le.0) RETURN

!     KEEP MAXIMUM EQUATION
      IF (I.GT.MXEQ)THEN
         MXEQ=I
         IF (I.GT.MXANEQ)STOP ' TOO MANY EQUATIONS - RECOMPILE'
      ENDIF

      IR=0
      IP=IFIRST(I)
      K=MATE(IP)
 10   IF (K.LE.J.AND.IP.NE.0) THEN
!      IF(np.le.30) print '(6i9)',np,k,j,ip,next(ip),mate(next(ip))
         IR=IP
         IP=NEXT(IP)
         K=MATE(IP)
         GOTO 10
      ENDIF

      IF(IFLAG.EQ.1) THEN
         MM=MM+1
         NUSED=MM
      ELSE
         NN=NN+1
         NUSED=NN
      ENDIF

      MATE(NUSED)=J
      PROG(NUSED)=NP
      IF (IP.EQ.0) THEN
         IF (IFIRST(I).EQ.0) THEN
            IFIRST(I)=NUSED
         ELSE
            NEXT(IR)=NUSED
         ENDIF
      ELSEIF (IR.EQ.0) THEN
         NEXT(NUSED)=IFIRST(I)
         IFIRST(I)=NUSED
      ELSE
         NEXT(NUSED)=NEXT(IR)
         NEXT(IR)=NUSED
      ENDIF

      END

!   ********************************************************************

SUBROUTINE FLIPPT
      use common_GP
!      integer prog(0:2*mxaneq)

      DO 40 K=1,MXEQ
         KC=IFIRST(K)
         NEWPT=0
 20      IF (KC.EQ.0) GOTO 30
         IOLDPT=NEXT(KC)
         NEXT(KC)=NEWPT
         NEWPT=KC
         KC=IOLDPT
         IF(KC.NE.0.AND.MATE(KC).EQ.MATE(NEWPT)) THEN
            DO I=1,3
               POST(I,IOLDPT)=POST(I,NEWPT)
            ENDDO
         ENDIF
         GOTO 20
 30      IFIRST(K)=NEWPT
 40   CONTINUE
      RETURN
      END

!  *********************************************************************

SUBROUTINE LOGADD(summ,n)
      real*8 summ(3),add,t

      if(n.eq.1) return
      if(n.eq.2) then
         summ(1)=add(summ(1),summ(2))
      else
         do i=1,2
            mm=i
            do j=i+1,3
               if(summ(j).lt.summ(mm)) then
                  t=summ(mm)
                  summ(mm)=summ(j)
                  summ(j)=t
               endif
            enddo
         enddo
         summ(1)=add(summ(1),summ(2))
         summ(1)=add(summ(1),summ(3))
      endif
      return
      end

! *****************************************************************
function add(x1,x2)
      real*8 add,x1,x2,diff,expmax
      expmax=300.d0
      diff=x1-x2
      if(diff.gt.expmax) then
         add=x1
      elseif(-diff.gt.expmax) then
         add=x2
      else
         add=x2 + dlog(dexp(diff) + 1.d0)
      endif
      return
      end


Subroutine info(phet, phom, phethw, phomhw, probindex)
	implicit none
	Real(Kind=8)	:: phet, phom, phethw, phomhw, probindex
	Real(Kind=8)	:: r, y, x, pi, d30, yobs, xobs, hwy, hwx, xi, xj, Bangle, a, b, c, xy

	r = 1.
	y = .5 * r
	x = sqrt(.75) * r
	pi = 3.141593
	d30 = 30 * pi / 180

	yobs = -y + phet * (r + y)
	xobs = -x + (phet * TAN(d30) + phom / COS(d30)) * (r + y)
	hwy = -y + phethw * (r + y)
	hwx = -x + (phethw * TAN(d30) + phomhw / COS(d30)) * (r + y)
	xi = (xobs + yobs * hwy / hwx) / (1 + (hwy / hwx) ** 2)
	xj = xi * hwy / hwx

	IF(abs(hwx - xobs).lt.0.00001.and.abs(hwy - yobs).lt.0.00001) THEN
	 probindex = 0.
	 return
	endif

	Bangle = atan(sqrt((xobs - xi) ** 2 + (yobs - xj) ** 2) / sqrt((hwx - xi) ** 2 + (hwy - xj) ** 2))
	IF(xi/hwx.gt.1.) Bangle = pi - Bangle

	a = sqrt(hwx ** 2 + hwy ** 2)
	b = r
	c = a * COS(Bangle) +  sqrt((a * COS(Bangle)) ** 2 - (a ** 2 - b ** 2))
	 !c2 = a * COS(Bangle) - sqrt((a * COS(Bangle)) ** 2 - (a ** 2 - b ** 2))
	 !c2 is other form of solution to quadratic formula - first seems always OK
	xy = sqrt((xobs - hwx) ** 2 + (yobs - hwy) ** 2)

	probindex = 100 * xy / c

END

! ********************************************************************

subroutine PVseq(mode)

USE commonbits

implicit none

character (LEN=lengan), ALLOCATABLE :: holdsireid(:), holddamid(:)
character (LEN=lengan), ALLOCATABLE :: holdid(:), SortedId(:), SortedSire(:), SortedDam(:)
character (LEN=lengan)              :: IDhold

integer, ALLOCATABLE                :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
integer, ALLOCATABLE                :: OldN(:), NewN(:), holdsire(:), holddam(:)

INTEGER :: mode    ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
INTEGER :: i, j, k, kk, newid, itth, itho, ihun, iten, iunit
integer :: nsires, ndams, newsires, newdams, nbisexuals, flag
INTEGER :: ns, nd, iextra, oldnobs, kn, kb, oldkn, ks, kd
INTEGER :: Noffset, Limit, Switch, ihold, ipoint

do j = 1, nobs
  If (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') Then
    dam(j) = '0'
    seqdam(j)=0
  endif
  If (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') Then
    sire(j) = '0'
    seqsire(j)=0
  endif
enddo !j

if(mode.eq.1) then
 PRINT*,  ' Inserting dummy IDs ... '
 newid=0
 do j = 1, nobs
   if(((sire(j) == '0').and.(dam(j).ne.'0'))  .or. ((sire(j).ne.'0').and.(dam(j) == '0'))) then
         newid=newid+1
         if(newid.gt.99999) then
           PRINT*, newid, ' ...'
           stop 'too many dummy single parent IDs'
         endif
         itth=int(newid/10000)
         itho=int(newid/1000)-10*itth
         ihun=int(newid/100)-10*itho-100*itth
         iten=int(newid/10)-10*ihun-100*itho-1000*itth
         iunit=newid-10*iten-100*ihun-1000*itho-10000*itth
         if(sire(j) == '0') sire(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
         if( dam(j) == '0')  dam(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
   endif
 enddo
endif

PRINT*,  ' Sorting Sires ... '

ALLOCATE  (SortedId(nobs), SortedIdIndex(nobs))

SortedId(1:nobs) = Sire(1:nobs)

  Noffset = INT(nobs/2)
  DO WHILE (Noffset>0)
      Limit = nobs - Noffset
      switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
               IDhold=SortedId(i)
               SortedId(i)=SortedId(i + Noffset)
               SortedId(i + Noffset)=IDhold

               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
  enddo

nsires=0
IF(SortedId(1) /= '0') nsires=1
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nsires=nsires+1
end do

ALLOCATE  (SortedSire(0:nsires), SortedSireIndex(nsires))
SortedSire(0) = '0'

nsires=0
IF(SortedId(1) /= '0') THEN
 nsires=1
 SortedSire(1) = SortedId(1)
ENDIF
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
   nsires=nsires+1
   SortedSire(nsires) = SortedId(i)
  ENDIF
end do

PRINT*,  ' Sorting Dams ... '

SortedId(1:nobs) = Dam(1:nobs)

  Noffset = INT(nobs/2)
  DO WHILE (Noffset>0)
      Limit = nobs - Noffset
      switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
               IDhold=SortedId(i)
               SortedId(i)=SortedId(i + Noffset)
               SortedId(i + Noffset)=IDhold

               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
  enddo

nDams=0
IF(SortedId(1) /= '0') nDams=1
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nDams=nDams+1
end do

ALLOCATE  (SortedDam(0:nDams), SortedDamIndex(ndams))
SortedDam(0)='0'

nDams=0
IF(SortedId(1) /= '0') THEN
 nDams=1
 SortedDam(1) = SortedId(1)
ENDIF
do i=2,nobs
  IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
   nDams=nDams+1
   SortedDam(nDams) = SortedId(i)
  ENDIF
end do


PRINT*,  ' Sorting IDs ... '

SortedId(1:nobs) = ID(1:nobs)
do i=1,nobs
 SortedIdIndex(i) = i
end do

  Noffset = INT(nobs/2)
  DO WHILE (Noffset>0)
      Limit = nobs - Noffset
      switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
               IDhold=SortedId(i)
               SortedId(i)=SortedId(i + Noffset)
               SortedId(i + Noffset)=IDhold

               ihold=SortedIdIndex(i)
               SortedIdIndex(i)=SortedIdIndex(i + Noffset)
               SortedIdIndex(i + Noffset)=ihold

               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
  enddo

PRINT*,  ' Check for duplicate IDs ... '
flag = -1
Do i = 2, nobs
  If (SortedID(i) == SortedID(i - 1)) Then
   If (flag == -1) Then
     open (1,FILE='ID_err.txt',STATUS = 'unknown')
     WRITE(1,*) 'Duplicated IDs ...'
     flag = 0
   End If
   WRITE(1,*) SortedID(i)
   flag = flag + 1
  End If
enddo

 If (flag > -1) Then
  Close (1)
  PRINT*, flag,' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
 End If

PRINT*,  ' Males ... '
PRINT*,  '  Find or set sire indices ... '

newsires = 0

do j=1,nsires

! check if already listed as an individual
   ipoint=INT(nobs/2)
   Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (SortedSire(j).lt.SortedId(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo

    kn=0
    if (SortedSire(j)==SortedId(ipoint)) kn=1
    do while (ipoint<nobs .and. kn==0 .and. SortedSire(j) > SortedId(ipoint))
     ipoint=ipoint+1
    enddo
    if (SortedSire(j)==SortedId(ipoint)) kn=1
    do while (ipoint>1 .and. kn==0 .and. SortedSire(j) < SortedId(ipoint))
     ipoint=ipoint-1
    enddo
    if (SortedSire(j)==SortedId(ipoint)) kn=1

    IF(kn==1) then
     SortedSireIndex(j) = SortedIdIndex(ipoint)
    else    ! sire is unlisted base sire
     newsires = newsires + 1
     SortedSireIndex(j) = nobs + newsires ! for now
    endif
end do !j

 ALLOCATE  (holdsireid(newsires))
 kn=0
 do j=1,nsires
  if (SortedSireIndex(j) > nobs) then
   kn=kn+1
   holdsireid(SortedSireIndex(j)-nobs) = SortedSire(j)
  end if
 enddo
 IF(kn /= newsires) stop'newsires error'


PRINT*,  '  Find seqsire ... '

do j = 1, nobs

  If (sire(j) == '0') Then
    seqsire(j)=0
  else

    ipoint=INT(nsires/2)
    Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (Sire(j).lt.SortedSire(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo

    kn=0
    if (Sire(j)==SortedSire(ipoint)) kn=1
    do while (ipoint<nsires .and. kn==0 .and. Sire(j) > SortedSire(ipoint))
     ipoint=ipoint+1
    enddo
    if (Sire(j)==SortedSire(ipoint)) kn=1
    do while (ipoint>1 .and. kn==0 .and. Sire(j) < SortedSire(ipoint))
     ipoint=ipoint-1
    enddo
    if (Sire(j)==SortedSire(ipoint)) kn=1
    IF(kn==1) then
     seqsire(j) = SortedSireIndex(ipoint)
    else
     PRINT*, ' Error: Sire missing: ', Sire(j)
     stop
    endif

  endif

ENDDO !j

PRINT*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'

PRINT*,  ' Females ... '
PRINT*,  '  Find or set dam indices ... '

newdams = 0
nbisexuals = 0

do j=1,ndams

! check if already listed as an individual
   ipoint=INT(nobs/2)
   Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (Sorteddam(j).lt.SortedId(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo

    kn=0
    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint  ! store ipoint here as ipoint can change with bisexuals
    do while (ipoint<nobs .and. kn==0 .and. Sorteddam(j) > SortedId(ipoint))
     ipoint=ipoint+1
    enddo
    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
    do while (ipoint>1 .and. kn==0 .and. Sorteddam(j) < SortedId(ipoint))
     ipoint=ipoint-1
    enddo
    if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
! check if already listed as a sire (and therefore bisexual)
   ipoint=INT(nsires/2)
   Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (SortedDam(j).lt.SortedSire(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo

    kb=0
    if (SortedDam(j)==SortedSire(ipoint)) kb=1
    do while (ipoint<nsires .and. kb==0 .and. SortedDam(j) > SortedSire(ipoint))
     ipoint=ipoint+1
    enddo
    if (SortedDam(j)==SortedSire(ipoint)) kb=1
    do while (ipoint>1 .and. kb==0 .and. SortedDam(j) < SortedSire(ipoint))
     ipoint=ipoint-1
    enddo
    if (SortedDam(j)==SortedSire(ipoint)) kb=1

    IF(kb==1) then
      nbisexuals = nbisexuals + 1
      open (1,FILE='bisex.txt',position = 'append')
       WRITE(1,*) SortedDam(j)
      close(1)
    endif

    if (kb==1) then
     SorteddamIndex(j) = SortedSireIndex(ipoint)
    elseif (kn>=1) then
     SorteddamIndex(j) = SortedIdIndex(kn)
    else    ! dam is unlisted base dam
     newdams = newdams + 1
     SorteddamIndex(j) = nobs + newsires + newdams ! for now
    endif

end do !j

If (nbisexuals > 0)  PRINT*, nbisexuals,' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'

 ALLOCATE  (holddamid(newdams))
 kn=0
 do j=1,ndams
  if (SortedDamIndex(j) > nobs+newsires) then
   kn=kn+1
   holddamid(SortedDamIndex(j)-nobs-newsires) = SortedDam(j)
  end if
 enddo
 IF(kn /= newdams) stop'newdams error'



PRINT*,  '  Find seqdam ... '

do j = 1, nobs

  If (dam(j) == '0') Then
    seqdam(j)=0
  else

    ipoint=INT(ndams/2)
    Noffset = INT(ipoint/2)
   do while (Noffset>1)
    IF (dam(j).lt.Sorteddam(ipoint)) THEN
     ipoint = ipoint - Noffset
     Noffset = INT(Noffset/2)
    else
     ipoint = ipoint + Noffset
     Noffset = INT(Noffset/2)
    endif
   enddo

    kn=0
    if (dam(j)==Sorteddam(ipoint)) kn=1
    do while (ipoint<ndams .and. kn==0 .and. dam(j) > Sorteddam(ipoint))
     ipoint=ipoint+1
    enddo
    if (dam(j)==Sorteddam(ipoint)) kn=1
    do while (ipoint>1 .and. kn==0 .and. dam(j) < Sorteddam(ipoint))
     ipoint=ipoint-1
    enddo
    if (dam(j)==Sorteddam(ipoint)) kn=1
    IF(kn==1) then
     seqdam(j) = SorteddamIndex(ipoint)
    else
     PRINT*, ' Error: dam missing: ', dam(j)
     stop
    endif

  endif

ENDDO !j

PRINT*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'

PRINT*,  ' Arranging unlisted base parents ... '


iextra = newsires + newdams

If (iextra > 0) then
      PRINT*, ' ', iextra, ' unlisted base parents found.'

 ! SortedId and SortedIdIndex just used as a holder while redimensioning

 SortedId(1:nobs)=id(1:nobs)
 deallocate (id)
 ALLOCATE(id(nobs+iextra))
 id(1+iextra:nobs+iextra)=SortedId(1:nobs)

 SortedId(1:nobs)=sire(1:nobs)
 deallocate (sire)
 ALLOCATE(sire(nobs+iextra))
 sire(1+iextra:nobs+iextra)=SortedId(1:nobs)

 SortedId(1:nobs)=dam(1:nobs)
 deallocate (dam)
 ALLOCATE(dam(nobs+iextra))
 dam(1+iextra:nobs+iextra)=SortedId(1:nobs)

 SortedIdIndex(1:nobs)=seqsire(1:nobs)
 deallocate (seqsire)
 ALLOCATE(seqsire(nobs+iextra))
 seqsire(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)

 SortedIdIndex(1:nobs)=seqdam(1:nobs)
 deallocate (seqdam)
 ALLOCATE(seqdam(nobs+iextra))
 seqdam(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)

endif

PRINT*, ' Inserting unlisted base parents ...'

oldnobs = nobs
nobs = nobs + iextra
PRINT*, ' Total number of animals = ',nobs

ALLOCATE (passedorder(nobs))
passedorder=0

do i = 1+iextra, nobs
 passedorder(i)= i-iextra

 If (sire(i) == '0')then
   seqsire(i) = 0
 Else
   seqsire(i) = iextra + seqsire(i)
   If (seqsire(i) > nobs)  seqsire(i) = seqsire(i) - nobs  ! for unlisted sires
 End If

 If (dam(i) == '0') Then
   seqdam(i) = 0
  Else
   seqdam(i) = iextra + seqdam(i)
   If (seqdam(i) > nobs)  seqdam(i) = seqdam(i) - nobs
  End If
ENDDO !i

do i = 1, newsires
 ID(i) = holdsireid(i)
 passedorder(i)=0
 seqsire(i) = 0
 seqdam(i) = 0
ENDDO !i
do i = newsires + 1, newsires + newdams
 ID(i) = holddamid(i - newsires)
 passedorder(i)=0
 seqsire(i) = 0
 seqdam(i) = 0
ENDDO !i

DEALLOCATE(holdsireid, holddamid, SortedIdIndex, SortedId)

flag = 0
Do i = 1, nobs
If (i <= seqsire(i) .Or. i <= seqdam(i) ) flag = 1
enddo !i
If (flag == 0) return  !PRINT*, 'not needed'!return

PRINT*, ' Re-Ordering pedigree ...'


Allocate ( OldN(0:nobs), NewN(0:nobs) )
ALLOCATE ( holdid(0:nobs), holdsire(nobs), holddam(nobs) )

OldN(0) = 0
NewN=0
!seqsire(0) = 0 !not needed !
!seqdam(0) = 0

holdid(1:nobs) = ID(1:nobs)
holdsire = seqsire
holddam = seqdam

!Find base ancestors ...
kn = 0
do i = 1, nobs
 If (seqsire(i) == 0 .And. seqdam(i) == 0) Then
      kn = kn + 1
      NewN(i) = kn
      OldN(kn) = i
 End If
ENDDO !i

!Re-order pedigree ...
NewN(0) = nobs + 1
flag = 0
Do While (kn < nobs)
 oldkn = kn
 do i = 1, nobs
  If (NewN(i) == 0) Then !And ID(i) <> 'UniqueNULL' Then
    Ks = seqsire(i)
    Kd = seqdam(i)
    If (NewN(Ks) > 0 .And. NewN(Kd) > 0) Then
      kn = kn + 1
      NewN(i) = kn
      OldN(kn) = i
    End If
  End If
 enddo !i

 ! to avoid hang on unexpected problem ...
 If (kn == oldkn) Then
  flag = flag + 1
 Else
  flag = 0
 endif

 If (flag > 10) Then
   open(1,file='ped_err.txt',status='unknown')
   write(1,*) 'Pedigree errors found involving two or more of the following relationships ...'
   write(1,*)
   write(1,*) '       Index numbers are followed by names.'
   write(1,*) '       Index number 0 means unknown, whence name is blank.'
   write(1,*)
   do i = 1, nobs
    If (NewN(i) == 0) Then
     write(1,*) 'Individual:',          i, ':  ', ID(i)
     write(1,*) '    Father:', seqsire(i), ':  ', ID(seqsire(i))
     write(1,*) '    Mother:',  seqdam(i), ':  ', ID(seqdam(i))
     write(1,*)
    End If
   ENDDO !i
   Close (1)
   PRINT*,  'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
   stop
 End If
ENDDO !while

NewN(0) = 0

do i = 1, nobs
 ID(i) = holdid(OldN(i))
enddo

do i = 1, nobs
seqsire(i) = NewN(holdsire(OldN(i)))
seqdam(i) = NewN(holddam(OldN(i)))
If (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
   PRINT*,  'out of order'
   stop
endif
ENDDO !i

DO i = 1, nobs
  holdsire(i) = passedorder(i)  ! holdsire just because it is free
enddo
DO i = 1, nobs
  passedorder(i) = holdsire(OldN(i))
enddo

deallocate ( OldN, NewN, holdid, holdsire, holddam) ! holdrec)

!do i = 1, nobs
! PRINT'(3i5,2x,3a4,i5)', i, seqsire(i), seqdam(i), id(i), sire(i), dam(i), passedorder(i)
!enddo

end subroutine



subroutine GetMaxFS!!(maxfs, maxmates, nfamilies)

USE commonbits
implicit none

INTEGER (KIND= 8), allocatable :: family(:)       ! max value is 9,223,372,036,854,775,807 allowing for plenty of space
INTEGER (KIND= 8)              :: holdfamily, multiplier
!!INTEGER :: maxfs, maxmates, maxmates1, maxmates2, nfamilies
INTEGER :: maxmates1, maxmates2
INTEGER :: i, Noffset, Limit, Switch, fs, mates, parent1, parent2, oldparent1, oldparent2

ALLOCATE (family(nobs))

PRINT*,  ' Finding maximum FS family size and maximum mates ... '

holdfamily=0
multiplier = 100000000   ! allows for up to 99,999,999 in number system.

do i=1,nobs
 family(i) = multiplier * seqsire(i) + seqdam(i)
! IF(family(i) /= 0) PRINT'(3i7,i15)', i, seqsire(i), seqdam(i), family(i)
end do

  Noffset = INT(nobs/2)
  DO WHILE (Noffset>0)
      Limit = nobs - Noffset
      switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (family(i).gt.family(i + Noffset)) THEN
               holdfamily=family(i)
               family(i)=family(i + Noffset)
               family(i + Noffset)=holdfamily

               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
  enddo

nfamilies=0
do i=2,nobs
  IF(family(i) /= family(i-1)) then
    nfamilies=nfamilies+1
  endif
end do


maxfs = 0
fs=1

do i=2,nobs
  parent1=INT(family(i)/multiplier)
  parent2=family(i)-multiplier*parent1
 IF(parent1 /= 0 .and. parent2 /= 0) then  ! note how this is handled
  IF(family(i) == family(i-1)) then
    fs = fs + 1
  else
    IF(fs > maxfs) then
     maxfs = fs
     holdfamily=family(i-1)
    endif
    fs = 1
  endif
 endif
end do

parent1=INT(holdfamily/multiplier)
parent2=holdfamily-multiplier*parent1
!PRINT*, '  Maximum FS family size ... ',maxfs, ' for parents: ',parent1,' ',parent2
PRINT*, '  Maximum FS family size ... ',maxfs, ' for parents: ',id(parent1),' ',id(parent2)

Maxmates1 = 0
mates=1
oldparent1=0
oldparent2=0

do i=2,nobs

  parent1=INT(family(i)/multiplier)
 IF(parent1 /= 0) then
  parent2=family(i)-multiplier*parent1
   IF(parent1 == oldparent1) then
     IF(parent2 /= oldparent2) mates=mates+1
   else
    IF(mates > maxmates1)then
     maxmates1 = mates
     Limit=oldparent1  ! to record the one with most mates
    endif
    mates=1
   endif
	If (i==nobs .and. mates > maxmates1) Then  ! need to cover the last observation
		maxmates1 = mates
		Limit = oldparent1 ! to record the one with most mates
	End If
   oldparent1=parent1
   oldparent2=parent2
 endif
enddo

if (Limit > 0) then
 PRINT*, '  Max female mates of males=', maxmates1, ', for male ID, seqID: ', id(Limit), Limit
else
 PRINT*, '  Max female mates of males=', maxmates1
endif


!Now max mates for females
do i=1,nobs
 family(i) = multiplier * seqdam(i) + seqsire(i)
end do

  Noffset = INT(nobs/2)
  DO WHILE (Noffset>0)
      Limit = nobs - Noffset
      switch=1
    DO WHILE (Switch.ne.0)
       Switch = 0
       do i = 1, Limit
          IF (family(i).gt.family(i + Noffset)) THEN
               holdfamily=family(i)
               family(i)=family(i + Noffset)
               family(i + Noffset)=holdfamily

               Switch = i
          endif
       enddo
       Limit = Switch - Noffset
    enddo
    Noffset = INT(Noffset/2)
  enddo


Maxmates2 = 0
mates=1
oldparent1=0
oldparent2=0

do i=2,nobs
  parent1=INT(family(i)/multiplier)
 IF(parent1 /= 0) then
  parent2=family(i)-multiplier*parent1
   IF(parent1 == oldparent1) then
     IF(parent2 /= oldparent2) mates=mates+1
   else
    IF(mates > maxmates2)then
     maxmates2 = mates
     Limit=oldparent1  ! to record the one with most mates
    endif
    mates=1
   endif
   	If (i==nobs .and. mates > maxmates2) Then  ! need to cover the last observation
   		maxmates2 = mates
   		Limit = oldparent1 ! to record the one with most mates
   	End If
   oldparent1=parent1
   oldparent2=parent2

 endif
enddo

if (Limit > 0) then
 PRINT*, '  Max male mates of females=', maxmates2, ', for female ID, seqID: ', id(Limit), Limit
else
 PRINT*, '  Max male mates of females=', maxmates2
endif

maxmates = MAX(maxmates1,maxmates2)

end subroutine

subroutine mkdir(tmpdir)

	character(len=*) :: tmpdir

	open(unit=1000,file=".tmpsh",status="unknown")
	write(1000,*) "if [ ! -d "// trim(tmpdir) //" ]"
	write(1000,*) "then mkdir " // trim(tmpdir)
	write(1000,*) "fi"
	close(1000)

	call system("chmod a+x .tmpsh")
	call system("./.tmpsh")
	call system("rm .tmpsh")

end subroutine mkdir




