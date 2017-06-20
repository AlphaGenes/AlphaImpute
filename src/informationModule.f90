! #ifndef _WIN32

! #define STRINGIFY(x) #x
! #define TOSTRING(x) STRINGIFY(x)

! #DEFINE DASH "/"
! #DEFINE COPY "cp"
! #DEFINE MD "mkdir"
! #DEFINE RMDIR "rm -r"
! #DEFINE RM "rm"
! #DEFINE RENAME "mv"
! #DEFINE SH "sh"
! #DEFINE EXE ""
! #DEFINE NULL ""

! #else

! #define STRINGIFY(x) #x
! #define TOSTRING(x) STRINGIFY(x)

! #DEFINE DASH "\"
! #DEFINE COPY "copy"
! #DEFINE MD "md"
! #DEFINE RMDIR "RMDIR /S /Q"
! #DEFINE RM "del"
! #DEFINE RENAME "MOVE /Y"
! #DEFINE SH "BAT"
! #DEFINE EXE ".exe"
! #DEFINE NULL " >NUL"
! #endif

! module informationModule

! contains

!     subroutine Checker
!         use Global

        
!         use AlphaImputeSpecFileModule
!         use alphahouseMod, only :countLines
!         implicit none

!         integer :: h,i,j,k,l,nAnisTest,CountCatTest(6)
!         integer :: SummaryStats(3,6),Div,CountLen,Counter
!         integer, dimension(:), allocatable :: Work,WorkTmp,GenoStratIndex
!         real :: SummaryProps(3,6),SumPat(6),SumMat(6)
!         character(len=300) :: Names(6),FileName,dumC
!         integer,allocatable,dimension(:) :: RecTestId,FinalSetter
!         integer,allocatable,dimension(:,:) :: TrueGenos,RawGenos,TestMat
!         real,allocatable,dimension(:) :: Correlations
!         real,allocatable,dimension(:,:) :: AnisSummary,RealTestGenos
!         character(len=lengan),allocatable,dimension(:) :: TrueGenosId
!         type(AlphaImputeInput), pointer :: inputParams

!         inputParams => defaultInput

!         allocate(Work(inputParams%nSnpRaw))
!         allocate(WorkTmp(inputParams%nSnpRaw))
!         allocate(GenoStratIndex(ped%pedigreeSize))


!         FileName=trim(inputParams%TrueGenotypeFile)
!         ! call CountLines(FileName,nAnisTest)
!         nAnisTest = CountLines(FileName)

!         call system(RMDIR // " TempTestAlphaImpute")
!         call system(MD // " TempTestAlphaImpute")

!         open (unit=35,file=trim(inputParams%TrueGenotypeFile),status="old")
!         open (unit=36,file=trim(inputParams%GenotypeFile),status="unknown")
!         ! open (unit=37,file="./TempTestAlphaImpute/IndividualAnimalAccuracy.txt",status="unknown")
!         ! open (unit=38,file="./TempTestAlphaImpute/SummaryAnimalAccuracy.txt",status="unknown")
!         ! open (unit=44,file="./TempTestAlphaImpute/IndividualSummaryAccuracy.txt",status="unknown")
!         ! open (unit=45,file="./TempTestAlphaImpute/IndividualSummaryYield.txt",status="unknown")
!         open (unit=37, file="." // DASH // "TempTestAlphaImpute" // DASH // "IndividualAnimalAccuracy.txt", status="unknown")
!         open (unit=38, file="." // DASH // "TempTestAlphaImpute" // DASH // "SummaryAnimalAccuracy.txt", status="unknown")
!         open (unit=44, file="." // DASH // "TempTestAlphaImpute" // DASH // "IndividualSummaryAccuracy.txt", status="unknown")
!         open (unit=45, file="." // DASH // "TempTestAlphaImpute" // DASH // "IndividualSummaryYield.txt", status="unknown")

!         Names(1)="Both Parents Genotyped"
!         Names(2)="Sire and Maternal GrandSire Genotyped"
!         Names(3)="Dam and Paternal Grandsire Genotyped"
!         Names(4)="Sire Genotyped"
!         Names(5)="Dam Genotyped"
!         Names(6)="Other Relatives Genotyped"

!         allocate(FinalSetter(0:ped%pedigreeSize))
!         FinalSetter=0
!         do i=1,ped%nGenotyped
!             read (36,*) dumC,WorkTmp(:)
!             Counter=0
!             do j=1,inputParams%nSnpRaw
!                 if ((WorkTmp(j)>=0).and.(WorkTmp(j)<=2)) Counter=Counter+1
!             enddo
!             if (float(Counter)>(float(inputParams%nSnpRaw)/2)) then
!                 block
!                     integer :: tmpID
!                     tmpID = ped%dictionary%getValue(dumC)
!                     if (tmpId /= DICT_NULL) then
!                         FinalSetter(tmpId)=1
!                     endif
!                 endblock
!             endif
!         enddo
!         rewind(36)

!         if (inputParams%outopt==0) then

!             allocate(TrueGenos(nAnisTest,inputParams%nsnp))
!             allocate(TrueGenosId(nAnisTest))
!             allocate(RawGenos(nAnisTest,inputParams%nsnp))
!             allocate(TestMat(nAnisTest,inputParams%nsnp))
!             allocate(RecTestId(nAnisTest))
!             allocate(AnisSummary(nAnisTest,5))
!             allocate(Correlations(6))
!             allocate(RealTestGenos(nAnisTest,inputParams%nsnp))
!             do i=1,nAnisTest
!                 read (35,*) TrueGenosId(i),Work(:)
!                 k=0
!                 do j=1,inputParams%nSnpRaw
!                     if (SnpIncluded(j)/=0) then
!                         k=k+1
!                         TrueGenos(i,k)=Work(j)
!                     endif
!                 enddo
!             enddo
!             block
!                 integer :: tmpID

!                 do i=1,nAnisTest
!                     tmpID = ped%dictionary%getValue(TrueGenosId(i))
!                     if (tmpID /= dict_null) then
!                         RecTestId(i)=tmpId

!                     endif
!                 enddo
!                 do i=1,ped%nGenotyped
!                     read (36,*) dumC,WorkTmp(:)
!                     do j=1,nAnisTest
!                         ! TODO this check can likely be avoided
!                         if (trim(TrueGenosId(j))==dumC) then
!                             k=0
!                             do l=1,inputParams%nSnpRaw
!                                 if (SnpIncluded(l)==1) then
!                                     k=k+1
!                                     RawGenos(j,k)=WorkTmp(l)
!                                 endif
!                             enddo
!                             exit
!                         endif
!                     enddo
!                 enddo

!             endblock
!             GenoStratIndex(:)=0
!             do i=1,ped%pedigreeSize
!                 if (FinalSetter(i)/=1) then
!                     GenoStratIndex(i)=6
!                     if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
!                         GenoStratIndex(i)=5
!                         if (FinalSetter(ped%Pedigree(i)%getPaternalGrandDamRecodedIndex())==1) then
!                             GenoStratIndex(i)=3
!                         endif
!                     endif
!                     if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
!                         GenoStratIndex(i)=4
!                         if (FinalSetter(ped%Pedigree(i)%getMaternalGrandDamRecodedIndex())==1) then
!                             GenoStratIndex(i)=2
!                         endif
!                     endif
!                     if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1)) then
!                         GenoStratIndex(i)=1
!                     endif
!                 endif
!             enddo
!             TestMat=4
!             CountCatTest=0
!             AnisSummary=0.0
!             do i=1,nAnisTest
!                 do j=1,inputParams%nsnp
!                     if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
!                         TestMat(i,j)=5
!                     else
!                         if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
!                             if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
!                             if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
!                             if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
!                         endif
!                     endif
!                 enddo
!                 write (37,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
!                 CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
!                 Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
!                 AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
!                 AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
!                 AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
!                 AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nsnp)
!                 AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nsnp)
!                 write (44,'(a20,i3,5f7.2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:)
!             enddo
!             SummaryStats=0
!             SumPat=0.0
!             SumMat=0.0
!             do i=1,nAnisTest
!                 do j=1,6
!                     if (GenoStratIndex(RecTestId(i))==j) then
!                         SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
!                         SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
!                         SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
!                         SumPat(j)=SumPat(j)+AnisSummary(i,4)
!                         SumMat(j)=SumMat(j)+AnisSummary(i,5)
!                     endif
!                 enddo
!             enddo

!             SummaryProps=0.0
!             do i=1,3
!                 do j=1,6
!                     if (CountCatTest(j)==0) then
!                         SummaryProps(i,j)=0.0
!                     else
!                         SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
!                         SummaryProps(i,j)=SummaryProps(i,j)*100
!                     endif
!                 enddo
!             enddo

!             do h=1,6
!                 CountLen=0
!                 do i=1,nAnisTest
!                     if(GenoStratIndex(RecTestId(i))==h) then
!                         do j=1,inputParams%nsnp
!                             if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
!                                 if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
!                                     CountLen=CountLen+1
!                                 endif
!                             endif
!                         enddo
!                     endif
!                 enddo
!             enddo

!             do i=1,6
!                 SumPat(i)=SumPat(i)/CountCatTest(i)
!                 SumMat(i)=SumMat(i)/CountCatTest(i)
!                 write (38,'(5f7.2,i7,a40)') SummaryProps(:,i),SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
!             enddo

!             print*, " "
!             do i=1,6
!                 if (CountCatTest(i)>0) write (*,'(3f7.2,a3,a3,2f7.2,i7,a40)') SummaryProps(:,i),"   "&
!                     ,"   ",SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
!             enddo

!             do i=1, ped%pedigreeSize
!                 if (ped%pedigree(i)%isDummy) then
!                     EXIT
!                 endif
!                 write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nsnp &
!                     ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nsnp
!             enddo
!         else
!             allocate(TrueGenos(nAnisTest,inputParams%nSnpRaw))
!             allocate(TrueGenosId(nAnisTest))
!             allocate(RawGenos(nAnisTest,inputParams%nSnpRaw))
!             allocate(TestMat(nAnisTest,inputParams%nSnpRaw))
!             allocate(RecTestId(nAnisTest))
!             allocate(AnisSummary(nAnisTest,5))

!             do i=1,nAnisTest
!                 read (35,*) TrueGenosId(i),TrueGenos(i,:)
!             enddo

!             do i=1,nAnisTest
!                 do j=1,ped%pedigreeSize
!                     if (trim(TrueGenosId(i))==trim((ped%pedigree(j)%originalID))) then
!                         RecTestId(i)=j
!                         exit
!                     endif
!                 enddo
!             enddo

!             do i=1,ped%nGenotyped
!                 read (36,*) dumC,WorkTmp(:)
!                 do j=1,nAnisTest
!                     if (trim(TrueGenosId(j))==dumC) then
!                         RawGenos(j,:)=WorkTmp(:)
!                         exit
!                     endif
!                 enddo
!             enddo
!             GenoStratIndex(:)=0
!             do i=1,ped%pedigreeSize
!                 if (FinalSetter(i)/=1) then
!                     GenoStratIndex(i)=6
!                     if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
!                         GenoStratIndex(i)=5
!                         if (FinalSetter(ped%pedigree(i)%getPaternalGrandSireRecodedIndex())==1) then
!                             GenoStratIndex(i)=3
!                         endif
!                     endif
!                     if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
!                         GenoStratIndex(i)=4
!                         if (FinalSetter(ped%pedigree(i)%getMaternalGrandSireRecodedIndex())==1) then
!                             GenoStratIndex(i)=2
!                         endif
!                     endif
!                     if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1)) then
!                         GenoStratIndex(i)=1
!                     endif
!                 endif
!             enddo
!             TestMat=4
!             CountCatTest=0
!             AnisSummary=0.0
!             do i=1,nAnisTest
!                 do j=1,inputParams%nSnpRaw
!                     if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
!                         TestMat(i,j)=5
!                     else
!                         if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
!                             if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
!                             if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
!                             if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
!                         endif
!                     endif
!                 enddo
!                 write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
!                 CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
!                 Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
!                 AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
!                 AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
!                 AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
!                 AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nSnpRaw)
!                 AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nSnpRaw)
!                 write (44,'(a20,i3,5f7.2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:)
!             enddo

!             SummaryStats=0
!             SumPat=0.0
!             SumMat=0.0
!             do i=1,nAnisTest
!                 do j=1,6
!                     if (GenoStratIndex(RecTestId(i))==j) then
!                         SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
!                         SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
!                         SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
!                         SumPat(j)=SumPat(j)+AnisSummary(i,4)
!                         SumMat(j)=SumMat(j)+AnisSummary(i,5)
!                     endif
!                 enddo
!             enddo

!             SummaryProps=0.0
!             do i=1,3
!                 do j=1,6
!                     if (CountCatTest(j)==0) then
!                         SummaryProps(i,j)=0.0
!                     else
!                         SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
!                         SummaryProps(i,j)=SummaryProps(i,j)*100
!                     endif
!                 enddo
!             enddo

!             do h=1,6
!                 CountLen=0
!                 do i=1,nAnisTest
!                     if(GenoStratIndex(RecTestId(i))==h) then
!                         do j=1,inputParams%nSnpRaw
!                             if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
!                                 if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
!                                     CountLen=CountLen+1
!                                 endif
!                             endif
!                         enddo
!                     endif
!                 enddo
!             enddo

!             do i=1,6
!                 SumPat(i)=SumPat(i)/CountCatTest(i)
!                 SumMat(i)=SumMat(i)/CountCatTest(i)
!                 write (38,'(5f7.2,i7,a40)') SummaryProps(:,i),SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
!             enddo

!             print*, " "
!             do i=1,6
!                 if (CountCatTest(i)>0) write (*,'(3f7.2,a3,a3,2f7.2,i7,a40)') SummaryProps(:,i),"   "&
!                     ,"   ",SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
!             enddo

!             do i=1, ped%pedigreeSize
!                 if (ped%pedigree(i)%isDummy) then
!                     exit
!                 endif
!                 write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nSnpRaw &
!                     ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nSnpRaw
!             enddo

!         endif
!         close(35)
!         close(36)
!         close(37)
!         close(38)
!         close(44)
!         close(45)

!         deallocate(Work)
!         deallocate(WorkTmp)
!         deallocate(GenoStratIndex)


!     end subroutine Checker

!     !#############################################################################################################################################################################################################################

!     subroutine CurrentYield
!         use Global
!         use AlphaImputeSpecFileModule
!         implicit none

!         integer :: CountPatAl,CountMatAl,CountGeno
!         real :: PropPatAl,PropMatAl,PropGeno,NotKnownStart,NotKnownEnd
!         type(AlphaImputeInput), pointer :: inputParams

!         inputParams => defaultInput
!         CountPatAl=count(ImputePhase(:,:,1)==9)
!         CountMatAl=count(ImputePhase(:,:,2)==9)
!         CountGeno=count(ImputeGenos(:,:)/=9)

!         PropPatAl=100*(float(CountPatAl)/(ped%pedigreeSize*inputParams%nsnp))
!         PropMatAl=100*(float(CountMatAl)/(ped%pedigreeSize*inputParams%nsnp))

!         NotKnownStart=(ped%pedigreeSize*inputParams%nsnp)-CountRawGenos
!         NotKnownEnd=(ped%pedigreeSize*inputParams%nsnp)-CountGeno
!         PropGeno=100*((NotKnownStart-NotKnownEnd)/NotKnownStart)

!         print*, " "
!         print*, "           ","Proportion not imputed:"
!         write(*,'(a10,1x,a15,1x,a15,1x,a33)') "           ","Paternal allele","Maternal allele","Proportion missing now genotyped"
!         write (*,'(a10,1x,2f15.2,1x,f33.2)') "           ",PropPatAl,PropMatAl,PropGeno

!     end subroutine CurrentYield



!     !######################################################################################################################################################################################



    ! !#############################################################################################################################################################################################################################

    ! subroutine FinalChecker
    ! ! TODO rewrite
    !     use Global

        
    !     use alphahouseMod, only : countLines
    !     use alphastatmod, only : pearsn, moment
    !     use AlphaImputeSpecFileModule
    !     implicit none

    !     integer :: h,i,j,k,l,nAnisTest,CountCatTest(6),tmpIDInt
    !     integer :: SummaryStats(3,6),Div,CountLen,Counter,Top1,Top2,Top3,Top4,Bot,ContSnpCor,CountValAnim(6)
    !     double precision :: SummaryProps(3,6),SumPat(6),SumMat(6),MeanCorPerInd(6),StdDevPerGrp(6),AveCategoryInformativeness(6,6)
    !     double precision :: Tmpave,Tmpadev,Tmpvar,Tmpskew,Tmpcurt
    !     character(len=300) :: Names(6),FileName,dumC
    !     integer,allocatable,dimension(:) :: RecTestId,FinalSetter,WorkTmp,GenoStratIndex,Work
    !     integer,allocatable,dimension(:,:) :: TrueGenos,RawGenos,TestMat,TestAnimInformativeness
    !     double precision,allocatable,dimension(:) :: Correlations,CorrelationPerAnimal,TmpVarPerGrp
    !     double precision,allocatable,dimension(:,:) :: AnisSummary,WorkVec,RealTestGenos,CalcCorPerAnimal
    !     type(AlphaImputeInput), pointer :: inputParams

    !     character(len=lengan),allocatable,dimension(:) :: TrueGenosId

    !     inputParams =>defaultInput

    !     allocate(Work(inputParams%nSnpRaw))
    !     allocate(WorkTmp(inputParams%nSnpRaw))
    !     allocate(GenoStratIndex(ped%pedigreeSize-ped%nDummys))

    !     FileName=trim(inputParams%TrueGenotypeFile)
    !     ! call CountLines(FileName,nAnisTest)
    !     nAnisTest = CountLines(FileName)


    !     call system(RMDIR // " TestAlphaImpute")
    !     call system(MD // " TestAlphaImpute")

    !     open (unit=35,file=trim(inputParams%TrueGenotypeFile),status="old")
    !     open (unit=36,file=trim(inputParams%GenotypeFile),status="unknown")
    !     ! open (unit=37,file="./TestAlphaImpute/IndividualAnimalAccuracy.txt",status="unknown")
    !     ! open (unit=38,file="./TestAlphaImpute/SummaryAnimalAccuracy.txt",status="unknown")
    !     ! open (unit=44,file="./TestAlphaImpute/IndividualSummaryAccuracy.txt",status="unknown")
    !     ! open (unit=45,file="./TestAlphaImpute/IndividualSummaryYield.txt",status="unknown")
    !     ! open (unit=48,file="./TestAlphaImpute/IndividualSnpAccuracy.txt",status="unknown")
    !     open (unit=37, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualAnimalAccuracy.txt", status="unknown")
    !     open (unit=38, file="." // DASH // "TestAlphaImpute" // DASH // "SummaryAnimalAccuracy.txt", status="unknown")
    !     open (unit=44, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualSummaryAccuracy.txt", status="unknown")
    !     open (unit=45, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualSummaryYield.txt", status="unknown")
    !     open (unit=48, file="." // DASH // "TestAlphaImpute" // DASH // "IndividualSnpAccuracy.txt", status="unknown")

    !     Names(1)="Both Parents Genotyped"
    !     Names(2)="Sire and Maternal GrandSire Genotyped"
    !     Names(3)="Dam and Paternal Grandsire Genotyped"
    !     Names(4)="Sire Genotyped"
    !     Names(5)="Dam Genotyped"
    !     Names(6)="Other Relatives Genotyped"

    !     if (allocated(GlobalTmpCountInf)==.FALSE.) then
    !         allocate(GlobalTmpCountInf(ped%pedigreeSize-ped%nDummys,6))
    !         GlobalTmpCountInf(:,:)=0
    !     endif

    !     allocate(FinalSetter(0:ped%pedigreeSize-ped%nDummys))
    !     FinalSetter=0
    !     do i=1,ped%nGenotyped
    !         read (36,*) dumC,WorkTmp(:)
    !         Counter=0
    !         do j=1,inputParams%nSnpRaw
    !             if ((WorkTmp(j)>=0).and.(WorkTmp(j)<=2)) Counter=Counter+1
    !         enddo
    !         if (float(Counter)>(float(inputParams%nSnpRaw)/2)) then
    !             tmpIDInt = ped%dictionary%getValue(dumC)
    !             if (tmpIDInt /= DICT_NULL) then
    !                 FinalSetter(tmpIDInt)=1
    !             endif

    !         endif
    !     enddo
    !     rewind(36)

    !     allocate(TestAnimInformativeness(nAnisTest,6))

    !     if (inputParams%outopt==0) then

    !         allocate(TrueGenos(nAnisTest,inputParams%nsnp))
    !         allocate(TrueGenosId(nAnisTest))
    !         allocate(RawGenos(nAnisTest,inputParams%nsnp))
    !         allocate(TestMat(nAnisTest,inputParams%nsnp))
    !         allocate(RecTestId(nAnisTest))
    !         allocate(AnisSummary(nAnisTest,5))
    !         allocate(Correlations(6))
    !         allocate(RealTestGenos(nAnisTest,inputParams%nsnp))
    !         allocate(CalcCorPerAnimal(inputParams%nsnp,2))
    !         allocate(CorrelationPerAnimal(nAnisTest))
    !         allocate(TmpVarPerGrp(nAnisTest))

    !         do i=1,nAnisTest
    !             read (35,*) TrueGenosId(i),Work(:)
    !             k=0
    !             do j=1,inputParams%nSnpRaw
    !                 if (SnpIncluded(j)/=0) then
    !                     k=k+1
    !                     TrueGenos(i,k)=Work(j)
    !                 endif
    !             enddo
    !         enddo

    !         RecTestId(:)=-99
    !         do i=1,nAnisTest
    !             do j=1,ped%pedigreeSize-ped%nDummys
    !                 if (trim(TrueGenosId(i))==trim((ped%pedigree(j)%originalID))) then
    !                     RecTestId(i)=j
    !                     TestAnimInformativeness(i,:)=GlobalTmpCountInf(j,1:6)
    !                     exit
    !                 endif
    !             enddo
    !         enddo
    !         if (count(RecTestId(:)==-99)>0) print*, "Error - There seems to be unidentifiablecount ",count(RecTestId(:)==-99)," individuals in the test file"

    !         do i=1,ped%nGenotyped
    !             read (36,*) dumC,WorkTmp(:)
    !             do j=1,nAnisTest
    !                 if (trim(TrueGenosId(j))==dumC) then
    !                     k=0
    !                     do l=1,inputParams%nSnpRaw
    !                         if (SnpIncluded(l)==1) then
    !                             k=k+1
    !                             RawGenos(j,k)=WorkTmp(l)
    !                         endif
    !                     enddo
    !                     exit
    !                 endif
    !             enddo
    !         enddo

    !         GenoStratIndex(:)=0
    !         do i=1,ped%pedigreeSize-ped%nDummys
    !             if (FinalSetter(i)/=1) then
    !                 GenoStratIndex(i)=6
    !                 if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
    !                     GenoStratIndex(i)=5
    !                     if (FinalSetter(ped%pedigree(i)%getPaternalGrandSireRecodedIndex())==1) then
    !                         GenoStratIndex(i)=3
    !                     endif
    !                 endif
    !                 if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
    !                     GenoStratIndex(i)=4
    !                     if (FinalSetter(ped%pedigree(i)%getMaternalGrandSireRecodedIndex())==1) then
    !                         GenoStratIndex(i)=2
    !                     endif
    !                 endif
    !                 if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1)) then
    !                     GenoStratIndex(i)=1
    !                 endif
    !             endif
    !         enddo

    !         TestMat=4
    !         CountCatTest=0
    !         AnisSummary=0.0

    !         do i=1,nAnisTest
    !             do j=1,inputParams%nsnp
    !                 if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
    !                     TestMat(i,j)=5
    !                 else
    !                     if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
    !                         if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
    !                         if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
    !                         if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
    !                     endif
    !                 endif
    !             enddo
    !             RealTestGenos(i,:)=ProbImputeGenos(RecTestId(i),:)
    !             write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
    !             CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
    !             Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
    !             AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
    !             AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
    !             AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
    !             AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nsnp)
    !             AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nsnp)
    !         enddo

    !         SummaryStats=0
    !         SumPat=0.0
    !         SumMat=0.0
    !         do i=1,nAnisTest
    !             do j=1,6
    !                 if (GenoStratIndex(RecTestId(i))==j) then
    !                     SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
    !                     SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
    !                     SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
    !                     SumPat(j)=SumPat(j)+AnisSummary(i,4)
    !                     SumMat(j)=SumMat(j)+AnisSummary(i,5)
    !                 endif
    !             enddo
    !         enddo

    !         SummaryProps=0.0
    !         do i=1,3
    !             do j=1,6
    !                 if (CountCatTest(j)==0) then
    !                     SummaryProps(i,j)=0.0
    !                 else
    !                     SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
    !                     SummaryProps(i,j)=SummaryProps(i,j)*100
    !                 endif
    !             enddo
    !         enddo

    !         MeanCorPerInd(:)=0.0
    !         CountValAnim(:)=0
    !         AveCategoryInformativeness(:,:)=0.0
    !         do h=1,6
    !             CountLen=0
    !             do i=1,nAnisTest
    !                 if(GenoStratIndex(RecTestId(i))==h) then
    !                     AveCategoryInformativeness(h,:)=AveCategoryInformativeness(h,:)+TestAnimInformativeness(i,:)
    !                     do j=1,inputParams%nsnp
    !                         if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
    !                             if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
    !                                 CountLen=CountLen+1
    !                             endif
    !                         endif
    !                     enddo
    !                 endif
    !             enddo
    !             allocate(WorkVec(CountLen,2))
    !             CountLen=0
    !             do i=1,nAnisTest
    !                 if(GenoStratIndex(RecTestId(i))==h) then
    !                     ContSnpCor=0
    !                     do j=1,inputParams%nsnp
    !                         if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
    !                             if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
    !                                 ContSnpCor=ContSnpCor+1                                         !#HereToday
    !                                 CountLen=CountLen+1
    !                                 WorkVec(CountLen,1)=float(TrueGenos(i,j))-(2*maf(j))
    !                                 WorkVec(CountLen,2)=RealTestGenos(i,j)-(2*maf(j))
    !                                 CalcCorPerAnimal(ContSnpCor,1)=float(TrueGenos(i,j))-(2*maf(j))
    !                                 CalcCorPerAnimal(ContSnpCor,2)=RealTestGenos(i,j)-(2*maf(j))
    !                             endif
    !                         endif
    !                     enddo
    !                     if (ContSnpCor>5) then
    !                         call Pearsn (CalcCorPerAnimal(1:ContSnpCor,1),CalcCorPerAnimal(1:ContSnpCor,2),ContSnpCor,CorrelationPerAnimal(i))
    !                         MeanCorPerInd(h)=MeanCorPerInd(h)+CorrelationPerAnimal(i)
    !                         CountValAnim(h)=CountValAnim(h)+1
    !                         TmpVarPerGrp(CountValAnim(h))=CorrelationPerAnimal(i)

    !                     else
    !                         CorrelationPerAnimal(i)=-99.0
    !                     endif

    !                 endif
    !             enddo
    !             if (CountLen>5) then
    !                 call Pearsn (WorkVec(:,1),WorkVec(:,2),CountLen,Correlations(h))
    !             else
    !                 Correlations(h)=-99.0
    !             endif
    !             deallocate(WorkVec)
    !             if (CountValAnim(h)>5) then
    !                 MeanCorPerInd(h)=MeanCorPerInd(h)/CountValAnim(h)

    !                 call moment(TmpVarPerGrp(1:CountValAnim(h)),CountValAnim(h),Tmpave,Tmpadev,StdDevPerGrp(h),Tmpvar,Tmpskew,Tmpcurt)
    !             else
    !                 MeanCorPerInd(h)=-99.0
    !             endif
    !         enddo

    !         do i=1,6
    !             SumPat(i)=SumPat(i)/CountCatTest(i)
    !             SumMat(i)=SumMat(i)/CountCatTest(i)
    !             AveCategoryInformativeness(i,:)=AveCategoryInformativeness(i,:)/CountCatTest(i)
    !             write (38,'(14f7.2,i7,a40)') SummaryProps(:,i),Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i),SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    !         enddo

    !         print*, " "
    !         do i=1,6
    !             if (CountCatTest(i)>0) write (*,'(3f7.2,a3,3f7.2,a3,8f7.2,i7,a40)') SummaryProps(:,i),"   ",Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i)&
    !                 ,"   ",SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    !         enddo

    !         do i=1, ped%pedigreeSize-ped%nDummys
    !             if (ped%pedigree(i)%isDummy) then
    !                 exit
    !             endif
    !             write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nsnp &
    !                 ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nsnp
    !         enddo

    !         do j=1,inputParams%nsnp
    !             Bot=count(TestMat(:,j)/=4)
    !             Top1=count(TestMat(:,j)==1)
    !             Top2=count(TestMat(:,j)==2)
    !             Top3=count(TestMat(:,j)==3)
    !             Top4=count(TestMat(:,j)==5)
    !             write (48,'(2i10,4f9.2)') j,Bot,100*(float(Top1)/Bot),100*(float(Top2)/Bot),100*(float(Top3)/Bot),100*(float(Top4)/Bot)
    !         enddo

    !         do i=1,nAnisTest
    !             write (44,'(a20,i3,6f7.2,6i10)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:),CorrelationPerAnimal(i),TestAnimInformativeness(i,:)
    !         enddo


    !     else
    !         allocate(TrueGenos(nAnisTest,inputParams%nSnpRaw))
    !         allocate(TrueGenosId(nAnisTest))
    !         allocate(RawGenos(nAnisTest,inputParams%nSnpRaw))
    !         allocate(TestMat(nAnisTest,inputParams%nSnpRaw))
    !         allocate(RecTestId(nAnisTest))
    !         allocate(AnisSummary(nAnisTest,5))
    !         allocate(Correlations(6))
    !         allocate(RealTestGenos(nAnisTest,inputParams%nSnpRaw))
    !         allocate(CalcCorPerAnimal(inputParams%nSnpRaw,2))
    !         allocate(CorrelationPerAnimal(nAnisTest))
    !         allocate(TmpVarPerGrp(nAnisTest))


    !         RecTestId(:)=-99
    !         do i=1,nAnisTest
    !             read (35,*) TrueGenosId(i),TrueGenos(i,:)
    !             tmpIDInt = ped%dictionary%getValue(TrueGenosId(i))

    !             if (tmpIDInt /= DICT_NULL) then
    !                 RecTestId(i)=tmpIDInt
    !                 testAnimInformativeness(i,:)=GlobalTmpCountInf(tmpIDInt,1:6)
    !             endif
    !         enddo

    !         if (any(RecTestId(:)==-99)) print*, "Error - There seems to be unidentifiablecount ",count(RecTestId(:)==-99)," individuals in the test file"

    !         do i=1,ped%nGenotyped
    !             read (36,*) dumC,WorkTmp(:)
    !             do j=1,nAnisTest
    !                 if (trim(TrueGenosId(j))==dumC) then
    !                     RawGenos(j,:)=WorkTmp(:)
    !                     exit
    !                 endif
    !             enddo
    !         enddo
    !         GenoStratIndex(:)=0
    !         do i=1,ped%pedigreeSize-ped%nDummys
    !             if (FinalSetter(i)/=1) then
    !                 GenoStratIndex(i)=6
    !                 if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
    !                     GenoStratIndex(i)=5
    !                     if (FinalSetter(ped%pedigree(i)%getPaternalGrandSireRecodedIndex())==1) then
    !                         GenoStratIndex(i)=3
    !                     endif
    !                 endif
    !                 if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
    !                     GenoStratIndex(i)=4
    !                     if (FinalSetter(ped%pedigree(i)%getMaternalGrandSireRecodedIndex())==1) then
    !                         GenoStratIndex(i)=2
    !                     endif
    !                 endif
    !                 if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1)) then
    !                     GenoStratIndex(i)=1
    !                 endif
    !             endif
    !         enddo
    !         TestMat=4
    !         CountCatTest=0
    !         AnisSummary=0.0
    !         do i=1,nAnisTest
    !             do j=1,inputParams%nSnpRaw
    !                 if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
    !                     TestMat(i,j)=5
    !                 else
    !                     if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
    !                         if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
    !                         if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
    !                         if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
    !                     endif
    !                 endif
    !             enddo
    !             RealTestGenos(i,:)=ProbImputeGenos(RecTestId(i),:)
    !             write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
    !             CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
    !             Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
    !             AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
    !             AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
    !             AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
    !             AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nSnpRaw)
    !             AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nSnpRaw)
    !         enddo

    !         SummaryStats=0
    !         SumPat=0.0
    !         SumMat=0.0
    !         do i=1,nAnisTest
    !             do j=1,6
    !                 if (GenoStratIndex(RecTestId(i))==j) then
    !                     SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
    !                     SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
    !                     SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
    !                     SumPat(j)=SumPat(j)+AnisSummary(i,4)
    !                     SumMat(j)=SumMat(j)+AnisSummary(i,5)
    !                 endif
    !             enddo
    !         enddo

    !         SummaryProps=0.0
    !         do i=1,3
    !             do j=1,6
    !                 if (CountCatTest(j)==0) then
    !                     SummaryProps(i,j)=0.0
    !                 else
    !                     SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
    !                     SummaryProps(i,j)=SummaryProps(i,j)*100
    !                 endif
    !             enddo
    !         enddo

    !         MeanCorPerInd(:)=0.0
    !         CountValAnim(:)=0
    !         AveCategoryInformativeness(:,:)=0.0
    !         do h=1,6
    !             CountLen=0
    !             do i=1,nAnisTest
    !                 if(GenoStratIndex(RecTestId(i))==h) then
    !                     AveCategoryInformativeness(h,:)=AveCategoryInformativeness(h,:)+TestAnimInformativeness(i,:)
    !                     do j=1,inputParams%nSnpRaw
    !                         if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
    !                             if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
    !                                 CountLen=CountLen+1
    !                             endif
    !                         endif
    !                     enddo
    !                 endif
    !             enddo
    !             allocate(WorkVec(CountLen,2))
    !             CountLen=0
    !             do i=1,nAnisTest
    !                 if(GenoStratIndex(RecTestId(i))==h) then
    !                     ContSnpCor=0
    !                     do j=1,inputParams%nSnpRaw
    !                         if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
    !                             if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
    !                                 ContSnpCor=ContSnpCor+1
    !                                 CountLen=CountLen+1
    !                                 WorkVec(CountLen,1)=float(TrueGenos(i,j))-(2*maf(j))
    !                                 WorkVec(CountLen,2)=RealTestGenos(i,j)-(2*maf(j))
    !                                 CalcCorPerAnimal(ContSnpCor,1)=float(TrueGenos(i,j))-(2*maf(j))
    !                                 CalcCorPerAnimal(ContSnpCor,2)=RealTestGenos(i,j)-(2*maf(j))
    !                             endif
    !                         endif
    !                     enddo
    !                     if (ContSnpCor>5) then
    !                         call Pearsn (CalcCorPerAnimal(1:ContSnpCor,1),CalcCorPerAnimal(1:ContSnpCor,2),ContSnpCor,CorrelationPerAnimal(i))
    !                         MeanCorPerInd(h)=MeanCorPerInd(h)+CorrelationPerAnimal(i)
    !                         CountValAnim(h)=CountValAnim(h)+1
    !                         TmpVarPerGrp(CountValAnim(h))=CorrelationPerAnimal(i)
    !                     else
    !                         CorrelationPerAnimal(i)=-99.0
    !                     endif
    !                 endif
    !             enddo

    !             if (CountLen>5) then
    !                 call Pearsn(WorkVec(:,1),WorkVec(:,2),CountLen,Correlations(h))
    !             else
    !                 Correlations(h)=-99.0
    !             endif

    !             deallocate(WorkVec)

    !             if (CountValAnim(h)>5) then
    !                 MeanCorPerInd(h)=MeanCorPerInd(h)/CountValAnim(h)
    !                 call moment(TmpVarPerGrp(1:CountValAnim(h)),CountValAnim(h),Tmpave,Tmpadev,StdDevPerGrp(h),Tmpvar,Tmpskew,Tmpcurt)
    !             else
    !                 MeanCorPerInd(h)=-99.0
    !             endif
    !         enddo

    !         do i=1,6
    !             SumPat(i)=SumPat(i)/CountCatTest(i)
    !             SumMat(i)=SumMat(i)/CountCatTest(i)
    !             AveCategoryInformativeness(i,:)=AveCategoryInformativeness(i,:)/CountCatTest(i)
    !             write (38,'(14f7.2,i7,a40)') SummaryProps(:,i),Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i),SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    !         enddo

    !         print*, " "
    !         do i=1,6
    !             if (CountCatTest(i)>0) write (*,'(3f7.2,a3,3f7.2,a3,8f7.2,i7,a40)') SummaryProps(:,i),"   ",Correlations(i),MeanCorPerInd(i),StdDevPerGrp(i)&
    !                 ,"   ",SumPat(i),SumMat(i),AveCategoryInformativeness(i,:),CountCatTest(i),trim(Names(i))
    !         enddo

    !         do i=1, ped%pedigreeSize-ped%nDummys
    !             if (ped%pedigree(i)%isDummy) then
    !                 exit
    !             endif
    !             write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nSnpRaw &
    !                 ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nSnpRaw
    !         enddo

    !         do j=1,inputParams%nsnp
    !             Bot=count(TestMat(:,j)/=4)
    !             Top1=count(TestMat(:,j)==1)
    !             Top2=count(TestMat(:,j)==2)
    !             Top3=count(TestMat(:,j)==3)
    !             Top4=count(TestMat(:,j)==5)
    !             write (48,'(2i10,4f9.2)') j,Bot,100*(float(Top1)/Bot),100*(float(Top2)/Bot),100*(float(Top3)/Bot),100*(float(Top4)/Bot)
    !         enddo

    !         do i=1,nAnisTest
    !             write (44,'(a20,i3,6f7.2,6i10)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:),CorrelationPerAnimal(i),TestAnimInformativeness(i,:)
    !         enddo
    !     endif

    !     deallocate(Work)
    !     deallocate(WorkTmp)
    !     deallocate(GenoStratIndex)
    ! end subroutine FinalChecker


! end module informationModule
