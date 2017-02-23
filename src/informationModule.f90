#ifdef OS_UNIX

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"
#DEFINE SH "sh"
#DEFINE EXE ""
#DEFINE NULL ""

#else

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S /Q"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#DEFINE SH "BAT"
#DEFINE EXE ".exe"
#DEFINE NULL " >NUL"
#endif

module informationModule

contains

    subroutine Checker
        use Global

        use Utils
        use alphaimputeinmod
        use alphahouseMod, only :countLines
        implicit none

        integer :: h,i,j,k,l,nAnisTest,CountCatTest(6)
        integer :: SummaryStats(3,6),Div,CountLen,Counter
        integer, dimension(:), allocatable :: Work,WorkTmp,GenoStratIndex
        real :: SummaryProps(3,6),SumPat(6),SumMat(6)
        character(len=300) :: Names(6),FileName,dumC
        integer,allocatable,dimension(:) :: RecTestId,FinalSetter
        integer,allocatable,dimension(:,:) :: TrueGenos,RawGenos,TestMat
        real,allocatable,dimension(:) :: Correlations
        real,allocatable,dimension(:,:) :: AnisSummary,RealTestGenos
        character(len=lengan),allocatable,dimension(:) :: TrueGenosId
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        allocate(Work(inputParams%nSnpRaw))
        allocate(WorkTmp(inputParams%nSnpRaw))
        allocate(GenoStratIndex(ped%pedigreeSize))


        FileName=trim(inputParams%TrueGenotypeFile)
        ! call CountLines(FileName,nAnisTest)
        nAnisTest = CountLines(FileName)

        call system(RMDIR // " TempTestAlphaImpute")
        call system(MD // " TempTestAlphaImpute")

        open (unit=35,file=trim(inputParams%TrueGenotypeFile),status="old")
        open (unit=36,file=trim(inputParams%GenotypeFile),status="unknown")
        ! open (unit=37,file="./TempTestAlphaImpute/IndividualAnimalAccuracy.txt",status="unknown")
        ! open (unit=38,file="./TempTestAlphaImpute/SummaryAnimalAccuracy.txt",status="unknown")
        ! open (unit=44,file="./TempTestAlphaImpute/IndividualSummaryAccuracy.txt",status="unknown")
        ! open (unit=45,file="./TempTestAlphaImpute/IndividualSummaryYield.txt",status="unknown")
        open (unit=37, file="." // DASH // "TempTestAlphaImpute" // DASH // "IndividualAnimalAccuracy.txt", status="unknown")
        open (unit=38, file="." // DASH // "TempTestAlphaImpute" // DASH // "SummaryAnimalAccuracy.txt", status="unknown")
        open (unit=44, file="." // DASH // "TempTestAlphaImpute" // DASH // "IndividualSummaryAccuracy.txt", status="unknown")
        open (unit=45, file="." // DASH // "TempTestAlphaImpute" // DASH // "IndividualSummaryYield.txt", status="unknown")

        Names(1)="Both Parents Genotyped"
        Names(2)="Sire and Maternal GrandSire Genotyped"
        Names(3)="Dam and Paternal Grandsire Genotyped"
        Names(4)="Sire Genotyped"
        Names(5)="Dam Genotyped"
        Names(6)="Other Relatives Genotyped"

        allocate(FinalSetter(0:ped%pedigreeSize))
        FinalSetter=0
        do i=1,ped%nGenotyped
            read (36,*) dumC,WorkTmp(:)
            Counter=0
            do j=1,inputParams%nSnpRaw
                if ((WorkTmp(j)>=0).and.(WorkTmp(j)<=2)) Counter=Counter+1
            enddo
            if (float(Counter)>(float(inputParams%nSnpRaw)/2)) then
                block
                    integer :: tmpID
                    tmpID = ped%dictionary%getValue(dumC)
                    if (tmpId /= DICT_NULL) then
                        FinalSetter(tmpId)=1
                    endif
                endblock
            endif
        enddo
        rewind(36)

        if (inputParams%outopt==0) then

            allocate(TrueGenos(nAnisTest,inputParams%nsnp))
            allocate(TrueGenosId(nAnisTest))
            allocate(RawGenos(nAnisTest,inputParams%nsnp))
            allocate(TestMat(nAnisTest,inputParams%nsnp))
            allocate(RecTestId(nAnisTest))
            allocate(AnisSummary(nAnisTest,5))
            allocate(Correlations(6))
            allocate(RealTestGenos(nAnisTest,inputParams%nsnp))
            do i=1,nAnisTest
                read (35,*) TrueGenosId(i),Work(:)
                k=0
                do j=1,inputParams%nSnpRaw
                    if (SnpIncluded(j)/=0) then
                        k=k+1
                        TrueGenos(i,k)=Work(j)
                    endif
                enddo
            enddo
            block
                integer :: tmpID

                do i=1,nAnisTest
                    tmpID = ped%dictionary%getValue(TrueGenosId(i))
                    if (tmpID /= dict_null) then
                        RecTestId(i)=tmpId

                    endif
                enddo
                do i=1,ped%nGenotyped
                    read (36,*) dumC,WorkTmp(:)
                    do j=1,nAnisTest
                        ! TODO this check can likely be avoided
                        if (trim(TrueGenosId(j))==dumC) then
                            k=0
                            do l=1,inputParams%nSnpRaw
                                if (SnpIncluded(l)==1) then
                                    k=k+1
                                    RawGenos(j,k)=WorkTmp(l)
                                endif
                            enddo
                            exit
                        endif
                    enddo
                enddo

            endblock
            GenoStratIndex(:)=0
            do i=1,ped%pedigreeSize
                if (FinalSetter(i)/=1) then
                    GenoStratIndex(i)=6
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
                        GenoStratIndex(i)=5
                        if (FinalSetter(ped%Pedigree(i)%getPaternalGrandDamRecodedIndex())==1) then
                            GenoStratIndex(i)=3
                        endif
                    endif
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
                        GenoStratIndex(i)=4
                        if (FinalSetter(ped%Pedigree(i)%getMaternalGrandDamRecodedIndex())==1) then
                            GenoStratIndex(i)=2
                        endif
                    endif
                    if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1)) then
                        GenoStratIndex(i)=1
                    endif
                endif
            enddo
            TestMat=4
            CountCatTest=0
            AnisSummary=0.0
            do i=1,nAnisTest
                do j=1,inputParams%nsnp
                    if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
                        TestMat(i,j)=5
                    else
                        if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                            if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
                            if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
                            if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
                        endif
                    endif
                enddo
                write (37,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
                CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
                Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
                AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
                AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
                AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
                AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nsnp)
                AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nsnp)
                write (44,'(a20,i3,5f7.2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:)
            enddo
            SummaryStats=0
            SumPat=0.0
            SumMat=0.0
            do i=1,nAnisTest
                do j=1,6
                    if (GenoStratIndex(RecTestId(i))==j) then
                        SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
                        SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
                        SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
                        SumPat(j)=SumPat(j)+AnisSummary(i,4)
                        SumMat(j)=SumMat(j)+AnisSummary(i,5)
                    endif
                enddo
            enddo

            SummaryProps=0.0
            do i=1,3
                do j=1,6
                    if (CountCatTest(j)==0) then
                        SummaryProps(i,j)=0.0
                    else
                        SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
                        SummaryProps(i,j)=SummaryProps(i,j)*100
                    endif
                enddo
            enddo

            do h=1,6
                CountLen=0
                do i=1,nAnisTest
                    if(GenoStratIndex(RecTestId(i))==h) then
                        do j=1,inputParams%nsnp
                            if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                                if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                                    CountLen=CountLen+1
                                endif
                            endif
                        enddo
                    endif
                enddo
            enddo

            do i=1,6
                SumPat(i)=SumPat(i)/CountCatTest(i)
                SumMat(i)=SumMat(i)/CountCatTest(i)
                write (38,'(5f7.2,i7,a40)') SummaryProps(:,i),SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
            enddo

            print*, " "
            do i=1,6
                if (CountCatTest(i)>0) write (*,'(3f7.2,a3,a3,2f7.2,i7,a40)') SummaryProps(:,i),"   "&
                    ,"   ",SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
            enddo

            do i=1, ped%pedigreeSize
                if (ped%pedigree(i)%isDummy) then
                    EXIT
                endif
                write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nsnp &
                    ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nsnp
            enddo
        else
            allocate(TrueGenos(nAnisTest,inputParams%nSnpRaw))
            allocate(TrueGenosId(nAnisTest))
            allocate(RawGenos(nAnisTest,inputParams%nSnpRaw))
            allocate(TestMat(nAnisTest,inputParams%nSnpRaw))
            allocate(RecTestId(nAnisTest))
            allocate(AnisSummary(nAnisTest,5))

            do i=1,nAnisTest
                read (35,*) TrueGenosId(i),TrueGenos(i,:)
            enddo

            do i=1,nAnisTest
                do j=1,ped%pedigreeSize
                    if (trim(TrueGenosId(i))==trim((ped%pedigree(j)%originalID))) then
                        RecTestId(i)=j
                        exit
                    endif
                enddo
            enddo

            do i=1,ped%nGenotyped
                read (36,*) dumC,WorkTmp(:)
                do j=1,nAnisTest
                    if (trim(TrueGenosId(j))==dumC) then
                        RawGenos(j,:)=WorkTmp(:)
                        exit
                    endif
                enddo
            enddo
            GenoStratIndex(:)=0
            do i=1,ped%pedigreeSize
                if (FinalSetter(i)/=1) then
                    GenoStratIndex(i)=6
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1) then
                        GenoStratIndex(i)=5
                        if (FinalSetter(ped%pedigree(i)%getPaternalGrandSireRecodedIndex())==1) then
                            GenoStratIndex(i)=3
                        endif
                    endif
                    if (FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1) then
                        GenoStratIndex(i)=4
                        if (FinalSetter(ped%pedigree(i)%getMaternalGrandSireRecodedIndex())==1) then
                            GenoStratIndex(i)=2
                        endif
                    endif
                    if ((FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(2))==1).and.(FinalSetter(ped%pedigree(i)%getSireDamNewIDByIndex(3))==1)) then
                        GenoStratIndex(i)=1
                    endif
                endif
            enddo
            TestMat=4
            CountCatTest=0
            AnisSummary=0.0
            do i=1,nAnisTest
                do j=1,inputParams%nSnpRaw
                    if ((TrueGenos(i,j)<0).or.(TrueGenos(i,j)>2)) then
                        TestMat(i,j)=5
                    else
                        if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                            if (TrueGenos(i,j)==ImputeGenos(RecTestId(i),j)) TestMat(i,j)=1
                            if (TrueGenos(i,j)/=ImputeGenos(RecTestId(i),j)) TestMat(i,j)=2
                            if (ImputeGenos(RecTestId(i),j)==9) TestMat(i,j)=3
                        endif
                    endif
                enddo
                write (37,'(a20,i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),TestMat(i,:)
                CountCatTest(GenoStratIndex(RecTestId(i)))=CountCatTest(GenoStratIndex(RecTestId(i)))+1
                Div=(count(TestMat(i,:)==1))+(count(TestMat(i,:)==2))+(count(TestMat(i,:)==3))
                AnisSummary(i,1)=100*(float(count(TestMat(i,:)==1))/Div)
                AnisSummary(i,2)=100*(float(count(TestMat(i,:)==2))/Div)
                AnisSummary(i,3)=100*(float(count(TestMat(i,:)==3))/Div)
                AnisSummary(i,4)=100*(float(count(ImputePhase(RecTestId(i),:,1)/=9))/inputParams%nSnpRaw)
                AnisSummary(i,5)=100*(float(count(ImputePhase(RecTestId(i),:,2)/=9))/inputParams%nSnpRaw)
                write (44,'(a20,i3,5f7.2)') TrueGenosId(i),GenoStratIndex(RecTestId(i)),AnisSummary(i,:)
            enddo

            SummaryStats=0
            SumPat=0.0
            SumMat=0.0
            do i=1,nAnisTest
                do j=1,6
                    if (GenoStratIndex(RecTestId(i))==j) then
                        SummaryStats(1,j)=SummaryStats(1,j)+(count(TestMat(i,:)==1))
                        SummaryStats(2,j)=SummaryStats(2,j)+(count(TestMat(i,:)==2))
                        SummaryStats(3,j)=SummaryStats(3,j)+(count(TestMat(i,:)==3))
                        SumPat(j)=SumPat(j)+AnisSummary(i,4)
                        SumMat(j)=SumMat(j)+AnisSummary(i,5)
                    endif
                enddo
            enddo

            SummaryProps=0.0
            do i=1,3
                do j=1,6
                    if (CountCatTest(j)==0) then
                        SummaryProps(i,j)=0.0
                    else
                        SummaryProps(i,j)=float(SummaryStats(i,j))/sum(SummaryStats(:,j))
                        SummaryProps(i,j)=SummaryProps(i,j)*100
                    endif
                enddo
            enddo

            do h=1,6
                CountLen=0
                do i=1,nAnisTest
                    if(GenoStratIndex(RecTestId(i))==h) then
                        do j=1,inputParams%nSnpRaw
                            if ((RawGenos(i,j)<0).or.(RawGenos(i,j)>2)) then
                                if ((TrueGenos(i,j)>=0).and.(TrueGenos(i,j)<=2)) then
                                    CountLen=CountLen+1
                                endif
                            endif
                        enddo
                    endif
                enddo
            enddo

            do i=1,6
                SumPat(i)=SumPat(i)/CountCatTest(i)
                SumMat(i)=SumMat(i)/CountCatTest(i)
                write (38,'(5f7.2,i7,a40)') SummaryProps(:,i),SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
            enddo

            print*, " "
            do i=1,6
                if (CountCatTest(i)>0) write (*,'(3f7.2,a3,a3,2f7.2,i7,a40)') SummaryProps(:,i),"   "&
                    ,"   ",SumPat(i),SumMat(i),CountCatTest(i),trim(Names(i))
            enddo

            do i=1, ped%pedigreeSize
                if (ped%pedigree(i)%isDummy) then
                    exit
                endif
                write (45,'(a25,i3,2f7.2)') ped%pedigree(i)%originalID,FinalSetter(i),float(count(ImputePhase(i,:,1)/=9))/inputParams%nSnpRaw &
                    ,float(count(ImputePhase(i,:,2)/=9))/inputParams%nSnpRaw
            enddo

        endif
        close(35)
        close(36)
        close(37)
        close(38)
        close(44)
        close(45)

        deallocate(Work)
        deallocate(WorkTmp)
        deallocate(GenoStratIndex)


    end subroutine Checker

    !#############################################################################################################################################################################################################################

    subroutine CurrentYield
        use Global
        use alphaimputeinmod
        implicit none

        integer :: CountPatAl,CountMatAl,CountGeno
        real :: PropPatAl,PropMatAl,PropGeno,NotKnownStart,NotKnownEnd
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput
        CountPatAl=count(ImputePhase(:,:,1)==9)
        CountMatAl=count(ImputePhase(:,:,2)==9)
        CountGeno=count(ImputeGenos(:,:)/=9)

        PropPatAl=100*(float(CountPatAl)/(ped%pedigreeSize*inputParams%nsnp))
        PropMatAl=100*(float(CountMatAl)/(ped%pedigreeSize*inputParams%nsnp))

        NotKnownStart=(ped%pedigreeSize*inputParams%nsnp)-CountRawGenos
        NotKnownEnd=(ped%pedigreeSize*inputParams%nsnp)-CountGeno
        PropGeno=100*((NotKnownStart-NotKnownEnd)/NotKnownStart)

        print*, " "
        print*, "           ","Proportion not imputed:"
        write(*,'(a10,1x,a15,1x,a15,1x,a33)') "           ","Paternal allele","Maternal allele","Proportion missing now genotyped"
        write (*,'(a10,1x,2f15.2,1x,f33.2)') "           ",PropPatAl,PropMatAl,PropGeno

    end subroutine CurrentYield


    subroutine FromHMM2ImputePhase
        ! Impute alleles from HMM dosage probabilities
        use Global
        use GlobalVariablesHmmMaCH

        use AlphaImputeInMod

        implicit none

        integer :: i,j,k
        type(AlphaImputeInput), pointer :: inputParams

        inputParams => defaultInput

        do i=1,ped%nGenotyped
            do j=1,inputParams%nsnp
                do k=1,2
                    if (FullH(i,j,k)<0.001.and.FullH(i,j,k)>=0.0) Then
                        ImputePhase(i,j,k)=0
                    elseif (FullH(i,j,k)>0.999.and.FullH(i,j,k)<=1.0) then
                        ImputePhase(i,j,k)=1
                    else
                        ImputePhase(i,j,k)=9
                    endif
                enddo
            enddo
        enddo

    end subroutine FromHMM2ImputePhase

    !######################################################################################################################################################################################

    subroutine InsteadOfReReadGeneProb
        ! Phase alleles in the SEX CHROMOSOME whenever it is possible (homozygous case).
        ! Phasing information is store in the variable GlobalWorkPhase
        use Global
        use AlphaImputeInMod
        implicit none

        type(AlphaImputeInput), pointer :: inputParams
        integer :: e,i,j,ParId
        integer, dimension(:,:) , allocatable :: Genos !  Temp variable


        genos = ped%getGenotypesAsArray()

        inputParams => defaultInput
        if (defaultInput%SexOpt==1) then                                         ! Sex chromosome
            deallocate(GlobalWorkPhase)
            allocate(GlobalWorkPhase(0:ped%pedigreeSize-ped%nDummys,inputParams%nsnp,2))
            GlobalWorkPhase=9
            do i=1,ped%pedigreeSize
                do j=1,inputParams%nsnp                                                     ! Phase alleles in the homozygous case
                    if (Genos(i,j)==0) GlobalWorkPhase(i,j,:)=0
                    if (Genos(i,j)==2) GlobalWorkPhase(i,j,:)=1
                enddo
                if (ped%pedigree(i)%gender/=inputParams%hetGameticStatus) then
                    do e=1,2                                                    ! Phase alleles for homogametic individuals in the homozygous case
                        ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
                        do j=1,inputParams%nsnp
                            if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,e)=0
                            if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,e)=1
                        enddo
                    enddo
                else
                    ParId=ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1)                          ! Phase alleles for heterogametic individuals in the homozygous case
                    do j=1,inputParams%nsnp
                        if (Genos(ParId,j)==0) GlobalWorkPhase(i,j,:)=0
                        if (Genos(ParId,j)==2) GlobalWorkPhase(i,j,:)=1
                    enddo
                endif
            enddo
            GlobalWorkPhase(0,:,:)=9
        else                                ! Nothing is done in other chromosomes
            !! WARNING: This should be some copied, pasted and erased stuff
        endif

    end subroutine InsteadOfReReadGeneProb


end module informationModule
