MODULE Recombination

CONTAINS

!#############################################################################################################################################################################################################################
subroutine ModelRecombSegment(nSnpBegin,nSnpFinal,nRecomb)
use Global

use AlphaImputeInMod
implicit none

INTEGER, INTENT(IN) :: nSnpBegin,nSnpFinal
INTEGER, INTENT(OUT) :: nRecomb(:)

! Local Variables
integer :: e,i,j,k,l,SuperJ,StartDisFound,EndDisFound,HetEnd,HetStart,RSide,LSide,PatMat,SireDamRL,dimSnps,nRec
integer :: StartDisPrev,EndDisPrev
integer :: GamA,GamB,StartDisOld,StartDisTmp
integer :: CountRightSwitch,CountLeftSwitch,PedId,StartDis,EndDis
integer(kind=1),allocatable,dimension(:,:,:) :: WorkPhase,TempWork
integer,allocatable,dimension(:) :: WorkLeft,WorkRight,TempVec,StR,EnR,StRNarrow,EnRNarrow
real,allocatable,dimension(:) :: LengthVec
real,allocatable,dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
character(len=7) :: cm
type(AlphaImputeInput), pointer :: inputParams

inputParams => defaultInput

write(cm,'(I7)') inputParams%nSnpRaw !for formatting
cm = adjustl(cm)

dimSnps=nSnpFinal-nSnpBegin+1

! allocate(GlobalWorkPhase(0:nAnisP,nSnp,2))
allocate(WorkPhase(0:nAnisP,dimSnps,2))
allocate(TempVec(dimSnps))
allocate(LengthVec(dimSnps))
allocate(WorkLeft(dimSnps))
allocate(WorkRight(dimSnps))
allocate(TempWork(0:nAnisP,dimSnps,2))
allocate(PatAlleleProb(dimSnps,2))
allocate(MatAlleleProb(dimSnps,2))
allocate(GeneProbWork(dimSnps,4))
allocate(StR(dimSnps))
allocate(EnR(dimSnps))
allocate(StRNarrow(dimSnps))
allocate(EnRNarrow(dimSnps))

WorkPhase=9

l=0
do j=nSnpBegin,nSnpFinal

    if (SnpIncluded(j)==1) then
        l=l+1
        WorkPhase(:,j,1)=GlobalWorkPhase(:,l,1)
        WorkPhase(:,j,2)=GlobalWorkPhase(:,l,2)
    endif
enddo

do i=1,nAnisP
    HetEnd=-1
    HetStart=-1
    WorkRight(:)=9
    WorkLeft(:)=9
    do e=1,2
        nRec=0
        nRecomb(i)=nRec
        StR=-9
        EnR=-9
        StRNarrow=-9
        EnRNarrow=-9

        PatMat=e
        SireDamRL=e+1
        CountLeftSwitch=0
        CountRightSwitch=0

        PedId= ped%pedigree(i)%getSireDamNewIDByIndex(SireDamRL) !get index of sire/dam
        if (ped%isDummy(pedId)) cycle
        if ((inputParams%SexOpt==1).and.(ped%pedigree(pedId)%gender==inputParams%HetGameticStatus)) cycle
        if ((PedId>0).and.((float(count(ImputePhase(PedId,:,:)==9))/(2*dimSnps))<0.30)) then          !(RecIdHDIndex(PedId)==1)
            WorkRight=9
            RSide=9
            do j=nSnpBegin,nSnpFinal
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.(ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetStart=j
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,1)) then
                        WorkRight(HetStart)=1
                        RSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetStart,PatMat)==ImputePhase(PedId,HetStart,2)) then
                        WorkRight(HetStart)=2
                        RSide=2
                        exit
                    endif
                endif
            enddo
            if (RSide/=9) then
                do j=HetStart+1,nSnpFinal
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,RSide)).and.(ImputePhase(PedId,j,RSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
                        RSide=abs((RSide-1)-1)+1
                        CountRightSwitch=CountRightSwitch+1
                    endif
                    WorkRight(j)=RSide
                enddo
            endif

            LSide=9
            do j=nSnpFinal,nSnpBegin,-1
                if ((ImputePhase(PedId,j,1)/=ImputePhase(PedId,j,2)).and.(ImputePhase(PedId,j,1)/=9).and.(ImputePhase(PedId,j,2)/=9))  then
                    HetEnd=j
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,1)) then
                        WorkLeft(HetEnd)=1  !£$$$$
                        LSide=1
                        exit
                    endif
                    if (ImputePhase(i,HetEnd,PatMat)==ImputePhase(PedId,HetEnd,2)) then
                        WorkLeft(HetEnd)=2  !£$$$$
                        LSide=2
                        exit
                    endif
                endif
            enddo
            if (LSide/=9) then
                do j=HetEnd-1,1,-1
                    if ((ImputePhase(i,j,PatMat)/=ImputePhase(PedId,j,LSide)).and.(ImputePhase(PedId,j,LSide)/=9).and.(ImputePhase(i,j,PatMat)/=9)) then
                        LSide=abs((LSide-1)-1)+1
                        CountLeftSwitch=CountLeftSwitch+1
                    endif
                    WorkLeft(j)=LSide
                enddo
            endif
            StartDis=-9
            EndDis=-9
            TempVec=9
            LengthVec=0.0

            !prototype start

            StartDis=nSnpBegin
            StartDisOld=nSnpBegin
            EndDis=nSnpFinal
            nRec=0
            nRecomb(i)=nRec
            SuperJ=nSnpFinal
            TempVec=9
            LengthVec=0.0

            StartDisPrev=StartDis
            EndDisPrev=-9

            do while (SuperJ<nSnpFinal)
                SuperJ=SuperJ+1

                !Finding StartDis and Moving it left and EndDis and Movie it Right

                !Find StartDis
                StartDisFound=0
                if (abs(WorkLeft(SuperJ)-WorkRight(SuperJ))==1) then
                    StartDisOld=StartDis
                    StartDis=SuperJ
                    StartDisFound=1
                endif

                if (StartDisFound==1) then
                    nRec=nRec+1
                    nRecomb(i)=nRec
                    print *, i,nRec

                    !Find EndDis
                    EndDisFound=0
                    do k=StartDis+1,nSnpFinal
                        if (WorkLeft(k)==WorkRight(k)) then
                            EndDis=k
                            SuperJ=k
                            EndDisFound=1
                            exit
                        endif
                    enddo
                    if (EndDisFound==0) then
                        EndDis=nSnpFinal
                        SuperJ=nSnpFinal
                    endif

                    !Move StartDis to the left informative marker
                    StartDisTmp=StartDis
                    do k=StartDis,nSnpBegin,-1
                        if ((WorkPhase(PedId,k,1)+WorkPhase(PedId,k,2))==1) then
                            if (WorkPhase(i,k,e)/=9) then
                                StartDis=k
                                exit
                            endif
                        endif
                    enddo
                    if (StartDis<=StartDisOld) StartDis=StartDisTmp
                    if (StartDis<1) StartDis=1
                    if (StartDis<EndDisPrev) StartDis=EndDisPrev

                    !Move EndDis to the right informative marker
                    do k=EndDis,nSnpFinal
                        if ((WorkPhase(PedId,k,1)+WorkPhase(PedId,k,2))==1) then
                            if (WorkPhase(i,k,e)/=9) then
                                EndDis=k
                                SuperJ=k
                                exit
                            endif
                        endif
                    enddo

                    StR(nRec)=StartDis
                    EnR(nRec)=EndDis
                    LengthVec(StartDis:EndDis)=1.0/(EndDis-StartDis)
                    StartDisPrev=StartDis
                    EndDisPrev=EndDis

                endif
            enddo
            !prototype end (temp)
            StR(nRec+1)=-9
            EnR(nRec+1)=-9

            GamA=9
            GamB=9
            k=1
            do j=nSnpBegin,nSnpFinal

                if (j==StR(k)) then
                    k=k+1
                    GamA=WorkRight(j)

                    If (GamA==1) GamB=2
                    If (GamA==2) GamB=1
                endif
                if (GamA==9) cycle
            enddo
        endif
    enddo
    WorkPhase(i,:,:)=ImputePhase(i,:,:)
    ! nRecomb(i)=nRec
enddo

deallocate(WorkPhase)
deallocate(TempVec)
deallocate(LengthVec)
deallocate(WorkLeft)
deallocate(WorkRight)
deallocate(TempWork)
deallocate(PatAlleleProb)
deallocate(MatAlleleProb)
deallocate(GeneProbWork)
deallocate(StR)
deallocate(EnR)
deallocate(StRNarrow)
deallocate(EnRNarrow)

end subroutine ModelRecombSegment

!######################################################################
subroutine RemoveRecombinationForSegment
use Global
use GlobalVariablesHmmMaCH
use Utils
use AlphaImputeInMod

implicit none
integer :: StartSnp, StopSnp, nSegments, SegmentSize
integer,allocatable :: nRecomb(:)
integer :: i,k,indv
type(AlphaImputeInput), pointer :: inputParams

inputParams => defaultInput

allocate(nRecomb(nAnisP))
StartSnp=1
StopSnp=nSnpHmm
nSegments=nSnpHmm/inputParams%windowLength
if (MOD(nSnpHmm,inputParams%windowLength)/=0) nSegments=nSegments+1
SegmentSize=inputParams%windowLength

do i=1,nSegments
    nRecomb=0
    StartSnp=SegmentSize*(i-1)+1
    StopSnp=SegmentSize*i
    if (StopSnp>nSegments) StopSnp=nSnpHmm

    call ModelRecombSegment(StartSnp, StopSnp, nRecomb)
    k=0
    do indv=1,nAnisP
        if (ped%pedigree(i)%genotyped==1 .AND. nRecomb(indv)>1) then
            k=k+1

            call RemoveGenotypeInformationIndividualSegment(k,StartSnp,StopSnp)
        endif
    enddo
enddo

end subroutine RemoveRecombinationForSegment

END MODULE