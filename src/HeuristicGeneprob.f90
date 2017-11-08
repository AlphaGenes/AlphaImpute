module HeuristicGeneprobModule

    use iso_fortran_env
    implicit none
    contains


   subroutine HeuristicGeneprob(inputParams, ped)

        use PedigreeModule
        use AlphaImputeSpecFileModule

        type(AlphaImputeInput) :: inputParams
        type(PedigreeHolder) :: ped
        real :: log5, loge, log1
        type(individual), pointer :: ind, off
        integer :: i, j, h
        
        integer :: g !< genotype at snp
        integer :: e !< if ind is childs sire (1) or dam (2)
        integer :: p !< offspring phase

        integer :: maxVal,maxPos, max2Val,gOtherParent

        real(kind=real64), dimension(0:2) :: genosProbs !< array for each animal, 0,1,2
        real, dimension(0:2, 0:2) :: heuristicTrace
        logical :: sire
        integer(kind=1) :: childSeg, childPhase 
        
        log1 = 0
        loge = log(1e-5)
        log5 = log(0.5)

        heuristicTrace(:, 0) = [log1, log5, loge]
        heuristicTrace(:, 1) = [loge, log5, log1]
        heuristicTrace(:, 2) = [loge, log5, log1]


        do i=1, ped%pedigreeSize
            ind => ped%pedigree(i)            

            do j=1, inputParams%nsnp
                genosProbs = 0

                if(ind%individualGenotype%getGenotype(j)/=9) cycle
                
                do h=1, ind%nOffs
                    off => ind%offsprings(h)%p


                    g = off%individualGenotype%getGenotype(j)
                    !< if snp is set at this pos, skip 
                    if (g ==9) cycle
                    ! figure out if this animal is sire or dam of this offspring
                    if (associated(off%sirePointer, ind)) then
                        e = 1
                        sire = .true.
                    else
                        e = 2
                        sire = .false.
                    endif

                    if (g == 1) then
                        
                        p = off%individualPhase(e)%getPhase(j)
                        if (p  /=9) then
                            genosProbs(:) = genosProbs(:) + heuristicTrace(:, p)
                        else
                            !See if we can implicitly phase them.
                            if (sire) then
                                gOtherParent = off%damPointer%individualGenotype%getgenotype(j)
                            else
                                gOtherParent = off%sirePointer%individualGenotype%getgenotype(j)
                            endif

                            if (gOtherParent ==0 .or. gOtherParent == 1) then
                                p = 1 - gOtherParent/2 !This takes 0 to 0 and 2 to 1
                                genosProbs(:) = genosProbs(:) + heuristicTrace(:, p)
                            endif

                        endif   
                    else if (g == 2 .or. g == 0) then
                        genosProbs(:) = genosProbs(:) + heuristicTrace(:, g)
                    else
                        print *, "Andrew says this never should get hit"     
                    endif                    
                enddo !< offsprings
            
                maxPos = 0
                maxVal = genosProbs(0) 
                max2Val = 0

                do h=1,2
                    if (genosProbs(h) > maxVal ) then
                        max2Val = maxVal
                        maxVal = genosProbs(h)
                        maxPos = h
                        
                    else if (genosProbs(h) > max2val ) then
                        max2Val = genosProbs(h)
                    endif
                enddo !< genotype probs

                    if ((maxVal - max2val) > 15) then
                        call ped%pedigree(i)%individualGenotype%setGenotype(j, maxPos)
                    endif
            enddo   !<SNPS!
        enddo !< animals

    end subroutine HeuristicGeneprob



    subroutine heuristicMLP(ped)

        use PedigreeModule

        type(PedigreeHolder), intent(inout) :: ped
        call peelDown(ped)
        call updateSeg(ped)
        call peelUp(ped)


        contains

            subroutine peelUp(ped)

                use PedigreeModule
                use AlphaImputeSpecFileModule

                type(AlphaImputeInput) :: inputParams
                type(PedigreeHolder) :: ped
                real :: log5, loge, log1,logd
                type(individual), pointer :: ind, off
                integer :: i, j, h
                
                integer(kind=1) :: g !< genotype at snp
                integer(kind=1) :: e !< if ind is childs sire (1) or dam (2)
                integer(kind=1) :: p !< offspring phase

                integer :: maxVal,maxPos, max2Val,max2Pos

                real(kind=real64), dimension(1:4) :: phaseProbs !< array for each animal, 0,1,2
                real(kind=real64), dimension(0:2) :: genosProbs !< array for each animal, 0,1,2
                real, dimension(1:4, 0:1, 1:2) :: heuristicTraceSeg !allele, phase, seg 
                real, dimension(1:4, 0:2) :: heuristicTraceNoSeg
                logical :: sire
                integer(kind=1) :: childSeg, childPhase
                
                log1 = 0
                loge = log(1e-5)
                log5 = log(0.5)
                logd = loge

                heuristicTraceNoSeg(:, 0) = [log1, log5, log5, loge]
                heuristicTraceNoSeg(:, 1) = [loge, log5, log5, log1]
                heuristicTraceNoSeg(:, 2) = [loge, log5, log5, log1]

                heuristicTraceSeg(:, 0, 1) = [log1, log1, logd, logd] !seg = 0
                heuristicTraceSeg(:, 1, 1) = [logd, logd, log1, log1]
                heuristicTraceSeg(:, 0, 2) = [log1, logd, log1, logd] !seg = 1
                heuristicTraceSeg(:, 1, 2) = [logd, log1, logd, log1]

                do i=1, ped%pedigreeSize
                    ind => ped%pedigree(i)
                    if (ind%Genotyped) cycle
                    

                   
                    do j=1, inputParams%nsnp
                        phaseProbs = 0

                        do h=1, ind%nOffs
                            off => ind%offsprings(h)%p



                            g = off%individualGenotype%getGenotype(j)
                            !< if snp is set at this pos, skip 
                            if (g /=9) cycle
                            ! figure out if this animal is sire or dam of this offspring
                            if (associated(off%sirePointer, ind)) then
                                e = 1
                                sire = .true.
                            else
                                e = 2
                                sire = .false.
                            endif
                            childSeg = ped%pedigree(i)%seg(j,e)
                            if(childSeg /=9) then
                                childPhase = ped%pedigree(i)%individualPhase(e)%getPhase(j)
                                if(childPhase /=9) then
                                    phaseProbs(:) = phaseProbs(:) + heuristicTraceSeg(:, childPhase, childSeg)
                                endif
                            else
                                if (g == 1) then
                                    
                                    p = off%individualPhase(e)%getPhase(j)
                                    if (p  /=9) then
                                        phaseProbs(:) = phaseProbs(:) + heuristicTraceNoSeg(:, g)
                                
                                    endif
                                else if (g == 2 .or. g == 0) then
                                    phaseProbs(:) = phaseProbs(:) + heuristicTraceNoSeg(:, g)

                                else
                                    if (sire) then
                                        g = off%damPointer%individualGenotype%getgenotype(j)

                                    else
                                        g = off%sirePointer%individualGenotype%getgenotype(j)
                                    endif

                                    if (g == 9) cycle !< don't care if animal is missing
                                    phaseProbs(:) = phaseProbs(:) + heuristicTraceNoSeg(:,g)

                                endif        
                            end if            
                        enddo !< offsprings
                    

                        !First can we phase?

                        maxPos = 1
                        maxVal = phaseProbs(1) 
                        max2Val = 0
                        max2Pos = 1

                        do h=2,4
                            if (phaseProbs(h) > maxVal ) then
                                max2Val = maxVal
                                max2Pos = maxPos
                                maxVal = phaseProbs(h)
                                maxPos = h
                                
                            else if (phaseProbs(h) > max2val ) then
                                max2Val = phaseProbs(h)
                                max2Pos = h
                            endif
                        enddo !< genotype probs

                        if ((maxVal - max2val) > 15-1) then
                            select case( maxPos )
                                case(1)
                                    call ped%pedigree(i)%individualGenotype%setGenotype(j, 0)
                                    call ped%pedigree(i)%individualPhase(1)%setPhase(j, 0)
                                    call ped%pedigree(i)%individualPhase(2)%setPhase(j, 0)
                                case(2)
                                    call ped%pedigree(i)%individualGenotype%setGenotype(j, 1)
                                    call ped%pedigree(i)%individualPhase(1)%setPhase(j, 0)
                                    call ped%pedigree(i)%individualPhase(2)%setPhase(j, 1)
                                case(3)
                                    call ped%pedigree(i)%individualGenotype%setGenotype(j, 1)
                                    call ped%pedigree(i)%individualPhase(1)%setPhase(j, 1)
                                    call ped%pedigree(i)%individualPhase(2)%setPhase(j, 0)
                                case(4)
                                    call ped%pedigree(i)%individualGenotype%setGenotype(j, 2)
                                    call ped%pedigree(i)%individualPhase(1)%setPhase(j, 1)
                                    call ped%pedigree(i)%individualPhase(2)%setPhase(j, 1)
                            end select
                        else
                            genosProbs(0:2) = [phaseProbs(1), max(phaseProbs(2),phaseProbs(3)), phaseProbs(4)]
                            maxPos = 0
                            maxVal = genosProbs(0) 
                            max2Val = 0

                            do h=1,2
                                if (genosProbs(h) > maxVal ) then
                                    max2Val = maxVal
                                    maxVal = genosProbs(h)
                                    maxPos = h
                                    
                                else if (genosProbs(h) > max2val ) then
                                    max2Val = genosProbs(h)
                                endif
                            enddo !< genotype probs

                            if ((maxVal - max2val) > 15) then
                                call ped%pedigree(i)%individualGenotype%setGenotype(j, maxPos)
                            endif                        
                        endif
                    enddo   !<SNPS!
                enddo !< animals
            end subroutine peelup


            subroutine updateSeg(ped)
                use ConstantModule
                use PedigreeModule

                implicit none
                type(PedigreeHolder) :: ped
                integer(kind=1) :: geno(2),childPhase, phase(2)
                integer :: startLoc 
                integer :: lastLoc
                integer :: currentSeg
                integer :: numConsistent
                integer :: threshold
                integer :: i,j,e, nsnps

                type(individual), pointer :: parent

                
                threshold = 3
                nsnps = ped%pedigree(ped%genotypeMap(1))%individualPhase(1)%length

                do i = 1, ped%pedigreeSize

                    if (ped%pedigree(i)%founder) cycle
                    do j=1, nSnps

                        geno(1) = ped%pedigree(i)%sirePointer%individualGenotype%getgenotype(j) 
                        geno(2) = ped%pedigree(i)%damPointer%individualGenotype%getgenotype(j) 

                        ! if parents are homozygous
                        if (MODULO(geno(1),2) == 0 .and. MODULO(geno(2),2) == 0) cycle

                        if (ped%pedigree(i)%individualPhase(1)%getPhase(j) /= 9 .and. ped%pedigree(i)%individualPhase(2)%getPHase(j) /= 9) then

                            do e=1,2
                                parent => ped%pedigree(i)%getSireDamObjectByIndex(e+1)
                                if (geno(e) == 1) then
                                    phase(1) = parent%individualPhase(1)%getPhase(j) 
                                    phase(2) = parent%individualPhase(2)%getPhase(j) 
                                    if (phase(1) /=9 .and. phase(2) /= 9) then

                                        childPhase = ped%pedigree(i)%individualPhase(e)%getPhase(j)
                                        if (childPhase == phase(1)) then
                                            call ped%pedigree(i)%setSeg(j,e,1)
                                        else if (childPhase == phase(2)) then
                                            call ped%pedigree(i)%setSeg(j,e,2)
                                        else
                                            write(error_unit, *) "WARNING: GENO and PHASE don't match: ", childphase, phase(1), phase(2), geno(e),parent%individualGenotype%getGenotype(j) 
                                        endif
                                    endif
                                endif
                            enddo
                        endif
                    enddo
                    do e =1,2
                         startLoc = 1
                        lastLoc = 1
                        currentSeg = -1 
                        numConsistent = 0 !< after seeing first allele, this will be 1
                        

                        do j =1, nSnps
                            if (ped%pedigree(i)%seg(j,e) /=9) then
                                if (ped%pedigree(i)%seg(j,e) == currentSeg .or. currentSeg == -1) then
                                    lastLoc = j
                                    numConsistent = numConsistent +1
                                else
                                    if (numConsistent > threshold) then
                                        ped%pedigree(i)%seg(startLoc:lastLoc,e) = currentSeg
                                    else 
                                        ped%pedigree(i)%seg(startLoc:lastLoc,e) = MISSINGGENOTYPECODE
                                    endif
                                    
                                    startLoc = j
                                    lastLoc = j
                                    currentSeg = ped%pedigree(i)%seg(j,e)
                                    numConsistent = 1
                                endif
                            endif
                        enddo

                        lastLoc = nSnps
                        if (numConsistent > threshold) then
                            ped%pedigree(i)%seg(startLoc:lastLoc,e) = currentSeg
                        else 
                            ped%pedigree(i)%seg(startLoc:lastLoc,e) = MISSINGGENOTYPECODE
                        endif
                    enddo
                enddo 

            end subroutine updateSeg


            subroutine peelDown(ped)
                use PedigreeModule
                use ConstantModule
                implicit none

                type(PedigreeHolder), intent(inout) :: ped
                integer :: segregation
                integer :: i,j,e,nsnps


                integer(kind=1) :: parentphase(2,2), parentGenotype(2)



                nsnps = ped%pedigree(ped%genotypeMap(1))%individualPhase(1)%length
                do i=1, ped%pedigreeSize

                    call ped%pedigree(i)%setSegToMissing
                    if (ped%pedigree(i)%founder) cycle

                    do j=1, nsnps

                        parentphase(1,1) = ped%pedigree(i)%sirePointer%individualPhase(1)%getPhase(j)
                        parentphase(1,2) = ped%pedigree(i)%sirePointer%individualPhase(2)%getPhase(j)
                        parentphase(2,1) = ped%pedigree(i)%damPointer%individualPhase(1)%getPhase(j)
                        parentphase(2,2) = ped%pedigree(i)%damPointer%individualPhase(2)%getPhase(j)
                        parentGenotype(1) = ped%pedigree(i)%sirePointer%individualGenotype%getgenotype(j)
                        parentGenotype(2) = ped%pedigree(i)%damPointer%individualGenotype%getgenotype(j)

                        do e=1,2

                            if (parentGenotype(e) == 2 .or. parentGenotype(e) == 0) then
                                call ped%pedigree(i)%individualPhase(e)%setPhase(j,parentGenotype(e) /2_ONEBYTEINT)
                            endif

                            if (parentGenotype(e) == 1) then

                                if (parentphase(e,1) /= 9 .and. parentphase(e,2) /= 9) then
                                    segregation = ped%pedigree(i)%getSeg(j,e)
                                    if (segregation /=9) then
                                       call ped%pedigree(i)%individualPhase(e)%setPhase(j,parentphase(e,segregation))
                                    endif
                                endif
                            endif
                        enddo


                    enddo

                    call ped%pedigree(i)%makeIndividualGenotypeFromPhase
                    call ped%pedigree(i)%makeIndividualPhaseCompliment
                     
                enddo
            end subroutine peelDown

    end subroutine heuristicMLP




end module HeuristicGeneprobModule 
