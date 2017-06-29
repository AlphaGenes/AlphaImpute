module HeuristicGeneprobModule
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

        integer :: maxVal,maxPos, max2Val, max2Pos

        real, dimension(:,:,:), allocatable :: genosProbs !< array for each animal, 0,1,2
        real, dimension(0:2, 0:2) :: heuristicTrace
        logical :: sire
        
        log1 = 0
        loge = log(1e-5)
        log5 = log(0.5)

        heuristicTrace(:, 0) = [log1, log5, loge]
        heuristicTrace(:, 1) = [loge, log5, log1]
        heuristicTrace(:, 2) = [loge, log5, log1]


        allocate(GenosProbs(ped%pedigreeSize,inputParams%nsnp, 0:2))

        do i=1, ped%pedigreeSize
            ind => ped%pedigree(i)
            if (ind%Genotyped) cycle
            

           
            do j=1, inputParams%nsnp
                

                do h=1, ind%nOffs
                    off => ind%offsprings(h)%p

                    if (off%genotyped) then


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

                        if (g == 1) then
                            
                            p = off%individualPhase(e)%getPhase(j)
                            if (p  /=9) then
                                genosProbs(i,j,:) = genosProbs(i,j,:) + heuristicTrace(p, :)
                        
                            endif
                        else if (g == 2 .or. g == 0) then
                            genosProbs(i,j,:) = genosProbs(i,j,:) + heuristicTrace(g, :)

                        else
                            if (sire) then
                                g = off%damPointer%individualGenotype%getgenotype(j)

                            else
                                g = off%sirePointer%individualGenotype%getgenotype(j)
                            endif

                            if (g == 9) cycle !< don't care if animal is missing
                            genosProbs(i,j,:) = genosProbs(i,j,:) + heuristicTrace(g, :)

                        endif                    
                    endif
                enddo !< offsprings
            
                maxPos = 0
                maxVal = genosProbs(i,j,0) 
                max2Val = 0

                do h=1,2
                    if (genosProbs(i,j,h) > maxVal ) then
                        max2Val = maxVal
                        maxVal = genosProbs(i,j,h)
                        maxPos = h
                        
                    else if (genosProbs(i,j,h) > max2val ) then
                        max2Val = genosProbs(i,j,h)
                    endif
                enddo !< genotype probs

                    if ((maxVal - max2val) > 15) then
                        call ped%pedigree(i)%individualGenotype%setGenotype(j, maxPos)
                    endif
            enddo   !<SNPS!
        enddo !< animals

    end subroutine HeuristicGeneprob



    subroutine outputGenosProbs(filename,GenosProbs, ped, nsnp )

        use PedigreeModule

        character(len=*), intent(in) :: filename

        real, dimension(:,:,:), allocatable, intent(in) :: genosProbs
        integer, intent(in) :: nsnp

        type(PedigreeHolder),intent(in) :: ped
        character(:), allocatable :: rowfmt

        integer :: i,j,unit

        open(newunit=unit, file=filename)
        WRITE(rowfmt,'(A,I9,A)') '(a,',nsnp+10,'f10.4)'
        do i=1, ped%pedigreeSize

            do j=1, nsnp                    
                write(unit,rowfmt) ped%pedigree(i)%originalID, genosProbs(i,j,1)
                write(unit,rowfmt) ped%pedigree(i)%originalID, genosProbs(i,j,2)
            enddo

        enddo

        close(unit)

    end subroutine outputGenosProbs



          subroutine outputGenotypesTestNew(filename, ped, nsnps )

          use PedigreeModule
         character(len=*), intent(in) :: filename

        integer, intent(in) :: nsnps

        integer :: consensusFile,i

        type(PedigreeHolder),intent(in) :: ped
        character(len=300) :: rowfmt

        open(newunit=consensusFile, file=filename, status="unknown")
        WRITE(rowfmt,'(A,I9,A)') '(a,',nSnps+10,'I2)'
        do i = 1, ped%pedigreeSize
            
            write(consensusFile, rowfmt) ped%pedigree(i)%originalID,ped%pedigree(i)%individualGenotype%toIntegerArray()
        enddo


        end subroutine outputGenotypesTestNew 


end module HeuristicGeneprobModule 
