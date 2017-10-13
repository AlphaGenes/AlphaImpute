#ifdef _WIN32

#define STRINGIFY(x)#x
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


#else

#define STRINGIFY(x)#x
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


#endif

!#############################################################################################################################################################################################################################

module AlphaImputeModule
	implicit none
	contains




		subroutine PhasingManagementNew(results)
			use AlphaImputeSpecFileModule
			use omp_lib
			use AlphaPhaseParametersModule
			use AlphaPhaseFunctions
			use AlphaPhaseResultsModule
			use Global, only : ped, OPT_RESTART_PHASING
			use inputoutput, only : MakeDirectories
			use OutputParametersModule
			implicit none

			type(AlphaImputeInput), pointer :: inputParams
			type(AlphaPhaseParameters) :: params
			type(AlphaPhaseResultsContainer), intent(out) :: results
			integer :: nCoreLengths,i, coreIndexes
			type(OutputParameters) :: oParams
			type(PedigreeHolder), allocatable :: hdPed

			inputParams=> defaultInput
			nCoreLengths = size(inputParams%CoreAndTailLengths)
			results%nResults = nCoreLengths*2

			allocate(results%results(nCoreLengths*2))
			params = newParameters()
			params%iterateType = inputParams%iterateMethod
			params%iterateNumber = inputParams%PhaseSubsetSize
			params%numIter = inputParams%PhaseNIterations
			params%minOverlap = inputparams%minoverlaphaplotype
			params%percGenoHaploDisagree = inputparams%GenotypeErrorPhase*0.01
			params%Offset = .true.

			if (inputparams%minoverlaphaplotype /= 0) then
				params%percMinPresent = 0
			endif
			call omp_set_nested(.true.)



			write(6,*) " "
			write(6,*) " ", "Running AlphaPhase"


			allocate(hdPed)
			call ped%getHDPedigree(hdPed)

			!$OMP parallel do schedule(dynamic)&
			!$OMP default(shared) &
			!$OMP FIRSTPRIVATE(params) &
			!$OMP PRIVATE(coreIndexes, i,oParams)
			do i= 1,nCoreLengths*2
				coreIndexes = i

				if (i > nCoreLengths) Then
					params%offset = .false.
					coreIndexes = i - nCoreLengths
				endif




				params%CoreAndTailLength = inputParams%CoreAndTailLengths(coreIndexes)
				params%jump = inputParams%CoreLengths(coreIndexes)
				params%numsurrdisagree = 1
				params%useSurrsN = 10
				results%results(i) = phaseAndCreateLibraries(hdPed, params, quiet=.true., updatePedigree=.false.)
				if (inputParams%restartOption==OPT_RESTART_PHASING) then
					oParams = newOutputParametersImpute()
					write(oParams%outputDirectory,'("."a"Phasing",a,"Phase"i0)') DASH,DASH, i
					call writeAlphaPhaseResults(results%results(i), hdPed, oParams)
				endif

			enddo
			!$omp end parallel do
			write(6,*) " ", "Finished Running AlphaPhase"

			deallocate(hdPed)
			if (inputParams%restartOption==OPT_RESTART_PHASING) then
				write(6,*) "Restart option 1 stops program after Phasing has been managed"
				stop
			endif

		end subroutine PhasingManagementNew




#ifdef MPIACTIVE
		subroutine phasingManagementCluster

			use AlphaImputeSpecFileModule
			use AlphaPhaseParametersModule
			use AlphaPhaseFunctions
			use AlphaPhaseResultsModule
			use Global, only : ped, OPT_RESTART_PHASING
			use OutputParametersModule
			use InputOutput
			use MPIUtilities

			implicit none

			integer :: totalToDo
			type(AlphaImputeInput), pointer :: inputParams
			type(AlphaPhaseParameters) :: params
			integer :: nCoreLengths,i, coreIndexes
			type(OutputParameters) :: oParams
			integer :: index
			type(AlphaPhaseResults) :: result
			inputParams=> defaultInput
			nCoreLengths = size(inputParams%CoreAndTailLengths)

			params = newParameters()
			params%iterateType = inputParams%iterateMethod
			params%iterateNumber = inputParams%PhaseSubsetSize
			params%numIter = inputParams%PhaseNIterations
			params%minOverlap = inputparams%minoverlaphaplotype
			params%percGenoHaploDisagree = inputparams%GenotypeErrorPhase*0.01
			params%Offset = .true.
			oParams = newOutputParametersImpute()

			if (inputparams%minoverlaphaplotype /= 0) then
				params%percMinPresent = 0
			endif

			call initialiseMPI
			print *,"MPI SIZE:", mpiSize
			totalToDo = (nCoreLengths*2)/mpiSize
			print *,"TOTALTODO", totalTOdo
			if (totalToDo <1) then
				if (mpiRank+1 < nCoreLengths*2) return
			endif
			write(6,*) " "
			write(6,*) " ", "Running AlphaPhase"

			! TODO can omp this loop
			do i=1,totalToDo

				index = (mpiRank+1)+((i-1) *nCoreLengths)
				if (index > nCoreLengths) Then
					params%offset = .false.
					coreIndexes = index - nCoreLengths
				else
					coreIndexes = index
				endif

				params%CoreAndTailLength = inputParams%CoreAndTailLengths(coreIndexes)
				params%jump = inputParams%CoreLengths(coreIndexes)
				params%numsurrdisagree = 1
				params%useSurrsN = 10
				result = phaseAndCreateLibraries(ped, params, quiet=.true., updatePedigree=.false.)

				print *,"Finished running AlphaPhase run ",index
				write(oParams%outputDirectory,'("."a"Phasing",a,"Phase"i0)') DASH,DASH, index
				call writeAlphaPhaseResults(result, ped, oParams)
			enddo


			call endMPI

			stop
		end subroutine phasingManagementCluster

#endif
		!#############################################################################################################################################################################################################################

		subroutine IterateInsteadOfGeneProbs
			use Global
			use Imputation
			use ConstantModule
			use AlphaImputeSpecFileModule
			implicit none

			integer :: e,i,j,k,Counter,ParId

			integer(kind=1) :: genos
			real,allocatable,dimension(:) :: TempAlleleFreq

			inputParams => defaultInput
			if (inputParams%SexOpt==1) then
				if (inputParams%outopt==0) nSnpIterate=inputParams%nsnp
				if (inputParams%outopt==1) nSnpIterate=inputParams%nSnpRaw


				if (.not. allocated(ProbImputeGenos)) allocate(ProbImputeGenos(0:ped%pedigreeSize,nSnpIterate))
				if (.not. allocated(ProbImputePhase)) allocate(ProbImputePhase(0:ped%pedigreeSize,nSnpIterate,2))
				if (.not. allocated(TempAlleleFreq)) allocate(TempAlleleFreq(nSnpIterate))

				TempAlleleFreq=0.0

				do j=1,size(SnpIncluded)
					Counter=0
					do i=1,ped%pedigreeSize-ped%nDummys

						genos = ped%pedigree(i)%individualGenotype%getGenotype(SnpIncluded(j))
						if (genos/=MISSINGGENOTYPECODE) then
							TempAlleleFreq(j)=TempAlleleFreq(j)+genos
							Counter=Counter+2
						endif

					enddo
					if (Counter/=0) then
						TempAlleleFreq(j)=TempAlleleFreq(j)/Counter
					else
						TempAlleleFreq(j)=9.0
					endif
				enddo

				ProbImputePhase(0,:,1)=TempAlleleFreq(:)
				ProbImputePhase(0,:,2)=TempAlleleFreq(:)
				ProbImputeGenos(0,:)=2*TempAlleleFreq(:)
				ProbImputeGenos(1:ped%pedigreeSize,:)=-9.0
				ProbImputePhase(1:ped%pedigreeSize,:,:)=-9.0

				call IterateParentHomoFill
				call ped%PhaseComplement
				call ped%MakeGenotype

				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,nSnpIterate
						do e=1,2
							block
								integer(kind=1) :: phase
								phase = ped%pedigree(i)%individualPhase(e)%getPhase(SnpIncluded(j))
								if (phase/=9) ProbImputePhase(i,j,e)=float(phase)
							end block
						enddo
					enddo
				enddo

				do i=1,ped%pedigreeSize-ped%nDummys
					do e=1,2
						parId = ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
						if (ParId==0) then
							do j=1,nSnpIterate
								block
									integer(kind=1) :: phase
									phase = ped%pedigree(i)%individualPhase(e)%getPhase(SnpIncluded(j))
									if (phase==9 .and. TempAlleleFreq(j)/=9) ProbImputePhase(i,j,e)=TempAlleleFreq(j)
								end block
							enddo
						endif
						if (ped%pedigree(i)%gender==inputParams%HomGameticStatus) then
							do j=1,nSnpIterate
								block
									integer(kind=1) :: phase
									phase = ped%pedigree(i)%individualPhase(e)%getPhase(SnpIncluded(j))
									if (phase==9) then
										ProbImputePhase(i,j,e)=(sum(ProbImputePhase(ParId,j,:))/2)
									endif
								end block
							enddo
						endif
					enddo
					if (ped%pedigree(i)%gender==inputParams%hetGameticStatus) then
						ParId=ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%hetGameticStatus+1)
						do j=1,nSnpIterate

							if (ped%pedigree(i)%individualPhase(1)%getPhase(SnpIncluded(j))==9) then
								ProbImputePhase(i,j,:)=(sum(ProbImputePhase(ParId,j,:))/2)
							endif
						enddo
					endif
				enddo

				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,nSnpIterate
						do k=1,2
							block
								integer(kind=1) :: phase
								phase = ped%pedigree(i)%individualPhase(k)%getPhase(SnpIncluded(j))
								if (phase/=9) ProbImputePhase(i,j,k)=float(phase)
							end block
						enddo
						block
							integer(kind=1) :: genos
							genos = ped%pedigree(i)%individualGenotype%getGenotype(SnpIncluded(j))
							if (genos/=9) then
								ProbImputeGenos(i,j)=float(genos)
							else
								ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
							endif
						end block
					enddo
				enddo

				if (inputParams%SexOpt==1) then
					do i=1,ped%pedigreeSize-ped%nDummys
						if (ped%pedigree(i)%gender==inputParams%hetGameticStatus) then
							! setFromOtherIfMissing
							call ped%pedigree(i)%individualPhase(1)%setFromOtherIfMissing(ped%pedigree(i)%individualPhase(2))
							call ped%pedigree(i)%individualPhase(2)%setFromOtherIfMissing(ped%pedigree(i)%individualPhase(1))
							do j=1,nSnpIterate
								if ((ProbImputePhase(i,j,1)==-9.0).and.(ProbImputePhase(i,j,2)/=-9.0)) ProbImputePhase(i,j,1)=ProbImputePhase(i,j,2)
								if ((ProbImputePhase(i,j,2)==-9.0).and.(ProbImputePhase(i,j,1)/=-9.0)) ProbImputePhase(i,j,2)=ProbImputePhase(i,j,1)
							enddo
						endif
					enddo
				endif

				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,nSnpIterate
						if (ProbImputeGenos(i,j)==-9.0) ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
						if (ProbImputeGenos(i,j)>1.999)  call ped%pedigree(i)%individualGenotype%setGenotype(SnpIncluded(j),2)
						if (ProbImputeGenos(i,j)<0.0001) call ped%pedigree(i)%individualGenotype%setGenotype(SnpIncluded(j),0)
						if ((ProbImputeGenos(i,j)>0.999).and.(ProbImputeGenos(i,j)<1.00001)) call ped%pedigree(i)%individualGenotype%setGenotype(SnpIncluded(j),1)

					enddo
				enddo


				do j=1,nSnpIterate
					if (TempAlleleFreq(j)==9.0) then
						ProbImputeGenos(:,j)=9.0
						ProbImputePhase(:,j,:)=9.0
					endif
				enddo

			else

				if (inputParams%outopt==0) nSnpIterate=inputParams%nsnp
				if (inputParams%outopt==1) nSnpIterate=inputParams%nSnpRaw

				if (.not. allocated(ProbImputeGenos)) allocate(ProbImputeGenos(0:ped%pedigreeSize,nSnpIterate))
				if (.not. allocated(ProbImputePhase)) allocate(ProbImputePhase(0:ped%pedigreeSize,nSnpIterate,2))
				if (.not. allocated(TempAlleleFreq)) allocate(TempAlleleFreq(nSnpIterate))

				TempAlleleFreq=0.0
				do j=1,nSnpIterate
					Counter=0
					do i=1,ped%pedigreeSize-ped%nDummys
						block

							integer(kind=1) :: geno
							geno = ped%pedigree(i)%individualGenotype%getGenotype(j)
							if (Geno/=9) then
								TempAlleleFreq(j)=TempAlleleFreq(j)+geno
								Counter=Counter+2
							endif
						end block
					enddo
					if (Counter/=0) then
						TempAlleleFreq(j)=TempAlleleFreq(j)/Counter
					else
						TempAlleleFreq(j)=9.0
					endif
				enddo

				ProbImputePhase(0,:,1)=TempAlleleFreq(:)
				ProbImputePhase(0,:,2)=TempAlleleFreq(:)
				ProbImputeGenos(0,:)=2*TempAlleleFreq(:)
				ProbImputeGenos(1:ped%pedigreeSize-ped%nDummys,:)=-9.0
				ProbImputePhase(1:ped%pedigreeSize-ped%nDummys,:,:)=-9.0

				call IterateParentHomoFill
				call ped%PhaseComplement
				call ped%MakeGenotype

				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,nSnpIterate
						do e=1,2
							block
								integer(kind=1) :: phase

								phase = ped%pedigree(i)%individualPhase(e)%getPhase(j)
								if (phase /=9) ProbImputePhase(i,j,e)=float(phase)
							end block
						enddo
					enddo
				enddo

				do i=1,ped%pedigreeSize-ped%nDummys
					do e=1,2
						parID=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
						if (ParId==0) then
							do j=1,nSnpIterate
								block
									integer :: phase

									phase = ped%pedigree(i)%individualPhase(e)%getPhase(j)
									if (phase==9) ProbImputePhase(i,j,e)=TempAlleleFreq(j)
								end block
							enddo
						endif
						do j=1,nSnpIterate
							if (ped%pedigree(i)%individualPhase(e)%isMissing(j)) then
								ProbImputePhase(i,j,e)=(sum(ProbImputePhase(ParId,j,:))/2)
							endif
						enddo
					enddo
				enddo

				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,nSnpIterate
						do k=1,2
							block
								integer(kind=1) :: phase
								phase = ped%pedigree(i)%individualPhase(k)%getPhase(j)
								if (phase/=9) ProbImputePhase(i,j,k)=float(phase)
							end block
						enddo

						block
							integer(kind=1) :: geno

							geno = ped%pedigree(i)%individualGenotype%getGenotype(j)
							if (geno/=9) then
								ProbImputeGenos(i,j)=float(geno)
							else
								ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
							endif
						end block
					enddo
				enddo

				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,nSnpIterate
						if (ProbImputeGenos(i,j)==-9.0) ProbImputeGenos(i,j)=sum(ProbImputePhase(i,j,:))
						if (ProbImputeGenos(i,j)>1.999) call ped%pedigree(i)%individualGenotype%setGenotype(j,2)
						if (ProbImputeGenos(i,j)<0.0001) call ped%pedigree(i)%individualGenotype%setGenotype(j,0)
						if ((ProbImputeGenos(i,j)>0.999).and.(ProbImputeGenos(i,j)<1.00001)) call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
					enddo
				enddo

				do j=1,nSnpIterate
					if (TempAlleleFreq(j)==9.0) then
						ProbImputeGenos(:,j)=9.0
						ProbImputePhase(:,j,:)=9.0
					endif
				enddo

			endif

		end subroutine IterateInsteadOfGeneProbs

		!#############################################################################################################################################################################################################################
		subroutine WriteOutResults
			use Global

			use AlphaImputeInputOutputModule
			use AlphaImputeSpecFileModule
			use HMMInputOutputModule

			implicit none

			integer :: i,j,l
			integer(kind=1),allocatable,dimension (:,:) :: ImputeGenosHMM
			integer(kind=1),allocatable,dimension (:,:,:) :: ImputePhaseHMM

			integer(kind=1),allocatable,dimension (:,:) :: tmpGenos
			integer(kind=1),allocatable,dimension (:,:,:) :: tmpPhase
			type(AlphaImputeInput), pointer :: inputParams
			integer :: nOutputSnps
			integer :: n0, n1, n2
			integer,allocatable,dimension(:):: WorkTmp
			character(len=300) :: TmpId
			real(kind=real64) :: ImputationQuality(ped%pedigreeSize-ped%nDummys,6)


			inputParams => defaultInput
			! if we are doing output

			if (inputParams%outopt==0) then


				nOutputSnps = nOutputSnps
				TmpGenos = ped%getGenotypesAsArrayWitHMissing()
				tmpPhase = ped%getPhaseAsArrayWithMissing()
			else if (inputParams%outopt==1) then
				nOutputSnps = inputParams%nsnpRaw
				allocate(WorkTmp(nOutputSnps))
				allocate(TmpGenos(ped%pedigreeSize,nOutputSnps))
				allocate(TmpPhase(ped%pedigreeSize,nOutputSnps,2))

				tmpGenos = MISSINGGENOTYPECODE
				tmpPhase = MISSINGPHASECODE

				if (inputParams%hmmoption==RUN_HMM_NGS) then
					do i=1,inputParams%nsnpRaw
						SnpIncluded(i)=i
					enddo
				endif

				l=0
				do j=1,nOutputSnps
					if (SnpIncluded(j)/=0) then
						l=l+1
						TmpGenos(:,j)=ped%getAllGenotypesAtPositionWithUngenotypedAnimals(l)
						TmpPhase(:,j,1)=ped%getPhaseAtPositionUngenotypedAnimals(l,1)
						TmpPhase(:,j,2)=ped%getPhaseAtPositionUngenotypedAnimals(l,2)
					endif
				enddo
				block
					integer :: tmpIDInt,f
					if (inputParams%inteditstat == 1) then
						open (unit=42,file=trim(inputParams%GenotypeFile),status='old')
						do
							read (42,*, iostat=f) TmpId,WorkTmp(:)

							if (f /= 0) then
								exit
							endif

							tmpIDInt = ped%dictionary%getValue(trim(tmpID))

							if (tmpIDInt /= DICT_NULL) then
								do j=1,inputParams%nSnpRaw
									if (SnpIncluded(j)==0) then
										if (WorkTmp(j)==1) TmpGenos(tmpIDInt,j)=1
										if (WorkTmp(j)==0) then
											TmpPhase(tmpIDInt,j,:)=0
											TmpGenos(tmpIDInt,j)=0
										endif
										if (WorkTmp(j)==2) then
											TmpPhase(tmpIDInt,j,:)=1
											TmpGenos(tmpIDInt,j)=2
										endif
									endif
								enddo
							endif

						enddo

						close(42)
					end if
				end block
				deallocate(WorkTmp)
				call ped%setGenotypeFromArray(tmpGenos)
				call ped%setPhaseFromArray(TmpPhase)
				call IterateInsteadOfGeneProbs
			endif


			if (.not. inputParams%modelRecomb) then
				if (inputParams%outputOnlyGenotypedAnimals) then
					call ped%writeOutGenotypes(trim(inputparams%resultFolderPath)// DASH // "ImputeGenotypes.txt")
					call ped%writeOutPhase(trim(inputparams%resultFolderPath)// DASH // "ImputePhase.txt")
				else
					call ped%writeOutGenotypesAll(trim(inputparams%resultFolderPath)// DASH // "ImputeGenotypes.txt")
					call ped%writeOutPhaseAll(trim(inputparams%resultFolderPath)// DASH // "ImputePhase.txt")
				endif
			endif





			if (inputParams%hmmoption /= RUN_HMM_NO) then
				if (.not. allocated(ImputeGenosHMM)) then
					allocate(ImputeGenosHMM(0:ped%pedigreeSize,nOutputSnps))
				end if
				if (.not. allocated(ImputePhaseHMM)) then
					allocate(ImputePhaseHMM(0:ped%pedigreeSize,nOutputSnps,2))
				end if
				if (.not. allocated(ProbImputeGenos)) then
					allocate(ProbImputeGenos(0:ped%pedigreeSize,nOutputSnps))
				end if
				if (.not. allocated(ProbImputePhase)) then
					allocate(ProbImputePhase(0:ped%pedigreeSize,nOutputSnps,2))
				end if
				if (.not. allocated(maf)) then
					allocate(Maf(nOutputSnps))
				end if

				probImputeGenos(1:ped%pedigreeSize-ped%nDummys,:) = 9.0
				ProbImputePhase(1:ped%pedigreeSize-ped%nDummys,:,:) = 9.0

				! Feed Impute and Phase probabilites
				l=0
				do j=1,nOutputSnps
					l=l+1
					do i=1,ped%nGenotyped
						ProbImputeGenos(ped%genotypeMap(i),j)   = ProbImputeGenosHmm(i,j)
						ProbImputePhase(ped%genotypeMap(i),j,1) = ProbImputePhaseHmm(i,j,1)
						ProbImputePhase(ped%genotypeMap(i),j,2) = ProbImputePhaseHmm(i,j,2)
					enddo
				enddo

				! Impute the most likely genotypes. (Most frequent genotype)
				do i=1,ped%nGenotyped
					do j=1,nOutputSnps
						n2 = GenosCounts(i,j,2)                           ! Homozygous: 2 case
						n1 = GenosCounts(i,j,1)                           ! Heterozygous
						n0 = (inputParams%nRoundsHmm-inputParams%HmmBurnInRound) - n1 - n2        ! Homozygous: 0 case
						if ((n0>n1).and.(n0>n2)) then
							ImputeGenosHMM(ped%genotypeMap(i),j)   = 0
							ImputePhaseHMM(ped%genotypeMap(i),j,:) = 0
						elseif (n1>n2) then
							! call ped%pedigree(ped%genotypeMap(i))%individualGenotype%setGenotype(j,1)
							ImputeGenosHMM(ped%genotypeMap(i),j)   = 1
							if (ProbImputePhaseHmm(ped%genotypeMap(i),j,1) > ProbImputePhaseHmm(ped%genotypeMap(i),j,2) ) then
								ImputePhaseHMM(ped%genotypeMap(i),j,1) = 1
								ImputePhaseHMM(ped%genotypeMap(i),j,2) = 0
							else
								ImputePhaseHMM(ped%genotypeMap(i),j,1) = 0
								ImputePhaseHMM(ped%genotypeMap(i),j,2) = 1
							endif
						else
							ImputeGenosHMM(ped%genotypeMap(i),j)   = 2
							ImputePhaseHMM(ped%genotypeMap(i),j,:) = 1
						endif
					enddo
				enddo

				BLOCK
					integer :: hmmID

					open (unit=53,file="." // DASH// trim(inputparams%resultFolderPath) // DASH // "ImputePhaseHMM.txt",status="unknown")
					open (unit=54,file="." // DASH// trim(inputparams%resultFolderPath) // DASH // "ImputeGenotypesHMM.txt",status="unknown")
					open (unit=40,file="." // DASH// trim(inputparams%resultFolderPath) // DASH // "ImputePhaseProbabilities.txt",status="unknown")
					open (unit=41,file="." // DASH// trim(inputparams%resultFolderPath)// DASH // "ImputeGenotypeProbabilities.txt",status="unknown")
					do i=1,ped%nGenotyped
						hmmID = ped%genotypeMap(i)
						write (53,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputePhaseHMM(hmmID,:,1)
						write (53,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputePhaseHMM(hmmID,:,2)
						write (54,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(hmmID)%originalID,ImputeGenosHMM(hmmID,:)

						if (.not. inputParams%modelRecomb) then
							write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputePhase(hmmID,:,1)
							write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputePhase(hmmID,:,2)
							write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(hmmID)%originalID,ProbImputeGenos(hmmID,:)
						endif
					enddo

					close(53)
					close(54)
					close(40)
					close(41)
				END BLOCK
			else
				if (.not. inputParams%modelRecomb) then

					open (unit=40,file="." // DASH// trim(inputparams%resultFolderPath) // DASH // "ImputePhaseProbabilities.txt",status="unknown")
					open (unit=41,file="." // DASH// trim(inputparams%resultFolderPath)// DASH // "ImputeGenotypeProbabilities.txt",status="unknown")
					do i=1, ped%pedigreeSize-ped%nDummys
						if (ped%pedigree(i)%isDummy) then
							exit
						endif
						write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(ped%inputmap(i))%originalID,ProbImputePhase(ped%inputmap(i),:,1)
						write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(ped%inputmap(i))%originalID,ProbImputePhase(ped%inputmap(i),:,2)
						write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(ped%inputmap(i))%originalID,ProbImputeGenos(ped%inputmap(i),:)
					enddo
					close(40)
					close(41)
				endif
			endif

			allocate(Maf(nOutputSnps))
			if (inputParams%SexOpt==1) then
				do j=1,nOutputSnps
					Maf(j)=sum(ProbImputeGenos(:,j))/(2*ped%pedigreeSize-ped%nDummys)
				enddo
				open(unit=111,file="." // DASH // "Miscellaneous" // DASH // "MinorAlleleFrequency.txt", status="unknown")

				do j=1,nOutputSnps
					write (111,*) j,Maf(j)
				enddo
				close(111)
			endif

			ImputationQuality(:,1)=sum(2*Maf(:))/nOutputSnps

			ImputationQuality(:,2)=0.0

			do i=1, ped%pedigreeSize-ped%nDummys
				if (ped%pedigree(i)%isDummy) then
					exit
				endif
				do j=1,nOutputSnps
					!ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-(2*Maf(j)))
					ImputationQuality(i,2)=ImputationQuality(i,2)+abs(ProbImputeGenos(i,j)-((Maf(j)**4)+(4*(Maf(j)**2)*((1.0-Maf(j))**2))+(1.0-Maf(j)**4)))
				enddo
				ImputationQuality(i,2)=ImputationQuality(i,2)/nOutputSnps
				ImputationQuality(i,3)=float(nOutputSnps-ped%pedigree(i)%IndividualPhase(1)%numberMissing())/nOutputSnps
				ImputationQuality(i,4)=float(nOutputSnps-ped%pedigree(i)%IndividualPhase(2)%numberMissing())/nOutputSnps
				ImputationQuality(i,5)=(ImputationQuality(i,3)+ImputationQuality(i,4))/2

				ImputationQuality(i,6)=float(nOutputSnps-ped%pedigree(i)%IndividualGenotype%numMissing())/nOutputSnps
				write (50,'(a20,20000f7.2)') ped%pedigree(i)%originalID,ImputationQuality(i,:)
			enddo

			open (unit=51,file="." // DASH// trim(inputparams%resultFolderPath) // DASH // "ImputationQualitySnp.txt",status="unknown")
			do j=1,nOutputSnps
				write (51,'(i10,20000f7.2)') j,float(((ped%pedigreeSize-(ped%nDummys+1))+1)-ped%countMissingGenotypesNoDummys())/((ped%pedigreeSize-(ped%nDummys+1))+1)
			enddo
			close(51)
			inputParams%WellPhasedThresh=inputParams%WellPhasedThresh/100

			open (unit=52,file="." // DASH// trim(inputparams%resultFolderPath) // DASH // "WellPhasedIndividuals.txt",status="unknown")
			do i=1, ped%pedigreeSize-ped%nDummys
				if (ImputationQuality(i,5)>=inputParams%WellPhasedThresh) then
					write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(ped%inputmap(i))%originalID,tmpPhase(ped%inputmap(i),:,1)
					write (52,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)') ped%pedigree(ped%inputmap(i))%originalID,tmpPhase(ped%inputmap(i),:,2)
				endif
			enddo
			close(52)



		end subroutine WriteOutResults


		!#############################################################################################################################################################################################################################

		subroutine ModelRecomb
			use Global
			use imputation, only : InsteadOfReReadGeneProb
			use AlphaImputeSpecFileModule
			use AlphaImputeInputOutputModule
			implicit none

			integer :: e,i,j,k,l,SuperJ,StartDisFound,EndDisFound,HetEnd,HetStart,RSide,LSide,PatMat,SireDamRL,nSnpFinal,Counter
			integer :: StartDisPrev,EndDisPrev,RecombOnOff
			integer :: GamA,GamB,StartDisOld,StartDisTmp
			integer :: CountRightSwitch,CountLeftSwitch,PedId,StartDis,EndDis,StartJ,nRec
			integer(kind=1),allocatable,dimension(:,:,:) :: WorkPhase,TempWork
			integer,allocatable,dimension(:) :: WorkLeft,WorkRight,TempVec,StR,EnR,StRNarrow,EnRNarrow
			real,allocatable,dimension(:) :: LengthVec
			real,allocatable,dimension(:,:) :: PatAlleleProb,MatAlleleProb,GeneProbWork
			character(len=7) :: cm
			type(AlphaImputeInput), pointer :: inputParams
			integer(kind=1) :: phase(2),tmpPhase
			integer :: tmp
			inputParams => defaultInput


			write(cm,'(I7)') inputParams%nSnpRaw !for formatting
			cm = adjustl(cm)

			open (unit=42, file=trim(inputparams%resultFolderPath) // DASH // "RecombinationInformation.txt")
			open (unit=43, file=trim(inputparams%resultFolderPath) // DASH // "RecombinationInformationNarrow.txt")
			open (unit=44, file=trim(inputparams%resultFolderPath) // DASH // "NumberRecombinations.txt")
			open (unit=45, file=trim(inputparams%resultFolderPath) // DASH // "RecombinationInformationR.txt")
			open (unit=46, file=trim(inputparams%resultFolderPath) // DASH // "RecombinationInformationNarrowR.txt")


			counter = 0

			! Check whether to consider all the raw snps or only the snps left after the edition procedure
			! If EditedSnpOut in Spec file
			if (inputParams%outopt==0) then
				nSnpFinal=inputParams%nsnp
				! If AllSnpOut in Spec file
			else
				nSnpFinal=inputParams%nSnpRaw
			endif
			if (allocated(GlobalWorkPhase)) then
				deallocate(GlobalWorkPhase)
			endif
			allocate(GlobalWorkPhase(0:ped%pedigreeSize,nSnpFinal,2))
			allocate(WorkPhase(0:ped%pedigreeSize,nSnpFinal,2))
			allocate(TempVec(nSnpFinal))
			allocate(LengthVec(nSnpFinal))
			allocate(WorkLeft(nSnpFinal))
			allocate(WorkRight(nSnpFinal))
			allocate(TempWork(0:ped%pedigreeSize,nSnpFinal,2))
			allocate(PatAlleleProb(nSnpFinal,2))
			allocate(MatAlleleProb(nSnpFinal,2))
			allocate(GeneProbWork(nSnpFinal,4))
			allocate(StR(nSnpFinal))
			allocate(EnR(nSnpFinal))
			allocate(StRNarrow(nSnpFinal))
			allocate(EnRNarrow(nSnpFinal))

			WorkPhase=ped%getPhaseAsArray()
			call InsteadOfReReadGeneProb

			l=0
			do j=1,nSnpFinal

				if (SnpIncluded(j)/=0) then
					l=l+1
					WorkPhase(:,j,1)=GlobalWorkPhase(:,l,1)
					WorkPhase(:,j,2)=GlobalWorkPhase(:,l,2)
				endif
			enddo

			do i=1,ped%pedigreeSize-ped%nDummys
				HetEnd=-1
				HetStart=-1
				WorkRight(:)=9
				WorkLeft(:)=9
				do e=1,2
					nRec=0
					StR=-9
					EnR=-9
					StRNarrow=-9
					EnRNarrow=-9

					PatMat=e
					SireDamRL=e+1
					CountLeftSwitch=0
					CountRightSwitch=0

					pedID=ped%pedigree(i)%getSireDamNewIDByIndexNoDummy(e+1)
					! means parent is a dummy or no parent
					if (pedId == 0) then
						cycle
					endif
					if ((inputParams%SexOpt==1).and.(ped%pedigree(PedId)%gender==inputParams%hetGameticStatus)) cycle
					if (((float(ped%pedigree(i)%individualPhase(1)%numberMissing()+ped%pedigree(i)%individualPhase(2)%numberMissing())/(2*nSnpFinal))<0.30)) then          !(RecIdHDIndex(PedId)==1)
						WorkRight=9
						RSide=9
						do j=1,nSnpFinal
							phase(1) = ped%pedigree(pedId)%individualPhase(1)%getPhase(j)
							phase(2)  = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
							if ((phase(1)/=phase(2)).and.(phase(1)/=9).and.(phase(2)/=9))  then
								HetStart=j
								tmpPhase = ped%pedigree(i)%individualPhase(PatMat)%getPhase(HetStart)
								if (tmpPhase==phase(1)) then
									WorkRight(HetStart)=1
									RSide=1
									exit
								endif
								if (tmpPhase==phase(2)) then
									WorkRight(HetStart)=2
									RSide=2
									exit
								endif
							endif
						enddo
						if (RSide/=9) then
							do j=HetStart+1,nSnpFinal
								phase(1) = ped%pedigree(pedId)%individualPhase(RSide)%getPhase(j)
								tmpPhase = ped%pedigree(i)%individualPhase(patMat)%getPhase(j)
								if ((tmpPhase /= phase(1)).and.(phase(1)/=9).and.(tmpPhase/=9)) then
									RSide=abs((RSide-1)-1)+1
									CountRightSwitch=CountRightSwitch+1
								endif
								WorkRight(j)=RSide
							enddo
						endif

						LSide=9
						do j=nSnpFinal,1,-1
							phase(1) = ped%pedigree(pedId)%individualPhase(1)%getPhase(j)
							phase(2) = ped%pedigree(pedId)%individualPhase(2)%getPhase(j)
							if ((phase(1)/=phase(2)).and.(phase(1)/=9).and.(phase(2)/=9))  then
								HetEnd=j
								tmpPhase = ped%pedigree(i)%individualPhase(patMat)%getPhase(HetEnd)
								if (tmpPhase==phase(1)) then
									WorkLeft(HetEnd)=1  !£$$$$
									LSide=1
									exit
								endif
								if (tmpPhase==phase(2)) then
									WorkLeft(HetEnd)=2  !£$$$$
									LSide=2
									exit
								endif
							endif
						enddo
						if (LSide/=9) then
							do j=HetEnd-1,1,-1

								phase(1) = ped%pedigree(i)%individualPhase(patMat)%getPhase(j)
								phase(2) = ped%pedigree(pedId)%individualPhase(LSide)%getPhase(j)
								if (phase(1)/=phase(2) .and.(phase(2)/=9).and.(phase(1)/=9)) then
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
						StartJ=1

						!prototype start

						StartDis=1
						StartDisOld=1
						EndDis=nSnpFinal
						nRec=0
						SuperJ=0
						TempVec=9
						LengthVec=0.0

						StartDisPrev=StartDis
						EndDisPrev=-9

						do while (SuperJ<inputParams%nsnp)
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
								do k=StartDis,1,-1
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

								tmp =(EndDis-StartDis)

								if (tmp==0 ) then
									LengthVec(StartDis:EndDis) = 0
								else
									LengthVec(StartDis:EndDis)=1.0/tmp
								endif
								StartDisPrev=StartDis
								EndDisPrev=EndDis

							endif
						enddo
						!prototype end (temp)
						StR(nRec+1)=-9
						EnR(nRec+1)=-9

						GamA=9
						GamB=9
						RecombOnOff=0
						k=1
						do j=1,nSnpFinal

							if (j==StR(k)) then
								k=k+1
								RecombOnOff=1
								Counter=0
								GamA=WorkRight(j)

								If (GamA==1) GamB=2
								If (GamA==2) GamB=1
							endif
							if (k>1) then
								if (j==EnR(k-1)) then
									RecombOnOff=0
								endif
							endif
							if (GamA==9) cycle
							if (RecombOnOff==1) then
								Counter=Counter+1
								if (Counter/=1) then
									block

										integer(kind=1) :: phase1, phase2

										phase1 = ped%pedigree(PedId)%individualPhase(1)%getPhase(j)
										phase2 = ped%pedigree(PedId)%individualPhase(2)%getPhase(j)
										if (phase1/= phase2) then
											if ((phase1/=9).and.(phase2/=9)) then

												block

													integer :: one, two
													call ped%pedigree(i)%individualPhase(e)%setPhase(j,9)
													call ped%pedigree(i)%individualGenotype%setGenotype(j,9)
													one = ((1.0-(LengthVec(j)*Counter))*ped%pedigree(pedid)%individualPhase(gamA)%getPhase(j))
													two = (LengthVec(j)*Counter*ped%pedigree(pedid)%individualPhase(gamB)%getPhase(j))


													ProbImputePhase(i,j,e) = one + two
													ProbImputeGenos(i,j)=ProbImputePhase(i,j,1)+ProbImputePhase(i,j,2)

												end block
											endif
										endif
									end block
								endif
							endif
						enddo
					endif
					write (44,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec
					write (42,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec,StR(1:nRec)
					write (42,'(a20,20000i20)') ped%pedigree(i)%originalID,nRec,EnR(1:nRec)
					do j=1,nRec
						write (45,'(a20,20000i20)') ped%pedigree(i)%originalID,e,StR(j),EnR(j)
					enddo

				enddo


			enddo



			open (unit=40,file=trim(inputparams%resultFolderPath) // DASH // "ImputePhaseProbabilities.txt",status="unknown")
			open (unit=41,file=trim(inputparams%resultFolderPath) // DASH // "ImputeGenotypeProbabilities.txt",status="unknown")


			call ped%writeOutGenotypesAll(trim(inputparams%resultFolderPath) // DASH // "ModelRecomb.txt")

			if (inputParams%outputOnlyGenotypedAnimals) then
				call ped%writeOutGenotypes(trim(inputparams%resultFolderPath)// DASH // "ImputeGenotypes.txt")
				call ped%writeOutPhase(trim(inputparams%resultFolderPath)// DASH // "ImputePhase.txt")
			else
				call ped%writeOutGenotypesAll(trim(inputparams%resultFolderPath)// DASH // "ImputeGenotypes.txt")
				call ped%writeOutPhaseAll(trim(inputparams%resultFolderPath)// DASH // "ImputePhase.txt")
			endif

			block
				integer :: tmpID, limit

				if (inputParams%outputonlygenotypedanimals) then
					limit = ped%nGenotyped
				else
					limit = ped%pedigreeSize - ped%nDummys
				endif

				do i=1, limit
					if (ped%pedigree(i)%isDummy) then
						exit
					endif
					if (inputParams%outputonlygenotypedanimals) then
						tmpId = ped%genotypeMap(i)
					else
						tmpId = ped%inputMap(i)
					endif
					write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(tmpId)%originalID,ProbImputePhase(tmpId,:,1)
					write (40,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(tmpId)%originalID,ProbImputePhase(tmpId,:,2)
					write (41,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') ped%pedigree(tmpId)%originalID,ProbImputeGenos(tmpId,:)
				enddo
			end block
			print*, " "
			print*, " ","Imputation by detection of recombination events completed"

		end subroutine ModelRecomb


		!#############################################################################################################################################################################################################################

		subroutine IterateParentHomoFill
			use Global
			use AlphaImputeSpecFileModule
			implicit none

			integer :: e,i,PedLoc, id
			type(AlphaImputeInput), pointer :: inputParams
			type(genotype) :: tmpGeno
			inputParams => defaultInput

			if (inputParams%SexOpt==0) then
				do i=1,ped%pedigreeSize-ped%nDummys
					do e=1,2
						PedLoc=e+1
						id = ped%pedigree(i)%getSireDamNewIDByIndex(pedLoc)
						if (id == 0) then
							cycle
						endif
						call tmpGeno%newGenotypeHap(ped%pedigree(id)%IndividualPhase(1),ped%pedigree(id)%individualPhase(2))
						call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(e))
					enddo
				enddo
			else
				do i=1,ped%pedigreeSize-ped%nDummys
					if (ped%pedigree(i)%gender==inputParams%HomGameticStatus) then
						do e=1,2
							id = ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
							if (id == 0) then
								cycle
							endif
							call tmpGeno%newGenotypeHap(ped%pedigree(id)%IndividualPhase(1),ped%pedigree(id)%individualPhase(2))
							call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(e))
						enddo
					else
						id = ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1)
						if (id == 0) then
							cycle
						endif
						call tmpGeno%newGenotypeHap(ped%pedigree(id)%IndividualPhase(1),ped%pedigree(id)%individualPhase(2))
						call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(1))
						call tmpGeno%setHaplotypeFromGenotypeIfMissing(ped%pedigree(i)%individualPhase(2))
					endif
				enddo
			endif


		end subroutine IterateParentHomoFill


		!######################################################################################################################################################################################

		subroutine InsteadOfGeneProb
			! Phase haplotypes whenever there is enough information from the parents (homozygous case)
			use Global
			use AlphaImputeSpecFileModule
			implicit none

			integer :: e,i,j,ParId
			type(AlphaImputeInput), pointer :: inputParams

			inputParams => defaultInput
			if (inputParams%SexOpt==1) then                                                     ! Sex chromosome

				GlobalWorkPhase=9
				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,inputParams%nsnp                                                     ! Phase in the homozygous case

						if (ped%pedigree(i)%individualGenotype%getGenotype(j)==0) then
							GlobalWorkPhase(i,j,:)=0
						else if (ped%pedigree(i)%individualGenotype%getGenotype(j)==2) then
							GlobalWorkPhase(i,j,:)=1
						endif
					enddo
					if (ped%pedigree(i)%gender/=inputParams%hetGameticStatus) then                        ! Am I homogametic?
						do e=1,2                                                    ! Phase a single haplotype whenever my parents are homozygous
							ParId=ped%pedigree(i)%getSireDamNewIDByIndex(e+1)
							do j=1,inputParams%nsnp
								if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==0) then
									GlobalWorkPhase(i,j,e)=0

								else if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==2) then
									GlobalWorkPhase(i,j,e)=1
								endif
							enddo
						enddo
					else                                                            ! Am I heterogametic?
						ParId=ped%pedigree(i)%getSireDamNewIDByIndex(inputParams%HomGameticStatus+1)
						do j=1,inputParams%nsnp                                                 ! Phase the two haplotypes whenever my homogametic parent is homozygous
							if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==0) then
								GlobalWorkPhase(i,j,:)=0
							else if (ped%pedigree(parId)%individualGenotype%getGenotype(j)==2) then
								GlobalWorkPhase(i,j,:)=1
							endif
						enddo
					endif
				enddo

				GlobalWorkPhase(0,:,:)=9

				! todo  -feel this can be optimised- do we need to do this?
				call ped%setPhaseFromArray(GlobalWorkPhase)

				allocate(GlobalTmpCountInf(ped%pedigreeSize-ped%nDummys,6))
				GlobalTmpCountInf(:,:)=0

			else                                                                    ! Other chromosome

				GlobalWorkPhase=9
				do i=1,ped%pedigreeSize-ped%nDummys
					do j=1,inputParams%nsnp                                                     ! Phase in the homozygous case
						if (ped%pedigree(i)%individualGenotype%getGenotype(j)==0) then
							GlobalWorkPhase(i,j,:)=0
						else if (ped%pedigree(i)%individualGenotype%getGenotype(j)==2) then
							GlobalWorkPhase(i,j,:)=1
						endif
					enddo
					if (associated(ped%pedigree(i)%sirePointer)) then
						if (ped%pedigree(i)%sirePointer%isDummy) exit
						ParId=ped%pedigree(i)%sirePointer%id

						do j=1,inputParams%nsnp                                                     ! Phase if my father is homozygous
							if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==0) then
								GlobalWorkPhase(i,j,1)=0
							else if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==2) then
								GlobalWorkPhase(i,j,1)=1
							endif
						enddo
					endif
					if (associated(ped%pedigree(i)%damPointer)) then
						if (ped%pedigree(i)%damPointer%isDummy) exit
						ParId=ped%pedigree(i)%damPointer%id
						do j=1,inputParams%nsnp                                                     ! Phase if my mother is homozygous
							if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==0) then
								GlobalWorkPhase(i,j,2)=0
							else if (ped%pedigree(ParId)%individualGenotype%getGenotype(j)==2) then
								GlobalWorkPhase(i,j,2)=1
							endif
						enddo
					endif

				enddo

				GlobalWorkPhase(0,:,:)=9

				call ped%setPhaseFromArray(GlobalWorkPhase)

			endif

		end subroutine InsteadOfGeneProb


		!#############################################################################################################################################################################################################################
		subroutine ClassifyAnimByChips
			! Classify animals according to the HD chip information with a margin of missing markers
			! The condition for an animal to be classify with a particular HD snp chip is:
			! If the missing markers are not above a threshold and the number of markers is below the
			! nominal number of markers for that chip
			! LD animals or HD animals with missing markers above a threshold are not assign to any
			! snp panel

			use Global

			use AlphaImputeSpecFileModule
			use ISO_Fortran_Env
			implicit none

			integer :: i, j, CountMiss, UOutputs
			logical, allocatable :: printed(:)
			type(AlphaImputeInput), pointer :: inputParams

			inputParams => defaultInput
			open(newunit=UOutputs, file="." // DASH // "Miscellaneous" // DASH // "SnpCallRateByAnimalByChip.txt",status='unknown')
			allocate(animChip(ped%pedigreeSize-ped%nDummys))
			animChip(:)=0

			allocate(printed(ped%pedigreeSize-ped%nDummys))
			printed=.FALSE.

			do i=1,ped%pedigreeSize-ped%nDummys
				CountMiss=ped%pedigree(i)%individualGenotype%numMissing()
				do j=1,inputParams%MultiHD
					if ( (CountMiss-(inputParams%nsnp-inputParams%nSnpByChip(j))) < (1.0-inputParams%PercGenoForHD)*inputParams%nSnpByChip(j)&
					.and. (inputParams%nsnp-CountMiss)<inputParams%nSnpByChip(j)&
					.and. ped%pedigree(i)%genotyped) then
					animChip(i)=j
					write(UOutputs,'(a20,6f5.1)') ped%pedigree(i)%originalID, (inputParams%nsnp-CountMiss)*100/real(inputParams%nSnpByChip(j))
					exit
				endif
				if ((CountMiss-(inputParams%nsnp-inputParams%nSnpByChip(j))) > (1.0-inputParams%PercGenoForHD)*inputParams%nSnpByChip(j)&
				! .and. animChip(i)/=0&
				.and. printed(i)==.false.&
				.and. ped%pedigree(i)%genotyped) Then
				write(UOutputs,'(a20,6f5.1)') ped%pedigree(i)%originalID, (inputParams%nsnp-CountMiss)*100/real(inputParams%nSnpByChip(j))
				printed(i)=.true.
			end if
		enddo
	enddo
	close(UOutputs)
end subroutine ClassifyAnimByChips

!#############################################################################################################################################################################################################################
subroutine InternalEdit
	use Global
	use AlphaImputeSpecFileModule
	implicit none

	integer :: i,j,k,CountMiss,CountHD,nSnpR,dum, phaseUnit
	real, allocatable, dimension(:) :: SnpSummary, TempFreq,Counter
	character (len=300) :: dumC,FileName
	type(AlphaImputeInput), pointer :: inputParams

	inputParams => defaultInput

	allocate(SnpSummary(inputParams%nsnp))
	allocate(TempFreq(inputParams%nsnp))
	allocate(Counter(inputParams%nsnp))
	allocate(SnpIncluded(inputParams%nsnp))
	allocate(Setter(0:ped%pedigreeSize-ped%nDummys))

	SnpIncluded(:)=0
	if ((inputParams%managephaseon1off0==0).and.(inputParams%NoPhasing==1)) then
		FileName = trim(inputParams%phasePath) // DASH // "EditingSnpSummary.txt"
		open(newunit=phaseUnit,file=trim(FileName),status="old")
		do i=1,inputParams%nsnp
			read (phaseUnit,*) dum,SnpSummary(i),SnpIncluded(i)
		enddo
		close(phaseUnit)
	endif

	if (inputParams%outopt==1 .or. inputParams%NoPhasing==1) then
		do i=1,inputParams%nsnp
			SnpIncluded(i)=i
		enddo

	endif
	! I user do not specify any file with HD individuals
	if (inputParams%UserDefinedHD==0) then
		Setter(0)=0
		Setter(1:ped%pedigreeSize-ped%nDummys)=1
		do i=1,ped%pedigreeSize-ped%nDummys
			CountMiss=ped%pedigree(i)%individualGenotype%numMissing()
			if (inputParams%MultiHD/=0) then
				! Disregard animals at LD or those HD animals with a number of markers missing
				if (animChip(i)==0) then
					Setter(i)=0
				endif
			else
				if ((float(CountMiss)/inputParams%nsnp)>(1.0-inputParams%PercGenoForHD)) then
					Setter(i)=0
				endif
			endif
		enddo
		CountHD=count(Setter(:)==1)
	else                                ! User has specified HD individuals
		Setter(0)=0
		Setter(1:ped%pedigreeSize-ped%nDummys)=0

		CountHD=0
		block
			integer :: tmpID
			do
				read (inputParams%AnimalFileUnit,*, iostat=k) dumC

				CountHD=CountHD+1
				if (k/=0) then
					CountHD=CountHD-1
					exit
				endif

				tmpID = ped%dictionary%getValue(dumC)
				if (tmpID /= DICT_NULL) then
					Setter(tmpID)=1
					exit
				endif
			enddo
		end block
		CountHD=count(Setter(:)==1)
		print*, " "
		print*, " ",CountHD," valid indiviudals in the user specified AlphaPhase file"
	endif

	if (inputParams%managephaseon1off0==1) then
		TempFreq(:)=0.0
		Counter(:)=0
		SnpSummary=0.0
		do j=1,inputParams%nsnp
			do i=1, ped%pedigreeSize-ped%nDummys
				if (.not. ped%pedigree(i)%individualGenotype%isMissing(j)) then
					TempFreq(j)=TempFreq(j)+float(ped%pedigree(i)%individualGenotype%getGenotype(j))
					Counter(j)=Counter(j)+2
				else
					if (Setter(i) == 1) then
						SnpSummary(j)=SnpSummary(j)+1.0
					endif
				end if
			end do
			if (Counter(j)>0.000000) TempFreq(j)=TempFreq(j)/Counter(j)
		end do
	endif

	if (inputParams%MultiHD/=0 .or. inputParams%IntEditStat==0) then
		nSnpR=inputParams%nsnp
		if (inputParams%managephaseon1off0==1) then
			do i=1,inputParams%nsnp
				SnpIncluded(i)=i
			enddo
		endif
	else
		if (inputParams%managephaseon1off0==1) then
			SnpSummary(:)=SnpSummary(:)/CountHD
			nSnpR=0
			do j=1,inputParams%nsnp
				if ((SnpSummary(j)<inputParams%PercSnpMiss).and.((TempFreq(j)>0.00000001).and.(TempFreq(j)<0.9999999))) nSnpR=nSnpR+1
			enddo

		else
			nSnpR=count(SnpIncluded(:)/=0)
		endif


		if (nSnpR==inputParams%nsnp) then
			do i=1,inputParams%nsnp
				SnpIncluded(i)=i
			enddo
		else
			if (inputParams%managephaseon1off0==1) then


				! Remove snps from individuals
				block
					integer(kind=1),dimension(:),allocatable :: old, temp ,tempphase1,tempphase2, oldphase1,oldphase2

					allocate(temp(nSNpR))
					allocate(tempphase1(nSNpR))
					allocate(tempphase2(nSNpR))
					do i=1, ped%pedigreeSize
						k=0
						old = ped%pedigree(i)%individualGenotype%toIntegerArray()
						oldphase1 = ped%pedigree(i)%individualPhase(1)%toIntegerArray()
						oldphase2 = ped%pedigree(i)%individualPhase(2)%toIntegerArray()
						do j=1,inputParams%nsnp
							if ((SnpSummary(j)<inputParams%PercSnpMiss).and.((TempFreq(j)>0.00000001).and.(TempFreq(j)<0.9999999))) then
								k=k+1
								SnpIncluded(j)=k
								temp(k) = old(j)
								tempphase1(k) =oldphase1(j)
								tempphase2(k) =oldphase2(j)
							else
								SnpIncluded(j)=0
							endif

						enddo
						call ped%pedigree(i)%individualGenotype%newGenotypeInt(temp)
						call ped%pedigree(i)%individualphase(1)%newhaplotypeInt(tempphase1)
						call ped%pedigree(i)%individualphase(2)%newhaplotypeInt(tempphase2)
					enddo
					deallocate(temp)
					deallocate(tempphase1)
					deallocate(tempphase2)

				end block
				inputParams%nsnp=nSnpR
			endif
		endif

		if (inputParams%UserDefinedHD==0) then
			Setter = 0
			do i=1,ped%nGenotyped
				setter(ped%genotypeMap(i)) =1
				CountMiss=ped%pedigree(ped%genotypeMap(i))%individualGenotype%numMissing()
				if ((float(CountMiss)/inputParams%nsnp)>(1.0-inputParams%SecondPercGenoForHD)) then
					Setter(ped%genotypeMap(i))=0
				endif
			enddo
			CountHD=count(Setter(:)==1)
		else
			do i=1,ped%nGenotyped
				if (Setter(ped%genotypeMap(i))==1) then
					CountMiss=ped%pedigree(ped%genotypeMap(i))%individualGenotype%numMissing()
					if ((float(CountMiss)/inputParams%nsnp)>(1.0-inputParams%SecondPercGenoForHD)) then
						Setter(ped%genotypeMap(i))=0
					endif
				endif
			enddo
			CountHD=count(Setter(:)==1)
		endif
	endif

	open (unit=102,file="." // DASH // "Miscellaneous" // DASH // "EditingSnpSummary.txt",status="unknown")

	if ((inputParams%ManagePhaseOn1Off0==0).and.(inputParams%NoPhasing==1)) then
		FileName = trim(inputParams%PhasePath) // DASH // "EditingSnpSummary.txt"
	else
		FileName = "." // DASH // "Phasing" // DASH // "EditingSnpSummary.txt"
	end if
	open (unit=112,file=FileName,status="unknown")


	do j=1,inputParams%nSnpRaw
		write (102,*) j,SnpSummary(j),SnpIncluded(j)        !'(i,1x,f5.3,1x,i)'
		write (112,*) j,SnpSummary(j),SnpIncluded(j)        !'(i,1x,f5.3,1x,i)'
	enddo
	close(102)
	close(112)
	print*, " "
	print*, " "
	print *,ped%nGenotyped
	print*, " ",CountHD," indiviudals passed to AlphaPhase"
	print*, " ",inputParams%nsnp," snp remain after editing"

	! TODO clean this whole subroutine


	! we add animals to hd list here
	do i=1, ped%nGenotyped
		if (setter(ped%GenotypeMap(i)) == 1) then
			call ped%setAnimalAsHD(ped%GenotypeMap(i))
		endif
	enddo
	deallocate(SnpSummary)
	deallocate(TempFreq)
	deallocate(Counter)
end subroutine InternalEdit

!#############################################################################################################################################################################################################################

subroutine FillInBasedOnOffspring
	! Genotype SNPs based on the genetic information of my offsprings
	use Global
	use AlphaImputeSpecFileModule
	implicit none

	integer :: i,j,k
	integer, allocatable, dimension(:) :: count0,count1,count2
	integer :: tmpVal
	type(AlphaImputeInput), pointer :: inputParams
	type(individual) ,pointer :: tmpOff
	inputParams => defaultInput

	allocate(Count0(inputParams%nsnp))
	allocate(Count1(inputParams%nsnp))
	allocate(Count2(inputParams%nsnp))
	do i=1,ped%pedigreeSize-ped%nDummys ! These are parents
		! This three variables will count the different number of genotypes of the offsprings
		Count0=0
		Count1=0
		Count2=0
		do j=1,ped%pedigree(i)%nOffs ! These are offsprings
			tmpOff => ped%pedigree(i)%offsprings(j)%p
			if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and.(tmpOff%gender==inputParams%hetGameticStatus)) cycle
			do k=1,inputParams%nsnp
				if (ped%pedigree(i)%individualGenotype%isMissing(k)) then ! If my parent is not genotyped
					tmpVal = tmpOff%individualGenotype%getGenotype(k)
					if (tmpVal==0) Count0(k)=Count0(k)+1    ! Number of offspring genotype as 0
					if (tmpVal==1) Count1(k)=Count1(k)+1    ! Number of offspring genotype as 1
					if (tmpVal==2) Count2(k)=Count2(k)+1    ! Number of offspring genotype as 2
				endif
			enddo

		enddo

		do k=1,inputParams%nsnp
			if ((Count0(k)+Count1(k)+Count2(k))>OffspringFillMin) then
				if (Count0(k)==0) call ped%pedigree(i)%individualGenotype%setGenotype(k,2)                       ! This is the most likely thing, but it might be not true
				if (Count2(k)==0) call ped%pedigree(i)%individualGenotype%setGenotype(k,0)                       ! This is the most likely thing, but it might be not true
				if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus)) cycle
				if ((Count0(k)>2).and.(Count2(k)>2)) call ped%pedigree(i)%individualGenotype%setGenotype(k,1)    ! This is the most likely thing, but it might be not true
			endif
		enddo
	enddo

end subroutine FillInBasedOnOffspring

!#############################################################################################################################################################################################################################

subroutine FillInSnp
	! Genotype SNPs based on the pedigree information
	use Global
	use AlphaImputeSpecFileModule
	implicit none

	integer :: i,j,k,TurnOn
	type(AlphaImputeInput), pointer :: inputParams
	integer :: tmpParentId
	type(individual), pointer :: tmpMother, tmpFather, tmpAnim
	inputParams => defaultInput

	do i=1,ped%pedigreeSize-ped%nDummys
		do k=2,3
			TurnOn=1
			tmpParentId = ped%pedigree(i)%getSireDamNewIDByIndexNoDummy(k)
			! tmpAnim => ped%pedigree(i)%getSireDamObjectByIndex(k)

			! if the proband is heterogametic, and
			! considering the heterogametic parent, then avoid!!
			if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and. ((k-1)==inputParams%hetGameticStatus)) TurnOn=0
			if (tmpParentId /= 0) then
				tmpAnim => ped%pedigree(tmpParentId)
				! if ((inputParams%SexOpt==1).and.(ped%pedigree(i)%gender==inputParams%hetGameticStatus).and. (pedigree%(i)%parent(k-1)%gender==inputParams%hetGameticStatus)) TurnOn=0
				! Homogametic individuals and the homogametic parent of a heterogametic individual
				if (TurnOn==1) then
					do j=1,inputParams%nsnp
						if ((ped%pedigree(i)%individualGenotype%getGenotype(j)==0).and.(ped%pedigree(tmpParentId)%individualGenotype%getGenotype(j)==2)) then
							call ped%pedigree(i)%individualGenotype%setGenotype(j,9)
							call tmpAnim%individualGenotype%setGenotype(j,9)
						endif
						if ((ped%pedigree(i)%individualGenotype%getGenotype(j)==2).and.(ped%pedigree(tmpParentId)%individualGenotype%getGenotype(j)==0)) then
							call ped%pedigree(i)%individualGenotype%setGenotype(j,9)
							call tmpAnim%individualGenotype%setGenotype(j,9)
						endif
					enddo
				endif
			endif
		enddo
	enddo



	! WARNING: This can be refactored
	do i=1,ped%pedigreeSize-ped%nDummys

		if(.not. ped%pedigree(i)%Founder) then
			! tmpFather =>ped%pedigree(i)%getSireDamObjectByIndex(2)
			! tmpMother =>ped%pedigree(i)%getSireDamObjectByIndex(3)

			tmpFather =>ped%pedigree(ped%pedigree(i)%getSireDamNewIDByIndex(2))
			tmpMother =>ped%pedigree(ped%pedigree(i)%getSireDamNewIDByIndex(3))

			do j=1,inputParams%nsnp


				if (ped%pedigree(i)%individualGenotype%isMissing(j)) then
					if ((tmpFather%individualGenotype%getGenotype(j)==0).and.(tmpMother%individualGenotype%getGenotype(j)==0)) then
						call ped%pedigree(i)%individualGenotype%setGenotype(j,0)
					else if ((tmpFather%individualGenotype%getGenotype(j)==2).and.(tmpMother%individualGenotype%getGenotype(j)==2)) then
						call ped%pedigree(i)%individualGenotype%setGenotype(j,2)
					endif

					if (inputParams%SexOpt==1) then
						if (ped%pedigree(i)%gender/=inputParams%hetGameticStatus) then
							if ((tmpFather%individualGenotype%getGenotype(j)==0).and.((tmpMother%individualGenotype%getGenotype(j)==2))) call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
							if ((tmpFather%individualGenotype%getGenotype(j)==2).and.((tmpMother%individualGenotype%getGenotype(j)==0))) call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
						else
							! HomGameticSatus(1 or 2) +1 = sire (2) or dam (3)
							if (inputParams%HomGameticStatus == 1) then
								if (tmpFather%individualGenotype%getGenotype(j)==0) call ped%pedigree(i)%individualGenotype%setGenotype(j,0)
								if (tmpFather%individualGenotype%getGenotype(j)==2) call ped%pedigree(i)%individualGenotype%setGenotype(j,2)
							else if (inputParams%HomGameticStatus == 2) then
								if (tmpMother%individualGenotype%getGenotype(j)==0) call ped%pedigree(i)%individualGenotype%setGenotype(j,0)
								if (tmpMother%individualGenotype%getGenotype(j)==2) call ped%pedigree(i)%individualGenotype%setGenotype(j,2)
							endif
						endif
					else
						if ((tmpFather%individualGenotype%getGenotype(j)==0).and.(tmpMother%individualGenotype%getGenotype(j)==2)) Then
							call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
						else if ((tmpFather%individualGenotype%getGenotype(j)==2).and.(tmpMother%individualGenotype%getGenotype(j)==0)) Then
							call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
						endif
					endif
				endif
			enddo
		endif
	enddo





end subroutine FillInSnp

!#############################################################################################################################################################################################################################

subroutine CheckParentage
	! Check which is the parentage relation between animals.
	! Set the vector baseline that specify whether an animal
	! has no parents or its parents have been pruned

	use AlphaImputeSpecFileModule
	use Global, only: ped, DisagreeThreshold

	implicit none

	type(AlphaImputeInput), pointer :: inputParams
	integer :: inconsistencies

	inputParams => defaultInput

	call ped%sortPedigreeAndOverwrite()
	inconsistencies = ped%findMendelianInconsistencies(DisagreeThreshold,"." // DASH // "Miscellaneous" // DASH // "PedigreeMistakes.txt","." // DASH // "Miscellaneous" // DASH // "snpMistakes.txt")
	call ped%outputSortedPedigreeInAlphaImputeFormat("." // DASH // "Miscellaneous" // DASH // "InternalDataRecoding.txt")

end subroutine CheckParentage


!#############################################################################################################################################################################################################################
SUBROUTINE SnpCallRate()
	use Global

	use AlphaImputeSpecFileModule
	use ISO_Fortran_Env

	implicit none

	integer :: i, CountMiss, UOutputs
	type(AlphaImputeInput), pointer :: inputParams

	inputParams => defaultInput


	open(newunit=UOutputs, file="." // DASH // "Miscellaneous" // DASH // "SnpCallRateByAnimal.txt",status='unknown')

	do i=1,ped%nGenotyped
		countMiss = ped%pedigree(ped%genotypeMap(i))%individualGenotype%numMissing()
		write(UOutputs,'(a20,6f5.1)') ped%pedigree(ped%genotypeMap(i))%originalID, (inputParams%nsnp-CountMiss)*100/real(inputParams%nsnp)
	end do
	close(UOutputs)

END SUBROUTINE SnpCallRate

!#############################################################################################################################################################################################################################

subroutine Titles

	call PrintVersion
	print *, ""
	print *, ""
	print *, ""

end subroutine Titles

!#############################################################################################################################################################################################################################

subroutine Header

	print *, ""
	print *, "                              ***********************                         "
	print *, "                              *                     *                         "
	print *, "                              *     AlphaImpute     *                         "
	print *, "                              *                     *                         "
	print *, "                              ***********************                         "
	print *, "                                                                              "
	print *, "                    Software For Phasing and Imputing Genotypes               "

end subroutine Header

!#############################################################################################################################################################################################################################

subroutine PrintVersion

	call Header
	print *, ""
	print *, "                              Commit:   "//TOSTRING(COMMIT),"                     "
	print *, "                              Compiled: "//__DATE__//", "//__TIME__
	print *, ""

end subroutine PrintVersion

!#############################################################################################################################################################################################################################

subroutine PrintTimerTitles
	use Global
	use iso_fortran_env
	implicit none

	real :: etime          ! Declare the type of etime()
	real :: elapsed(2)     ! For receiving user and system time
	real :: total,Minutes,Hours,Seconds

	print *, ""
	print *, ""
	call Header
	print*, ""
	print*, "                                  No Liability"
	print*, ""
	print*, "                Analysis Finished                         "

	total=etime(elapsed)
	Minutes=total/60
	Seconds=Total-(INT(Minutes)*60)
	Hours=Minutes/60
	Minutes=INT(Minutes)-(INT(Hours)*60)
	print '(A107,A7,I3,A9,I3,A9,F6.2)', "Time Elapsed","Hours", INT(Hours),"Minutes",INT(Minutes),"Seconds",Seconds

	open (unit=32,file="." // DASH // "Miscellaneous" // DASH // "Timer.txt",status="unknown")

	write(32,'(A27,A7,I3,A9,I3,A9,F6.2)') "Time Elapsed","Hours", INT(Hours),"Minutes",INT(Minutes),"Seconds",Seconds
	close(32)
end subroutine PrintTimerTitles



subroutine runAlphaImpute(in, pedIn)

	use Global
	use AlphaImputeInputOutputModule
	use AlphaImputeSpecFileModule
	use Imputation
	use ModuleRunFerdosi
	use AlphaPhaseResultsModule

	class(baseSpecFile), target :: in
	type(pedigreeHolder), optional :: pedIn






	select type(in)

	type is (AlphaImputeInput)
	defaultInput => in
	class default
	write(error_unit, *) "ERROR: AlphaImpute given correct object type as input"
	call abort()
end select



inputParams => defaultInput


if (inputParams%hmmoption /= RUN_HMM_NGS) then
	if (inputParams%restartOption<OPT_RESTART_IMPUTATION) call MakeDirectories(RUN_HMM_NULL)

	if (.not. present(pedIn)) then
		call ReadInData
	else
		ped = pedIn
	endif

	inputParams%nSnpRaw = inputParams%nsnp
	call SnpCallRate

			
	call CheckParentage
	if (inputParams%MultiHD/=0) then
		call ClassifyAnimByChips
	endif

	call FillInSnp

	call FillInBasedOnOffspring
	call InternalEdit

	if (inputParams%PreProcess==.true.) then
		print*, "Data preprocessed"
		stop
	endif
	allocate(GlobalWorkPhase(0:ped%pedigreeSize,inputParams%nsnpraw,2))
	call InitialiseArrays

	if (inputParams%hmmoption==RUN_HMM_ONLY) then

		print*, ""
		print*, "Bypass calculation of probabilities and phasing"

	else ! if hmm option is not ngs or only
		write(6,*) " "
		write(6,*) " ","Data editing completed"

		if (inputParams%useFerdosi) then

			call doFerdosi(ped)
		endif



		if (inputParams%managephaseon1off0==1) then


			if (inputParams%restartOption<OPT_RESTART_IMPUTATION) Then

				if (inputParams%cluster) then
#ifdef MPIACTIVE
					call phasingManagementCluster
#else
					write(error_unit,*) "WARNING: CLUSTER HAS BEEN SPECIFIED BUT MPI NOT ENABLED. Falling back on OpenMP version"
					call PhasingManagementNew(APResults)
#endif
				else
					call PhasingManagementNew(APResults)
				endif

			endif
		endif
	endif

	if (inputParams%restartOption> OPT_RESTART_PHASING) Then
		print *,"Reading in Phasing information"
		! Read back in geneprob data

		if (inputParams%managephaseon1off0==1) then
			inputParams%phasePath = "." // DASH //"Phasing"
		endif
		block

			use OutputParametersModule
			use InputOutput
			use AlphaPhaseResultsModule
			integer :: i
			type(OutputParameters) :: oParams
			oParams = newOutputParametersImpute()
			ApResults%nResults = inputparams%nPhaseInternal
			allocate(ApResults%results(ApResults%nResults))
			do i=1, ApResults%nResults
				write(oParams%outputDirectory,'(a,a,a,"Phase"i0)') trim(inputParams%phasePath),DASH,DASH, i
				call readAlphaPhaseResults(ApResults%results(i), oParams, ped)
			enddo
		end block
	endif
	print *, "Phasing Completed"

	! If we only want to phase data, then skip all the imputation steps
	if (inputParams%PhaseTheDataOnly==0) Then
		call ImputationManagement
		call WriteOutResults

		! WARNING: Skip the modelling the recombination because it interferes with HMM propabilites
		! TODO:
		if (.not. inputparams%ModelRecomb .or. inputParams%hmmoption /= RUN_HMM_NO) then
			write(*,*) "ModelRecomb has been Bypassed"
		else
			call ModelRecomb
		endif


		if (inputParams%TrueGenos1None0==1) then
			block
				use informationModule

				print *,""
				print *,"**************************************************************************************************"
				print *, "Yield", checkYield(ped)
				print *,"Accuracy per animal:",calculateaccuracyPerAnimal(ped,inputParams%TrueGenotypeFile, "perAnimal.txt", "Miscellaneous"// DASH// "ImputationErrors.txt")
			end block
		endif

	endif

else if (inputParams%hmmoption == RUN_HMM_NGS) then
	call MakeDirectories(RUN_HMM_NGS)
	call ReadInData
	call SnpCallRate
	allocate(SnpIncluded(inputParams%nsnp))
	call CheckParentage
	call ped%addSequenceFromFile(inputparams%GenotypeFile, nsnps=inputParams%nsnpRaw, maximumReads=MAX_READS_COUNT)


	block
		use AlphaHmmInMod
		use ExternalHMMWrappers
		type (AlphaHMMinput) :: inputParamsHMM
		integer(kind=1) ,dimension(:,:), allocatable :: res

		inputParamsHMM%nsnp = inputParams%nsnp
		inputParamsHMM%nHapInSubH = inputParams%nHapInSubH
		inputParamsHMM%HmmBurnInRound = inputParams%HmmBurnInRound
		inputParamsHMM%nRoundsHmm = inputParams%nRoundsHmm
		inputParamsHMM%useProcs = inputParams%useProcs
		inputParamsHMM%imputedThreshold = inputParams%imputedThreshold
		inputParamsHMM%phasedThreshold = inputParams%phasedThreshold
		inputParamsHMM%HapList = inputParams%HapList

		res = ped%getGenotypesAsArray()
		call AlphaImputeHMMRunner(inputParamsHMM, ped, ProbImputeGenosHmm, ProbImputePhaseHmm, GenosCounts, FullH)


	end block

	call FromHMM2ImputePhase
	call WriteOutResults

endif

! call ped%destroyPedigree()
call PrintTimerTitles



end subroutine runAlphaImpute

end module AlphaImputeModule





