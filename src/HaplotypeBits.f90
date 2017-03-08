!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: PhaseRounds
!
!> @file        HaplotypeBits.f90
!
! DESCRIPTION: 
!> @brief       Module to bit-wise operations
!>
!> @details     This MODULE includes routines to work with haplotypes in a bit-wise way
!
!> @author      Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date        Aug 01, 2016
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016.08.01  RAntolin - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------
MODULE HaplotypeBits
    USE, INTRINSIC :: ISO_FORTRAN_ENV
    implicit none

    PRIVATE

    ! PUBLIC compareHaplotype, compareHaplotypeAllowMissing, compareHaplotypeAllowSimultaneousMissing
    ! PUBLIC BitCountAllelesImputed, BitCountAllelesPhased, BitCountAllelesMissing
    ! PUBLIC BitCompletePhased, BitCompleteMissing
    !  PUBLIC BitCountRefAlleles, BitCountAltAlleles

    !  PUBLIC compareHaplotype2AllowMissingThreshold



    !   INTERFACE compareHaplotype
    !   MODULE PROCEDURE compareHaplotypeThreshold, compareHaplotypeExtrict, compareHaplotype2Threshold, compareHaplotype2Extrict
    ! END INTERFACE compareHaplotype

    ! INTERFACE compareHaplotypeAllowMissing
    !   MODULE PROCEDURE compareHaplotypeAllowMissingThreshold, compareHaplotypeAllowMissingExtrict
    ! END INTERFACE compareHaplotypeAllowMissing

    ! INTERFACE compareHaplotypeAllowSimultaneousMissing
    !   MODULE PROCEDURE compareHaplotypeAllowSimultaneousMissingThreshold, compareHaplotypeAllowSimultaneousMissingExtrict
    ! END INTERFACE compareHaplotypeAllowSimultaneousMissing


    TYPE, PUBLIC :: BitSection
        ! PUBLIC
        integer :: numSections
        integer :: overhang

    contains
        private
        procedure :: compareHaplotypeThreshold
        procedure :: compareHaplotypeExtrict
        procedure :: compareHaplotype2Threshold
        procedure :: compareHaplotype2Extrict
        generic,public :: compareHaplotype => compareHaplotypeThreshold,compareHaplotypeExtrict, compareHaplotype2Threshold, compareHaplotype2Extrict

        procedure :: compareHaplotypeAllowMissingThreshold
        procedure :: compareHaplotypeAllowMissingExtrict
        generic,public ::  compareHaplotypeAllowMissing => compareHaplotypeAllowMissingThreshold, compareHaplotypeAllowMissingExtrict

        procedure :: compareHaplotypeAllowSimultaneousMissingThreshold
        procedure :: compareHaplotypeAllowSimultaneousMissingExtrict
        generic, public :: compareHaplotypeAllowSimultaneousMissing => compareHaplotypeAllowSimultaneousMissingThreshold, compareHaplotypeAllowSimultaneousMissingExtrict
        procedure, public :: BitCompletePhased
        procedure, public :: BitCompleteMissing
        procedure, public :: BitCountAllelesPhased
        procedure, public :: BitCountAllelesImputed
        procedure, public :: BitCountAllelesMissing
    END TYPE BitSection




    INTERFACE BitSection
        MODULE PROCEDURE newBitSection
    END INTERFACE BitSection

    ! TYPE, PUBLIC :: HaplotypeBit
    !   integer :: numSections
    !   integer(kind = int64), allocatable, dimension(:,:) :: alt
    !   integer(kind = int64), allocatable, dimension(:,:) :: missing
    ! END TYPE

    ! INTERFACE HaplotypeBit
    !   MODULE PROCEDURE newHaplotypeBit
    ! END INTERFACE HaplotypeBit

CONTAINS

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Initialise new Bit Haplotype
    !
    !> @details    Initialise new Bit Haplotype
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       August 09, 2016
    !
    ! PARAMETERS:
    !> @param[inout]  this  HaplotypeBit
    !---------------------------------------------------------------------------
    ! FUNCTION newHaplotypeBit(nsecs, ids) result(this)
    !   integer, intent(in) :: nsecs, ids
    !   type(HaplotypeBit)    :: this

    !   this%numSections = nsecs
    !   allocate(this%alt(nsecs,ids))
    !   allocate(this%missing(nsecs,ids))

    !   this%alt = 0
    !   this%missing = 0

    ! END FUNCTION newHaplotypeBit

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Initialise new bit section
    !
    !> @details    Initialise new bit section
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       August 02, 2016
    !
    ! PARAMETERS:
    !> @param[inout]  this  BitSection
    !---------------------------------------------------------------------------
    FUNCTION newBitSection(snps, bits) result(this)
        integer, intent(in) :: snps, bits
        type(BitSection)    :: this

        integer :: nsecs

        nsecs = snps / bits
        if (MOD(snps, bits)/=0) then
            nsecs = nsecs + 1
        end if

        this%numSections = nsecs
        this%overhang = bits - (snps - (nsecs - 1) * bits)

    END FUNCTION newBitSection

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least a number of alleles (thres) of hap1 is
    !>             different from the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles of both hap1 and hap2
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1  Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2  Bit-wise array with the second haplotype to compare
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @param[in]  thres  Number of alleles allow to differ
    !> @return     .TRUE. if haplotypes are the same, .FALSE. otherwise 
    !---------------------------------------------------------------------------  
    FUNCTION compareHaplotypeThreshold(this,hap1, hap2, miss1, miss2, thres) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
        integer, intent(in) :: thres
        logical :: same

        integer :: c, i

        c = 0
        same = .TRUE.
        do i = 1, this%numSections
            ! either hap(1) OR hap(2) (not both) and not miss(1) | miss(2)
            ! count the occurences of this
            c = c + POPCNT( IAND(IEOR(hap1(i), hap2(i)),&
                NOT(IOR(miss1(i), miss2(i)))))
            if ( c >= thres ) then
                same = .FALSE.
                exit
            end if
        end do

    END FUNCTION compareHaplotypeThreshold

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least one allele of hap1 is different from
    !>             the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles of both hap1 and hap2
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1  Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2  Bit-wise array with the second haplotype to compare
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     .TRUE. if haplotypes are the extrictly the same, .FALSE. otherwise 
    !---------------------------------------------------------------------------  
    FUNCTION compareHaplotypeExtrict(this,hap1, hap2, miss1, miss2) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
        logical :: same

        same = .TRUE.
        same = this%compareHaplotypeThreshold(hap1, hap2, miss1, miss2, 1)

    END FUNCTION compareHaplotypeExtrict

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least one allele of hap1 is different from
    !>             the homologous imputed allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles of hap2
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1  Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2  Bit-wise array with the second haplotype to compare
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     .TRUE. if haplotypes are the extrictly the same, .FALSE. otherwise 
    !---------------------------------------------------------------------------  
    FUNCTION compareHaplotypeAllowMissingExtrict(this,hap1, hap2, miss1, miss2) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
        logical :: same

        same = .TRUE.
        same = this%compareHaplotypeAllowMissingThreshold(hap1, hap2, miss1, miss2, 1)

    END FUNCTION compareHaplotypeAllowMissingExtrict

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least a number of alleles (thres) of hap1 is
    !>             different from the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles of hap2
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1  Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2  Bit-wise array with the second haplotype to compare
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @param[in]  thres  Number of alleles allow to differ
    !> @return     .TRUE. if haplotypes are the same, .FALSE. otherwise 
    !---------------------------------------------------------------------------  
    FUNCTION compareHaplotypeAllowMissingThreshold(this,hap1, hap2, miss1, miss2, thres) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
        integer, intent(in) :: thres
        logical :: same

        integer :: c
        integer :: i

        c = 0
        same = .TRUE.
        do i = 1, this%numSections
            c = c + POPCNT( IOR(IAND(IEOR(hap1(i), hap2(i)), NOT(IOR(miss1(i), miss2(i)))),&
                IAND(NOT(miss2(i)), miss1(i))                              &
                ) )

            if ( c >= thres ) then
                same = .FALSE.
                exit
            end if
        end do

    END FUNCTION compareHaplotypeAllowMissingThreshold

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least a number of alleles (thres) of hap1 is
    !>             different from the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles from alleles given by miss
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1   Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2   Bit-wise array with the second haplotype to compare
    !> @param[in]  miss   Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @param[in]  thres  Number of alleles allow to differ
    !> @return     .TRUE. if haplotypes are the same, .FALSE. otherwise
    !---------------------------------------------------------------------------
    FUNCTION compareHaplotype2Threshold(this,hap1, hap2, miss, thres) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss
        integer, intent(in) :: thres
        logical :: same

        integer :: i, c

        c = 0
        same = .TRUE.
        do i = 1, this%numSections
            c = c + POPCNT( IAND( IEOR(hap1(i), hap2(i)),&
                NOT(miss(i))))
            if ( c >= thres ) then
                same = .FALSE.
                exit
            end if
        end do

    END FUNCTION compareHaplotype2Threshold

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least an allele (thres) of hap1 is
    !>             different from the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles from alleles given by miss
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1   Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2   Bit-wise array with the second haplotype to compare
    !> @param[in]  miss   Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     .TRUE. if haplotypes are extrictly the same, .FALSE. otherwise
    !---------------------------------------------------------------------------
    FUNCTION compareHaplotype2Extrict(this,hap1, hap2, miss) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss
        logical :: same

        same = .TRUE.
        same = this%compareHaplotype2Threshold(hap1, hap2, miss, 1)

    END FUNCTION compareHaplotype2Extrict

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least a number of alleles (thres) of hap1 is
    !>             different from the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles, but not simultaneous
    !>                   missing alleles of both hap1 and hap2
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1  Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2  Bit-wise array with the second haplotype to compare
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @param[in]  thres  Number of alleles allow to differ
    !> @return     .TRUE. if haplotypes are the same, .FALSE. otherwise
    !---------------------------------------------------------------------------
    FUNCTION compareHaplotypeAllowSimultaneousMissingThreshold(this,hap1, hap2, miss1, miss2, thres) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
        integer, intent(in) :: thres
        logical :: same

        integer :: c, i

        c = 0
        same = .TRUE.
        do i = 1, this%numSections
            c = c + POPCNT( IAND(IEOR(hap1(i), hap2(i)),&
                NOT(IAND(miss1(i), miss2(i)))))
            if ( c >= thres ) then
                same = .FALSE.
                exit
            end if
        end do

    END FUNCTION compareHaplotypeAllowSimultaneousMissingThreshold

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Compare two haplotypes
    !
    !> @details    Return .FALSE. if at least one allele of hap1 is different from
    !>             the homologous allele of hap2, and .TRUE. otherwise.
    !>
    !>             Note: This function skips missing alleles, but not simultaneous
    !>                   missing alleles of both hap1 and hap2
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap1  Bit-wise array with the first haplotype to compare
    !> @param[in]  hap2  Bit-wise array with the second haplotype to compare
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     .TRUE. if haplotypes are the extrictly the same, .FALSE. otherwise
    !---------------------------------------------------------------------------
    FUNCTION compareHaplotypeAllowSimultaneousMissingExtrict(this,hap1, hap2, miss1, miss2) result(same)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
        logical :: same

        same = .TRUE.
        same = this%compareHaplotypeAllowSimultaneousMissingThreshold(hap1, hap2, miss1, miss2, 1)

    END FUNCTION compareHaplotypeAllowSimultaneousMissingExtrict


    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the total number of markers completely phased
    !
    !> @details    Count the total number of alleles completely phased
    !>             for a pair of haplotyes
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  miss1  Bit-wise array with the missing alleles of haplotype 1
    !> @param[in]  miss2  Bit-wise array with the missing alleles of haplotype 2
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     Number of markers phased
    !---------------------------------------------------------------------------
    FUNCTION BitCountAllelesPhased(this,miss1, miss2) result(c)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: miss1, miss2
        integer :: c

        integer :: i

        c = 0
        do i = 1, this%numSections
            c = c + POPCNT(NOT(IOR(miss1(i),miss2(i)))) - this%overhang
        end do

    END FUNCTION BitCountAllelesPhased

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the total number of alleles completely imputed
    !
    !> @details    Count the total number of alleles completely imputed in a haplotype
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 01, 2016
    !
    ! PARAMETERS:
    !> @param[in]  miss  Bit-wise array with the missing alleles of haplotype
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     Number of alleles imputed
    !---------------------------------------------------------------------------
    FUNCTION BitCountAllelesImputed(this,miss) result(c)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: miss
        integer :: c

        integer :: i

        c = 0
        do i = 1, this%numSections
            ! DW subtracted overhang fir not
            c = c + POPCNT(NOT(miss(i))) - this%overhang
        end do

    END FUNCTION BitCountAllelesImputed

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the total number of alleles missing
    !
    !> @details    Count the total number of alleles missing in a haplotype
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 02, 2016
    !
    ! PARAMETERS:
    !> @param[in]  miss  Bit-wise array with the missing alleles of haplotype
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     Number of alleles missing
    !---------------------------------------------------------------------------
    FUNCTION BitCountAllelesMissing(this,miss) result(c)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: miss
        integer :: c

        integer :: i

        c = 0
        do i = 1, this%numSections
            c = c + POPCNT(miss(i))
        end do

    END FUNCTION BitCountAllelesMissing

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Determine if the haplotype is not fully phased
    !
    !> @details    Determine if the haplotype is phased
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 02, 2016
    !
    ! PARAMETERS:
    !> @param[in]  miss  Bit-wise array with the missing alleles of haplotype
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     .TRUE. if haplotype is phased, .FALSE. otherwise
    !---------------------------------------------------------------------------
    FUNCTION BitCompletePhased(this,miss) result(phased)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: miss
        logical :: phased

        integer :: i

        phased = .TRUE.
        do i = 1, this%numSections
            if ( POPCNT(miss(i)) > 0) then
                phased = .FALSE.
                exit
            end if
        end do

    END FUNCTION BitCompletePhased

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Determine if the haplotype is fully phased
    !
    !> @details    Determine if the haplotype is fully phased
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 02, 2016
    !
    ! PARAMETERS:
    !> @param[in]  miss  Bit-wise array with the missing alleles of haplotype
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     .TRUE. if haplotype is phased, .FALSE. otherwise
    !---------------------------------------------------------------------------
    FUNCTION BitCompleteMissing(this,miss) result(phased)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: miss
        logical :: phased

        integer :: i

        phased = .TRUE.
        do i = 1, this%numSections
            ! DW added subtracting of overhang to deal with false population
            if ( (POPCNT(NOT(miss(i))) - this%overhang) > 0) then
                phased = .FALSE.
                exit
            end if
        end do
        !    phased = BitCompletePhased(NOT(miss), nSecs)

    END FUNCTION BitCompleteMissing

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the number of reference alleles
    !
    !> @details    Count the number of reference alleles (# 0-bits - # missing)
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 04, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap   Bit-wise array with the alternative alleles of haplotype
    !> @param[in]  miss  Bit-wise array with the missing alleles of haplotype
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     Number of reference alleles
    !---------------------------------------------------------------------------
    FUNCTION BitCountRefAlleles(this,hap, miss) result(RefA)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap, miss
        integer :: RefA

        integer :: i

        RefA = 0
        do i = 1, this%numSections
            ! DW added subtracting overhang to deal with notting of hap bits
            RefA = RefA + POPCNT(NOT(hap(i))) - POPCNT(miss(i)) - this%overhang 
        end do

    END FUNCTION BitCountRefAlleles

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the number of alternative alleles
    !
    !> @details    Count the number of alternative alleles (# 1-bits)
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date        Aug 04, 2016
    !
    ! PARAMETERS:
    !> @param[in]  hap   Bit-wise array with the alternative alleles of haplotype
    !> @param[in]  nSecs  Number of elements of the bit-wise arrays
    !> @return     Number of reference alleles
    !---------------------------------------------------------------------------
    FUNCTION BitCountAltAlleles(this,hap) result(AltA)
    class(BitSection) :: this
        integer(kind=int64), dimension(:), intent(in) :: hap
        integer :: AltA

        integer :: i

        AltA = 0
        do i = 1, this%numSections
            AltA = AltA + POPCNT(hap(i))
        end do

    END FUNCTION BitCountAltAlleles

END MODULE HaplotypeBits
