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

  PUBLIC compareHaplotype, compareHaplotypeAllowMissing, compareHaplotypeAllowSimultaneousMissing
  PUBLIC BitCountAllelesImputed, BitCountAllelesPhased, BitCountAllelesMissing
  PUBLIC BitCompletePhased, BitCompleteMissing
!  PUBLIC BitCountRefAlleles, BitCountAltAlleles

!  PUBLIC compareHaplotype2AllowMissingThreshold
  
  INTERFACE compareHaplotype
    MODULE PROCEDURE compareHaplotypeThreshold, compareHaplotypeExtrict, compareHaplotype2Threshold, compareHaplotype2Extrict
  END INTERFACE compareHaplotype

  INTERFACE compareHaplotypeAllowMissing
    MODULE PROCEDURE compareHaplotypeAllowMissingThreshold, compareHaplotypeAllowMissingExtrict
  END INTERFACE compareHaplotypeAllowMissing

  INTERFACE compareHaplotypeAllowSimultaneousMissing
    MODULE PROCEDURE compareHaplotypeAllowSimultaneousMissingThreshold, compareHaplotypeAllowSimultaneousMissingExtrict
  END INTERFACE compareHaplotypeAllowSimultaneousMissing

  TYPE, PUBLIC :: BitSection
    ! PUBLIC
    integer :: numSections
    integer :: overhang
  END TYPE BitSection
  
  INTERFACE BitSection
    MODULE PROCEDURE newBitSection
  END INTERFACE BitSection

  TYPE, PUBLIC :: HaplotypeBit
    integer :: numSections
    integer(kind = int64), allocatable, dimension(:,:) :: alt
    integer(kind = int64), allocatable, dimension(:,:) :: missing
  END TYPE

  INTERFACE HaplotypeBit
    MODULE PROCEDURE newHaplotypeBit
  END INTERFACE HaplotypeBit

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
  FUNCTION newHaplotypeBit(nsecs, ids) result(this)
    integer, intent(in) :: nsecs, ids
    type(HaplotypeBit)    :: this

    this%numSections = nsecs
    allocate(this%alt(nsecs,ids))
    allocate(this%missing(nsecs,ids))

    this%alt = 0
    this%missing = 0

  END FUNCTION newHaplotypeBit

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
  FUNCTION compareHaplotypeThreshold(hap1, hap2, miss1, miss2, nSecs, thres) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
    integer, intent(in) :: nSecs, thres
    logical :: same

    integer :: c, i

    c = 0
    same = .TRUE.
    do i = 1, nSecs
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
  FUNCTION compareHaplotypeExtrict(hap1, hap2, miss1, miss2, nSecs) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
    integer, intent(in) :: nSecs
    logical :: same

    same = .TRUE.
    same = compareHaplotypeThreshold(hap1, hap2, miss1, miss2, nSecs, 1)

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
  FUNCTION compareHaplotypeAllowMissingExtrict(hap1, hap2, miss1, miss2, nSecs) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
    integer, intent(in) :: nSecs
    logical :: same

    same = .TRUE.
    same = compareHaplotypeAllowMissingThreshold(hap1, hap2, miss1, miss2, nSecs, 1)

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
  FUNCTION compareHaplotypeAllowMissingThreshold(hap1, hap2, miss1, miss2, nSecs, thres) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
    integer, intent(in) :: nSecs, thres
    logical :: same

    integer :: c
    integer :: i

    c = 0
    same = .TRUE.
    do i = 1, nSecs
      c = c + POPCNT( IOR(                                                           &
                          IAND(IEOR(hap1(i), hap2(i)), NOT(IOR(miss1(i), miss2(i)))),&
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
  FUNCTION compareHaplotype2Threshold(hap1, hap2, miss, nSecs, thres) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss
    integer, intent(in) :: nSecs, thres
    logical :: same

    integer :: i, c

    c = 0
    same = .TRUE.
    do i = 1, nSecs
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
  FUNCTION compareHaplotype2Extrict(hap1, hap2, miss, nSecs) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss
    integer, intent(in) :: nSecs
    logical :: same

    same = .TRUE.
    same = compareHaplotype2Threshold(hap1, hap2, miss, nSecs, 1)

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
  FUNCTION compareHaplotypeAllowSimultaneousMissingThreshold(hap1, hap2, miss1, miss2, nSecs, thres) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
    integer, intent(in) :: nSecs, thres
    logical :: same

    integer :: c, i

    c = 0
    same = .TRUE.
    do i = 1, nSecs
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
  FUNCTION compareHaplotypeAllowSimultaneousMissingExtrict(hap1, hap2, miss1, miss2, nSecs) result(same)
    integer(kind=8), dimension(:), intent(in) :: hap1, hap2, miss1, miss2
    integer, intent(in) :: nSecs
    logical :: same

    same = .TRUE.
    same = compareHaplotypeAllowSimultaneousMissingThreshold(hap1, hap2, miss1, miss2, nSecs, 1)

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
  FUNCTION BitCountAllelesPhased(miss1, miss2, nSecs) result(c)
    integer(kind=8), dimension(:), intent(in) :: miss1, miss2
    integer, intent(in) :: nSecs
    integer :: c

    integer :: i

    c = 0
    do i = 1, nSecs
      c = c + POPCNT(NOT(IOR(miss1(i),miss2(i))))
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
  FUNCTION BitCountAllelesImputed(miss, nSecs) result(c)
    integer(kind=8), dimension(:), intent(in) :: miss
    integer, intent(in) :: nSecs
    integer :: c

    integer :: i

    c = 0
    do i = 1, nSecs
      c = c + POPCNT(NOT(miss(i)))
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
  FUNCTION BitCountAllelesMissing(miss, nSecs) result(c)
    integer(kind=8), dimension(:), intent(in) :: miss
    integer, intent(in) :: nSecs
    integer :: c

    integer :: i

    c = 0
    do i = 1, nSecs
      c = c + POPCNT(miss(i))
    end do

  END FUNCTION BitCountAllelesMissing

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Determine if the haplotype is phased
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
  FUNCTION BitCompletePhased(miss, nSecs) result(phased)
    integer(kind=8), dimension(:), intent(in) :: miss
    integer, intent(in) :: nSecs
    logical :: phased

    integer :: i

    phased = .TRUE.
    do i = 1, nSecs
      if ( POPCNT(miss(i)) > 0) then
        phased = .FALSE.
        exit
      end if
    end do

  END FUNCTION BitCompletePhased

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Determine if the haplotype is phased
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
  FUNCTION BitCompleteMissing(miss, nSecs) result(phased)
    integer(kind=8), dimension(:), intent(in) :: miss
    integer, intent(in) :: nSecs
    logical :: phased

    integer :: i

    phased = .TRUE.
    do i = 1, nSecs
      if ( POPCNT(NOT(miss(i))) > 0) then
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
  FUNCTION BitCountRefAlleles(hap, miss, nSecs) result(RefA)
    integer(kind=8), dimension(:), intent(in) :: hap, miss
    integer, intent(in) :: nSecs
    integer :: RefA

    integer :: i

    RefA = 0
    do i = 1, nSecs
      RefA = RefA + POPCNT(NOT(hap(i))) - POPCNT(miss(i))
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
  FUNCTION BitCountAltAlleles(hap, nSecs) result(AltA)
    integer(kind=8), dimension(:), intent(in) :: hap
    integer, intent(in) :: nSecs
    integer :: AltA

    integer :: i

    AltA = 0
    do i = 1, nSecs
      AltA = AltA + POPCNT(hap(i))
    end do

  END FUNCTION BitCountAltAlleles

END MODULE HaplotypeBits
