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
  implicit none
  PRIVATE

  PUBLIC compareHaplotype, compareHaplotypeAllowMissing
  PUBLIC BitCountAlleleImputed, BitCountAllelePhased
  
  INTERFACE compareHaplotype
    MODULE PROCEDURE compareHaplotypeThreshold, compareHaplotypeExtrict, compareHaplotype2Threshold, compareHaplotype2Extrict
  END INTERFACE compareHaplotype
!
  INTERFACE compareHaplotypeAllowMissing
    MODULE PROCEDURE compareHaplotypeAllowMissingThreshold, compareHaplotypeAllowMissingExtrict
  END INTERFACE compareHaplotypeAllowMissing

  
CONTAINS

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
  FUNCTION BitCountAllelePhased(miss1, miss2, nSecs) result(c)
    integer(kind=8), dimension(:), intent(in) :: miss1, miss2
    integer, intent(in) :: nSecs
    integer :: c

    integer :: i

    c = 0
    do i = 1, nSecs
      c = c + POPCNT(NOT(IOR(miss1(i),miss2(i))))
    end do

  END FUNCTION BitCountAllelePhased

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
  FUNCTION BitCountAlleleImputed(miss, nSecs) result(c)
    integer(kind=8), dimension(:), intent(in) :: miss
    integer, intent(in) :: nSecs
    integer :: c

    integer :: i

    c = 0
    do i = 1, nSecs
      c = c + POPCNT(NOT(miss(i)))
    end do

  END FUNCTION BitCountAlleleImputed

END MODULE HaplotypeBits
