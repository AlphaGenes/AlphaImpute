#!/usr/bin/python2.7
# encoding: utf-8
'''
alphaimpute -- Generate the Spec file for AlphaImpute
Different
@author:     Roberto Antolín
@copyright:  2015 Roberto Antolín. All rights reserved.
@license:    license
@contact:    roberto dot antolin at roslin dot ed dot ac dot uk
@deffield    updated: Updated
'''

import sys, os
import subprocess

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.2
__date__ = '2014-12-07'
__updated__ = '2015-03-10'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main(argv=None):    # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
  Created by rantolin on %s.
  Copyright 2015 Roslin Institute. All rights reserved.
  Licensed under the GPL 3.0v
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument(dest="output", help="Output file", metavar="file", type=str)
        parser.add_argument("-P", "--pedigree", dest="pedigree", help="File containing the pedigree information", metavar="file", type=file)
        parser.add_argument("-G", "--genotype", dest="genotype", help="File containing the genotypes", metavar="file", type=file)
        parser.add_argument("-X", "--sexchrom", help="Sex Chromosome File", metavar="file")
        parser.add_argument("-f", "--female", help="Set Female as the heterogametic sex. [Default: Male]", action="store_true")
        parser.add_argument("-S", "--snp", help="Number of SNP in the genotype file", type=int, metavar="nSNP", required=True)
        parser.add_argument("-E", "--edit", help="Edit data internally", action="store_true", dest="edit")
        parser.add_argument("-e", "--edit_param", help="Parameters to edit data internally [Default: %(default)s]", metavar="int", nargs=3, default=[95.0,2.0,98.0], dest='editParam')
        parser.add_argument("-o", "--edit-output", help="Output of the editing phase [Default: %(default)s]", choices=["AllSnpOut", "EditedSnpOut"], default="AllSnpOut", metavar="str")
        parser.add_argument("--phased", help="Specify if phasing rounds have been done previously", choices=["PhaseDone","NoPhase"], metavar="str")
        parser.add_argument("--phasepath", help="Path where the phasing rounds are store", metavar="path")
        parser.add_argument("-r", "--phasing_runs", help="Number of phasing runs", default=20, type=int, metavar="int")
        parser.add_argument("-t", "--tiles", help="Core and Tail lengths [Default: %(default)s]", nargs="*", default=[600,700,800,900,1000,1100,1200,1300,1600,1800], type=int, metavar="int")
        parser.add_argument("-c", "--cores", help="Core lengths [Default: %(default)s]", nargs="*", type=int, default=[500,600,700,800,900,1000,1100,1200,1400,1600], metavar="int")
        parser.add_argument("-F", "--freephasing", help="Pedigree free phasing", action="store_true")
        parser.add_argument("-g", "--genotype-error", help="Genotype Error [Default: %(default)3.1f]", type=float, default=0.0, metavar="float")
        parser.add_argument("-n", "--processors", help="Number of Processors Available [Default: %(default)d]", default=2, type=int, metavar="int")
        parser.add_argument("-i", "--iterations", help="Internal Iterations [Default: %(default)d]", type=int, default=3, metavar="int")
        parser.add_argument("--data-only", help="Pre process data only", action="store_true", dest="data")
        parser.add_argument("--phase-only", help="Phase Only", action="store_true", dest="phase")
        parser.add_argument("-l", "--library", help="Conservative Haplotype Library use", action="store_true")
        parser.add_argument("-w", "--well_phased_thres", help="Well phase threshold [Default: %(default)3.1f]", type=float, default=99.0, dest="wellthres", metavar="float")
        parser.add_argument("-U", "--user_phase", help="User defined AlphaPhase animals file", type=file, metavar="file", dest="userfile")
        parser.add_argument("-p", "--prephase", help="Pre-phased file", type=file, metavar="file")
        parser.add_argument("-b", "--bypass", help="Bypass GeneProb", action="store_true")
        parser.add_argument("-R", "--restart", help="Restart Option [Default: %(default)d]", type=int, default=0, required=True, metavar="int")
        parser.add_argument("-M", "--hmm", help="Use Hidden Markov Model [Default: %(default)s]", dest="hmm", choices=["No","Only", "Yes"], metavar="str")
        parser.add_argument("-m", "--hmm_param", help="Hidden Markov Model parameters[Default: %(default)s]", nargs=5, default=[300,19,20,4,-123456788], metavar="int", dest="hmmParam")

        parser.add_argument("-T", "--truegenotype", help="True Genotype File", metavar="file", type=file)
        parser.add_argument("-V", "--version", action="version", version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        outputFile = args.output
        pedigreeFile = args.pedigree
        genotypeFile = args.genotype
        sexChromosome = args.sexchrom
        female = args.female
        nSnps = args.snp
        edit = args.edit
        editParameters = args.editParam
        editOutput = args.edit_output
        phased = args.phased
        phasePath = args.phasepath
        phaseRuns = args.phasing_runs
        tilesLength = args.tiles
        coresLength = args.cores
        freePhasing = args.freephasing
        genotypeError = args.genotype_error
        nProcessors = args.processors
        nIterations = args.iterations
        dataOnly = args.data
        phaseOnly = args.phase
        library = args.library
        wellPhasedThres = args.wellthres
        userPhaseFile = args.userfile
        prePhasedFile = args.prephase
        bypass = args.bypass
        restartOption = args.restart
        hmm = args.hmm
        hmmParameters = args.hmmParam
        trueGenotypeFile = args.truegenotype


        # Construct file
        spec= 'PedigreeFile\t\t\t\t,{0}\n'.format(pedigreeFile.name)
        spec+= 'GenotypeFile\t\t\t\t,{0}\n'.format(genotypeFile.name)
        if sexChromosome is None:
            spec+= 'SexChrom\t\t\t\t,No\n'
        elif os.path.isfile(sexChromosome):
            heterogamete='Male'
            if female: heterogamete='Female'
            spec+= 'SexChrom\t\t\t\t,Yes,{0},{1}\n'.format(sexChromosome,heterogamete)
        else:
            parser.error(program_name + ": " + '<{0}> is not a valid sex chromosome file'.format(sexChromosome))
        spec+= 'NumberSnp\t\t\t\t,{0}\n'.format(nSnps)
        if edit:
            spec+= 'InternalEdit\t\t\t\t,Yes\n'
            spec+= 'EditingParameters\t\t\t,'
            for param in editParameters:
                spec+= str(param) + ','
            spec+= '{0}\n'.format(editOutput)
        else:
            spec+= 'InternalEdit\t\t\t\t,No\n'
            spec+= 'EditingParameters\t\t\t,0.0,0.0,0.0,AllSnpOut\n'
        if phased=='PhaseDone':
            if phasePath is not None and os.path.isfile(phasePath):
                spec+= 'NumberPhasingRuns\t\t\t,{0},{1},{2}\n'.format(phased,phasePath,phaseRuns)
            else:
                parser.error(program_name + ": " + '<{0}> is not a valid file name for phased data'.format(phasePath))
        elif phased=='NoPhase':
            spec+= 'NumberPhasingRuns\t\t\t,{0}\n'.format(phased)
        else:
            spec+= 'NumberPhasingRuns\t\t\t,{0}\n'.format(phaseRuns)

        spec+= 'CoreAndTailLengths\t\t\t'
        for tile in tilesLength:
            spec+= ',' + str(tile)
        spec+= '\n'

        spec+= 'CoreLengths\t\t\t\t'
        for core in coresLength:
            spec+= ',' + str(core)
        spec+= '\n'

        if freePhasing:
            spec+= 'PedigreeFreePhasing\t\t\t,Yes\n'
        else:
            spec+= 'PedigreeFreePhasing\t\t\t,No\n'

        spec+= 'GenotypeError\t\t\t\t,{0}\n'.format(genotypeError)

        spec+= 'NumberOfProcessorsAvailable\t\t,{0}\n'.format(nProcessors)
        spec+= 'InternalIterations\t\t\t,{0}\n'.format(nIterations)

        if dataOnly:
            spec+= 'PreprocessDataOnly\t\t\t,Yes\n'
        else:
            spec+= 'PreprocessDataOnly\t\t\t,No\n'

        if phaseOnly:
            spec+= 'PhasingOnly\t\t\t\t,Yes\n'
        else:
            spec+= 'PhasingOnly\t\t\t\t,No\n'

        if library:
            spec+= 'ConservativeHaplotypeLibraryUse\t\t,Yes\n'
        else:
            spec+= 'ConservativeHaplotypeLibraryUse\t\t,No\n'

        spec+= 'WellPhasedThreshold\t\t\t,{0}\n'.format(wellPhasedThres)

        if userPhaseFile is not None and os.path.isfile(userPhaseFile):
            spec+= 'UserDefinedAlphaPhaseAnimalsFile\t,{0}\n'.format(userPhaseFile)
        else:
            spec+= 'UserDefinedAlphaPhaseAnimalsFile\t,None\n'

        if prePhasedFile is not None and os.path.isfile(prePhasedFile):
            spec+= 'PrePhasedFile\t\t\t\t,{0}\n'.format(prePhasedFile)
        else:
            spec+= 'PrePhasedFile\t\t\t\t,None\n'

        if bypass:
            spec+= 'BypassGeneProb\t\t\t\t,Yes\n'
        else:
            spec+= 'BypassGeneProb\t\t\t\t,No\n'


        spec+= 'RestartOption\t\t\t\t,{0}\n'.format(restartOption)

        if hmm=="Yes":
            spec+= 'HMMOption\t\t\t\t,Yes\n'
        elif hmm=="Only":
            spec+= 'HMMOption\t\t\t\t,Only\n'
        else:
            spec+= 'HMMOption\t\t\t\t,No\n'

        spec+= 'HmmParameters\t\t\t\t'
        for param in hmmParameters:
            spec+= ',' + str(param)
        spec+= '\n'

        if trueGenotypeFile is not None and os.path.isfile(trueGenotypeFile.name):
            spec+= 'TrueGenotypeFile\t\t\t,{0}\n'.format(trueGenotypeFile.name)
        else:
            spec+= 'TrueGenotypeFile\t\t\t,None\n'


        # Write and close file
        # specFile = open(outputFile, 'w')
        with open(outputFile, 'w') as specFile:
            specFile.write(spec)

        specFile.close()

        # Bye bye
        print outputFile + " has been created successfully.\n"

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n")
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'argparse_module_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
