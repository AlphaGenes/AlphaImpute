#!/usr/bin/python2.7
# encoding: utf-8
'''
alphaimpute -- 
Different
@author:     Roberto Antolín
@copyright:  2014 Roberto Antolín. All rights reserved.
@license:    license
@contact:    roberto dot antolin at roslin dot ed dot ac dot uk
@deffield    updated: Updated
'''

import sys, os
import subprocess

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2014-12-07'
__updated__ = '2014-12-07'

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


# def runCommand(verb, cmd):
#     proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=False)
#     if verb == 2:
#         for line in iter(proc.stdout.readline, ''):
#             line = line.replace('\r', '').replace('\n', '')
#             print line
#             sys.stdout.flush()
#     elif verb == 1:
#         print cmd
#         proc.wait()
#     else:
#         proc.wait()


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
  Copyright 2014 Roslin Institute. All rights reserved.
  Licensed under the GPL 3.0v
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-P", "--pedigree", dest="pedigree", help="File containing the pedigree information", metavar="file", type=file)
        parser.add_argument("-G", "--genotype", dest='genotype', help="File containing the genotypes", metavar="file", type=file)
        parser.add_argument("-x", "--sexchrom", help="Sex Chromosome File", action="store_true")
        parser.add_argument("-f", "--heterogamete", help="The heterogametic sex [Default: %(default)s]", choices=['Female', 'Male'], default='Female')
        parser.add_argument("-S", "--snp", help="Number of SNP in the genotype file", type=int, metavar="nSNP", required=True)
        parser.add_argument("-e", "--edit", help="Parameters to edit data internally. First three parameters have to be float representing percentages. The last one can be one option between {AllSnpOut, EditedSnpOut} [Default: %(default)s]", metavar="", nargs=4, default=[95.0,2.0,98.0,'AllSnpOut']) # First three parameters are numerical (need to cast to float). The last one is a boolean
        parser.add_argument("--phased", help="Specify if phasing rounds have been done previously", choices=['PhaseDone','NoPhase'], metavar="")
        parser.add_argument("--phasepath", help="Path where the phasing rounds are store", metavar='PATH')
        parser.add_argument("-r", "--phasing_runs", help="Number of phasing runs", default=20, type=int, metavar='int')
        parser.add_argument("-t", "--tiles", help="Core and Tail lengths [Default: %(default)s]", nargs='*', default=[600,700,800,900,1000,1100,1200,1300,1600,1800], type=int, metavar='int')
        parser.add_argument("-c", "--cores", help="Core lengths [Default: %(default)s]", nargs='*', type=int, default=[500,600,700,800,900,1000,1100,1200,1400,1600], metavar='int')
        parser.add_argument("-E", "--genotype_error", help="Genotype Error [Default: %(default)3.1f]", type=float, default=0.0)
        parser.add_argument("-n", "--processors", help="Number of Processors Available [Default: %(default)d]", default=2, type=int)
        parser.add_argument("-i", "--iterations", help="Internal Iterations [Default: %(default)d]", type=int, default=3)
        parser.add_argument("--data-only", help="Pre process data only", action="store_true", dest="data")
        parser.add_argument("--phase-only", help="Phase Only", action="store_true", dest="phase")
        parser.add_argument("-l", "--library", help="Conservative Haplotype Library Use", action="store_true")
        parser.add_argument("-w", "--well_phased_thres", help="Well phase threshold [Default: %(default)3.1f]", type=float, default=99.0, dest="wellthres")
        parser.add_argument("-u", "--user_phase", help="User defined AlphaPhase animals file", type=file, metavar="file", dest="userfile")
        parser.add_argument("-p", "--prephase", help="Pre-phased file", type=file, metavar="file")
        parser.add_argument("-b", "--bypass", help="Bypass GeneProb", action="store_true")
        parser.add_argument("-R", "--restart", help="Restart Option [Default: %(default)d; ]", type=int, default=0, required=True)
        parser.add_argument("-H", "--HMM", help="Hidden Markov Model parameters", nargs=4, default=[2,3,5,-123456788], metavar="")
        parser.add_argument("-T", "--truegenotype", help="True Genotype File", metavar="file", type=file)
        parser.add_argument("-V", "--version", action="version", version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        pedigreeFile = args.pedigree
        genotypeFile = args.genotype
        sexChromosome = args.sexchrom
        heterogamete = args.heterogamete
        nSnps = args.snp
        editParameters = args.edit
        phased = args.phased
        phasePath = args.phasepath
        phaseRuns = args.phasing_runs
        tilesLength = args.tiles
        coresLength = args.cores
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
        hmm = args.HMM
        trueGenotypeFile = args.truegenotype

        with open('AlphaImputeSpec.txt', 'w') as spec:
            spec.write('PedigreeFile\t\t,{0}\n'.format(pedigreeFile.name))
            spec.write('GenotypeFile\t\t,{0}\n'.format(genotypeFile.name))


        print "[DONE]"

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
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