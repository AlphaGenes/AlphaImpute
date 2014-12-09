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
        parser.add_argument("-i", "--pedigree" dest="predigree", help="Pedigree file", metavar="file", required=True)
        parser.add_argument("-g", "--genotype", dest='genotype', help="Genotype file", metavar="file", required=True)
        parser.add_argument("-x", "--sexchrom", help="Sex Chromosome File")   
        parser.add_argument("-f", "--hetero-female", help="Female is the heterogametic sex", choices=['Male', 'Female'])
        parser.add_argument("-s", "--snp", help="Number of SNPs", type=int, metavar="nSNP", required=True)
        parser.add_argument("-E", "--edit", help="Internal Edition of data")
        parser.add_argument("-n", "--phasing-run", help="Number of phasing runs")
        parser.add_argument("-t", "--tiles", help="Core and Tail lengths", type=list)
        parser.add_argument("-c", "--cores", help="Core lengths", type=list)
        parser.add_argument("-e", "--generror", help="Genotype Error [Default: %(default)3.1f]", type=float, default=0.0)
        parser.add_argument("-p", "--processors", help="Number of Processors Available [Default: %(default)d]", default=4, type=int)
        parser.add_argument("-I", "--iterations", help="Internal Iterations [Default: %(default)d]", type=int, default=5)
        parser.add_argument("-D", "--data-only", help="Process Data Only", action="store_true")
        parser.add_argument("-P", "--phase-only", help="Phase Only", action="store_true")
        parser.add_argument("-l", "--library-use", help="Conservative Haplotype Library Use", action="store_true")
        parser.add_argument("-w", "--well-phase", help="Well phase Threshold [Default: %(default)3.1f]", type=float, default=99.0)
        parser.add_argument("-u", "--user-phase", help="User Defined Alpha Phase Animals File")
        parser.add_argument("-R", "--prephase", help="Pre-phased file")
        parser.add_argument("-b", "--bypass", help="Bypass GeneProb", action="store_true")
        parser.add_argument("-r", "--restart", help="Restart Option", type=int)
        parser.add_argument("-H", "--HMM", help="Hidden Markov Model parameters", nargs=4)
        parser.add_argument("-T", "--truegenotype", help="True Genotype File", metavar="file")
        # parser.add_argument("-h", "--help", help="Display this help", action="help")
        parser.add_argument("-V", "--version", action="version", version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        pedigreeFile = args.pedigree


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