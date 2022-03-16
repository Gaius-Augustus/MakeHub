#!/usr/bin/env python3

# WARNINGS:
#
# 1) Display of hints has been adapted to BRAKER1 hints format
# 2) This script retrieves the binaries for linux and os x64 bit systems. It
#    will not work on other architectures and systems unless the required
#    UCSC tool binaries are already present.


import os
import errno
import os.path
import argparse
import re
import shutil
import platform
import urllib.request
import subprocess
from difflib import SequenceMatcher
from inspect import currentframe, getframeinfo
try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    frameinfo = getframeinfo(currentframe())
    raise ImportError('In file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' +
                      'Failed to import biophython modules. ' +
                      'Try installing with \"pip3 install biopython\"')


__author__ = "Kathairna J. Hoff"
__copyright__ = "Copyright 2019-2022. All rights reserved."
__credits__ = "Mario Stanke"
__license__ = "Artistic Licsense"
__version__ = "1.0.6"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

ucsc_tools = {'bedToBigBed': '', 'genePredCheck': '', 'faToTwoBit': '',
              'gtfToGenePred': '', 'hgGcPercent': '', 'ixIxx': '',
              'twoBitInfo': '', 'wigToBigWig': '', 'genePredToBed': '',
              'genePredToBigGenePred': ''}

ucsc_as_files = {'bigGenePred.as': None, 'cytoBand.as': None}

augustus_tools = {'bam2wig': ''}

parser = argparse.ArgumentParser(
    description='Generate UCSC assembly hub (e.g. from BRAKER or MAKER output).')
parser.add_argument('-p', '--printUsageExamples', action='store_true',
                    help="Print usage examples for make_hub.py")
parser.add_argument('-e', '--email', required=False, type=str,
                    help='Contact e-mail adress for assembly hub')
parser.add_argument('-g', '--genome', required=False, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('-L', '--long_label', required=False, type=str,
                    help='Long label for hub, e.g. species name in english ' +
                    'and latin, pass in single quotation marks, e.g. ' +
                    '--long_label \'Dorosphila melanogster (fruit fly)\'')
parser.add_argument('-l', '--short_label', required=False, type=str,
                    help='Short label for hub, will also be used as ' +
                    'directory name for hub, should not contain spaces or ' +
                    'special characters, e.g. --short_label fly')
parser.add_argument('-b', '--bam', required=False, type=str, nargs="+",
                    help='BAM file(s) - space separated - with RNA-Seq ' +
                    'information, by default will be displayed as bigWig')
parser.add_argument('-c', '--cores', required=False, type=int, default=1,
                    help='Number of cores for samtools sort processes')
parser.add_argument('-d', '--display_bam_as_bam', action='store_true',
                    help="Display BAM file(s) as bam tracks")
parser.add_argument('-E', '--gemoma_filtered_predictions', required=False, type=str,
                    help="GFF3 output file of Gemoma")
parser.add_argument('-X', '--braker_out_dir', required=False, type=str,
                    help="BRAKER output directory with GTF files")
parser.add_argument('-M', '--maker_gff', required=False, type=str,
                    help='MAKER2 output file in GFF3 format')
parser.add_argument('-I', '--glimmer_gff', required= False, type=str, 
                    help="GFF3 output file of GlimmerHMM")
parser.add_argument('-S', '--snap_gff', required=False, type = str,
                    help='SNAP output file in GFF3 format')
parser.add_argument('-a', '--annot', required=False, type=str,
                    help='GTF file with reference annotation')
parser.add_argument('-G', '--gene_track', required=False, nargs='+',
                    help="Gene track with user specified label, argument " +
                    "must be formatted as  follows: --gene_track file.gtf tracklabel")
# The following argument is for adding a track to an existing hub, i.e.
# producing the files required for the track and writing into existing
# configuration files; requires the directory tmp to be present
parser.add_argument('-A', '--add_track', action='store_true',
                    help="Add track(s) to existing hub")
parser.add_argument('-o', '--outdir', required=False, type=str, default='.',
                    help="output directory to write hub to")
parser.add_argument('-n', '--no_repeats', action='store_true',
                    help="Disable repeat track generation from softmasked " +
                    "genome sequence (saves time)")
parser.add_argument('-s', '--SAMTOOLS_PATH', required=False, type=str,
                    help="Path to samtools executable")
parser.add_argument('-B', '--BAM2WIG_PATH', required=False, type=str,
                    help="Path to bam2wig executable")
parser.add_argument('-i', '--hints', required=False, type=str,
                    help='GFF file with AUGUSTUS hints')
parser.add_argument('-t', '--traingenes', required=False, type=str,
                    help='GTF file with training genes')
parser.add_argument('-m', '--genemark', required=False, type=str,
                    help='GTF file with GeneMark predictions')
parser.add_argument('-w', '--aug_ab_initio', required=False, type=str,
                    help='GTF file with ab initio AUGUSTUS predictions')
parser.add_argument('-x', '--aug_hints', required=False, type=str,
                    help='GTF file with AUGUSTUS predictions with hints')
parser.add_argument('-y', '--aug_ab_initio_utr', required=False, type=str,
                    help='GTF file with ab initio AUGUSTUS predictions with ' +
                    'UTRs')
parser.add_argument('-z', '--aug_hints_utr', required=False, type=str,
                    help='GTF file with AUGUSTUS predictions with hints with ' +
                    'UTRs')
parser.add_argument('-N', '--latin_name', required=False, type=str,
                    help="Latin species name, e.g. \"Drosophila melanogaster\". This " +
                    "argument must be provided if the hub is supposed to be added to the " +
                    "public UCSC list.")
parser.add_argument('-V', '--assembly_version', required=False, type=str,
                    help="Assembly version, e.g. \"BDGP R4/dm3\". This " +
                    "argument must be provided if the hub is supposed to be added to the " +
                    "public UCSC list.")
parser.add_argument('-r', '--no_tmp_rm', action='store_true',
                    help="Do not delete temporary files ")
parser.add_argument('-P', '--no_genePredToBigGenePred', action='store_true',
                    help='Option for the special case in which the precompiled' +
                    ' UCSC binaries are not working on your system, and you installed ' +
                    'kentutils from the older ENCODE github repository; if activated, ' +
                    'gene prediction tracks will be output to bigBed instead of bigGenePred ' +
                    'format and amino acid display will not be possible in gene tracks.')
parser.add_argument('-u', '--verbosity', required=False, type=int, default=0,
                    help="If INT>0 verbose output log is produced")
parser.add_argument('-v', '--version', action='version',
                    version='%(prog)s ' + __version__)
args = parser.parse_args()


if args.verbosity > 3:
    print(args)

''' Print usage examples if that is desired '''

if args.printUsageExamples:
    print("\nUsage example for generating a novel hub:\n")
    print("make_hub.py -l hmi2 -L \"Rodent tapeworm\" -g data/genome.fa -e " +
          "anonymous@anonymous.de -a data/annot.gtf -b data/rnaseq.bam -d\n\n")
    print("Usage example for adding a gene prediction track to an existing " +
          "hub (hub resides in the directory where this command is executed):\n")
    print("make_hub.py -l hmi2 -e anonymous@anonymous.de -i data/hintsfile.gff " +
          "-A -M data/maker.gff -X data\n")
    exit(1)
else:
    if args.email is None or args.short_label is None:
        print("Error: the following arguments are required: -e/--email, -l/--short_label")
        exit(1)


''' Check whether sufficient options have been provided '''

if (args.long_label is None) and (args.short_label is not None) and (args.add_track is False):
    print("Warning: no long label specified for creating novel track hub, " +
          "will use short label \"" + args.short_label + "\" as long label!")
    args.long_label = args.short_label

if ((args.email is None) or (args.genome is None) or (args.short_label is None)) and (args.add_track is False) and (args.printUsageExamples is False):
    frameinfo = getframeinfo(currentframe())
    print('Usage error in file ' + frameinfo.filename + ' at line ' + str(frameinfo.lineno) + ': '
          + 'If a novel track is created, the following arguments are ' +
          'required: -e/--email, -g/--genome, -l/--short_label')
    exit(1)


''' Set args in case args.braker_out_dir is specified and they are otherwise unset '''

if args.braker_out_dir:
    if re.search(r'[^/]$', args.braker_out_dir):
        args.braker_out_dir = args.braker_out_dir + "/"
    if os.path.isdir(args.braker_out_dir):
        # set args.hints
        if (args.hints is None) and os.path.isfile(args.braker_out_dir + "hintsfile.gff"):
            args.hints = args.braker_out_dir + "hintsfile.gff"
        # set args.traingenes
        if args.traingenes is None:
            if os.path.isfile(args.braker_out_dir + "GeneMark-ES/genemark.d.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-ES/genemark.d.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ES/genemark.f.good.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-ES/genemark.f.good.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ET/genemark.d.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-ET/genemark.d.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ET/genemark.f.good.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-ET/genemark.f.good.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-EP/genemark.d.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-EP/genemark.d.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-EP/genemark.f.good.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-EP/genemark.f.good.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ETP/genemark.d.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-ETP/genemark.d.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ETP/genemark.f.good.gtf"):
                args.traingenes = args.braker_out_dir + "GeneMark-ETP/genemark.f.good.gtf"
            elif os.path.isfile(args.braker_out_dir + "gthTrainGenes.gtf"):
                args.traingenes = args.braker_out_dir + "gthTrainGenes.gtf"
            else:
                print("Warning: did not set args.traingenes because none of the " +
                      "expected training gene structure files exists in " +
                      args.braker_out_dir)
        # set args.genemark
        if args.genemark is None:
            if os.path.isfile(args.braker_out_dir + "GeneMark-ES/genemark.gtf"):
                args.genemark = args.braker_out_dir + "GeneMark-ES/genemark.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ET/genemark.gtf"):
                args.genemark = args.braker_out_dir + "GeneMark-ET/genemark.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-EP/genemark.gtf"):
                args.genemark = args.braker_out_dir + "GeneMark-EP/genemark.gtf"
            elif os.path.isfile(args.braker_out_dir + "GeneMark-ETP/genemark.gtf"):
                args.genemark = args.braker_out_dir + "GeneMark-ETP/genemark.gtf"
            else:
                print("Warning: did not set args.genemark because none of the " +
                      "expected genemark.gtf files exists in " +
                      args.braker_out_dir)
        # set args.aug_ab_initio
        if args.aug_ab_initio is None:
            if os.path.isfile(args.braker_out_dir + "augustus.ab_initio.gtf"):
                args.aug_ab_initio = args.braker_out_dir + "augustus.ab_initio.gtf"
            else:
                print("Warning: did not set args.aug_ab_initio because the " +
                      "file augustus.ab_initio.gtf does not exist in " +
                      args.braker_out_dir)
        # set args.aug_hints
        if args.aug_hints is None:
            if os.path.isfile(args.braker_out_dir + "augustus.hints.gtf"):
                args.aug_hints = args.braker_out_dir + "augustus.hints.gtf"
            else:
                print("Warning: did not set args.aug_hints because the " +
                      "file augustus.hints.gtf does not exist in " +
                      args.braker_out_dir)
        # set args.aug_ab_initio_utr
        if args.aug_ab_initio_utr is None:
            if os.path.isfile(args.braker_out_dir + "augustus.ab_initio_utr.gtf"):
                args.aug_ab_initio_utr = args.braker_out_dir + "augustus.ab_initio_utr.gtf"
            else:
                print("Warning: did not set args.aug_ab_initio_utr because the " +
                      "file augustus.ab_initio_utr.gtf does not exist in " +
                      args.braker_out_dir)
        # set args.aug_hints_utr
        if args.aug_hints_utr is None:
            if os.path.isfile(args.braker_out_dir + "augustus.hints_utr.gtf"):
                args.aug_hints_utr = args.braker_out_dir + "augustus.hints_utr.gtf"
            else:
                print("Warning: did not set args.aug_hints_utr because the " +
                      "file augustus.hints_utr.gtf does not exist in " +
                      args.braker_out_dir)


''' Check whether args.outdir exists and is writable '''

if not os.path.isdir(args.outdir):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --outdir ' +
          "(" + args.outdir + ") is not a directory!")
    exit(1)
elif not os.access(args.outdir, os.W_OK):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --outdir ' +
          "(" + args.outdir + ") is not writable!")
    exit(1)

tmp_dir = args.outdir + "/tmp-" + args.short_label + "/"
hub_dir = args.outdir + "/" + args.short_label + \
    "/" + args.short_label + "/"


''' Find samtools (if bam file provided) '''

samtools = ""
if args.bam:
    if args.verbosity > 0:
        print("Searching for samtools:")
if args.bam and args.SAMTOOLS_PATH:
    samtools = args.SAMTOOLS_PATH + "/samtools"
    if not(os.access(samtools, os.X_OK)):
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + samtools + " is not executable!")
        exit(1)
    else:
        if args.verbosity > 0:
            print("Will use " + samtools)
elif args.bam:
    if shutil.which("samtools") is not None:
        samtools = shutil.which("samtools")
        if args.verbosity > 0:
            print("Will use " + samtools)
    else:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': '
              + "Unable to locate samtools binary!")
        print("samtools is available as package in many Linux distributions.")
        print("For example, on Ubuntu, try installing with:")
        print("\"sudo apt install samtools\"")
        print("If samtools is unavailable as a package, you can obtain it " +
              "from github at:")
        print("https://github.com/samtools/samtools")
        exit(1)


''' Find either Augustus auxprog bam2wig (optional dependency) '''

bam2wig_aug = False
if args.bam:
    if args.verbosity > 0:
        print("Searching for bam2wig from AUGUSTUS auxprogs:")
    if args.BAM2WIG_PATH is not None:
        if os.access(args.BAM2WIG_PATH + "/bam2wig", os.X_OK):
            augustus_tools['bam2wig'] = args.BAM2WIG_PATH + "/bam2wig"
            bam2wig_aug = True
    elif shutil.which("bam2wig") is not None:
        augustus_tools['bam2wig'] = shutil.which("bam2wig")
        bam2wig_aug = True
    else:
        frameinfo = getframeinfo(currentframe())
        print('Was unable to locate bam2wig from ' +
              'AUGUSTUS auxprogs (available at ' +
              'https://github.com/Gaius-Augustus/Augustus). ' +
              'Will convert bam to wig using samtools pileup ' +
              'and line-by-line processing, instead.')


''' Find gzip '''

gzip_tool = ""
if (not args.add_track) or args.bam:
    gzip_tool = shutil.which('gzip')
    if gzip_tool is None:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' +
              "Unable to locate gzip")
        print("gzip is available as package in many Linux distributions.")
        print("For example, pn Ubuntu, try installing with:")
        print("\"sudo apt install gzip\"")
        print("If gzip is unavailable as a package, you can obtain it " +
              "from:")
        print("https://ftp.gnu.org/gnu/gzip/")
        quit(1)


''' Find bash sort '''

sort_tool = shutil.which('sort')
if sort_tool is None:
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': ' + "Unable to locate bash tool 'sort'")
    print('sort is part of most Linux distributions. On Ubuntu, it is ' +
          'part of the package coreutils. Try re-installing your bash if sort' +
          ' is missing on your system.')
    quit(1)


''' Find or obtain UCSC tools '''
# the URLs of UCSC tool download are hardcoded for linux.x84_64

arch = platform.machine()
plat_sys = platform.system()

if args.verbosity > 0:
    print("Searching for required UCSC tools:")
for key, val in ucsc_tools.items():
    # genePredToBigBed is optional because some users might not have easy access
    if key == 'genePredToBigBed' and no_genePredToBigBed:
        continue
    if shutil.which(key) is not None:
        ucsc_tools[key] = shutil.which(key)
    elif os.path.isfile(os.getcwd() + "/" + key):
        ucsc_tools[key] = os.getcwd() + "/" + key
        if not(os.access(ucsc_tools[key], os.X_OK)):
            os.chmod(ucsc_tools[key], 0o777)
    else:
        if not(arch == 'x86_64'):
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "This script depends on binaries that are available for " +
                  "x86_64 architecture for linux and MacOX. " +
                  "We have determined that your system architecture is " +
                  arch + "." +
                  " Please try downloading " + key +
                  " for your architecture from: " +
                  "http://hgdownload.soe.ucsc.edu/admin/exe")
            exit(1)
        elif not(plat_sys == 'Linux') and not(plat_sys == 'Darwin'):
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "This script depends on binaries that are available for " +
                  "x86_64 architecture for linux and MacOX. " +
                  "We have determined that your system is " +
                  plat_sys + "." +
                  " Please try downloading " + key +
                  " for your operating system from: " +
                  "http://hgdownload.soe.ucsc.edu/admin/exe")
            exit(1)
        else:
            if plat_sys == 'Linux':
                tool_url = "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/" + key
            elif plat_sys == 'Darwin':
                tool_url = "http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/" + key
            print("Was unable to locate " + key +
                  " on your system, will try to download it from " + tool_url + "...")
            with urllib.request.urlopen(tool_url) as response, open(key, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
            ucsc_tools[key] = os.getcwd() + "/" + key
            os.chmod(ucsc_tools[key], 0o777)
    if args.verbosity > 0:
        print("Will use " + ucsc_tools[key])


''' track color defintion '''

col_idx = 0
rgb_cols = ['0,0,0',       '255,0,0',   '0,255,0',     '0,0,255',     '176,196,222',
            '0,255,255',   '255,0,255', '192,192,192', '128,128,128', '128,0,0',
            '128,128,0',   '0,128,0',   '128,0,128',   '0,128,128',   '0,0,128',
            '139,0,0',     '220,20,60', '233,150,122', '255,140,0',   '218,165,32',
            '154,205,50',  '34,139,34', '152,251,152', '72,209,203',  '176,224,230', 
            '138,43,226']  # 0 - 25, number: 26

''' Create hub directory structure '''

try:
    os.makedirs(hub_dir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

try:
    os.makedirs(tmp_dir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


''' When UCSC tools are present, also check for the *.as files or obtain them;

store the as files either make_hub.py directory if it is writable, otherwise
store in temporary hub directory '''

as_url = 'http://genome.ucsc.edu/goldenPath/help/examples/'
for key, val in ucsc_as_files.items():
    if ((not args.no_genePredToBigGenePred) and (key == 'BigGenePred.as')) or (key != 'BigGenePred.as'):
        as_dir = os.path.dirname(os.path.realpath(__file__))
        if os.path.isfile(as_dir + '/' + key):
            ucsc_as_files[key] = as_dir + '/' + key
        elif os.path.isfile(tmp_dir + '/' + key):
            ucsc_as_files[key] = tmp_dir + '/' + key
        else:
            print("Was unable to locate file " + key + " in " + as_dir + " and " +
                  tmp_dir + ", will try to download it from " + as_url + "...")
            if (not os.access(as_dir, os.W_OK)) and os.access(tmp_dir, os.W_OK):
                as_dir = tmp_dir
            elif (not os.access(as_dir, os.W_OK)) and (not os.access(tmp_dir, os.W_OK)):
                print('Error: neither ' + as_dir + " nor " +
                      tmp_dir + " are writable directories.")
                exit(1)
            ucsc_as_files[key] = as_dir + '/' + key
            print("Trying to download " + as_url + key +
                  " and store into " + ucsc_as_files[key])
            try:
                with urllib.request.urlopen(as_url + key) as response, open(ucsc_as_files[key], 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
                print("Stored file at " + ucsc_as_files[key] + ".")
            except:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ": Downloading file " + key +
                      " failed. If your " +
                      "machine does not have internet access, you may " +
                      "manually download the file from " + as_url + key +
                      " and store it in the directory where make_hub.py " +
                      "resides.")


''' ******************* BEGIN FUNCTIONS *************************************'''


def set_color(c, col_lst):
    if c <= len(col_lst)-2:
        c = c+1
    else:
        c = 0
    return c, col_lst[c]


def set_visibility(v):
    if v < 10:
        print_vis = "dense"
        v = v + 1
    else:
        print_vis = "hide"
    return v, print_vis


''' Function that runs a subprocess with arguments '''


def run_simple_process(args_lst):
    try:
        # bed files need sorting with LC_COLLATE=C
        myenv = os.environ.copy()
        myenv['LC_COLLATE'] = 'C'
        if args.verbosity > 0:
            print("Trying to execute the following command:")
            print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=myenv)
        if args.verbosity > 0:
            print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " +
                  str(result.returncode) + str(result.args))
            quit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


''' Function that runs a subprocess with input from STDIN '''


def run_process_stdinput(args_lst, byte_obj):
    try:
        # bed files need sorting with LC_COLLATE=C
        myenv = os.environ.copy()
        myenv['LC_COLLATE'] = 'C'
        if args.verbosity > 0:
            print("Trying to execute the following command with input from " +
                  "STDIN:")
            print(" ".join(args_lst))
        result = subprocess.run(args_lst, input=byte_obj,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=myenv)
        if args.verbosity > 0:
            print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': '
                  + "run_process_stdinput: return code of subprocess was "
                  + str(result.returncode))
            quit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' +
              "Failed executing: ", " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


''' Function that sanity checks and fixes a gtf file, e.g. for strand problems '''


def make_gtf_sane(annot_file, ucsc_file):
    try:
        with open(annot_file, "r") as annot_handle:
            txs = {}
            for line in annot_handle:
                if ('transcript_id' in line):
                    seq, first_part, strand, second_part, txid = re.search(
                        r'(\S+)(\t\S+\t\S+\t\d+\t\d+\t\S+\t)(\S+)(\t\S+\t).*transcript_id\s\"(\S+)\"', line).groups()
                    if txid not in txs:
                        txs[txid] = {'strand': strand,
                                     'lines': [], 'sane': True}
                    else:
                        if not(strand == txs[txid]['strand']):
                            txs[txid]['sane'] = False
                    txs[txid]['lines'].append(line)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              annot_file + " for reading!")
        quit(1)
    try:
        with open(ucsc_file, "w") as ucsc_handle:
            for key, value in txs.items():
                if value['sane'] == True:
                    for line in value['lines']:
                        ucsc_handle.write(line)
                else:
                    for line in value['lines']:
                        seq, first_part, strand, second_part, txid = re.search(
                            r'^(\S+)(\t\S+\t\S+\t\d+\t\d+\t\S+\t)(\S+)(\t\S+\t).*transcript_id\s\"(\S+)\"', line).groups()
                        new_gid = txid + "_" + seq + "_" + strand
                        new_txid = new_gid + ".t1"
                        ucsc_handle.write(seq + first_part + strand + second_part +
                                          "gene_id \"" + new_gid + "\"; transcript_id \"" + new_txid + "\";\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              ucsc_file + " for writing!")
        quit(1)


''' Function that converts Gemoma filtered predictions in GFF3 format
    to GTF format that is similar to AUGUSTUS GTF format '''


def gemoma2aug_gtf(gemoma_file, aug_like_file):
    try:
        with open(aug_like_file, "w") as aug_handle:
            try:
                with open(gemoma_file, "r") as gemoma_handle:
                    for line in gemoma_handle:
                        line = line.strip()
                        if re.search(r'\tCDS\t', line) or re.search(r'\tfive_prime_UTR\t', line) or re.search(r'\tthree_prime_UTR\t', line):
                            f0, f1, f2, f3, f4, f5, f6, f7, f8 = line.split()
                            if re.match(r"ID=", f8):
                                gff3_part = f8.split(";")
                                tid_part = gff3_part[1].split("=")
                                gid = tid = tid_part[1]
                            else:
                                thismatch = re.search(r'Parent=([^;]+);?', f8)
                                gid = tid = thismatch.group(1)
                            if re.search(r'_R\d+', gid):
                                gid_groups = re.search(
                                    r'(^\S+)_R\d+', gid).groups()
                                gid = gid_groups[0]
                            if f2 == 'five_prime_UTR':
                                f2 = '5\'-UTR'
                            elif f2 == 'three_prime_UTR':
                                f2 = '3\'-UTR'
                            aug_handle.write(f0 + "\t" + f1 + "\t" + f2 + "\t" +
                                             f3 + "\t" + f4 + "\t" + f5 + "\t" +
                                             f6 + "\t" + f7 + "\t" +
                                             "gene_id \"" + gid + "\"; " +
                                             "transcript_id \"" + tid +
                                             "\";\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " + gemoma_file +
                      " for reading!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + aug_like_file +
              " for writing!")
        quit(1)


''' Function that converts GlimmerHMM predictions in GFF3 format to GTF format
    that is similar to AUGUSTUS GTF format '''


def glimmer2aug_gtf(glimmer_file, aug_like_file):
    try:
        with open(aug_like_file, "w") as aug_handle:
            try:
                with open(glimmer_file, "r") as glimmer_handle:
                    txid = None
                    for line in glimmer_handle:
                        line = line.strip()
                        if re.search(r'\tCDS\t', line):
                            f0, f1, f2, f3, f4, f5, f6, f7, f8 = line.split()
                            gff3_part = f8.split(";")
                            tid_part = gff3_part[1].split("=")
                            gid = tid_part[1]
                            aug_handle.write(f0 + "\t" + f1 + "\t" + f2 + "\t" +
                                             f3 + "\t" + f4 + "\t" + f5 + "\t" +
                                             f6 + "\t" + f7 + "\t" +
                                             "gene_id \"" + gid + ".1\"; " +
                                             "transcript_id \"" + gid +
                                             "\";\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " + glimmer_file +
                      " for reading!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + aug_like_file +
              " for writing!")
        quit(1)


''' Function that converts SNAP predictions in GFF3 format to GTF format
    that is similar to AUGUSTUS GTF format '''


def snap2aug_gtf(snap_file, aug_like_file):
    try:
        with open(aug_like_file, "w") as aug_handle:
            try:
                with open(snap_file, "r") as snap_handle:
                    txid = None
                    for line in snap_handle:
                        line = line.strip()
                        if re.search(r'\tCDS\t', line):
                            f0, f1, f2, f3, f4, f5, f6, f7 = line.split()
                            tid_part = f7.split("=")
                            gid = tid_part[1]
                            aug_handle.write(f0 + "\t" + f1 + "\t" + f2 + "\t" +
                                             f3 + "\t" + f4 + "\t" + f5 + "\t" +
                                             f6 + "\t.\t" +
                                             "gene_id \"" + gid + ".1\"; " +
                                             "transcript_id \"" + gid +
                                             "\";\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " + snap_file +
                      " for reading!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + aug_like_file +
              " for writing!")
        quit(1)


''' Function that reformats AUGUSTUS gtf format to UCSC gtf format '''


def aug2ucsc_gtf(augustus_file, ucsc_file):
    try:
        with open(augustus_file, "r") as aug_handle:
            try:
                with open(ucsc_file, "w") as ucsc_handle:
                    for line in aug_handle:
                        if re.search(r'\t\S+\tCDS\t', line) or re.search(r'\t\S+\texon\t', line) or re.search(r'\t\S+\tstart_codon\t', line) or re.search(r'\t\S+\t\d+\'-UTR\t', line):
                            if re.search(r'\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\ttranscript_id\s\"\S+\";\sgene_id\s\"\S+\";', line):
                                first_part, feature, second_part, gid, txid = re.search(
                                    r'(\S+\t\S+\t)(\S+)(\t\d+\t\d+\t\S+\t\S+\t\S+\t)(transcript_id\s\"\S+\";)\s(gene_id\s\"\S+\";)', line).groups()
                            elif re.search(r'\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\tgene_id\s\"\S+\";\stranscript_id\s\"\S+\";', line):
                                first_part, feature, second_part, txid, gid = re.search(
                                    r'(\S+\t\S+\t)(\S+)(\t\d+\t\d+\t\S+\t\S+\t\S+\t)(gene_id\s\"\S+\";)\s(transcript_id\s\"\S+\";)', line).groups()
                            if re.search(r'5\'-UTR', feature):
                                feature = "5UTR"
                            elif re.search(r'3\'-UTR', feature):
                                feature = "3UTR"
                            ucsc_handle.write(
                                first_part + feature + second_part + gid + " " +
                                txid + "\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      ucsc_file + " for writing!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + augustus_file +
              " for reading!")
        quit(1)


''' Function that writes subprocess byte object to flat file '''


def write_byteobj(byte_obj, outfile):
    try:
        with open(outfile, 'w') as byteobj_handle:
            byteobj_handle.write(byte_obj.decode('utf-8'))
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + outfile +
              " for writing!")
        quit(1)


''' Function that determines which regions in genome are softmasked '''


def find_masked_intervals(genome_file, bed3_file):
    try:
        masked_intervals = []
        with open(genome_file, "r") as genome_handle:
            for record in SeqIO.parse(genome_handle, "fasta"):
                masked_seq = str(record.seq)
                inMasked = False
                mmEnd = 0
                for x in range(0, len(masked_seq)):
                    if masked_seq[x].islower() and inMasked == False:
                        mmStart = x + 1
                        inMasked = True
                    elif not(masked_seq[x].islower()) and inMasked == True:
                        mmEnd = x
                        inMasked = False
                        if x > 0 and mmEnd > 0:
                            masked_intervals.append(
                                {'id': record.id, 'start': mmStart, 'end': mmEnd})
                if inMasked == True:
                    masked_intervals.append(
                        {'id': record.id, 'start': mmStart, 'end': len(masked_seq)})
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + genome_file +
              " for reading!")
        quit(1)

    try:
        with open(bed3_file, "w+") as bed2_handle:
            for x in masked_intervals:
                bed2_handle.write(
                    x['id'] + '\t' + str(x['start']) + '\t' + str(x['end']) +
                    '\n')
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + bed3_file +
              " for writing!")
        quit(1)


''' Function that sorts a bed3 file with LC_COLLATE=C via custom subprocess function'''


def sort_bed3(bed3_file, bed3_sorted_file):
    subprcs_args = [sort_tool, '-k1,1', '-k2,2n', bed3_file]
    result = run_simple_process(subprcs_args)
    write_byteobj(result.stdout, bed3_sorted_file)


''' Function that converts bed to bigBed '''


def bed2bigBed(btype, bed_file, chrom_size_file, bigBedFile, as_file):
    print('Generating bigBed file for ' + bed_file + '...')
    if as_file is None:
        subprcs_args = [ucsc_tools['bedToBigBed'], '-type=bed' +
                        str(btype), bed_file, chrom_size_file, bigBedFile]
    else:
        subprcs_args = [ucsc_tools['bedToBigBed'], '-type=bed' +
                        str(btype), bed_file, '-as=' + as_file, chrom_size_file, bigBedFile]
    run_simple_process(subprcs_args)


''' Function that converts bigGenePred bed-format to special bigBed fromat
for gene predictions '''


def bigGenePredToBigBed(as_file, bigGenePred_file, chrom_size_file, bb_file):
    print('Generating bigBed file from bigGenePred format for ' +
          bigGenePred_file + '...')
    subprcs_args = [ucsc_tools['bedToBigBed'], '-as=' + as_file,
                    '-type=bed12+8', bigGenePred_file, chrom_size_file,
                    bb_file]
    run_simple_process(subprcs_args)


''' Function that converts gtf to genePred format, used in gtf2bb and gtf2bgpbb '''


def gtf2genePred(gtf_file, gp_file, info_out_file, frameGiven):
    if frameGiven is True:
        subprcs_args = [ucsc_tools['gtfToGenePred'], '-infoOut=' +
                        info_out_file, '-genePredExt', gtf_file, gp_file]
    else:
        subprcs_args = [ucsc_tools['gtfToGenePred'], '-infoOut=' +
                        info_out_file, gtf_file, gp_file]
    run_simple_process(subprcs_args)
    subprcs_args = [ucsc_tools['genePredCheck'], gp_file]
    result = run_simple_process(subprcs_args)
    # parse result for failed annotations
    annotation_validation_result = result.stderr.decode('utf-8')
    if args.verbosity > 0:
        print(annotation_validation_result)
    regex_result = re.match(
        r'checked: \d+ failed: (\d+)', annotation_validation_result)
    if int(regex_result.group(1)) > 0:
        print("Warning: " + regex_result.group(1) +
              " annotations did not pass validation!")


''' Function that converts gtf to bb format '''


def gtf2bb(gtf_file, gp_file, bed_file, bb_file, info_out_file, chrom_size_file,
           sort_tool, frameGiven):
    gtf2genePred(gtf_file, gp_file, info_out_file, frameGiven)
    subprcs_args = [ucsc_tools['genePredToBed'], gp_file, 'stdout']
    result = run_simple_process(subprcs_args)
    subprcs_args = [sort_tool, '-k1,1', '-k2,2n']
    result = run_process_stdinput(subprcs_args, result.stdout)
    write_byteobj(result.stdout, bed_file)
    bed2bigBed(12, bed_file, chrom_size_file, bb_file, None)


''' Function that converts gtf to bigGenePred bb format (allows display of
amino acids in browser), alternative to gtf2bb '''


def gtf2bgpbb(gtf_file, gp_file, bigGenePred_file, info_out_file, chrom_size_file,
              as_file, sort_tool, bb_file, short_label):
    gtf2genePred(gtf_file, gp_file, info_out_file, True)
    subprcs_args = [ucsc_tools['genePredToBigGenePred'], gp_file, 'stdout']
    result = run_simple_process(subprcs_args)
    subprcs_args = [sort_tool, '-k1,1', '-k2,2n']
    result = run_process_stdinput(subprcs_args, result.stdout)
    write_byteobj(result.stdout, bigGenePred_file)
    # type and geneType field are empty in resulting file, fix it:
    bigGenePredFixed_file = tmp_dir + short_label + ".fixed.bed"
    try:
        with open(bigGenePred_file, "r") as bigGenePred_handle:
            try:
                with open(bigGenePredFixed_file, "w") as fixed_handle:
                    for line in bigGenePred_handle:
                        line_items = re.split(r'\t+', line)
                        for i in range(0, 16):
                            fixed_handle.write(line_items[i] + '\t')
                        fixed_handle.write(
                            'none\t' + line_items[16] + '\tnone\t' + line_items[17] + '\n')
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      bigGenePredFixed_file + " for writing!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              bigGenePred_file + " for reading!")
        quit(1)
    bigGenePredToBigBed(as_file, bigGenePredFixed_file,
                        chrom_size_file, bb_file)


''' Function that writes info about gene pred or hints to trackDb file '''


def info_to_trackDB(trackDb_file, short_label, long_label, rgb_color, group,
                    bed_no, visibility, ttype):
    # if maker track, use separate name for label...
    track_name = short_label
    if(re.search(r'_maker', track_name)):
        track_name = re.sub(r'_maker', r'', track_name)
    try:
        with open(trackDb_file, "a") as trackDb_handle:
            trackDb_handle.write("track " + short_label + "\n" +
                                 "longLabel " + long_label + "\n" +
                                 "shortLabel " + track_name + "\n" +
                                 "group " + group + "\ntype " + ttype)
            if bed_no is not None:
                trackDb_handle.write(" " + str(bed_no) + " .")
            trackDb_handle.write("\n" +
                                 "bigDataUrl " + short_label + ".bb\n" +
                                 "color " + rgb_color + "\n" +
                                 "visibility " + visibility +
                                 "\nhtml " + short_label + ".html\n\n" +
                                 "group " + group + "\ntype " + ttype)
            if bed_no is not None:
                trackDb_handle.write(" " + str(bed_no) + " .")
            trackDb_handle.write("\n" +
                                 "bigDataUrl " + short_label + ".bb\n" +
                                 "color " + rgb_color + "\n\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              trackDb_file + " for writing!")
        quit(1)


''' Function that converts gtf to gene prediction track '''


def make_gtf_track(trackDb_file, gtf_file, chrom_size_file, short_label, long_label, rgb_color, visibility, no_genePredToBigGenePredFlag, frameGiven):
    gp_file = tmp_dir + short_label + ".gp"
    info_out_file = tmp_dir + short_label + ".infoOut.txt"
    bed_file = tmp_dir + short_label + ".bed"
    bb_file = hub_dir + short_label + ".bb"
    if (no_genePredToBigGenePredFlag is True) or (frameGiven is False):
        gtf2bb(gtf_file, gp_file, bed_file, bb_file,
               info_out_file, chrom_size_file, sort_tool, frameGiven)
        info_to_trackDB(trackDb_file, short_label, long_label,
                        rgb_color, "genePreds", 12, visibility, 'bigBed')
    else:
        gtf2bgpbb(gtf_file, gp_file, bed_file, info_out_file, chrom_size_file,
                  ucsc_as_files['bigGenePred.as'], sort_tool, bb_file, short_label)
        print("Next step should be writing trackDb entry...")
        info_to_trackDB(trackDb_file, short_label, long_label,
                        rgb_color, "genePreds", None, visibility, 'bigGenePred')
    # parse info_out_file to produce txt file for creating nameIndex files
    name_index_txt_file = tmp_dir + short_label + ".nameIndex.txt"
    try:
        with open(name_index_txt_file, "w") as name_index_handle:
            try:
                with open(info_out_file, "r") as info_out_handle:
                    for line in info_out_handle:
                        if not(re.match(r'^\#', line)) and (re.match(r'(\S+)\s+(\S+)\s+(\S+)\s+', line)):
                            tx_id, g_id, src = re.match(
                                r'(\S+)\s+(\S+)\s+(\S+)\s+', line).groups()
                            name_index_handle.write(
                                tx_id + "\t" + g_id + "," + src + ",,,\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      info_out_file + " for reading!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              name_index_txt_file + " for writing!")
        quit(1)
    ix_file = hub_dir + short_label + ".nameIndex.ix"
    ixx_file = hub_dir + short_label + ".nameIndex.ixx"
    subprcs_args = [ucsc_tools['ixIxx'],
                    name_index_txt_file, ix_file, ixx_file]
    run_simple_process(subprcs_args)



''' Function that cleans wiggle files '''


def clean_wig(tmp_wig_file, wig_file, chrom_size_file):
    # observed that in rare cases a wig file might contain coverage for one
    # more base than present in the sequence; seems to be an alignment
    # error, not a bam2wig error, because the same problem arises if I
    # convert bam to wig in python in a different way
    # therefore check wiggle file for sanity and modify if required
    # (i.e. cleave coverage at sequence end)
    if args.verbosity > 0:
        print("Sanity checking wig file...")
    chrom_sizes = {}
    try:
        with open(chrom_size_file, "r") as size_handle:
            for line in size_handle:
                split_list = re.split(r'\t', line.rstrip('\n'))
                chrom_sizes[split_list[0]] = int(split_list[1])
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Could not open file " +
              chrom_size_file + " for reading!")
    try:
        with open(tmp_wig_file, "r") as tmp_handle:
            try:
                with open(wig_file, "w") as wig_handle:
                    seqname = ""
                    for line in tmp_handle:
                        match = re.search(r'chrom=(\S+)', line)
                        if match:
                            seqname = match.group(1)
                        match = re.search(r'(\d+) \d+', line)
                        if match:
                            if int(match.group(1)) <= chrom_sizes[seqname]:
                                wig_handle.write(line)
                        else:
                            wig_handle.write(line)
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Could not open file " +
                      wig_file + " for writing!")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Could not open file " +
              tmp_wig_file + " for reading!")


''' Function that converts bam file to wig file using AUGUSTUS auxprog bam2wig'''


def bam2wig(bam_file, wig_file, size_file):
    tmp_wig_file = wig_file + ".tmp"
    try:
        with open(tmp_wig_file, "w") as wig_handle:
            subprcs_args = [augustus_tools['bam2wig'],
                            "-t", bam_file, bam_file]
            if args.verbosity > 0:
                print("Trying to execute the following command:")
                print(" ".join(subprcs_args))
            result = subprocess.run(
                subprcs_args, stdout=wig_handle, stderr=subprocess.PIPE)
            if args.verbosity > 0:
                print("Suceeded in executing command.")
            if(result.returncode != 0):
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' +
                      "Return code of subprocess was " + str(result.returncode))
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Could not open file " +
              tmp_wig_file + " for writing!")
    clean_wig(tmp_wig_file, wig_file, size_file)


''' Function that converts bam file to wig file with RSeQC '''


def bamToWig(bam_file, wig_file, size_file):
    pileup_file = wig_file + ".mpileup"
    subprcs_args = [samtools, "mpileup", "-o", pileup_file, bam_file]
    run_simple_process(subprcs_args)
    # observed that in rare cases a wig file might contain coverage for one
    # more base than present in the sequence; seems to be an alignment
    # error, not a bam2wig error, because the same problem arises if I
    # convert bam to wig in python in a different way
    # therefore check wiggle file for sanity and modify if required
    # (i.e. cleave coverage at sequence end)
    if args.verbosity > 0:
        print("Reading chrom_sizes for wig sanity check...")
    chrom_sizes = {}
    try:
        with open(size_file, "r") as size_handle:
            for line in size_handle:
                split_list = re.split(r'\t', line.rstrip('\n'))
                chrom_sizes[split_list[0]] = int(split_list[1])
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Could not open file " +
              size_file + " for reading!")
    try:
        with open(pileup_file, "r") as pileup_handle:
            try:
                with open(wig_file, "w") as wig_handle:
                    wig_handle.write(
                        "track name=" + bam_file + " type=wiggle_0\n")
                    lastSeq = ""
                    for line in pileup_handle:
                        seq, start, t1, depth, t2, t3 = line.split()
                        if (seq != lastSeq):
                            wig_handle.write(
                                "variableStep chrom=" + seq + "\n")
                        if(int(start) <= chrom_sizes[seq]):
                            wig_handle.write(start + "\t" + depth + "\n")
                        lastSeq = seq
            except IOError:
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' +
                      "Could not open file " + wig_file + " for writing!")
    except IOError:
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' +
              "Could not open file " + pilup_file + " for reading!")


''' Function that checks whether a track already exists in trackDb '''


def check_trackDb(trackDb, label):
    if os.path.isfile(trackDb):
        my_regex = r"" + re.escape(label) + r""
        try:
            with open(trackDb, "r") as trackDb_handle:
                for line in trackDb_handle:
                    if re.search(my_regex, line):
                        frameinfo = getframeinfo(currentframe())
                        print('Error in file ' + frameinfo.filename + ' at line ' +
                              str(frameinfo.lineno) + ': ' +
                              "the label \"" + label + "\" already exists in " +
                              "trackDb_file " + trackDb + ". Gene prediction " +
                              "tracks with pre-defined labels can only be added " +
                              "to an assembly hub, once.\n\nEither create a new hub\n\n" +
                              "or manually edit the file " + args.short_label +
                              "/" + args.short_label + "/trackDb.txt to rename " +
                              "the already existing track with that label and " +
                              "try re-running make_hub.py in --add_track mode, " +
                              "subsequently. \n\nRemember, that in case of manual " +
                              "editing the trackDb.txt file, not only do you " +
                              "have to change the label \"" + label + "\", but " +
                              "you also have to change the string after \"track \"" +
                              "above that term to be unique after adding another " +
                              "track with that same label.")
                        exit(1)
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  trackDb_file + " for reading!")


''' Function that converts MAKER GFF3 to Bed '''


def maker_to_bed(line_lst, outfile, chrom_size_file):
    chrom_sizes = {}
    try:
        with open(chrom_size_file, "r") as size_handle:
            for line in size_handle:
                split_list = re.split(r'\t', line.rstrip('\n'))
                chrom_sizes[split_list[0]] = int(split_list[1])
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Could not open file " +
              chrom_size_file + " for reading!")
    bed_items = {}
    for line in line_lst:
        if re.search(r'\tmatch_part\t', line):
            seq, start, end, score, strand, parent = re.search(
                r'(\S+)\t\S+\tmatch_part\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t\S+\t(\S+)', line).groups()
            parent = re.sub(r'\S+Parent=', r'', parent)
            parent = re.sub(r';Target=\S+', r'', parent)
            if not parent in bed_items:
                bed_items[parent] = {}
                bed_items[parent]['seq'] = seq
                bed_items[parent]['chromStart'] = int(start)
                # + 1 because the segment is drawn until the next position
                bed_items[parent]['chromEnd'] = int(end)
                bed_items[parent]['name'] = parent
                bed_items[parent]['score'] = 0
                bed_items[parent]['strand'] = strand
                bed_items[parent]['itemRgb'] = 0
                bed_items[parent]['blockCount'] = 1
                bed_items[parent]['blockSizes'] = []
                bed_items[parent]['blockSizes'].append(
                    (int(end) - int(start)) + 1)
                bed_items[parent]['blockStarts'] = []
                bed_items[parent]['blockStarts'].append(int(start))
                bed_items[parent]['lines'] = []
                # store lines in case we later need to resolve overlap issues
                bed_items[parent]['lines'].append(line)
            else:
                if int(start) <= bed_items[parent]['chromStart']:
                    bed_items[parent]['chromStart'] = int(start)
                if int(end) >= bed_items[parent]['chromEnd']:
                    bed_items[parent]['chromEnd'] = int(end)
                bed_items[parent]['blockCount'] = bed_items[parent]['blockCount'] + 1
                bed_items[parent]['blockSizes'].append(
                    (int(end) - int(start)) + 1)
                bed_items[parent]['blockStarts'].append(int(start))
                bed_items[parent]['lines'].append(line)
    # check whether blocks of match_part overlap in one parent and separate
    # safe to assume that they are sorted by MAKER
    # typical case are BLAST hits that often overlap
    new_items = {}
    for parent in bed_items:
        overlap = False
        bed_items[parent]['blockStarts'][:] = [x - bed_items[parent]
                                               ['chromStart'] for x in bed_items[parent]['blockStarts']]
        if bed_items[parent]['strand'] == '-':
            bed_items[parent]['blockSizes'].reverse()
            bed_items[parent]['blockStarts'].reverse()
        for i in range(0, len(bed_items[parent]['blockStarts']) - 1):
            if (bed_items[parent]['blockStarts'][i] + bed_items[parent]['blockSizes'][i]) >= bed_items[parent]['blockStarts'][i + 1]:
                overlap = True
        if overlap is True:
            c = 1
            for line in bed_items[parent]['lines']:
                seq, start, end, score, strand, p = re.search(
                    r'(\S+)\t\S+\tmatch_part\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t\S+\t(\S+)', line).groups()
                p = re.sub(r'\S+Parent=', r'', p)
                p = re.sub(r';Target=\S+', r'', p)
                new_parent = parent + "_" + str(c)
                new_items[new_parent] = {}
                new_items[new_parent]['seq'] = seq
                new_items[new_parent]['chromStart'] = int(start)
                # + 1 because the segment is drawn until the next position, if we don't do it, bed format is inconsistent in many lines?!
                # note all lines, fix this, again, below
                new_items[new_parent]['chromEnd'] = int(end) + 1
                new_items[new_parent]['name'] = new_parent
                new_items[new_parent]['score'] = 0
                new_items[new_parent]['strand'] = strand
                new_items[new_parent]['itemRgb'] = 0
                new_items[new_parent]['blockCount'] = 1
                new_items[new_parent]['blockSizes'] = []
                new_items[new_parent]['blockSizes'].append(
                    (int(end) - int(start)) + 1)
                new_items[new_parent]['blockStarts'] = []
                new_items[new_parent]['blockStarts'].append(
                    int(start) - new_items[new_parent]['chromStart'])
                c = c + 1
        else:
            new_items[parent] = bed_items[parent]
    bed_items = new_items
    del new_items
    # check whether items don't exceed lengths of sequences:
    for parent in bed_items:
        # I am not sure why sometimes the bed format is inconsistent, usually by one, fix it:
        if (bed_items[parent]['chromStart'] + bed_items[parent]['blockStarts'][-1] + bed_items[parent]['blockSizes'][-1]) != bed_items[parent]['chromEnd']:
            bed_items[parent]['chromEnd'] = bed_items[parent]['chromStart'] + \
                bed_items[parent]['blockStarts'][-1] + \
                bed_items[parent]['blockSizes'][-1]
        if bed_items[parent]['chromEnd'] > chrom_sizes[bed_items[parent]['seq']]:
            the_diff = bed_items[parent]['chromEnd'] - \
                chrom_sizes[bed_items[parent]['seq']]
            bed_items[parent]['chromEnd'] = bed_items[parent]['chromEnd'] - the_diff
            bed_items[parent]['blockSizes'][-1] = bed_items[parent]['blockSizes'][-1] - the_diff
    # write bed file
    try:
        with open(outfile, "w") as out_handle:
            for parent in bed_items:
                out_handle.write(bed_items[parent]['seq'] + "\t" +
                                 str(bed_items[parent]['chromStart']) + "\t" +
                                 str(bed_items[parent]['chromEnd']) + "\t" +
                                 bed_items[parent]['name'] + "\t" +
                                 str(bed_items[parent]['score']) + "\t" +
                                 bed_items[parent]['strand'] + "\t" +
                                 str(bed_items[parent]['chromStart']) + "\t" +
                                 str(bed_items[parent]['chromStart']) + "\t" +
                                 str(bed_items[parent]['itemRgb']) + "\t" +
                                 str(bed_items[parent]['blockCount']) + "\t")

                out_handle.write(','.join(str(x)
                                          for x in bed_items[parent]['blockSizes']))
                out_handle.write('\t')

                out_handle.write(','.join(str(x)
                                          for x in bed_items[parent]['blockStarts']))
                out_handle.write("\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              outfile + " for writing!")


''' Function that converts MAKER GFF3 to GTF '''


def maker_to_gtf(line_lst, outfile):
    # establish tx->gid relationship
    gids = {}
    for line in line_lst:
        if re.search(r'\tmRNA\t', line):
            f0, f1, f2, f3, f4, f5, f6, f7, f8 = line.split()
            gff3_part = f8.split(";")
            parent = re.search(r'Parent=(\S+)', gff3_part[1]).groups()
            child = re.search(r'ID=(\S+)', gff3_part[0]).groups()
            if child[0] not in gids:
                gids[child[0]] = parent[0]
    try:
        with open(outfile, "w") as aug_handle:
            for line in line_lst:
                line = line.strip()
                if re.search(r'\tCDS\t', line) or re.search(r'\tfive_prime_UTR\t', line) or re.search(r'\tthree_prime_UTR\t', line):
                    f0, f1, f2, f3, f4, f5, f6, f7, f8 = line.split()
                    gff3_part = f8.split(";")
                    tid_part = gff3_part[1].split("=")
                    tx = tid_part[1]
                    gid = gids[tx]
                    if f2 == 'five_prime_UTR':
                        f2 = '5\'-UTR'
                    elif f2 == 'three_prime_UTR':
                        f2 = '3\'-UTR'
                    aug_handle.write(f0 + "\t" + f1 + "\t" + f2 + "\t" +
                                     f3 + "\t" + f4 + "\t" + f5 + "\t" +
                                     f6 + "\t" + f7 + "\t" +
                                     "gene_id \"" + gid + "\"; " +
                                     "transcript_id \"" + tx +
                                     "\";\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + outfile +
              " for writing!")
        quit(1)


''' Function that writes aboutHub.html file '''


def write_aboutHub(outfile, title, latin, assembly, email):
    try:
        with open(outfile, "w") as handle:
            handle.write(
                "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd\">\n")
            handle.write("<html>\n<head>\n<title>" + title +
                         "</title>\n</head>\n<body>\n<hr>\n")
            handle.write("<h1>" + title +
                         "</h1>\n<hr><p>This is an Assembly Hub for display with the UCSC Genome Browser that was automatically generated for the single species " +
                         title + " (" + latin + "), assembly version " +
                         assembly + " with make_hub.py version " + __version__ +
                         ".</p>")
            handle.write("<p>Contact: " + email + "</p>")
            handle.write("<hr>\n</body>\n</html>\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + outfile +
              " for writing!")
        quit(1)


''' Function that writes track description page '''


def write_trackDescription(outfile, track_name, short_method, email, track_color):
    try:
        with open(outfile, "w") as handle:
            handle.write(
                "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd\">\n")
            handle.write("<html>\n<head>\n<title>" + track_name +
                         "</title>\n</head>\n<body>\n")
            handle.write("<h2>Description</h2>\n" +
                         "The track " + track_name + " was automatically generated using <a href=\"https://github.com/Gaius-Augustus/MakeHub\">make_hub.py</a>. " +
                         "We strongly recommend that the creator of this hub modifies this page to be more informative!</p>" +
                         "<h2>Display Conventions and Configuration</h2>\n" +
                         "<p>This section describes track configuration controls, or any special display conventions such as the meaning of different colors in your tracks.</p>" +
                         "<p>By default, all items in this track have the RGB color " +
                         track_color + ".</p>\n" +
                         "<h2>Methods</h2>\n" +
                         "<p>This section describes the methods used to generate and analyze the data and helps users understand how the track data was produced and sometimes has subsections if useful.</p>" +
                         "<p>The default methods description for this track provided by make_hub.py is:</p>\n" +
                         "<p>" + short_method + "</p>\n" +
                         "<h2>Credits</h2>\n" +
                         "<p>Please list invidividuals and/or organizations who contributed to the collection and analysis of the data. Be sure to include a preferred contact email address for users who have questions concerning the data.</p>" +
                         "<p>Please direct questions about this assembly hub to " + email + ".</p>\n" +
                         "<h2>References</h2>" +
                         "<p>Hoff KJ (2019) <a href=\"https://github.com/Gaius-Augustus/MakeHub\">MakeHub</a> (manuscript in preparation)</p>\n" +
                         "</body>\n</html>\n")

    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + outfile +
              " for writing!")
        quit(1)


''' Function that determines index of bam entries in trackDb file (if present) 
    will use index of bw files, not bam files, in case both are present '''


def find_bam_index(trackDb):
    bam_index = 1
    if os.path.isfile(trackDb):
        try:
            with open(trackDb, "r") as db_handle:
                for line in db_handle:
                    if re.search(r'rnaseq_\d+.bw', line):
                        re_result = re.search(
                            r'rnaseq_(\d+).bw', line).groups()
                        bam_index = int(re_result[0]) + 1
                        pass
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  trackDb + " for reading!")
            quit(1)
    return bam_index


''' ******************* END FUNCTIONS ***************************************'''

''' Globally required files '''

ChromSizes_file = hub_dir + args.short_label + \
    ".chrom.sizes"  # for making tracks
trackDb_file = hub_dir + "trackDb.txt"  # main UCSC hub configuration file


''' Generate essential files for genome display '''

visibility_counter = 0  # not more than 10 tracks should be unhidden

if not args.add_track:
    TwoBit_file = hub_dir + args.short_label + ".2bit"
    print('Generating genome 2bit file ' + TwoBit_file + '...')
    subprcs_args = [ucsc_tools['faToTwoBit'], args.genome, TwoBit_file]
    run_simple_process(subprcs_args)

    print('Generating chromsome size info file ' + ChromSizes_file + '...')
    subprcs_args = [ucsc_tools['twoBitInfo'], TwoBit_file, 'stdout']
    result = run_simple_process(subprcs_args)
    subprcs_args = [sort_tool, '-k2rn']
    result = run_process_stdinput(subprcs_args, result.stdout)
    write_byteobj(result.stdout, ChromSizes_file)

    WigVarStep_file = tmp_dir + args.short_label + ".gc5Base.wigVarStep"
    WigVarStep_file_compr = WigVarStep_file + ".gz"
    print('Generating variable step wiggle file for GC content ' +
          WigVarStep_file_compr + '...')
    subprcs_args = [ucsc_tools['hgGcPercent'], '-wigOut', '-doGaps',
                    '-file=stdout', '-win=5', '-verbose=0', args.short_label,
                    TwoBit_file]
    result = run_simple_process(subprcs_args)
    write_byteobj(result.stdout, WigVarStep_file)
    if os.path.isfile(WigVarStep_file_compr):
        os.unlink(WigVarStep_file_compr)
    subprcs_args = [gzip_tool, WigVarStep_file]
    run_simple_process(subprcs_args)

    BigWigGC_file = hub_dir + args.short_label + ".gc5Base.bw"
    print('Generating bigWig file for GC content ' + BigWigGC_file + '...')
    subprcs_args = [ucsc_tools['wigToBigWig'],
                    WigVarStep_file_compr, ChromSizes_file, BigWigGC_file]
    run_simple_process(subprcs_args)

    hub_txt_file = args.outdir + "/" + args.short_label + "/hub.txt"
    try:
        with open(hub_txt_file, "w+") as hub_txt_handle:
            hub_txt_handle.write("hub " + args.short_label + "\n")
            hub_txt_handle.write("shortLabel " + args.short_label + "\n")
            hub_txt_handle.write("longLabel " + args.long_label + "\n")
            hub_txt_handle.write("genomesFile genomes.txt\n")
            hub_txt_handle.write("email " + args.email + "\n")
            hub_txt_handle.write("descriptionUrl aboutHub.html\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + hub_txt_file +
              " for writing!")
        quit(1)
    about_file = args.outdir + "/" + args.short_label + "/aboutHub.html"
    ass_vers = ""
    if args.assembly_version is not None:
        ass_vers = args.assembly_version
    if args.latin_name is None:
        write_aboutHub(about_file, args.long_label,
                       args.long_label, ass_vers, args.email)
    else:
        write_aboutHub(about_file, args.long_label,
                       args.latin_name, ass_vers, args.email)

    default_seq_id = ""
    default_seq_end = 0
    try:
        with open(args.genome, "r") as genome_handle:
            nSeq = 0
            for record in SeqIO.parse(genome_handle, "fasta"):
                default_seq_id = record.id
                if len(record.seq) > 15000:
                    default_seq_end = 15000
                else:
                    default_seq_end = len(record.seq)
                nSeq = nSeq + 1
                if nSeq > 0:
                    continue
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + args.genome +
              " for reading!")
        quit(1)

    genomes_txt_file = args.outdir + "/" + args.short_label + "/genomes.txt"
    try:
        with open(genomes_txt_file, "w+") as genomes_txt_handle:
            genomes_txt_handle.write("genome " + args.short_label + "\n")
            genomes_txt_handle.write(
                "trackDb " + args.short_label + "/trackDb.txt\n")
            genomes_txt_handle.write(
                "groups " + args.short_label + "/groups.txt\n")
            genomes_txt_handle.write(
                "description ")
            if args.latin_name is not None:
                genomes_txt_handle.write("(" + args.latin_name + ")")
            genomes_txt_handle.write("[generated by MakeHub")
            if args.assembly_version is not None:
                genomes_txt_handle.write(" with " + args.assembly_version)
            genomes_txt_handle.write("]\n")
            genomes_txt_handle.write(
                "twoBitPath " + args.short_label + "/" + args.short_label +
                ".2bit\n")
            genomes_txt_handle.write("organism " + args.long_label + "\n")
            genomes_txt_handle.write(
                "defaultPos " + default_seq_id + ":1-" + str(default_seq_end) +
                "\n")
            genomes_txt_handle.write("orderKey 4800\n")
            genomes_txt_handle.write(
                "scientificName " + args.long_label + "\n")
            genomes_txt_handle.write("htmlPath hmi/description.html\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              genomes_txt_file + " for writing!")
        quit(1)

    try:
        with open(trackDb_file, "w+") as trackDb_handle:
            visibility_counter, visibility = set_visibility(visibility_counter)
            trackDb_handle.write("track gcPercent\nlongLabel GC Percent in 5-base " +
                                 "Window\nshortLabel GC Percent\n" +
                                 "type bigWig 0 100\ngroup map\nvisibility " + visibility +
                                 "\nwindowingFunction Mean\nbigDataUrl " +
                                 args.short_label + ".gc5Base.bw\npriority 2\nautoScale Off\n" +
                                 "maxHeightPixels 128:36:16\ngraphTypeDefault Bar\ngridDefault OFF\n" +
                                 "ncolor 0,0,0\naltColor 128,128,128\nviewLimits 30:70\nhtml gcPercent.html\n\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              trackDb_file + " for writing!")
        quit(1)
    html_file = hub_dir + "gcPercent.html"
    write_trackDescription(html_file, "gcPercent", "This GC-content track " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, "128,128,128")
    groups_txt_file = hub_dir + "groups.txt"
    try:
        with open(groups_txt_file, "w+") as groups_handle:
            groups_handle.write(
                "name genePreds\nlabel Gene Predictions\npriority 2\ndefaultIsClosed 0\n\n")
            groups_handle.write(
                "name reps\nlabel Repeats\npriority 2\ndefaultIsClosed 0\n\n")
            groups_handle.write(
                "name hints\nlabel Hints\npriority 2\ndefaultIsClosed 0\n\n")
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              groups_txt_file + " for writing!")
        quit(1)

    ''' Generate cytoband track '''

    if ucsc_as_files['cytoBand.as'] is not None:
        print('Generating cytoband track...')
        cytoBand_sort_file = tmp_dir + args.short_label + ".cytoBand.sorted"
        cytoBand_bed = tmp_dir + args.short_label + ".cytoband.bed"
        cytoBand_file = hub_dir + "cytoBandIdeo.bb"
        try:
            with open(ChromSizes_file, "r") as chrom_sizes_handle:
                chrom_sizes_data = chrom_sizes_handle.read()
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  ChromSizes_file + " for reading!")
            quit(1)
        subprcs_args = [sort_tool, '-k1,1', '-k2,2n']
        result = run_process_stdinput(
            subprcs_args, chrom_sizes_data.encode('utf-8'))
        write_byteobj(result.stdout, cytoBand_sort_file)
        try:
            with open(cytoBand_bed, "w") as cyto_bed_handle:
                try:
                    with open(cytoBand_sort_file, "r") as cyto_sort_handle:
                        for line in cyto_sort_handle:
                            line = line.strip()
                            line_entries = re.split(r'\t', line)
                            cyto_bed_handle.write(
                                line_entries[0] + "\t0\t" + line_entries[1] + "\t" +
                                line_entries[1] + "\tgneg\n")
                except IOError:
                    frameinfo = getframeinfo(currentframe())
                    print('Error in file ' + frameinfo.filename + ' at line ' +
                          str(frameinfo.lineno) + ': ' + "Failed to open file " +
                          cytoBand_sort_file + " for reading!")
                    quit(1)
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  cytoBand_bed + " for writing!")
            quit(1)
        bed2bigBed(4, cytoBand_bed, ChromSizes_file,
                   cytoBand_file, ucsc_as_files['cytoBand.as'])
        try:
            with open(trackDb_file, "a+") as trackDb_handle:
                trackDb_handle.write("track cytoBandIdeo\nlongLabel Chromosome " +
                                     "ideogram\nshortLabel cytoBandIdeo\n" + 
                                     "bigDataUrl cytoBandIdeo.bb\ntype bigBed\n\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                trackDb_file + " for writing!")
            quit(1)

    ''' Generate repeat masking track '''

    if not args.no_repeats:
        col_idx, this_color = set_color(col_idx, rgb_cols)
        softmaskedBed3_file = tmp_dir + args.short_label + ".RMsoft.bed3"
        print('Generating softmasking information bed3 file ' +
              softmaskedBed3_file + " from genome data (this may take a while)...")
        find_masked_intervals(args.genome, softmaskedBed3_file)
        softmaskedBed3_sorted_file = tmp_dir + args.short_label + ".RMsoft.s.bed3"
        print('Sorting file ' + softmaskedBed3_file + '...')
        sort_bed3(softmaskedBed3_file, softmaskedBed3_sorted_file)
        softmaskedBigBed_file = hub_dir + args.short_label + "_softmasking.bb"
        bed2bigBed(3, softmaskedBed3_sorted_file,
                   ChromSizes_file, softmaskedBigBed_file, None)
        try:
            with open(trackDb_file, "a") as trackDb_handle:
                visibility_counter, visibility = set_visibility(
                    visibility_counter)
                trackDb_handle.write("track RMsoft\nlongLabel Softmasked Repeats\n" +
                                     "shortLabel Repeats\ngroup reps\ntype bigBed 3 .\n" +
                                     "bigDataUrl " + args.short_label +
                                     "_softmasking.bb\nvisibility " + visibility + "\n" +
                                     "color " + this_color + "\n" +
                                     "html RMsoft.html" +
                                     "\n\ngroup reps\ntype bigBed 3 .\n" +
                                     "bigDataUrl " + args.short_label +
                                     "_softmasking.bb\n" +
                                     "color " + this_color + "\n\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  trackDb_file + " for writing!")
            quit(1)
        html_file = hub_dir + "RMsoft.html"
        write_trackDescription(html_file, "RMsoft", "This repeat masking track " +
                               "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                               "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>." +
                               "Repeat information was automatically extracted from softmasked genome file.</p>",
                               args.email, this_color)


''' check how many tracks are already visible in existing hub '''

if args.add_track:
    try:
        with open(trackDb_file, "r") as trackDb_handle:
            for line in trackDb_handle:
                if re.search(r'visibility', line):
                    if not(re.search(r'hide', line)):
                        visibility_counter = visibility_counter + 1
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              trackDb_file + " for reading!")
        quit(1)


''' Creating RNA-Seq bam track(s) '''

if args.bam and args.display_bam_as_bam:
    print('Generating BAM track(s)...')
    bam_index = find_bam_index(trackDb_file)
    for bam_file in args.bam:
        bam_sorted_file = hub_dir + "rnaseq_" + str(bam_index) + \
            ".s.bam"
        subprcs_args = [samtools, "sort", "-@",
                        str(args.cores), bam_file, "-o", bam_sorted_file]
        run_simple_process(subprcs_args)
        bam_index_file = hub_dir + "rnaseq_" + str(bam_index) + ".s.bam.bai"
        subprcs_args = [samtools, "index", "-@", str(args.cores),
                        bam_sorted_file, bam_index_file]
        run_simple_process(subprcs_args)
        try:
            with open(trackDb_file, "a") as trackDb_handle:
                visibility_counter, visibility = set_visibility(
                    visibility_counter)
                trackDb_handle.write("track RNASeq_" + str(bam_index) + "\n" +
                                     "bigDataUrl rnaseq_" + str(bam_index) +
                                     ".s.bam\n" +
                                     "shortLabel RNASeq_" + str(bam_index) +
                                     "\n" +
                                     "longLabel RNASeq bam file "
                                     + str(bam_index) + " from file " +
                                     bam_file +
                                     "\ngroup hints\nvisibility ")
                trackDb_handle.write(visibility)
                trackDb_handle.write("\ntype bam\nhtml " + "RNASeq_" +
                                     str(bam_index) + ".html\n\n" +
                                     "group hints\nbigDataUrl rnaseq_" +
                                     str(bam_index) + ".s.bam\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  trackDb_file + " for writing!")
            quit(1)
        html_file = hub_dir + "RNASeq_" + str(bam_index) + ".html"
        write_trackDescription(html_file, "RNASeq_" + str(bam_index), "This RNA-Seq BAM track " +
                               "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                               "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                               args.email, "None")
        bam_index = bam_index + 1


''' Creating bigWig RNA-Seq track(s) from bam '''

if args.bam:
    print('Generating bigWig RNA-Seq track(s) from BAM...')
    bam_index = find_bam_index(trackDb_file)
    for bam_file in args.bam:
        bam_sorted_file = hub_dir + "rnaseq_" + str(bam_index) + ".s.bam"
        if not args.display_bam_as_bam:
            subprcs_args = [samtools, "sort", "-@",
                            str(args.cores), bam_file, "-o", bam_sorted_file]
            run_simple_process(subprcs_args)
        wig_file = tmp_dir + "rnaseq_" + str(bam_index) + ".wig"
        if bam2wig_aug is True:
            bam2wig(bam_sorted_file, wig_file, ChromSizes_file)
        else:
            bamToWig(bam_sorted_file, wig_file, ChromSizes_file)
        wig_compr_file = tmp_dir + "rnaseq_" + str(bam_index) + ".wig.gz"
        if os.path.isfile(wig_compr_file):
            os.unlink(wig_compr_file)
        subprcs_args = [gzip_tool, wig_file]
        run_simple_process(subprcs_args)
        print("Converting compressed wig to bigWig")
        big_wig_file = hub_dir + "rnaseq_" + str(bam_index) + ".bw"
        subprcs_args = [ucsc_tools['wigToBigWig'],
                        wig_compr_file, ChromSizes_file, big_wig_file]
        run_simple_process(subprcs_args)
        print("Writing to track db")
        col_idx, this_color = set_color(col_idx, rgb_cols)
        try:
            with open(trackDb_file, "a") as trackDb_handle:
                visibility_counter, visibility = set_visibility(
                    visibility_counter)
                trackDb_handle.write("track RNASeq_wig_" + str(bam_index)
                                     + "\n" + "type bigWig\n" +
                                     "bigDataUrl rnaseq_" + str(bam_index) +
                                     ".bw\n" +
                                     "shortLabel RNASeq_" + str(bam_index) +
                                     "\n" +
                                     "longLabel RNASeq Wiggle " +
                                     str(bam_index) + " from bam file " +
                                     bam_file + "\ncolor " + this_color +
                                     "\nyLineOnOff on\nyLineMark 0\ngridDefault on\nvisibility ")
                trackDb_handle.write(
                    visibility + "\nhtml " + "RNASeq_wig_" + str(bam_index) + ".html")
                trackDb_handle.write("\n\ngroup hints\ntype bigWig\nbigDataUrl rnaseq_" +
                                     str(bam_index) + ".bw\ncolor " + this_color + "\n\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  trackDb_file + " for writing!")
            quit(1)
        html_file = hub_dir + "RNASeq_wig_" + str(bam_index) + ".html"
        write_trackDescription(html_file, "RNASeq_wig_" + str(bam_index), "This RNA-Seq WIG track " +
                               "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                               "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                               args.email, this_color)
        bam_index = bam_index + 1


''' Creating reference annotation gene track '''

if args.annot:
    print('Generating reference annotation gene track...')
    check_trackDb(trackDb_file, "Reference Annotation")
    ucsc_file = tmp_dir + "annot_ucsc.gtf"
    make_gtf_sane(args.annot, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "Annot",
                   "Reference Annotation", this_color, visibility,
                   args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "Annot.html"
    write_trackDescription(html_file, "annot", "This gene prediction track with reference annotation " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating genemark prediction track '''

if args.genemark:
    print('Generating GeneMark prediction track...')
    check_trackDb(trackDb_file, "GeneMark predicitons")
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, args.genemark, ChromSizes_file, "Genemark",
                   "GeneMark predictions", this_color, visibility,
                   args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "Genemark.html"
    write_trackDescription(html_file, "Genemark", "This gene prediction track with GeneMark-ES/ET predictions " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating AUGUSTUS ab initio track '''

if args.aug_ab_initio:
    print('Generating AUGUSTUS ab initio prediction (no UTRs) track...')
    check_trackDb(trackDb_file, "AUGUSTUS ab initio predictions without UTRs")
    ucsc_file = tmp_dir + "aug_ab_initio_ucsc.gtf"
    aug2ucsc_gtf(args.aug_ab_initio, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file,
                   "AUG_ab_initio_no_utr",
                   "AUGUSTUS ab initio predictions without UTRs",
                   this_color, visibility, args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "AUG_ab_initio_no_utr.html"
    write_trackDescription(html_file, "AUG_ab_initio_no_utr", "This gene prediction track with AUGUSTUS <i>ab initio</i> predictions (without UTRs) " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating AUGUSTUS with hints track '''

if args.aug_hints:
    print('Generating AUGUSTUS prediction with hints (no UTRs) track...')
    check_trackDb(trackDb_file, "AUGUSTUS predictions with hints without UTRs")
    ucsc_file = tmp_dir + "aug_hints_ucsc.gtf"
    aug2ucsc_gtf(args.aug_hints, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file,
                   "AUG_hints_no_utr",
                   "AUGUSTUS predictions with hints without UTRs",
                   this_color, visibility, args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "AUG_hints_no_utr.html"
    write_trackDescription(html_file, "AUG_hints_no_utr", "This gene prediction track with AUGUSTUS predictions with hints (without UTRs) " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating AUGUSTUS ab initio with UTRs track '''

if args.aug_ab_initio_utr:
    print('Generating AUGUSTUS ab initio prediction with UTRs track...')
    ucsc_file = tmp_dir + "aug_ab_initio_utr_ucsc.gtf"
    check_trackDb(trackDb_file, "AUGUSTUS ab initio predictions with UTRs")
    aug2ucsc_gtf(args.aug_ab_initio_utr, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file,
                   "AUG_ab_initio_utr",
                   "AUGUSTUS ab initio predictions with UTRs",
                   this_color, visibility, args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "AUG_ab_initio_utr.html"
    write_trackDescription(html_file, "AUG_ab_initio_utr", "This gene prediction track with AUGUSTUS <i>ab initio</i> predictions with UTRs " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating AUGUSTUS with hints and UTRs track '''

if args.aug_hints_utr:
    print('Generating AUGUSTUS prediction with hints and with UTRs track...')
    check_trackDb(trackDb_file, "AUGUSTUS predictions with hints and UTRs")
    ucsc_file = tmp_dir + "aug_hints_utr_ucsc.gtf"
    aug2ucsc_gtf(args.aug_hints_utr, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "AUG_hints_utr",
                   "AUGUSTUS predictions with hints and UTRs",
                   this_color, visibility, args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "AUG_hints_utr.html"
    write_trackDescription(html_file, "AUG_hints_utr", "This gene prediction track with AUGUSTUS predictions with hints and UTRs" +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating Gemoma filtered predictions track '''

if args.gemoma_filtered_predictions:
    print('Generating Gemoma filtered predictions track...')
    check_trackDb(trackDb_file, "Gemoma filtered predictions")
    aug_file = tmp_dir + "gemoma_aug.gtf"
    gemoma2aug_gtf(args.gemoma_filtered_predictions, aug_file)
    ucsc_file = tmp_dir + "gemoma_ucsc.gtf"
    aug2ucsc_gtf(aug_file, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "GeMoMa",
                   "Gemoma filtered predictions",
                   this_color, visibility, args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "GeMoMa.html"
    write_trackDescription(html_file, "GeMoMa", "This gene prediction track with Gemoma predictions " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating GlimmerHMM predictions track '''

if args.glimmer_gff:
    print('Generating GlimmerHMM predictions track...')
    check_trackDb(trackDb_file, "GlimmerHMM predictions")
    aug_file = tmp_dir + "glimmer_aug.gtf"
    glimmer2aug_gtf(args.glimmer_gff, aug_file)
    ucsc_file = tmp_dir + "glimmer_ucsc.gtf"
    aug2ucsc_gtf(aug_file, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "GlimmerHMM",
                   "GlimmerHMM predictions",
                   this_color, visibility, args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "GlimmerHMM.html"
    write_trackDescription(html_file, "GlimmerHMM", "This gene prediction track with GlimmerHMM predictions " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating SNAP predictions track (native SNAP, not MAKER-SNAP) '''

if args.snap_gff:
    print('Generating SNAP predictions track...')
    check_trackDb(trackDb_file, "SNAP predictions")
    aug_file = tmp_dir + "snap_aug.gtf"
    snap2aug_gtf(args.snap_gff, aug_file)
    ucsc_file = tmp_dir + "snap_ucsc.gtf"
    aug2ucsc_gtf(aug_file, ucsc_file)
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, "SNAP",
                   "SNAP predictions",
                   this_color, visibility, args.no_genePredToBigGenePred, False)
    html_file = hub_dir + "SNAP.html"
    write_trackDescription(html_file, "SNAP", "This gene prediction track with SNAP predictions " +
                           "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating general gene track with custom label '''

if args.gene_track:
    print('Generating gene track from file ' +
          args.gene_track[0] + ' with label ' + args.gene_track[1])
    check_trackDb(trackDb_file, args.gene_track[1])
    try:
        with open(args.gene_track[0], "r") as gene_file_handle:
            lineNo = 1
            for line in gene_file_handle:
                if(re.search(r'\s+AUGUSTUS\s+', line)):
                    ucsc_file = tmp_dir + "gene_track_ucsc.gtf"
                    if(os.path.isfile(ucsc_file)):
                        os.remove(ucsc_file)
                    aug2ucsc_gtf(args.gene_track[0], ucsc_file)
                    break
                if lineNo > 30:
                    ucsc_file = args.gene_track[0]
                    break
                lineNo = lineNo + 1
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              args.gene_track[0] + " for reading!")
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, args.gene_track[1],
                   args.gene_track[1], this_color, visibility, 
                   args.no_genePredToBigGenePred, True)
    html_file = hub_dir + args.gene_track[1] + ".html"
    write_trackDescription(html_file, args.gene_track[1], "This gene prediction track was added to the assembly hub as <i>general gene track with custom label</i>. " +
                           "It was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating traingenes track '''

if args.traingenes:
    print('Generating training gene track...')
    check_trackDb(trackDb_file, "Training genes")
    visibility_counter, visibility = set_visibility(visibility_counter)
    col_idx, this_color = set_color(col_idx, rgb_cols)
    make_gtf_track(trackDb_file, args.traingenes, ChromSizes_file, "Traingenes",
                   "Training genes", this_color, visibility,
                   args.no_genePredToBigGenePred, True)
    html_file = hub_dir + "Traingenes.html"
    write_trackDescription(html_file, args.traingenes, "This training gene track was added to the assembly hub by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)


''' Creating MAKER tracks '''

if args.maker_gff:
    print('Generating tracks from MAKER output file...')
    groups_txt_file = hub_dir + "groups.txt"
    extend_groups = True
    try:
        with open(groups_txt_file, "r") as groups_handle:
            for line in groups_handle:
                if re.search(r'Evidence generated by MAKER', line):
                    extend_groups = False
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " +
              groups_txt_file + " for reading!")
        quit(1)
    if extend_groups:
        try:
            with open(groups_txt_file, "a") as groups_handle:
                groups_handle.write(
                    "name makerEvidence\nlabel Evidence generated by MAKER\npriority 2\ndefaultIsClosed 0\n\n")
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  groups_txt_file + " for writing!")
            quit(1)
    else:
        print("Warning: MAKER evidence group configuration already exists in file " +
              groups_txt_file + "! You may not have intended to add further MAKER tracks " +
              "prior modifying existing tracks in trackDb.txt file...")
    maker_categ = {}
    try:
        with open(args.maker_gff, "r") as maker_handle:
            for line in maker_handle:
                line = line.strip()
                if re.search(r'^\S+\t\S+\t', line):
                    m_type = re.search(r'^\S+\t(\S+)\t', line).groups()
                    if (m_type[0] not in maker_categ) and (m_type[0] != '.'):
                        maker_categ[m_type[0]] = {}
                        maker_categ[m_type[0]]['lines'] = []
                        maker_categ[m_type[0]]['gene'] = False
                    if (m_type[0] != '.') and re.search(r'\tgene\t', line):
                        maker_categ[m_type[0]]['gene'] = True
                    if m_type[0] != '.':
                        maker_categ[m_type[0]]['lines'].append(line)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + args.maker_gff +
              " for reading!")
        quit(1)
    for feature in maker_categ:
        if maker_categ[feature]['gene'] is False:
            print("Creating BED MAKER track for feature " + feature + "...")
            this_maker_file = tmp_dir + feature + "_" + "maker.bed"
            maker_to_bed(maker_categ[feature]['lines'],
                         this_maker_file, ChromSizes_file)
            this_sorted_maker_file = this_maker_file + ".s"
            sort_bed3(this_maker_file, this_sorted_maker_file)
            this_maker_bb_file = hub_dir + feature + "_" + "maker.bb"
            bed2bigBed(12, this_sorted_maker_file,
                       ChromSizes_file, this_maker_bb_file, None)
            visibility_counter, visibility = set_visibility(visibility_counter)
            col_idx, this_color = set_color(col_idx, rgb_cols)
            info_to_trackDB(trackDb_file, feature + "_maker", "Evidence by MAKER from source " + feature, this_color,
                            "makerEvidence", 12, visibility, "bigBed")
            html_file = hub_dir + feature + "_maker" + ".html"
            write_trackDescription(html_file, feature + "_maker", "This gene prediction track was added to the assembly hub with the argument <i>--maker_gff MAKER_GFF</i>. " +
                                   "It was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                                   "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                                   args.email, this_color)
        else:
            print("Creating GTF MAKER track for feature " + feature + "...")
            this_maker_file = tmp_dir + feature + "_" + "maker.gtf"
            maker_to_gtf(maker_categ[feature]['lines'], this_maker_file)
            ucsc_file = this_maker_file + ".f"
            aug2ucsc_gtf(this_maker_file, ucsc_file)
            visibility_counter, visibility = set_visibility(visibility_counter)
            col_idx, this_color = set_color(col_idx, rgb_cols)
            make_gtf_track(trackDb_file, ucsc_file, ChromSizes_file, feature,
                           "MAKER gene predictions", this_color, visibility,
                           args.no_genePredToBigGenePred, True)
            html_file = hub_dir + "maker.html"
            write_trackDescription(html_file, "MAKER", "This MAKER gene track was added to the assembly hub by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                           "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                           args.email, this_color)
    # delete because it consumes a lot of RAM
    del maker_categ


''' Creating hints tracks '''

if args.hints:
    print('Separating hints into separate tracks...')
    hint_categ = {}
    try:
        with open(args.hints, "r") as hints_handle:
            for line in hints_handle:
                h_src, h_type = re.search(
                    r'^\S+\t(\S+)\t(\S+)\t', line).groups()
                if h_src not in hint_categ:
                    hint_categ[h_src] = {}
                    hint_categ[h_src][h_type] = []
                else:
                    if h_type not in hint_categ[h_src]:
                        hint_categ[h_src][h_type] = []
                hint_categ[h_src][h_type].append(line)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + args.hints +
              " for reading!")
        quit(1)
    # determine whether there have been hints tracks before, in this hub
    hint_file_no = 1
    if os.path.isfile(trackDb_file):
        try:
            with open(trackDb_file, "r") as trackDb_handle:
                for line in trackDb_handle:
                    if re.search(r'_hints_\d+$', line):
                        re_result = re.search(r'_hints_(\d+)$', line).groups()
                        hint_file_no = int(re_result[0]) + 1
        except IOError:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Failed to open file " +
                  trackDb_file + " for reading!")
    for h_src in hint_categ:
        for h_type in hint_categ[h_src]:
            this_hints_file = tmp_dir + h_src + "_" + \
                h_type + ".hints_" + str(hint_file_no)
            try:
                with open(this_hints_file, "w") as hints_handle:
                    hint_no = 0
                    for hint in hint_categ[h_src][h_type]:
                        hint_no = hint_no + 1
                        if (h_type == 'CDS') or (h_type == 'CDSpart') or (h_type == 'exon'):
                            # GTF like hints (all treated as exons)
                            if re.match(r'(\S+\t\S+\t)\S+(\t\S+\t\S+\t\S+\t\S+\t\S+\t)', hint):
                                first_part, second_part = re.match(
                                    r'(\S+\t\S+\t)\S+(\t\S+\t\S+\t\S+\t\S+\t\S+\t)', hint).groups()
                                if re.search(r'grp=', hint):
                                    seq_id, strand, grp = re.search(
                                        r'(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t\S+\t\S*grp=([^;]+)', hint).groups()
                                    tx_id = seq_id + "_" + grp + "_" + strand
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id" +
                                                       "\" " + tx_id + "\"; transcript_id \"" + tx_id + ".t1\"\n")
                                elif re.search(r'group=', hint):
                                    seq_id, strand, grp = re.search(
                                        r'(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t\S+\t\S*group=([^;]+)', hint).groups()
                                    tx_id = seq_id + "_" + grp + "_" + strand
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id" +
                                                       "\" " + grp + "\"; transcript_id \"" + grp + ".t1\"\n")
                                elif re.search(r'mult=', hint):
                                    mult = re.search(
                                        r'mult=(\d+)', hint).group(1)
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id \"hint_" + str(hint_no) +
                                                       "_mult_" + mult + "\"; transcript_id \"hint_" + str(hint_no) + "_mult_" + mult + ".t1\"\n")
                                else:
                                    hints_handle.write(first_part + "exon" + second_part + "gene_id \"hint_" + str(hint_no) +
                                                       "_mult_1\"; transcript_id \"hint_" + str(hint_no) + "_mult_1.t1\"\n")

                        else:
                            # segmental hints, i.e. introns
                            if re.match(
                                    r'(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t', hint):
                                seq, start, end, strand = re.match(
                                    r'^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t', hint).groups()
                                if re.search(r'mult=', hint):
                                    mult = re.search(
                                        r'mult=(\d+)', hint).group(1)
                                else:
                                    mult = 1
                                hints_handle.write(seq + "\t" + start + "\t" +
                                                   end + "\t" +
                                                   "mult_" +
                                                   str(mult) + "_" +
                                                   strand + "id_" + str(hint_no) + "\t1\t" + strand +
                                                   "\t" + start + "\t" + end +
                                                   "\t0\t1\t" +
                                                   str(int(end) - int(start)) +
                                                   "\t0\n")

            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed top open file " +
                      this_hints_file + " for writing!")
                quit(1)
            bb_file = hub_dir + h_src + "_" + h_type + \
                "_hints_" + str(hint_file_no) + ".bb"
            if (h_type == 'CDS') or (h_type == 'CDSpart') or (h_type == 'exon'):
                gp_file = this_hints_file + ".gp"
                bed_file = this_hints_file + ".bed"
                info_out_file = this_hints_file + ".infoOut.txt"
                gtf2bb(this_hints_file, gp_file, bed_file,
                       bb_file, info_out_file, ChromSizes_file, sort_tool, False) # check whether sometimes CDS hints have frame?
            else:
                this_sorted_hints_file = this_hints_file + ".sorted"
                sort_bed3(this_hints_file, this_sorted_hints_file)
                bed2bigBed(12, this_sorted_hints_file,
                           ChromSizes_file, bb_file, None)
            visibility_counter, visibility = set_visibility(visibility_counter)
            col_idx, this_color = set_color(col_idx, rgb_cols)
            info_to_trackDB(trackDb_file, h_src + "_" + h_type + "_hints_" +
                            str(hint_file_no), "Hints of type " + h_type +
                            " from source " + h_src, this_color,
                            "hints", 12, visibility, "bigBed")
            html_file = hub_dir + h_src + "_" + h_type + \
                "_hints_" + str(hint_file_no) + ".html"
            write_trackDescription(html_file, h_src + "_" + h_type + "_hints_" +
                                   str(hint_file_no), "This AUGUSTUS hints track was added to the assembly hub " +
                                   "was generated by <a href=\"https://github.com/Gaius-Augustus/MakeHub\">" +
                                   "make_hub.py</a> using <a href=\"http://hgdownload.soe.ucsc.edu/admin/exe\">UCSC tools</a>.</p>",
                                   args.email, this_color)


''' Delete temporary directory '''

if not args.no_tmp_rm:
    if args.verbosity > 0:
        print("Deleting temporary files...")
    shutil.rmtree(tmp_dir)

print('\nHub is ready, please copy to a web server, e.g.')
print('\"scp -r ' + args.outdir + "/" +
      args.short_label + " user@server:/target/location\"")
print('Feed the link to hub.txt to UCSC genome browser:')
print('\tMy Data -> Track Hubs -> MyHubs')
