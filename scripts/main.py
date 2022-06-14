# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import logging
import argparse

from region import GRegion
from ensembl_rest import getENSEMBLGeneLoc, getENSEMBLGFF, getENSEMBLFASTA, getGeneTranscriptLists, getMANETranscripts
from make_genbank import updateGFFRegions, renamingRecord, joinFeature, flattenFeatures, check_gff
from safeports import createSAFEPORTS
from recoding import recode
import os
import io
import os.path as path
import sys
from BCBio import GFF
from Bio import SeqIO

# Get the version:
version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

# A function to quit with an error:
def error(msg, exit_code=1):
    logging.error(msg)
    sys.exit(exit_code)

# A function to write data to a given file, with correct logging:
def writeDataFile(x, filename, datatype):
    try:
        with open(filename, 'w') as output_handle: output_handle.write(x)
    except Exception as err: error('failed to write {} data to {} ({})'.format(datatype, filename, err))


# Handle unix pipes properly:
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

# This is the main entrypoint to the safeports command script:
def main():
    # Define the CLI defaults:
    defaults = {'verbosity':'info', 'flank':0, 'output_dir':'./',\
                'sp_filter': ['ALL'], 'sp_length': 3,\
                'endprotect':0, 'transcriptpriority': 'max',
                'join':True,
                'primedesign': False}
    
    # Create the command line interface:
    parser = argparse.ArgumentParser(description='Annotate genome safe regions for gene editing')
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {0}'.format(version['__version__']))
    parser.add_argument('-v', '--verbose', dest='verbosity', default=defaults['verbosity'], choices=['error', 'warning', 'info', 'debug'], help='Set logging level (default {verbosity})'.format(**defaults))
    parser.add_argument('-f', '--gene-flank', dest='flank', metavar='n', default=defaults['flank'], type=int, help='Gene flanking region size (default {flank})'.format(**defaults))
    # The output options:
    parser.add_argument('-o', '--output-dir', dest='output_dir', metavar='dir', default=defaults['output_dir'], help='Output directory (default {output_dir})'.format(**defaults))
    parser.add_argument('-n', '--output-prefix', dest='output_prefix', metavar='prefix', help='Output file prefix (default uses input region)')
    # The input region:
    parser.add_argument(dest='region', metavar='region', help='Either gene name/symbol or region of interest')
    # The safeport options
    parser.add_argument('-p', '--safeport-filter', nargs='+', dest = 'sp_filter', default=defaults['sp_filter'], help='List of annotations that are NOT safeports (default {sp_filter})'.format(**defaults))
    parser.add_argument('-l', '--safeport-length', dest = 'sp_length', default=defaults['sp_length'], type = int, help='Minimum length for safeports (default {sp_length}nt)'.format(**defaults))
    # The recoding options
    parser.add_argument('-c','--coding-file', dest='codon_map', metavar='FILE', help = 'Optional codon mapping file')
    parser.add_argument('-e','--endprotect',default = defaults['endprotect'], type = int, help = 'Do you want to protect the end of each coding sequence by not modifing the last n nt?')
    parser.add_argument("-tp",'--transcriptpriority',type = str, default = defaults['transcriptpriority'], choices = ['max','soft','maneonly'], help = "Choose either 'max' (if none of the overlapping cds is a MANE transcript, the longest cds will be taken),'soft' (if none of the overlapping transcripts is a MANE, none will be recoded) or 'maneonly' (only MANE transcripts will be recoded).")
    # Optional
    parser.add_argument('-j', '--join', type = bool, default=defaults['join'], help = 'If True, a joined CompoundLocation annotations will be created for each transcript.')
    parser.add_argument('-pd', '--primedesign', type = bool, default=defaults['primedesign'], help = 'Outputs file input for primedesign tool if True.')
    # Parse the CLI:
    args = parser.parse_args()
    
    # Set up logging based on the verbosity level set by the command line arguments:
    logging.basicConfig(format='%(levelname)s: %(message)s', level=args.verbosity.upper())
    
    #Load the coding file:
    if args.codon_map is None:
        codonmap = []
        logging.warning('No codon map provided. {} will be used as recoding scheme!'.format(codonmap))
    else:
        codon_data = {}
        codonmap = []
        try:
            with open(args.codon_map, 'rt') as codon_handle:
                for line in codon_handle:
                    line = line.strip()
                    if len(line) == 0 or line.startswith('#'):
                        continue
                    line_data = line.split('\t')
                    if len(line_data) != 2:
                        logging.warning('invalid line in codon mapping file')
                        continue
                    if line_data[1] not in codon_data.keys():
                        codon_data[line_data[1]] = list()
                    codon_data[line_data[1]].append(line_data[0])
        except: error('failed to read codon mapping file {}'.format(args.codon_map))
        codonmap = []
        for k, v in codon_data.items():
            codonmap.append((v, k))
    
    logging.debug('\n ----------ADJUST REGION \n')
    # Check that the output directory exists, and make it if not:
    output_dir = os.path.expanduser(os.path.abspath(args.output_dir))
    if not os.path.exists(output_dir):
        try:
            logging.debug('creating output directory {}'.format(output_dir))
            os.makedirs(output_dir)
        except Exception as err: error('failed to create output directory {}'.format(output_dir))
    if not os.path.isdir(output_dir): 
        error('specified output {} is not a directory'.format(output_dir))        
    
    # If we've been given a gene ID, convert it to a region:
    try: 
        query_region = GRegion.fromString(args.region)
    except:
        logging.debug('mapping ENSEMBL identifier {} to region'.format(args.region))
        try:
            query_region = getENSEMBLGeneLoc(args.region)
            logging.debug('ENSEMBL identifier {} maps to {}'.format(args.region, query_region))
        except Exception as err: error('failed to map identifier ({})'.format(err))
    
    # Apply region flanking:
    if args.flank > 0:
        logging.debug('expanding query region by {} nt'.format(args.flank))
        query_region = query_region.flank(args.flank)
    
    # Log the query region:
    logging.info('using query region {}'.format(query_region))

    logging.debug('\n ----------BUILD SETTINGS \n')   
    # Build the output file prefix:
    if args.output_prefix is None: 
        args.output_prefix = 'chr{}_{}_{}'.format(query_region.chromosome, query_region.start, query_region.end)
    #args.output_prefix = path.join(output_dir, args.output_prefix) # commented out since I want to use the output_prefix without path at one point
    logging.debug('output file prefix is {}'.format(path.join(output_dir, args.output_prefix)))

    # Build safeport filter settings
    if args.sp_filter == ['ALL']:
        args.sp_filter = ['exon','CDS','CTCF_binding_site','promoter','sequence_variant','repeat_region','TF_binding_site']
    elif args.sp_filter == ['noSNP']:
        args.sp_filter = ['exon','CDS','CTCF_binding_site','promoter','repeat_region','TF_binding_site']
    ### else: the user input is taken. How to make sure the user is not typing **wrong** things?
    logging.info('The following list of annotations indicates where a safeport should NOT be created: {} \n'.format(args.sp_filter))

    logging.debug('\n ----------IMPORT DATA \n')
    
    # Get the ENSEMBL GFF data and save it to file:
    try: 
        gff_data = getENSEMBLGFF(query_region)
        logging.debug('  GFF3 data was downloaded.')
    except Exception as err: 
        error('failed to pull GFF data from ENSEMBL ({})'.format(err))
    writeDataFile(gff_data, '{}.gff'.format(path.join(output_dir, args.output_prefix)), 'GFF')
    
    # Get the ENSEMBL FASTA data and save it to file:
    try: 
        fasta_data = getENSEMBLFASTA(query_region)
        logging.debug('  FASTA data was downloaded.')
    except Exception as err: 
        error('failed to pull FASTA data from ENSEMBL ({})'.format(err))
    writeDataFile(fasta_data, '{}.fasta'.format(path.join(output_dir, args.output_prefix)), 'FASTA')
    
    # Get list of genes and transcripts in defined region
    try:
        genelist, transcriptlist = getGeneTranscriptLists(GFF.parse('{}.gff'.format(path.join(output_dir, args.output_prefix))))
        logging.debug('  Following genes were identified: {}'.format(genelist))
        logging.debug('  Following transcripts were identified: {}'.format(transcriptlist))
    except:
        genelist = []
        transcriptlist = []
        logging.warning('Genes and transcripts could not be extracted. This will impact the creation of joined features and recoding.')
    
    # Get the MANE list
    try:
        manelist = getMANETranscripts(genelist)
    except:
        manelist = []
        logging.warning('Downloading the list of MANE transcripts from Tark failed. The recoding will be done without taking into account MANE transcripts.')
    set_mane = [args.transcriptpriority,manelist]
    logging.info('For recoding, the overlapping transcript is chosen based on transcriptpriority {}. \n The identified MANE transcripts include {}'.format(args.transcriptpriority, manelist))
    
    # Write out the updated GFF file:
    try:
        with open('{}-fixed.gff'.format(path.join(output_dir, args.output_prefix)), 'wt') as gff_out_handle:
            with io.StringIO(gff_data) as gff_in_handle:
                GFF.write(updateGFFRegions(GFF.parse(gff_in_handle), query_region), gff_out_handle)
    except Exception as err: 
        error('failed to updated GFF coordinates ({})'.format(err))

    logging.debug('\n ----------START PARSING \n')
    # Create parser for combined GFF3 and FASTA
    try:
        fasta_input = SeqIO.to_dict(SeqIO.parse('{}.fasta'.format(path.join(output_dir, args.output_prefix)), "fasta"))
        data_iter = GFF.parse('{}-fixed.gff'.format(path.join(output_dir, args.output_prefix)),fasta_input)
    except Exception as err:
        error('Failed to combine GFF3 and FASTA data.')

    # Flatten features and rename record
    try:
        # data_handle = flattenFeatures(data_handle)
        data_handle = check_gff(data_iter) # this is replacing flattenFeatures for now, since flattenFeatures doesn't work
        data_handle = renamingRecord(data_handle,query_region,args.output_prefix) # this works
    except:
        logging.warning('Renaming record failed.')
    
    # Create safeports
    try:
        data_handle = createSAFEPORTS(data_handle,len(query_region),args.sp_filter,args.sp_length)
    except:
        logging.warning('Creating safeports failed.')

    # Recode
    if len(codonmap) > 0:
        try:
            logging.info('Recoding map used: {}'.format(codonmap))
            data_handle = recode(data_handle, codonmap, args.endprotect, set_mane, output_dir, args.output_prefix, args.primedesign)
        except Exception as err:
            error('Recoding failed.')
    else:
        logging.info('Recoding is skipped due to empty codonmap.')

    # Create join-annotations for better visualization of the data
    if args.join == True:
        try:
            data_handle = joinFeature(data_handle, transcriptlist)
        except:
            logging.warning('Joined feature annotations were not created because of an error. GenBank file has been created anyway.')
 
    # Write GenBank file
    try:
        SeqIO.write(data_handle,'{}.gb'.format(path.join(output_dir, args.output_prefix)),'genbank')
        logging.info('Data has been written into a GenBank file.')
    except Exception as err:
        error('Writing the GenBank file failed.')

if __name__ == '__main__':
    main()
