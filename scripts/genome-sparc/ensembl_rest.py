#from safeports.region import GRegion
from region import GRegion
import logging
import urllib.parse
import requests
import numpy as np

def getENSEMBLGeneLoc(gene, ensembl_netloc='rest.ensembl.org'):
    # Build the ENSEMBL REST request URI:
    if gene.startswith('ENS'): 
        request_path = 'lookup/id/{}'.format(gene)
    else: 
        request_path = 'lookup/symbol/homo_sapiens/{}'.format(gene)
    request_uri = urllib.parse.urlunsplit(('https', ensembl_netloc, request_path, '', ''))
    # Query ENSEMBL:
    res = requests.get(request_uri, headers={"Content-Type" : "application/json"})
    if not res.ok: 
        res.raise_for_status()
    # Capture the output:
    res_json = res.json()
    return GRegion(res_json['seq_region_name'], int(res_json['start']), int(res_json['end']))

def getENSEMBLGFF(region, ensembl_netloc='rest.ensembl.org'):
    """ Import GFF3 for region from Ensemble. """
    # Build the ENSEMBL REST request URI:
    request_path = 'overlap/region/human/{}'.format(region)
    request_uri = urllib.parse.urlunsplit(('https', ensembl_netloc, request_path, 'feature=regulatory;feature=cds;feature=exon;feature=variation;feature=repeat;variant_set=1kg_3_com', ''))
    # Query ENSEMBL:
    res = requests.get(request_uri, headers={"Content-Type" : "text/x-gff3"})
    if not res.ok: 
        res.raise_for_status()
    return res.text

def getENSEMBLFASTA(region, ensembl_netloc='rest.ensembl.org'):
    """ Import FASTA for region from Ensemble. """
    # Build the ENSEMBL REST request URI:
    request_path = 'sequence/region/human/{}:1'.format(region)
    request_uri = urllib.parse.urlunsplit(('https', ensembl_netloc, request_path, '', ''))
    # Query ENSEMBL:
    res = requests.get(request_uri, headers={"Content-Type" : "text/x-fasta"})
    if not res.ok: 
        res.raise_for_status()
    # Return with an updated FASTA header:
    return '>{}\n{}'.format(region.chromosome, res.text.partition('\n')[2])

def getMANETranscripts(genelist, tark_netloc='dev-tark.ensembl.org'):
    """ Import list of MANE transcripts from Tark. """
    tarkdict_list = []
    mane_list = []
    # Query Tark for every gene in the genelist
    for x in genelist:
        # request_path = '/api/transcript/search/?identifier_field={}'.format(x)
        # request_uri = urllib.parse.urlunsplit(('http', tark_netloc, request_path,'expand=%transcript_release_set%2Cgenes%'))
        ext = '/api/transcript/search/?identifier_field={}&expand=%transcript_release_set%2Cgenes%'.format(x)
        request_uri = 'http://' + tark_netloc + ext
        logging.debug(request_uri)
        # Query Tark
        res = requests.get(request_uri)
        if not res.ok:
            res.raise_for_status()

        tarkdict_list.append(res.json())
        logging.debug('  Appended to Tark-Dict list: {}'.format(res.json()))

    # Iterate over transcripts and get MANE transcript ID
    for entrydict_list in tarkdict_list:
        for transcript in entrydict_list:
            if 'mane_transcript_type' in transcript.keys() and transcript['stable_id'].startswith('EN') and transcript['assembly'] == 'GRCh38':
                    mane_list.append(transcript['stable_id'])
    return mane_list

def getGeneTranscriptLists(gff_handle):
    """ Extracts a list of all genes and transcripts in GFF file. """
    gene_list = []
    transcript_list = []
    for record in gff_handle:
        # (Recursively) iterate over features and extract genes and transcripts
        for feature in record.features:
            # Genes are in first level of annotations and don't have subfeatures
            if feature.type == 'gene':
                gene_list.append(feature.qualifiers['gene_id'][0])
            else:
                gene_list = getCDS(feature, transcript_list)
    return np.unique(gene_list), np.unique(transcript_list)

            
def getCDS(feature, transcript_list, log_depth = 0):
    if feature.type == 'CDS':
        tID = feature.qualifiers['Parent'][0].split(':')[1]
        transcript_list.append(tID)
    for subfeature in feature.sub_features:
        getCDS(subfeature, transcript_list, log_depth = log_depth + 1)
    return transcript_list

