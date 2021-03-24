from BCBio import GFF
from Bio import Seq, SeqFeature, SeqRecord
import logging
import time

# A function to clip a feature to a specified region, and returns a BioPython FeatureLocation.
# If required, it translates the location first
def translateFeatureLocation(location, region, translation=0):
    location2 = location + translation + 1
    if location2.end < 0:
        logging.debug('Error-prone feature detected: {}'.format(location2))
        return SeqFeature.FeatureLocation(
            start = 0,
            end = 0,
            strand = 0
        )
    else:
        return SeqFeature.FeatureLocation(
            start = max(0, location2.start),
            end = min(location2.end, region.end),
            strand = location2.strand
        )

# A function to update a feature's location:
def updateFeature(feature, region, translation=0, log_depth=0):
    location = feature.location
    feature.location = translateFeatureLocation(location, region, translation=translation)
    logging.debug('{}updated {}feature {} {} -> {}'.format(' ' * (2 * log_depth), 'sub-' * log_depth,feature.type, location, feature.location))
    for subfeature in feature.sub_features:
        updateFeature(subfeature, region, translation, log_depth = log_depth + 1)

def updateGFFRegions(gff_iterator,region):    
    # Iterate through all positions in the GFF datastream:
    for record in gff_iterator:
        for feature in record.features:
            updateFeature(feature, region, translation=-region.start)
        # Update the record sequence:
        record.seq = record.seq[region.start:region.end]
        # Return the updated record:
        yield record

def addFeature(feature, feature_list):
    feature.qualifiers['note'] = feature.id # changes name for visualization
    if feature.type == 'remark':
        logging.debug('Following remark has been removed: \n {}'.format(feature))
    else: 
        feature_list.append(feature, feature_list)
        for subfeature in feature.sub_features:
            addFeature(subfeature, feature_list)
    return feature_list

def flattenFeatures(gff_iterator):
    """ Brings all the subfeatures and subsubfeatures to the same level with the features for GenBank format"""
    for record in gff_iterator:
        record.annotations["molecule_type"] = "DNA"
        feature_list = []
        for feature in record.features:
            feature_list = addFeature(feature, feature_list)
        record.features = feature_list
        logging.info('{} features have been flattened.'.format(len(feature_list)))
        yield record


def joinFeature(gff_iterator, transcriptlist):
    """Creates joint feature annotations for CDS """
    if len(transcriptlist) == 0:
        return ValueError('Transcriptlist is empty.')

    for rec in gff_iterator:
        transcripts = {k:[] for k in transcriptlist}
        # Get transcript information and locations
        for i in range(len(rec.features)):
            if rec.features[i].type == 'CDS':
                loc = SeqFeature.FeatureLocation(rec.features[i].location.start,rec.features[i].location.end,rec.features[i].strand)
                tID = rec.features[i].qualifiers['Parent'][0].split(":")[1]
                if len(transcripts[tID]) == 0:
                    tQual = {'ID': rec.features[i].qualifiers['ID'],\
                            'assembly_name': rec.features[i].qualifiers['assembly_name'],\
                            'protein_id': rec.features[i].qualifiers['protein_id'],\
                            'source': rec.features[i].qualifiers['source']}
                    transcripts[tID].extend([tQual,loc])
                else:
                    transcripts[tID].append(loc)
        # Create annotation for each transcript
        for tKey,tVal in transcripts.items():
            new_feature = SeqFeature.SeqFeature(SeqFeature.CompoundLocation(tVal[1:]),type = 'joint CDS',id = tKey, qualifiers = tVal[0])
            rec.features.append(new_feature)
            logging.debug('Joined feature has been created for {}.'.format(tKey))
        yield rec

#### old functions since some of the above don't work


def check_gff(gff_iterator):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    Modified from https://github.com/chapmanb/bcbb since Bio.Aplhabet was removed from BioPython.
    """
    # should be included partially into flattenFeatures (the isinstance is not necessary)
    for rec in gff_iterator:
        rec.annotations["molecule_type"] = "DNA"
        if isinstance(rec.seq, Seq.UnknownSeq):
            logging.warning("Warning: FASTA sequence not found for {} in GFF file".format(rec.id))
        yield flatten_features(rec)

def flatten_features(rec):
    # should be replaced by flattenFeatures
    # Modified from https://github.com/chapmanb/bcbb
    """Make sub_features in an input rec flat for output.

    GenBank does not handle nested features, so we want to make
    everything top level.

    Modified by JW Nov 2020 to include naming of features.
    Modified by JW March 2021 to remove features that have [0:0] location.
    """
    out = []
    for f in rec.features:
        #logging.info('current feature: {}'.format(f))
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                if curf.type == 'remark':
                    logging.debug('Following remark annotation removed: \n {}'.format(curf))
                    continue
                
                elif curf.location.start == 0 and curf.location.end == 0:
                    logging.debug('Location was previously of of range and is now removed: {}'.format(curf))
                
                else:
                    out.append(curf)

                curf.qualifiers['note'] = curf.id

                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
      
            cur = nextf
  
    rec.features = out
    logging.info('{} features have been flattened.'.format(len(out)))
    return rec

def renamingRecord(gff_iterator,region,prefix):
    for record in gff_iterator:
        record.description = 'Recoded Chr{} from {} to {}'.format(region.chromosome, region.start, region.end)
        record.name = '{}'.format(prefix)
        record.annotations['taxonomy'] = ["Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa","Bilateria", "Deuterostomia","Chordata","Craniata","Vertebrata","Gnathostomata","Teleostomi","Euteleostomi","Sarcopterygii","Dipnotetrapodomorpha","Tetrapoda","Amniota","Mammalia","Theria","Eutheria","Boreoeutheria","Euarchontoglires","Primates","Haplorrhini","Simiiformes","Catarrhini","Hominoidea","Hominidae"]
        record.annotations['date'] = time.strftime("%d-%b-%Y").upper()
        record.annotations['organism'] =  "Homo sapiens"
        record.id = 'chromosome:GRCh38:{}:{}:{}'.format(region.chromosome, region.start, region.end)
        logging.debug('Renaming of record to {} (name) and {} (id)'.format(record.name,record.id))
        yield record
