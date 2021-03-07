import logging
from Bio import SeqFeature

def createSAFEPORTS(gff_iter,g_len,safeport_filter,min_len):
    """ Create safeport annotations based on existing annotations in the SeqRecord. """
    for rec in gff_iter:
        # Identify where a safeport should be created
        saferegions = getSAFERegions(rec,g_len,safeport_filter,min_len)
        logging.info('{} new safeport locations identified. \n'.format(len(saferegions)))
        
        # create annotation for safeports
        for x in saferegions:
            new_feature = createFEATUREannot(x,"safeport", s = 0)
            rec.features.append(new_feature)
            logging.debug('new feature created: {}'.format(new_feature))
        yield rec

def getSAFERegions(seq_record,g_len,filter_criteria,min_len):
    """ Identify regions that have no relevant annotation."""
    try: danger_regions = getDANGERRegions(seq_record,g_len,filter_criteria)
    except:
        logging.info('No relevant annotations detected - whole region is safe')
        return([[0,g_len]])
    try: safe_regions = filterRegions(invertRegions(danger_regions,g_len), min_len)
    except: logging.error('Getting saferegions failed.')
    return safe_regions

def getDANGERRegions(seq_record,g_len,filter_criteria):
    """ Identify regions that should be avoided for safeports."""
    danger_regions = []
    logging.debug('Filter criteria: {}'.format(filter_criteria))
    for i in range(len(seq_record.features)):
        # determine which annotations should be taken into account
        if seq_record.features[i].type in filter_criteria:
            #logging.debug('annotation: type {}'.format(seq_record.features[i].type))
            loc_list = [int(seq_record.features[i].location.start),int(seq_record.features[i].location.end)]
            danger_regions.append(loc_list)
    if len(danger_regions) == 0:
        raise Exception('No relevant annotations for danger regions detected.')
        # return([[0,g_len]]) # whole region is returned a safe region

    logging.info('{} annotations taken into account for creating safeports'.format(len(danger_regions)))
    danger_regions = mergeRegions(danger_regions)
    logging.debug('Excluded danger regions after merging: {}'.format(danger_regions))
    return(danger_regions)

def mergeRegions(regionlist):
    regionlist.sort(key=lambda interval: interval[0])
    merged = [regionlist[0]]
    if len(regionlist)>1:
        for curr in regionlist:
            previous = merged[-1]
            if curr[0] <= previous[1]: # if current start position is smaller or equal to the previous end position, change previous end position to new max.
                previous[1] = max(curr[1], previous[1])
            else:
                merged.append(curr) # if start position is higher, add a new danger region 
    return merged

def filterRegions(regionlist, minlength):
    filtered = []
    for x in regionlist:
        if x[1] - x[0] >= minlength:
            filtered.append(x)
    return(filtered)

def createFEATUREannot(loc_range, featuretype,s):
    """ Creates a new SeqFeature with ExactPositions based on range."""
    location = SeqFeature.FeatureLocation(SeqFeature.ExactPosition(loc_range[0]),SeqFeature.ExactPosition(loc_range[1]))
    new_feature = SeqFeature.SeqFeature(location,type=featuretype,strand=s)
    return(new_feature)

def invertRegions(regionlist,end, start = 0):
    inverted = []
    for i in range(len(regionlist)):
        if i == 0 and regionlist[i][0] > start:
            # if the first few bases are free of an annotation, we want to define a safe region before the first annotation
            inverted.append([start,regionlist[i][0]]) # before annotation starts
            if len(regionlist) == 1: # if this is the only danger region, the end position is at the end of the region of interest
                inverted.append([regionlist[i][1],end]) # add region after until the end
            else: 
                inverted.append([regionlist[i][1],regionlist[i+1][0]]) # add region after until next annotation
        elif i == len(regionlist)-1:
            # last danger region
            if regionlist[i][1] >= end: # region of interests already ends, so no more safeport
                continue 
            else:
                inverted.append([regionlist[i][1],end])
        else:
            inverted.append([regionlist[i][1],regionlist[i+1][0]])   
    return inverted