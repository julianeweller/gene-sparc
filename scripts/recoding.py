import re
import logging
from Bio import SeqFeature
from operator import itemgetter
from library_design import export_primedesign_format

def recode(gff_iterator, codonmap, endprotect, set_mane, outdir, outprefix, primedesign):
    """ Replaces the codons based on codonmap in the Sequence Record, 
    but protects nucleotides at the beginning or end of the coding sequence
    and handles overlaps based on MANE transcript list and rules definet by set_mane.
    """
    complement_dict = {'A':'T','T':'A','C':'G','G':'C'}
    
    for record in gff_iterator:
        
        cds_ranges = getAnnotRange(record,target = 'CDS') # takes into account transcript ID, returns [[start,end,phase,transcript_id, strand]]
        
        # Handle overlaps --> select the cds that should be taken into account. Alternativly, in the long term, a duplication of the overlapping sequence could be introduced here.
        # Use set_mane to only choose the well supported transcripts / default transcript per gene.
        selected_cds_ranges, conflict_cds_ranges = selectCDS(cds_ranges,set_mane)
        logging.debug('Selected cds ranges: {}'.format(selected_cds_ranges))

        # Create annotations for conflict cds ranges --> user should look at them manually
        logging.info('Conflict cds ranges: {}'.format(conflict_cds_ranges))
        for x in conflict_cds_ranges:
            create_feature_annot([x[0],x[1]], "recoding-conflict",x[4])
            logging.debug('  Recoding conflict annotation created for {}'.format(x))
        
        # Selected CDS: Create a dictionary for every transcript (key) and the corresponding ranges that have been selected (value: list of ranges (list) with 0 - start, 1 - end, 2 - phase, 3 - transcript id, 4 - strand)
        cds_ranges_for, cds_ranges_rev = createRangeDict(selected_cds_ranges) 
        logging.debug('Ranges of relevance: \n  Forward: {} \n  Reverse: {}'.format(cds_ranges_for, cds_ranges_rev))
        
        # All CDS: Create a dictionary for every transcript (key) and the corresponding ranges (value: list of ranges (list) with 0 - start, 1 - end, 2 - phase, 3 - transcript id, 4 - strand)
        cds_ranges_for_all, cds_ranges_rev_all = createRangeDict(cds_ranges) # we need all ranges to take into account splitted ORFs
        logging.debug('all ranges: \n  Forward: {} \n  Reverse: {}'.format(cds_ranges_for_all,cds_ranges_rev_all))
        
        # Create a dictionary with positions where the codon should be inserted, key: codon that should be insterted, value: list of positions where it should be inserted  
        replace_map_for = getReplaceMap(record,'forward',cds_ranges_for,cds_ranges_for_all,codonmap,endprotect)
        replace_map_rev = getReplaceMap(record,'reverse',cds_ranges_rev,cds_ranges_rev_all,codonmap,endprotect)
        logging.debug('Following codons will we inserted: \n  Forward: {} \n  Reverse: {}'.format(replace_map_for, replace_map_rev))


        # output for primeDesign software if wanted
        if primedesign == True:
            export_primedesign_format(record, replace_map_for,replace_map_rev, outdir, outprefix)

        # replace codons
        mutable_seq = record.seq.tomutable() # makes sequence mutable
        for x in replace_map_for.keys():
            codon = list(x)
            for pos in replace_map_for[x]:
                # create a new annotation for each recoded nucleotide with original nucleotide
                for k in pos:
                    new_feature = create_feature_annot([k,k+1],'replaced-'+mutable_seq[k], 1)
                    record.features.append(new_feature)
                # Recode    
                mutable_seq = replaceCODON(pos,codon,mutable_seq)
        for x in replace_map_rev.keys():
            codon = [complement_dict[y] for y in list(x)]
            for pos in replace_map_rev[x]:
                # create a new annotation for each recoded cnucleotide with original nucleotide
                for k in pos:
                    new_feature = create_feature_annot([k,k+1],'replaced-'+complement_dict[mutable_seq[k]], 0)
                    record.features.append(new_feature)
                mutable_seq = replaceCODON(pos,codon,mutable_seq)
 
        # finish editing
        record.seq = mutable_seq.toseq()

        yield record


def replaceCODON(positionlist, codon, sequence):
    """ Replacing codon in sequence at position """
    sequence[positionlist[0]] = codon[0]
    sequence[positionlist[1]] = codon[1]
    sequence[positionlist[2]] = codon[2]
    return(sequence)

def getReplaceMap(record, strand, rangedict, all_ranges_dict, codonmap, endprotect):
    """ Find the positions for the codons of interest in the + strand
    """
    replace_map = {}
    
    for tID, regionlist in rangedict.items():
        i = 0
        for region in regionlist:
            # determine shift based on region information in gff file
            shift = region[2]
            # check if ORF is complete or if part of next CDS has to be included
            curr_remainder = (region[1] - region[0] - shift) % 3
            # Create sub-region with full codons (include/exclude partial codons that are splitted between CDSs)
            sub_reg = getSubRegion(record, region, tID, regionlist, all_ranges_dict,curr_remainder, shift, strand)
                
            # find the positions of the nucleotides that should be replaced for each codon
            for x in codonmap:
                # Find occurences in substring
                if strand == 'forward':
                    occ_list = find_occurence(x[0],str(sub_reg.seq),endprotect) # [[pos of 1st codon, 1 for being last codon / 0 for every codon within CDS]]
                    # add start of region to have the real position
                    occ_list = [[x+region[0]+shift,n] for x,n in occ_list]
                elif strand == 'reverse':
                    occ_list = find_occurence(x[0],str(sub_reg.seq.reverse_complement()),endprotect) # [[pos of 1st codon, 1 for being last codon / 0 for every codon within CDS]]
                    # get the real position
                    occ_list = [[region[1]-shift-x-1,n] for x,n in occ_list]
                else:
                    logging.error('Strand direction not clearly defined. strand = {}'.format(strand))
                # get list of nucleotides that should be replaced
                nt_list = getNTReplace(occ_list,regionlist,curr_remainder,strand,i)
                # add the nucleotide list to the dictionary
                replace_map[x[1]] = replace_map.get(x[1], []) + nt_list
            
            # update i
            i += 1

    return(replace_map)

def getSubRegion(record, region, tID, regionlist, all_ranges_dict,curr_remainder, shift, strand):
    # Find target codon nucleotides in target sequence
    if curr_remainder != 0:
        # Find current region in all region list
        m = all_ranges_dict[tID].index(region) # position of this region in all ranges for this transcript
        # Find next region in all region list
        try:
            next_reg = all_ranges_dict[tID][m+1]
        except:
            logging.warning('next region not found for current region {}. It seems like the transcripts stops without finishing the reading frame. The last {}nt were not considered for recoding.'.format(region, (3-curr_remainder)))
        # If next region is in regionlist of interest, all of them should be taken into account
        if next_reg in regionlist:
            if strand == 'forward':
                next_sub_reg = record[next_reg[0] : next_reg[0] + (3-curr_remainder)] # get the first partial codon from next range
                curr_sub_reg = record[region[0]+shift : region[1]]
            elif strand == 'reverse':
                next_sub_reg = record[next_reg[1]-(3-curr_remainder) : next_reg[1]] # get the last partial codon from next range
                curr_sub_reg = record[region[0]:region[1]-shift] # region stops earlier since strand is read reversly
            else:
                logging.error('Strand direction is not clearly indicated')
            # Get target subset of rec but without remainder from previous region 
            sub_reg = next_sub_reg + curr_sub_reg # join to create a subregion with complete codons 

        # else, the splitted codon won't be taken into account
        else:
            if strand == 'forward':
                sub_reg = record[region[0]+shift : region[1] - curr_remainder] # region starts later, since its for stand and end earlier due to not taking into account split codon
            elif strand == 'reverse':
                sub_reg = record[region[0]+curr_remainder:region[1]-shift] # region stops earlier since strand is read reversly, and starts later due to reverse strand and not taking into account split codon
            else:     
                logging.error('Strand direction is not clearly indicated')
    else:
        # Get target subset of rec but without remainder from previous region
        if strand == 'forward':
            sub_reg = record[region[0]+shift : region[1]]
        elif strand == 'reverse':
            sub_reg = record[region[0] : region[1]-shift]
        else:
            logging.error('Strand direction is not clearly indicated')
    
    return sub_reg

def getNTReplace(occ_list,regionlist,curr_remainder,strand, i):
    nt_list = []
    for occ in occ_list:
        if occ[1] == 1:
            if curr_remainder == 2:
                if strand == 'forward':
                    second_nt = occ[0]+1 # still part of current region
                    third_nt = regionlist[i+1][0]+1 # part of next region
                elif strand == 'reverse':
                    second_nt = occ[0]-1 # still part of current region
                    third_nt = regionlist[i+1][1]-1 # part of next region starting from the end
                else:
                    logging.error('Strand direction is not clearly indicated')
                logging.debug('      last codon positions: {}, remainder 2'.format([occ[0],second_nt,third_nt]))
            elif curr_remainder == 1: #current remainder should be one
                if strand == 'forward':
                    second_nt = regionlist[i+1][0]+1 # part of next region
                    third_nt = regionlist[i+1][0]+2 # also part of next region
                elif strand == 'reverse':
                    second_nt = regionlist[i+1][1]-1 # part of next region starting from the end
                    third_nt = regionlist[i+1][1]-2 # part of next region starting from the end
                logging.debug('      last codon positions: {}, remainder 1'.format([occ[0],second_nt,third_nt]))
            elif curr_remainder == 0:
                if strand == 'forward':
                    second_nt = occ[0]+1 # still part of current region
                    third_nt = occ[0]+2 # part of next region
                elif strand == 'reverse':
                    second_nt = occ[0]-1
                    third_nt = occ[0]-2
                logging.debug('      last codon positions: {}, remainder 0'.format([occ[0],second_nt,third_nt]))
            else:
                logging.warning("Problem: remainder should be 1 or 2, but is {}".format(curr_remainder))
                logging.warning('The codon at {} should be replaced, but failed.'.format(occ))
                continue
        else: 
            if strand == 'forward':
                second_nt = occ[0]+1 # still part of current region
                third_nt = occ[0]+2 # part of next region
            elif strand == 'reverse':
                second_nt = occ[0]-1
                third_nt = occ[0]-2
            else:
                logging.error('Strand direction not clearly indicated.')
        nt_list.append([occ[0],second_nt,third_nt])
    return(nt_list)

def selectCDS(rangelist, set_mane):
    # sort the list based on start
    rangelist.sort(key=lambda interval: interval[0])
    # Get settings for transcript priority and list of MANE transcripts
    transcriptpriority = set_mane[0]
    mane_list = set_mane[1]
    # Store selected ranges in list
    keep_cds = []
    conflict_cds = []
    j = 0 # to skip items that have already been taken into account
    for i in range(len(rangelist)):
        # skip the cds that have already been taken into account
        if i < j and i != 0: continue
        # for last cds, no comparison with next cds can be made. If it's not skipped, it needs to be recoded
        elif i == len(rangelist)-1:
            keep_cds.append(rangelist[i])
            continue
        # if we are only interested in MANE transcripts for recoding: loop to the next range if it's not a MANE transcript
        elif transcriptpriority == 'maneonly' and rangelist[i][3] not in mane_list: continue
        
        # check if current start is before the end of the last appended cds --> overlap with added
        elif len(keep_cds) != 0 and rangelist[i][0] < keep_cds[-1][1]:
            # check the orf: end of last orf + shift - end of current of + shift needs to be dividable by 3
            subremainder = (keep_cds[-1][1] + keep_cds[-1][2] - rangelist[i][0] + rangelist[i][2]) % 3
            # if the end is also before the end of the last appended cds and they are in the reading frame
            if rangelist[i][1] <= keep_cds[-1][1] and subremainder == 0: continue        
            # if soft, last one has to be mane. So if current is mane add to conflict
            elif transcriptpriority == 'soft':
                if rangelist[i][3] in mane_list: conflict_cds.append(rangelist[i]) # both are MANE, so add to conflict
                else: continue # it overlaps with a MANE, so it's not recoded
            # if maximize, last one doesn't have to be MANE
            elif transcriptpriority == 'max':
                # if both are MANE add the current one to conflict
                if rangelist[i][3] and keep_cds[-1][3] in mane_list: conflict_cds.append(rangelist[i])
                # if only last one is MANE, don't recode current one
                elif keep_cds[-1][3] in mane_list: continue
                # if only current one is MANE, don't recode last one
                elif rangelist[i][3] in mane_list: keep_cds[-1] = rangelist[i]
                else: 
                    conflict_cds.append(rangelist[i])
                    logging.warning("It looks like you have a special case that has not been considered. Current cds is added to conflict")
            else:
                conflict_cds.append(rangelist[i])
                logging.warning("It looks like you have a special case that has not been considered. The current cds will be added to conflict")

        # check if end of current cds is higher than the start of the next one --> overlap with next
        elif rangelist[i][1] > rangelist[i+1][0]:
            # check if more transcripts start at the same position
            n = i + 1
            subrangelist = [rangelist[i]]
            while n <= len(rangelist)-1 and rangelist[n][0] == rangelist[i+1][0]:
                subrangelist.append(rangelist[n])
                n += 1
            j = n    
            # check which transcript is a MANE transcript
            submanelist = []

            for k in range(len(subrangelist)):
                if subrangelist[k][3] in mane_list:
                    submanelist.append(k) # k is the position of transcript in the subrangelist
            # if there is only one mane transcript, everything is ok and this will be used for recoding
            if len(submanelist) == 1:
                keep_cds.append(subrangelist[k])
            # if no transcript is the mane transcript, check if they are in the same reading frame
            elif len(submanelist) == 0:
                if transcriptpriority == 'soft':
                    # in the soft regime, in case of an overlap of two non-MANE transcripts, nothing is recoded
                    continue
                elif transcriptpriority == 'max':
                    # decide which one to recode
                    orf_shiftlist = [item[2] for item in subrangelist]
                    # if they have the same orf shift
                    if (len(set(orf_shiftlist))==1):
                        # if they have also the same end:
                        if len(set([item[1] for item in subrangelist])):
                            keep_cds.append(subrangelist[0])
                        # they have different end positions --> take the longer one
                        else:
                            cds_id, cds_maxvalue = max(enumerate(map(itemgetter(1), subrangelist)),key=itemgetter(1))
                            keep_cds.append(subrangelist[cds_id])
                    # if they don't have the same orf, create annotation that this hasn't been recoded --> add to conflict list
                    else:
                        conflict_cds = conflict_cds + subrangelist
                else:
                    logging.warning("This looks like an error. You are neither following the maximize nor the soft regime for recoding.")

            # if there are too many mane transcripts, only the longest will be taken into account, but this should never happen
            else:
                submanes = []
                for l in submanelist:
                    submanes.append(subrangelist[l])
                orf_shiftlist = [item[2] for item in submanes]
                # if they have the same orf shift
                if (len(set(orf_shiftlist))==1):
                    # if they have also the same end:
                    if len(set([item[1] for item in submanes])):
                        keep_cds.append(submanes[0])
                        logging.info("Too many MANE transcript identified, but they have same length and orf. Appended first transcript {}.".format(submanes[0]))
                    else:
                        cds_id, cds_maxvalue = max(enumerate(map(itemgetter(1), submanes)),key=itemgetter(1))
                        keep_cds.append(submanes[cds_id])
                        logging.info("Several MANE transcripts identified. Appended longest transcript {} with end at {}".format(submanes[cds_id][3],cds_maxvalue))
                # if they don't have the same orf, create annotation that this hasn't been recoded --> add to conflict list
                else:
                    conflict_cds = conflict_cds + submanes                          
        # no overlap --> recode
        else:
            keep_cds.append(rangelist[i])       
    return(keep_cds, conflict_cds)

def createRangeDict(ranges):
    # Every range consist of 0 - start, 1 - end, 2 - phase, 3 - transcript id, 4 - strand
    for_dict = {}
    rev_dict = {}
    for cds in ranges:
        if cds[4] == 1: for_dict[cds [3]] = for_dict.get(cds [3], []) + [[cds[0],cds[1],cds[2]]]
        elif cds[4] == -1: rev_dict[cds [3]] = rev_dict.get(cds [3], []) + [[cds[0],cds[1],cds[2]]]
        else: logging.warning("Problem finding the strand of {} - not included".format(cds))  
    # sort them
    for transcript in for_dict: for_dict[transcript].sort(key=lambda interval: interval[0])
    for transcript in rev_dict: rev_dict[transcript].sort(key=lambda interval: interval[0],reverse = True)
    return for_dict, rev_dict

def find_occurence(find_strings, in_string,protect):
    """ Find startposition of find_strings in in_strings, exclude n = endprotect nt from start and end """    
    occ_pos = []
    occ_list = [y.start() for x in find_strings for y in re.finditer('(?=' + x +')', in_string)]
    occ_list = list(filter(lambda x: (x % 3 == 0), occ_list)) # filter by ORF (divisible by 3)
    for x in occ_list:
        # don't recode codons if it's at the beginning or end
        if x < protect or x >= len(in_string)-protect: continue
        if x >= len(in_string)-3: 
            occ_pos.append([x,1]) # the 1 indicates that it is the last codon in the substring
        else: occ_pos.append([x,0]) # the 0 indicates that it is not the last codon in the substring
    return(occ_pos)

def create_feature_annot(loc_range, featuretype,s):
    """ Create a new feature annotation at loc_range with featuretype on strand s. """
    location = SeqFeature.FeatureLocation(SeqFeature.ExactPosition(loc_range[0]),SeqFeature.ExactPosition(loc_range[1]))
    new_feature = SeqFeature.SeqFeature(location,type=featuretype,strand=s)
    return(new_feature)

def getAnnotRange(record,target = "CDS"):
    """ Get position range for annotation type of interest."""
    annot_range = []
    # Get location ranges of interesting features
    for i in range(len(record.features)):
        # determine which annotations should be taken into account
        if record.features[i].type == target:
            # get their start and end location and phase: [[[start,end,phase,transcript_id, strand]]
            loc_list = [int(record.features[i].location.start),\
                        int(record.features[i].location.end),\
                        int(record.features[i].qualifiers['phase'][0]),\
                        record.features[i].qualifiers['Parent'][0].split(":")[1],\
                        record.features[i].location.strand] 
            annot_range.append(loc_list)
    if len(annot_range) > 0 and annot_range[0] == '': 
        del annot_range[0]
    return(annot_range)