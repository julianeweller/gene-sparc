import pandas as pd
import os
from Bio import SeqFeature
import logging

def export_primedesign_format(record, replace_map_for,replace_map_rev, outdir, outprefix):     
    complement_dict = {'A':'T','T':'A','C':'G','G':'C'}
    # output file for PrimeDesign
    rev_brackets = {}
    non_brackets = {}
    for x in replace_map_for.keys():
        codon = x
        print("codon: ",codon)
        for pos in replace_map_for[x]:
            primedesign_id = f"for_{x}_{pos[0]}-{pos[1]}-{pos[2]}"
            codon_origin = record.seq[pos[0]] + record.seq[pos[1]] + record.seq[pos[2]]
            print("for: ",codon_origin)
            # if all three nucleotides are consecutive
            if (pos[0] == pos[1]-1) and (pos[1] == pos[2]-1):
                if (codon_origin[2] == codon[2]) and ((codon_origin[1] == codon[1]) and (codon_origin[0] == codon[0])):
                    # in this case, no gene editing is required
                    continue
                elif codon_origin[0] == codon[0]:
                    if codon_origin[1] == codon[1]:
                        bracket = f"{codon_origin[0]}{codon_origin[1]}({codon_origin[2]}/{codon[2]})"
                    elif codon_origin[2] == codon[2]:
                        bracket = f"{codon_origin[0]}({codon_origin[1]}/{codon[1]}){codon_origin[2]}"
                    else:
                        bracket = f"{codon_origin[0]}({codon_origin[1]}{codon_origin[2]}/{codon[1]}{codon[2]})"
                elif (codon_origin[1] == codon[1]) and (codon_origin[2] == codon[2]):
                    bracket = f"({codon_origin[0]}/{codon[0]}){codon_origin[1]}{codon_origin[2]}"
                else:
                    bracket = f"({codon_origin}/{codon})"
                primedesign_seq = str(record.seq[(pos[0]-30):pos[0]]) + bracket + str(record.seq[(pos[2]+1):pos[2]+31])
                rev_brackets[primedesign_id] = primedesign_seq
            else:
                print(pos)
                print(f"codon: {codon}, codon_origin {codon_origin}")
                # in this case, the codon spreads across a splicing site. Check if the edit can be made with an individual pegRNA
                if (pos[0] == pos[1] - 1):
                    # first two are consecutive --> third nt needs to be the same to make it happen
                    if codon_origin[2] == codon[2]:
                        bracket = f"({codon_origin[0]}{codon_origin[1]}/{codon[0]}{codon[1]})"
                        primedesign_seq = str(record.seq[(pos[2]-30):pos[2]]) + bracket + str(record.seq[(pos[1]+1):pos[1]+32])
                        rev_brackets[primedesign_id] = primedesign_seq
                    else:
                        non_brackets[primedesign_id] = [x,codon_origin, pos]
                else:
                    # second two are consecutive --> first nt needs to be the same to make it happen
                    if codon_origin[0] == codon[0]:
                        bracket = f"({codon_origin[1]}{codon_origin[2]}/{codon[1]}{codon[2]})"
                        primedesign_seq = str(record.seq[(pos[1]-31):pos[1]]) + bracket + str(record.seq[(pos[2]+1):pos[2]+31])
                        rev_brackets[primedesign_id] = primedesign_seq
                    else:
                        non_brackets[primedesign_id] = [x,codon_origin, pos]

    for x in replace_map_rev.keys():
        codon = [complement_dict[y] for y in list(x)]
        codon = "".join(str(x) for x in codon)
        print("codon:", codon)
        for pos in replace_map_rev[x]:
            primedesign_id = f"rev_{x}_{pos[2]}-{pos[1]}-{pos[0]}"
            codon_origin = [record.seq[pos[2]], record.seq[pos[1]],record.seq[pos[0]]]
            codon_origin = "".join(str(x) for x in codon_origin)
            print("rev: ", codon_origin)
            # if all three nucleotides are consecutive
            if (pos[0] == pos[1] + 1) and (pos[1] == pos[2]+1):
                if (codon_origin[2] == codon[2]) and ((codon_origin[1] == codon[1]) and (codon_origin[0] == codon[0])):
                    # in this case, no gene editing is required
                    continue
                elif codon_origin[2] == codon[2]:
                    if codon_origin[1] == codon[1]:
                        bracket = f"{codon_origin[2]}{codon_origin[1]}({codon_origin[0]}/{codon[0]})"
                    elif codon_origin[0] == codon[0]:
                        bracket = f"{codon_origin[2]}({codon_origin[1]}/{codon[1]}){codon_origin[0]}"
                    else:
                        bracket = f"{codon_origin[2]}({codon_origin[1]}{codon_origin[0]}/{codon[1]}{codon[0]})"
                elif (codon_origin[1] == codon[1]) and (codon_origin[0] == codon[0]):
                    bracket = f"({codon_origin[2]}/{codon[2]}){codon_origin[1]}{codon_origin[0]}"
                else:
                    bracket = f"({codon_origin}/{codon})"
                primedesign_seq = str(record.seq[(pos[2]-30):pos[2]]) + bracket + str(record.seq[(pos[0]+1):pos[0]+31])
                rev_brackets[primedesign_id] = primedesign_seq
            else:
                print(pos)
                print(f"codon: {codon}, codon_origin {codon_origin}")
                # in this case, the codon spreads across a splicing site. Check if the edit can be made with an individual pegRNA
                if (pos[1] == pos[2] + 1):
                    # first two are consecutive --> third nt needs to be the same to make it happen
                    if codon_origin[0] == codon[0]:
                        bracket = f"({codon_origin[2]}{codon_origin[1]}/{codon[2]}{codon[1]})"
                        primedesign_seq = str(record.seq[(pos[2]-30):pos[2]]) + bracket + str(record.seq[(pos[1]+1):pos[1]+32])
                        rev_brackets[primedesign_id] = primedesign_seq
                    else:
                        non_brackets[primedesign_id] = [x, pos]
                else:
                    # second two are consecutive --> first nt needs to be the same to make it happen
                    if codon_origin[0] == codon[0]:
                        bracket = f"({codon_origin[2]}{codon_origin[1]}/{codon[2]}{codon[1]})"
                        primedesign_seq = str(record.seq[(pos[2]-31):pos[2]]) + bracket + str(record.seq[(pos[1]+1):pos[1]+31])
                        rev_brackets[primedesign_id] = primedesign_seq
                    else:
                        non_brackets[primedesign_id] = [x, pos]
    
    
    pd.DataFrame.from_dict(rev_brackets, orient="index").to_csv(os.path.join(outdir, outprefix + "_primedesign.csv"))
    pd.DataFrame.from_dict(replace_map_for, orient="index").to_csv(os.path.join(outdir, outprefix + "_primedesign_replacemap_for.csv"))

    if non_brackets:
        failedpath = os.path.join(outdir, outprefix + "_primedesign_failed.csv")
        pd.DataFrame.from_dict(non_brackets, orient="index").to_csv(failedpath)
        logging.info(f"Edits that could not be prepared for prime editing are saved in {failedpath}")
