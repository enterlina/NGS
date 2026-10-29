import numpy as np
import pandas as pd
from IPython.display import display, HTML
import pysam
from Bio import SeqIO
import warnings

warnings.filterwarnings("ignore")


def error_stat(reference, sam_err, sam_cor):
    with open(reference) as reference_handle:
        file = SeqIO.parse(reference_handle, "fasta")
        reference_dict = SeqIO.to_dict(file)
        reference_seq = ''.join([str(chr_data.seq) for chr, chr_data in reference_dict.items()])

        samfile_err = pysam.AlignmentFile(sam_err, "rb")
        samfile_corr = pysam.AlignmentFile(sam_cor, "rb")

        undetected, detected_and_corrected, detected_and_removed, false_corrected, correctly_unmodified, incorrectly_removed = 0, 0, 0, 0, 0, 0
        read_counter = 0
        for raw_read, corr_read in zip(samfile_err, samfile_corr):
            read_counter = read_counter + 1
            if not raw_read.is_unmapped and not corr_read.is_unmapped:
                raw_aligned_bases = raw_read.get_aligned_pairs(with_seq=True)
                raw_dict = dict([(k[0], (k[1], k[2])) for k in raw_aligned_bases])
                corr_aligned_bases = corr_read.get_aligned_pairs(with_seq=True)
                corr_dict = dict([(k[0], (k[1], k[2])) for k in corr_aligned_bases])
                corr_dict2 = dict([(k[1], (k[0], k[2])) for k in corr_aligned_bases])

                for raw_idx, (ref_idx, raw_ref_seq) in raw_dict.items():
                    if raw_idx is None:
                        if ref_idx in corr_dict2.keys():

                            if corr_dict2[ref_idx][1] is None:
                                detected_and_removed += 1
                            elif corr_dict2[ref_idx][1].islower():
                                undetected += 1
                            else:
                                detected_and_corrected += 1
                        else:
                            detected_and_removed += 1
                    else:
                        corr_ref_seq = corr_dict[raw_idx][1]
                        corr_seq = corr_read.seq[raw_idx]

                        if raw_ref_seq is None or raw_ref_seq.islower():
                            if corr_ref_seq is None or corr_seq == 'n' or corr_seq == 'N':
                                detected_and_removed += 1
                            elif corr_ref_seq.islower():
                                undetected += 1
                            else:
                                detected_and_corrected += 1
                        else:
                            if corr_ref_seq is None or corr_seq == 'n' or corr_seq == 'N':
                                incorrectly_removed += 1
                            elif corr_ref_seq.islower():
                                false_corrected += 1
                            else:  # correct in corrected
                                correctly_unmodified += 1
        stats = pd.DataFrame({'Error in corrected reads': [0, 1],
                              'Correct base in corrected reads': [2, 3],
                              'Base is absent in corrected reads': [4, 5]
                              })
        stats = stats.set_index(pd.Series(['Error in raw data', 'Correct base in raw data']))
        stats['Error in corrected reads'][0] = undetected
        stats['Error in corrected reads'][1] = false_corrected

        stats['Correct base in corrected reads'][0] = detected_and_corrected
        stats['Correct base in corrected reads'][1] = correctly_unmodified

        stats['Base is absent in corrected reads'][0] = detected_and_removed
        stats['Base is absent in corrected reads'][1] = incorrectly_removed

    return stats