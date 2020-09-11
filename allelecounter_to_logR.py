import pandas as pd
import numpy as np

def allele_count_to_logR(tumor_file,
                         normal_file,
                         out_dir,
                         sample):
    ''' Takes allele count data from AlleleCounter and SNP loci and writes
    logR file formatted for use with ASCAT.

    - tumor logR file is computed as:

        tumorR = Tumor Count / Normal Count
        tumor_logR = log2(tumorR/mean(tumorR))

    - normal LogR file is computed as:

        normal_logR = 0

    Parameters
    ----------
    tumor_file : str
        File name for tumor AlleleCounter output counts file.
    normal_file : str
        File name for normal AlleleCounter output counts file.
    out_dir : str
        Output directory for output file for logR values.
    sample_name : str
        Sample name.

    '''
    tumor = pd.read_csv(tumor_file,
                         sep='\t',
                         skiprows=1,
                         names=['chr',
                                'pos',
                                'A',
                                'C',
                                'G',
                                'T',
                                'DP'])

    normal = pd.read_csv(normal_file,
                         sep='\t',
                         skiprows=1,
                         names=['chr',
                                'pos',
                                'A',
                                'C',
                                'G',
                                'T',
                                'DP'])

    merged = tumor.merge(normal,on=['chr','pos'],
                         suffixes=('_tumor','_normal'))

    merged['tumor_R'] = np.divide(merged['DP_tumor'],
                                    merged['DP_normal']).replace([np.inf, -np.inf], np.nan)
    merged = merged.dropna(subset=['tumor_R'])
    merged['tumor_logR'] = np.log2(merged['tumor_R']/np.average(merged['tumor_R']))
    merged['normal_logR'] = 0

    merged.index = merged['chr'] + merged['pos'].astype(str)

    tumor_out_file = merged[['chr','pos','tumor_logR']]
    tumor_out_file.columns = ['chrs','pos',sample]
    tumor_out_file.to_csv('{}/{}_Tumor_LogR.txt'.format(out_dir,sample),sep='\t')

    normal_out_file = merged[['chr','pos','normal_logR']]
    normal_out_file.columns = ['chrs','pos',sample]
    normal_out_file.to_csv('{}/{}_Normal_LogR.txt'.format(out_dir,sample),sep='\t')
