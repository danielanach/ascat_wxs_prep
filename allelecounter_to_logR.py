import pandas as pd
import numpy as np
import argparse

def allele_count_to_logR(tumor_file,
                         normal_file,
                         out_dir,
                         sample,
                         ref_snp_csv):
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
    print('Reading tumor')
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
    print('Reading normal')
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
    print('Merging files')
    merged = tumor.merge(normal,on=['chr','pos'],
                         suffixes=('_tumor','_normal'))

    print('Loading reference SNPs')
    ref_snp_df = pd.read_csv(ref_snp_csv)
    merged = merged.merge(ref_snp_df,how='inner')


    merged['tumor_R'] = np.divide(merged['DP_tumor'],
                                    merged['DP_normal']).replace([np.inf, -np.inf], np.nan)
    merged = merged.dropna(subset=['tumor_R'])
    merged['tumor_logR'] = np.log2(merged['tumor_R']/np.average(merged['tumor_R']))
    merged['normal_logR'] = 0

    merged['chr'] = [c.split('chr')[-1] for c in merged['chr']]

    merged.index = merged['chr'] + '_' + merged['pos'].astype(str)
    merged = merged.sort_values(by=['chr','pos','tumor_logR'])
    merged = merged[~merged.index.duplicated(keep='last')]

    tumor_out_file = merged[['chr','pos','tumor_logR']]
    tumor_out_file.columns = ['chrs','pos',sample]
    tumor_file_name = '{}/{}_Tumor_LogR.txt'.format(out_dir,sample)
    tumor_out_file.to_csv(tumor_file_name,sep='\t')
    print('Wrote tumor file {}'.format(tumor_file_name))

    normal_out_file = merged[['chr','pos','normal_logR']]
    normal_out_file.columns = ['chrs','pos',sample]
    normal_file_name = '{}/{}_Normal_LogR.txt'.format(out_dir,sample)
    normal_out_file.to_csv(normal_file_name ,sep='\t')
    print('Wrote normal file {}'.format(normal_file_name))

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--ref_snp_csv", help="File name for file containing the 1000G SNPs to profile")
parser.add_argument("-t", "--tumor", help="Tumor AlleleCounter File")
parser.add_argument("-n", "--normal", help="Normal AlleleCounter File")
parser.add_argument("-o", "--out_dir", help="Output directory")
parser.add_argument("-s", "--sample", help="Sample name")

args = parser.parse_args()

print( "Sample: {}\nTumor: {}\nNormal: {}\nOutput directory: {}".format(
        args.sample,
        args.tumor,
        args.normal,
        args.out_dir
        ))

allele_count_to_logR(tumor_file=args.tumor,
                     normal_file=args.normal,
                     out_dir=args.out_dir,
                     sample=args.sample,
                     ref_snp_csv = args.ref_snp_csv)
