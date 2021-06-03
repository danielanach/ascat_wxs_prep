import pandas as pd
import argparse

def allele_count_to_BAF(ref_snp_csv,
                        ac_file,
                        out_dir,
                        sample,
                        sample_type):

    ''' Takes allele count data from AlleleCounter and SNP loci and writes
    BAF file formatted for use with ASCAT

    Parameters
    ----------
    ref_snp_csv : str
        File name for file containing the 1000G SNPs to profile.
    ac_file : str
        File name for AlleleCounter output counts file.
    out_dir : str
        Output directory for output file for BAF values.
    sample_name : str
        Sample name.

    '''
    print('Loading reference SNPs')
    ref_snp_df = pd.read_csv(ref_snp_csv)
    print('Loading allelecounts file')
    df = pd.read_csv(ac_file,
                     sep='\t',
                     skiprows=1,
                     names=['chr',
                            'pos',
                            'A',
                            'C',
                            'G',
                            'T',
                            'DP'])
    print('Intersecting reference SNPs and allele counts')
    annot_df = df.merge(ref_snp_df,how='inner')

    print('Computing BAF')
    annot_df['ref_count'] = [annot_df[annot_df['ref'][i]][i] for i in annot_df.index]
    annot_df['alt_count'] = [annot_df[annot_df['alt'][i]][i] for i in annot_df.index]
    annot_df['baf'] = annot_df['alt_count'] / (annot_df['ref_count'] + annot_df['alt_count'])

    out_df = annot_df[['chr','pos','baf']].reset_index(drop=True).dropna()
    out_df = out_df.sort_values(by=['chr','pos','baf'])

    out_df['chr'] = [c.split('chr')[-1] for c in out_df['chr']]
    out_df.index = out_df['chr'] + '_' + out_df['pos'].astype(str)
    out_df.columns = ['chrs','pos',sample]
    out_df = out_df[~out_df.index.duplicated(keep='last')]

    out_file_name = '{}/{}_{}_BAF.txt'.format(out_dir,sample,sample_type)
    print('Writing to file {}'.format(out_file_name))
    out_df.to_csv(out_file_name, sep='\t')

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--ref_snp_csv", help="File name for file containing the 1000G SNPs to profile")
parser.add_argument("-a", "--allele_count", help="AlleleCounter File")
parser.add_argument("-o", "--out_dir", help="Output directory")
parser.add_argument("-s", "--sample", help="Sample name")
parser.add_argument("-t", "--type", help="Sample type")

args = parser.parse_args()

print( "Sample: {}\nSample type: {}\nAlleleCounter file: {}\nRef file: {}\nOutput directory: {}".format(
        args.sample,
        args.type,
        args.allele_count,
        args.ref_snp_csv,
        args.out_dir
        ))

allele_count_to_BAF(ref_snp_csv = args.ref_snp_csv,
                        ac_file=args.allele_count,
                        out_dir=args.out_dir,
                        sample=args.sample,
                        sample_type=args.type)
