import pandas as pd

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

    ref_snp_df = pd.read_csv(ref_snp_csv)
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

    annot_df = df.merge(ref_snp_df)

    annot_df['ref_count'] = [annot_df[annot_df['ref'][i]][i] for i in annot_df.index]
    annot_df['alt_count'] = [annot_df[annot_df['alt'][i]][i] for i in annot_df.index]
    annot_df['baf'] = annot_df['alt_count'] / (annot_df['ref_count'] + annot_df['alt_count'])

    out_df = annot_df[['chr','pos','baf']].reset_index(drop=True).dropna()
    out_df = out_df.sort_values(by=['chr','pos','baf'])

    out_df['chr'] = [c.split('chr')[-1] for c in out_df['chr']]
    out_df.index = out_df['chr'] + '_' + out_df['pos'].astype(str)
    out_df.columns = ['chrs','pos',sample]
    out_df = out_df[~out_df.index.duplicated(keep='last')]
    out_df.to_csv('{}/{}_{}_BAF.txt'.format(out_dir,sample,sample_type),
                   sep='\t')
