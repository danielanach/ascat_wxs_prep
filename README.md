# ascat_wxs_prep

Small scripts to prepare ASCAT input files from alleleCount on WXS data 

These scripts run on outputs from allelecount:
https://github.com/cancerit/alleleCount

Prior to using these scripts, use alellecount on BAM files using 1000G SNP loci intersected with your exome panel

## Prepare log R files 

Using allele count output files from paired normal / tumor samples:

```
python allelecounter_to_logR.py --tumor A047-tumor-allelecount.txt --normal A047-normal-allelecount.txt --out_dir output_dir --sample A047
```

## Prepare BAF files 

Using allele count output files separately on paired normal / tumor samples and 1000G bi-allelic SNPs:
(the algorithm is a little slow and may take a few minutes)

```
python allelecounter_to_baf.py --allele_count A047-tumor-allelecount.txt  --ref_snp_csv 1000G_Exome_baSNPs.csv  -o ../ --sample A047 --type Tumor
python allelecounter_to_baf.py --allele_count A047-normal-allelecount.txt  --ref_snp_csv 1000G_Exome_baSNPs.csv  -o ../ --sample A047 --type Normal
```

## Dependencies

Tested with:
* numpy 1.18.3
* pandas 1.0.3
* python 3.7.6 
