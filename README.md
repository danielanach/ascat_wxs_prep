# ascat_wxs_prep

Small scripts to prepare ASCAT input files from alleleCount on WXS data 

These scripts run on outputs from allelecount:
https://github.com/cancerit/alleleCount

Prior to using these scripts, use alellecount on BAM files using 1000G SNP loci intersected with your exome panel

## Prepare log R files 

Using allele count output files from paired normal / tumor samples:

```
python allelecounter_to_logR.py -t A047-tumor-allelecount.txt -n A047-normal-allelecount.txt -o output_dir -s A047
```

## Prepare BAF files 

Using allele count output files separately on paired normal / tumor samples and 1000G bi-allelic SNPs:
(the algorithm is a little slow and may take a few minutes)

```
python allelecounter_to_baf.py -a A047-tumor-allelecount.txt  -r 1000G_Exome_baSNPs.csv  -o ../ -s A047 -t Tumor
python allelecounter_to_baf.py -a A047-normal-allelecount.txt  -r 1000G_Exome_baSNPs.csv  -o ../ -s A047 -t Normal
```
