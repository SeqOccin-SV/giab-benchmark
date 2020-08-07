Bovine SV reference set construction
--------------------------------------

#### Available sequence data sets

Trio2 offspring variants were detected with several technologies

| Technology   |   runs |  depth of Coverage | 
|:-------------|-----------------:|---------------------:|
| Illumina     |  2 | 120X  | 
| 10x Chromium |  1 ?  | 90X |
| Nanopore     |  5 |  55X |
| PacBio CLR   |  2 |  90X  | 
| PacBio CCS   | 4  | 15X |


#####  Available variants data sets

Different methologies that produced the results.  
All vcf files are located under the following dir
```
/work2/project/seqoccin/svdetection/results/Trio2_offspring
```
The vcf files in the following table are relative to this root dir


| Technology |  caller   | vcf files |
|:----------|:-----------|:------|
| Illumina   | manta      |  illumina_manta/diploidSV.vcf.gz| 
| Illumina   | cnvpipeline  | illumina_cnvpipeline/DEL/concat_DEL_final.vcf.gz|
| Nanopore   | svim        |  nanopore_svim/final_results.vcf.gz |
| Nanopore   | sniffles    |  nanopore_sniffles/sniffles_FullvsARS.vcf.gz | 
| PacBio CLR | pbsv        |  pacbio_CLR_pbsv/Trio2-offspring-CLR_ARS_pbmm2.vcf.gz <br> pacbio_CLR-GP_pbsv/Trio2-offspring-CLR-GetP_ARSonlyChr_pbmm2.vcf.gz | 
| PacBio CCS | pbsv        |  pacbio_CCS_pbsv/Trio2-offspring-CCS15kb_ARSonlyChr_pbmm2.vcf.gz |
| 10X        | longranger  | 10x_longranger/dels.vcf.gz|
| 10X        | linkedSV    | ? |

Best policy should be SV discovered by at least 2 methods/technologies.  
For next version, reads subsampling and cross technology local assembly could be a good idea.

# Basic steps

1. Pooling results
2. Renaming variants IDs per method/technology
3. Merging with svimmer (! svtypes INV/DUP/INS !)
4. Filtering variants discover by >= 2 methods
5. Filtering regions with overlapping variants
6. Vcf polishing for truvari usage
7. Statistics comparison on v0.1 set


