
### Running NAIBR on HG002

#### Running the pipeline

##### Retrieving alignment from LongRanger

You need to produce the `phased_possorted_bam.bam` file by using the LongRanger pipeline
or you can use the `10x/giabftp/bams/LR-phased_possorted_bam.bam` available in the datadir

##### Setting the environment

```bash
# module load system/Python-2.7.15
# module load bioinfo/bcftools-1.9
# module load bioinfo/samtools-1.9
# python -m venv naibrenv
# source naibrenv/bin/activate
# pip install -r requirements.txt
conda activate /work/project/seqoccin/tools/miniconda/miniconda3/envs/naibr/
git clone https://github.com/raphael-group/NAIBR.git
```

##### Running NAIBR

...

#### Comparing with the truthset

##### Reformating NAIBR output

```bash
rm -f naibr_del*.vcf* noheader.bedpe
tail -n+2 NAIBR_SVs.bedpe > noheader.bedpe
python naibr2vcf.py -i nohead.bedpe.bedpe -o naibr_del.vcf -r ../../../data/genome/hs37d5_hsa10.fa
cat naibr_del.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > naibr_del_sorted.vcf
bgzip naibr_del_sorted.vcf
tabix naibr.vcf.gz
```

##### Running comparison
```bash
truthset_dir="../../../truthset"
rm -rf truvari_del
truvari -b $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10_DEL.vcf.gz -c naibr_del_sorted.vcf.gz --passonly --includebed $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10.bed -o truvari_del --pctsim 0
```
