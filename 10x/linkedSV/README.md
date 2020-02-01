
#### Running the linkedSV pipeline on HG002

#### Running the pipeline


....



#### Comparing with the truthset


Setting the environment
```bash
module load system/Python-3.6.3
module load bioinfo/bcftools-1.9
python -m venv linkedsvenv
source linkedsvenv/bin/activate
pip install -r requirements.txt
```

Reformating linkedSV output

```bash
rm -f linkedSV.vcf.gz
python linkedsv2vcf.py -i phased_possorted_bam.bam.small_deletions.bedpe -o linkedSV.vcf -r ../../../data/genome/hs37d5_hsa10.fa
bgzip linkedSV.vcf
tabix linkedSV.vcf.gz
```



##### Comparing with the truthset
```bash
truthset_dir="../../../truthset"
rm -fr truvari_del
truvari -b $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10_DEL.vcf.gz -c linkedSV.vcf.gz --passonly --includebed $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10.bed -o truvari_del --pctsim 0
```
