
#### Running the manta pipeline on HG002

##### genologin execution

Setting the environment
```bash
conda activate cnvpipeline"
```

Snakemake command

```bash
python cnvpipelines/cnvpipelines.py run detection -r hs37d5.fa -s bamlist.txt -t delly lumpy pindel --chromosomes 10 -w svdetection -p -n
```

##### Comparing with the truthset
```bash
truthset_dir="../../../truthset"
truvari -b $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10_DEL.vcf.gz -c mantasv/manta_DEL.vcf.gz --passonly --includebed $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10.bed -o truvari_del --pctsim 0
```
