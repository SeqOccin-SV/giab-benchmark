

#### Snakemake pipeline for manta detection

```bash
snakemake --configfile config.yaml \
          --cluster-config cluster.yaml \
          --drmaa " --mem={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" \
          --jobs 4 -p -n
```

#### Comparing with the truthset
```bash
truthset_dir="../../truthset"
truvari -b $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10_DEL.vcf.gz -c detection/manta_DEL.vcf.gz --passonly --includebed $truthset_dir/HG002_SVs_Tier1_v0.6_hsa10.bed -o truvari_del --pctsim 0
```
