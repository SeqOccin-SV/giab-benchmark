### Running PacBio SV detection pipeline on HG002

#### Settin python env with PacBio helper script

# check for module load on genotoul
conda env create -p ./pbsv -f environment.yaml
conda activate ./pbsv/

#### Running the pipeline

##### Using pbmm2 to run minimap2 with preset for CCS
```bash
pbmm2 align {ref} {reads} {output.bam} --sort --preset CCS --sample {sample} --rg '@RG\tID:movie{sample}'
```

##### Running Pacbio Detection
```bash
pbsv discover {output.bam} {output.svsig.gz}
pbsv call {ref} {output.svsig.gz} {output.vcf}
```

