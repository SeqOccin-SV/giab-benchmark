### Running LongRanger on HG002

#### Running the pipeline

##### Setting reference genome for LongRanger
```bash
dataset_dir="../../../data"
module load bioinfo/longranger-2.2.2
srun -c 1 -p workq -J lr-mkref --export=ALL longranger mkref $dataset_dir/genome/hs37d5_hsa10.fa
```
Adding picard .dict on previously created reference
```bash
module load bioinfo/picard-2.20.7
srun -c 1 -p workq -J PicardLR --export=ALL java -Xmx4g -jar $PICARD CreateSequenceDictionary R=./refdata-hs37d5_hsa10/fasta/genome.fa O=./refdata-hs37d5_hsa10.fa/fasta/genome.dict
```

##### Setting symbolic link with proper Illumina naming
```bash
mkdir data
ln -s $datadir/10x/giabftp/fastq/NA24385_giabftp-noflag_lariat-chr10_1.fq.gz data/NA24385_giabftp_S5_L004_R1_001.fastq.gz
ln -s $datadir/10x/giabftp/fastq/NA24385_giabftp-noflag_lariat-chr10_2.fq.gz data/NA24385_giabftp_S5_L004_R2_001.fastq.gz
```

##### Running LongRanger full pipeline
```bash
create.launch.longranger.wgs.sh NA24385_giabftp refdata-hs37d5_hsa10/ data/ 10 60G
sbatch --mem=60G --cpus-per-task=10 -J lr-wgs_NA24385chr10 --workdir=$PWD --export=ALL -p workq ./launch.longranger.wgs.sh
```
output will be in NA24385_giabftp/outs/
