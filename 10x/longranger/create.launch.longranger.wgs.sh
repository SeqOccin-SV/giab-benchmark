#!/bin/bash

if [ $# -ne 5 ]
then
	echo 'USAGE: create.launch.longranger.wgs.sh <id> <ref> <fastq_path> <ncore> <mem>'
else
	echo "#!/bin/bash

cd $PWD
module load bioinfo/GATK-3.8-1-0
module load bioinfo/longranger-2.2.2
longranger wgs --id=$1 --fastq=$3 --reference=$2 --vcmode=gatk:/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar --localcores=$4 --localmem=$5 > longranger-$1.out 2> longranger-$1.err
" > launch.longranger.wgs.$1.sh
fi