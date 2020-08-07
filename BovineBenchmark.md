Bovine SV reference set construction
--------------------------------------

#### Available sequence data sets


| Technology   |   runs |  depth of Coverage | 
|:-------------|-----------------:|---------------------:|
| Illumina     |  2 | 120X  | 
| Nanopore     |  5 |  55X |
| PacBio CLR   |  2 |  90X  | 
| PacBio CCS   | 4  | 15X |


#####  Available variants data sets

All vcf files are located under the following dir
```
/work2/project/seqoccin/svdetection/results/Trio2_offspring
```
The directories detailed in the following table are relative to this root dir


| Technology |  caller   | sv dir |
|:----------|:-----------|:------|
| Illumina   | manta      |  illumina_manta/diploid | 
| Illumina   | cnvpipeline  | illumina_cnvpipeline  |
| Nanopore   | svim        |  nanopore_svim |
| Nanopore   | sniffles    |  nanopore_sniffles | 
| PacBio CLR | pbsv        |  pacbio_CLR_pbsv + pacbio_CLR-GP_pbsv | 
| PacBio CCS | pbsv        |  pacbio_CCS_pbsv |
| 10X        | longranger  | 10x_longranger |
| 10X        | linkedSV    | 10x_longranger |
 

_Should the exact file names be given instead ?_

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

## Pooling results

Actually, Trio2 offspring variants were detected with several technologies

1. Illumina technology
2. 10x chromium library
3. Oxford Nanopore (ONT Ultralong read)
4. PacBio CLR
5. PacBio CCS

At the moment, We did not produce sufficient read coverage on PacBio CCS. It is unlikely that new data will arrive in time.

Here are the various methodology that produced results

- Illumina CNV pipeline
- Illumina Manta
- Illumina DeepVariant (check if indels)
- 10X LongRanger
- 10X LinkedSV
- 10X GROCSVs
- ONT Sniffles
- ONT SVIM
- PB CLR pbsv
- PB CCS pbsv

Most of the following test and scripts dev was done in ~/test/VCFModules/Pysam

# Summary

	for f in *.vcf.gz; do bcftools view  -i '(SVTYPE=="DEL" || SVTYPE=="INS" || SVTYPE=="DUP")' $f | bgzip -c > $(basename $f .vcf.gz)_filtered.vcf.gz; done
	# SHOULD I RENAME THE FILE BEFORE END ?
	for f in  *_filtered.vcf.gz; do NAME=$(basename $f _filtered.vcf.gz); zcat $f | awk -F "\t" -v NAME="${NAME}" 'BEGIN {OFS=FS}; $1~/^#/ {print $0; next} {$3=NAME"-"$3; print $0}' | bgzip -c > ${NAME}_filtered_renamed.vcf.gz; done
	for f in *_filtered_renamed.vcf.gz; do echo $f; done > svimmer_full.idl
	for f in *_filtered_renamed.vcf.gz; do tabix $f; done
	svimmer --ids --output svimmer_full.vcf --threads 6 svimmer_full.idl $(seq -s ' ' 1 29) X MT
	python ~/scripts/svimmerParse4RefSet.py svimmer_full.vcf | bgzip -c > svimmer_full_Filtered.vcf.gz
	bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%ID\t%INFO/MERGED_IDS\n' svimmer_full_Filtered.vcf.gz > test.bed
	bedtools cluster -d 1000 -i test.bed > test.bed.cluster
	python ~/seqoccin/scripts/createBed.py test.bed.cluster teloCentro_ARS1.3.bed
	bedtools subtract -a createBed.bed -b teloCentro_ARS1.3.bed > ReferenceSet_v0.2.bed
	python ~/seqoccin/scripts/testPysamVcf.py svimmer_full_Filtered.vcf.gz | bgzip -c > modif.vcf.gz
	zcat modif.vcf.gz | java -jar /home/adf/git_rep/jvarkit/dist/fixvcf.jar | bgzip -c > modif_fix.vcf.gz
	python ~/seqoccin/scripts/filterSVbyID.py SVFiltered.idl modif_fix.vcf.gz | bgzip -c > ReferenceSet_v0.2.vcf.gz
	

## Filtering variants type

	for f in *.vcf.gz; do bcftools view  -i '(SVTYPE=="DEL" || SVTYPE=="INS" || SVTYPE=="DUP")' $f | bgzip -c > $(basename $f .vcf.gz)_filtered.vcf.gz; done

## Renaming variants IDS

	zcat *_filtered.vcf.gz | awk -F "\t" 'BEGIN {OFS=FS}; $1~/^#/ {next} { print $3}' | head
	zcat trio2.offspring.run1.ARS-UCD1.2.minimap2.vcf.gz | awk -F "\t" 'BEGIN {OFS=FS}; $1~/^#/ {print $0; next} {$3=$3"_"NR; print $0}' | bgzip -c > sniffles_fix.vcf.gz
	for f in  *_filtered.vcf.gz; do NAME=$(basename $f _filtered.vcf.gz); zcat $f | awk -F "\t" -v NAME="${NAME}" 'BEGIN {OFS=FS}; $1~/^#/ {print $0; next} {$3=NAME"-"$3; print $0}' | bgzip -c > ${NAME}_filtered_renamed.vcf.gz; done

## Merging with svimmer

### Issues

	for f in *.gz; do tabix $f; done
	for f in *.gz; do echo $f; done > svimmer_full.idl
	svimmer --ids --output svimmer_fullTest.vcf --threads 6 svimmer_full.idl 1 2 3 4 5 6
	svimmer --ids --output svimmer_fullTest.vcf --threads 6 svimmer_full.idl $(seq -s ' ' 1 29) X MT 
	
LinkedSv Results are not even ok for tabix...
svimmer assertion error ... which files have an unknown svtype ? should modify sv.py directly...

	zcat *.vcf.gz | grep -v '#' | perl -nle 'if (m/SVTYPE/) { s/.*?SVTYPE=(.*?);.*/$1/; print;} else { print "NO\n";}' | sort | uniq -c > type.idl
		22695 BND
		627 cnv
		94952 DEL
		11150 DUP
		80734 INS
		 2204 INV
	        3 TRA
		  547 UNK

Lots of stupid svtype. Will filter before hand for DEL,DUP,INS

	zcat ILL_manta.vcf.gz | perl -nle 'if (substr($_,0,1) eq "#" ) { print;} elsif (m/SVTYPE=[DEL|INS|DUP]/) { print;}' | bgzip -c > prout
	bcftools view  -i '(SVTYPE=="DEL" || SVTYPE=="INS" || SVTYPE=="DUP")' ILL_manta.vcf.gz
	for f in *.vcf.gz; do bcftools view  -i '(SVTYPE=="DEL" || SVTYPE=="INS" || SVTYPE=="DUP")' $f | bgzip -c > $(basename $f .vcf.gz)_filtered.vcf.gz; done
	for f in *_filtered.vcf.gz; do tabix $f; done
	for f in *_filtered.vcf.gz; do echo $f; done > svimmer_full_filtered.idl
	svimmer --ids --output svimmer_fullTest.vcf --threads 6 svimmer_full_filtered.idl $(seq -s ' ' 1 29) X MT
	
## Filtering variants for 2+ methods

First filter by methodology
	python ~/scripts/svimmerParse4RefSet.py svimmer_fullTest.vcf | bgzip -c > svimmer_fullTest_Filtered.vcf.gz

Should also filter indels but better to keep them to determine regions

## Filtering regions with overlapping variants

Get the regions
	zcat svimmer_fullTest_Filtered.vcf.gz | grep -v '#' | perl -nle 'my @F = split( /\t/ , $_ ); (my $end = $F[7]) =~ s/.*?END=(\d+)\;.*/$1/; print join("\t",($F[0],$F[1],$end));' > test.bed
	bcftools query -f '%CHROM\t%POS\t%INFO/END\n' svimmer_fullTest_Filtered.vcf.gz > test2.bed
	bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/MERGED_IDS\n' svimmer_fullTest_Filtered.vcf.gz > test2.bed

	bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%ID\t%INFO/MERGED_IDS\n' svimmer_fullTest_Filtered.vcf.gz > test3.bed

Regions with telomeric and centromeric regions
	tail -n+3 ~/Téléchargements/giaa021_supplemental_data/TableS3_centromeric_and_telomeric_repeats.csv | cut -f1-3 | perl -nle 'print unless (m/^\s+$/ || m/^\d+\s+$/)' > teloCentro_ARS1.3.bed

Filtering indels
	bcftools query -f '%INFO/SVLEN\n' svimmer_fullTest_Filtered.vcf.gz | awk '{ if (sqrt($1^2)>49) print $0 }' | wc -l

Clustering on 1k window	
	bedtools cluster -d 1000 -i test2.bed > test2.bed.cluster
	bedtools cluster -d 1000 -i test3.bed > test3.bed.cluster

Defining regions around Good variants within 1k	
	python createBed.py test3.bed.cluster teloCentro_ARS1.3.bed
produce createBed.bed , SVFiltered.idl , SVLenIssue.idl
Removing include regions inside telomeric, centromeric regions
	bedtools subtract -a createBed.bed -b teloCentro_ARS1.3.bed > ReferenceSet_v0.1.bed
	
## Vcf polishing for truvari usage

Obviously stupid truvari depending on Pysam and VCF module cannot read svimmer uninformative VCF
Problems could range from header definition to absence of specific format/sample information

Quick-Fix for Format and sample 
	python testPysamVcf.py svimmer_fullTest_Filtered.vcf.gz | bgzip -c > modif.vcf.gz
	
next issue is about bed file for the construction of intervalTree.
It does not accept bed with same start and end
Done previously

incoherent content between header and INFO field
Testing fixvcf from jvarkit

	zcat modif.vcf.gz | java -jar /home/adf/git_rep/jvarkit/dist/fixvcf.jar | bgzip -c > modif_fix.vcf.gz

Need to filter SV based on createBed.py results
Not what I want, Looking for SV ID filtering, here it is sample filtering
	#bcftools view -S SVFiltered.idl -Oz -o ReferenceSet_v0.1.vcf.gz svimmer_fullTest_Filtered.vcf.gz
	python ./filterSVbyID.py SVFiltered.idl modif_fix.vcf.gz | bgzip -c > ReferenceSet_v0.1.vcf.gz

## Truvari

	truvari -b ReferenceSet_v0.1.vcf.gz -c variants_svimmer_q15.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_test

	truvari -b ReferenceSet_v0.1.vcf.gz -c variants_svimmer_q15.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_ONTSvimq15
	truvari -b ReferenceSet_v0.1.vcf.gz -c trio2.offspring.run1.ARS-UCD1.2.minimap2.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_ONTSniffles
	truvari -b ReferenceSet_v0.1.vcf.gz -c diploidSVconverted.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_ILLmanta
	truvari -b ReferenceSet_v0.1.vcf.gz -c merge_sort.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_TenXLR
	truvari -b ReferenceSet_v0.1.vcf.gz -c Trio2-offspring-CLR_ARS_pbmm2.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_PBCLR
	truvari -b ReferenceSet_v0.1.vcf.gz -c Trio2-offspring-CLR-GetP_ARSonlyChr_pbmm2.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_PDCLRGP
	truvari -b ReferenceSet_v0.1.vcf.gz -c Trio2-offspring-CCS15kb_ARSonlyChr_pbmm2.vcf.gz --passonly --pctsim 0 --includebed ReferenceSet_v0.1.bed -o truvari_Ref01_PBCCS



# Size graph TP

	for f in truvari_Ref01_*; do bcftools query -f '%INFO/SVLEN\n' ${f}/tp-base.vcf > ${f}/size_$(echo $f | sed 's/truvari_Ref01_//').txt; done
	bcftools query -f '%INFO/SVLEN\n' ReferenceSet_v0.1.vcf.gz > size_ReferenceSet.txt
