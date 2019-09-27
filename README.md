# UNIX Assignment
Axelle Weeger

## Data Inspection

### Attributes of `fang_et_al_genotypes`

```
wc fang_et_al_genotypes.txt 
du -h fang_et_al_genotypes.txt 
awk -F "\t" '{print NF; exit}' fang_et_al_genotypes.txt 
file fang_et_al_genotypes.txt 
```

By inspecting this file I learned that:

1. I have 2783 lines
2. I have 2744038 words
3. I have 11051939 bytes
4. It is a large 6.1M file
5. I have 986 columns
6. This file is encoded in ASCII, with very long lines


### Attributes of `snp_position`

```
wc snp_position.txt 
du -h snp_position.txt 
awk -F "\t" '{print NF; exit}' snp_position.txt 
file snp_position.txt 
```

By inspecting this file I learned that:

1. I have 984 lines
2. I have 13198 words
3. I have 82763 bytes
4. This file is 38K
5. This file has 15 columns
6. This file is ASCII text


## Data Processing

### Prepping the files for joining


Get maize and teosinte datasets out of `fang_et_al_genotypes` first.

```
grep -E "(ZMMIL|ZMMLR|ZMMMR)" fang_et_al_genotypes.txt > maize_genotype.txt
grep -E "(ZMPBA|ZMPIL|ZMPJA)" fang_et_al_genotypes.txt > teosinte_genotype.txt
```

Make header file from `fang_et_al_genotypes`, \(non transposed\) and add it back to each subset. 

```
head -n 1  fang_et_al_genotypes.txt  >header.txt
cat header.txt maize_genotype.txt > maize_genotypeHead.txt
cat header.txt teosinte_genotype.txt > teosinte_genotypeHead.txt
```

Transpose the subsets with header.

```
awk -f transpose.awk maize_genotypeHead.txt > maize_transposed_genotype.txt
awk -f transpose.awk teosinte_genotypeHead.txt > teosinte_transposed_genotype.txt
```

Make a header from `maize_transposed_genotype` \(it will also be used for `teosinte_transposed_genotype`\), and `snp_position`. 

```
head -n 1 maize_transposed_genotype.txt > fang_transposed_head.txt
head -n 1 snp_position.txt > snp_head.txt
```

We use column 1, the SNP_ID as our matching column for join. First, we sort each file and add the correct header back. 

```
sort -k1,1 maize_transposed_genotype.txt > maize_sorted_transposed_genotype.txt
cat fang_transposed_head.txt maize_sorted_transposed_genotype.txt > maize_stgh.txt

sort -k1,1 teosinte_transposed_genotype.txt > teosinte_sorted_transposed_genotype.txt
cat fang_transposed_head.txt teosinte_sorted_transposed_genotype.txt > teosinte_stgh.txt

sort -k1,1 snp_position.txt > snp_position_sorted.txt   
cat snp_head.txt snp_position_sorted.txt > snp_sh.txt
```

We extract the 3 columns of `snp_position_sorted` that we are interested in. 

```
cut -f 1,3,4 snp_sh.txt > snp_sh_redux.txt
```

### Joining the Maize file

```
awk -F "\t" '{print NF; exit}' maize_stgh.txt
awk -F "\t" '{print NF; exit}' snp_sh.txt
join -1 1 -2 1 -t $'\t' snp_sh_redux.txt maize_stgh.txt --header > maize_snp.txt		
awk -F "\t" '{print NF; exit}' maize_snp.txt
cut -f 1-10 maize_snp.txt | head | column -t
```
Check number of columns of each file. 
Join both files, with option -t to make sure the file is tab delimitted, and --header to ignore header when joining
Check that columns add up, and the files were successfully joined. 
Visualize the file. 


### Joining the Teosinte file

```
awk -F "\t" '{print NF; exit}' teosinte_stgh.txt
awk -F "\t" '{print NF; exit}' snp_sh.txt
join -1 1 -2 1 -t $'\t' snp_sh_redux.txt teosinte_stgh.txt --header > teosinte_snp.txt		
awk -F "\t" '{print NF; exit}' teosinte_snp.txt
cut -f 1-10 teosinte_snp.txt | head | column -t
```

# Files are now joined! 

##M aize Data Extraction

#### Extract unknown and multiple SNP locations. Visualize file for quality control. 

Make header from `maize_snp` file. 

```
head -n 1 maize_snp.txt > maize_snpHead.txt
```

```
awk '$2 ~ /unknown/' maize_snp.txt > maize_snp_unknownHL.txt
cut -f 1-10 maize_snp_unknownHL.txt | head | column -t
cat maize_snpHead.txt maize_snp_unknownHL.txt > maize_snp_unknown.txt

awk '$2 ~ /multiple/' maize_snp.txt > maize_snp_multiHL.txt
cut -f 1-10 maize_snp_multiHL.txt | head | column -t
cat maize_snpHead.txt maize_snp_multiHL.txt > maize_snp_multi.txt
```

#### Extract Chr1 1 data, sorted in ascending order. Add header back onto file. 
```
awk '$2 == 1' maize_snp.txt > maize_snp_01hl.txt
sort -k3 maize_snp_01hl.txt > maize_snp_01hls.txt   
cut -f 1-10 maize_snp_01hls.txt | head | column -t
cat maize_snpHead.txt maize_snp_01hls.txt > maize_snp_01.txt
cut -f 1-10 maize_snp_01.txt | head | column -t
```
#### Sort Chr1 in descending order. Add header back onto file. Replace all instance of '?' with '-'
```
sort -k3 -r maize_snp_01hl.txt > maize_snp_01hlsr.txt   
cut -f 1-10 maize_snp_01hlsr.txt | head | column -t
cat maize_snpHead.txt maize_snp_01hlsr.txt > maize_snp_01r0.txt
cut -f 1-10 maize_snp_01r0.txt | head | column -t
sed 's/?/-/g' maize_snp_01r0.txt > maize_snp_01r.txt
cut -f 1-10 maize_snp_01r.txt | head | column -t
```

#### Repeat file generation for all Chromosomes

Chr2

```
awk '$2 == 2' maize_snp.txt > maize_snp_02hl.txt | sort -k3 maize_snp_02hl.txt > maize_snp_02hls.txt | cat maize_snpHead.txt maize_snp_02hls.txt > maize_snp_02.txt

sort -k3 -r maize_snp_02hl.txt > maize_snp_02hlsr.txt | cat maize_snpHead.txt maize_snp_02hlsr.txt > maize_snp_02r0.txt | sed 's/?/-/g' maize_snp_02r0.txt > maize_snp_02r.txt
```

Chr3

```
awk '$2 == 3' maize_snp.txt > maize_snp_03hl.txt | sort -k3 maize_snp_03hl.txt > maize_snp_03hls.txt | cat maize_snpHead.txt maize_snp_03hls.txt > maize_snp_03.txt

sort -k3 -r maize_snp_03hl.txt > maize_snp_03hlsr.txt | cat maize_snpHead.txt maize_snp_03hlsr.txt > maize_snp_03r0.txt | sed 's/?/-/g' maize_snp_03r0.txt > maize_snp_03r.txt

```

Chr4

```
awk '$2 == 4' maize_snp.txt > maize_snp_04hl.txt | sort -k3 maize_snp_04hl.txt > maize_snp_04hls.txt | cat maize_snpHead.txt maize_snp_04hls.txt > maize_snp_04.txt

sort -k3 -r maize_snp_04hl.txt > maize_snp_04hlsr.txt | cat maize_snpHead.txt maize_snp_04hlsr.txt > maize_snp_04r0.txt | sed 's/?/-/g' maize_snp_04r0.txt > maize_snp_04r.txt
```

Chr5

```
awk '$2 == 5' maize_snp.txt > maize_snp_05hl.txt | sort -k3 maize_snp_05hl.txt > maize_snp_05hls.txt | cat maize_snpHead.txt maize_snp_05hls.txt > maize_snp_05.txt

sort -k3 -r maize_snp_05hl.txt > maize_snp_05hlsr.txt | cat maize_snpHead.txt maize_snp_05hlsr.txt > maize_snp_05r0.txt | sed 's/?/-/g' maize_snp_05r0.txt > maize_snp_05r.txt

```

Chr6

```
awk '$2 == 6' maize_snp.txt > maize_snp_06hl.txt | sort -k3 maize_snp_06hl.txt > maize_snp_06hls.txt | cat maize_snpHead.txt maize_snp_06hls.txt > maize_snp_06.txt

sort -k3 -r maize_snp_06hl.txt > maize_snp_06hlsr.txt | cat maize_snpHead.txt maize_snp_06hlsr.txt > maize_snp_06r0.txt | sed 's/?/-/g' maize_snp_06r0.txt > maize_snp_06r.txt

```

Chr7

```
awk '$2 == 7' maize_snp.txt > maize_snp_07hl.txt | sort -k3 maize_snp_07hl.txt > maize_snp_07hls.txt | cat maize_snpHead.txt maize_snp_07hls.txt > maize_snp_07.txt

sort -k3 -r maize_snp_07hl.txt > maize_snp_07hlsr.txt | cat maize_snpHead.txt maize_snp_07hlsr.txt > maize_snp_07r0.txt | sed 's/?/-/g' maize_snp_07r0.txt > maize_snp_07r.txt

```

Chr8

```
awk '$2 == 8' maize_snp.txt > maize_snp_08hl.txt | sort -k3 maize_snp_08hl.txt > maize_snp_08hls.txt | cat maize_snpHead.txt maize_snp_08hls.txt > maize_snp_08.txt

sort -k3 -r maize_snp_08hl.txt > maize_snp_08hlsr.txt | cat maize_snpHead.txt maize_snp_08hlsr.txt > maize_snp_08r0.txt | sed 's/?/-/g' maize_snp_08r0.txt > maize_snp_08r.txt

```

Chr9

```
awk '$2 == 9' maize_snp.txt > maize_snp_09hl.txt | sort -k3 maize_snp_09hl.txt > maize_snp_09hls.txt | cat maize_snpHead.txt maize_snp_09hls.txt > maize_snp_09.txt

sort -k3 -r maize_snp_09hl.txt > maize_snp_09hlsr.txt | cat maize_snpHead.txt maize_snp_09hlsr.txt > maize_snp_09r0.txt | sed 's/?/-/g' maize_snp_09r0.txt > maize_snp_09r.txt

```

Chr10

```
awk '$2 == 10' maize_snp.txt > maize_snp_010hl.txt | sort -k3 maize_snp_010hl.txt > maize_snp_010hls.txt | cat maize_snpHead.txt maize_snp_010hls.txt > maize_snp_010.txt

sort -k3 -r maize_snp_010hl.txt > maize_snp_010hlsr.txt | cat maize_snpHead.txt maize_snp_010hlsr.txt > maize_snp_010r0.txt | sed 's/?/-/g' maize_snp_010r0.txt > maize_snp_010r.txt


```

### Teosinte Data Extraction


#### Extract unknown and multiple SNP locations. Visualize file for quality control. 

Make header from `teosinte_snp` file

```
head -n 1 teosinte_snp.txt > teosinte_snpHead.txt
```

```
awk '$2 ~ /unknown/' teosinte_snp.txt > teosinte_snp_unknownHL.txt | cat teosinte_snpHead.txt teosinte_snp_unknownHL.txt > teosinte_snp_unknown.txt
cut -f 1-10 teosinte_snp_unknown.txt | head | column -t

awk '$2 ~ /multiple/' teosinte_snp.txt > teosinte_snp_multiHL.txt | cat teosinte_snpHead.txt teosinte_snp_multiHL.txt> teosinte_snp_multi.txt
cut -f 1-10 teosinte_snp_multi.txt | head | column -t
```

#### Extract Chromosome 1 data, sorted in ascending order. Add header back onto file. 
```
awk '$2 == 1' teosinte_snp.txt > teosinte_snp_01hl.txt 
sort -k3 teosinte_snp_01hl.txt > teosinte_snp_01hls.txt   
cut -f 1-10 teosinte_snp_01hls.txt | head | column -t
cat teosinte_snpHead.txt teosinte_snp_01hls.txt > teosinte_snp_01.txt
cut -f 1-10 teosinte_snp_01.txt | head | column -t
```
#### Sort Chr1 in descending order. Add header back onto file. replace all instance of '?' with '-'
```

sort -k3 -r teosinte_snp_01hl.txt > teosinte_snp_01hlsr.txt   
cut -f 1-10 teosinte_snp_01hlsr.txt | head | column -t
cat teosinte_snpHead.txt teosinte_snp_01hlsr.txt > teosinte_snp_01r0.txt
cut -f 1-10 teosinte_snp_01r0.txt | head | column -t
sed 's/?/-/g' teosinte_snp_01r0.txt > teosinte_snp_01r.txt
cut -f 1-10 teosinte_snp_01r.txt | head | column -t
```

#### Repeat file generation for all Chromosomes

Chr2

```
awk '$2 == 2' teosinte_snp.txt > teosinte_snp_02hl.txt | sort -k3 teosinte_snp_02hl.txt > teosinte_snp_02hls.txt | cat teosinte_snpHead.txt teosinte_snp_02hls.txt > teosinte_snp_02.txt

sort -k3 -r teosinte_snp_02hl.txt > teosinte_snp_02hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_02hlsr.txt > teosinte_snp_02r0.txt | sed 's/?/-/g' teosinte_snp_02r0.txt > teosinte_snp_02r.txt
```

Chr3

```
awk '$2 == 3' teosinte_snp.txt > teosinte_snp_03hl.txt | sort -k3 teosinte_snp_03hl.txt > teosinte_snp_03hls.txt | cat teosinte_snpHead.txt teosinte_snp_03hls.txt > teosinte_snp_03.txt

sort -k3 -r teosinte_snp_03hl.txt > teosinte_snp_03hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_03hlsr.txt > teosinte_snp_03r0.txt | sed 's/?/-/g' teosinte_snp_03r0.txt > teosinte_snp_03r.txt

```

Chr4

```
awk '$2 == 4' teosinte_snp.txt > teosinte_snp_04hl.txt | sort -k3 teosinte_snp_04hl.txt > teosinte_snp_04hls.txt | cat teosinte_snpHead.txt teosinte_snp_04hls.txt > teosinte_snp_04.txt

sort -k3 -r teosinte_snp_04hl.txt > teosinte_snp_04hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_04hlsr.txt > teosinte_snp_04r0.txt | sed 's/?/-/g' teosinte_snp_04r0.txt > teosinte_snp_04r.txt

```

Chr5

```
awk '$2 == 5' teosinte_snp.txt > teosinte_snp_05hl.txt | sort -k3 teosinte_snp_05hl.txt > teosinte_snp_05hls.txt | cat teosinte_snpHead.txt teosinte_snp_05hls.txt > teosinte_snp_05.txt

sort -k3 -r teosinte_snp_05hl.txt > teosinte_snp_05hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_05hlsr.txt > teosinte_snp_05r0.txt | sed 's/?/-/g' teosinte_snp_05r0.txt > teosinte_snp_05r.txt

```

Chr6

```
awk '$2 == 6' teosinte_snp.txt > teosinte_snp_06hl.txt | sort -k3 teosinte_snp_06hl.txt > teosinte_snp_06hls.txt | cat teosinte_snpHead.txt teosinte_snp_06hls.txt > teosinte_snp_06.txt

sort -k3 -r teosinte_snp_06hl.txt > teosinte_snp_06hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_06hlsr.txt > teosinte_snp_06r0.txt | sed 's/?/-/g' teosinte_snp_06r0.txt > teosinte_snp_06r.txt

```

Chr7

```
awk '$2 == 7' teosinte_snp.txt > teosinte_snp_07hl.txt | sort -k3 teosinte_snp_07hl.txt > teosinte_snp_07hls.txt | cat teosinte_snpHead.txt teosinte_snp_07hls.txt > teosinte_snp_07.txt

sort -k3 -r teosinte_snp_07hl.txt > teosinte_snp_07hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_07hlsr.txt > teosinte_snp_07r0.txt | sed 's/?/-/g' teosinte_snp_07r0.txt > teosinte_snp_07r.txt

```

Chr8

```
awk '$2 == 8' teosinte_snp.txt > teosinte_snp_08hl.txt | sort -k3 teosinte_snp_08hl.txt > teosinte_snp_08hls.txt | cat teosinte_snpHead.txt teosinte_snp_08hls.txt > teosinte_snp_08.txt

sort -k3 -r teosinte_snp_08hl.txt > teosinte_snp_08hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_08hlsr.txt > teosinte_snp_08r0.txt | sed 's/?/-/g' teosinte_snp_08r0.txt > teosinte_snp_08r.txt

```

Chr9

```
awk '$2 == 9' teosinte_snp.txt > teosinte_snp_09hl.txt | sort -k3 teosinte_snp_09hl.txt > teosinte_snp_09hls.txt | cat teosinte_snpHead.txt teosinte_snp_09hls.txt > teosinte_snp_09.txt

sort -k3 -r teosinte_snp_09hl.txt > teosinte_snp_09hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_09hlsr.txt > teosinte_snp_09r0.txt | sed 's/?/-/g' teosinte_snp_09r0.txt > teosinte_snp_09r.txt

```

Chr10

```
awk '$2 == 10' teosinte_snp.txt > teosinte_snp_010hl.txt | sort -k3 teosinte_snp_010hl.txt > teosinte_snp_010hls.txt | cat teosinte_snpHead.txt teosinte_snp_010hls.txt > teosinte_snp_010.txt

sort -k3 -r teosinte_snp_010hl.txt > teosinte_snp_010hlsr.txt | cat teosinte_snpHead.txt teosinte_snp_010hlsr.txt > teosinte_snp_010r0.txt | sed 's/?/-/g' teosinte_snp_010r0.txt > teosinte_snp_010r.tx
```
