##  1.准备数据
### 重命名样品名
#### 由于样品测序文件名中包含LANE等信息，因此需要清洗
rename 's/_.*_1//g' data/*/*.fq.gz
### 1.2 解压
#### 解压样品文件，usearch不识别压缩文件
for i in `ls -d data/*/`;do gzip -d ${i}/*.gz;done
# 准备临时文件夹
mkdir temp
##  2. 合并双端序列并按样品重命名
while read old new
do
	usearch --fastq_mergepairs data/${old}/${old}.R1.fq \
	        --reverse data/${old}/${old}.R2.fq \
	        --fastqout temp/${new}.merged.fq \
	        --relabel ${new}.  # 千万注意最后这个点
done<raw/metadata.tsv
cat temp/*.merged.fq > temp/all.fq
##  3. 切除引物与质控
	# 左边19右边20是338F和806R的长度，并且已确定合并后的文件中没有切除相应序列。
	# ACTCCTACGGGAGGCAGCA 338F
	# GGACTACHVGGGTWTCTAAT  806R
time usearch --fastq_filter temp/all.fq \
	--fastq_stripleft 19 \
	--fastq_stripright 20 \
	--fastq_maxee_rate 0.01 \
	--fastaout temp/filtered.fa
##  4. 序列去冗余与去噪
### 4.1 去除低丰度噪音并增加计算速度
time usearch --fastx_uniques temp/filtered.fa \
	--fastaout temp/uniques.fa \
	--relabel Uni \
	--minuniquesize 10 \
	--sizeout
### 4.2 去噪
### 按97%做的软件就叫聚类，按99%做的软件就叫降噪（因为这个没聚类，只是去掉一些“错”的。unoise3就是降噪
usearch -unoise3 temp/uniques.fa -zotus temp/zotus.fa
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/asv.fa # 将header改成ASV更符合习惯
##  5. 特征表和筛选 Feature table
mkdir -p result/01.raw
cp temp/otus.fa result/01.raw
### 按97%将reads回贴，计算丰度。这步用vsearch更快
vsearch --usearch_global temp/filtered.fa \
        --db result/01.raw/asv.fa \
        --otutabout result/01.raw/asvtab.txt \
        --id 0.97 \
        --threads 8
##  6.使用SINTAX算法注释
### 记得用新版的vsearch
### 由于 SINTAX算法不基于Bayes，因此可以不用特地将V3-V4区域剪切出来
### 6.1 用qiime整理好的转换为sintax friendly的开头
vsearch --sintax result/01.raw/otus.fa \
	--db /home/database/micro/SILVA138_RESCRIPt.fasta \
	--tabbedout result/01.raw/otus.sintax \
	--sintax_cutoff 0.8 # 参考美格，使用0.8
### 6.2 过滤叶绿体、线粒体等序列
Rscript ~/software/EasyMicrobiomeScript/otutab_filter_nonBac.R \
	--input result/01.raw/asvtab.txt \
	--taxonomy result/01.raw/asv.sintax \
	--output result/asvtab.txt \
	--stat result/01.raw/asvtab_nonBac.stat \
	--discard result/01.raw/asv.sintax.discard
### 6.3 只保留过滤后的ASV
cut -f 1 result/asvtab.txt | tail -n+2 > result/asvtab.id
usearch -fastx_getseqs result/01.raw/asv.fa -labels result/asvtab.id -fastaout result/asv.fa
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' \
result/01.raw/asv.sintax result/asvtab.id \
> result/asv.sintax
sed -i 's/\t$/\td:Unassigned/' result/asv.sintax
## 7.构建进化树
### 使用qiime fastme默认参数流程
conda activate qiime2-2021.4
mkdir result/tree
qiime tools import --type FeatureData[Sequence] --input-path result/asv.fa --output-path result/asv.fa.qza
qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences result/asv.fa.qza \
	--o-alignment result/tree/asv.aln.fa \
	--o-masked-alignment result/tree/asv.masked.aln.fa \
	--o-tree result/tree/asv.nwk \
	--o-rooted-tree result/tree/asv.root.nwk \
	--p-n-threads 10
qiime tools export --input-path asv.root.nwk.qza --output-path asv.root.nwk
mv asv.root.nwk/tree.nwk ../asv.rooted.nwk