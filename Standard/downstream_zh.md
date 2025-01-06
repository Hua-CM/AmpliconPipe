# 16S下游分析流程

## 内容

本流程包括了如下常规分析内容

> [!Important]
>
> 目前该流程只支持单分组分析

- 稀释曲线
- α多样性
- β多样性
- 物种组成
- 核心菌群
- 物种差异分析
- 网络分析



## 使用方法

使用R studio打开rmd文件，修改其中如下所示的输入文件路径（**必要时可能需要手动处理这些文件，以使最终读取的数据符合输入对象的要求**），然后运行每个模块。



## 输入

```{r set_data}
table_sam_path <- "sample.meta.tsv"
table_tax_path <- "asv.sintax"
table_otu_path <- "asvtab.tsv"
tree_path <- "tree.root.nwk"
output_dir_path <- "demo16S"

# Analysis settings
group1 = "group"
order1 = c('CK', 'Soil', 'GE', 'P', 'PGE')
color1 =  color_scheme("Plan2")
```



### 文件路径

- `table_sam_path` 样本元信息表路径，第一列为样品名
- `table_tax_path` OTU物种分类信息文件路径。默认使用**sintax分类算法默认的输出文件**
- `table_otu_path` OTU丰度表路径。行为OTU名，列为样品名。
- `tree_path` 进化树路径。

- `output_dir_path`：输出文件夹路径。



### 分析设置

- `group1`: 分组名，table_sample中的一个列名。例如`"group"`
- `order1`: （ 非必选项）画图时展示的分组顺序。例如` c('CK', 'Soil', 'GE', 'P', 'PGE')`
- `color1`: （ 非必选项）每个分组的颜色，和`order1`中记录的顺序相对应。使用microeco中的默认配色。



### 输入对象

根据上述路径读取数据之后，会将其转换为数据对象。这是直接用于构建分析对象的输入对象。**因此，如果输入文件格式和上述不同的话，手动处理数据，保证以下输入对象的格式符合下述描述即可**

- `table_otu`：行名为*ASV*，列名为*样品名*。 data.frame对象
- `table_sample`：行名为*样品名*，列名为*变量名*。
- `table_tax`：行名为*ASV*，列名为*分类层级名*。包含，且只能包含："Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
- `tree_otu`：使用`ape::read.tree`读取的nwk格式的树



## 输出

输出文件夹目录结构如下：

~~~shell
Output
├─01_rarefication
│      rarefication.pdf
│      rarefied_otu_table.tsv
│
├─02_Alpha
│      alpha_group_summary_result.tsv
│      alpha_result.tsv
│      Observed.pdf
│      Pielou.pdf
│      Shannon.pdf
│
├─03_Beta
│  ├─bray
│  │      ANOSIM.tsv
│  │      NMDS.pdf
│  │      NMDS.tsv
│  │      PCoA.pdf
│  │      PCoA.tsv
│  │      PERANOVA.tsv
│  ├─jaccard
│  ├─unwei_unifrac
│  └─wei_unifrac
│
├─04_Composition
│      Genus_abun.pdf
│      Genus_abun.tsv
│      Genus_group_abun.pdf
│      Genus_group_abun.tsv
│      petal.pdf
│      Phylum_abun.pdf
│      Phylum_abun.tsv
│      Phylum_group_abun.pdf
│      Phylum_group_abun.tsv
│      ASV_category.tsv
│
├─05_Diff
│      LEfSe.tsv
│      LEfSe.pdf
│
└─06_Network
        edge.tsv
        key_node.pdf
        network.gexf
        network_attribute.tsv
        node.tsv
~~~

`01_rarefication`：抽平相关结果

- `rarefication.pdf`稀释曲线图

- `rarefied_otu_table.tsv` 抽平后的OTU丰度表。**后续分析均基于此表进行。**

`02_Alpha`：Alpha多样性相关结果

- `alpha_group_summary_result.tsv`：分组Alpha多样性指标统计值
- `alpha_result.tsv`：各样品Alpha多样性指标值
- `Observed.pdf` ：观测物种数分组柱状图
- ``Pielou.pdf` ：Pielou指数分组柱状图
- `Shannon.pdf` ：Shannon指数分组柱状图

`03_Beta`： Beta多样性相关结果

该文件夹下分为四个子文件夹，分别是基于`Bray-Curtis`、`jaccard`、`unweighted unifrac`、`weighted unifrac`矩阵计算的Beta多样性结果。每个子文件下的文件一样，各文件的内容如下：

- `ANOSIM.tsv`:ANOSIM检验结果
- `NMDS.pdf`: NMDS分析图
- `NMDS.tsv` : 样品-NMDS值矩阵
- `PCoA.pdf`: 主坐标分析图
- `PCoA.tsv`: 样品-主坐标值矩阵
- `PERMANOVA.tsv`: PERMANOVA检验结果

`04_Composition` ：样品物种组成相关结果

- `Genus_abun.pdf`：属水平相对丰度图

- `Genus_abun.tsv`：属水平相对丰度表
- `Genus_group_abun.pdf`：属水平分组相对丰度图
- `Genus_group_abun.tsv`：属水平分组相对丰度表
- `venn.pdf` 维恩图。**当分组小于5组时会输出维恩图**
- `petal.pdf`：花瓣图
- `Phylum_abun.pdf`：门水平相对丰度图
- `Phylum_abun.tsv`：门水平相对丰度表
- `Phylum_group_abun.pdf`：门水平分组相对丰度图
- `Phylum_group_abun.tsv`：门水平分组相对丰度表
- `ASV_category.tsv`：各分组间共有ASV结果

> [!Note]
>
> 有关维恩图和花瓣图的相关计算标准，请参阅关键方法参数

`05_Diff`：差异物种检验结果

- `LEfSe.tsv`：LEfSe判别结果
- `LEfSe.pdf` ：LEfSe柱状图

`06_Network`：网络分析结果

- `edge.tsv`：边文件
- `key_node.pdf`：关键节点图
- `network.gexf`：网络文件（Gephi软件使用）
- `network_attribute.tsv`：网络属性统计结果
- `node.tsv`：节点文件



## 关键方法参数

> [!Important] 
>
> microeco中所有函数使用相对丰度筛选原则：将每个ASV在所有样本中的reads数相加，再除以所有样本的reads总和。这种方法非常不合理。

1. alpha多样性的分组比较使用Kruskal-Wallis检验配合Dunn后置检验进行

2. 各组之间的共有ASV的筛选过程如下：

   step1：先筛选各组的ASV：在每组50%样品中出现，且在每个样品出现样品中的相对丰度都高于0.02%的ASV，被认为是该组的ASV。

   step2：使用每组的ASV计算交集，获得各组之间的共有ASV。

3. 用于构建网络的ASV：

   step1：先筛选各组的ASV：在每组50%样品中出现，且在每个样品出现样品中的相对丰度都高于0.02%的ASV，被认为是该组的ASV。

   step2：进一步筛选出在至少50%的全部样本中出现的ASV。
