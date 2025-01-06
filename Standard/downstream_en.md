# 16S Downstream Analysis Pipeline

## Content

This pipeline includes the following routine analysis contents:

> [!Important]
>
> Currently, this pipeline only supports single group analysis.

- **Rarefaction Curve**
- **Alpha Diversity**
- **Beta Diversity**
- **Species Composition**
- **Core Microbiome**
- **Species Differential Analysis**
- **Network Analysis**

## Usage

Open the Rmd file in R studio, modify the input file paths as shown below (it may be necessary to manually process these files to ensure that the final read data meets the requirements of the input objects), and then run each module.

## Input

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

### File Paths

- `table_sam_path`: Path to the sample metadata table, with the first column as the sample name.
- `table_tax_path`: Path to the OTU species classification information file. Default output file of the **sintax classification algorithm**.
- `table_otu_path`: Path to the OTU abundance table. Rows are OTU names, columns are sample names.
- `tree_path`: Path to the evolutionary tree.
- `output_dir_path`: Path to the output directory.

### Analysis Settings

- `group1`: Group name, a column name in `table_sample`. For example, `"group"`.
- `order1`: (Optional) Order of groups to be displayed in the plot. For example, `c('CK', 'Soil', 'GE', 'P', 'PGE')`.
- `color1`: (Optional) Colors for each group, corresponding to the order in `order1`. Use the default color scheme in `microeco`.

### Input Objects

After reading the data according to the above paths, it will be converted into data objects. These are the input objects directly used to build the analysis objects. Therefore, if the input file format is different from the above, manually process the data to ensure that the format of the following input objects meets the following description:

- `table_otu`: Row names are *ASV*, column names are *sample names*. `data.frame` object.
- `table_sample`: Row names are *sample names*, column names are *variable names*.
- `table_tax`: Row names are *ASV*, column names are *classification level names*. Contains, and only contains: "Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species"
- `tree_otu`: A tree read using `ape::read.tree` in nwk format.

## Output

The output directory structure is as follows:

shell

```shell
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
```

`01_rarefication`: Rarefaction related results

- `rarefication.pdf`: Rarefaction curve plot.
- `rarefied_otu_table.tsv`: Rarefied OTU abundance table. **Subsequent analyses are based on this table.**

`02_Alpha`: Alpha diversity related results

- `alpha_group_summary_result.tsv`: Group summary of alpha diversity indices.
- `alpha_result.tsv`: Alpha diversity indices for each sample.
- `Observed.pdf`: Bar plot of observed species numbers by group.
- `Pielou.pdf`: Bar plot of Pielou's index by group.
- `Shannon.pdf`: Bar plot of Shannon's index by group.

`03_Beta`: Beta diversity related results

This directory is divided into four subdirectories, calculated based on the `Bray-Curtis`, `jaccard`, `unweighted unifrac`, and `weighted unifrac` matrices. The files in each subdirectory are the same, and the contents are as follows:

- `ANOSIM.tsv`: ANOSIM test results.
- `NMDS.pdf`: NMDS analysis plot.
- `NMDS.tsv`: Sample-NMDS value matrix.
- `PCoA.pdf`: Principal coordinate analysis plot.
- `PCoA.tsv`: Sample-principal coordinate value matrix.
- `PERMANOVA.tsv`: PERMANOVA test results.

`04_Composition`: Sample species composition related results

- `Genus_abun.pdf`: Relative abundance plot at the genus level.
- `Genus_abun.tsv`: Relative abundance table at the genus level.
- `Genus_group_abun.pdf`: Relative abundance plot at the genus level by group.
- `Genus_group_abun.tsv`: Relative abundance table at the genus level by group.
- `venn.pdf`: Venn diagram. **When the number of groups is less than 5, a Venn diagram will be output.**
- `petal.pdf`: Petal plot.
- `Phylum_abun.pdf`: Relative abundance plot at the phylum level.
- `Phylum_abun.tsv`: Relative abundance table at the phylum level.
- `Phylum_group_abun.pdf`: Relative abundance plot at the phylum level by group.
- `Phylum_group_abun.tsv`: Relative abundance table at the phylum level by group.
- `ASV_category.tsv`: ASV results common to each group.

> [!Note]
>
> For the calculation standards of Venn diagrams and petal plots, please refer to the key method parameters.

`05_Diff`: Differential species test results

- `LEfSe.tsv`: LEfSe discrimination results.
- `LEfSe.pdf`: LEfSe bar plot.

`06_Network`: Network analysis results

- `edge.tsv`: Edge file.
- `key_node.pdf`: Key node plot.
- `network.gexf`: Network file (for use with Gephi software).
- `network_attribute.tsv`: Network attribute statistics results.
- `node.tsv`: Node file.



## Key Method Parameters

> [!Important]
>
> In `microeco`, all functions use relative abundance filtering principles: sum the reads of each ASV across all samples, then divide by the total reads of all samples. This method is very unreasonable.

1. **Group comparison of alpha diversity** uses the Kruskal-Wallis test with Dunn post-hoc test.

2. **Selection of ASV common to each group** is as follows:

   Step 1: First, select ASV for each group: ASV that appears in 50% of the samples in each group and has a relative abundance greater than 0.02% in each sample is considered ASV of that group.

   Step 2: Calculate the intersection using the ASV of each group to obtain ASV common to each group.

3. **ASV for network construction**:

   Step 1: First, select ASV for each group: ASV that appears in 50% of the samples in each group and has a relative abundance greater than 0.02% in each sample is considered ASV of that group.

   Step 2: Further select ASV that appear in at least 50% of all samples.

