# Format of input data for each step 
## 1. Multi-omics Integration by MixOmics 
### 1.1. Nomalised data :shipit:

```
head(Strep_trans_data)
```
|             |   g000004   |    g000007   |   g000009_1  | ...          | 
| ----------- | ----------- |  ----------- |  ----------- | -----------  | 
| 5448_1_RPMI |  12.70983   |   8.779187   |   6.675756   | ...          | 
| 5448_1_SERA |  12.30743   |   8.512199   |   7.697677   | ...          | 
| 5448_2_RPMI |  12.32878   |   8.786838   |   7.606702   | ...          | 
| 5448_2_SERA |  11.49004   |   9.012668   |   8.913019   | ...          | 
| ...         |  ...        |  ...         |  ...         | ...          | 

```
head(Strep_prot_data)
```
|             |   g000009_1 |    g000106   |   g000226    | ...          | 
| ----------- | ----------- |  ----------- |  ----------- | -----------  | 
| 5448_1_RPMI |  29.84732   |   26.85443   |   32.74102   | ...          | 
| 5448_1_SERA |  28.29469   |   28.59965   |   32.50166   | ...          | 
| 5448_2_RPMI |  29.24547   |   28.68960	 |   32.70440   | ...          | 
| 5448_2_SERA |  29.42946   |   27.61975   |   32.00997   | ...          | 
| ...         |  ...        |  ...         |  ...         | ...          | 


```
head(Strep_met_GC_data)
```
|             |  METAB_3.Aminoglutaric.acid |  METAB_HMDB0000034 |   METAB_HMDB0000050  | ...          | 
| ----------- | --------------------------- |  ----------------- |  ------------------- | -----------  | 
| 5448_1_RPMI |  -4.8902845                 |   5.2105774        |   3.423466           | ...          | 
| 5448_1_SERA |  1.7011588                  |   4.0033932	       |   4.016438           | ...          | 
| 5448_2_RPMI |  -6.9542200                 |   4.3488065	       |   2.682624	          | ...          | 
| 5448_2_SERA |  -3.1877999                 |   2.6740439        |   3.488940           | ...          | 
| ...         |  ...                        |  ...               |  ...                 | ...          | 


## 2. Network Analysis for Pathway Enrichment Outputs
### 2.1. Results from DIABLO :shipit:

```
loadings_data <- final.plsda.diablo$loadings 

loadings_data_trans <- loadings_data$Transcript   %>%  as.data.frame(.) %>% 
                 rownames_to_column(., var = "core_entity_id") %>% 
                 rename(loadings = comp1)

head(loadings_data_trans)
```
|core_entity_id |   loadings    |
| ------------- | ------------- |
| g000004       |-3.735591e-02	|
| g000007       |-6.332548e-03	|
| g000009_1     |3.790403e-03	  |
| ...           |  ...          |


### 2.2. Results from Differentially Expression (DE) Analysis :shipit:

```
head(all_DE.df.strep)
```
| strain |    entity_id    |  core_entity_id |   logFC      |   FDR        | gene_name    | core_consensus_annotation                      |  
| ------ | --------------- | --------------- |  ----------- | ------------ | ---------    | ---------------------------------------------- |
| 5448   |  EW021_RS00005	 |   g010876       |   -3.1449886 | 3.507962e-20 | dnaA         | chromosomal replication initiator protein DnaA |
| 5448   |  EW021_RS00010  |   g009819	     |   -3.0151278 | 9.302070e-26 | dnaN	        | DNA polymerase III subunit beta                |
| 5448   |  EW021_RS00015  |   g008944       |   -0.7921787 | 2.949191e-01 | EW021_RS00015| DUF951 family protein                          |
| 5448   |  EW021_RS00020	 |   g000450       |   -1.6880047 | 5.128776e-06 | ychF         | redox-regulated ATPase YchF                    |
| ...    |  ...            |  ...            |  ...         | ...          | ...          | ...                                            |

