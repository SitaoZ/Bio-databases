## VFDB
[Virulence factor database](https://www.mgc.ac.cn/VFs/main.htm) 整理了细菌病原体毒力因子的信息

- `致病性 Pathogenicity`  
  能够导致疾病发生的细菌称为细菌病原物(`bacterial pathogen`)。它的治病能力称为致病性(`pathogenicity`)。
- `毒力 Virulence`  
  毒力提供了致病性或引起疾病的可能性的定量测量
- `毒力因子 Virulence factor`  
  毒力因子是指使微生物能够在特定物种的宿主上或宿主内定居并增强其致病潜力的特性(即基因产物)。
  毒力因子包括细菌毒素、介导细菌附着的细胞表面蛋白、保护细菌的细胞表面碳水化合物和蛋白质以及可能导致细菌致病性的水解酶。


## 数据库下载
`核心数据集core dataset`仅包括与实验验证的 VF 相关的代表性基因，而`完整数据集(full dataset)`涵盖全部已知和预测的 VF 相关的所有基因。

```bash
# full dataset 
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz 
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
zcat VFDB_setB_pro.fas.gz | grep ">" | wc -l
# 27991


# core dataset
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz
zcat VFDB_setA_pro.fas.gz | grep ">" | wc -l
# 4261

# stat
wget https://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz
wget https://www.mgc.ac.cn/VFs/Down/VFs.xls.gz

```

## 建立索引
diamond比对需要
```bash
diamond makedb --in VFDB_setA_pro.fas --db VFDB_setA_pro

```
