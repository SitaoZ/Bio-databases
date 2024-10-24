## KEGG


### HSA
生物通路是细胞内分子之间的一系列相互作用，可导致特定产物或细胞变化。通路涉及一系列发挥不同作用的基因，这些基因构成了“通路基因集”。
KEGG 通路是最常用的通路数据库。它通过 REST API（https://rest.kegg.jp/）提供数据。有几种命令可以检索特定类型的数据。要检索通路基因集，
我们可以使用“link”命令，如以下 URL 所示（“link”表示将基因链接到通路）。在 Web 浏览器中输入 URL 时：
```r
keggGeneSets = read.table(url("https://rest.kegg.jp/link/pathway/hsa"), sep = "\t")
head(keggGeneSets)
#         V1            V2
# 1 hsa:10327 path:hsa00010
# 2   hsa:124 path:hsa00010
# 3   hsa:125 path:hsa00010
# 4   hsa:126 path:hsa00010
# 5   hsa:127 path:hsa00010
# 6   hsa:128 path:hsa00010

# 在这个两列表格中，第一列包含 Entrez ID 类型的基因。我们删除“hsa:”前缀，同时删除第二列通路 ID 的“path:”前缀。
#      V1       V2
# 1 10327 hsa00010
# 2   124 hsa00010
# 3   125 hsa00010
# 4   126 hsa00010
# 5   127 hsa00010
# 6   128 hsa00010

keggGeneSets[, 1] = gsub("hsa:", "", keggGeneSets[, 1])
keggGeneSets[, 2] = gsub("path:", "", keggGeneSets[, 2])
head(keggGeneSets)

# 可以通过“list”命令获取完整的通路名称。
keggNames = read.table(url("https://rest.kegg.jp/list/pathway/hsa"), sep = "\t")
head(keggNames)

#         V1                                                     V2
# 1 hsa01100              Metabolic pathways - Homo sapiens (human)
# 2 hsa01200               Carbon metabolism - Homo sapiens (human)
# 3 hsa01210 2-Oxocarboxylic acid metabolism - Homo sapiens (human)
# 4 hsa01212           Fatty acid metabolism - Homo sapiens (human)
# 5 hsa01230     Biosynthesis of amino acids - Homo sapiens (human)
# 6 hsa01232           Nucleotide metabolism - Homo sapiens (human)

```

在这两个命令中，我们获取了人类的数据，其中对应的 KEGG 代码为“hsa”。其他生物的代码可以从 KEGG 网站（例如小鼠的代码为“mmu”）或通过 https://rest.kegg.jp/list/organism 找到。

请记住，KEGG 通路仅对学术用户免费。如果您将其用于商业目的，请联系 KEGG 团队获取许可。


除了直接从 URL 读取数据外，还有一些 R 包可以帮助从 KEGG 数据库获取数据，例如 KEGGREST 包或 clusterProfiler 包中的 download_KEGG() 函数。但本质上，它们都是使用 REST API 获取 KEGG 数据。

