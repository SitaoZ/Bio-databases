
## PubMed 
PubMed® 包含来自 MEDLINE、生命科学期刊和在线书籍的 3800 多万条生物医学文献引文。引文可能包括来自 PubMed Central 和出版商网站的全文内容链接。

## 文献查找

```bash
# pumbed 搜索框直接输入PMID,注意冒号后面需要空格
PMID: 40078721
```

## BioPython 

```python
from Bio import Entrez
import time

Entrez.email = "your_email@example.com"  # 必须提供邮箱

def search_pubmed(keyword, max_results=10):
    """搜索包含关键词的文献，返回PMID列表"""
    handle = Entrez.esearch(
        db="pubmed",
        term=keyword,
        retmax=max_results,
        usehistory="y"
    )
    result = Entrez.read(handle)
    handle.close()
    return result["IdList"]

# 示例：搜索全文中包含"CRISPR"的文献
pmids = search_pubmed("CRISPR", max_results=50)
print("找到PMID列表:", pmids)


```
