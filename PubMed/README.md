
## PubMed 
PubMed® 包含来自 MEDLINE、生命科学期刊和在线书籍的 3800 多万条生物医学文献引文。引文可能包括来自 PubMed Central 和出版商网站的全文内容链接。

## 文献查找

```bash
# pumbed 搜索框直接输入PMID,注意冒号后面需要空格
PMID: 40078721
```

## BioPython 

- 获取总数
```bash
def get_all_pmids(keyword, max_total=20000):
    # 初始化搜索，获取总文献数和会话参数
    handle = Entrez.esearch(
        db="pubmed",
        term=keyword,
        retmax=0,          # 不返回具体ID，仅获取总数
        usehistory="y"      # 启用历史会话
    )
    result = Entrez.read(handle)
    total = int(result["Count"]) # 总数
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]
    print(f"总文献数: {total}")

# 示例：搜索全文中包含"CRISPR"的文献
get_all_pmids("CRISPR", max_total=20000)
```

- 获取PMID
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
