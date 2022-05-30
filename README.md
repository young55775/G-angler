
# G-angler 0.0.1
A novel way to detect target gene in large scale screening




```pip install Gangler```
This package can help find candidate genes in high throughput mutagenesis and suppressor screening experiments
without mapping.
Please call variants via freebayes and annotate with snpeff before using this package.
Any advise is welcomed, please contact
 ### [e-mail]:guozhengyang980525@yahoo.co.jp


```gangler.prepare```

```a = Gangler.prepare.pool(filepath)```

   multiple renamed vcf files must be included in the filepath.
   examples:
```
filepath
|---juz113.vcf
|---juz114.vcf
|---juz115.vcf
|---juz116.vcf
|---juz117.vcf
|---juz118.vcg
```

```a = Gangler.prepare.getpool(filepath,filename)```

   multiple renamed folders must be included in the filepath, and a vcf file with filename must be included
   in the folders, the name of the folder must be splited with '_' to divide the folder name into strain name and
   WGS order name:
examples:
```
filepath
|---juz113_20221011jxskaosdosh---filename.vcf
|---juz114_20221011jxskaosdosh---filename.vcf
|---juz115_20221011jxskaosdosh---filename.vcf
|---juz116_20221011jxskaosdosh---filename.vcf
|---juz117_20221011jxskaosdosh---filename.vcf
|---juz118_20221011jxskaosdosh---filename.vcf
```


a temp folder will be automatically created in the filepath including renamed vcf files:
  ```
filepath
|----temp
|     |---juz113.vcf
|     |---juz114.vcf
|     |---juz115.vcf
|     |---juz116.vcf
|     |---juz117.vcf
|     |---juz118.vcg
|---juz113_20221011jxskaosdosh---filename.vcf
|---juz114_20221011jxskaosdosh---filename.vcf
|---juz115_20221011jxskaosdosh---filename.vcf
|---juz116_20221011jxskaosdosh---filename.vcf
|---juz117_20221011jxskaosdosh---filename.vcf
|---juz118_20221011jxskaosdosh---filename.vcf
```



```a = Gangler.prepare.pool()```
         ```a = Gangler.prepare.getpool()```

   a.taglist : a list of Dataframes which included columns: 'gene', 'ID', 'type', 'base', 'protein'，'tag'
   column 'tag' will be filled with strain name

example:
   a.taglist[1]:

|   gene |    ID      |   type   | base | protein |     tag |
| ---- | -------- |  ------ |  -------  |  ----  | -----  |
| ttn-1 | WBGenexxxx  | missense | C<G | Asp666Asn | juz113 |
| cla-1 | WBGenexxxx  | missense | C<G | Asp223Asn | juz114 |


```a = Gangler.pool.snpool(poollist,targetlist)```

a.result will contain all the result you need, small m_value indicate that there is high possibility that this gene is the target gene in this screening. Details will be explained in bioRxiv paper. 

# Examples

```

import Gangler as gl
a = gl.prepare.pool(r"C:\Users\YOUNG\Desktop\geneA")
b = gl.prepare.pool(r"C:\Users\YOUNG\Desktop\geneB")
c = gl.prepare.pool(r"C:\Users\YOUNG\Desktop\geneC")
d = gl.prepare.pool(r"C:\Users\YOUNG\Desktop\geneD")
e = gl.prepare.pool(r"C:\Users\YOUNG\Desktop\geneE")
f = gl.prepare.pool(r"C:\Users\YOUNG\Desktop\geneF")
j = gl.pool.snpool([a,b,c,d,e,f],['geneA','geneB',geneC','geneD','geneE','geneF'])
j.result.to_csv(r"C:\Users\YOUNG\Desktop\temp.csv")

```
