# Introduce

Up to now, many studies have demonstrated that the analysis of disease similarity is a useful way to understand the causes of human diseases and promote new use of existing drugs. Besides, the quantification of disease similarity has been applied to the discovery of different types of disease-related biological entity associations, such as disease-related genes, disease-related miRNAs. Moreover, identification of similar diseases offers a better knowledge to improve the existing systems of disease classification, since the more similar two diseases are, the more likely they belong to the same category. Furthermore, the disease networks not only provides a rapid visual reference of links among diseases, but also offers a valuable global perspective for physicians, genetic counselors and biomedical researchers.  

In the past several decades, many computational approaches have been proposed to quantify the similarity between two diseases from a single and isolated view to a systematic level. To facilitate the computation and analysis of disease similarity, we developed a python package named **NetBasedDSim**.

Besides, theare are two public resources for disease similarity implemented by [our group](http://bioinformatics.csu.edu.cn/).
+ [DSNet](http://bioinformatics.csu.edu.cn/DSNet/)
+ [dSimer](http://www.bioconductor.org/packages/release/bioc/html/dSimer.html)

### Integrated methods

- Semantic-based methods:
	+ ResnikSim
	+ XuanSim

- Function-based methods:
	+ NetSim
	+ ModuleSim
	+ FunSim

- Text mining-based methods:
	+ MimMiner
	+ CosineDFV
	+ MicrbeSim

- Information fusion-based methods:
	+ IDN
	+ RADAR
	+ mpDisNet


-----------------------------------
## Installation
-----------------------------------
* software
 > python 3
* python package
 > gensim--3.6.0  
 > network--2.1  
 > numpy--1.17.4  
 > fastdtw--0.3.4  
 > python-igraph--0.8.0  




