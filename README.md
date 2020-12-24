# kegg数据整理

br08201下载地址：

```shell
https://www.genome.jp/kegg-bin/get_htext?br08201.keg
```



从kegg所有反应的list中提取出反应物。

```shell
cat br08201.keg |awk '$1 == "E"'|awk -F'<=>' '{print $1}'|tr -s ' '|awk '{ for(i=3; i<NF; i++) printf "%s",$i OFS; if(NF) printf "%s",$NF; printf ORS}'|sed 's/ + /*/g'|awk -F "*" '{print $1}'|sort|uniq|sed 's/[0-9][0-9] \|[0-9] \|^n \|^n-1 \|^[0-9]n //'|sort|uniq > reactants1.list
```



反应物最多有5个

最终得到的5个文件为：

```shell
reactants1.list
reactants2.list
reactants3.list
reactants4.list
reactants5.list
```



整合：

```shell
cat reactants*|sort|uniq > Allreactants.list
```



提取产物：

```shell
cat br08201.keg |awk '$1 == "E"'|awk -F'<=>' '{print $2}'|tr -s ' '|sed 's/^ //'|sed 's/ + /*/g'|awk -F "*" '{print $1}'|sed 's/^[0-9] \|^n \|^[0-9][0-9] \|^n-[0-9] \|[0-9]n //g'|sort|uniq > product1.list
```



产物文件一共7个：

```shell
product1.list
product2.list
product3.list
product4.list
product5.list
product6.list
product7.list
```



整合产物文件：

```shell
cat product*|sort|uniq > allproduct.list
```



整合反应物和产物：

```shell
cat allproduct.list allreactant.list |sort|uniq > allmetabolites.list
```



写一个爬虫爬取

所有kegg compounds的id对应的name

然后用allmetabolites.list map到爬下来的文件得到所有代谢物的cid，但实际上一些代谢物是属于glycan类的，这些是没有结构的。



爬虫内容：

```python
import urllib3
from bs4 import BeautifulSoup
from urllib.error import HTTPError,URLError
import threading
import sys
import io
import json
import socket
import logging

kegg_cids = []
kegg_names = []
non_exist_compounds = []
connection_error = []
cids_names_dic = {}
headers = {'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.89 Safari/537.36'}


with open('./kegg_cids_names/redownload.list','r') as f:
    for i in f:
        i = i.strip('\n')
        kegg_cids.append(i)

class CrawlerThread(threading.Thread):
    def __init__(self,kegg_cid,url,tid):
        threading.Thread.__init__(self)
        self.url = url
        self.tid = tid
        self.kegg_cid = kegg_cid

    def run(self):
        global cids_names_dic
        try:
            http = urllib3.PoolManager(timeout=20)
            resp = http.request('GET',self.url,headers=headers)
        except Exception as e:
            print('kegg compound {}:timeout!'.format(self.kegg_cid))
            logging.basicConfig(filename='./kegg_cids_names/logger.log', level=logging.INFO)
            logging.error('error:timeout failed compound cid:%s',self.kegg_cid)
        try:
            data = BeautifulSoup(resp.data,'html.parser')
            name = data.find(class_ = 'td21').div.div.text.strip('\n').replace("\n","")
            cids_names_dic[self.kegg_cid] = name
            print(self.kegg_cid + " : " + name)
        except Exception as e:
            print(e)
            print("no such compound:" + self.kegg_cid)
            non_exist_compounds.append(self.kegg_cid)
            pass


class Crawler():
    def __init__(self,name,thread_num):
        self.name = name
        self.thread_num = thread_num
        self.thread_pool = []

    def getAllName(self):
        url_index = 0
        while url_index < len(kegg_cids):
            thread_index = 0
            print("processing:")
            while thread_index < self.thread_num and url_index + thread_index < len(kegg_cids):
                kegg_cid = kegg_cids[url_index + thread_index]
                getNameResult = self.getName(kegg_cid,thread_index)
                thread_index += 1
            url_index = url_index + thread_index
            for thread in self.thread_pool:
                thread.join()
            with open('./kegg_cids_names/redownload3Compounds{}.list'.format(str(int((url_index + 1)/thread_num))),'w') as f:
                f.write(json.dumps(cids_names_dic))


    def getName(self,kegg_cid,tid):
        url = 'https://www.genome.jp/dbget-bin/www_bget?cpd:' + kegg_cid
        print(url)
        crawler_thread = CrawlerThread(kegg_cid,url,tid)
        self.thread_pool.append(crawler_thread)
        crawler_thread.start()


if __name__ == "__main__":
    thread_num = 200
    name = 'kegg_crawler'
    kegg_crawler = Crawler(name,thread_num)
    kegg_crawler.getAllName()
```



如何mapping?

```shell
cat allkegg_cids_names.list |sed 's/  /;/g' > allkegg_cids_names2.list

#这里涉及到一个常用的awk处理两个文件的方法FNR == NR,在处理第一个文件时只执行第一个brace的内容
#第二个文件执行了第二个花括号的内容，先是指定了OFS然后根据column是否在第一个文件产生的数组中来拼接字符串，然后输出到文件中。我们对所有column都进行了
#类似的操作。
awk -F';' 'FNR == NR{a[$1] ;next} {OFS=";" ;$0=$1 OFS ($2 in a ? $2 : "")}1' \ 
allmetabolites_sort.list allkegg_cids_names2.list > column2.list

gawk -F';' 'ARGIND == 1 {if($2) print $0;next} ARGIND ==2 {if($2) print $0 ;next} ARGIND == 3 {if($2) print $0 ;next} ARGIND == 4 {if ($2) print $0 ;next} ARGIND == 5 {if ($2) print $0 ;next} ARGIND == 6 {if ($2) print $0 ;next} ARGIND == 7 {if ($2) print $0 ;next} ARGIND == 8 {if ($2) print $0 ;next} ARGIND == 9 {if ($2) print $0 ;next} ARGIND == 10 {if ($2) print $0 ;next} ARGIND == 11 {if ($2) print $0 ;next} ARGIND == 12 {if ($2) print $0 ;next} ARGIND == 13 {if ($2) print $0 ;next} ARGIND == 14 {if ($2) print $0 ;next} ARGIND == 15 {if ($2) print $0 ;next} ARGIND == 16 {if ($2) print $0 ;next} ARGIND == 17 {if ($2) print $0 ;next} ARGIND == 18 {if ($2) print $0 ;next} ARGIND == 19 {if ($2) print $0 ;next} ARGIND == 20 {if ($2) print $0 ;next} ARGIND == 21 {if ($2) print $0 ;next} ARGIND == 22 {if ($2) print $0 ;next} ARGIND == 23 {if ($2) print $0 ;next} ARGIND == 24 {if ($2) print $0 ;next} ARGIND == 25 {if ($2) print $0 ;next}' column2.list column3.list column4.list column5.list column6.list column7.list column8.list column9.list column10.list column11.list column12.list column13.list column14.list column15.list column16.list column17.list column18.list column19.list column20.list column21.list column22.list column23.list column24.list column25.list column26.list column27.list|sort|uniq > extract_all_metabolites.list
```



提取相应的pubchem id

```shell
grep -Ff extract_all_metabolites2.list kegg_compound_pubchem.list |cut -f2|cut -d':' -f2 > extract_all_metabolites3.list
```



