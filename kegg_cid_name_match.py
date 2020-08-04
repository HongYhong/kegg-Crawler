import urllib3
from bs4 import BeautifulSoup
from urllib.error import HTTPError,URLError
import threading
import sys
import io
import json

kegg_cids = []
kegg_names = []
non_exist_compounds = []
connection_error = []
cids_names_dic = {}
headers = {'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.89 Safari/537.36'}


with open('./kegg_cid.list','r') as f:
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
            data = BeautifulSoup(resp.data,'html.parser')
            try:
                name = data.find(class_ = 'td21').div.div.text.strip('\n').replace("\n","")
                cids_names_dic[self.kegg_cid] = name
                print(self.kegg_cid + " : " + name)
            except Exception as e:
                print(e)
                print("no such compound:" + self.kegg_cid)
                non_exist_compounds.append(self.kegg_cid)
                pass

        except HTTPError as e:
            connection_error.append(self.kegg_cid)
            print('The server couldn\'t fulfill the request.')
            print('Error code: ', e.code)   
            logging.basicConfig(filename='./kegg_cids_names/logger.log', level=logging.INFO)
            logging.error('error:timeout failed compound cid:%s',self.kegg_cid)
        except URLError as e:
            connection_error.append(self.kegg_cid)
            print('We failed to reach a server.')
            print('Reason: ', e.reason)

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
            with open('./kegg_cids_names/cids_names_startAt{}.list'.format(url_index),'w') as f:
                f.write(json.dumps(cids_names_dic))
            for thread in self.thread_pool:
                thread.join()
                

    def getName(self,kegg_cid,tid):
        url = 'https://www.genome.jp/dbget-bin/www_bget?cpd:' + kegg_cid
        print(url)
        crawler_thread = CrawlerThread(kegg_cid,url,tid)
        self.thread_pool.append(crawler_thread)
        crawler_thread.start()


if __name__ == "__main__":
    thread_num = 100
    name = 'kegg_crawler'
    kegg_crawler = Crawler(name,thread_num)
    kegg_crawler.getAllName()
    with open('./non_exist_compounds.list','w') as f:
        f.write(non_exist_compounds)
    with open('./connection_error','w') as f:
        f.write(connection_error)

