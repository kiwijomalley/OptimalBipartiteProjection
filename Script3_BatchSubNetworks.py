import csv
import time
import sys

fp=open("../HRR-PHN-list/PHN-index.txt",'r')
lines=fp.readlines()
N=len(lines)
fp.close()
print N

def subPHN(year, start, end):
    year = str(year)
    path00 = "/projects/active/54054/idata/jomalley"
    os.chdir(path00)
    path0 = "./PHN-subnetwork/" + year # results
    if not os.path.exists(path0):
        os.makedirs(path0)

    path1 = "/projects/active/54054/idata/jomalley/PHNSplit/PHNphy/"

    for i in range(start, end):  #N
        print('processing %d out of %d items...'%(i+1, end), '\\r',)
        time.sleep(0.1)
        batch_id=int(i/160)+1
        
        # batch directory
        fo_dir = path0 + "/batch-" + str(batch_id)
        if not os.path.exists(fo_dir):
            os.mkdir(fo_dir)
            print("Making new directory" + fo_dir)
        else:\
            print("Directory: " + fo_dir + " already exists")

        #name ="./batch-"+str(batch_id)+"/PHN"+str(i)

        phyname = path1 + "phylist-PHN" + str(i) + "sound.txt"
         
        fpr = open(phyname,'r')
        phylist=\{\}
        while 1:
            line = fpr.readline()
            if not line:
                break
            words=line.split(",")
            phylist[words[0]]=1
        fpr.close()

        print(str(i)+"PHNlist "+str(len(phylist)))
        #print(phylist)

        andname = fo_dir + "/PHN" + str(i) + "-anda.txt" #Line to change if forming different networks
        #name+"-and.txt"
        fp1 = open(andname,'w')
        #orname=name+"-or.txt"
        #fp2 = open(orname,'w')

        phnedges=\{\}
        recordname = "./directed_network_2018a.txt"
        f=open(recordname,'r')
        while 1:
            line=f.readline()
            if not line:
                break
            words=line.split(";")
            nodes=words[0].split("_") #Splits dyad ID into NPIs
            if nodes[0] in phylist and nodes[1] in phylist:
                if words[0] not in phnedges:
                    phnedges[words[0]]=[0 for col in range(80)]
                for j in range(1,len(words)-1):
                    phnedges[words[0]][j-1]=round(float(words[j]),3)
            
            if nodes[0] in phylist or nodes[1] in phylist:
                output=nodes[0]+","+nodes[1]+","
                for j in range(1,len(words)-1):
                    output=output+words[j]+","
                output=output+"\\n"
                fp2.write(output)
        f.close()
        
        recordname = "./undirected_network_2018a.txt"
        
        f=open(recordname,'r')
        while 1:
            line=f.readline()
            if not line:
                break
            words=line.split(";")
            nodes=words[0].split("_")
            if nodes[0] in phylist and nodes[1] in phylist:
                if words[0] not in phnedges:
                    phnedges[words[0]]=[0 for col in range(80)]
                for j in range(1,21):
                    phnedges[words[0]][j+19]=round(float(words[j]),3)
        f.close()
        
        recordname = "./binary_dir_network_2018a.txt"
        f=open(recordname,'r')
        while 1:
            line=f.readline()
            if not line:
                break
            words=line.split(",")
            nodes=words[0].split("_")
            if nodes[0] in phylist and nodes[1] in phylist:
                if words[0] not in phnedges:
                    phnedges[words[0]]=[0 for col in range(80)]
                for j in range(1,21):
                    phnedges[words[0]][j+39]=round(float(words[j]),3)
		    print(words[j]) # 01/07/2020
        f.close()
        
        recordname = "./binary_undir_network_2018a.txt"
        
        f=open(recordname,'r')
        while 1:
            line=f.readline()
            if not line:
                break
            words=line.split(",") 
            nodes=words[0].split("_")
            if nodes[0] in phylist and nodes[1] in phylist:
                if words[0] not in phnedges:
                    phnedges[words[0]]=[0 for col in range(80)]
                for j in range(1,21):
                    phnedges[words[0]][j+59]=round(float(words[j]),3)
        f.close()
        
        
        for item in phnedges:
            res=str(item)+";"
            tmp=phnedges[item]
            for i in range(80):
                res=res+str(tmp[i])+";"
            fp1.write(res+"\\n")
        fp1.close()

if __name__ == "__main__":
    year = int(sys.argv[1])
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    subPHN(year, start, end)}