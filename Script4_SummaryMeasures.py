import os
import sys
import numpy as np
import random
import networkx as nx
import math
import decimal
network_measure_dim=31 

# 16 kinds of triads in paper Analysis of referral network
def judge(p,q,s,g):
    res=["003","012","102","021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300"]
    dic=\{\}
    for i in range(0,16):
        dic[res[i]] = i+1
\
    #print "id"+str(p)+","+str(q)+","+str(s)
    mu=0
    asy=0
    nul=0
    extra=""
    if p in g[q] and q in g[p]:
        mu=mu+1
    elif p in g[q] or q in g[p]:
        asy=asy+1
    else:
        nul=nul+1

    if p in g[s] and s in g[p]:
        mu=mu+1
    elif p in g[s] or s in g[p]:
        asy=asy+1
    else:
        nul=nul+1

    if s in g[q] and q in g[s]:
        mu=mu+1
    elif s in g[q] or q in g[s]:
        asy=asy+1
    else:
        nul=nul+1
    key = str(mu)+str(asy)+str(nul)

    if(key=="021"):
        if ((q in g[p] and s in g[p]) or (p in g[q] and s in g[q]) or (q in g[s] and p in g[s])):
            extra="D"
        elif ((p in g[q] and p in g[s]) or (q in g[p] and q in g[s]) or (s in g[q] and s in g[p])):
            extra = "U"
        else:
            extra = "C"
    if(key=="111"):
        if ((q in g[p] or s in g[p]) and (p in g[q] or s in g[q]) and (q in g[s] or p in g[s])):
            extra="D"
        else:
            extra = "U"
    if(key=="030"):
        if ((q in g[p] or s in g[p]) and (p in g[q] or s in g[q]) and (q in g[s] or p in g[s])):
            extra="C"
        else:\
            extra="T"

    if(key=="120"):
        if ((p not in g[q] and p not in g[s]) or(q not in g[s] and q not in g[p]) or (s not in g[q] and s not in g[p])):
            extra="D"
        elif (q not in g[p] and s not in g[p]) or (p not in g[q] and s not in g[q]) or (p not in g[s] and q not in g[s]):
            extra="U"
        else:
            extra="C"

    key = key+extra
    #print key
    return dic[key]

def findmax_dia(g, s,curmax,n):
    local_visited = [-1] * n
    local_visited[s]=0
    queue =[s]
    node_count=0
    while (len(queue)!=0):
        v = queue.pop(0)
        for item in g[v]:
            if local_visited[item]==-1:
                local_visited[item]=local_visited[v]+1
                node_count+=1
                queue.append(item)
    #print max(visited)
    return max(local_visited)

def avg_distance(s,g,n):
    local_visited = [-1] * n
    local_visited[s]=0
    queue =[s]
    node_count=0
    while (len(queue)!=0):
        v = queue.pop(0)
        for item in g[v]:
            if local_visited[item]==-1:
                local_visited[item]=local_visited[v]+1
                node_count+=1
                #if node_count%10000==0:
                 #   print "now bfs nodes "+str(node_count)
                queue.append(item)
    #print max(visited)
    valid_num=0
    total_sum=0
    for i in range(0, len(local_visited)):
        if local_visited[i]>-1 and i!=s:
            valid_num=valid_num+1
            total_sum=total_sum+local_visited[i]
    if valid_num==0:
        return [-1,0]
    else:
        return [1.0/total_sum, max(local_visited)]

def localC(alpha,beta,loc,n):
    if loc < n*beta:
        res = loc*(1-alpha)/2/beta
    else:
        res = (loc-beta)*(1-alpha)/2/(n-beta)+0.5*(1+alpha)
    return res

def TransP(maxE,curE):
    delta = maxE-curE
    if delta <0:
        return 1.0
    elif maxE>0:
        return math.exp(-delta*100/maxE)
    else:
        return 0.0

# reciprocity for directed and binary network
def dir_bi_recip(t1, t2, numNodes, den):
    mutual = sum([t1[i] * t2[i] for i in range(len(t1))])
    return (2*mutual/numNodes/(numNodes-1) - den*den)
    
def compute_measures(index,edges,w_id):
    print str(w_id)+" "+str(len(index))+" "+str(len(edges))
    node_p_dict=\{\} #node position dictionary

    local_net_measure=[0 for col in range(network_measure_dim)] #Initialize all networks measures to 0
    local_net_measure[0]=len(index) # number of non-isolated nodes
    if len(index)<=20 or len(edges)<20: # number of non-isolated nodes less than 20 was suppressed, 0120
        return [local_net_measure,node_p_dict] 
    
    tmp_sum=0
    for item in edges:
        tmp_sum+=edges[item] #Total number of edges, numerator of strength
        if w_id<20 or (w_id>=40 and w_id<60):
        local_net_measure[1]=round(tmp_sum/local_net_measure[0]/(local_net_measure[0]-1), 3) 
    elif (w_id>=20 and w_id<40):
        #print(local_net_measure[0])
	#print(local_net_measure[0]-1)
        local_net_measure[1]=round(tmp_sum*2/local_net_measure[0]/(local_net_measure[0]-1), 3)
    else:
        undir_total_edges=\{\} #Computes a set with undirected edges included?
        for item in edges:
            rev_item=(item[1],item[0])
            if (item not in undir_total_edges) and (rev_item not in undir_total_edges):
                undir_total_edges[item]=1
        #print(len(undir_total_edges))
	#print(local_net_measure[0])
        local_net_measure[1]=round(len(undir_total_edges)*2.0/local_net_measure[0]/(local_net_measure[0]-1), 3)
	#print(local_net_measure[1])
	#print('len of index is', len(index))
	#print('len of edges is', len(edges))

    node_strength=\{\}
    outdeg=\{\}
    indeg=\{\}
    for item in edges:
        if item[0] not in node_strength:
            node_strength[item[0]]=0
        if item[0] not in outdeg:
            outdeg[item[0]]=0
        if item[1] not in indeg:
            indeg[item[1]]=0
        if item[1] not in node_strength:
            node_strength[item[1]]=0
        node_strength[item[0]]+=edges[item]
        node_strength[item[1]]+=edges[item]
        outdeg[item[0]]+=edges[item]
        indeg[item[1]]+=edges[item]

    local_net_measure[2]=round(np.var(node_strength.values()),3) #Centralization
    outdeg_list=[]
    indeg_list=[]
    for node in outdeg:
        outdeg_list.append(outdeg[node])
        if node in indeg:
            indeg_list.append(indeg[node])
        else:
            indeg_list.append(0)
    for node in indeg:
        if node not in outdeg:
            outdeg_list.append(0)
            indeg_list.append(indeg[node])
    tmp_res=np.corrcoef(outdeg_list,indeg_list)
    local_net_measure[3]=round(tmp_res[0][1],3) #Correlation between in-degree and out-degree

    tmp_list1=[]
    tmp_list2=[]
    for item in edges:
        if item[0] in indeg:
            tmp_list1.append(indeg[item[0]])
        else:
            tmp_list1.append(0)
        if item[1] in indeg:
            tmp_list2.append(indeg[item[1]])
        else:
            tmp_list2.append(0)
    tmp_res=np.corrcoef(tmp_list1,tmp_list2)
    local_net_measure[4]=round(tmp_res[0][1],3) #in-in Assortativity

    tmp_list1=[]
    tmp_list2=[]
    for item in edges:
        if item[0] in indeg:
            tmp_list1.append(indeg[item[0]])
        else:
            tmp_list1.append(0)
        if item[1] in outdeg:
            tmp_list2.append(outdeg[item[1]])
        else:
            tmp_list2.append(0)
    tmp_res=np.corrcoef(tmp_list1,tmp_list2)
    local_net_measure[5]=round(tmp_res[0][1],3) #in-out Assortativity

    tmp_list1=[]
    tmp_list2=[]
    for item in edges:
        if item[0] in outdeg:
            tmp_list1.append(outdeg[item[0]])
        else:
            tmp_list1.append(0)
        if item[1] in indeg:
            tmp_list2.append(indeg[item[1]])
        else:
            tmp_list2.append(0)
    tmp_res=np.corrcoef(tmp_list1,tmp_list2)
    local_net_measure[6]=round(tmp_res[0][1],3) #out-in Assortativity

    tmp_list1=[]
    tmp_list2=[]
    for item in edges:
        if item[0] in outdeg:
            tmp_list1.append(outdeg[item[0]])
        else:
            tmp_list1.append(0)
        if item[1] in outdeg:
            tmp_list2.append(outdeg[item[1]])
        else:
            tmp_list2.append(0)
    tmp_res=np.corrcoef(tmp_list1,tmp_list2)\
    local_net_measure[7]=round(tmp_res[0][1],3) #out-out Assortativity
    
    # reciprocity
    tmp_list1=[]
    tmp_list2=[]
    for item in edges:
        rev_key=(item[1], item[0])
        tmp_list1.append(edges[item])
        if rev_key in edges:
            tmp_list2.append(edges[rev_key])
        else:
            tmp_list2.append(0)
    if (w_id >= 40 and w_id <= 59): # edited on 12/09/2019
        tmp_res = dir_bi_recip(tmp_list1, tmp_list2, local_net_measure[0], local_net_measure[1]) #Calls function above
	local_net_measure[8] = round(tmp_res, 3)
	print("reciprocity: " + str(local_net_measure[8]))
    else: #Might want to alter so that only compute for directed network
        tmp_res=np.corrcoef(tmp_list1,tmp_list2)
        local_net_measure[8]=round(tmp_res[0][1],3)
	print("reciprocity: " + str(local_net_measure[8]))
    #for binary existence, reciprocity may be nan


    n=len(index) #index = nodes
    g = \{x:[] for x in xrange(n)\}
    for item in edges:
        x=item[0]
        y=item[1]
        g[x].append(y) #. Adds to queue or array
        g[y].append(x)
    for i in range(0,n):
        g[i]=list(set(g[i]))

    fenmu=0
    fenzi=0
    csum=0
    for i in range(0,n):
        deg = len(g[i])
        localfenzi=0
        for j in range(0,deg-1):
            for k in range(j+1,deg):
                if g[i][k] in g[g[i][j]]:
                    fenzi=fenzi+1
                    localfenzi=localfenzi+1
        fenmu=fenmu+deg*(deg-1)/2
        if deg>1:
            csum=csum+localfenzi*2.0/(deg-1)/deg

            if i not in node_p_dict:
                node_p_dict[i]=[-99 for col in range(6)]
            node_p_dict[i][1]=round(localfenzi*2.0/(deg-1)/deg,3)

    if fenmu>0:
        local_net_measure[9]=round(fenzi*1.0/fenmu,3) #Transitivity (undirected)
        local_net_measure[10]=round(csum*1.0/n,3) #Clustering

    if w_id <=19 or (w_id>=40 and w_id<=59):
        #triads
        g=\{\}
        for i in range(0,n):
            g[i]=[]
        for item in edges:
            x=item[0]
            y=item[1]
            g[x].append(y)
            
        kind = [0 for col in range(17)]
        for j in range(0,100000):
            #if(j%1000==0):
            #    print j
            p = random.randint(0,n-1)
            while(len(g[p])==0):
                p = random.randint(0,n-1)
            #print g[p]
            q = random.randint(0,len(g[p])-1)
            #while(q==p):
            #print q
            q = g[p][q] #random.randint(0,n-1)
            s = random.randint(0,n-1)
            while(s==p or s==q):
                s =random.randint(0,n-1)

            tmp_index=int(judge(p,q,s,g))
            kind[tmp_index]=kind[tmp_index]+1

        #11-25
        for i in range(2,17):
            local_net_measure[11+i-2]=kind[i] #Relative frequency of each type of triad


    g = \{x:[] for x in xrange(n)\}
    for item in edges:
        x=item[0]
        y=item[1]
        g[x].append(y)
        g[y].append(x)
    for i in range(0,n):
        g[i]=list(set(g[i]))
    visited = [0] * n
    ret = 0
    for i in xrange(n):
        queue = [i]
        ret += 1 if i in g else 0
            #print queue
        for j in queue:
            if j in g:
                queue += g[j]
                visited[j]=ret
                for item in g[j]:
                    visited[item]=ret
                del g[j]

    local_net_measure[26]=ret #Number of components

    g = \{x:[] for x in xrange(n)\}
    for item in edges:
        x=item[0]
        y=item[1]
        g[x].append(y)
        g[y].append(x)
    max_dia=-1

    for j in range(0,1000):
        i = random.randint(0,n-1)
        temp = findmax_dia(g, i,max_dia,n)
        if temp > max_dia:
            max_dia = temp
    local_net_measure[27]=max_dia #Diameter

    local_net_measure[28]=local_net_measure[18]+2*local_net_measure[21]+2*local_net_measure[22]+local_net_measure[23]+4*local_net_measure[24]+6*local_net_measure[25] #Transitive triads (directed)

    local_net_measure[29]=local_net_measure[19]+local_net_measure[23]+local_net_measure[24]+2*local_net_measure[25] #Three-cycles

    #Node-level measures
    for item in node_strength:
        if item not in node_p_dict:
            node_p_dict[item]=[0 for col in range(6)]
        node_p_dict[item][0]=node_strength[item] #Node strength

    g = \{x:[] for x in xrange(n)\}
    for item in edges:
        x=item[0]
        y=item[1]
        g[x].append(y)
        g[y].append(x)
    for i in range(0,n):
        g[i]=list(set(g[i]))

    for i in range(0,n): #n
        node_p_dict[i][2]=round(avg_distance(i,g,n)[0],6) #?
        node_p_dict[i][5]=round(avg_distance(i,g,n)[1],3) #?

    #find the main component for eigen vector
    distri=\{\}
    for item in visited:
        if item not in distri:
            distri[item]=1
        else:
            distri[item]=distri[item]+1

    max_freq=-1
    max_key =-1
    for item in distri:
        if distri[item]>max_freq:
            max_freq=distri[item]
            max_key=item

    if max_key!=-1:
        #print max_key
        #print str(max_freq)+" nodes in main component"
        count=0
        component_list=\{\}
        for kk in range(n):
            if visited[kk]==max_key:
                component_list[kk]=1
        del visited

        main_n=0
        main_index=\{\}
        inverted_index=\{\}
        main_edges=[]
        for item in edges:
            x=item[0]
            y=item[1]
            if x in component_list and y in component_list:
                if x not in main_index:
                    main_index[x]=main_n
                    inverted_index[main_n]=x
                    main_n=main_n+1
                if y not in main_index:
                    main_index[y]=main_n
                    inverted_index[main_n]=y
                    main_n=main_n+1
                main_edges.append([main_index[x], main_index[y]])

        G=nx.Graph()
        for x,y in main_edges:
            G.add_edge(x,y)
        try:
            tmp_eigen=nx.eigenvector_centrality(G)
            for tmp_id in tmp_eigen:
                raw_id=inverted_index[tmp_id]
                node_p_dict[raw_id][3]=round(tmp_eigen[tmp_id],3) #Eigenvector centrality
        except:
            print "eigenvector error"

    
    n=len(index)
    g = \{x:[] for x in xrange(n)\}
    for item in edges:
        x=item[0]
        y=item[1]
        g[x].append(y)
        g[y].append(x)
    for i in range(0,n):
        g[i]=list(set(g[i]))

    CScore = [0 for col in range(n)]
    unit =0.05
    ubound = int(1/unit)
    looptime = 2
    for jj in range(0,ubound):
        for kk in range(0,ubound):
            alpha = jj*unit+unit
            beta = kk*unit+unit
            maxE = -1

            Cl = [0 for col in range(n)]

            Seq0 =[0 for col in range(n)]
            Seq=[0 for col in range(n)]
            for i in range(0,n):
                Seq0[i]=i              #Seq[i] is the index of node i
            #random.shuffle(Seq)
            i=0
            while (i<looptime):       #simulated annealing
                j=0
                for p in range(0,n):
                    Seq[p]=Seq0[p]
                while(j<looptime):        #each time swap 100 pairs location
                    ida = random.randint(0,n-1)
                    idb = random.randint(0,n-1)
                    tmp = Seq[ida]
                    Seq[ida] = Seq[idb]
                    Seq[idb] = tmp
                    j=j+1

                for p in range(0,n):
                    Cl[p] = localC(alpha,beta,Seq[p],n)
                curE = 0  # core quality in (10)
                for p in range(0,n-1):
                    for q in range(p+1,n):
                        if p in g[q]:
                            curE = curE+ Cl[p]*Cl[q]
                prob = TransP(maxE,curE)
                if random.random() < prob:
                    maxE = curE
                    for p in range(0,n):
                        Seq0[p] = Seq[p]
                    #change, update
                #print str(i)+","+str(maxE)
                #print Seq0
                i=i+1

            #compute score
            maxlc=-1
            maxlcid=0
            for p in range(0,n):
                tmp=localC(alpha,beta,Seq0[p],n)
                CScore[p]= CScore[p]+maxE*tmp
                if tmp>maxlc:
                    maxlcid=p
                    maxlc=tmp
    maxcs=-1
    for p in range(0,n):
        if CScore[p]>maxcs:
            maxcs = CScore[p]
    for p in range(0,n):
        node_p_dict[p][4]=round(CScore[p]*1.0/maxcs,3) #Local clustering coefficient?
    
    return [local_net_measure,node_p_dict]

if __name__ == "__main__":
    batchid=int(sys.argv[1])\
    #phnstart = int(sys.argv[2]) # edited on 02/10/2020
    #phnend = int(sys.argv[3])
    batch_size=160
    os.chdir("/projects/active/54054/idata/jomalley/PHN-subnetwork/")
    for PHN_id in range((batchid-1)*batch_size, batchid*batch_size): 
        network_measures=[]
        node_position=\{\}

        filename="./2018/batch-"+str(batchid)+"/PHN"+str(PHN_id)+"-anda.txt" # edited on 01/10/2020
        fpw_prefix="./2018/output-"+str(batchid)+"/"

        print filename
	if not os.path.exists(fpw_prefix):
	    os.makedirs(fpw_prefix)
        #print fpw_prefix
        #continue
        fpw_name=fpw_prefix+"PHN-"+str(PHN_id)+"-network-measures_a.txt"
        #fpw_name="./newoutput-"+str(batchid)+"/"+"PHN-"+str(PHN_id)+"-network-measures.txt"
        if os.path.isfile(fpw_name):
            continue
        
        fp=open(filename,'r')
        lines=fp.readlines()
	par_non_iso = [0 for i in range(80)] 
        for weight_id in range(0,80):
            index=\{\}
            edges=\{\}
            n=0
            count=0
	    allNodes = []
            for line in lines:
                words=line.split(";")
                pair=words[0].split("_")
		if pair[0] not in allNodes: 
		    allNodes.append(pair[0]) 
		if pair[1] not in allNodes: 
		    allNodes.append(pair[1]) 
                if float(words[weight_id+1])!=0: # note on 01/20, already teased out isolated nodes
		    if pair[0] == pair[1]: 
		        continue 
                    if pair[0] not in index:
                        index[pair[0]]=n
                        n=n+1
                    if pair[1] not in index:
                        index[pair[1]]=n
                        n=n+1
                    edges[(index[pair[0]], index[pair[1]])]=round(float(words[weight_id+1]),3)
		    #print(edges) 
                    count=count+1
                    #if count%10000==0:
                    #    print count
            #local_net_measure=[0 for col in range(network_measure_dim)]
            [a,tmp_nodes]=compute_measures(index,edges,weight_id) #Key function as generates results for all measures
	    a[30]=len(allNodes) - len(index) #number of isolated nodes
	    print('len of non iso nodes', len(index))
	    print('len of all nodes', len(allNodes))
            if weight_id>=20 and weight_id<=39: #undir and cts network
                for tmpi in range(3,9):
                    a[tmpi]=-99
                for tmpi in range(11,26):
                    a[tmpi]=-99
                a[28]=-99
                a[29]=-99
            #if weight_id>=40 and weight_id<=59: # edited on 12/09
            #    a[8]=-99 # recalculated for dir and bi network
            if weight_id>=60 and weight_id<=79: #undir and bi network
                for tmpi in range(3,9):
                    a[tmpi]=-99
                for tmpi in range(11,26):
                    a[tmpi]=-99
                a[28]=-99
                a[29]=-99

            network_measures.append(a)
                        
            node_inverse_dict=\{\}
            for item in index:
                node_inverse_dict[index[item]]=item

            #print node_inverse_dict

            for item in tmp_nodes:
                tmp_key=node_inverse_dict[item]

                if tmp_key not in node_position:
                    node_position[tmp_key]=[[0 for col in range(6)] for row in range(80)]
                node_position[tmp_key][weight_id]=tmp_nodes[item]
            
	# generate a list of indices to put the num of non iso and iso next to each other
	listOfIndex = [i for i in range(0,30)] # 0-29
	listOfIndex.insert(1,30)
	
        fpw_name=fpw_prefix+"PHN-"+str(PHN_id)+"-network-measures_a.txt"
        #fpw_name="./newoutput-"+str(batchid)+"/PHN-"+str(PHN_id)+"-network-measures.txt"
        fpw=open(fpw_name,'w')
        for tmpi in range(len(network_measures)):
            output=""
            #for tmpj in range(network_measure_dim):
            for tmpj in listOfIndex:
                output=output+str(network_measures[tmpi][tmpj])+","
            output=output+"\\n"
            fpw.write(output)
        fpw.close()

        fpw_name=fpw_prefix+"PHN-"+str(PHN_id)+"-node-measures_a.txt"
        
        #fpw_name="./newoutput-"+str(batchid)+"/PHN-"+str(PHN_id)+"-node-measures.txt"
        fpw=open(fpw_name,'w')
        for item in node_position:
            for tmpi in range(len(node_position[item])):
                output=str(item)+"\\t"+str(tmpi)+"\\t"
                for tmpj in range(6):
                    output=output+str(node_position[item][tmpi][tmpj])+"\\t"
                output=output+"\\n"
                fpw.write(output)
        fpw.close()
        
        #save results

        fp.close()}