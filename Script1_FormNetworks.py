#First generates referral paths and then calculates the edge weights by considering intermediate or return 
#    visits or not, to retain as directed or reduce to undirected, and to leave as weighted or reduce to binary.
#The output files from this script include the referral paths and the network measures for the
# 20 ways of forming directed and the undirected networks (dyad ID in column 1, measures 1-20 in cols 2-21)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

@author: created by Chuankai An, adapted by Xin Ran and by James O'Malley

"""
import math
#from sets import Set
from datetime import datetime, date, time

# generate referral paths
def referral_paths(visit_info, curbene):
    global fpw # need to comment out for unit testing
    # handle the first line in raw data
    # now handle this in scan_visit_info()
    #if len(visit_info)==0 or "SFROMDT" in visit_info:
    #    return
    #visit_info=date;NPI;RVU>
    prephy="-1"
    path_list=[]
    node_list=[]
    visit_num=0
    rvu_sum=0

    previsit_date=sorted(visit_info)[0]
    for key in sorted(visit_info):
        single_day_visit = visit_info[key]

        for item in single_day_visit:
            curphy = item.split(";")[0]
            curRVU = str(item.split(";")[1])
            if curRVU=="None":
                curRVU=1.0
            curvisit_date = key
            pre_date_format = datetime.strptime(previsit_date,"%Y-%m-%d") 
            cur_date_format = datetime.strptime(curvisit_date,"%Y-%m-%d") 
            tmp_gap =  str(cur_date_format-pre_date_format)\
            if tmp_gap=="0:00:00":
                date_gap = 0
            else:
                date_gap = int(tmp_gap.split(" ")[0])
            # when gap>30, cut the node list to form a new referral 
            # starting from curphy
            if date_gap>365: #30, 365 (the latter allows all patient encounters within the year to build network)
                # wrap up the current referral path
                # by adding the number of visits and RVU 
                # of the last visit on the path
                node_list[len(node_list)-1]=node_list[len(node_list)-1]+str(visit_num)+";"+str(rvu_sum)+";" #Inserts number of visits into path record
                # store the current complete path by adding the arrow
                if len(node_list)>1:
                    tmp_path =""
                    for node_unit in node_list:
                        tmp_path=tmp_path+str(node_unit)+">"
                    path_list.append(tmp_path)
                    # where fpw fits in 
                    fpw.write(str(curbene)+":"+tmp_path+"\\n")
                    print("saving to referral paths")
                # empty it for next path, and start the next new path
                node_list = []
                one_node = str(key)+";"+str(curphy)+";"
                node_list.append(one_node) # same-day visits are already taken account for as referrals

                visit_num=1
                rvu_sum=float(curRVU)


            else:
                if curphy!=prephy:# and prephy !="-1":
                    if len(node_list)>0:
                        node_list[len(node_list)-1]=node_list[len(node_list)-1]+str(visit_num)+";"+str(rvu_sum)+";"
                    visit_num=1

                    one_node = str(key)+";"+str(curphy)+";"
                    node_list.append(one_node)
                    rvu_sum=float(curRVU)

                else:
                    visit_num=visit_num+1
                    rvu_sum = rvu_sum+float(curRVU)

            prephy = curphy
            previsit_date=curvisit_date

    # deal with the last phy visit entry, add its info
    node_list[len(node_list)-1]=node_list[len(node_list)-1]+str(visit_num)+";"+str(rvu_sum)+";"
    # store the last path into path_list
    if len(node_list)>1:
        tmp_path =""
        for node_unit in node_list:
            tmp_path=tmp_path+str(node_unit)+">"
        fpw.write(str(curbene)+":"+tmp_path+"\\n") #Writes: BeneID:date,NPI,visits,rvus>...
        path_list.append(tmp_path)
    return path_list
    
# calculate different edge weights
def edge_weights(path_list):
    global pairs
    for unit in path_list:
        nodes = unit.split(">") #Splits path into segments separated by using >
        phy_count=\{\}
        for i in range(0, len(nodes)-1):
            source=nodes[i].split(";") #source=beneid:date (or just date after first entry);npi;sum_visits;sum_rvu
            if source[1] not in phy_count:
                phy_count[source[1]]=0
            phy_count[source[1]]+=float(source[2]) #sum number of visits to npi across referral path
        ABk_no_others=\{\} # note!
        ABk_no_others = no_inter_no_return(nodes, phy_count, ABk_no_others) #Sets up how many visits with physician
        ABk_inter=\{\}
        ABk_inter = inter_no_return(nodes, phy_count, ABk_inter)
        return_no_inter(nodes, ABk_no_others)
        phy_on_path_list=list(phy_count.keys()) # edited
        inter_return(nodes, phy_on_path_list, ABk_inter)
    return

#nodes: beneid:date (or just date after first entry);npi;sum_visits;sum_rvu for each segment of referral path
#phy_count: number of visits to given physician across entire referral path

# no intermediate, no return
def no_inter_no_return(nodes, phy_count, ABk_no_others):
    global pairs\
    for i in range(0,len(nodes)-2): #-2 as look two positions ahead
        source=nodes[i].split(";") #Current physician
        target=nodes[i+1].split(";") #Next physician
        key = source[1]+"_"+target[1] #Dyad identifier
        if key not in pairs:
            # initiate 20 weights for this edge
            pairs[key]=[0 for col in range(20)]
        if key not in ABk_no_others:
            ABk_no_others[key]=0
            pairs[key][0]+=1 # existence; at least one visit to directed dyad
        # total count
        pairs[key][1]=pairs[key][1]+1 #+ instance of visit by source followed by visit to target
        # total ordering
        pairs[key][2]=pairs[key][2]+float(source[2])*float(target[2]) #+ number consec visits to source * consec visits to target (phys may repeat)
        ABk_no_others[key]+=float(source[2])*float(target[2]) #Prep for revisit scenario: return_no_inter
        #pairs[key][3]=pairs[key][3]+math.sqrt(float(source[3])*float(target[3]))
        # norm1
        pairs[key][3]=pairs[key][3]+1.0/(phy_count[source[1]]+phy_count[target[1]]) #+ [1] divided by summed total visits to source and target
        # norm2\
        pairs[key][4]=pairs[key][4]+float(source[2])*float(target[2])/phy_count[source[1]]/phy_count[target[1]]
        #+ [2] divided by maximum possible total count for given patient's referral path
    
    return ABk_no_others

# with intermediate, no return\
def inter_no_return(nodes, phy_count, ABk_inter):
    global pairs
    for i in range(0,len(nodes)-2):
        for j in range(i+1,len(nodes)-1):
            source=nodes[i].split(";")
            target=nodes[j].split(";")
            key = source[1]+"_"+target[1]

            if key not in pairs:
                pairs[key]=[0 for col in range(20)]
            if key not in ABk_inter:
                ABk_inter[key]=0
                pairs[key][5]+=1 # existence
            # total count and total ordering  
            pairs[key][6]=pairs[key][6]+1
            pairs[key][7]=pairs[key][7]+float(source[2])*float(target[2])
            ABk_inter[key]+=float(source[2])*float(target[2]) #Prep for revisit scenario: inter_return
            #pairs[key][3]=pairs[key][3]+math.sqrt(float(source[3])*float(target[3]))
            # norm1 and norm2
            pairs[key][8]=pairs[key][8]+1.0/(phy_count[source[1]]+phy_count[target[1]])
            pairs[key][9]=pairs[key][9]+float(source[2])*float(target[2])/phy_count[source[1]]/phy_count[target[1]]
    return ABk_inter

# with return, no intermediate
def return_no_inter(nodes, ABk_no_others):
    global pairs
    ABAk=\{\}
    ABA_visit=\{\}
    for i in range(0,len(nodes)-3): #-3 as look two positions ahead
        source=nodes[i].split(";")
        target=nodes[i+1].split(";")
        back=nodes[i+2].split(";") #Looks ahead 2 places to see who next visit is to 
        key = source[1]+"_"+target[1]
        if source[1]==back[1]: #Return to source physician
            if key not in pairs:
                pairs[key]=[0 for col in range(20)]
            if key not in ABAk:
                ABAk[key]=0
                ABA_visit[key]=0
                pairs[key][10]+=1 # existence
            pairs[key][11]+=1 #Total count
            pairs[key][12]+=float(source[2])*float(target[2])*float(back[2]) #Total ordering
            ABAk[key]+=1
            ABA_visit[key]+=float(source[2])*float(target[2])*float(back[2]) 

    for key in ABAk:
        phy_pair=key.split("_")
        reverse_key=phy_pair[1]+"_"+phy_pair[0]
        reverse_value=0
        if reverse_key in ABAk:
            reverse_value=ABAk[reverse_key] #Appears to just split weight between the two directions?
        pairs[key][13]+=ABAk[key]*1.0/(ABAk[key]+reverse_value) #Scale factor differs from non-revisit as condition on ABAk[key] and reverse
        reverse_value=0
        if reverse_key in ABk_no_others:
            reverse_value=ABk_no_others[reverse_key]
        pairs[key][14]+=ABA_visit[key]*1.0/(ABk_no_others[key]*reverse_value) #Scale factor differs from non-revisit as condition on ABk[key] and reverse
    return

# with both intermediate and return
def inter_return(nodes, phy_on_path_list, ABk_inter):
    global pairs
    ABAk_inter=\{\}
    ABAk_visit_inter=\{\}
    for i in range(len(phy_on_path_list)):
        for j in range(len(phy_on_path_list)):
            if j!=i:
                source_phy=phy_on_path_list[i]
                target_phy=phy_on_path_list[j]
                key=source_phy+"_"+target_phy

                sub_node_list=[]
                for k in range(len(nodes)-1):
                    cur_phy=nodes[k].split(";")
                    if cur_phy[1]==source_phy or cur_phy[1]==target_phy:
                        sub_node_list.append(nodes[k])
                #ABBK, ABBBABA
                node_merge_list=[]

                cur_phy=sub_node_list[0].split(";")[1]
                cur_visit=int(sub_node_list[0].split(";")[2])
                for k in range(1,len(sub_node_list)):
                    tmp=sub_node_list[k].split(";")
                    if tmp[1]==cur_phy:
                        cur_visit+=int(tmp[2])
                    else:
                        node_merge_list.append(cur_phy+";"+str(cur_visit))
                        cur_phy=tmp[1]
                        cur_visit=int(tmp[2])

                node_merge_list.append(cur_phy+";"+str(cur_visit))

                for k in range(0, len(node_merge_list)-2):
                    first_phy=node_merge_list[k].split(";")
                    second_phy=node_merge_list[k+1].split(";")
                    third_phy=node_merge_list[k+2].split(";")
                    if first_phy[0]==source_phy and second_phy[0]==target_phy and third_phy[0]==source_phy:
                        if key not in pairs:
                            pairs[key]=[0 for col in range(20)]
                        if key not in ABAk_inter:
                            pairs[key][15]+=1
                            ABAk_inter[key]=0
                            ABAk_visit_inter[key]=0
                        pairs[key][16]+=1
                        ABAk_inter[key]+=1
                        pairs[key][17]+=float(first_phy[1])*float(second_phy[1])*float(third_phy[1])
                        ABAk_visit_inter[key]+=float(first_phy[1])*float(second_phy[1])*float(third_phy[1])

    for key in ABAk_inter:
        phy_pair=key.split("_")
        reverse_key=phy_pair[1]+"_"+phy_pair[0]
        reverse_value=0
        if reverse_key in ABAk_inter:
            reverse_value=ABAk_inter[reverse_key] #Why not ABk_inter[reverse_key]??? Appears to just split weight between the two directions
        pairs[key][18]+=ABAk_inter[key]*1.0/(ABAk_inter[key]+reverse_value) #Why not divide by ABk_inter[key]+reverse_value??? Equivalent?
        reverse_value=0
        if reverse_key in ABk_inter:
            reverse_value=ABk_inter[reverse_key]
        pairs[key][19]+= ABAk_visit_inter[key]*1.0/(ABk_inter[key]*reverse_value)
    return

 
def scan_visit_info(visit_info, curbene):
    if len(visit_info)==0 or "SFROMDT" in visit_info: #If null or not a new record
        return
    path_list = referral_paths(visit_info, curbene) #Calls referral paths function
    edge_weights(path_list) #Calls edge-weights function
    return
    

def main():
    # fileName = input("Please enter the file name of the data you want to process: ")
    # print("You entered: " + fileName)
    # #fileName = "raw_cohort2016.txt"
    
    # year = fileName[10:14]
    # print("For the cohort of year " + year)
    
    # state = input("Please enter the state abbreviation you want to use to build the subnetwork: ")
    # print("You entered: " + state)
    
    filePath1 = "/projects/active/54054/idata/jomalley/raw-data_Sound.txt"
      
    filePath2 = '/projects/active/54054/idata/jomalley/'

    # save the HRR-specific results
    # filePath3 = '/projects/active/14593/Project3_HospClosure_VQI/programs/xran/networkResults/existence_updated/HRR-specific/'

    # for a specific HRR
    # hrrFile = input("Please enter the HRR file name: ")
    # fpw4 = open(filePath2 + hrrFile, "r")
    # hrrs = fpw4.readlines()
    
    # for h in hrrs:
    #     hrrID = h.strip('\\n')
        
    global fpw
    # underscore = '_'
    fpw=open(filePath2 + 'referral_path_2018a' + '.txt','w')
    fpw2=open(filePath2 + 'directed_network_2018a' + ".txt",'w')
    fpw3=open(filePath2 + 'undirected_network_2018a' + ".txt",'w')
      
    global pairs
    pairs = \{\}
    #undir_patient_num=\{\}
    #dir_patient_ranking_weight=\{\}
    #undir_patient_ranking_weight=\{\}
    visit_info=\{\}
    #output=\{\}
    #distinct_patients=\{\}

    count =0
    curbene="nobody" 
    
    # generate and save referral paths
    #fp=open("../../patient-phy-visit/2011/raw-data.txt",'r')
    
    fp=open(filePath1,'r')
    while 1:
        line = fp.readline()       
        if not line:
            break
        count = count+1
        if count%10000==0:
            print(count)

        words = line.split(";") #Fields within each line: [0] = BeneID, [2] = SFROMDT, [5] = PRFNPI, [20] = WorkRVU
        if words[0]!=curbene: 
            # processing the collected info of previous patient
            scan_visit_info(visit_info, curbene) #Call that leads to calls to referral path and edge_weight functions
            del visit_info
            visit_info=\{\}
        if len(str(words[20]))<1:
            words[20]=1.0 #Make RVU not less than 1?
        if words[2] not in visit_info:
            visit_info[words[2]]=[]
            visit_info[words[2]].append(str(words[5])+";"+str(words[20])) #Augment visit_info with date, NPI, RVU
        else:
            visit_info[words[2]].append(str(words[5])+";"+str(words[20])) 
        curbene = words[0]
    fp.close()
    # processing for the last patient 
    scan_visit_info(visit_info,curbene) 
    fpw.close() # finish generating referral paths and write results to file
    
    # save weights to directed networks
    for item in pairs:
        tmp=pairs[item]
        res=str(item)+";"
        for i in range(20): 
             res=res+str(round(tmp[i],4))+";"
        res=res+"\\n"
        fpw2.write(res)
        print("saving to directed network")
    fpw2.close()
    
    # save weights to undirected networks
    added_pairs=\{\}
    for item in pairs:
        tmp=pairs[item]
        phy_pair=item.split("_")
        reverse_key=phy_pair[1]+"_"+phy_pair[0]
        if item in added_pairs or reverse_key in added_pairs:
            continue
        added_pairs[item]=0
        added_pairs[reverse_key]=0
        reverse_tmp=[0 for col in range(20)]
        if reverse_key in pairs:
            reverse_tmp=pairs[reverse_key]
        res=str(item)+";"
        for i in range(20): #Creates sequence from 0 to 19
            if i%5==0:
                indicator=tmp[i]+reverse_tmp[i] #Existence case
                # if indicator>1:
                #     indicator=1 # do not binarize
                res=res+str(indicator)+";"
            else:
                res=res+str(round(tmp[i]+reverse_tmp[i], 4))+";" # add up weights in each direction
        res=res+"\\n"
        fpw3.write(res)
        print("saving to undirected network")
    fpw3.close()


if __name__ == "__main__":
    main()

# script for calculating 30 nw measures
# reciprocity for binary networkResults
# testing cases}