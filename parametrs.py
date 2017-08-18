NaN="NaN"
inf=100
ro=1
sigma=1

def delta(elems):
    if(len(elems)==0):
            return 0
    if(len(elems)==1):
            return 1
    if(len(elems)==2):
            if(elems[0]==elems[1]):
                    return 1
            else:
                    return -1
    ##################
    dic={}
    n=len(elems)
    for p in elems:
        if(dic.has_key(p)):
            dic[p]+=1
        else:
            dic[p]=1
    #calculation aleternativs:
    #1) #chars/count of different chars .. + normalization
    #return 1.0*n/len(dic) - n/2.0
    #2) same = 1, different = -1
    #3)
    if(len(dic)==1):
    	return 1
    else:
    	return -1
