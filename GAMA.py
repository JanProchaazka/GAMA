#################################################
#       Global Affine Mutiple Alignment         #
#             (for more strings)                #
#                                               #
#                           Jan Prochazka, 2007 #
#################################################


import copy
import parameters
import grey
#################################################

def buildArray(dim, size, value=None):
# creates an array of the size size in each of the dim dimensions
    if (dim!=len(size)) or (dim<=0):
        print "!!! Error in params !!!"
        return parameters.NaN
    ret=[]
    if dim>1:
        item=buildArray(dim-1, size[1:],value)
        for i in range(size[0]):
            ret=ret+[copy.deepcopy(item)]
        return ret
    else:   #dim==1
        for i in range(size[0]):
            ret=ret+copy.deepcopy([value])
        return ret

def sort_by_len(ar=[]):
    if len(ar)<=1:
        return ar
    pivot=ar[0]
    lower=sort_by_len(filter(lambda x: len(x) < len(pivot), ar[1:]))
    bigger=sort_by_len(filter(lambda x: len(x) >= len(pivot), ar[1:]))
    return lower+[pivot]+bigger

#################################################

def score2_opt(v,w): #score of global affinne mulitple alignment of 2 strings (with memory saving)
    if(len(v)<len(w)):  #A - shorter, B - longer
        A=v
        B=w
    else:
        A=w
        B=v

    L0=buildArray(1,[len(A)+1],-parameters.inf)   # main
    L1=buildArray(1,[len(A)+1],-parameters.inf)   # insert B
    L2=buildArray(1,[len(A)+1],-parameters.inf)   # insert A
    L0[0]=0

    for i in xrange(len(B)+1):     #for each row
        prevL0=-parameters.inf
        for j in xrange(len(A)+1):    #for each col
            #L1
            L1[j]=max(L1[j]-parameters.sigma, L0[j]-(parameters.ro+parameters.sigma))
            #L2
            if(j==0):
                L2[j]=-parameters.inf
            else:
                L2[j]=max(L2[j-1]-parameters.sigma,L0[j-1]-(parameters.ro+parameters.sigma))
            #L0
            old=L0[j]
            if(j*i==0):
                val=-parameters.inf
            else:
                val=prevL0+parameters.delta([A[j-1],B[i-1]])
            if(i+j!=0):
                L0[j]=max(val,L1[j],L2[j])
            prevL0=old
#        print L0

    if(len(L0)>0):
        return L0[-1]
    else:
        if(len(B)==0):
            return 0
        else:
            return -(parameters.ro+len(B)*parameters.sigma)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

def score(strings): #score of global afinne multiple alignment of strings

    n=len(strings)  #dimensions

    if(n==0):
        return 0
    #end if

    if(n==1):
        ret=0
        for i in xrange(len(strings[0])):
            ret+=parameters.delta([strings[0][i]])
        #next i
        return ret
    #end if

    #so n>=2 ..
################################################################
    #the longest set as the last one (saves memory)
    sort_by_len(strings) ;

    #initialize
    lenghts=[]    #lenghts of strings w/o the last one, i.e. searching object size (L)
    dimensions=[]
    for i in xrange(n):
        lenghts+=[len(strings[i])]
        dimensions+=[lenghts[-1]+1]
    #next i

    layers=range(1,(2**n)-1)
    layersCount=len(layers)

    #layers:
    L0=buildArray(n,dimensions,-parameters.inf) #init all -oo

    zeros=buildArray(1,[len(dimensions)],0)    #init origin zero
    writeTo(L0,zeros,0)

    L_I=buildArray(1+n,[layersCount]+dimensions,-parameters.inf)
    L_D=buildArray(1+n,[layersCount]+dimensions,-parameters.inf)

    ones=buildArray(1,[n],1)

    x=buildArray(1,[len(dimensions)],0)
    x=next_pos(x,dimensions)
    while(x[-1]>=0):  #for each point then returns [0,...0,-1]
        #INSERT LAYERS:
        for lay in xrange(layersCount):         #for each insert layer..
            lay_kod=grey.dec2bin(layers[lay])
            lay_kod=([0]*(n-len(lay_kod)))+lay_kod #insert leading zeros

            #for each element x (coord in L)
            prev=read(L_I[lay],diff(x,lay_kod))-parameters.sigma
            main=read(L0,diff(x,lay_kod))-(ones(lay_kod)*parameters.ro+parameters.sigma)
            val=maximum([prev, main])
            writeTo(L_I[lay], x, val)
        #next lay

        #DELETE LAYERS:
        for lay in xrange(layersCount):         #for each delete layer..
            lay_kod=grey.dec2bin(layers[lay])
            lay_kod=([0]*(n-len(lay_kod)))+lay_kod #insert leading zeros
            ley_kod=grey.invert(lay_kod)

            #for each element x (coord in L)
            prev=read(L_D[lay],diff(x,lay_kod))-parameters.sigma
            main=read(L0,diff(x,lay_kod))-(ones(lay_kod)*parameters.ro+parameters.sigma)
            val=maximum([prev, main])
            writeTo(L_D[lay], x, val)
        #next lay

        #MAIN LAYER:
        #which chars we try to MATCH / REPLACE(MISMATCH)
        chars=[]
        outOfRange=False
        for j in xrange(n):
            if((0<x[j]) and (x[j]<=len(strings[j]))):
                chars+=[strings[j][x[j]-1] ]
            else:
                outOfRange=True
            #end if
        #next j
        if(outOfRange):
            chars=[]
        #end if

        #"end of insert" possibility
        neighbours=[]
        for j in xrange(layersCount):
            neighbours+=[read(L_I[j],x)]
            neighbours+=[read(L_D[j],x)]
        #next j

        if(outOfRange):
            prev=-parameters.inf
        else:
            prev=read(L0,diff(x,ones))+parameters.delta(chars)
        #end if

        val=maximum([prev]+neighbours)
        writeTo(L0,x,val)

        x=next_pos(x,dimensions)
    #wend

    ret=L0
    for j in xrange(n):
        ret=ret[-1]

#    print "L0:\n---\n",L0,"\n---"
#    print "L_I:\n----\n",L_I,"\n----"
#    print "L_D:\n----\n",L_D,"\n----"
    return ret
#\def

#########################
def maximum(ar):
    if(len(ar)==0):
        print "Error: maximumm of empty list!"
        return 0

    x=ar[0]
    for i in xrange(1,len(ar)):
        if(x<ar[i]):
            x=ar[i]
    return x

def read(ar,coords):
    if(coords[0]<0): #negative coord -> predecessor of 0 -> impossible ~ -oo
        return -parameters.inf
    if(len(coords)>1):
        return read(ar[coords[0]],coords[1:])
    if(len(coords)==1):
        return ar[coords[0]]
#    print "Bad coord"
    return ar

def writeTo(ar,coords,value):
    if(len(coords)>1):
        return writeTo(ar[coords[0]],coords[1:],value)
    if(len(coords)==1):
        ar[coords[0]]=value
        return True
    print "Bad coord"
    return False

def next_pos(x, lenghts):
    if(x==[]):
        return [-1]
    val=x[0]+1
    if(val >= lenghts[0]):
        return [0]+next_pos(x[1:],lenghts[1:])
    else:
        return [val]+x[1:]

def diff(a,b):    #a-b for a,b vectors
    ret=[]
    if(len(a)!=len(b)):
        print "Subtraction of vectors of different lenghts!!!"
        return ret
    for i in xrange(len(a)):
        ret+=[a[i]-b[i]]
    return ret

def p(ar):    # max(#0,#1)
    zeros=0
    ones=0
    for i in ar:
        if(i==0):
            zeros+=1
        else:
            ones+=1 //TODO: sizeof(ar)-zeros
    return max(zeros, ones)
#    return max(sum(ar),len(ar)-sum(ar)) //equivalent notation (for 0,1 vectors)

def ones(ar):
    return sum(ar)
#    return 1


################################################################################

print "debug:"
#print score([])
#print score(["aabbaa"])
#print score(["aabbaa","aabaa"])
#print score(["aaaa","aabaa","aabbaa"])
#print score2_opt("bcbaca","baa"),"==", score(["bcbaca","baa"])
#print score2_opt("abca","aabbaca"),"==",score(["abca","aabbaca"])
#print score2_opt("aba","abcca"),"==",score(["aba","abcca"])
#print score(["abbba","abba","aa"])

strs=["abcd","abcd","xxxx"]
print "score(",strs,")=", score(strs)
strs=["aaaa","bbbb","xxxx"]
print "score(",strs,")=", score(strs)
