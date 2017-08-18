#########################################################
# Grey-code manipulation library
# Jan Prochazka, 2007
#########################################################
## DEC vs. GREY

def dec2grey(c):
    #converts decimal number to Grey code
    rad=2
    ret=[]
    tmp=2
    while(tmp>=1.5):
        tmp=(c*1.0)/rad + 0.5
        if(int(tmp)%2 == 0):
            ret=[0]+ret
        else:
            ret=[1]+ret
        rad*=2
    if(c==0):
        ret=[0]
    return ret

def grey2dec(g):
    #converts Grey code to decimal number
    ret=0
    tmp=1
    sgn=1
    for i in range(len(g)):
        if(g[-1-i]==1): #from the end
            ret+=sgn*tmp
            sgn=-sgn
        tmp=((tmp+1)*2)-1 #1,3,7,15,..,[(2**n)-1]
    return abs(ret)

#########################################################
## BIN vs. GREY

def bin2grey(b):
    #converts binary notation of a number to its Grey code
    arr=[0]+b
    ret =[]
    for i in range(1,len(arr)):
        if(arr[i-1]==arr[i]):
            ret+=[0]
        else:
            ret+=[1]
    return ret

def grey2bin(g):
    #converts Grey code to its binary notated number
    arr=g
    ret=[0]
    for i in range(len(arr)):
        if(arr[i]==ret[-1]): #ret[-1] ~ last. subresult
            ret+=[0]
        else:
            ret+=[1]
    return ret[1:]

#########################################################
## BIN vs. DEC

def bin2dec(b):
    #converts binary number to decimal
    ret=0
    tmp=1
    for i in xrange(len(b)):
        if(b[-1-i]==1):
            ret+=tmp
        tmp*=2
    return ret

def dec2bin(d):
    #converts decimal number to binary
    ret=[]
    while (d>0):
        ret=[d%2]+ret
        d=d/2
    if (ret==[]):
        ret=[0]
    return ret

#########################################################

def invert(b):
	for i in xrange(len(b)):
		if(b[i]==0):
			b[i]=1
		else:
			b[i]=0
	return b
#########################################################
## debug:

#def test(od,po):
#    err=False
#    for i in xrange(od,po):
#        if(i!=grey2dec(dec2grey(i))):
#            err=True
#            print "Doesn't work for", i
#    return err
#
#a=0
#b=100000
#print "TEST of interval [",a,";",b,"]:"
#if (test(a,b)):
#    print "TEST finished unsuccessfully!"
#else:
#    print "TEST finished successfully!"
#
#####

#####
#for i in xrange(40000):
#    a=bin2grey(dec2bin(i))
#    b=dec2grey(i)
#    if(a!=b):
#        print i, a==b
#        print "  ", a, "I have:",b
#
#########################################################
