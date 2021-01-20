from pymongo import MongoClient                #creating DB
democlient = MongoClient()
myclient = MongoClient('localhost',27017)
mydb = myclient["geneinfo"]
mycoll=mydb["genetable"]
import sys
import copy
import re
mylist=[]
elementno=0
fin=open('genes.txt','r')
line=fin.readline()
while len(line)>0:
    mylist.append({'Start':0,'End':0,'Strand':'','Length':0,'PID':'','Gene':'','Synonym':'','Code':'','COG':'','Product':''})
    line=re.sub(r'\s+',' ',line).strip()
    #print(line)
    #line=' '.join(line.split())
    i=0
    start=0
    for i in range(len(line)):
        if line[i]>='0' and line[i]<='9':
            start=start*10+int(line[i],10)
        elif line[i]=='.':
            if line[i+1]=='.':
                i+=2
                break
            else:
                print('STARTERROR!')
                sys.exit()
        else:
            print('STARTERROR!')
            sys.exit()
    mylist[elementno]['Start']=start
    end=0
    while(1):
        if line[i]>='0' and line[i]<='9':
            end=end*10+int(line[i],10)
        elif line[i]==' ':
            break
        else:
            print('ENDERROR!')
            sys.exit()
        i+=1
    mylist[elementno]['End']=end
    while(1):
        if line[i]==' ':
            pass
        elif line[i]=='+' or line[i]=='-':
            mylist[elementno]['Strand']=line[i]
            i+=1
            break
        else:
            print('STRANDERROR!')
            sys.exit()
        i+=1
    length=0
    while(1):
        if line[i]==' ' and length==0:
            pass
        elif line[i]>='0' and line[i]<='9':
            length=length*10+int(line[i],10)
        elif length>0 and line[i]==' ':
            break
        else:
            print(length)
            print('LENGTHERROR!')
            sys.exit()
        i+=1
    mylist[elementno]['Length']=length
    pid=''
    while(1):
        if line[i]==' ' and pid=='':
            pass
        elif line[i]>='0' and line[i]<='9':
            pid+=line[i]
        elif pid!='' and line[i]==' ':
            break
        else:
            print('PIDERROR!')
            sys.exit()
        i+=1
    mylist[elementno]['PID']+=pid
    gene=''
    while(1):
        if line[i]==' ' and gene=='':
            pass
        elif gene!='' and line[i]==' ':
            break
        else:
            gene+=line[i]
        i+=1
    mylist[elementno]['Gene']+=gene
    synonym=''
    while(1):
        if line[i]==' ' and synonym=='':
            pass
        elif synonym!='' and line[i]==' ':
            break
        else:
            synonym+=line[i]
        i+=1
    mylist[elementno]['Synonym']+=synonym
    code=''
    while(1):
        if line[i]==' ' and code=='':
            pass
        elif code!='' and line[i]==' ':
            break
        else:
            code+=line[i]
        i+=1
    mylist[elementno]['Code']+=code
    cog=''
    while(1):
        if line[i]==' ' and cog=='':
            pass
        elif cog!='' and line[i]==' ':
            break
        else:
            cog+=line[i]
        i+=1
    mylist[elementno]['COG']+=cog
    product=''
    while i<len(line):
        if line[i]==' ' and product=='':
            pass
        else:
            product+=line[i]
        i+=1
    mylist[elementno]['Product']+=product
    line=fin.readline()
    elementno+=1

x = mycoll.insert_many(mylist)

dblist = myclient.list_database_names()
if input('Enter DB') in dblist:
    print('The database exists.')
else:
    print('Not Present')

print(dblist)
