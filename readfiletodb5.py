from pymongo import MongoClient                #creating DB
democlient = MongoClient()
myclient = MongoClient('localhost',27017)
mydb = myclient["geneinfo"]
mycoll=mydb["maintable"]
import sys
import copy
mylist=[]
fin=open('Eco2.txt','r')
firstline=fin.readline()
geneseq=[]
geneseq2=''
start=0
end=0
secondhalf=0
sino=0
#print('SI_No.          Information          Gene\n')

if firstline[:4]=='>gi|' and firstline[13:27]=='|gb|U00096.3|:':
  for i in firstline[4:13]:
    if i not in ['0','1','2','3','4','5','6','7','8','9']:
      print('ERROR in firstline position ' + i)
      sys.exit()
  j=27
  for i in firstline[27:]:
    if i=='-':
      j+=1
      break
    start=(start*10)+int(i,10)
    j+=1
  for i in firstline[j:]:
    if i==' ':
      break
    end=(end*10)+int(i,10)
    j+=1
  secondhalf=j

mylist.append({'SI No':sino+1,'Information':firstline,'Gene':'','A Count':0,'T Count':0,'C Count':0,'G Count':0,'Start':start,'End':end,'Length':0,'G+C percent':0,'AA Sequence':'','Effective codon no':0})
listelementno=0
sino+=1
#print(sino,firstline[1:].strip())
line=fin.readline()
if '\n' in line and len(line)==1:
  print('ERROR!')
  sys.exit()
linenumber=2
k=0
comma=0

while len(line)>0:
  if line[0]=='>':
    complement=0
    comma=0
    mylist[listelementno]['Gene']=geneseq2
    mylist[listelementno]['G+C percent']=100*(mylist[listelementno]['G Count']+mylist[listelementno]['C Count'])/(mylist[listelementno]['A Count']+mylist[listelementno]['T Count']+mylist[listelementno]['C Count']+mylist[listelementno]['G Count'])
    geneseq2=''
    if line[:27]==firstline[:27]:
      start=end=0
      j=27
      if line[j]=='c':
        complement=1
        j+=1
      for i in line[j:]:
        if i==',':
          comma=1
          break
        if i=='-':
          j+=1
          break
        start=(start*10)+int(i,10)
        j+=1
      for i in line[j:]:
        if i==',':
          comma=1
          break
        elif i==' ' or i=='-':
          break
        else:
          #print('endERROR!')
          pass
        end=(end*10)+int(i,10)
        j+=1
      sino+=1
      #print(sino,line[1:].strip())
    else:
      #print(line[:27])
      #print(firstline[:27])
      #print('notfirstlineERROR!')
      sys.exit()
    if comma==0:
      listelementno=listelementno+1
      mylist.append({'SI No':sino,'Information':line,'Gene':'','A Count':0,'T Count':0,'C Count':0,'G Count':0,'Start':0,'End':0,'Length':0,'G+C percent':0,'AA Sequence':'','Effective codon no':0})
      if complement==1:
        mylist[listelementno]['Start']=end
        mylist[listelementno]['End']=start
      else:
        mylist[listelementno]['Start']=start
        mylist[listelementno]['End']=end
  else:
    if comma==0:
      for i in line.strip():
        if i in ['A','T','C','G']:
          if i=='A':
            mylist[listelementno]['A Count']+=1
          if i=='T':
            mylist[listelementno]['T Count']+=1
          if i=='C':
            mylist[listelementno]['C Count']+=1
          if i=='G':
            mylist[listelementno]['G Count']+=1
          geneseq2=geneseq2+i
          geneseq.append(i)
          mylist[listelementno]['Length']+=1
          #print(i),
        else:
          print('ERROR!')
          print(i)
          print(geneseq)
          sys.exit()
    k+=1
  #print('\n')
  line=fin.readline()
  if '\n' in line and len(line)==1:
    #print('ERROR!')
    sys.exit()
  linenumber+=1

mylist[listelementno]['Gene']=geneseq2
mylist[listelementno]['G+C percent']=100*(mylist[listelementno]['G Count']+mylist[listelementno]['C Count'])/(mylist[listelementno]['A Count']+mylist[listelementno]['T Count']+mylist[listelementno]['C Count']+mylist[listelementno]['G Count'])

i=0
triplet=''
TTT=TTC=TTA=TTG=TCT=TCC=TCA=TCG=0
TAT=TAC=TGT=TGC=TGG=CTT=CTC=CTA=0
CTG=CCT=CCC=CCA=CCG=CAT=CAC=CAA=0
CAG=CGT=CGC=CGA=CGG=AGA=AGG=ATT=0
ATC=ATA=ATG=ACT=ACC=ACA=ACG=AAT=0
AAC=AAA=AAG=AGT=AGC=GTT=GTC=GTA=0
GTG=GCT=GCC=GCA=GCG=GAT=GAC=GAA=0
GAG=GGT=GGC=GGA=GGG=0
F=L=S=Y=C=W=L=P=H=Q=R=0
I=M=T=N=K=S=V=A=D=E=G=0
NC=0
FAA=0
while(i<len(mylist)):
  j=0
  while(j<len(mylist[i]['Gene'])):
    if len(triplet)<3:
      triplet=triplet+mylist[i]['Gene'][j]
      j+=1
    else:
      #print(triplet)
      if triplet in ['TTT','TTC']:
        mylist[i]['AA Sequence']+='F'
        F+=1
        if triplet=='TTT':
          TTT+=1
        if triplet=='TTC':
          TTC+=1
      elif triplet in ['TCT','TCC','TCA','TCG','AGT','AGC']:
        mylist[i]['AA Sequence']+='S'
        S+=1
        if triplet=='TCT':
          TCT+=1
        if triplet=='TCC':
          TCC+=1
        if triplet=='TCA':
          TCA+=1
        if triplet=='TCG':
          TCG+=1
        if triplet=='AGT':
          AGT+=1
        if triplet=='AGC':
          AGC+=1
      elif triplet in ['TAT','TAC']:
        mylist[i]['AA Sequence']+='Y'
        Y+=1
        if triplet=='TAT':
          TAT+=1
        if triplet=='TAC':
          TAC+=1
      elif triplet in ['TGT','TGC']:
        mylist[i]['AA Sequence']+='C'
        C+=1
        if triplet=='TGT':
          TGT+=1
        if triplet=='TGC':
          TGC+=1
      elif triplet in ['TGG']:
        mylist[i]['AA Sequence']+='W'
        W+=1
        if triplet=='TGG':
          TGG+=1
      elif triplet in ['CTT','CTC','CTA','CTG','TTA','TTG']:
        mylist[i]['AA Sequence']+='L'
        L+=1
        if triplet=='CTT':
          CTT+=1
        if triplet=='CTC':
          CTC+=1
        if triplet=='CTA':
          CTA+=1
        if triplet=='CTG':
          CTG+=1
        if triplet=='TTA':
          TTA+=1
        if triplet=='TTG':
          TTG+=1
      elif triplet in ['CCT','CCC','CCA','CCG']:
        mylist[i]['AA Sequence']+='P'
        P+=1
        if triplet=='CCT':
          CCT+=1
        if triplet=='CCC':
          CCC+=1
        if triplet=='CCA':
          CCA+=1
        if triplet=='CCG':
          CCG+=1
      elif triplet in ['CAT','CAC']:
        mylist[i]['AA Sequence']+='H'
        H+=1
        if triplet=='CAT':
          CAT+=1
        if triplet=='CAC':
          CAC+=1
      elif triplet in ['CAA','CAG']:
        mylist[i]['AA Sequence']+='Q'
        Q+=1
        if triplet=='CAA':
          CAA+=1
        if triplet=='CAG':
          CAG+=1
      elif triplet in ['CGT','CGC','CGA','CGG','AGA','AGG']:
        mylist[i]['AA Sequence']+='R'
        R+=1
        if triplet=='CGT':
          CGT+=1
        if triplet=='CGC':
          CGC+=1
        if triplet=='CGA':
          CGA+=1
        if triplet=='CGG':
          CGG+=1
        if triplet=='AGA':
          AGA+=1
        if triplet=='AGG':
          AGG+=1
      elif triplet in ['ATT','ATC','ATA']:
        mylist[i]['AA Sequence']+='I'
        I+=1
        if triplet=='ATT':
          ATT+=1
        if triplet=='ATC':
          ATC+=1
        if triplet=='ATA':
          ATA+=1
      elif triplet in ['ATG']:
        mylist[i]['AA Sequence']+='M'
        M+=1
        if triplet=='ATG':
          ATG+=1
      elif triplet in ['ACT','ACC','ACA','ACG']:
        mylist[i]['AA Sequence']+='T'
        T+=1
        if triplet=='ACT':
          ACT+=1
        if triplet=='ACC':
          ACC+=1
        if triplet=='ACA':
          ACA+=1
        if triplet=='ACG':
          ACG+=1
      elif triplet in ['AAT','AAC']:
        mylist[i]['AA Sequence']+='N'
        N+=1
        if triplet=='AAT':
          AAT+=1
        if triplet=='AAC':
          AAC+=1
      elif triplet in ['AAA','AAG']:
        mylist[i]['AA Sequence']+='K'
        K+=1
        if triplet=='AAA':
          AAA+=1
        if triplet=='AAG':
          AAG+=1
      elif triplet in ['GTT','GTC','GTA','GTG']:
        mylist[i]['AA Sequence']+='V'
        V+=1
        if triplet=='GTT':
          GTT+=1
        if triplet=='GTC':
          GTC+=1
        if triplet=='GTA':
          GTA+=1
        if triplet=='GTG':
          GTG+=1
      elif triplet in ['GCT','GCC','GCA','GCG']:
        mylist[i]['AA Sequence']+='A'
        A+=1
        if triplet=='GCT':
          GCT+=1
        if triplet=='GCC':
          GCC+=1
        if triplet=='GCA':
          GCA+=1
        if triplet=='GCG':
          GCG+=1
      elif triplet in ['GAT','GAC']:
        mylist[i]['AA Sequence']+='D'
        D+=1
        if triplet=='GAT':
          GAT+=1
        if triplet=='GAC':
          GAC+=1
      elif triplet in ['GAA','GAG']:
        mylist[i]['AA Sequence']+='E'
        E+=1
        if triplet=='GAA':
          GAA+=1
        if triplet=='GAG':
          GAG+=1
      elif triplet in ['GGT','GGC','GGA','GGG']:
        mylist[i]['AA Sequence']+='G'
        G+=1
        if triplet=='GGT':
          GGT+=1
        if triplet=='GGC':
          GGC+=1
        if triplet=='GGA':
          GGA+=1
        if triplet=='GGG':
          GGG+=1
      elif triplet in ['TAA','TAG','TGA']:
        pass
      else:
        print('ERROR!')
        sys.exit()
      triplet=''
  if triplet in ['TAA','TAG','TGA']:
    mylist[i]['AA Sequence']+='Stop'
  if F>0:
    FAA=(TTT/F)*(TTT/F)+(TTC/F)*(TTC/F)
    if FAA>0:
      NC+=(1/FAA)
  if S>0:
    FAA=(TCT/S)*(TCT/S)+(TCC/S)*(TCC/S)+(TCA/S)*(TCA/S)+(TCG/S)*(TCG/S)+(AGT/S)*(AGT/S)+(AGC/S)*(AGC/S)
    if FAA>0:
      NC+=(1/FAA)
  if Y>0:
    FAA=(TAT/Y)*(TAT/Y)+(TAC/Y)*(TAC/Y)
    if FAA>0:
      NC+=(1/FAA)
  if C>0:
    FAA=(TGT/C)*(TGT/C)+(TGC/C)*(TGC/C)
    if FAA>0:
      NC+=(1/FAA)
  if W>0:
    FAA=(TGG/W)*(TGG/W)
    if FAA>0:
      NC+=(1/FAA)
  if L>0:
    FAA=(CTT/L)*(CTT/L)+(CTC/L)*(CTC/L)+(CTA/L)*(CTA/L)+(CTG/L)*(CTG/L)+(TTA/L)*(TTA/L)+(TTG/L)*(TTG/L)
    if FAA>0:
      NC+=(1/FAA)
  if H>0:
    FAA=(CAT/H)*(CAT/H)+(CAC/H)*(CAC/H)
    if FAA>0:
      NC+=(1/FAA)
  if P>0:
    FAA=(CCT/P)*(CCT/P)+(CCC/P)*(CCC/P)+(CCA/P)*(CCA/P)+(CCG/P)*(CCG/P)
    if FAA>0:
      NC+=(1/FAA)
  if Q>0:
    FAA=(CAA/Q)*(CAA/Q)+(CAG/Q)*(CAG/Q)
    if FAA>0:
      NC+=(1/FAA)
  if R>0:
    FAA=(CGT/R)*(CGT/R)+(CGC/R)*(CGC/R)+(CGA/R)*(CGA/R)+(CGG/R)*(CGG/R)+(AGA/R)*(AGA/R)+(AGG/R)*(AGG/R)
    if FAA>0:
      NC+=(1/FAA)
  if I>0:
    FAA=(ATT/I)*(ATT/I)+(ATC/I)*(ATC/I)+(ATA/I)*(ATA/I)
    if FAA>0:
      NC+=(1/FAA)
  if M>0:
    FAA=(ATG/M)*(ATG/M)
    if FAA>0:
      NC+=(1/FAA)
  if T>0:
    FAA=(ACT/T)*(ACT/T)+(ACC/T)*(ACC/T)+(ACA/T)*(ACA/T)+(ACG/T)*(ACG/T)
    if FAA>0:
      NC+=(1/FAA)
  if N>0:
    FAA=(AAT/N)*(AAT/N)+(AAC/N)*(AAC/N)
    if FAA>0:
      NC+=(1/FAA)
  if K>0:
    FAA=(AAA/K)*(AAA/K)+(AAG/K)*(AAG/K)
    if FAA>0:
      NC+=(1/FAA)
  if V>0:
    FAA=(GTT/V)*(GTT/V)+(GTC/V)*(GTC/V)+(GTA/V)*(GTA/V)+(GTG/V)*(GTG/V)
    if FAA>0:
      NC+=(1/FAA)
  if A>0:
    FAA=(GCT/A)*(GCT/A)+(GCC/A)*(GCC/A)+(GCA/A)*(GCA/A)+(GCG/A)*(GCG/A)
    if FAA>0:
      NC+=(1/FAA)
  if D>0:
    FAA=(GAT/D)*(GAT/D)+(GAC/D)*(GAC/D)
    if FAA>0:
      NC+=(1/FAA)
  if E>0:
    FAA=(GAA/E)*(GAA/E)+(GAG/E)*(GAG/E)
    if FAA>0:
      NC+=(1/FAA)
  if G>0:
    FAA=(GGT/G)*(GGT/G)+(GGC/G)*(GGC/G)+(GGA/G)*(GGA/G)+(GGG/G)*(GGG/G)
    if FAA>0:
      NC+=(1/FAA)
  mylist[i]['Effective codon no']=NC
  i+=1
  FAA=0
  NC=0
  TTT=TTC=TTA=TTG=TCT=TCC=TCA=TCG=0
  TAT=TAC=TGT=TGC=TGG=CTT=CTC=CTA=0
  CTG=CCT=CCC=CCA=CCG=CAT=CAC=CAA=0
  CAG=CGT=CGC=CGA=CGG=AGA=AGG=ATT=0
  ATC=ATA=ATG=ACT=ACC=ACA=ACG=AAT=0
  AAC=AAA=AAG=AGT=AGC=GTT=GTC=GTA=0
  GTG=GCT=GCC=GCA=GCG=GAT=GAC=GAA=0
  GAG=GGT=GGC=GGA=GGG=0
  F=L=S=Y=C=W=L=P=H=Q=R=0
  I=M=T=N=K=S=V=A=D=E=G=0

x = mycoll.insert_many(mylist)

dblist = myclient.list_database_names()
if input('Enter DB') in dblist:
    print('The database exists.')
else:
    print('Not Present')

print(dblist)
