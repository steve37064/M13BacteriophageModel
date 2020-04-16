import sys
import numpy as np
def degrade_DEF_eq(mRNA_str,D):
  eq_1 = mRNA_str + '() + '+ ' + '.join(["RBS" + str(i)+'()' for i in D]) +' -> 0 0.3*C12_'+mRNA_str+'*'
  eq_2 = '*'.join(["RBS" + str(i)+ 'Removal()' for i in D])
  return '  '+eq_1+eq_2+'\n'

def degrade_DEF_eq2(mRNA_str,D):
  eq_1 = mRNA_str + '() + '+ ' + '.join(["RBS" + str(i)+'()' for i in D]) +' -> 0 C12_'+mRNA_str+'*'
  eq_2 = '*'.join(["RBS" + str(i)+ 'Removal()' for i in D])
  return '  '+eq_1+eq_2+'\n'

def updateELRate(mRNA,species,lengths):
  #find which  position corresponds to given mrna
  length=np.sum(lengths[1,np.where(np.in1d(lengths[0], species))[0]])
  rate=67/(length-50)
  rate= round(rate,5)
  rate = '  C10_' +mRNA +'   '+ str(rate)+'   #Units: [(Second)^{-1}] -------------- RNAP Elongation of mRNA '+mRNA+'\n' 
  return rate
def updateELeq(mRNA,species):
  EL = ["RBS" + str(i)+'()' for i in species]
  EL = ' + '.join(EL)
  eq = '  EL'+mRNA+'() -> RNAP() + ' + mRNA +'() + '+ EL + ' C10_'+ mRNA+'\n'
  return eq
def degrade_mRNA_eq(mRNA,species):
 
  removed_genes = species[:2] 
  eq = '  '+mRNA+'() + RBS' + str(removed_genes[0])+'() -> D() C12_A*RBS' + str(removed_genes[0])+'Removal()\n'  
  return eq
#for lines 23: 32 
def updateDegradationRate(mRNA,species,lengths):
  length=np.sum(lengths[1,np.where(np.in1d(lengths[0], species))[0]])
  coeff = [ 1.70684029e-03, -6.96366259e+00]
  rate = np.exp(coeff[1])*np.exp(coeff[0]*length)
  rate = round(rate, 5) 
  rate = '  C12_'+ mRNA+ '   '+ str(rate)+'     #Units: [(Second)^{-1}] -------------- Degradation Rate of mRNA '+mRNA+'\n'
  return rate

def createModel(input_file,swap1=0,swap2=0):
  with open(input_file) as file:  
    model = file.readlines() 
  A=[2,5,7,9,8]
  B=[10,5,7,9,8]
  H=[9,8]
  Z=[3,6]
  Y=[3,6,1,11]
  W=[4]
  D=[5,7,9,8]
  species=[A,B,H,Z,Y,W,D]
  if(swap1 !=0 and swap2!=0):
    for i in range(len(species)):
      try:
        idx1=species[i].index(swap1)
      except ValueError:
        idx1=-1
      try:
        idx2=species[i].index(swap2)
      except ValueError:
        idx2=-1
      if(idx1 !=-1 and idx2 !=-1):
        species[i][idx1]=swap2
        species[i][idx2]=swap1
      elif(idx1 != -1 and idx2==-1):
        species[i][idx1]=swap2
      elif(idx1==-1 and idx2!=-1):
        species[i][idx2]=swap1 
      else: 
        None
  A=np.array(species[0])
  B=np.array(species[1])
  H=np.array(species[2])
  Z=np.array(species[3])
  Y=np.array(species[4])
  W=np.array(species[5])
  D=np.array(species[6])
  lengths=np.array([(2,10,5,7,9,8,3,6,1,11,4),(1200,350,270,101,98,220,1300,340,1040,310,1270)])
  species=np.array([A,B,H,Z,Y,W])
  E=D
  F=D
  G=H
  ###########write elongation rates (lines 17:21)
  mRNAs=['A','B','H','W','Z']
  lengths=np.array([(2,10,5,7,9,8,3,6,1,11,4),(1200,350,270,101,98,220,1300,340,1040,310,1270)])
  species=np.array([A,B,H,W,Z])
  #print('Elongation rates:')
  idx=16
  for m in range(len(mRNAs)):
    model[idx]=updateELRate(mRNAs[m],species[m],lengths)
    idx=idx+1
  ##########write degradation rates (lines 23:32):
  mRNAs=['A','B','D','E','F','G','H','W','Y','Z']
  species=np.array([A,B,D,E,F,G,H,H,H,H]) #what to do about e f g ???
  idx=22
  for m in range(len(mRNAs)):
    #print(updateDegradationRate(mRNAs[m],species[m],lengths))
    model[idx]=updateDegradationRate(mRNAs[m],species[m],lengths)
    idx=idx+1
  ###########rewrite active transcription eqs: 394:400

  mRNAs=['A','B','H','Z','Y','W']
  species=np.array([A,B,H,Z,Y,W])
  idx=393
  for m in range(len(mRNAs)):
    species_temp=np.delete(species[m], np.argwhere(species[m] == 7))
    if(mRNAs[m]=='Z'):
      #print('EL'+mRNAs[m]+'() -> RNAP() + ' + ' + '.join(["RBS" + str(i)+'()' for i in species[m]])+ ' C10_Z')
      #print('0    -> Z() C11*C10_Z*ELZ()')
    
      model[idx]='EL'+mRNAs[m]+'() -> RNAP() + ' + ' + '.join(["RBS" + str(i)+'()' for i in species_temp])+ ' C10_Z\n'
      idx=idx+1
      model[idx]='0    -> Z() C11*C10_Z*ELZ()\n'
    elif(mRNAs[m]=='Y'):
      #print('0    -> Y() + '+ ' + '.join(["RBS"+str(i)+'()' for i in species[m][2::]]) +' (1-C11)*C10_Z*ELZ()')
      model[idx]='0    -> Y() + '+ ' + '.join(["RBS"+str(i)+'()' for i in species_temp[2::]]) +' (1-C11)*C10_Z*ELZ()\n'
    else:
      #print(updateELeq(mRNAs[m],species[m]))
      model[idx]=updateELeq(mRNAs[m],species_temp)
    idx=idx+1

  #############write degredation eqs: 406:420

  #removed_genes = species[0][:2] 
  #D=species[0][2:]
  D=np.delete(D, np.argwhere(D == 7))
  #print('A() + RBS' + str(removed_genes[0])+'() -> D() C12_A*RBS' + str(removed_genes[0])+'Removal()' ) #line 406
  model[405]='A() + RBS' + str(A[0])+'() -> D() C12_A*RBS' + str(A[0])+'Removal()\n'
  #print('B() + RBS' + str(removed_genes[1])+'() -> D() C12_B*RBS' + str(removed_genes[1])+'Removal()' )#line 407
  model[406]='B() + RBS' + str(B[0])+'() -> D() C12_B*RBS' + str(B[0])+'Removal()\n'
  decay_products=['D','E','F','G','H']
  species=np.array([D,E,F,G,H])
  idx=407
  for i in range(len(decay_products)):
    s=np.delete(species[i], np.argwhere(species[i] == 7))
    if decay_products[i]=='F':
      #print(decay_products[i]+'() RBS'+str(D[0])+'() -> '+decay_products[i+1]+ '() 0.7*C12_'+decay_products[i]+'*RBS'+str(D[0])+'Removal()')
      model[idx]=decay_products[i]+'() + RBS'+str(D[0])+'() -> '+decay_products[i+1]+ '() 0.7*C12_'+decay_products[i]+'*RBS'+str(D[0])+'Removal()\n'
      idx=idx+1
    if decay_products[i]=='H':
      model[idx]=degrade_DEF_eq2(decay_products[i],s)
    if decay_products[i]!='H' and decay_products[i]!='F':
      model[idx]= decay_products[i]+'()   -> '+decay_products[i+1]+ '() 0.7*C12_'+decay_products[i]+'\n'
      idx=idx+1
      #print(decay_products[i]+'()   -> '+decay_products[i+1]+ '() 0.7*C12_'+decay_products[i])
    if decay_products[i]!='H':
      model[idx]=degrade_DEF_eq(decay_products[i],s)
      idx=idx+1
  #now do W:
  #print('W() + RBS'+str(W[0])+'() -> 0 C12_W*RBS'+str(W[0])+'Removal()')
  model[417]='W() + RBS'+str(W[0])+'() -> 0 C12_W*RBS'+str(W[0])+'Removal()\n'
  #do Y:
  model[418]='Y() + '+ ' + '.join(["RBS" + str(i)+'()' for i in Y])+' -> 0 C12_Y*'+ '*'.join(["RBS" + str(i)+ 'Removal()' for i in Y])+'\n'
  #print('Y() + '+ ' + '.join(["RBS" + str(i)+'()' for i in Y])+' -> 0 C12_Y*'+ '*'.join(["RBS" + str(i)+ 'Removal()' for i in Y]))
  #do Z:
  model[419]='Z() + ' +' + '.join(["RBS" + str(i)+'()' for i in Z])+' -> 0 C12_Z*'+ '*'.join(["RBS" + str(i)+ 'Removal()' for i in Z])+'\n'
  #print('Z() + ' +' + '.join(["RBS" + str(i)+'()' for i in Z])+' -> 0 C12_Z*'+ '*'.join(["RBS" + str(i)+ 'Removal()' for i in Z]))
  f= open(input_file,"w+")
  for i in range(len(model)):
    f.write(model[i])
  f.close() 

if __name__=="__main__":
    input_file=sys.argv[1]

    if(len(sys.argv)>2):
      createModel(input_file,swap1=int(sys.argv[2]),swap2=int(sys.argv[3]))
    else:
      createModel(input_file)

