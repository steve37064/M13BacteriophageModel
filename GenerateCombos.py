from createModel import * 
from itertools import combinations 
import numpy as np 
import os

input_file = "Model.bngl"

Folder = "GeneratedModels/"


Name = "RemakeOrginal"
BNGLOut = Folder+Name
createModel(input_file,output_file=BNGLOut+".bngl")
os.system("perl /Applications/RuleBender.app/Contents/eclipse/BioNetGen-2.5.0/BNG2.pl "+BNGLOut+".bngl")
os.system("rm *.cdat")
os.system("rm *.gdat")
os.system("rm *.net")
os.system("mv  "+Name+".m "+ Folder+"MatlabVersions/" +"." )

ProteinSwaps = np.arange(1,12)
for i_swap,j_swap in combinations(ProteinSwaps,2): 
    Name = "Swap_" + str(i_swap)+"_"+str(j_swap)
    BNGLOut = Folder+Name
    createModel(input_file,i_swap,j_swap,output_file=BNGLOut+".bngl")
    os.system("perl /Applications/RuleBender.app/Contents/eclipse/BioNetGen-2.5.0/BNG2.pl "+BNGLOut+".bngl")
    os.system("rm *.cdat")
    os.system("rm *.gdat")
    os.system("rm *.net")
    os.system("mv  "+Name+".m "+ Folder+"MatlabVersions/" +"." )

os.system("perl -pi -w -e 's/   0.0001,   .../   1e-8,   .../g;' "+Folder+"MatlabVersions/*.m")