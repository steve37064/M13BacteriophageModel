from createModel import * 
import os

input_file = "Model.bngl"

Folder = "GeneratedModels/"


Name = "RemakeOrginal"
BNGLOut = Folder+Name
createModel(input_file,output_file=BNGLOut+".bngl")

Name = "Swap5_6"
BNGLOut = Folder+Name
createModel(input_file,5,6,output_file=BNGLOut+".bngl")

os.system("perl /Applications/RuleBender.app/Contents/eclipse/BioNetGen-2.5.0/BNG2.pl "+BNGLOut+".bngl")
os.system("rm *.cdat")
os.system("rm *.gdat")
os.system("rm *.net")
os.system("mv  "+Name+".m "+ Folder+"MatlabVersions/" +"." )