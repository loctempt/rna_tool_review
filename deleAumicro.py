import os 
import sys
inputpath=sys.argv
dirlist=os.listdir(inputpath)
for dir in dirlist:
    filenum = len(os.listdir(dir))
    if len(os.listdir(dir))<2:
        os.removedirs(inputpath+'/'+dir)
