import os


## SYMM JUDGE

for filename in os.listdir('./symm'):
    filenam = filename[:-4].split('_')
    symm = filenam[0]
    molname = filenam[1]
    
 
