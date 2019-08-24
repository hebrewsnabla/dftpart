#!/usr/bin/env python
# -* - coding: UTF-8 -* -
#Filename:

def labc_parser(labc):
    flag = -1
    number = 0
    sys = []
    LAB_tem = []
    Lab = []
    with open(labc,'r') as f:
        for line in f:
            line = line.strip()
            if 'Table: Lab for subsystems' in line:
                flag = 0
            if flag ==0:
                number = number + 1
            if number==5:
                flag = 1
            if flag ==1:
                if '=' in line:
                    break
                line = line.split(" ")
                line = [e for e in line if e!='']
                sys.append(line[0])
                line[1] = line[1].split(",")
                Lab = []
                for i in range(len(line[1])):
                  if '-' in line[1][i]:
                      t = line[1][i].index('-')
                      delt = int(line[1][i][t+1:])-int(line[1][i][0:t])
                      delt = delt + 1
                      s = int(line[1][i][0:t])
                      line[1][i] = ''
                      for m in range(delt):
                        line[1][i] = line[1][i] + str(int(s+m))+ ' '
                      line[1][i]=line[1][i].split(' ')
                      line[1][i] = [e for e in line[1][i] if e!= '']
                      for j in range(len(line[1][i])):
                        Lab.append(int(line[1][i][j]))
                  elif '0' in line[1][i] :
                    if '*' in line[1][i] :
                      t=line[1][i].index('*')
                      line[1][i] = '0 '* int(line[1][i][0:t])
                      line[1][i]=line[1][i].split(' ')
                      line[1][i] = [e for e in line[1][i] if e!= '']
                      for j in range(len(line[1][i])):
                        Lab.append(int(line[1][i][j]))
                    else :
                       Lab.append(int(line[1][i]))

                  else :
                      Lab.append(int(line[1][i]))
            LAB_tem.append(Lab)
        LAB = []
        for item in LAB_tem:
            if item==[]:
                pass
                #LAB.pop(item)
            else:
                LAB.append(item)
                #print("item",item)
            #LAB.append(Lab)
        #print(LAB)
        nsys = len(sys)
        #LAB = np.hstack(LAB).reshape(nsys,-1)
        return LAB

#LAB = labc_parser(labcnam)
#print(LAB)
#print('Hello World')
