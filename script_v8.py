# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:12:48 2023

@author: armen
"""

from os import listdir,chdir
import subprocess
from random import randrange,sample
from shutil import copy,copytree,rmtree
from os import mkdir,remove
import csv
import time
import shutil
import math
import random

start_time = time.time()# start timing

while True:
    exptype = int(input("Which experiment do you want to run? (1/2)"))
    if exptype != 1 and exptype !=2:
        exptype = int(input("Which experiment do you want to run? (1/2)"))
    else:
        if exptype == 2:
            expnum = 5
        if exptype == 1:
            expnum = int(input("Number of Experiments: "))
            while True:
                dist = int(input("How will the the experiments be distributed in sequence (1) or geometrically (2)?"))
                if dist != 1 and dist !=2:
                    dist = int(input("How will the the experiments be distributed in sequence (1) or geometrically (2)?"))
                else:
                    if dist ==2:
                        start = int(input("How many reasources to start with?"))
                    break
        break
    

repnum= int(input("How many replicates do you want?"))
print(f"You want {repnum} replicates, each replicate will contain {expnum} for a total of {repnum*expnum}")
confirm = input("Is this correct? (y/n)")
if confirm == "y":
    for rep in range(repnum):
        copy('./default_files/avida.cfg','./avida.cfg')
        with open("avida.cfg", "r") as seedfile:
            contents = seedfile.read()
            
        random_seed = random.randint(1,1000000) 
        print(f"Starting replicate{rep+1}")
        
        new_contents = contents.replace("RANDOM_SEED -1","RANDOM_SEED "+str(random_seed))
        
        with open("avida.cfg", "w") as f:
            f.write(new_contents)
        if exptype == 2:
            for ex in range (expnum):
                num= ex + 1
                mkdir('experiment_'+str(num)) #create a directory named cdir
                copy('./default_files/food_webs/fw_'+str(num)+'/environment.cfg','./experiment_'+str(num)+'/environment.cfg') #copy one file from one place to another
                copy('./experiment_' + str(num) + '/environment.cfg', './')
                print("starting experiment_" + str(num))
                    
                proc = subprocess.Popen(['./avida'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                # Wait for the subprocess to finish
                output, error = proc.communicate()
                    
                # Close the subprocess
                proc.terminate()
                    
                shutil.move('./data','./experiment_' +str(num))
                remove('./environment.cfg')
            
            
            
            
        if exptype == 1:

            if dist == 1:
                for ex in range (expnum):
                    num= ex + 1   
                    mkdir('experiment_'+str(num)) #create a directory named cdir
                    copy('./default_files/environment.cfg','./experiment_'+str(num)+'/environment.cfg') #copy one file from one place to another
            
                    #WSIZE
                    env_file = open('./experiment_'+str(num)+'/environment.cfg','r')
                    env_file_content = []
                    for i in env_file:
                        env_file_content.append(i.split(' '))
            
                    #env_file_content = [['a','b'],['c','d']]
                    #access first line: env_file_content[0] ; note that index in python start from 0
                    #access first element in second line: env_file_content[1][0] ; note that index in python start from 0
            
                    n = num #number of resources
                    var0 = '100' #resource inflow
                    var1 = '0.01' #resource outflow
                    var3 = '1' # The minimum amount of resource required
                    reactiontype = ['not','nand','and','orn','or','andn','nor','xor','equ']
                    reward = ["1.0","1.0","2.0","2.0","4.0","4.0","8.0","8.0","16.0"]
            
            
                    #n = sample(range(10),1)[0]
            
                    out = open('./experiment_'+str(num)+'/environment.cfg','w')
                    for i in range(n):
                        out.write('RESOURCE res'+str(i)+':inflow='+var0+':outflow='+var1+'\n')
            
                    sc=0
                    for i in range(n):
                        out.write('REACTION reaction'+str(i)+' '+reactiontype[sc]+' process:resource=res'+str(i)+':value='+reward[sc]+':min='+var3+   '\n')
                        sc+=1
                        if sc==len(reactiontype):
                            sc = 0 
                    out.close()
                    env_file.close()
                
                    ##RUN Avida from python
                    copy('./experiment_' + str(num) + '/environment.cfg', './')
                    print("starting experiment_" + str(num))
                    
                    proc = subprocess.Popen(['./avida'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                    # Wait for the subprocess to finish
                    output, error = proc.communicate()
                    
                    # Close the subprocess
                    proc.terminate()
                    
                    shutil.move('./data','./experiment_' +str(num))
                    #copytree('./data', './experiment_' + str(num) + '/data')
                    #rmtree('./data')
                    remove('./environment.cfg')
                
            if dist ==2:
                
                limit = start
                ratio = 2
                for lim in range((expnum)-1):
                    limit*=2
                n = start    
                for ex in range (expnum):
                    num= ex + 1   
                    mkdir('experiment_'+str(num)) #create a directory named cdir
                    copy('./default_files/environment.cfg','./experiment_'+str(num)+'/environment.cfg') #copy one file from one place to another
            
                    #WSIZE
                    env_file = open('./experiment_'+str(num)+'/environment.cfg','r')
                    env_file_content = []
                    for i in env_file:
                        env_file_content.append(i.split(' '))
            
                    #env_file_content = [['a','b'],['c','d']]
                    #access first line: env_file_content[0] ; note that index in python start from 0
                    #access first element in second line: env_file_content[1][0] ; note that index in python start from 0
            

                    var0 = '100' #resource inflow
                    var1 = '0.01' #resource outflow
                    var3 = '1' # The minimum amount of resource required
                    reactiontype = ['not','nand','and','orn','or','andn','nor','xor','equ']
                    reward = ["1.0","1.0","2.0","2.0","4.0","4.0","8.0","8.0","16.0"]
            
            
                    #n = sample(range(10),1)[0]
            
                    out = open('./experiment_'+str(num)+'/environment.cfg','w')
                    for i in range(n):
                        out.write('RESOURCE res'+str(i)+':inflow='+var0+':outflow='+var1+'\n')
            
                    sc=0
                    for i in range(n):
                        out.write('REACTION reaction'+str(i)+' '+reactiontype[sc]+' process:resource=res'+str(i)+':value='+reward[sc]+':min='+var3+   '\n')
                        sc+=1
                        if sc==len(reactiontype):
                            sc = 0 
                    out.close()
                    env_file.close()
                
                    ##RUN Avida from python
                    copy('./experiment_' + str(num) + '/environment.cfg', './')
                    print("starting experiment_" + str(num))
                    
                    proc = subprocess.Popen(['./avida'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                    # Wait for the subprocess to finish
                    output, error = proc.communicate()
                    
                    # Close the subprocess
                    proc.terminate()
                    
                    shutil.move('./data','./experiment_' +str(num))
                    #copytree('./data', './experiment_' + str(num) + '/data')
                    #rmtree('./data')
                    remove('./environment.cfg')
                    n *=2

        replicatenum = rep + 1
        mkdir('replicate_'+str(replicatenum)) #create a directory named replicate_
        shutil.move('./avida.cfg','replicate_'+str(replicatenum))
        for repl in range(expnum):
            numb = repl + 1  
            source = './experiment_'+str(numb)
            dest = './replicate_' + str(replicatenum) + '/experiment_' + str(numb)
            shutil.move(source, dest)
  
    
#end timing
end_time = time.time()

elapsed_time = end_time - start_time
with open("elapsed_time.txt", "w") as file:
    file.write("Elapsed time: {} seconds".format(elapsed_time))
    