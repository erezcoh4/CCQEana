'''
    1. do 
        > jobsub_q --user=ecohen > held_jobs.txt
        > scp $uboone:/uboone/app/users/ecohen/ErezCCQE/held_jobs.txt ~/Desktop/
    2. copy to jobs.txt
    3. run 
        > python mac/pealoff_jobid_for_jobsub_rm.py
    4. copy and excecute result in gpvm
    '''

filename = "/Users/erezcohen/Desktop/held_jobs.txt"
string = ''
with open(filename,'r') as file:
    lines = file.readlines()
    for line in lines:
#        list=line.split("            ")
        list=line.split("          ")
        jobid = list[0]
        print 'adding',jobid
        string = string+jobid+','
print 'jobsub_rm --jobid='+string