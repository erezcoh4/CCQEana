filename = "/Users/erezcohen/Desktop/held.c"
string = ''
with open(filename,'r') as file:
    lines = file.readlines()
    for line in lines:
        list=line.split("            ")
        jobid = list[0]
        print 'adding',jobid
        string = string+jobid+','
print string