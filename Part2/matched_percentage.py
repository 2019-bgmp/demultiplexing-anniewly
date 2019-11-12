#!/usr/bin/env python3


with open("/Users/annie/Documents/BGMP/Bi622/Assignments/demultiplexing-anniewly/Part2/matched_R1.txt") as fh :        
#    while True:
    for i in range(24):
        line = fh.readline().split()
        bar = (line[1].split("_")[1].split(".")[0])
#        bar = str(bar)
#        print(bar)
#        print(line[1])
        print("The percentage of " + str(bar) + " is " + str("%.2f"%(int(line[0])/1275357552*100)))
#        x=([int(i)/1275357552 if i.isdigit() for i in line.split() if i.isdigit()])
#
#        if line == '':
#            break