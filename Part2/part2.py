import numpy as np
import gzip

def convert_phred(string):
    """Converts a string of characters into phred scores"""
    return [ord(x)-33 for x in string]

complement={"A":"T","T":"A","G":"C","C":"G","N":"N"}
def rev_complement(string):
    return ''.join([complement[x] for x in string[::-1]])

def cutoff(q1,q2,q3,q4):
    """Return True if the average score of the biological reads is below 25 or if any base of the indexes is below 10"""
    if np.mean(convert_phred(q1))<25:
        return True
    if np.mean(convert_phred(q4))<25:
        return True
    if min(convert_phred(q2)) <10:
        return True
    if min(convert_phred(q3)) <10:
        return True
    return False

def check_match(r2,r3):
    """Return True if two things are the same otherwise return false"""
    if r2 == r3:
        return True
    return False

 
## Read in indexes
barcodes=[]
with open('/projects/bgmp/shared/2017_sequencing/indexes.txt') as fb:
    [barcodes.append(x.split()[4])for x in fb.readlines()[1:]]
print(barcodes)

#open all file that is going to be written in and make a dictionary
#dictionary that has barcodes as key and file as value
## start a dictionary to count reads for each index
filedict_R1={}
filedict_R2={}
barcounts = {}
for i in barcodes:
    filedict_R1.update({i:open("R1_"+str(i)+".fq","a")})
    filedict_R2.update({i:open("R2_"+str(i)+".fq","a")})
    barcounts.update({i:0})
f1= open("unknown_R1.fq","a")
f2= open("unknown_R2.fq","a")
f3= open("hopped_R1.fq","a")
f4= open("hopped_R2.fq","a")


## start a counter for unknow, matched and hopped reads
counter_unknow = 0
counter_matched = 0
counter_hopped = 0

## Read in files
filenames = ['/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz','/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz','/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz','/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz']
with gzip.open(filenames[0],'rt') as R1, gzip.open(filenames[1],'rt') as R2, gzip.open(filenames[2],'rt') as R3, gzip.open(filenames[3],'rt') as R4:
## 
    while True:
        record_R1=[]
        record_R2=[]
        record_R3=[]
        record_R4=[]
        for i in range(4):
            record_R1 = np.append(record_R1,R1.readline().splitlines() )          
            record_R2 = np.append(record_R2,R2.readline().splitlines() )
            record_R3 = np.append(record_R3,R3.readline().splitlines() )
            record_R4 = np.append(record_R4,R4.readline().splitlines() )
        if record_R1[0] =='':
            break
        record_R3[1] = rev_complement(record_R3[1])
## Check if indexes are one of our indexes, otherwise output in unknown
        if record_R2[1] in barcodes and record_R3[1] in barcodes:
            ## quality check
            if cutoff(record_R1[3],record_R2[3],record_R3[3],record_R4[3]):
                record_R1[0] += str(record_R2[1])+ '_' + str(record_R3[1])
                record_R4[0] += str(record_R2[1])+ '_' + str(record_R3[1])
                # f1= open("unknown_R1.fq","a")
                f1.write('\n'.join((record_R1))+'\n')
                # f1.close()
                # f2= open("unknown_R2.fq","a")
                f2.write('\n'.join((record_R4))+'\n')
                # f2.close()
                counter_unknow+=1
            ## match check
            elif check_match(record_R2[1],record_R3[1]):
                record_R1[0] += str(record_R2[1])+ '_' + str(record_R3[1])
                record_R4[0] += str(record_R2[1])+ '_' + str(record_R3[1])
                # f1= open("R1_"+str(record_R2[1])+".fq","a")
                filedict_R1[record_R2[1]].write('\n'.join(record_R1)+'\n')
                # f1.close()
                # f2= open("R2_"+str(record_R3[1])+".fq","a")
                filedict_R2[record_R3[1]].write('\n'.join(record_R4)+'\n')
                # f2.close()
                counter_matched+=1
                barcounts[record_R2[1]]+=1
            ## the rest is hopped reads
            else:
                record_R1[0] += str(record_R2[1])+ '_' + str(record_R3[1])
                record_R4[0] += str(record_R2[1])+ '_' + str(record_R3[1])
                # f1= open("hopped_R1.fq","a")
                f3.write('\n'.join((record_R1))+'\n')
                # f3.close()
                # f2= open("hopped_R2.fq","a")
                f4.write('\n'.join((record_R4))+'\n')
                # f2.close()
                counter_hopped+=1
        else:
            record_R1[0] += str(record_R2[1])+ '_' + str(record_R3[1])
            record_R4[0] += str(record_R2[1])+ '_' + str(record_R3[1])
            # print(record_R1)
            # f1= open("unknown_R1.fq","a")
            f1.write('\n'.join(record_R1)+'\n')
            # f1.close()
            # f2= open("unknown_R2.fq","a")
            f2.write('\n'.join(record_R4)+'\n')
            # f2.close()
            counter_unknow+=1

            
print("The number of unknow reads is " + str(counter_unknow))
print("The number of matched reads is " + str(counter_matched))
print("The number of hopped reads is " + str(counter_hopped))

for i in barcodes:
    print("The percentage of " + str(i) + "is" + str(barcounts[i]/counter_matched))
#Percentage of reads from each sample
#Overall amount of index swapping
#Any figures/any other relevant data your code output
