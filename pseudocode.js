Part 1:
1. See "part1_1"
2. a) See inside the repo
   b) We should choose 10 as the cutoff. 
	A phred score of 10 means a 90% accuracy and 20 means a 99% accuracy.
	Therefore, we choose 10 to preseve most data and ensure at least 90% avarage.
   c) See screenshot. I first created symlink for R2 and R3 files. 

Part 2:
1. Define the problem:
	We indentify the index hopping and cut off those with unknown barcode or low quality. 
	Then we bin the duel-mathced records into their each individual file based on the barcodes.
	For our example, we will have 48 different files for each matched barcode.

2. Determine/describe what output would be informative:
	A heat map with each indexes as x and y label would be informative.
	Hopefully, the heat map will show the highest at diagnol and smaller everwhere else.

3. High level funtions:
def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter)-33

def rev_compliment(R3):
	"""Convert A,C,T,G into T,G,A,C then reverse the order"""
	Change A to T
	Change C to G
	Change T to A
	Change G to C
	Lastly, reveser it //can use string.reverse()
	return reverse compliment of the barcodes in R3

4. Pseudo code:
Read all four files line by line
record_R1(_R2,_R3,_R4) = store every four lines into array
	rev_R3 = use function rev_compliment on record_R3[1]
	so now record_R3 should have the reverse compliment read.
	if record_R2[1] in barcode file & record_R3[1] in barcode file
		if <cutoff
			output to unknow file 
			counter_unknow +=1
		elif record_R2[1] == record_R3[1]
			//record_R1[0] += str(record_R2[1])_str(record_R3[1])
			//record_R4[0] += str(record_R2[1])_str(record_R3[1])
			output record_R1 to file (str(record_R2[1])_R1.fq) and record_R3[1] to file (str(rev_R3))
			counter_matched +=1
		else
			output to hopped file 
			counter_hopped +=1
	else
		Add the barcodes in the end of the heade line
		// record_R1[0] += str(record_R2[1])_str(record_R3[1]) 
		// record_R4[0] += str(record_R2[1])_str(record_R3[1])
		record_R1 -> file(unknown_R1.fq) 
		record_R4 -> file (unknown_R2.fq)
		counter_unknow +=1


