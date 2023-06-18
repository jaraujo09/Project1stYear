from sys import argv 
import re
import pandas as pd
    
def format(string, bps:int = 10, group:int = 66)->list[str]:
    """
    @brief Function to sort nucleotide sequence in the right format
    @param bps: number of base pair per group
    @param group: number of characters (10 bp x 6 groups + 6 invisible characters)
    @return list of string of the sequence in the right format
    """           
    string1= " ".join([string[i:i+bps] for i in range(0, len(string), bps)])
    return "\n".join([string1[j:j+group]for j in range(0,len(string1), group)])


def lines_to_line(list0:list[str])->str:
    """
    @brief Function to join lines of FASTA file in one line (string)
    @param list0: list of strings form the above function
    @return Sequence of the FASTA file in one line
    """  
    a=''
    for line in list0:
        a += line.strip()
        b = a.split()
    list = ''.join(b)
    return(list)

#Reads in all required file, and create a output file for report.
#command line: python /Users/user/Desktop/UNI/MESTRADO/2º semestre/Projeto/new_test_web.py fasta_seq.fasta enzymes.txt resultado.txt

sequence_file= open(argv[1], "r") #reads in FASTA file with nucleotide sequence
enzyme_file= open(argv[2], "r") #reads in text file with restriction enzymes
output_file= open(argv[3], "w") #write the output to new file (must indicate the type, like .txt)

output_file.write(f"Restriction enzyme analysis of nucleotide sequence from file {argv[1]}.\n Cutting with enzymes from file {argv[2]}.")
output_file.write("\n" + '-' * 80 + "\n")


first_line = sequence_file.readline()  #reads the first line of FASTA file
if first_line.startswith(">"):      #If the first line contains the header ">", meaning its a FASTA file
    sequence_name = first_line.strip('>')   #use as sequence name
    sequence1 = sequence_file.readlines() #takes the rest lines as nucleotide sequence
    sequence = lines_to_line(sequence1) #uses lines_to_line function to combine lines into one string line
    
else:
    sequence_name = "N/A"         #if there is no '>', use N/A as sequence name
    sequence1 = [first_line] + sequence_file.readlines() #take the first and rest lines as whole sequence
    sequence= lines_to_line(sequence1)

output_file.write(f"Sequence name: {sequence_name}\n")
output_file.write(f"Sequence length: {len(sequence)} bp\n")


def get_pattern(s:str)->list[str]:
    """
    @brief Function that gets the pattern of ambiguous nucleotides
    @param s: letter to get that pattern 
    @return list of patterns
    """
    D = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
         'R': '[AG]', 'Y': '[CT]', 'N': '.',
         'S': '[GC]', 'W': '[AT]', 'K': '[GT]',
         'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
         'H': '[ACT]', 'V': '[ACG]',
         '^': ''}
    rL = [D.get(c, c) for c in s]
    return ''.join(rL)


no_cut = []
enzyme_counts = {}
for line in enzyme_file.readlines(): #read in each restriction enzyme, their sequence and their cutting site
    enzyme_name, enzyme_pattern = line.strip().split(";")  #strip by ";" (as in the enzyme_file) to get the enzyme name then the enzyme pattern
    enzyme_pattern = get_pattern(enzyme_pattern)  #apply get_pattern function to get the pattern of ambiguous nucleotides
    print(line)
    before_cutting_site, after_cutting_site = enzyme_pattern.split("|") #the '|' indicates the cleavage site
    
    cutting_sites = []   #list to cut sites


    
    for seq in re.finditer(before_cutting_site + after_cutting_site, sequence):  #find all hits that match the enzyme pattern in the sequence
        match_site = seq.start()      #get the index of start matching position
        cutting_sites.append(match_site+len(before_cutting_site))    #record all cutting sites into a list


    if len(cutting_sites) != 0:  #if there are matches for enzyme in sequence, print the count of cutting sites and sequence fragments
        output_file.write('-' * 80 + f"\nThere are {len(cutting_sites)} cutting sites for {enzyme_name}, cutting at {enzyme_pattern}\n")
        output_file.write(f"There are {len(cutting_sites)+1} fragments:\n\n")

        #length of each fragment and index info
        cutting_sites.append(len(sequence)) 
        start_cutting_site = 0
        for site in cutting_sites:
            fragment = sequence[start_cutting_site:site]
            output_file.write(f"Length: {len(fragment)} range: {start_cutting_site+1}-{site}\n")
            output_file.write(format(fragment) + "\n") #use formating function to put fragment sequence into the right format
            start_cutting_site = site
        enzyme_counts[enzyme_name] = len(cutting_sites)

    if len(cutting_sites) == 0:   #if no cutting site found, print the result at the end
        no_cut.append(enzyme_name)
output_file.write('-' * 80 + f"\nThere are no cutting sites for {no_cut} in this nucleotide sequence.\n")

# Visualization
most_cut_enzyme = max(sorted(enzyme_counts.values(), reverse = True))
least_cut_enzyme = min(sorted(enzyme_counts.values(), reverse = True))
order = sorted(enzyme_counts.values(), reverse = True)
print("Sorted", order)
most = [key for key, val in enzyme_counts.items() if val == most_cut_enzyme]
least = [key for key, val in enzyme_counts.items() if val == least_cut_enzyme]
print("Most Cut Enzyme:" , most)
print("Least Cut Enzyme: ", least)
print("No cuts: ", no_cut)

less_ten = []
ten_to_hundred = []
more_hundred = []
val = map(int, enzyme_counts.values())
print(type(val))

for enzyme, value in enzyme_counts.items():
    if 1 <= value < 10:
        less_ten.append(enzyme)
    elif value > 100:
        more_hundred.append(enzyme)
    else:
        ten_to_hundred.append(enzyme)

print("Less than 10:", less_ten)
print("Between 10 and 100:", ten_to_hundred)
print("More than 100:", more_hundred)

a = len(ten_to_hundred)
b = len(less_ten)
c = len(more_hundred)
d = a+b+c
print(d)

def create_result_file(filename)->None:
    """
    @brief Function to create a result file
    @param filename: name of the file to create
    """
    less_ten_count = len(less_ten)
    ten_to_hundred_count = len(ten_to_hundred)
    more_hundred_count = len(more_hundred)
    no_cut_count = len(no_cut)

    with open(filename, 'w') as file:
        file.write("Counting the number of fragments generated by each enzyme...\n")
        file.write("We can verify that, in the total of 404 restriction enzymes:\n")
        file.write("{} enzymes cut 1-10 times\n".format(less_ten_count))
        file.write("{} enzymes cut 10 to 100 times\n".format(ten_to_hundred_count))
        file.write("{} enzymes cut more than 100 times\n".format(more_hundred_count))
        file.write("{} enzymes have no cuts in this phage\n".format(no_cut_count))

create_result_file(filename="new_counting_T66.txt")


  
sequence_file.close()  #Close all files
enzyme_file.close()
output_file.close()
