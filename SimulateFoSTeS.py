## SimulateFoSTeS v1.0
## Copyright 2020, Samuel Tremblay-Belzile

import random

#FASTA file of reference genome from which you want to simulate FoSTeS
refgenome = "genome.fasta"
#Base name of the FastQ files
fastqfile = "FoSTeSSim"
#File that will contain the detailed positions of the FoSTeS events
FoSTeSfile = "Junctions" + fastqfile + ".txt"
#Files that will contain the positions of all the reads generated
reffileR1 = fastqfile + "RefR1.txt"
reffileR2 = fastqfile + "RefR2.txt"
#Number of reads to simulate
reads = 50000000
#Length of each read
readlen = 150
#Minimum and maximum values of the DNA fragment between R1 and R2 reads
mininsert = 200
maxinsert = 300
#Rate at which FoSTeS events will occur
tsrate = 0.001
#Ratio of FoSTeS events that will switch to another position on the same chromosome
samechrrate = 0.5

#Make a list of chromosome sequences in the reference genome

f = open(refgenome, 'r')
totalseqs = 0
chromosomes = []
seqname = ""
seqlen = 0
seq = ""
totalbases = 0

for line in f:
    line = line.strip('\n')
    if line=="":
        pass
    elif line[0] == ">":
        if seqname != "":
            chromosomes.append([seqname,seqlen,seq])
        totalseqs += 1
        totalbases += seqlen
        seqname = line[1:]
        seqlen = 0
        seq = ""
    else:
        seqlen += len(line)
        seq += line

chromosomes.append([seqname,seqlen,seq])
totalbases += seqlen

f.close()

g1 = open(fastqfile + "R1.fastq", "w")
g2 = open(fastqfile + "R2.fastq", "w")
h = open(FoSTeSfile, "w")
h2 = open(reffileR1, "w")
h3 = open(reffileR2, "w")

for i in range(reads):
    
    #Determine if there is a template switch in R1 or R2
    switchR1 = False
    switchR2 = False
    switch = random.randrange(0,1000000,1)/1000000.0
    if switch <= tsrate:
        if random.randrange(0,2,1)==0:
            switchR1 = True
        else:
            switchR2 = True

    #Find the starting position for sequencing
    while True:
        #Determine initial position of sequencing
        posstart = random.randrange(1,totalbases, 1)
        chrpos = posstart
        chrstart = 0
        while chrpos >= chromosomes[chrstart][1]:
            chrpos -= chromosomes[chrstart][1]
            chrstart += 1
            
        #Determine orientation of sequencing
        direction = random.randrange(0,2,1) #0 is forward, 1 is reverse

        #Look for new position until beginning and end bases are not N
        #and until the reads will not go beyond the chromosome end
        if direction == 0:
            if ((chromosomes[chrstart][2][chrpos] != 'N') and
                (chromosomes[chrstart][1] > (chrpos+(2*readlen)-1+maxinsert))):
                if (chromosomes[chrstart][2][chrpos+(2*readlen)-1+maxinsert] != 'N'):
                    break
        else:
            if ((chromosomes[chrstart][2][chrpos] != 'N') and
                ((chrpos-(2*readlen)-maxinsert+1) >= 0)):
                if (chromosomes[chrstart][2][chrpos-(2*readlen)-maxinsert+1] != 'N'):
                    break

    #Write the R1 file
    g1.write("@" + fastqfile + "." + str(i+1) + "\n")
    h2.write("@" + fastqfile + "." + str(i+1) + "\t")
    h2.write(chromosomes[chrstart][0] + "\t" + str(chrpos) + "\t" + str(direction) + "\n")
    if switchR1 == False:
        if direction == 0:
            for j in range(readlen):
                g1.write(chromosomes[chrstart][2][chrpos+j])
        if direction == 1:
            for j in range(readlen):
                letter = chromosomes[chrstart][2][chrpos-j]
                if letter == 'A' or letter == 'a':
                    g1.write('T')
                elif letter == 'T' or letter == 't':
                    g1.write('A')
                elif letter == 'C' or letter == 'c':
                    g1.write('G')
                elif letter == 'G' or letter == 'g':
                    g1.write('C')
                else:
                    g1.write(letter)
    else:
        #Determine how far into the read the switch occurs
        switchpoint = random.randrange(1,readlen,1)

        #Write the beginning of the sequence   
        if direction == 0:
            for j in range(switchpoint):
                g1.write(chromosomes[chrstart][2][chrpos+j])
            h.write(fastqfile + "." + str(i+1) + "1\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write("1\t" + str(switchpoint) + "\t")
            h.write(str(chrpos+1) + "\t" + str(chrpos+switchpoint) + "\n")
        if direction == 1:
            for j in range(switchpoint):
                letter = chromosomes[chrstart][2][chrpos-j]
                if letter == 'A' or letter == 'a':
                    g1.write('T')
                elif letter == 'T' or letter == 't':
                    g1.write('A')
                elif letter == 'C' or letter == 'c':
                    g1.write('G')
                elif letter == 'G' or letter == 'g':
                    g1.write('C')
                else:
                    g1.write(letter)
            h.write(fastqfile + "." + str(i+1) + "1\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write("1\t" + str(switchpoint) + "\t")
            h.write(str(chrpos+1) + "\t" + str(chrpos-switchpoint+2) + "\n")

        #Determine whether the switch is on same chromosome
        if random.random() >= samechrrate:
            samechr = False
        else:
            samechr = True

        #Determine the new genome position
        while True:
            if samechr:
                chrpos = random.randrange(0,chromosomes[chrstart][1],1)
                #Determine orientation of sequencing
                direction = random.randrange(0,2,1) #0 is forward, 1 is reverse

                #Look for new position until beginning and end bases are not N
                #and until the reads will not go beyond the chromosome end
                if direction == 0:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        (chromosomes[chrstart][1] > (chrpos+(2*readlen)-1+maxinsert))):
                        if (chromosomes[chrstart][2][chrpos+(2*readlen)-1+maxinsert] != 'N'):
                            break
                else:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        ((chrstart-(2*readlen)-maxinsert+1) >= 0)):
                        if (chromosomes[chrstart][2][chrstart-(2*readlen)-maxinsert+1] != 'N'):
                            break
            else:
                #Determine initial position of sequencing
                posstart = random.randrange(1,totalbases, 1)
                chrpos = posstart
                chrstart = 0
                while chrpos >= chromosomes[chrstart][1]:
                    chrpos -= chromosomes[chrstart][1]
                    chrstart += 1
                    
                #Determine orientation of sequencing
                direction = random.randrange(0,2,1) #0 is forward, 1 is reverse

                #Look for new position until beginning and end bases are not N
                #and until the reads will not go beyond the chromosome end
                if direction == 0:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        (chromosomes[chrstart][1] > (chrpos+(2*readlen)-1+maxinsert))):
                        if (chromosomes[chrstart][2][chrpos+(2*readlen)-1+maxinsert] != 'N'):
                            break
                else:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        ((chrstart-(2*readlen)-maxinsert+1) >= 0)):
                        if (chromosomes[chrstart][2][chrstart-(2*readlen)-maxinsert+1] != 'N'):
                            break
        #Write the rest of the sequence
        if direction == 0:
            for j in range(readlen-switchpoint):
                g1.write(chromosomes[chrstart][2][chrpos+j])
            h.write(fastqfile + "." + str(i+1) + "1\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write(str(switchpoint+1) + "\t" + str(readlen) + "\t")
            h.write(str(chrpos+1)+"\t"+str(chrpos+readlen-switchpoint) + "\n")
        if direction == 1:
            for j in range(readlen-switchpoint):
                letter = chromosomes[chrstart][2][chrpos-j]
                if letter == 'A' or letter == 'a':
                    g1.write('T')
                elif letter == 'T' or letter == 't':
                    g1.write('A')
                elif letter == 'C' or letter == 'c':
                    g1.write('G')
                elif letter == 'G' or letter == 'g':
                    g1.write('C')
                else:
                    g1.write(letter)
            h.write(fastqfile + "." + str(i+1) + "1\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write(str(switchpoint+1) + "\t" + str(readlen) + "\t")
            h.write(str(chrpos+1)+"\t"+str(chrpos-readlen+switchpoint+2) + "\n")

    #Write the rest of the FastQ info
    g1.write("\n+\n")
    for k in range(readlen):
        g1.write("I")
    g1.write("\n")

    #Write the R2 file
    g2.write("@" + fastqfile + "." + str(i+1) + "\n")
    h3.write("@" + fastqfile + "." + str(i+1) + "\t")
    h3.write(chromosomes[chrstart][0] + "\t" + str(chrpos) + "\t" + str(direction) + "\n")
    if switchR2 == False:
        if direction == 1:
            r2pos = (chrpos - (2*readlen) + 1 -
                     random.randrange(mininsert, maxinsert+1,1))
            for j in range(readlen):
                g2.write(chromosomes[chrstart][2][r2pos+j])
        else:
            r2pos = (chrpos + (2*readlen) - 1 +
                     random.randrange(mininsert, maxinsert+1,1))
            for j in range(readlen):
                letter = chromosomes[chrstart][2][r2pos-j]
                if letter == 'A' or letter == 'a':
                    g2.write('T')
                elif letter == 'T' or letter == 't':
                    g2.write('A')
                elif letter == 'C' or letter == 'c':
                    g2.write('G')
                elif letter == 'G' or letter == 'g':
                    g2.write('C')
                else:
                    g2.write(letter)
    else:
        #Determine how far into the read the switch occurs
        switchpoint = random.randrange(1,readlen,1)

        #Write the beginning of the sequence   
        if direction == 1:
            r2pos = (chrpos - (2*readlen) + 1 -
                     random.randrange(mininsert, maxinsert+1,1))
            for j in range(switchpoint):
                g2.write(chromosomes[chrstart][2][r2pos+j])
            h.write(fastqfile + "." + str(i+1) + "2\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write("1\t" + str(switchpoint) + "\t")
            h.write(str(r2pos+1) + "\t" + str(r2pos+switchpoint) + "\n")
        else:
            r2pos = (chrpos + (2*readlen) - 1 +
                     random.randrange(mininsert, maxinsert+1,1))
            for j in range(switchpoint):
                letter = chromosomes[chrstart][2][r2pos-j]
                if letter == 'A' or letter == 'a':
                    g2.write('T')
                elif letter == 'T' or letter == 't':
                    g2.write('A')
                elif letter == 'C' or letter == 'c':
                    g2.write('G')
                elif letter == 'G' or letter == 'g':
                    g2.write('C')
                else:
                    g2.write(letter)
            h.write(fastqfile + "." + str(i+1) + "2\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write("1\t" + str(switchpoint) + "\t")
            h.write(str(r2pos+1) + "\t" + str(r2pos-switchpoint+2) + "\n")

        #Determine whether the switch is on same chromosome
        if random.random() >= samechrrate:
            samechr = False
        else:
            samechr = True

        #Determine the new genome position
        while True:
            if samechr:
                chrpos = random.randrange(0,chromosomes[chrstart][1],1)
                #Determine orientation of sequencing
                direction = random.randrange(0,2,1) #0 is forward, 1 is reverse

                #Look for new position until beginning and end bases are not N
                #and until the reads will not go beyond the chromosome end
                if direction == 0:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        (chromosomes[chrstart][1] > (chrpos+readlen-switchpoint))):
                        if (chromosomes[chrstart][2][chrpos+readlen-switchpoint] != 'N'):
                            break
                else:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        ((chrstart-readlen+switchpoint) >= 0)):
                        if (chromosomes[chrstart][2][chrstart-readlen+switchpoint] != 'N'):
                            break
            else:
                #Determine initial position of sequencing
                posstart = random.randrange(1,totalbases, 1)
                chrpos = posstart
                chrstart = 0
                while chrpos >= chromosomes[chrstart][1]:
                    chrpos -= chromosomes[chrstart][1]
                    chrstart += 1
                    
                #Determine orientation of sequencing
                direction = random.randrange(0,2,1) #0 is forward, 1 is reverse

                #Look for new position until beginning and end bases are not N
                #and until the reads will not go beyond the chromosome end
                if direction == 0:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        (chromosomes[chrstart][1] > (chrpos+readlen-switchpoint))):
                        if (chromosomes[chrstart][2][chrpos+readlen-switchpoint] != 'N'):
                            break
                else:
                    if ((chromosomes[chrstart][2][chrpos] != 'N') and
                        ((chrstart-readlen+switchpoint) >= 0)):
                        if (chromosomes[chrstart][2][chrstart-readlen+switchpoint] != 'N'):
                            break
        #Write the rest of the sequence
        if direction == 0:
            for j in range(readlen-switchpoint):
                g2.write(chromosomes[chrstart][2][chrpos+j])
            h.write(fastqfile + "." + str(i+1) + "2\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write(str(switchpoint+1) + "\t" + str(readlen) + "\t")
            h.write(str(chrpos+1)+"\t"+str(chrpos+readlen-switchpoint) + "\n")
        if direction == 1:
            for j in range(readlen-switchpoint):
                letter = chromosomes[chrstart][2][chrpos-j]
                if letter == 'A' or letter == 'a':
                    g2.write('T')
                elif letter == 'T' or letter == 't':
                    g2.write('A')
                elif letter == 'C' or letter == 'c':
                    g2.write('G')
                elif letter == 'G' or letter == 'g':
                    g2.write('C')
                else:
                    g2.write(letter)
            h.write(fastqfile + "." + str(i+1) + "2\t")
            h.write(chromosomes[chrstart][0] + "\t")
            h.write(str(switchpoint+1) + "\t" + str(readlen) + "\t")
            h.write(str(chrpos+1)+"\t"+str(chrpos-readlen+switchpoint+2) + "\n")
                        
    #Write the rest of the FastQ info
    g2.write("\n+\n")
    for k in range(readlen):
        g2.write("I")
    g2.write("\n")

g1.close()              
g2.close()


    
