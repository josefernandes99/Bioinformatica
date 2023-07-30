def readDNA(f):
    DNA = [] 
    with open(f, 'r') as file:
        lines = file.readlines()[1:]
        for line in lines:
            for letter in line:  
                if((letter != '\n') and (len(DNA) < 30000)):
                    DNA.append(letter)
        return DNA
    
def calculateFrequency(DNA):
    numA = 0
    numT = 0
    numG = 0
    numC = 0
    
    for i in range(len(DNA)):
        if (DNA[i] == 'A'):
            numA += 1
        elif (DNA[i] == 'T'):
            numT += 1 
        elif (DNA[i] == 'G'):
            numG += 1
        elif (DNA[i] == 'C'):
            numC += 1
    numA = round((numA * 100)/30000,2)
    numT = round((numT * 100)/30000,2)
    numG = round((numG * 100)/30000,2)
    numC = round((numC * 100)/30000,2)
    
    print("2. Frequency (in %): A = {}, C = {}, G = {}, T = {}".format(numA, numC, numG, numT))
    print("3. GC Content = {}".format(numG + numC))
    
def searchCodons(DNA, DNA_reversed):
    startCodon1 = 0
    codonTAG1 = 0
    codonTAA1 = 0
    codonTGA1 = 0
    
    for i in range(len(DNA)-2):
        if(DNA[i] == "A"):
            if(DNA[i+1] == "T"):  
                if(DNA[i+2] == "G"):
                    startCodon1 += 1    
        if(DNA[i] == "T"):
            if(DNA[i+1] == "A"):
                if(DNA[i+2] == "G"):
                    codonTAG1 += 1 
                elif(DNA[i+2] == "A"):
                    codonTAA1 += 1
            elif(DNA[i+1] == "G"):
                if(DNA[i+2] == "A"):
                    codonTGA1 += 1
                    
    startCodon2 = 0
    codonTAG2 = 0
    codonTAA2 = 0
    codonTGA2 = 0
    
    for i in range(len(DNA_reversed)-2):
        if(DNA_reversed[i] == "A"):
            if(DNA_reversed[i+1] == "T"):
                if(DNA_reversed[i+2] == "G"):
                    startCodon2 += 1 
        if(DNA_reversed[i] == "T"):
            if(DNA_reversed[i+1] == "A"):
                if(DNA_reversed[i+2] == "G"):
                    codonTAG2 += 1
                elif(DNA_reversed[i+2] == "A"):
                    codonTAA2 += 1
            elif(DNA_reversed[i+1] == "G"):
                
                if(DNA_reversed[i+2] == "A"):
                    codonTGA2 += 1
                    
    temp1 = [[startCodon1, codonTAG1, codonTAA1, codonTGA1],["ATG","TAG","TAA","TGA"]]
    temp2 = [[startCodon2, codonTAG2, codonTAA2, codonTGA2],["ATG","TAG","TAA","TGA"]]
                
    print("4. Number of Start Codons: PositiveStrand = {}, NegativeStrand = {}".format(startCodon1, startCodon2))
    print("5. Number of Stop Codons: PositiveStrand = {}, NegativeStrand = {}".format((codonTAG1 + codonTAA1 + codonTGA1), (codonTAG2 + codonTAA2 + codonTGA2)))
    print("6. Most Frequent (+) = {} ({}), Least Frequent (+) = {} ({}), Most Frequent (-) = {} ({}), Least Frequent (-) = {} ({})".format(temp1[1][temp1[0].index(max(temp1[0]))], max(temp1[0]), temp1[1][temp1[0].index(min(temp1[0]))], min(temp1[0]), temp2[1][temp2[0].index(max(temp2[0]))], max(temp2[0]), temp2[1][temp2[0].index(min(temp2[0]))], min(temp2[0]))) 

def getORFs(DNA, DNA_reversed):
    ORFs = []
    inORF = 0
    tempString1 = ""
    tempString2 = ""
    startCoord = 0
    stopCoord = 0
    countORF = 1
    
    outF1 = open("all_potential_proteins.txt", "w")
    outF2 = open("orf_coordinates.txt", "w")
    
    for i in range(0, len(DNA)-2, 3):
        if(DNA[i] == 'A' and DNA[i+1] == 'T' and DNA[i+2] == 'G'): #found ATG
            inORF = 1
            startCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'G' and DNA[i+2] == 'A'): #found TGA
            inORF = 0
            stopCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'A' and DNA[i+2] == 'A'): #found TAA
            inORF = 0
            stopCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'A' and DNA[i+2] == 'G'): #found TAG
            inORF = 0
            stopCoord = i
        if(inORF == 1):
            tempString1 += DNA[i] + DNA[i+1] + DNA[i+2]
        elif(inORF == 0 and tempString1 != ""):
            if(len(tempString1) >= 150):
                tempString1 += "\n"
                tempString2 += str(startCoord) + ", " + str(stopCoord) + ", ORF" + str(countORF) + "\n"
                outF1.write(tempString1)
                ORFs.append(tempString1)
                outF2.write(tempString2)
                countORF += 1
            tempString1 = ""
            tempString2 = ""
            startCoord = 0
            stopCoord = 0
    
    tempString1 = ""
    tempString2 = ""
    startCoord = 0
    stopCoord = 0
    
    for i in range(1, len(DNA)-2, 3):
        if(DNA[i] == 'A' and DNA[i+1] == 'T' and DNA[i+2] == 'G'): #found ATG
            inORF = 1
            startCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'G' and DNA[i+2] == 'A'): #found TGA
            inORF = 0
            stopCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'A' and DNA[i+2] == 'A'): #found TAA
            inORF = 0
            stopCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'A' and DNA[i+2] == 'G'): #found TAG
            inORF = 0
            stopCoord = i
        if(inORF == 1):
            tempString1 += DNA[i] + DNA[i+1] + DNA[i+2]
        elif(inORF == 0 and tempString1 != ""):
            if(len(tempString1) >= 150):
                tempString1 += "\n"
                tempString2 += str(startCoord) + ", " + str(stopCoord) + ", ORF" + str(countORF) + "\n"
                outF1.write(tempString1)
                ORFs.append(tempString1)
                outF2.write(tempString2)
                countORF += 1
            tempString1 = ""
            tempString2 = ""
            startCoord = 0
            stopCoord = 0
            
    tempString1 = ""
    tempString2 = ""
    startCoord = 0
    stopCoord = 0
    
    for i in range(2, len(DNA)-2, 3):
        if(DNA[i] == 'A' and DNA[i+1] == 'T' and DNA[i+2] == 'G'): #found ATG
            inORF = 1
            startCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'G' and DNA[i+2] == 'A'): #found TGA
            inORF = 0
            stopCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'A' and DNA[i+2] == 'A'): #found TAA
            inORF = 0
            stopCoord = i
        elif(DNA[i] == 'T' and DNA[i+1] == 'A' and DNA[i+2] == 'G'): #found TAG
            inORF = 0
            stopCoord = i
        if(inORF == 1):
            tempString1 += DNA[i] + DNA[i+1] + DNA[i+2]
        elif(inORF == 0 and tempString1 != ""):
            if(len(tempString1) >= 150):
                tempString1 += "\n"
                tempString2 += str(startCoord) + ", " + str(stopCoord) + ", ORF" + str(countORF) + "\n"
                outF1.write(tempString1)
                ORFs.append(tempString1)
                outF2.write(tempString2)
                countORF += 1
            tempString1 = ""
            tempString2 = ""
            startCoord = 0
            stopCoord = 0
    
    tempString1 = ""
    tempString2 = ""
    startCoord = 0
    stopCoord = 0
            
    for i in range(0, len(DNA_reversed)-2, 3):
        if(DNA_reversed[i] == 'A' and DNA_reversed[i+1] == 'T' and DNA_reversed[i+2] == 'G'): #found ATG
            inORF = 1
            startCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'G' and DNA_reversed[i+2] == 'A'): #found TGA
            inORF = 0
            stopCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'A' and DNA_reversed[i+2] == 'A'): #found TAA
            inORF = 0
            stopCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'A' and DNA_reversed[i+2] == 'G'): #found TAG
            inORF = 0
            stopCoord = i
        if(inORF == 1):
            tempString1 += DNA_reversed[i] + DNA_reversed[i+1] + DNA_reversed[i+2]
        elif(inORF == 0 and tempString1 != ""):
            if(len(tempString1) >= 150):
                tempString1 += "\n"
                tempString2 += str(startCoord) + ", " + str(stopCoord) + ", ORF" + str(countORF) + "\n"
                outF1.write(tempString1)
                ORFs.append(tempString1)
                outF2.write(tempString2)
                countORF += 1
            tempString1 = ""
            tempString2 = ""
            startCoord = 0
            stopCoord = 0
    
    tempString1 = ""
    tempString2 = ""
    startCoord = 0
    stopCoord = 0
    
    for i in range(1, len(DNA_reversed)-2, 3):
        if(DNA_reversed[i] == 'A' and DNA_reversed[i+1] == 'T' and DNA_reversed[i+2] == 'G'): #found ATG
            inORF = 1
            startCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'G' and DNA_reversed[i+2] == 'A'): #found TGA
            inORF = 0
            stopCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'A' and DNA_reversed[i+2] == 'A'): #found TAA
            inORF = 0
            stopCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'A' and DNA_reversed[i+2] == 'G'): #found TAG
            inORF = 0
            stopCoord = i
        if(inORF == 1):
            tempString1 += DNA_reversed[i] + DNA_reversed[i+1] + DNA_reversed[i+2]
        elif(inORF == 0 and tempString1 != ""):
            if(len(tempString1) >= 150):
                tempString1 += "\n"
                tempString2 += str(startCoord) + ", " + str(stopCoord) + ", ORF" + str(countORF) + "\n"
                outF1.write(tempString1)
                ORFs.append(tempString1)
                outF2.write(tempString2)
                countORF += 1
            tempString1 = ""
            tempString2 = ""
            startCoord = 0
            stopCoord = 0
            
    tempString1 = ""
    tempString2 = ""
    startCoord = 0
    stopCoord = 0
    
    for i in range(2, len(DNA_reversed)-2, 3):
        if(DNA_reversed[i] == 'A' and DNA_reversed[i+1] == 'T' and DNA_reversed[i+2] == 'G'): #found ATG
            inORF = 1
            startCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'G' and DNA_reversed[i+2] == 'A'): #found TGA
            inORF = 0
            stopCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'A' and DNA_reversed[i+2] == 'A'): #found TAA
            inORF = 0
            stopCoord = i
        elif(DNA_reversed[i] == 'T' and DNA_reversed[i+1] == 'A' and DNA_reversed[i+2] == 'G'): #found TAG
            inORF = 0
            stopCoord = i
        if(inORF == 1):
            tempString1 += DNA_reversed[i] + DNA_reversed[i+1] + DNA_reversed[i+2]
        elif(inORF == 0 and tempString1 != ""):
            if(len(tempString1) >= 150):
                tempString1 += "\n"
                tempString2 += str(startCoord) + ", " + str(stopCoord) + ", ORF" + str(countORF) + "\n"
                outF1.write(tempString1)
                ORFs.append(tempString1)
                outF2.write(tempString2)
                countORF += 1
            tempString1 = ""
            tempString2 = ""
            startCoord = 0
            stopCoord = 0
            
    return ORFs

def readGenesFile(f):
    genes = []
    with open(f, "r") as arquivo:
        for linha in arquivo:
            if "exon" in linha:
                palavra = linha.split()
                gene = [palavra[3],palavra[4],palavra[9].replace(';', '').replace('"', '')]
                genes.append(gene)
    return genes

def calculateOverlap(A1, A2, B1, B2):
    intersection = min(A2, B2) - max(A1, B1)
    if intersection <= 0:
        return 0
    return (intersection / (A2-A1))

def OverlapWithAnnotation(f, genes): 
    orfs = []
    with open(f, "r") as arquivo:
        for linha in arquivo:
            palavras = linha.split()
            orf = [palavras[0].replace(',', ''),palavras[1].replace(',', ''),palavras[2]]
            orfs.append(orf)
               
    for i in range(len(genes)):
        maxIntersection = 0
        for j in orfs:
            intersection = calculateOverlap(int(genes[i][0]), int(genes[i][1]), int(j[0]), int(j[1]))
            if maxIntersection < intersection:
                maxIntersection = intersection
        print(genes[i][2], " ", maxIntersection*100, "%")

DNA = readDNA('sequence_chr1.fasta')
DNA_reversed = list(reversed(DNA))

print("1. Length of the sequence:", len(DNA))
calculateFrequency(DNA)
searchCodons(DNA, DNA_reversed)
ORFs = getORFs(DNA, DNA_reversed)

genes = readGenesFile('genes_chr1.gtf') #lê o ficheiro dos ORF fornecido e coloca num array genes, só com as coordenadas
outF1 = open("genes.txt", "w")
outF1.write(str(genes))
OverlapWithAnnotation('orf_coordinates.txt', genes)
                     
""" def separateGenesFile(f):
    outF1 = open("positive_strand.txt", "w")
    outF2 = open("negative_strand.txt", "w")
    with open(f, 'r') as file:
        lines = file.readlines()
        for line in lines:
            curLine = line.split()
            if(curLine[6] == "+"):
                outF1.write(line)
            elif(curLine[6] == "-"):
                outF2.write(line)  
                
def getStrands(f, DNA):
    startCoord = 0
    endCoord = 0
    strands = []
    pos = 0
    with open(f, 'r') as file:
        lines = file.readlines()
        for line in lines:
            curStrand = ""
            curLine = line.split()
            if(curLine[2] == "CDS"):
                startCoord = int(curLine[3])
                endCoord = int(curLine[4])
                for i in range(startCoord, endCoord+1):
                    curStrand += str(DNA[i])
                strands.append([curStrand])
                pos += 1
    return strands """
    
""" def convert(s):
 
    # initialization of string to ""
    new = ""
 
    # traverse in the string
    for x in s:
        new += x
 
    # return string
    return new """
    


#print(convert(DNA))
""" separateGenesFile('genes_chr1.gtf')
positiveStrands = getStrands('positive_strand.txt', DNA)
negativeStrands = getStrands('negative_strand.txt', DNA)
print(positiveStrands)
 """