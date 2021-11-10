"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename)
    text = f.read()
    x = ""
    for i in text.splitlines():
        x += i
   
    return x
       

    


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    a = [] 
    for i in range(startIndex,len(dna),3):
        a.append(dna[i:i+3])
        if dna[i:i+3] == "TAG" or dna[i:i+3] == "TAA" or dna[i:i+3] == "TGA":
            break
    b = [ ]
    for i in a:
        x = i.replace("T","U")
        b.append(x)
    return b
    
    
    
    
'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
import json
def makeCodonDictionary(filename):
    x = { }
    f = open(filename)
    j = json.load(f)
    for a,b in j.items():
        for i in b:
            y = i.replace("T","U")

            x[y] = a
    return x


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    a = []
    if codons[0] == "AUG":
        a.append("Start")
    for i in range(1,len(codons)):
        if codons[i] in codonD.keys():
            a.append(codonD[codons[i]])
    
    return a
    

        
    
       
        



   


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    x = readFile(dnaFilename)
    y = makeCodonDictionary(codonFilename)

    i = 0
    a = 0
    b = [ ]
    while i < len(x):
        if x[i:i+3] == "ATG":
            c = dnaToRna(x,i)
            
            d = generateProtein(c,y)
           
            b.append(d)
            i = i+3*len(c)
        else:
            i += 1
            a +=1

    return b


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    x = []

    for i in proteinList1 :
        for j in proteinList2:
            if i == j and i not in x:
                x.append(j)

    return x


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    x = []
    for i in proteinList:
        for j in i:
                x.append(j)
    return x

   



'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    x = { }
    for i in aaList:
        if i not in x:
            x[i] = 1
        else:
            x[i] += 1
    

    return x


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    combine1,combine2 = combineProteins(proteinList1), combineProteins(proteinList2)
    dict1,dict2 = aminoAcidDictionary(combine1),aminoAcidDictionary(combine2)
    temp,result=[],[]         
    freq_dict1,freq_dict2={},{}    
    for i in dict1:
        freq_dict1[i] = dict1[i]/len(combine1)
        if i not in temp and i !="Start" and i !="Stop":
            temp.append(i)
    for a in dict2:
        freq_dict2[a] = dict2[a]/len(combine2)
        if a not in temp and a !="Start" and a!="Stop":
            temp.append(a)
    for ac in temp:
        freq1,freq2=0,0
        if ac in freq_dict1:
            freq1= freq_dict1[ac]
        if ac in freq_dict2:
            freq2= freq_dict2[ac]
        difference = freq2-freq1
        if difference < -cutoff or difference > cutoff  :
            result.append([ac , freq1, freq2])
    return result


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    for i in commonalities:
        commonProteins = ""
        x = i[1:len(i)-1]
        count = 0
        for j in x:
            commonProteins+=j
            count+=1
            if count !=len(x):
                commonProteins+="-"
    for a in differences:
        print(a[0],":", round(a[1]*100,2),"% in seq1,",round(a[2]*100,2),"% in seq2")
    



    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    x = []
    for i in proteinList1:
        for j in i:
            if j not in x:
                x.append(j)
    for a in proteinList2:
        for y in a:
            if y not in x:
                x.append(y)
                
    return sorted(x)
    


    


    '''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    

    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testCombineProteins()

    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # test.runWeek2()
    test.testMakeAminoAcidLabels()
    

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
