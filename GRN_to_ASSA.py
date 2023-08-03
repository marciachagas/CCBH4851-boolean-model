
# Reads the GRN file (exactly as published CCBH-2022.csv) and returns the input file for 
# ASSA-PBN software https://satoss.uni.lu/software/ASSA-PBN/#quick for further calculation 
# of GRN atractors in dynamics. 
# 
# Parameter "perturbation" is set to 0.000. 
# 
# Nodes with at least an output edge are considered in the sub-network
# Boolean rules are written for those nodes having at least an input edge
#
# Use this script parsing the .csv file containing the GRN and file name .pbn for writing in ASSA-PBN input format
#
# Also generates the file conserved_genes.csv containing the gene names of regulators that are not targets.
# 
# (Runs in Python 3.8.10).

# by Marcelo (MTS), 18/09/2022
import sys
import csv
import re

def main():
    #error handling
    if len(sys.argv) != 3:
        sys.exit("Use this script parsing the .csv file containing the GRN and file name .pbn for writing in ASSA-PBN input format")
    elif not sys.argv[1].endswith(".csv"):
        sys.exit("Not a csv file!")
    elif not sys.argv[2].endswith(".pbn"):
        sys.exit("Not a pbn file!")

    # this function reads the grn file and returns grn data structure:
    grn = read_grn(sys.argv[1])
    print("Number of nodes in atractors sub-network = ", len(grn[0]))
    # grn[0] #list of nodes in sub-network
    # grn[1] #dict: key = node;  value = list of  positive regulators of node in key
    # grn[2] #dict: key = node;  value = list of negative regulators of node in key
    # grn[3] #dict: key = node;  value = list  unknown mode of regulators of node in key 
    
    # this function prints ASSA-PBN input file:
    print_assa_format(grn, sys.argv[2])
    print("generated file: ", sys.argv[2])
    print( "generated file: ", "conserved_genes.csv")


 
    
def read_grn(file_name):
    try:
        with open(file_name, "r") as f_name:
            reader = csv.DictReader(f_name, delimiter = ';')
            tfs = []
            targets = []

            for row in reader:
                #remove space, and parenthesis since ASSA-PBN does not admit them. Keep only gene names as ID
                row['Regulatory gene'] = re.sub(" \(.+\)", "", row['Regulatory gene'])
                row['Target gene'] = re.sub(" \(.+\)", "", row['Target gene'])
                
                #building TFs list (nodes with at least one output)
                if row['Regulatory gene'] not in tfs:
                    tfs.append(row['Regulatory gene'])
                
                #building targets list (nodes with at least one input)
                if row['Target gene'] not in targets:
                    targets.append(row['Target gene'])

        # a node to be relevant for attractors landscape calculations needs to have at least one output edge 
        nodes = tfs

        # Identify regulators that are not targets (conserved genes in the dynamics)
        with open("conserved_genes.csv", "w") as file_genes:
            for g in tfs:
                if g not in targets:
                    print(g, file = file_genes)
        
        # Defining sub-network interactions data structures 
        with open(file_name, "r") as f_name:
            reader = csv.DictReader(f_name,  delimiter=';')
            activators = {}
            repressors = {}
            unknown = {}
 
            for row in reader:
                #remove space, and parenthesis since ASSA-PBN does not admit them. Keep only gene names as ID
                row['Regulatory gene'] = re.sub(" \(.+\)", "", row['Regulatory gene'])
                row['Target gene'] = re.sub(" \(.+\)", "", row['Target gene'])
  
                if row['Target gene'] in nodes and row['Regulatory gene'] in nodes:
                    #building dictionaties for each type of interaction
                    if row['Mode of regulation'] == "+":
                        try:
                            activators[row['Target gene']].append(row['Regulatory gene'])
                        except KeyError:
                            activators[row['Target gene']]=[]
                            activators[row['Target gene']].append(row['Regulatory gene'])
                        
                    if row['Mode of regulation'] == "-":
                        try:
                            repressors[row['Target gene']].append(row['Regulatory gene'])
                        except KeyError:
                            repressors[row['Target gene']] = []
                            repressors[row['Target gene']].append(row['Regulatory gene'])

                    if row['Mode of regulation'] == "?" or row['Mode of regulation'] == "d":
                        try:
                            unknown[row['Target gene']].append(row['Regulatory gene'])
                        except KeyError:
                            unknown[row['Target gene']] =[]
                            unknown[row['Target gene']].append(row['Regulatory gene'])      
    
        return([nodes,activators,repressors,unknown])
    except FileNotFoundError:
        sys.exit("Could not open file!")



def print_assa_format(grn, file_name):

    nodes = grn[0] #list
    activators = grn[1] #dict: key = node;  value = list of regulators
    repressors = grn[2] #dict: key = node;  value = list of regulators
    unknown = grn[3] #dict: key = node;  value = list of regulators
    
    with open(file_name, "w") as f_name:
    
        #printing header
        print("//HEADER", file = f_name)
        print("type=synchronous", file = f_name)
        number = len(nodes)
        print(f"n={number}", file = f_name)
        print("perturbation=0.000" ,file = f_name)
        print("",  file = f_name)
        
        #printing nodes names
        print("//NODES NAMES",  file = f_name)
        print("nodeNames",  file = f_name)
        for node in nodes:
            print(node,  file = f_name)
        print("endNodeNames",  file = f_name)
        print(file = f_name)
        # Boolean rules are defined as in Crespo: STEM CELLS 2013;31:2127–2135 :
        # "The logic rule
        # applied by default is the following: if none of its inhibitors and at
        # least one of its activators is active, then a gene becomes active; 
        # otherwise the gene is inactive".

        # If more the one interaction with unknown mode is found, we group their positive interaction 
        # (with probability = 0.5) in one rule,  and their negative interactions 
        # (with probability = 0.5) in another one.

        # Defining Boolean rule:
        print("//BOOLEAN RULES",  file = f_name)
        for node in nodes:
            rule_activ = ""
            rule_repre = ""
            rule_activ_unk = ""
            rule_repre_unk = ""
            # Concatenates sub-strings for activators and repressors:
            try:
                for act in activators[node]:
                    rule_activ = rule_activ + act + "|"
            except KeyError:
                rule_activ
            rule_activ = re.sub("\|$", "", rule_activ) # rule substring with activators
            try:
                for rep in repressors[node]:
                    rule_repre = rule_repre + rep + "|"
            except KeyError:
                rule_repre
            rule_repre = re.sub("\|$", "", rule_repre) # rule substring with repressors

            # Concatenates sub-strings for Unknown interactions:
            try:
                unknown[node] 
                rule_activ_unk = rule_activ
                rule_repre_unk = rule_repre                
                for unk in unknown[node]:
                    rule_activ_unk = rule_activ_unk  + "|" +  unk 
                    rule_repre_unk  = rule_repre_unk  + "|" +  unk
            except KeyError:
                    rule_activ_unk
                    rule_repre_unk
            rule_activ_unk = re.sub("^\|", "", rule_activ_unk)  
            rule_repre_unk = re.sub("^\|", "", rule_repre_unk)
            
            # Printing Boolean rules
            
            
            # Concatenates and print rule if there are no unknown mode interactions (? ou d):
            if rule_activ_unk == "" and rule_repre_unk == "":
                if rule_activ == "" and rule_repre != "":
                    rule = "(!(" + rule_repre + "))"
                elif rule_activ != "" and rule_repre == "":
                    rule = rule_activ
                elif rule_activ != "" and rule_repre != "":
                    rule = "(" + rule_activ + ")" + "&" + "(!(" + rule_repre + "))"
                else:
                    rule = node
                #####
                print("node", node, file = f_name )
                print("1.0 : ", rule, file = f_name)
                #####
            # Concatenates and print rules if there are at least one unknown mode interaction (? ou d):
            else:
                if rule_activ == "" and rule_repre == "":
                    rule1 = "(!(" + rule_repre_unk + "))"
                    rule2 = rule_activ_unk
                elif rule_activ != "" and rule_repre != "":
                    rule1 = "(" + rule_activ_unk + ")" + "&" + "(!(" + rule_repre + "))"
                    rule2 = "(" + rule_activ + ")" + "&" + "(!(" + rule_repre_unk + "))"
                elif rule_activ != "" and rule_repre == "":
                    rule1 = "(" + rule_activ + ")" + "&" + "(!(" + rule_repre_unk + "))"
                    rule2 = rule_activ_unk 
                elif rule_activ == "" and rule_repre != "":
                    rule1 = "(" + rule_activ_unk + ")" + "&" + "(!(" + rule_repre + "))"
                    rule2 = "(!(" + rule_repre_unk + "))"

                ####
                print("node", node, file = f_name )
                print("0.5 :", rule1, file = f_name )
                print("0.5 :", rule2, file = f_name )
                #### 
            #
            print("endNode", file = f_name)
            #

####
main()
####
