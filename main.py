import networkx as nx
import matplotlib.pyplot as plt
import os
import time

def add_edge(MBG, gene1, gene2, origin):
    if MBG.has_edge(gene1, gene2):
        MBG[gene1][gene2]["value"] = MBG[gene1][gene2]["value"] + 1
        MBG[gene1][gene2]["genoma"].append(origin)
    else:
        MBG.add_edge(gene1, gene2, value = 1, genoma = [origin])

def create_MBG(MBG): # Read input file ("input.txt") to create the Multiple Breakpoint Graph
                     # The input needs to have a number N at the first line, representing the number of genes and each following line a different genoma
                     # each genoma consist of a sequence from 1 to N separated by a space and if a gene is reversed the number should be negative

    f_in = open("input.txt", "r")
    lines = f_in.readlines()
    size = 0
    genoma = -1

    # For each line it will add all edges and by consequence all nodes to the MBG
    for line in lines:
        f.write(str(line))
        genoma = genoma + 1
        line = line.split()
        if len(line) == 1:
            size = int(line[0])
            
        else:
            vet = []
            for i in range(len(line)):

                if line[i] == "(":
                    vet = []
                
                elif line[i] == ")":
                    for i in range(len(vet)):
                        num1 = vet[i]
                        num2 = vet[0]
                        if i != len(vet) - 1:
                            num2 = vet[i + 1]
                    
                        if num1 > 0:
                            if num2 > 0:
                                add_edge(MBG, str(num1) + "H", str(num2) + "T", genoma)
                            else:
                                add_edge(MBG, str(num1) + "H", str(-num2) + "H", genoma)
                        else:
                            if num2 > 0:
                                add_edge(MBG, str(-num1) + "T", str(num2) + "T", genoma)
                            else:
                                add_edge(MBG, str(-num1) + "T", str(-num2) + "H", genoma)

                else:
                    vet.append(int(line[i]))

    f_in.close()
    return size

def create_QG(QG, size): # Create Quadruple Graph, all nodes are ordered from lower to higher ie. 1T 2H or 1H 1T
    for i in range(size - 1):
        QG.add_node(str(i + 1)  + "H " + str(i + 1) + "T", value = 0)
        for j in range(i + 1, size):
            QG.add_node(str(i + 1)  + "H " + str(j + 1) + "H", value = 0)
            QG.add_node(str(i + 1)  + "T " + str(j + 1) + "H", value = 0)
            QG.add_node(str(i + 1)  + "H " + str(j + 1) + "T", value = 0)
            QG.add_node(str(i + 1)  + "T " + str(j + 1) + "T", value = 0)
    QG.add_node(str(size)  + "H " + str(size) + "T", value = 0)

    # Add the value to all nodes that represents an edge at the MBG, in this case a 2-cycle
    nodes = MBG.edges()
    node1 = 0
    node2 = 0

    for i in (nodes):
        node1 = i[0]
        node2 = i[1]

        if node1 > node2:
            node1, node2 = node2, node1
        
        name = node1 + " " + node2

        QG.nodes[name]["value"] = QG.nodes[name]["value"] + 1
        QG.nodes[name]["genoma"] = MBG.edges[str(node1), str(node2)]["genoma"]

def add_edge_QG(QG, MBG): # Add all edges to the Quadruple Graph
    all_nodes = list(QG.nodes)

    # This for will look for all nodes and if the sum of its value and the amount of edges connected to that node is less then 3,
    # Than there is at least of edge that needs to be added to the graph
    for i in range(len(all_nodes)):

        neighbors = list(QG[all_nodes[i]])
        if QG.nodes[all_nodes[i]]["value"] + len(neighbors) == 3:
            continue

        colors = [1, 2, 3]
        if QG.nodes[all_nodes[i]]["value"] > 0:
            for j in range(QG.nodes[all_nodes[i]]["value"]):
                colors.remove(QG.nodes[all_nodes[i]]["genoma"][j])
        
        for j in range (len(neighbors)):
            for k in QG.edges[all_nodes[i], neighbors[j]]["genoma"]:
                colors.remove(k)
        
        nodes = all_nodes[i].split()
        node1 = nodes[0]
        node2 = nodes[1]

        node_to_add1 = 0
        node_to_add2 = 0

        for j in colors:
            for k in list(MBG[node1]):
               if j in MBG[node1][k]["genoma"]:
                   node_to_add1 = k
                   break

            for k in list(MBG[node2]):
               if j in MBG[node2][k]["genoma"]:
                   node_to_add2 = k
                   break
               
            if node_to_add1 > node_to_add2:
                node_to_add1, node_to_add2 = node_to_add2, node_to_add1

            if QG.has_edge(node1 + " " + node2, node_to_add1 + " " + node_to_add2):
                QG[node1 + " " + node2][node_to_add1 + " " + node_to_add2]["value"] = QG[node1 + " " + node2][node_to_add1 + " " + node_to_add2]["value"] + 1
                QG[node1 + " " + node2][node_to_add1 + " " + node_to_add2]["genoma"].append(j)
            else:
                QG.add_edge(node1 + " " + node2, node_to_add1 + " " + node_to_add2, genoma = [j], value = 1)

def possible_solutions(QG, size): # Brute force attempt to select all possible solutions of max weight
                                  # If there are more than one optimal solution, this will return all viable ones.
    all_nodes = list(QG.nodes)
    total_weight = 0
    aux = []
    divided_subset = []
    subset = []
    solution = []
    i = 0
    while len(aux) < size:
        splited_i = all_nodes[i].split()
        if splited_i[0] not in divided_subset and splited_i[1] not in divided_subset:
            aux.append(i)
            divided_subset.append(splited_i[0])
            divided_subset.append(splited_i[1])
            subset.append(all_nodes[i])
        i += 1
    
    len_all_nodes = len(all_nodes)

    while aux[0] + size < len_all_nodes:
        i = size - 1
        while i < size:
            if aux[0] + size >= len_all_nodes:
                break
            if aux[i] == -1:
                aux[i] = aux[i- 1] + 1
                while True:
                    if aux[i] < len_all_nodes:
                        splited_i = all_nodes[aux[i]].split()
                        if splited_i[0] in divided_subset or splited_i[1] in divided_subset:
                            aux[i] += 1
                        else:
                            break
                    else:
                        break
            else:
                divided_subset.pop(i*2)
                divided_subset.pop(i*2)
                subset.pop(i)
                aux[i] += 1
                while True:
                    if aux[i] < len_all_nodes:
                        splited_i = all_nodes[aux[i]].split()
                        if splited_i[0] in divided_subset or splited_i[1] in divided_subset:
                            aux[i] += 1
                        else:
                            break
                    else:
                        break

            if aux[i] == len_all_nodes:
                aux[i] = -1
                i -= 1
            else:
                splited_i = all_nodes[aux[i]].split()
                divided_subset.append(splited_i[0])
                divided_subset.append(splited_i[1])
                subset.append(all_nodes[aux[i]])
                i += 1

        if len(divided_subset) == size + size:
            current_weight = 0
            for i in range (len(subset)):
                current_weight = QG.nodes[subset[i]]["value"] + current_weight
                for j in range(i + 1, len(subset)):
                    if QG.has_edge(subset[i], subset[j]):
                        current_weight = QG[subset[i]][subset[j]]["value"] + current_weight
                    
            if current_weight > total_weight:
                total_weight = current_weight
                solution = [subset.copy()]
            elif current_weight == total_weight:
                solution.append(subset.copy())

    return total_weight, solution


# Create directory to add the solutions
try:
    os.mkdir("solutions")
except FileExistsError:
    pass

instance = 0
while True:
    try:
        os.mkdir("solutions/instance_" + str(instance))
        break
    except FileExistsError:
        instance += 1

start_time = time.time()
f = open("solutions/instance_" + str(instance) + "/output.txt", "w")

# MBG refers to the multiple breakpoint graph, each edge has two attributes, called genoma, a list of all genomas it belongs (in case there are parallel edges)
# value, to represent the the quantity of genomas that a single edge belongs to.
MBG = nx.Graph()
size = create_MBG(MBG)

# Add colors for edges
edges = MBG.edges()
color_map = []
for i in edges:
    if MBG[i[0]][i[1]]["genoma"][0] == 1:
        color_map.append("red")
    elif MBG[i[0]][i[1]]["genoma"][0] == 2:
        color_map.append("blue")
    else:
        color_map.append("green")

# A Draw of the MBG
nx.draw_networkx(MBG, edge_color = color_map)
plt.savefig("solutions/instance_" + str(instance) + "/MBG.png")
plt.close()

QG = nx.Graph()
create_QG(QG, size)

add_edge_QG(QG, MBG)

weight, solution = possible_solutions(QG, size)


f.write("Weight = " + str(weight) + "\n")

f.write("--- %s seconds ---" % (time.time() - start_time) + "\n")
f.write("Optimal solutions: \n")
# Each solution will draw the QG, the nodes is red represents the solution while the ones in yellow all other nodes
quant = 0
for i in solution:
    f.write(str(i) + "\n")
    quant += 1
    name = "solutions/instance_" + str(instance) + "/solution_" + str(quant)+  ".png"
    color_map = ["red" if node in i else "yellow" for node in QG]
    edges = QG.edges()
    color_map_edges = []
    for j in edges:
        if QG[j[0]][j[1]]["value"] == 1:
            color_map_edges.append("red")
        else:
            color_map_edges.append("green")
    nx.draw_networkx(QG, node_color = color_map, edge_color = color_map_edges)
    plt.savefig(name)
    plt.show()

f.close()