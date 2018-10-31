#! /usr/bin/env python


import time

import networkx as nx
import numpy as np
import copy
import scipy.stats
from collections import defaultdict
import csv
import sys

'''
# =============================================================================
def print_usage():

    print (' ')
    print ('        usage: ./DIAMOnD network_file seed_file n alpha(optional) outfile_name (optional)')
    print ('        -----------------------------------------------------------------')
    print ('        network_file : The edgelist must be provided as any delimiter-separated')
    print ('                       table. Make sure the delimiter does not exit in gene IDs')
    print ('                       and is consistent across the file.' )
    print ('                       The first two columns of the table will be')
    print ('                       interpreted as an interaction gene1 <==> gene2')
    print ('        seed_file    : table containing the seed genes (if table contains')
    print ('                       more than one column they must be tab-separated;')
    print ('                       the first column will be used only)')
    print ('        n            : desired number of DIAMOnD genes, 200 is a reasonable')
    print ('                       starting point.')
    print ('        alpha        : an integer representing weight of the seeds,default')
    print ('                       value is set to 1')
    print ('        outfile_name : results will be saved under this file name')
    print ('                       by default the outfile_name is set to "first_n_added_nodes_weight_alpha.txt"')
    print (' ')


# =============================================================================
'''
def check_input_style(input_list):
    try:
        network_edgelist_file = input_list[1]
        seeds_file = input_list[2]
        max_number_of_added_nodes = int(input_list[3])
    # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return

    alpha = 1
    outfile_name = 'first_%d_added_nodes_weight_%d.txt'%(max_number_of_added_nodes,alpha)

    if len(input_list)==5:
        try:
            alpha = int(input_list[4])
            outfile_name = 'first_%d_added_weight_%d.txt'%(max_number_of_added_nodes,alpha)
        except:
            outfile_name = input_list[4]

    if len(input_list)==6:
        try:
            alpha = int(input_list[4])
            outfile_name = input_list[5]
        except:
            print_usage()
            sys.exit(0)
            return
    return network_edgelist_file,seeds_file,max_number_of_added_nodes,alpha,outfile_name

# =============================================================================
#def read_input(network_file,seed_file):
"""

    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file,'r'):
        if line[0]=='#':
            continue
        else:
            dialect = sniffer.sniff(line)
            line_delimiter = dialect.delimiter
	        break
    if line_delimiter == None:
        print ('network_file format not correct')
        sys.exit(0)


    # read the network:
    G = nx.Graph()
    for line in open(network_file,'r'):
        # lines starting with '#' will be ignored
        if line[0]=='#':
            continue
        # The first two columns in the line will be interpreted as an
        # interaction gene1 <=> gene2
        #line_data   = line.strip().split('\t')
        line_data = line.strip().split(line_delimiter)
        node1 = line_data[0]
        node2 = line_data[1]
        G.add_edge(node1,node2)

    # read the seed genes:
    seed_genes = set()
    for line in open(seed_file,'r'):
        # lines starting with '#' will be ignored
        if line[0]=='#':
            continue
        # the first column in the line will be interpreted as a seed
        # gene:
        line_data = line.strip().split('\t')
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)

    return G,seed_genes
"""

# ================================================================================
def compute_all_gamma_ln(N):
    """
    precomputes all logarithmic gammas
    """
    gamma_ln = {}
    for i in range(1,N+1):
        gamma_ln[i] = scipy.special.gammaln(i)

    return gamma_ln

# =============================================================================
def logchoose(n, k, gamma_ln):
    if n-k+1 <= 0:
        return scipy.infty
    lgn1  = gamma_ln[n+1]
    lgk1  = gamma_ln[k+1]
    lgnk1 = gamma_ln[n-k+1]
    return lgn1 - [lgnk1 + lgk1]

# =============================================================================
def gauss_hypergeom(x, r, b, n, gamma_ln):
    return np.exp(logchoose(r, x, gamma_ln) +
                  logchoose(b, n-x, gamma_ln) -
                  logchoose(r+b, n, gamma_ln))

# =============================================================================
def pvalue(kb, k, N, s, gamma_ln):
    """
    -------------------------------------------------------------------
    Computes the p-value for a node that has kb out of k links to
    seeds, given that there's a total of s sees in a network of N nodes.

    p-val = \sum_{n=kb}^{k} HypergemetricPDF(n,k,N,s)
    -------------------------------------------------------------------
    """
    p = 0.0
    for n in range(kb,k+1):
        if n > s:
            break
        prob = gauss_hypergeom(n, s, N-s, k, gamma_ln)
        # print prob
        p += prob


    if p > 1:
        return 1
    else:
        return float(p)

# =============================================================================
def get_neighbors_and_degrees(G):

    neighbors,all_degrees = {},{}
    for node in G.nodes():
        nn = set(G.neighbors(node))
        neighbors[node] = nn
        all_degrees[node] = G.degree(node)

    return neighbors,all_degrees

# =============================================================================
# Reduce number of calculations
# =============================================================================
def reduce_not_in_cluster_nodes(all_degrees,neighbors,G,not_in_cluster,cluster_nodes,alpha):
    reduced_not_in_cluster = {}
    kb2k = defaultdict(dict)
    for node in not_in_cluster:

        k = all_degrees[node]
        kb = 0
        # Going through all neighbors and counting the number of module neighbors
        for neighbor in neighbors[node]:
            if neighbor in cluster_nodes:
                kb += 1

        #adding wights to the the edges connected to seeds
        k += (alpha-1)*kb
        kb += (alpha-1)*kb
        kb2k[kb][k] =node

    # Going to choose the node with largest kb, given k
    k2kb = defaultdict(dict)
    for kb,k2node in kb2k.items():
        min_k = min(k2node.keys())
        node = k2node[min_k]
        k2kb[min_k][kb] = node

    for k,kb2node in k2kb.items():
        max_kb = max(kb2node.keys())
        node = kb2node[max_kb]
        reduced_not_in_cluster[node] =(max_kb,k)

    return reduced_not_in_cluster

#======================================================================================
#   C O R E    A L G O R I T H M
#======================================================================================
def diamond_iteration_of_first_X_nodes(G,S,X,alpha):

    """

    Parameters:
    ----------
    - G:     graph
    - S:     seeds
    - X:     the number of iterations, i.e only the first X gened will be
             pulled in
    - alpha: seeds weight

    Returns:
    --------

    - added_nodes: ordered list of nodes in the order by which they
      are agglomerated. Each entry has 4 info:

      * name : dito
      * k    : degree of the node
      * kb   : number of +1 neighbors
      * p    : p-value at agglomeration

    """

    N = G.number_of_nodes()

    added_nodes = []


    # ------------------------------------------------------------------
    # Setting up dictionaries with all neighbor lists
    # and all degrees
    # ------------------------------------------------------------------
    neighbors,all_degrees = get_neighbors_and_degrees(G)

    # ------------------------------------------------------------------
    # Setting up initial set of nodes in cluster
    # ------------------------------------------------------------------

    cluster_nodes = set(S)
    not_in_cluster = set()
    s0 = len(cluster_nodes)

    s0 += (alpha-1)*s0
    N +=(alpha-1)*s0

    # ------------------------------------------------------------------
    # precompute the logarithmic gamma functions
    # ------------------------------------------------------------------
    gamma_ln = compute_all_gamma_ln(N+1)

    # ------------------------------------------------------------------
    # Setting initial set of nodes not in cluster
    # ------------------------------------------------------------------
    for node in cluster_nodes:
        not_in_cluster |= neighbors[node]
    not_in_cluster -= cluster_nodes


    # ------------------------------------------------------------------
    #
    # M A I N     L O O P
    #
    # ------------------------------------------------------------------

    all_p = {}

    while len(added_nodes) < X:

        # ------------------------------------------------------------------
        #
        # Going through all nodes that are not in the cluster yet and
        # record k, kb and p
        #
        # ------------------------------------------------------------------

        info = {}

        pmin = 10
        next_node = 'nix'
        reduced_not_in_cluster = reduce_not_in_cluster_nodes(all_degrees,
                                                             neighbors,G,
                                                             not_in_cluster,
                                                             cluster_nodes,alpha)

        for node,kbk in reduced_not_in_cluster.items():
            # Getting the p-value of this kb,k
            # combination and save it in all_p, so computing it only once!
            kb,k = kbk
            try:
                p = all_p[(k,kb,s0)]
            except KeyError:
                p = pvalue(kb, k, N, s0, gamma_ln)
                all_p[(k,kb,s0)] = p

            # recording the node with smallest p-value
            if p < pmin:
                pmin = p
                next_node = node

            info[node] = (k,kb,p)

        # ---------------------------------------------------------------------
        # Adding node with smallest p-value to the list of aaglomerated nodes
        # ---------------------------------------------------------------------
        added_nodes.append((next_node,
                            info[next_node][0],
                            info[next_node][1],
                            info[next_node][2]))

        # Updating the list of cluster nodes and s0
        cluster_nodes.add(next_node)
        s0 = len(cluster_nodes)
        not_in_cluster |= ( neighbors[next_node] - cluster_nodes )
        not_in_cluster.remove(next_node)

    return added_nodes

# ===========================================================================
#
#   M A I N    D I A M O n D    A L G O R I T H M
#
# ===========================================================================
def DIAMOnD(G_original,seed_genes,max_number_of_added_nodes,alpha,outfile = None):

    """
    Runs the DIAMOnD algorithm

    Input:
    ------
     - G_original :
             The network
     - seed_genes :
             a set of seed genes
     - max_number_of_added_nodes:
             after how many added nodes should the algorithm stop
     - alpha:
             given weight to the sees
     - outfile:
             filename for the output generates by the algorithm,
             if not given the program will name it 'first_x_added_nodes.txt'

     Returns:
     --------
      - added_nodes: A list with 4 entries at each element:
            * name : name of the node
            * k    : degree of the node
            * kb   : number of neighbors that are part of the module (at agglomeration)
            * p    : connectivity p-value at agglomeration
      -
    """

    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network
    ignored_genes = seed_genes.difference(all_genes_in_network)

    #if len(disease_genes) != len(seed_genes):
        #print "DIAMOnD(): ignoring %s of %s seed genes that are not in the network" %(
            #len(seed_genes - all_genes_in_network), len(seed_genes))
    # 2. agglomeration algorithm.
    added_nodes = diamond_iteration_of_first_X_nodes(G_original,
                                                     disease_genes,
                                                     max_number_of_added_nodes,alpha)
    # 3. saving the results #Changed to fit the format of PASCAL
    #with open(outfile,'a') as fout:
        #print>>fout,'\t'.join(['DIAMOnD']) + '\t' + '\t', #'#rank','DIAMOnD_node'])
        #rank = 0
    #for DIAMOnD_node_info in added_nodes:
        #rank += 1
        #DIAMOnD_node = DIAMOnD_node_info[0]
        #p = float(DIAMOnD_node_info[3])
            #print>>fout,'\t'.join(map(str,([DIAMOnD_node]))) + '\t',
        #print map(str,([DIAMOnD_node]))


    nodes_name = [str(item[0]) for item in added_nodes]
    nodes_degree = [str(item[1]) for item in added_nodes]
    nodes_n_neighbors = [str(item[2]) for item in added_nodes]
    nodes_p = [str(item[3]) for item in added_nodes]
    #print (', '.join(disease_genes))
    #print (', '.join(ignored_genes))
    #print (', '.join(nodes_name))
    #print (', '.join(nodes_degree))
    #print (', '.join(nodes_n_neighbors))
    #print (', '.join(nodes_p))

	#print>>fout,'\n',

    #final_result = list[added_nodes, disease_genes, ignored_genes]

    return [added_nodes, list(disease_genes), list(ignored_genes)]


# ===========================================================================
#
# "Hey Ho, Let's go!" -- The Ramones (1976)
#
# ===========================================================================

def to_set(list_r):
	pyset = set(list_r)
	return pyset

def to_graph(array_r):
	G = nx.Graph()
	for row in array_r:
		node1 = row[0]
		node2 = row[1]
		G.add_edge(node1,node2)
	
	
	return G
