'''
Created on 3/10/2014

@author: zbdou
'''
import networkx as nx
from collections import deque
import matplotlib.pyplot as plt

def build_virtual_network():
    """
    build the virtual network.
    test topo:
        see the visio graph.
    """
    g = nx.Graph(done = deque())
    g.add_nodes_from(range(1, 4))
    g.add_edges_from([(1, 2), (2, 3), (1, 3)])
    
    # b(ev(i, j) )
    g[1][2]['bw'] = 30
    g[1][3]['bw'] = 50
    g[2][3]['bw'] = 10
    
    # cpu(nv), map(nv), phy_nodes(nv)
    cpus = [20, 60, 30]
    map_nv = [[5], [1, 4], [2, 3, 4]]
    for i in xrange(1, len(cpus) + 1):
        g.node[i]['cpu'] = cpus[i-1]
        g.node[i]['candidates'] = map_nv[i-1]
        g.node[i]['ns'] = None
        
    # dump topology
    print "------------------------------------------------"
    print g.graph
    print g.nodes(data = True)
    print g.edges(data = True)
    
    return g

def build_substrate_network():
    """
    build the physical network, return this physical network.
    test topo:
        see the visio graph.
    """
    g = nx.Graph()
    g.add_nodes_from(range(1, 6))
    g.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 1), (2, 5), (3, 5)])
    
    # b(es(i,j))
    g[1][2]['bw'] = 20
    g[2][3]['bw'] = 70
    g[3][4]['bw'] = 30
    g[1][4]['bw'] = 40
    g[2][5]['bw'] = 100
    g[3][5]['bw'] = 80
    
    # cpu(ns), allocated(ns)
    cpus = [100, 50, 30, 80, 20]
    for i in xrange(1, len(cpus) + 1):
        g.node[i]['cpu'] = cpus[i-1]
        g.node[i]['nv'] = deque()
        
    # dump topology
    print "------------------------------------------------"
    print g.graph
    print g.nodes(data = True)
    print g.edges(data = True)
    
    return g



if __name__ == '__main__':
    Gs = build_substrate_network()
    Gv = build_virtual_network()
    
    

