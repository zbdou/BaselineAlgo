#-*- coding: UTF-8 -*-
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
    g = nx.Graph(not_done = deque())
    g.add_nodes_from(range(1, 4))
    g.add_edges_from([(1, 2), (2, 3), (1, 3)])
    g.graph['not_done'] = g.nodes()
    
    # print "not_done", g.graph['not_done']
    
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
    """
    print "------------------------------------------------"
    print g.graph
    print g.nodes(data = True)
    print g.edges(data = True)
    """
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
        
    # dump topology
    """
    print "------------------------------------------------"
    print g.graph
    print g.nodes(data = True)
    print g.edges(data = True)
    """
    return g

class BaseLineAlgo:
    """Virtual Network Mapping, baseline algorithm"""
    def __init__(self, gv, gs):
        self.gv = gv
        self.gs = gs
        assert self.gv, self.gs
        
    def run(self, all_solutions = True,  verbose_output = True):
        """
        if all_solutions == False, then return if find the first solution
        if verbose_output == True, then output 1) the total number of iterations, 2) output the solution.
        
        return: all solution if all_solutions = True, otherwise, return the first or None, if no feasible solution found.
        """
        virtual_nodes = deque(self.gv.nodes())
          
        # find the total number of iterations
        total_iterations = reduce(lambda x, y: x+y, map(lambda x: len(self.gv.node[x]['candidates']), virtual_nodes))
        print "The total number of iterations is: {total_iterations}".format(total_iterations=total_iterations)
        
        self._run()
        
    def _run(self):
        '''recursive function'''
        # 0. not_done is empty
        if not self.gv.graph['not_done']:
            # print self.gs.nodes(data=True)
            print self.gv.nodes(data=True)
            return
        # 1. get the first nv that is in not_done
        nv = self.gv.graph['not_done'].pop(0)
        
        # 2. get the map(nv)
        candidates = self.gv.node[nv]['candidates']
        candidates_mapping = dict(zip(candidates, map(lambda x: self.gs.node[x]['cpu'], candidates)))
        
        # print candidates
        # print candidates_mapping
        desc_candi_list = sorted(candidates_mapping.iteritems(), key = lambda d:d[1], reverse = True)
        #print desc_candi_list
        assert desc_candi_list
        
        nv_cpu = self.gv.node[nv]['cpu']
        
        # step 3.        
        for ns, ns_cpu in desc_candi_list:
            # step 4.
            if nv_cpu > ns_cpu: 
                self.gv.graph['not_done'].insert(0, nv);
                return
            
            # step 5. update the resource allocation database
            self.gs.node[ns]['cpu'] -= nv_cpu
            self.gv.node[nv]['ns'] = ns
            
            # step 6. 
            self._run()
            
            # step 7. restore the resource allocation database
            self.gs.node[ns]['cpu'] += nv_cpu
            self.gv.node[nv]['ns'] = None
        
        self.gv.graph['not_done'].insert(0, nv);                             
        return

if __name__ == '__main__':
    Gv = build_virtual_network()    
    Gs = build_substrate_network()
    
    BaseLineAlgo(Gv, Gs).run()
    print Gv.graph
    

