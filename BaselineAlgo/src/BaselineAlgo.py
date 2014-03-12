#-*- coding: UTF-8 -*-
'''
Created on 3/10/2014
       KSP: https://github.com/Pent00/YenKSP
@author: zbdou
'''
import networkx as nx
from collections import deque
import matplotlib.pyplot as plt

from YenKSP import *

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

class BaselineAlgo:
    """Virtual Network Mapping, baseline algorithm"""
    def __init__(self, gv, gs, all_solutions = True, verbose = True):
        self.gv = gv
        self.gs = gs
        self.all_solutions = all_solutions
        self.verbose = verbose
        self.solution_found = False
        assert self.gv, self.gs
        
    def run(self):
        """
        if all_solutions == False, then return if find the first solution
        if verbose_output == True, then output 1) the total number of iterations, 2) output the solution.
        
        return: all solution if all_solutions = True, otherwise, return the first or None, if no feasible solution found.
        """
        virtual_nodes = deque(self.gv.nodes())
          
        # find the total number of iterations
        total_iterations = reduce(lambda x, y: x+y, map(lambda x: len(self.gv.node[x]['candidates']), virtual_nodes))
        print "The total number of iterations is: {total_iterations}".format(total_iterations=total_iterations)
        print "Find all solutions: {0}, verbose output: {1}".format(self.all_solutions, self.verbose)
        
        self._run()
        
    def _run(self):
        '''recursive function'''
        # 0. not_done is empty
        if not self.gv.graph['not_done']:
            self.solution_found = True
            print self.gv.nodes(data=True)

            # FIXME: how to define the cost?            
            # build the DiGraph
            digraph = graph.DiGraph()
            for ns in self.gs.nodes():
                digraph.add_node(ns)
                for neighbor in nx.all_neighbors(self.gs, ns):
                    # add edges
                    digraph.add_edge(ns, neighbor, 1.0/self.gs[ns][neighbor]['bw'])
                    digraph.add_edge(neighbor, ns, 1.0/self.gs[neighbor][ns]['bw'])
                    
            # digraph.save()
            # digraph.export()
            
            # K-shortest path algo.
            # 1. for all ev in Ev
            for src, dst in self.gv.edges():
                ns_src, ns_dst = self.gv.node[src]['ns'], self.gv.node[dst]['ns']
                paint = digraph.painter()
                paint.set_source_sink(str(ns_src), str(ns_dst))
                digraph.set_name("{0}_{1}".format(ns_src, ns_dst))
                digraph.export(False, paint)
                
                items = algorithms.ksp_yen(digraph, src, dst, 3)
                for path in items:
                    #print path
                    print "Cost:{0}\t{1}".format(path['cost'], "->".join(str(i) for i in (path['path'])))
            
            
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
            
            # check if we have found one solution and do not expect all solutions.
            if self.solution_found and not self.all_solutions:
                self.gv.graph['not_done'].insert(0, nv);                             
                return
            # step 7. restore the resource allocation database
            self.gs.node[ns]['cpu'] += nv_cpu
            self.gv.node[nv]['ns'] = None
        
        self.gv.graph['not_done'].insert(0, nv);                             
        return

def main():
    Gv = build_virtual_network()    
    Gs = build_substrate_network()
    
    BaselineAlgo(Gv, Gs, all_solutions = False).run()
    
    print Gv.graph
    
        
if __name__ == '__main__':
    main()

    

