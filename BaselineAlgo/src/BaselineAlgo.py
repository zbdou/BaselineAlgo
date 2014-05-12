#-*- coding: UTF-8 -*-
'''
Created on 3/10/2014
       KSP: https://github.com/Pent00/YenKSP
@author: zbdou
'''
import networkx as nx
from collections import deque
import ast
import random
import os
from YenKSP import graph, graphviz

# 2. for link bandwidth, use all simple paths algo.
#     2.1 split path ok
#     2.2 no split path
# 3. for link delay, use networkx shortest path algo.
# 4. 如果求解失败，则要回退，重新转到 节点映射 阶段

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
    nxGraph2KSPDiGraph(g, "virtual").generate()
    
    return g

def build_substrate_network():
    """
    build the physical network, return this physical network.
    test topo:
        see the visio graph.
    """
    MAX_NODES_NO = 26
    g = nx.Graph();
    g.add_nodes_from(range(1, MAX_NODES_NO))
    
    total = random.randint(20, 40)
    while True:
        if total == 0: break
        src = random.randint(1, MAX_NODES_NO-1)    
        dst = random.randint(1, MAX_NODES_NO-1)
        if src == dst or g.has_edge(src, dst): continue
        bw = 10 * random.randint(1, 10)
        g.add_edge(src, dst, {'bw':bw})            
        total -= 1
    
    # cpu(ns), allocated(ns)
    cpus = map(lambda x: 10 * random.randint(1, 10), xrange(1, MAX_NODES_NO))
    print "cpus = ", cpus
    
    for i in xrange(1, len(cpus) + 1):
        g.node[i]['cpu'] = cpus[i-1]
        
    # dump topology
    nxGraph2KSPDiGraph(g, "substrate").generate()
    
    return g

class nxGraph2KSPDiGraph:
    """ a wrapper class, Convert networkx Graph to KSP-like DiGraph for PNG output."""
    def __init__(self, nxGraph, name):
        self.name = name
        self.nxgraph = nxGraph
    
    def build(self):
        digraph = graph.DiGraph()
        for ns in self.nxgraph.nodes():
            node = "{0}({1})".format(ns, self.nxgraph.node[ns]['cpu'])
            digraph.add_node(node)
            for neighbor in nx.all_neighbors(self.nxgraph, ns):
                # add edges
                nei_node = "{0}({1})".format(neighbor, self.nxgraph.node[neighbor]['cpu'])
                digraph.add_edge(node, nei_node, self.nxgraph[ns][neighbor]['bw'])
                digraph.add_edge(nei_node, node, self.nxgraph[neighbor][ns]['bw']) 
        return digraph
    
    def generate(self):
        """
        generate the PNG graph
        """
        assert self.nxgraph, self.name
        digraph = self.build()
        digraph.set_name(self.name)
        digraph.export()

class BaselineAlgo:
    MAX_BW = 1e5
    """Virtual Network Mapping, baseline algorithm"""
    def __init__(self, gv, gs, all_solutions = True, verbose = True, isBW_SP = True):
        """
        @isBW_SP 资源分配方式，按照链路带宽，是否允许split path
        """
        self.gv = gv
        self.gs = gs
        self.all_solutions = all_solutions
        self.verbose = verbose
        self.solution_found = False
        self.split = isBW_SP
        self.solution_counter = 0
        assert self.gv, self.gs
        
    def run(self):
        """
        if all_solutions == False, then return if find the first solution
        if verbose_output == True, then output 1) the total number of iterations, 2) output the solution.
        
        return: all solution if all_solutions = True, otherwise, return the first or None, if no feasible solution found.
        """
        virtual_nodes = deque(self.gv.nodes())
          
        # find the total number of iterations
        total_iterations = reduce(lambda x, y: x*y, map(lambda x: len(self.gv.node[x]['candidates']), virtual_nodes))
        print "The total number of iterations is: {total_iterations}".format(total_iterations=total_iterations)
        print "Find all solutions: {0}, verbose output: {1}".format(self.all_solutions, self.verbose)
        
        self._run()
    
    def _virtual_path_mapping(self, vsrc, vdst, gv, gs):
        """
        @vsrc path's virtaul source node
        @vdst path's virtual sink node
        @gv virtual network
        @gs substrate network
        
        @return [path1:bandwidth1, path2:bandwidth2, ..., ],
                in desc order
        """
        # 1. get the physical psrc, pdst
        psrc, pdst = gv.node[vsrc]['ns'], gv.node[vdst]['ns']
        
        # 2. get all paths from psrc to pdst of the substrate network, P
        paths_mapping = dict()
        paths = nx.all_simple_paths(gs, psrc, pdst, None)
        for path in paths:
            # calculate the path bandwidth based on the link bandwidth
            path_bw = self.MAX_BW
            for link_src, link_dst in zip(path[:-1], path[1:]):
                # print gs[link_src][link_dst]['bw'],
                path_bw = min(gs[link_src][link_dst]['bw'], path_bw)
            # store the path bw into this path
            paths_mapping[str(path)] = path_bw
        
        sorted_paths = sorted(paths_mapping.iteritems(), key = lambda d: d[1], reverse = True)
        # print sorted_paths
        # 3. sorted_paths is sorted by the path's bandwidth
        return sorted_paths
    
    def _run(self):
        '''recursive function'''
        # 0. not_done is empty
        if not self.gv.graph['not_done']:
            self.solution_found = True
            print self.gv.nodes(data=True)

            # 1. for all ev in Ev
            gs = nx.Graph(self.gs)
            gv = nx.Graph(self.gv)
            failed = False
            for src, dst in gv.edges():
                if failed == True: break
                # print "src, dst", src, dst
                
                # 取该虚拟路径(src, dst)对应的所有物理路径，sp
                # get the virtual path's required bandwidth
                vbw = gv[src][dst]['bw']                
                gv[src][dst]['allocated_path'] = list()
                done = False
                while not done:
                    sp = self._virtual_path_mapping(src, dst, gv, gs)
                    if len(sp) == 0:
                        failed = True
                        break
                    # 计算总带宽
                    total_bw = sum(map(lambda d:d[1], sp))
                    path, first_bw = sp[0]
                    
                    if self.split == False and first_bw < vbw:
                        print "no split, and the vbw > first_bw", vbw, first_bw
                        failed = True
                        break
                    elif self.split == True and total_bw < vbw:
                        print "split, but no enough bw, total < vbw", total_bw, vbw
                        failed = True
                        break
                    # 对该虚拟路径(src, dst)依次分配物理路径和带宽，并更新资源视图
                    # print "vbw, first_bw", vbw, first_bw
                    bwmin = min(vbw, first_bw)
                    vbw -= bwmin
                    
                    # gv allocated path 
                    # gs resource update
                    path = ast.literal_eval(path)
                    gv[src][dst]['allocated_path'].append((path, bwmin))                        
                    for link_src, link_dst in zip(path[:-1], path[1:]):
                        gs[link_src][link_dst]['bw'] -= bwmin

                    if vbw == 0:
                        done = True
            if failed == True:
                self.solution_found = False
            else: 
                # a solution found, output the result
                self.solution_counter += 1
                print '-------------------------------------------------'
                print "A solution found: "
                # for each virtual path, draw a plot
                for src, dst in gv.edges():
                    # for each physical path
                    no = 1
                    for path, bw in gv[src][dst]['allocated_path']:
                        g = graph.DiGraph()
                        for node in path:
                            g.add_node(node)
                        for psrc, pdst in zip(path[:-1], path[1:]):
                            g.add_edge(psrc, pdst, bw)
                            g.add_edge(pdst, psrc, bw)
                        name = "solution({0})_{1}-{2}_{3}".format(self.solution_counter, src, dst, no)
                        g.set_name(name)
                        g.export()
                        no += 1
                print '-------------------------------------------------'                        
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

def remove_images(image_dir, dot_dir):
    for d in [image_dir, dot_dir]:
        for f in os.listdir(d):
            tf = os.path.join(d, f)
            if os.path.isfile(tf):
                os.remove(tf)
            
def main():
    remove_images(graphviz.Graphviz._directory_images, graphviz.Graphviz._directory_data)
    Gv = build_virtual_network()    
    Gs = build_substrate_network()
    BaselineAlgo(Gv, Gs, all_solutions = True, isBW_SP = True).run()

if __name__ == '__main__':
    main()

    

