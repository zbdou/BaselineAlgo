#-*- coding: UTF-8 -*-
'''
Created on 3/10/2014
       KSP: https://github.com/Pent00/YenKSP
@author: zbdou
'''
import networkx as nx
from collections import deque
import random, os, math, string
import matplotlib.pyplot as plt
from YenKSP import algorithms, graph

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

    # assign edge bandwidth (link resource)
    g[1][2]['bw'] = 30
    g[1][3]['bw'] = 26
    g[2][3]['bw'] = 50
    
    # assign node computation resource
    cpus = [20, 60, 30]
    map_nv = [[1,2], [1,3,4,5], [2,5]]
    for i in xrange(1, len(cpus) + 1):
        g.node[i]['cpu'] = cpus[i-1]
        g.node[i]['candidates'] = map_nv[i-1]
        g.node[i]['ns'] = None
        
    # dump topology
    pos = nx.shell_layout(g)
    nx.draw_networkx_nodes(g, pos, node_size=1500)
    nx.draw_networkx_edges(g, pos)
    
    node_labels = {n: "{0}({1})".format(n, g.node[n]['cpu']) for n in g.nodes()}    
    nx.draw_networkx_labels(g, pos,labels=node_labels)
    
    edge_labels = {(u,v): bw['bw'] for u,v,bw in g.edges(data=True)}
    nx.draw_networkx_edge_labels(g, pos, edge_labels, label_pos=0.3)
    
    plt.axis('off')
    plt.savefig("virtual", transparent=True, dpi=300)
    plt.clf()
        
    return g

def bandwidth_generator(bmin = 1, bmax = 100, times10 = True):
    assert bmin <= bmax
    r = random.randint(bmin, bmax)
    if times10:
        return r - r % 10 + 10
    return r

def delay_generator(dmin = 100, dmax = 100, times10 = True):
    assert dmin <= dmax
    r = random.randint(dmin, dmax)    
    if times10:
        return r - r % 10 + 10
    return r
    
def build_substrate_network(nodes_no = 10, max_edge_degree = 1):
    """
    build the substrate network which is of circle-like.
    @nodes_no 物理网络的节点总数
    @max_edge_degree 每个节点的最大度数，（不包含与前驱、后继节点相连的两条链路）
    """
    MAX_NODES_NO = nodes_no + 1
    NO_BASE = MAX_NODES_NO - 1
    g = nx.Graph()
    g.add_nodes_from(range(1, MAX_NODES_NO))
    
    # cpu resource, = 100, by default
    cpus = map(lambda x: 100, xrange(1, MAX_NODES_NO))
    #print "cpus = ", cpus
    
    for i in xrange(1, len(cpus) + 1):
        g.node[i]['cpu'] = cpus[i-1]

    # build edge and weight (bandwidth & delay)
    for i in xrange(1, MAX_NODES_NO):
        # excluding the direct (pre, i), (i, next) links
        dst_no = random.randint(0, max_edge_degree)
        pre = (i - 1 + NO_BASE - 1) % NO_BASE + 1
        nxt = (i - 1 + NO_BASE + 1) % NO_BASE + 1
        
        if not g.has_edge(pre, i):
            g.add_edge(pre, i, {'bw': bandwidth_generator(), 'dy': delay_generator()})
        if not g.has_edge(i, nxt):
            g.add_edge(i, nxt, {'bw': bandwidth_generator(), 'dy': delay_generator()})
            
        while dst_no > 0:
            dst = random.randint(1, NO_BASE)
            if not g.has_edge(i, dst):
                g.add_edge(i, dst, {'bw': bandwidth_generator(), 'dy': delay_generator()})
            dst_no -= 1

    pos = nx.circular_layout(g, scale=2)

    node_labels = {n: "{0}({1})".format(n, g.node[n]['cpu']) for n in g.nodes()}
    nx.draw_networkx_nodes(g, pos, node_sizes=5000)
    nx.draw_networkx_edges(g, pos)
    nx.draw_networkx_labels(g, pos, labels=node_labels, font_size=5)
    
    edge_labels = {(u,v): bw['bw'] for u, v, bw in g.edges(data=True)}
    nx.draw_networkx_edge_labels(g, pos, edge_labels, label_pos=0.3, font_size=8)
    
    plt.axis('off')       
    plt.savefig("substrate", transparent=True, dpi=300)
    plt.clf()

    return g

def graph2KSPGraph(g):
    g2 = graph.DiGraph()
    def add_edge(s, d, cost):
        g2.add_edge(s, d, cost)
        g2.add_edge(d, s, cost)
        
    map(lambda x: g2.add_node(x), g.nodes_iter())
    map(lambda x: add_edge(x[0], x[1], x[2]['dy']), g.edges_iter(data=True))
        
    #g2.set_name("aaa")
    #g2.export()
    return g2
        
class BaselineAlgo:
    MAX_BW = 1e5
    """Virtual Network Mapping, baseline algorithm"""
    def __init__(self, gv, gs, all_solutions = True, max_split_no = 1, max_k = 5):
        """
        @all_solutions: 是否求解所有的节点映射方案
        @max_split_no: 是否允许路径分割映射，= 1 表示不允许，> 1 表示一条虚拟路径可以映射到 @max_split_no 条物理路径
        @max_k: 仅在 @max_split_no == 1 时有意义，表示按照带宽作为物理路径的cost时，搜索空间，
                                       即最多搜索 @max_k 条物理cost最短的路径（KSP，K Shortest Paths）
        """
        self.gv = gv
        self.gs = gs
        self.all_solutions = all_solutions
        self.solution_found = False
        self.solution_counter = 0
        self.max_split_no = max_split_no
        self.max_k = max_k
        self.gs_snapshot = None
        assert self.gv, self.gs
        
    def run(self):
        """
        if all_solutions == False, then return if find the first solution
        return: all solution if all_solutions = True, otherwise, return the first or None, if no feasible solution found.
        """
        virtual_nodes = deque(self.gv.nodes())
          
        # find the total number of iterations
        total_iterations = reduce(lambda x, y: x*y, map(lambda x: len(self.gv.node[x]['candidates']), virtual_nodes))
        print "The total number of iterations is: {0}".format(total_iterations)
        print "Find all solutions: {0}".format(self.all_solutions)
        
        self._run()
    
    def _virtual_path_mapping(self, vsrc, vdst, gv, gs):
        """
        @vsrc path's virtaul source node
        @vdst path's virtual sink node
        @gv virtual network
        @gs substrate network
        
        @return list of paths:
            [{'path': [a,b,c], 'cost': cost, 'bw': path bandwidth}, {}, {}, ...]
        """
        # 1. get the physical psrc, pdst
        psrc, pdst = gv.node[vsrc]['ns'], gv.node[vdst]['ns']
        
        if nx.has_path(gs, psrc, pdst) == False: return None
        
        # check the abnormal case: psrc = pdst
        if psrc == pdst:
            return [{'path':[psrc, pdst], 'cost': 0, 'bw': self.MAX_BW}]
        
        if self.max_split_no != 1: 
            max_k = 1
        else:
            max_k = self.max_k
            
        paths= algorithms.ksp_yen(graph2KSPGraph(gs), psrc, pdst, max_k)

        # calculate the path bandwidth based on the link bandwidth
        for item in paths:
            path = item['path']
            bws = map(lambda src, dst: gs[src][dst]['bw'], path[:-1], path[1:])
            item['bw'] = min(bws)
        
        #print "path2 = ", paths    
        return paths

    def get_gs_snapshot(self):
        return self.gs_snapshot
    
    def dump_virtual_path(self, gs, gv, src, dst, path):
        # for each physical path
        g = nx.Graph()
                 
        nodes = list()
        map(lambda x: nodes.extend(x[0]), gv[src][dst]['allocated_path'])
        g.add_nodes_from(list(set(nodes)))
        
        pos = nx.circular_layout(g)
        nx.draw_networkx_nodes(g, node_color='r', pos=pos)
        colors = "bgrcmyk"
        for index, (path, bw, delay) in enumerate(gv[src][dst]['allocated_path']):
            edges = zip(path[:-1], path[1:])
            edgelabels = {edge:bw for edge in edges}
            psrc, pdst = path[0], path[-1]
            nx.draw_networkx_nodes(g, nodelist=[psrc, pdst], node_color='y', pos=pos)
            nx.draw_networkx_edges(g, edgelist=edges, pos=pos, edge_color=colors[index % len(colors)], label=str(path) + ", " + str(bw) + ", " + str(int(delay)))
            nx.draw_networkx_edge_labels(g, pos, edgelabels, label_pos=0.3)
            
            
            print "{0}->{1}->{2}, {3}, {4}".format(src, "->".join([str(x) for x in path]), dst, bw, delay)
        
        nx.draw_networkx_labels(g, pos)
        name = "solution({0})_{1}-{2}".format(self.solution_counter, src, dst)
        
        plt.axis('off')
        plt.legend(loc='best')
        plt.tight_layout(0)
        plt.savefig(name, transparent=True, dpi=300)
        plt.clf()

    def _run(self):
        '''recursive function'''
        # 0. not_done is empty
        if not self.gv.graph['not_done']:
            self.solution_found = True
            #print self.gv.nodes(data=True)

            # 1. for all ev in Ev
            gs = nx.Graph(self.gs)
            gv = nx.Graph(self.gv)
            failed = False
            for src, dst in gv.edges():
                if failed == True: break
                max_split = self.max_split_no
                
                # get the virtual path's required bandwidth
                vbw = gv[src][dst]['bw']                
                gv[src][dst]['allocated_path'] = list()
                done = False
                while not done:
                    if max_split == 0:
                        # FIXME: 不是最优解，最优解的搜索空间是： Cn^m,其中 n = max_k，m = max_split_no
                        print "max_split = 0"
                        failed = True
                        break
                    
                    sp = self._virtual_path_mapping(src, dst, gv, gs)
                    
                    if sp is None: 
                        print "sp is None"
                        failed = True
                        break
                    
                    path, delay, pbw = sp[0]['path'], sp[0]['cost'], sp[0]['bw'] 
                    if self.max_split_no == 1:
                        allowed_path = [(item['path'], item['cost'], item['bw']) for item in sp if item['bw'] >= vbw]
                        if len(allowed_path) == 0: 
                            print "no split and bw is not enough"
                            failed = True
                            break
                        else:
                            (path, delay, pbw)= allowed_path[0]
                    
                    # 对该虚拟路径(src, dst)依次分配物理路径和带宽，并更新资源视图
                    bwmin = min(vbw, pbw)
                    vbw -= bwmin
                    
                    # update gs resource view
                    gv[src][dst]['allocated_path'].append((path, bwmin, delay))                        
                    for link_src, link_dst in zip(path[:-1], path[1:]):
                        if link_src == link_dst: continue
                        gs[link_src][link_dst]['bw'] -= bwmin
                        if gs[link_src][link_dst]['bw'] == 0:
                            gs.remove_edge(link_src, link_dst)

                    if vbw == 0: done = True
                    max_split -= 1
                    
            if failed == True:
                self.solution_found = False
            else: 
                # a solution found, output the result
                self.gs_snapshot = nx.Graph(gs)
                self.solution_counter += 1
                #print '-------------------------------------------------'
                #print "A solution found: "
                # for each virtual path, draw a plot
                for src, dst in gv.edges():
                    self.dump_virtual_path(gs, gv, src, dst, path)
                #print '-------------------------------------------------'                        
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

def remove_images(image_dir):
    for d in [image_dir]:
        for f in os.listdir(d):
            tf = os.path.join(d, f)
            if os.path.splitext(tf)[1] == ".png" and os.path.isfile(tf):
                os.remove(tf)

def build_v_s_from_file(fname):
    """
    build the virtual network and substrate network from the file @fname
    return tuple (gv, gs)
    """
    
    # 1. build substrate network
    C = 3e5 # light speed, in km
    R = 42164.169637
    gs = nx.Graph()
    
    sats = dict()
    with open(fname) as f:
        for line in f:
            line = line[:-1] # skip \n
            infos = string.split(line, ",")
            (satname, loc, cpu) = int(infos[0][1:]), infos[1], infos[2]
            sats[satname] = (satname, float(loc), int(cpu))
            gs.add_node(satname, attr_dict = {'cpu': int(cpu)})
    
    total_nodes_no = gs.number_of_nodes()
    for sat_no in sats.iterkeys():
        pre = (sat_no - 1 + total_nodes_no - 1) % total_nodes_no + 1
        nxt = (sat_no - 1 + total_nodes_no + 1) % total_nodes_no + 1
        oth = ( sat_no - 1 + total_nodes_no + 17) % total_nodes_no + 1
        
        satlist = [pre, nxt, oth]
        
        def dist(n):
            the = abs(sats[n][1] - sats[sat_no][1])
            the = min(360 - the, the)
            
            return 2 * R * math.sin( math.radians(the/2))
        
        s = map(lambda x: dist(x), satlist)
        for sato, dis in zip(satlist, s):
            gs.add_edge(sat_no, sato, {"dy": int(dis)*1e3/C, "bw":100, "dis": dis})
    
    for n in xrange(1, gs.number_of_nodes() + 1):
        satinfos = sats[n]
        t = ""
        for nei in gs.neighbors(n):
            t += ", {0}/{1}".format(nei, int(gs[n][nei]['dis']) )
        #print "{0}, {1}, {2} {3}".format(n, satinfos[1], satinfos[2], t)
    
    pos = nx.circular_layout(gs, scale=2)
    node_labels = {n: "{0}({1})".format(n, gs.node[n]['cpu']) for n in gs.nodes()}
    nx.draw_networkx_nodes(gs, pos, node_sizes=5000)
    nx.draw_networkx_edges(gs, pos)
    nx.draw_networkx_labels(gs, pos, labels=node_labels, font_size=5)
    
    edge_labels = {(u,v): int(bw['dy']) for u, v, bw in gs.edges(data=True)}
    nx.draw_networkx_edge_labels(gs, pos, edge_labels, label_pos=0.3, font_size=8)
    
    #for n in gs.nodes():
        #print n, gs.neighbors(n)
    
    plt.axis('off')       
    plt.savefig("substrate", transparent=True, dpi=300)
    plt.clf()            
    
    # 2. build virtual network no.1
    gv1 = nx.Graph(not_done = deque())
    gv1.add_nodes_from(['CZ1', 'BJ'])
    gv1.graph['not_done'] = gv1.nodes()
    
    # cpu
    map_nv = {'CZ1':[x for x in xrange(26, 46+1)], 
              'BJ':[x for x in (range(1, 23+1) + range(48, 50+1))],
              'BER':[x for x in range(17, 34+1)],
              'BRSL':[x for x in range(24, 43+1)],
              'CPT':[x for x in range(14, 34+1)],
              'LON':[x for x in range(18, 36+1)],
              'MOW':[x for x in range(11, 30+1)],
              'NYC':[x for x in range(27, 45+1)],
              'SIN':[x for x in (range(1, 24+1) + range(49, 50+1))],
              'TYO':[x for x in (range(1, 20+1) + range(49, 50+1))],
              'WAS':[x for x in range(27, 46+1)]}
    map_nv['JC3'] = map_nv['CZ1']
    map_nv['CZ4'] = map_nv['CZ1']
        
    #print map_nv
    for n in gv1.nodes():
        gv1.node[n]['cpu'] = 20
        gv1.node[n]['candidates'] = [x for x in map_nv[n] if gs.node[x]['cpu'] >= gv1.node[n]['cpu']]
        gv1.node[n]['ns'] = None   
         
        
    # bw, dy
    gv1.add_edge('CZ1', 'BJ', {'bw':20})
    
    pos = nx.shell_layout(gv1)
    nx.draw_networkx_nodes(gv1, pos, node_size=1500)
    nx.draw_networkx_edges(gv1, pos)
    
    node_labels = {n: "{0}({1})".format(n, gv1.node[n]['cpu']) for n in gv1.nodes()}    
    nx.draw_networkx_labels(gv1, pos,labels=node_labels)
    
    edge_labels = {(u,v): bw['bw'] for u,v,bw in gv1.edges(data=True)}
    nx.draw_networkx_edge_labels(gv1, pos, edge_labels, label_pos=0.3)
    
    plt.axis('off')
    plt.savefig("virtual1", transparent=True, dpi=300)
    plt.clf()   
        
    # 3. build virtual network no.2
    gv2 = nx.Graph(not_done = deque())
    gv2.add_nodes_from(['JC3','CZ4', 'BJ'])
    gv2.graph['not_done'] = gv2.nodes()
    
    # cpu
    for n in gv2.nodes():
        gv2.node[n]['cpu'] = 20
        gv2.node[n]['candidates'] =  [x for x in map_nv[n] if gs.node[x]['cpu'] >= gv2.node[n]['cpu']]
        gv2.node[n]['ns'] = None 
        
    # bw, dy
    gv2.add_edge('JC3', 'BJ', {'bw':20})
    gv2.add_edge('CZ4', 'BJ', {'bw':20})    

    pos = nx.shell_layout(gv2)
    nx.draw_networkx_nodes(gv2, pos, node_size=1500)
    nx.draw_networkx_edges(gv2, pos)
    
    node_labels = {n: "{0}({1})".format(n, gv2.node[n]['cpu']) for n in gv2.nodes()}    
    nx.draw_networkx_labels(gv2, pos,labels=node_labels)
    
    edge_labels = {(u,v): bw['bw'] for u,v,bw in gv2.edges(data=True)}
    nx.draw_networkx_edge_labels(gv2, pos, edge_labels, label_pos=0.3)
    
    plt.axis('off')
    plt.savefig("virtual2", transparent=True, dpi=300)
    plt.clf() 
    
    
    # 4. build virtual network no. 3
    gv3 = nx.Graph(not_done = deque())
    gv3.add_nodes_from(['BJ','BER', 'BRSL', 'CPT', 'LON', 'MOW', 'NYC', 'SIN', 'TYO', 'WAS'])
    gv3.graph['not_done'] = gv3.nodes()   
    
    for n in gv3.nodes():
        if n == 'BJ':
            gv3.node[n]['cpu'] = 180
        else:
            gv3.node[n]['cpu'] = 20
        gv3.node[n]['candidates'] =  [x for x in map_nv[n] if gs.node[x]['cpu'] >= gv3.node[n]['cpu']]
        gv3.node[n]['ns'] = None 
    
    gv3.node['BJ']['cpu'] = 180    
    
    # bw, dy
    gv3.add_edge('BER', 'BJ', {'bw':20})
    gv3.add_edge('BRSL', 'BJ', {'bw':20})    
    gv3.add_edge('CPT', 'BJ', {'bw':20})    
    gv3.add_edge('LON', 'BJ', {'bw':20})    
    gv3.add_edge('MOW', 'BJ', {'bw':20})    
    gv3.add_edge('NYC', 'BJ', {'bw':20})    
    gv3.add_edge('SIN', 'BJ', {'bw':20})    
    gv3.add_edge('TYO', 'BJ', {'bw':20})    
    gv3.add_edge('WAS', 'BJ', {'bw':20})        

    pos = nx.shell_layout(gv3)
    nx.draw_networkx_nodes(gv3, pos, node_size=1500)
    nx.draw_networkx_edges(gv3, pos)
    
    node_labels = {n: "{0}({1})".format(n, gv3.node[n]['cpu']) for n in gv3.nodes()}    
    nx.draw_networkx_labels(gv3, pos,labels=node_labels)
    
    edge_labels = {(u,v): bw['bw'] for u,v,bw in gv3.edges(data=True)}
    nx.draw_networkx_edge_labels(gv3, pos, edge_labels, label_pos=0.3)
    
    plt.axis('off')
    plt.savefig("virtual3", transparent=True, dpi=300)
    plt.clf() 
        
    return (gv1, gv2, gv3, gs)

def main():
    remove_images("./")
    gv1, gv2, gv3, gs = build_v_s_from_file("sat.csv")
    #Gv = build_virtual_network()    
    #Gs = build_substrate_network(5, 2)
    algo = BaselineAlgo(gv2, gs, all_solutions = False, max_split_no = 1, max_k = 5)
    algo.run()
    #print algo.get_gs_snapshot().nodes(True)
    #print algo.get_gs_snapshot().edges(data=True)
    BaselineAlgo(gv3, algo.get_gs_snapshot(), all_solutions = False, max_split_no = 1, max_k = 5).run()    

if __name__ == '__main__':
    main()
