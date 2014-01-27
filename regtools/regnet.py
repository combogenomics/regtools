#!/usr/bin/python

import numpy as np
import networkx as nx
import itertools

regstates = ['plugged', 'unplugged',
             'ready', 'notready', 'absent',
             'missing']
             
def getRegStatesCombinations():
    for a, b in itertools.combinations(sorted(regstates), 2):
        yield a, b
        
def getRegStatesPermutations():
    for a, b in itertools.permutations(sorted(regstates), 2):
        yield a, b

def getRegState(org, reg, prom, gene):
    o = org
    
    if o in reg and o in prom and o in gene:
        return 'plugged'
    elif o in reg and o not in prom and o in gene:
        return 'unplugged'
    elif o not in reg and o in prom and o in gene:
        return 'ready'
    elif o not in reg and o not in prom and o in gene:
        return 'notready'
    elif o in reg and o not in prom and o not in gene:
        return 'absent'
    elif o not in reg and o not in prom and o not in gene:
        return 'missing'
    else:
        raise ValueError('Unattended case! (%s // %s - %s - %s)'%(org,
                        bool(org in reg), bool(org in prom),
                        bool(org in gene)))

class RegLink(object):
    header = '\t'.join( ('Regulator', 'Gene', '# orgs',
                        'Plugged', 'Unplugged', 'Ready',
                        'Not ready', 'Absent', 'Missing', 'Plug length') )

    def __init__(self, source, target, norgs):
        self.regulator = source
        self.gene = target
        self.norgs = norgs
        self.plugged = 0.
        self.pluglen = 0.
        self.unplugged = 0.
        self.ready = 0.
        self.notready = 0.
        self.absent = 0.
        self.missing = 0.
        
    def _sumAll(self):
        return sum( (self.plugged,
                   self.unplugged,
                   self.ready,
                   self.notready,
                   self.absent,
                   self.missing) )
        
    def __str__(self):
        return '\t'.join( [str(x)
                           for x in ( self.regulator, self.gene,
                                       self.norgs, round(self.plugged,2),
                                       round(self.unplugged,2),
                                       round(self.ready,2),
                                       round(self.notready,2),
                                       round(self.absent,2),
                                       round(self.missing,2),
                                       round(self.pluglen,2))])

def getMean(objs, attr):
    return np.array( [getattr(x, attr) for x in objs] ).mean()
    
def getMad(objs, attr):
    data = np.array([getattr(x, attr) for x in objs])
    return np.median( np.absolute(data - np.median(data)))

def getPlug(net, node, org):
    # Extract the DFS tree
    tree = nx.depth_first_search.dfs_tree(net, node)
    
    # Prune the DFS tree
    for a, b in tree.edges():
        orgs = net[a][b]['orgs'].split()
        
        if org not in orgs:
            del tree[a][b]
    
    # Do another DFS on the pruned tree
    tree = nx.depth_first_search.dfs_tree(tree, node)
    
    return tree.edges()
    
def getDownstream(net, node, org):
    # Extract the DFS tree
    tree = nx.depth_first_search.dfs_tree(net, node)
    
    # Prune the DFS tree
    for a, b in tree.edges():
        orgs = net[a][b]['orgs'].split()
        
        if org not in orgs:
            del tree[a][b]
    
    # Do another DFS on the pruned tree
    tree = nx.depth_first_search.dfs_tree(tree, node)
    
    return tree.nodes()

def getPlugLen(net, node, org):
    return len(getPlug(net, node, org))
