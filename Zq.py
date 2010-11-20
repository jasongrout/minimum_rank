# -*- coding: utf-8 -*-
"""
Zero forcing

This module implements zero forcing using fast bitsets and a
brute-force approach to trying various bitsets.
"""


#######################################################################
#
# Copyright (C) 2008 Steve Butler, Jason Grout.
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.
#######################################################################


try:
   from Zq_c import can_push, neighbors_connected_components
except ImportError:
   # assume everything is in the global space
   # this happens when, for example, someone "load"s the right files
   # in the Sage notebook
   pass

from itertools import combinations,chain
def subsets(s,r=None):
    """ Returns subsets of size r, or if r is None, return nonempty subsets of any size"""
    if r is not None:
        return combinations(s,r)
    else:
        s = list(s)
        return chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))


from sage.all import Bitset,FrozenBitset
def zero_forcing_sets(G=None,neighbors=None):
    """
    Calculate all zero forcing sets

    Returns the zero forcing number, a minimal zero forcing set, and a
    set of all zero forcing sets
    """
    if neighbors is None:
        n=G.order()
        neighbors=[Bitset(G.neighbors(i),capacity=n) for i in range(n)]
    n=len(neighbors)
    subgraph=Bitset(range(n),capacity=n)
    V=set(range(n))
    
       
    L=set(map(FrozenBitset,subsets(V,n)))
    L|=set(map(FrozenBitset,subsets(V,n-1)))
    zero_forcing_number=n-1
    #print L
    lastZ=set(list(V)[:-1])
    debug=False
    for sizeZ in range(n-2,-1,-1):
        #if debug: print "Exploring size",sizeZ
        for Z in subsets(V, sizeZ):
            Z=FrozenBitset(Z,capacity=n)
            
            pushed, W=can_push(neighbors, subgraph, Z)
            #if debug: print "W",W
            if pushed and W in L:
                #if debug: print "W in L"
                L.add(Z)
                lastZ=Z
                zero_forcing_number=sizeZ
        if zero_forcing_number>sizeZ: # we are done
            return zero_forcing_number, lastZ, L



def Z_pythonBitset(G,q,zfs_sets=None):
    """
    Assumes the graph has vertices 0,1,2,...,n-1
    """
    G=G.copy()

    if zfs_sets is None:
        relabel=G.relabel(return_map=True)
        reverse_map=dict( (v,k) for k,v in relabel.items())
    elif G.vertices()==range(G.order()):
        reverse_map=dict( (i,i) for i in range(G.order()))
    else:
        raise ValueError("If zfs_sets is supplied, G must have vertices 0,1,...,n-1, corresponding to the vertices in zfs_sets")

    n=G.order()
    V=FrozenBitset(G.vertices(),capacity=n)
    neighbors=[FrozenBitset(G.neighbors(i),capacity=n) for i in range(n)]
    if not G.is_connected():
        raise ValueError("G needs to be connected")
    if n<2:
        raise ValueError("G needs to have 2 or more vertices")
       
    if zfs_sets is None: 
       zfs_sets=zero_forcing_sets(neighbors=neighbors) 
    else: 
       zfs_sets=zfs_sets[0], zfs_sets[1], copy(zfs_sets[2]) 
    zero_forcing_number, lastZ,L=zfs_sets

    #print L
    debug=False
    for sizeZ in range(n-2,0,-1):
        #if debug: print "Exploring size",sizeZ
        for Z in subsets(V, sizeZ):
            Z=FrozenBitset(Z,capacity=n)
            if Z not in L:
                #H=[tuple(i) for i in G.subgraph(V.difference(Z)).connected_components()]
                #if debug: print "calling H",neighbors, V.difference(Z)
                H=neighbors_connected_components(neighbors, V.difference(Z))
                #if debug: print "H",H
                #if debug: print "looping through subsets:",list(subsets(H,q+1))
                for J in subsets(H,q+1):
                    #if debug: print "Hand to opponent: ",J
                    for K in subsets(J):
                        if len(K)==0:
                            continue # ignore the empty set
                        new_vertices=Bitset(Z,capacity=n)
                        for s in K:
                            new_vertices.update(FrozenBitset(s))
                        #new_vertices=FrozenBitset(Z.union(*K),capacity=n)
                        #if debug: print "Opponent hands back: ",K, "so we have vertices",new_vertices
                        # we don't have to do the empty subset, but it is simpler to just include it
                        
                        pushed, X=can_push(neighbors, new_vertices, Z)
                        #if debug: print "X, ZFS on opponent's return: ",X
                        if X not in L:
                            #if debug: print "X not in L"
                            break # go to next J
                        #if debug: print "X is in L"
                    else:
                        L.add(Z)
                        if sizeZ<zero_forcing_number:
                            lastZ=Z
                            zero_forcing_number=sizeZ
                        break # next Z
        #print zero_forcing_number, sizeZ
        if zero_forcing_number>sizeZ: # we are done
            return zero_forcing_number, set(reverse_map[i] for i in lastZ)



def Zq_inertia_lower_bound(G):
    """
    Run the Zq iteratively until we get a good lower bound.
    (n-q-Z(G,q)-1,q) and (q,n-q-Z(G,q)-1)is not in the inertia set of G
    By the Southwest lemma, everything south and west of the point is not in the inertia
    
    G is assumed to be connected
    """
    G=G.relabel(inplace=False)
    n=G.order()
    not_in_inertia=set()
    zfs_sets=zero_forcing_sets(G)
    zero_forcing_number=zfs_sets[0]
    compute_Zq=True
    for q in range(ceil(n/2)):
        print "calculating Z%s"%q
        if compute_Zq is True:
            Zq=Z_pythonBitset(G,q,zfs_sets=zfs_sets)[0]
        else:
            Zq=zero_forcing_number
        if Zq==zero_forcing_number:
            compute_Zq=False
        if n-q-Zq-1==q:
            not_in_inertia.update([(n-q-Zq-1,q),(q,n-q-Zq-1)])
        if n-q-Zq-1<=q:
            break
        not_in_inertia.update([(n-q-Zq-1,q),(q,n-q-Zq-1)])
    # Find the closure to give a point in each row and each column
    sorted_boundary=sorted(not_in_inertia,reverse=True)
    curr_x,curr_y=sorted_boundary[0]    
    for pt in sorted_boundary[1:]:
        new_x,new_y=pt
        if new_x<curr_x-1:
            not_in_inertia.update([(i,curr_y) for i in range(new_x+1,curr_x)])
            not_in_inertia.update([(curr_y,i) for i in range(new_x+1,curr_x)])
        curr_x,curr_y=new_x,new_y
    return not_in_inertia

def Zplus(G):
   return Z(G,q=0)

from sage.all import Graph, graphs
G=Graph()
G.add_edges([[1,2],[2,3],[3,4],[4,5],[5,6],[6,1],[1,4],[2,5],[3,6],[7,1],[7,2],[7,3]])

G2=graphs.CompleteGraph(4)
G2.subdivide_edges(G2.edges(),1)

from sage.all import points
def plot_inertia_lower_bound(g):
   return points(list(Zq_inertia_lower_bound(g)),
                pointsize=40,gridlines=True,
                ticks=[range(g.order()),range(g.order())],
                aspect_ratio=1)

"""
import cProfile as cp
cp.run('Z_pythonBitset(graphs.HeawoodGraph(),q=2)',sort='time')
"""

###################################################################
##   OLD CODE SUPERSEDED BY THE ABOVE CODE ########################
###################################################################

old_code='''

def playzerosgame(graph, initial_set=[], one_step=False):
   """
   Apply the color-change rule to a given graph given an optional
   initial set.

   :param graph: the graph on which to apply the rule
   :param initial_set: the set of "zero" (black) vertices in the graph
   
   :return: the list of zero (black) vertices in the resulting derived
        coloring

   EXAMPLES:: 

        sage: from sage.graphs.minrank import zerosgame
        sage: zerosgame(graphs.PathGraph(5))
        []
        sage: zerosgame(graphs.PathGraph(5),[0])
        [0, 1, 2, 3, 4]
   """
   new_zero_set=set(initial_set)
   zero_set=set([])
   zero_neighbors={}
   active_zero_set = set([])
   inactive_zero_set = set([])
   another_run=True
   while another_run:
       another_run=False
       # Add the new zero vertices
       zero_set.update(new_zero_set)
       active_zero_set.update(new_zero_set)
       active_zero_set.difference_update(inactive_zero_set)
       zero_neighbors.update([[i, 
                   set(graph.neighbors(i)).difference(zero_set)] 
                              for i in new_zero_set])
       # Find the next set of zero vertices
       new_zero_set.clear()
       inactive_zero_set.clear()
       for v in active_zero_set:
           zero_neighbors[v].difference_update(zero_set)
           if len(zero_neighbors[v])==1:
               new_zero_set.add(zero_neighbors[v].pop())
               inactive_zero_set.add(v)
               another_run=True
               if one_step is True:
                   return zero_set.union(new_zero_set)
       if one_step is True and another_run is False:
           return zero_set
   return list(zero_set)

def Z_python(G,q,zfs_sets=None):
    n=G.order()
    V=set(G.vertices())
    if not G.is_connected():
        raise ValueError("G needs to be connected")
    if n<2:
        raise ValueError("G needs to have 2 or more vertices")
       
    if zfs_sets is None:
        zfs_sets=zero_forcing_sets(G)
    zero_forcing_number, lastZ,L=zfs_sets

    #print L
    debug=False
    for sizeZ in range(n-2,0,-1):
        #if debug: print "Exploring size",sizeZ
        for Z in subsets(V, sizeZ):
            Z=frozenset(Z)
            if Z not in L:
                H=[tuple(i) for i in G.subgraph(V.difference(Z)).connected_components()]
                #if debug: print "H",H
                #if debug: print "looping through subsets:",list(subsets(H,q+1))
                for J in subsets(H,q+1):
                    #if debug: print "Hand to opponent: ",J
                    for K in subsets(J):
                        if len(K)==0:
                            continue # ignore the empty set
                        new_vertices=Z.union(*K)
                        #if debug: print "Opponent hands back: ",K, "so we have vertices",new_vertices
                        # we don't have to do the empty subset, but it is simpler to just include it
                        X=set(playzerosgame(G.subgraph(new_vertices),initial_set=Z))
                        #if debug: print "X, ZFS on opponent's return: ",X
                        if X not in L:
                            #if debug: print "X not in L"
                            break # go to next J
                        #if debug: print "X is in L"
                    else:
                        L.add(Z)
                        if sizeZ<zero_forcing_number:
                            lastZ=Z
                            zero_forcing_number=sizeZ
                        break # next Z
        #print zero_forcing_number, sizeZ
        if zero_forcing_number>sizeZ: # we are done
            return zero_forcing_number, lastZ

from sage.all import Subsets, Set
def Z_sage(G,q):

    n=G.order()
    V=Set(G.vertices())
    if not G.is_connected():
        raise ValueError("G needs to be connected")
    if n<2:
        raise ValueError("G needs to have 2 or more vertices")
       
    L=set(Subsets(V,n))
    L|=set(Subsets(V,n-1))
    lastZ=V
    debug=False
    for sizeZ in range(n-2,0,-1):
        #print "Exploring size",sizeZ
        for Z in Subsets(V, sizeZ):
            W=Set(playzerosgame(G,initial_set=Z)) # find right command
            #if debug: print "W",W
            if W in L:
                #if debug: print "W in L"
                L.add(Z)
                lastZ=Z
                zero_forcing_number=sizeZ
            else:
                #if debug: print "W not in L"
                H=[tuple(i) for i in G.subgraph(V-Z).connected_components()]
                #if debug: print "H",H
                #if debug: print "looping through subsets:",Subsets(H,q+1).list()
                for J in Subsets(H,q+1):
                    #if debug: print "Hand to opponent: ",J
                    for K in Subsets(J):
                        if len(K)==0:
                            continue # ignore the empty set
                        #if debug: print "Opponent hands back: ",K, "so we have vertices",Z+sum([Set(i) for i in K],Set([]))
                        # we don't have to do the empty subset, but it is simpler to just include it
                        X=Set(playzerosgame(G.subgraph(Z+Set(set([]).union(*K))),initial_set=Z))
                        #if debug: print "X, ZFS on opponent's return: ",X
                        if X not in L:
                            #if debug: print "X not in L"
                            break # go to next J
                        #if debug: print "X is in L"
                    else:
                        L.add(Z)
                        lastZ=Z
                        zero_forcing_number=sizeZ
                        break # next Z
        #print zero_forcing_number, sizeZ
        if zero_forcing_number>sizeZ: # we are done
            return zero_forcing_number, lastZ
'''
