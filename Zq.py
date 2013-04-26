# -*- coding: utf-8 -*-

#######################################################################
#
# Copyright (C) 2011 Steve Butler, Jason Grout.
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
    from Zq_c import push_zeros, push_zeros_looped, neighbors_connected_components
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


def Z_pythonBitsetold(G,q,zfs_sets=None):
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
       zfs_sets=zfs_sets[0], zfs_sets[1], zfs_sets[2].copy() 
    zero_forcing_number, lastZ,L=zfs_sets

    #print L
    debug=False
    # TODO: Should this be range(n-2, -1, -1) (so that sizeZ could be 1!)
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




def Zq_inertia_lower_bound(G,verbose=False):
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
    for q in range(n//2+1): # ceil(n/2)
        if verbose: print "calculating Z%s"%q
        if compute_Zq is True:
            Zq=Zq_bitset(G,q,zfs_sets=zfs_sets)[0]
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


# if in library mode, we need to import this.
# if it is just included in a sage notebook, then it is in the global namespace
# so we try importing it if we can.  If we can't import it, then just trust that
# everything is in the global namespace.
try:
    from zero_forcing_wavefront import zero_forcing_set_wavefront
    from inertia import InertiaSet
    from Zq_c import push_zeros
except ImportError:
    pass


def Zq_inertia_lower_bound(G, zero_forcing_function=None, verbose=False):
    """
    Run the Zq iteratively until we get a good lower bound.
    (n-q-Z(G,q)-1,q) and (q,n-q-Z(G,q)-1)is not in the inertia set of G
    By the Southwest lemma, everything south and west of the point is not in the inertia
    
    G is assumed to be connected
    """
    global Zq_compute
    if zero_forcing_function is None:
        zero_forcing_function=Zq_compute
    G=G.relabel(inplace=False)
    n=G.order()
    I = InertiaSet([(G.order(), G.order())])
    zero_forcing_number=zero_forcing_function(G,n)
    compute_Zq=True
    for q in range(n//2+1): # ceil(n/2)
        if verbose: print "calculating Z%s"%q
        if compute_Zq is True:
            Zq=zero_forcing_function(G,q)
        else:
            Zq=zero_forcing_number
        if Zq==zero_forcing_number:
            compute_Zq=False
        if n-q-Zq==q:
            I|=[(n-q-Zq,q)]
        if n-q-Zq<=q:
            break
        I|=[(n-q-Zq,q)]
    return I



from sage.all import points
def plot_inertia_lower_bound(g, Zq_args={}):
    n=g.order()
    return plot(Zq_inertia_lower_bound(g, **Zq_args),
                pointsize=40,gridlines=True,
                ticks=[range(n+1),range(n+1)],
                aspect_ratio=1)+line([(0,n),(n,0)],linestyle=':')



#################################################################
######  Better Algorithm
#################################################################


def Zq_graph_info(G):
    """
    Extract (and cache) necessary graph information for the Zq_bitset function.

    We've separated this out so that this busy work can be easily cached.
    """
    G=G.copy()

    relabel=G.relabel(return_map=True)
    reverse_map=dict( (v,k) for k,v in relabel.items())
    R=lambda x: reverse_map[x]

    n=G.order()
    V=FrozenBitset(G.vertices(),capacity=n)
    neighbors=[FrozenBitset(G.neighbors(i),capacity=n) for i in range(n)]
    return reverse_map, R, n, V, neighbors

def Zq_bitset(G,q, push_zeros, push_zeros_kwargs=dict(), return_track=False):
    """
    Calculate Zq, where you can have arbitrary color rules encoded in a push_zeros function.

    INPUT:
    :param G: a simple undirected graph
    :param q: the :math:`q` for the algorithm
    :param push_zeros: a function with the signature ``push_zeros(neighbors, subgraph, filled_set, return_bitset=False, **kwargs)``, where the extra ``kwargs`` are another parameter.
    :param push_zeros_kwargs: extra arguments to the push_zeros function
    :param return_track: (bool) whether to return a sequence of actions that obtain the Zq value.
    """
    # We aggressively cache the graph information
    if not isinstance(G, tuple):
        G = Zq_graph_info(G)
    reverse_map, R, n, V, neighbors = G
    # TODO: Why is this important?
    if n<2:
        raise ValueError("G needs to have 2 or more vertices")
    cost={V: 0}
    if return_track: 
        track={V:'done'}
    
    #print L
    debug=True
    for sizeZ in range(n-1,-1,-1):
        for Z in subsets(V, sizeZ):
            Z=FrozenBitset(Z,capacity=n)
            if push_zeros(neighbors, subgraph=V, filled_set=Z, return_bitset=False,
                          **push_zeros_kwargs):
                #print "can push, so skipping", Z
                continue
            b=n
            c=n
            cost[Z]=n
            if return_track:
                track[Z]='sentinal'
            H=neighbors_connected_components(neighbors, V.difference(Z))
            #if debug: print "H",H
            #if debug: print "looping through subsets:",list(subsets(H,q+1))
            for J in subsets(H,q+1):
                #if debug: print "Hand to opponent: ",J
                bb=-1
                for K in subsets(J):
                    if len(K)==0:
                        continue # ignore the empty set
                        
                    subgraph=Bitset(Z,capacity=n)
                    for s in K:
                        subgraph.update(FrozenBitset(s))
                    closed_Z=push_zeros(neighbors, subgraph=subgraph, 
                                filled_set=Z, return_bitset=True, **push_zeros_kwargs)
                    closed_Z=push_zeros(neighbors, subgraph=V, 
                                filled_set=closed_Z, return_bitset=True, **push_zeros_kwargs)
                    if cost[closed_Z]>bb:
                        bb=cost[closed_Z] #max(bb,cost[closed_Z])
                        if return_track:
                            bb_set=(J,K,closed_Z)
                if bb<b:
                    b=bb #min(b,bb)
                    if return_track:
                        b_set=bb_set
            for v in V-Z:
                closed_Z=Z.union(FrozenBitset([v],capacity=n))
                closed_Z=push_zeros(neighbors, subgraph=V, filled_set=closed_Z,
                                    return_bitset=True, **push_zeros_kwargs)
                if cost[closed_Z]+1<c:
                    c=cost[closed_Z]+1 #min(c, cost[closed_Z]+1)
                    if return_track:
                        c_vertex=v
                        c_closed=closed_Z
            if b<c:
                cost[Z]=b #min(b,c)
                if return_track:
                    track[Z]=('set: hand %s to adversary; adversary hands back %s, push to get %s'%([map(R,i) for i in b_set[0]], [map(R,i) for i in b_set[1]], map(R,b_set[2])), b_set[2])
            else:
                cost[Z]=c #min(b,c)
                if return_track:
                    track[Z]=('spend vertex %s, get %s'%(reverse_map[c_vertex],map(R,c_closed)), c_closed)
    if return_track:
        trail=''
        Z=FrozenBitset([], capacity=n)
        while Z!=FrozenBitset(V,capacity=n):
            trail+='%s\n'%(track[Z],)
            Z=track[Z][1]
        return cost[push_zeros(neighbors, subgraph=V, 
                               filled_set=FrozenBitset([], capacity=n), 
                               return_bitset=True, **push_zeros_kwargs)], trail
    else:
        return cost[push_zeros(neighbors, subgraph=V, 
                               filled_set=FrozenBitset([], capacity=n), 
                               return_bitset=True, **push_zeros_kwargs)]



def Zqhat_recurse(G,q,looped,unlooped, BEST_LOWER_BOUND, BEST_LOOPS, CACHE, G_info):
    """
    We construct a tree of possibilities of looping and unlooping
    vertices, where the root of the tree is all vertices unspecified,
    and the leaves are all possibilities of specifying loops and
    unloops.  Each node has two children, corresponding to making a
    fixed vertex looped or unlooped.

    This function returns the maximum Zqhat of all leaves.  We know
    that Zqhat decreases as we go down from the root (a looped/unlooped vertex has more
    options to force than an unspecified vertex), so if we ever
    find a leaf that is the same value as a node in the tree, we don't
    have to explore that subtree any more.

    BEST_LOWER_BOUND holds the global maximum lower bound for Zqhat.  If a node's Zq
    is lower than this, then we have no hope of increasing the best lower bound, so we
    shortcut return.

    BEST_LOOPS is a list (in which case we will fill it with the loopsets that give us the lower bound)
    or it is False, in which case we can shortcut operations even more to run faster.

    """
    n=G.order()
    unmarked = Bitset(range(n))-looped-unlooped

    canonical_label = (G.canonical_label(partition=[list(looped), list(unlooped), list(unmarked)]).graph6_string(), len(looped), len(unlooped))
    if canonical_label in CACHE:
        # already investigated this node (up to permutation respecting loops)
        #print "skipping looped: %s, unlooped %s"%(looped, unlooped)
        return
    else:
        CACHE.add(canonical_label)

    # Zq = how many vertices is Black forced to use
    Zq=Zq_bitset(G_info,q,push_zeros=push_zeros_looped, 
                 push_zeros_kwargs=dict(looped=looped,unlooped=unlooped))

    if Zq<BEST_LOWER_BOUND[0]:
        # Some leaf in another branch of the tree is already doing better than this entire branch
        # could possibly do, so prune here.
        return
    elif BEST_LOOPS is False and Zq==BEST_LOWER_BOUND[0]:
        # if we don't care about calculating the best loops, then we can shortcut even equality
        # since we only care about improving the bound
        return

    if len(unmarked)==0:
        # we are at a leaf and ready to evaluate
        # Furthermore, because of the above, we know we have a better (or equal if BEST_LOOPS is a list)
        # BEST_LOWER_BOUND than we have found so far, so store it.
        if BEST_LOOPS is not False:
            if Zq == BEST_LOWER_BOUND[0]:
                BEST_LOOPS.append(dict(looped=looped, unlooped=unlooped))
            else:
                BEST_LOWER_BOUND[0]=Zq
                BEST_LOOPS[:] = [dict(looped=looped, unlooped=unlooped)]
        else:
            # We must have a better bound
            BEST_LOWER_BOUND[0]=Zq
        return
    else:
        # we need to choose an unmarked vertex to go further down the branches
        v=unmarked.pop()

        new_looped=looped.union(FrozenBitset([v],capacity=n))
        Zqhat_recurse(G,q,looped=new_looped, unlooped=unlooped,
                      BEST_LOWER_BOUND=BEST_LOWER_BOUND, BEST_LOOPS=BEST_LOOPS, CACHE=CACHE, G_info=G_info)


        if BEST_LOOPS is False and Zq<=BEST_LOWER_BOUND[0]:
            # We've done as good as we possibly can on this branch,
            # and we don't care about exploring everything, so shortcut and return
            return

        new_unlooped=unlooped.union(FrozenBitset([v],capacity=n))
        Zqhat_recurse(G,q,looped=looped,unlooped=new_unlooped,
                      BEST_LOWER_BOUND=BEST_LOWER_BOUND, BEST_LOOPS=BEST_LOOPS, CACHE=CACHE, G_info=G_info)
        return

def Zqhat(G, q, return_loops=False):
    n=G.order()
    full_set=FrozenBitset(range(n))
    empty_set=FrozenBitset([],capacity=n)
    # calculate Zqhat for a few loopsets to get a trivial lower bound for Zqhat
    # BEST_* are mutable lists so that the recursive calls can change these
    # variables.  In a sense, this is a global variable for the Zqhat
    # recursive calls.
    BEST_LOWER_BOUND = [0]
    if return_loops:
	BEST_LOOPS = []
    else:
        BEST_LOOPS = False
    G_info = Zq_graph_info(G)
    for loopset in [dict(looped=full_set,unlooped=empty_set), dict(looped=empty_set,unlooped=full_set)]:
        Zq = Zq_bitset(G_info,q,push_zeros=push_zeros_looped,
                       push_zeros_kwargs=loopset)
        if Zq < BEST_LOWER_BOUND[0]:
            BEST_LOWER_BOUND[0] = Zq
            if BEST_LOOPS is not False:
                BEST_LOOPS[:] = [loopset]
        elif BEST_LOOPS is not False and Zq == BEST_LOWER_BOUND[0]:
            BEST_LOOPS.append(loopset)

    Zqhat_recurse(G, q, FrozenBitset([], capacity=n),
                  FrozenBitset([], capacity=n), BEST_LOWER_BOUND=BEST_LOWER_BOUND, BEST_LOOPS=BEST_LOOPS, CACHE=set(), G_info=G_info)
    if return_loops:
        return BEST_LOWER_BOUND[0], BEST_LOOPS
    else:
        return BEST_LOWER_BOUND[0]

def Zq_compute(G,q):
    return Zq_bitset(G,q,push_zeros=push_zeros)

def Zplus(G):
   return Zq_compute(G,0)

from sage.all import Graph, graphs
G=Graph()
G.add_edges([[1,2],[2,3],[3,4],[4,5],[5,6],[6,1],[1,4],[2,5],[3,6],[7,1],[7,2],[7,3]])

G2=graphs.CompleteGraph(4)
G2.subdivide_edges(G2.edges(),1)


def check_trees(start,end):
    for i in range(start,end):
        print "working on %s vertices"%i
        list(check_tree(list(graphs.trees(i))))

@parallel            
def check_tree(g):
    if not inertia_set(g,f)==Zq_inertia_lower_bound(g):
        if not inertia_set(g,f)==Zq_inertia_lower_bound(g, zero_forcing_function=Zqhat):
            print '\n%s'%g.graph6_string()
    return None
    

"""
import cProfile as cp
cp.run('Z_pythonBitset(graphs.HeawoodGraph(),q=2)',sort='time')
"""
