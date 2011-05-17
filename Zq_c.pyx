# -*- coding: utf-8 -*-
"""
Zero forcing (Z_q)

This module implements zero forcing (Z_q version) using fast bitsets
and a brute-force approach to trying various bitsets. 
"""

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


include "sage/misc/bitset_pxd.pxi"
include "sage/misc/bitset.pxi"
from sage.misc.bitset cimport FrozenBitset, Bitset    
    
cpdef push_zeros(list neighbors, FrozenBitset subgraph, FrozenBitset filled_set, bint return_bitset=True):
    """
    Run zero forcing as much as possible
    
    INPUT
    :param: neighbors (list of FrozenBitsets) -- the neighbors of each vertex
    :param: subgraph (FrozenBitset) -- the subgraph we are forcing on
    :param: filled_set (FrozenBitset) -- the initial filled vertices
    :param: return_bitset (bool) -- if True, return the set of filled vertices after playing the game,
        if False, return a boolean can_push, where can_push is True iff a force could happen.
        
        The returned filled set contains all of the initially filled vertices, even if they are outside
        of the subgraph.
    """
    cdef bitset_s *filled = filled_set._bitset
    cdef bitset_s *subgraph_bitset = subgraph._bitset

    cdef FrozenBitset v
    cdef bitset_s *current_neighbors
    cdef bitset_t unfilled_neighbors
    cdef bitset_t filled_active
    cdef bitset_t filled_active_copy
    cdef bitset_t unfilled # unfilled in the subgraph
    
    # Initialize filled_active to the complement of unfilled
    bitset_init(filled_active, filled.size)
    bitset_copy(filled_active, filled)
    
    bitset_init(unfilled, filled.size)
    bitset_complement(unfilled, filled)
    bitset_intersection(unfilled, unfilled, subgraph_bitset)

    bitset_init(filled_active_copy, filled.size)
    bitset_init(unfilled_neighbors, filled.size)
    
    cdef FrozenBitset ret = FrozenBitset(None, capacity=filled.size)
    cdef bitset_s *ret_bitset=ret._bitset

    cdef int new_filled, n
    cdef bint done = False
    cdef bint can_push = False
    
    while not done:
        done = True
        bitset_copy(filled_active_copy, filled_active)
        n = bitset_first(filled_active_copy)
        while n>=0:
            v = neighbors[n]
            current_neighbors = v._bitset
            bitset_intersection(unfilled_neighbors, current_neighbors, unfilled)
            new_filled = bitset_first(unfilled_neighbors)
            if new_filled < 0:
                # no unfilled neighbors
                bitset_discard(filled_active, n)
            else:
                # look for second unfilled neighbor
                if bitset_next(unfilled_neighbors, new_filled+1) < 0:
                    # No more unfilled neighbors
                    # push to the new_filled vertex
                    bitset_add(filled_active, new_filled)
                    bitset_remove(unfilled, new_filled)
                    bitset_remove(filled_active, n)
                    if return_bitset:
                        done = False
                    else:
                        can_push=True
                        # done is True, so the break leads to breaking out of both while loops
                        break
            n = bitset_next(filled_active_copy, n+1)

    if return_bitset:
        bitset_complement(ret_bitset, unfilled)
        bitset_intersection(ret_bitset, ret_bitset, subgraph_bitset)
        bitset_union(ret_bitset, ret_bitset, filled)

    # Free all memory used:
    bitset_free(filled_active)
    bitset_free(filled_active_copy)
    bitset_free(unfilled_neighbors)
    bitset_free(unfilled)
    
    if return_bitset:
        return ret
    else:
        return can_push
    

cpdef neighbors_connected_components(list neighbors, FrozenBitset subgraph):
    cdef int n=len(neighbors)
    cdef bitset_s *subneighbors = <bitset_s *> sage_malloc(n*sizeof(bitset_s))
    cdef int i, visit
    cdef bitset_s *b
    cdef FrozenBitset FB
    cdef bitset_s *subgraph_set = subgraph._bitset
    for i in range(n):
        FB=neighbors[i]
        b=FB._bitset        
        bitset_init(&subneighbors[i], n)
        bitset_intersection(&subneighbors[i], b, subgraph_set)
            
    components=set()
    cdef bitset_t seen, queue, component
    bitset_init(seen, n)
    bitset_init(queue, n)
    bitset_init(component, n)

    for i in range(n):
        if bitset_in(subgraph_set,i) and bitset_not_in(seen, i):
            bitset_clear(queue)
            bitset_clear(component)
            
            # do breadth/depth-first search from i
            bitset_add(queue,i)
            visit=i
            while visit>=0:
                bitset_remove(queue, visit)
                if bitset_not_in(component, visit):
                    bitset_add(component, visit)
                    bitset_union(queue, queue, &subneighbors[visit])
                visit=bitset_first(queue)

            components.add(tuple(bitset_list(component)))
            bitset_union(seen, seen, component)
            
            
    for i in range(n):
        bitset_free(&subneighbors[i])
    sage_free(subneighbors)

    bitset_free(seen)
    bitset_free(queue)
    bitset_free(component)
    return components
