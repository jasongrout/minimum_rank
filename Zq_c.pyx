# -*- coding: utf-8 -*-
"""
Zero forcing (Z_q)

This module implements zero forcing (Z_q version) using fast bitsets
and a brute-force approach to trying various bitsets. 
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


include "sage/misc/bitset_pxd.pxi"
include "sage/misc/bitset.pxi"
from sage.misc.bitset cimport FrozenBitset, Bitset

cpdef can_push(list neighbors, FrozenBitset subgraph, FrozenBitset filled):
    """
    INPUT: an array of neighbors
    a set of subgraph vertices (we are then assuming we are working on an induced subgraph)
    a set of filled vertices
    
    OUTPUT: True, if a push can happen, plus the new filled set
    False, if a push can't happen, and a copy of the original filled set
    """
    cdef bitset_s *filled_bitset = filled._bitset
    cdef bitset_s *subgraph_bitset = subgraph._bitset
    cdef int num_vertices = filled_bitset.size

    cdef bitset_t unfilled_neighbors
    bitset_init(unfilled_neighbors, num_vertices)

    cdef bitset_s *current_neighbors

    cdef bint can_push=False

    cdef FrozenBitset v
 
    # Make a copy of the filled bitset
    cdef FrozenBitset ret = FrozenBitset(None, capacity=num_vertices)
    cdef bitset_s *ret_bitset=ret._bitset
    bitset_copy(ret_bitset, filled_bitset)

    # Go through each neighborhood.
    cdef int new_filled
    

    cdef int n = bitset_first(filled_bitset)
    while n>=0:
        v = neighbors[n]
        current_neighbors = v._bitset
        bitset_intersection(unfilled_neighbors, current_neighbors, subgraph_bitset)
        bitset_difference(unfilled_neighbors, unfilled_neighbors, filled_bitset)
        new_filled = bitset_first(unfilled_neighbors)
        if new_filled >= 0:
            # look for second unfilled neighbor
            if bitset_next(unfilled_neighbors, new_filled+1) < 0:
                # No more unfilled neighbors
                can_push=True
                break
        n = bitset_next(filled_bitset, n+1)

    bitset_free(unfilled_neighbors)

    if can_push:
        bitset_add(ret_bitset, new_filled)
    return can_push, ret
    
    
    

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
