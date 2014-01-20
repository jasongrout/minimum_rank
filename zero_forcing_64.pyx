# -*- coding: utf-8 -*-
"""
Zero forcing

This module implements zero forcing using fast bitsets, automorphisms
of the graph, and a brute-force approach to trying various bitsets.
"""


#######################################################################
#
# Copyright (C) 2008 Jason Grout.
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

from sage.misc.misc import verbose

ctypedef unsigned long long bitset_t
cdef int BITSET_SIZE = 64

# If you change BITSET_SIZE and bitset_t, you must search for BITSET_SIZE 
# in the comments below and change the corresponding numbers.

cdef inline bitset_t bitset_set(bitset_t bitset, int pos):
    return bitset | (<bitset_t>1<<pos)

cdef inline bitset_t bitset_clear(bitset_t bitset, int pos):
    return bitset & ~(<bitset_t>1<<pos)

cdef inline bitset_t bitset_union(bitset_t bitset, bitset_t bitset2):
    return (bitset) | (bitset2)

cdef inline bitset_t bitset_intersection(bitset_t bitset, bitset_t bitset2):
    return (bitset) & (bitset2)

cdef inline bitset_t bitset_difference(bitset_t bitset, bitset_t bitset2):
    return (bitset) & ~(bitset2)

cdef inline int bitset_check(bitset_t bitset, int pos):
    return (bitset>>pos)&1

cdef bitset_t BITSET_EMPTY = <bitset_t>0

cdef inline bitset_t bitset_full(int length):
    return(<bitset_t>1<<length)-1

cdef int bitset_find(bitset_t bitset, int pos):
    cdef int i
    for i from pos<=i<BITSET_SIZE:
        if bitset_check(bitset,i):
            return i
    return -1

cdef int bitset_next(bitset_t bitset, int pos):
    return bitset_find(bitset, pos+1)

cdef bitset_t zeros_game(bitset_t neighbor_list[64], bitset_t initial_set): # BITSET_SIZE
    cdef int i,j, active_pos, first_nonzero_neighbor, second_nonzero_neighbor

    cdef bitset_t new_zero_set = initial_set
    cdef bitset_t zero_set = BITSET_EMPTY
    cdef bitset_t active_zero_set = BITSET_EMPTY
    cdef bitset_t inactive_zero_set = BITSET_EMPTY
    cdef bint another_run=1
    cdef bitset_t nonzero_neighbors[64] #BITSET_SIZE
    for i from 0<=i<BITSET_SIZE:
        nonzero_neighbors[i] = neighbor_list[i]
    while another_run:
        another_run=0

        # Update the sets
        zero_set = bitset_union(zero_set, new_zero_set)
        active_zero_set = bitset_union(active_zero_set, new_zero_set)
        active_zero_set = bitset_difference(active_zero_set, inactive_zero_set)
        new_zero_set = BITSET_EMPTY
        inactive_zero_set = BITSET_EMPTY

        active_pos = bitset_find(active_zero_set, 0)

        while active_pos >= 0:
            nonzero_neighbors[active_pos] = bitset_difference( nonzero_neighbors[active_pos], zero_set)
            first_nonzero_neighbor = bitset_find(nonzero_neighbors[active_pos],0)
            if first_nonzero_neighbor == -1:
                inactive_zero_set = bitset_set(inactive_zero_set, active_pos)
            else:
                second_nonzero_neighbor = bitset_next(nonzero_neighbors[active_pos],first_nonzero_neighbor)
                if second_nonzero_neighbor == -1:
                    new_zero_set = bitset_set(new_zero_set, first_nonzero_neighbor)
                    inactive_zero_set = bitset_set(inactive_zero_set,active_pos)
                    another_run=1
            active_pos = bitset_next(active_zero_set, active_pos)
    return zero_set


cdef bitset_t positions_to_bitset(int positions[64], int length): # BITSET_SIZE
    cdef bitset_t tmp = BITSET_EMPTY
    cdef int i
    for i from 0<=i<length:
#        print positions[i], " "
        tmp = bitset_set(tmp, positions[i])
#    print "\n"
    return tmp

cdef bitset_t binary_digits_to_bitset(digits, int length):
    cdef bitset_t tmp = BITSET_EMPTY
    cdef int i
    for i from 0<=i<length:
        if digits[i] != 0:
            tmp = bitset_set(tmp, i)
    return tmp


cdef inline bint earlier( bitset_t tuple1, bitset_t tuple2, int n):
    cdef int i
    for i from 1 <= i < n+1:
        if bitset_check(tuple1,n-i) < bitset_check(tuple2,n-i):
#            verbose("%s < %s"%([bitset_check(tuple1,i) for i in xrange(n)], [bitset_check(tuple2,i) for i in xrange(n)]))

            return True
        if bitset_check(tuple1,n-i) > bitset_check(tuple2,n-i):
#            verbose("%s > %s"%([bitset_check(tuple1,i) for i in xrange(n)], [bitset_check(tuple2,i) for i in xrange(n)]))

            return False 
#    verbose("%s = %s"%([bitset_check(tuple1,i) for i in xrange(n)], [bitset_check(tuple2,i) for i in xrange(n)]))
    return False



cdef inline bitset_t permute_diag ( bitset_t d, int *p, int n):
    """
    d = original diag
    p = permutation list (which is [i%n for i in (~p).list])
    returns the permuted diagonal
    """
    cdef int i
    cdef bitset_t dp = BITSET_EMPTY
    if bitset_check(d, p[n-1]):
        dp = bitset_set(dp,0)
    for i from 1<=i<n:
        if bitset_check(d, p[i-1]):
            dp = bitset_set(dp, i)
    return dp



cpdef zero_forcing_set_bruteforce_cython_connected(graph, int upper_bound=-1):
    cdef int n=len(graph.vertices())
    cdef int comb[64] # BITSET_SIZE
    cdef bitset_t adjacency[64] # BITSET_SIZE
    cdef bitset_t diag
    cdef int i, j, k
    cdef bitset_t result
    cdef int mindegree
    cdef int zero_degree_vertices
    

    # We assume that the graph is connected and that
    # the size of the graph is <= BITSET_SIZE
    if n > BITSET_SIZE:
        raise ValueError, "Graph is too large; maximum size is %s"%BITSET_SIZE
    if set(graph.vertices())!=set(xrange(n)):
        raise ValueError, "Graph vertices must be labeled 0 through n-1; use the graph.relabel() command"
    if n == 1:
        return 1, [0], 0, 0, 1
    for i from 0<=i<n:
            adjacency[i] = binary_digits_to_bitset(graph.am().row(i), n)
    for i from n<=i<BITSET_SIZE:
            adjacency[i] = BITSET_EMPTY

    if upper_bound == -1:
        upper_bound = n-1


    # START permutation code
    cdef int saved = 0
    cdef bint done
    cdef bitset_t diag_perm # BITSET_SIZE


    # For the shortcutting, we need to be able to test if some legal permutation
    # has been seen before.  So if an automorphism of the graph takes vertex 1 to 3 to 2
    # I.e., the permutation (1,3,2) is an automorphism
    # and we have the diagonal (0,1,1), then that is equivalent to (1,1,0), which has already
    # been taken care of before.  We can see it in this case by taking the list 
    # representation of (1,3,2) (namely, 
    # We want to compare the last element of the permuted diag with the original diag
    # Then the second-to-last element, etc.
    gp = graph.automorphism_group()
    # Get the generators
    gens = set(gp.gens())
    # and the inverses of the generators
    gens.update(set([~p for p in gens]))
    # Get more of the group by taking products of generators
    gens.update(set([x*y*z for x in gens for y in gens for z in gens]))
    # We really want ~p, but since we have all elements and their inverses
    # it makes no difference if we ask for p.list() or (~p).list()
    perm_p = [p.domain() for p in gens]
    verbose("%s"%perm_p)
    cdef int perms[128][32]
    cdef int num_perms = min(128, len(perm_p))


    for i from 0<=i<num_perms:
        for k from 0<=k<n:
            try:
                perms[i][k] = perm_p[i][k]
            except IndexError:
                perms[i][k] = (k+1) % n

    print "%d permutations"%num_perms



    mindegree=min(graph.degree())
    print "Min degree is %d, so starting from there"%mindegree
    for k from mindegree <= k <= upper_bound:
        print "Investigating subsets of size %s"%k

        # Some code to generate all combinations of n things
        # taken k at a time
        # initialize
        for i from 0<=i<k:
            comb[i]=i
        
        diag = positions_to_bitset(comb,k)
#        verbose("Trying %s"%[bitset_check(diag,i) for i in xrange(n)])
        # Check to see if we've seen this diag before
        done=0
        for p_i from 0<=p_i<num_perms:
            diag_perm = permute_diag(diag, perms[p_i], n)
            if earlier(diag_perm, diag,n):
#                verbose("%s -> %s by %s"%([bitset_check(diag,i) for i in xrange(n)], 
#                            [bitset_check(diag_perm,i) for i in xrange(n)],
#                            [perms[p_i][i] for i in xrange(n)]))

                done=1
                break
        if done:
            saved +=1
        else:
            # diag[j] contains the bit for the jth vertex
            result=zeros_game(adjacency, diag)
#            verbose("%s gives %s"%([bitset_check(diag,i) for i in xrange(n)], [bitset_check(result,i) for i in xrange(n)]))
            if result == bitset_full(n):
                return k,[comb[j] for j in xrange(k)], saved, num_perms, result
    
        while 1:        
            i = k-1
            comb[i] += 1
            while i>=0 and comb[i]>=n-k+1+i:
                i -= 1
                comb[i] += 1
            if comb[0] > n-k:
                break
            
            for j from i+1<=j<k:
                comb[j]=comb[j-1]+1

            diag = positions_to_bitset(comb,k)
#            verbose("Trying %s"%[bitset_check(diag,i) for i in xrange(n)])

            # Check to see if we've seen this diag before
            done=0
            for p_i from 0<=p_i<num_perms:
                diag_perm = permute_diag(diag, perms[p_i], n)
                if earlier(diag_perm, diag,n):
#                    verbose("%s -> %s by %s"%([bitset_check(diag,i) for i in xrange(n)], 
#                                [bitset_check(diag_perm,i) for i in xrange(n)],
#                                [perms[p_i][i] for i in xrange(n)]))
    
                    done=1
                    break
            if done:
                saved +=1
            else:
                # diag[j] contains the bit for the jth vertex


                result=zeros_game(adjacency, diag)
#                verbose("%s gives %s"%([bitset_check(diag,i) for i in xrange(n)], [bitset_check(result,i) for i in xrange(n)]))
                if result == bitset_full(n):
                    return k,[comb[j] for j in xrange(k)], saved, num_perms, result
#            if saved%10000==0:
#                print "Saved %d"%saved
    
    return False

cpdef zero_forcing_set_bruteforce_cython(graph, upper_bound=-1):
    graph = graph.copy()
    relabeling = graph.relabel(return_map=True)
    labeling = dict([(v,k) for k,v in relabeling.iteritems()])
    connected_components = graph.connected_components_subgraphs()
    n = graph.order()
    if upper_bound == -1:
        upper_bound = n-1
    current_zfs = []
    num_perms = 0
    saved_calculations = 0
    for g in connected_components:
        size, zfs, saved, perms, _ = zero_forcing_set_bruteforce_cython_connected(g, upper_bound - len(current_zfs))
        if zfs:
            current_zfs.extend([labeling[i] for i in zfs])
            saved_calculations += saved
            num_perms += perms
        else:
            return False
    return len(current_zfs), current_zfs, saved_calculations, num_perms

# zfs_size, zfs_set, saved_calculations, num_automorphisms = zero_forcing_set_bruteforce_cython(g)
# print "\nZFS minimum size is %d, given by the set:\n %s\n%d calculations were skipped using %d automorphisms of the graph."%(zfs_size, zfs_set, saved_calculations, num_automorphisms)

#min_ranks = [(1,0), (2,0), (3,1), (4,0), (5,1), (6,2), (7,1), (8,0), (9,1), (10,2), (11,2), (12,1), (13,2), (14,3), (15,2),    (16,2),(17,2), (18,1), (19,0), (20,1), (21,2), (22,2), (23,1),    (24,2), (25,3), (26,3), (27,2), (28,2), (29,2), (30,3), (31,4),    (32,2), (33,2), (34,3),(35,3), (36,3), (37,3), (38,3), (39,1),    (40,3), (41,3), (42,2), (43,3), (44,2), (45,2), (46,2), (47,3),    (48,2), (49,2), (50,2), (51,2), (52,1), (53,0), (54,1), (55,2),    (56,2), (57,1), (58,2), (59,3), (60,3), (61,3), (62,2), (63,2),    (64,2), (65,3), (66,4), (67,2), (68,3), (69,4), (70,4),(71,2),    (72,3), (73,3), (74,3), (75,3), (76,3), (77,2), (78,3), (79,4),    (80,4), (81,4), (82,3), (83,5), (84,3), (85,3), (86,1), (87,3),    (88,3),(89,2), (90,3), (91,2), (92,3), (93,4), (94,4), (95,4),    (96,4), (97,4), (98,4), (99,4), (100,3), (101,3), (102,4), (103,4),    (104,4), (105,4),(106,2), (107,2), (108,2), (109,3), (110,2),    (111,4), (112,4), (113,4), (114,3), (115,4), (116,2), (117,3),    (118,4), (119,3), (120,4), (121,3),(122,4), (123,4), (124,4),    (125,3), (126,3), (127,4), (128,4), (129,3), (130,3), (131,2),    (132,2), (133,3), (134,3), (135,3), (136,4), (137,4),(138,3),    (139,4), (140,3), (141,3), (142,3), (143,3), (144,3), (145,3),    (146,2), (147,4), (148,4), (149,3), (150,3), (151,3), (152,4),    (153,3),(154,3), (155,2), (156,3), (157,3), (158,3), (159,3),    (160,3), (161,2), (162,3), (163,3), (164,4), (165,2), (166,3),    (167,4), (168,3), (169,3),(170,3), (171,3), (172,3), (173,3),    (174,3), (175,2), (176,1), (177,3), (178,3), (179,3), (180,3),    (181,3), (182,3), (183,3), (184,3), (185,3),(186,3), (187,3),    (188,3), (189,2), (190,2), (191,2), (192,3), (193,3), (194,2),    (195,2), (196,3), (197,2), (198,3), (199,2), (200,2), (201,2),(202,3), (203,2), (204,2), (205,2), (206,2), (207,2), (208,1),    (209,0), (210,1), (211,2), (212,2), (213,1), (214,2), (215,3),    (216,3), (217,3),(218,2), (219,2), (220,2), (221,3), (222,4),    (223,2), (224,3), (225,4), (226,4), (227,4), (228,2), (229,3),    (230,3), (231,3), (232,3), (233,3),(234,2), (235,3), (236,4),    (237,4), (238,4), (239,3), (240,5), (241,3), (242,3), (243,3),    (244,4), (245,4), (246,5), (247,5), (248,3), (249,1),(250,3),    (251,3), (252,2), (253,3), (254,2), (255,3), (256,4), (257,4),    (258,4), (259,4), (260,4), (261,4), (262,3), (263,4), (264,3),    (265,4),(266,4), (267,4), (268,4), (269,2), (270,2), (271,3),    (272,4), (273,4), (274,4), (275,4), (276,5), (277,4), (278,4),    (279,5), (280,5), (281,4),(282,4), (283,4), (284,5), (285,3),    (286,6), (287,4), (288,4), (289,4), (290,2), (291,2), (292,3),    (293,2), (294,4), (295,4), (296,4), (297,3),(298,4), (299,2),    (300,3), (301,4), (302,3), (303,4), (304,3), (305,4), (306,4),    (307,4), (308,3), (309,3), (310,4), (311,3), (312,4), (313,3),   (314,3), (315,4), (316,4), (317,5), (318,4), (319,4), (320,5),    (321,5), (322,5), (323,4), (324,4), (325,4), (326,3), (327,5),    (328,5), (329,4),(330,4), (331,5), (332,5), (333,5), (334,5),    (335,3), (336,5), (337,5), (338,5), (339,4), (340,5), (341,5),    (342,5), (343,4), (344,4), (345,4),(346,4), (347,3), (348,5),    (349,5), (350,5), (351,5), (352,3), (353,5), (354,3), (355,2),    (356,2), (357,3), (358,3), (359,3), (360,4), (361,4),(362,3),    (363,4), (364,3), (365,3), (366,3), (367,3), (368,3), (369,3),    (370,2), (371,4), (372,4), (373,3), (374,3), (375,3), (376,3),    (377,4),(378,3), (379,4), (380,4), (381,5), (382,3), (383,5),    (384,4), (385,5), (386,4), (387,3), (388,4), (389,4), (390,5),    (391,5), (392,4), (393,5),(394,5), (395,4), (396,4), (397,3),    (398,5), (399,5), (400,5), (401,5), (402,5), (403,4), (404,4),    (405,4), (406,4), (407,4), (408,4), (409,4),(410,4), (411,4),    (412,5), (413,5), (414,5), (415,4), (416,4), (417,3), (418,3),    (419,4), (420,4), (421,5), (422,5), (423,5), (424,4), (425,4),   (426,4), (427,5), (428,4), (429,4), (430,4), (431,4), (432,5),    (433,5), (434,5), (435,5), (436,4), (437,5), (438,5), (439,5),    (440,4), (441,4),(442,4), (443,4), (444,4), (445,5), (446,5),    (447,4), (448,4), (449,4), (450,4), (451,3), (452,2), (453,3),    (454,3), (455,3), (456,3), (457,3),(458,2), (459,3), (460,3),    (461,4), (462,2), (463,3), (464,4), (465,3), (466,3), (467,3),    (468,3), (469,3), (470,3), (471,3), (472,2), (473,3),(474,4),    (475,4), (476,4), (477,4), (478,5), (479,5), (480,4), (481,4),    (482,5), (483,4), (484,5), (485,4), (486,4), (487,4), (488,5),    (489,5),(490,4), (491,4), (492,4), (493,4), (494,4), (495,4),    (496,3), (497,5), (498,4), (499,4), (500,4), (501,3), (502,3),    (503,4), (504,4), (505,4),(506,4), (507,3), (508,5), (509,5),    (510,4), (511,4), (512,5), (513,4), (514,4), (515,5), (516,5),    (517,5), (518,5), (519,4), (520,4), (521,4),(522,4), (523,4),    (524,4), (525,3), (526,5), (527,5), (528,5), (529,5), (530,5),    (531,4), (532,4), (533,5), (534,4), (535,4), (536,4), (537,4),   (538,4), (539,4), (540,4), (541,4), (542,4), (543,4), (544,4),    (545,4), (546,4), (547,4), (548,5), (549,4), (550,4), (551,3),    (552,4), (553,4),(554,3), (555,4), (556,4), (557,3), (558,3),    (559,5), (560,4), (561,5), (562,4), (563,5), (564,4), (565,4),    (566,5), (567,4), (568,4), (569,4),(570,3), (571,4), (572,4),    (573,4), (574,5), (575,5), (576,4), (577,4), (578,4), (579,4),    (580,4), (581,4), (582,2), (583,1), (584,3), (585,3),(586,3),    (587,3), (588,3), (589,3), (590,3), (591,3), (592,3), (593,3),    (594,3), (595,3), (596,2), (597,2), (598,4), (599,4), (600,4),    (601,4),(602,4), (603,4), (604,4), (605,4), (606,4), (607,4),    (608,4), (609,4), (610,3), (611,3), (612,3), (613,4), (614,4),    (615,3), (616,4), (617,4),(618,5), (619,3), (620,4), (621,4),    (622,5), (623,5), (624,4), (625,4), (626,4), (627,4), (628,4),    (629,4), (630,3), (631,4), (632,5), (633,4),(634,4), (635,4),    (636,4), (637,4), (638,4), (639,4), (640,5), (641,4), (642,4),    (643,4), (644,3), (645,4), (646,5), (647,4), (648,4), (649,4),(650,4), (651,4), (652,4), (653,4), (654,4), (655,4), (656,4),    (657,4), (658,4), (659,4), (660,4), (661,4), (662,4), (663,4),    (664,4), (665,4),(666,4), (667,3), (668,3), (669,3), (670,2),    (671,4), (672,4), (673,4), (674,4), (675,4), (676,4), (677,4),    (678,3), (679,4), (680,4), (681,3),(682,5), (683,5), (684,4),    (685,4), (686,3), (687,3), (688,4), (689,3), (690,5), (691,4),    (692,4), (693,4), (694,5), (695,4), (696,4), (697,4),(698,4),    (699,4), (700,4), (701,4), (702,4), (703,4), (704,4), (705,4),    (706,4), (707,4), (708,4), (709,4), (710,5), (711,4), (712,4),    (713,4),(714,4), (715,4), (716,4), (717,4), (718,4), (719,4),    (720,4), (721,3), (722,4), (723,4), (724,4), (725,3), (726,3),    (727,4), (728,4), (729,4),(730,3), (731,2), (732,3), (733,3),    (734,2), (735,2), (736,3), (737,2), (738,3), (739,2), (740,4),    (741,4), (742,4), (743,3), (744,4), (745,2),(746,4), (747,4),    (748,4), (749,4), (750,4), (751,4), (752,4), (753,4), (754,4),    (755,4), (756,4), (757,4), (758,4), (759,4), (760,4), (761,4),  (762,4), (763,4), (764,4), (765,4), (766,4), (767,4), (768,4),    (769,4), (770,4), (771,4), (772,4), (773,4), (774,4), (775,4),    (776,4), (777,4),(778,4), (779,4), (780,3), (781,3), (782,4),    (783,4), (784,4), (785,4), (786,3), (787,3), (788,4), (789,3),    (790,2), (791,3), (792,4), (793,4),(794,4), (795,4), (796,3),    (797,4), (798,4), (799,4), (800,4), (801,3), (802,4), (803,4),    (804,4), (805,4), (806,4), (807,4), (808,4), (809,3),(810,4),    (811,4), (812,3), (813,5), (814,3), (815,3), (816,4), (817,4),    (818,4), (819,4), (820,4), (821,5), (822,4), (823,4), (824,4),    (825,4),(826,4), (827,4), (828,4), (829,3), (830,4), (831,3),    (832,3), (833,4), (834,4), (835,4), (836,4), (837,4), (838,4),    (839,4), (840,5), (841,4),(842,4), (843,4), (844,4), (845,4),    (846,3), (847,4), (848,4), (849,4), (850,4), (851,3), (852,4),    (853,4), (854,4), (855,4), (856,3), (857,4),(858,4), (859,4),    (860,4), (861,4), (862,4), (863,3), (864,4), (865,3), (866,4),    (867,4), (868,4), (869,3), (870,4), (871,4), (872,3), (873,3),  (874,3), (875,4), (876,3), (877,4), (878,3), (879,2), (880,2),    (881,3), (882,2), (883,2), (884,3), (885,3), (886,4), (887,4),    (888,4), (889,4),(890,4), (891,4), (892,4), (893,3), (894,3),    (895,3), (896,3), (897,4), (898,3), (899,4), (900,4), (901,3),    (902,3), (903,4), (904,4), (905,4),(906,3), (907,4), (908,4),    (909,3), (910,4), (911,3), (912,3), (913,3), (914,4), (915,4),    (916,4), (917,4), (918,3), (919,4), (920,4), (921,4),(922,4),    (923,4), (924,3), (925,3), (926,4), (927,4), (928,4), (929,4),    (930,4), (931,3), (932,3), (933,4), (934,4), (935,4), (936,4),    (937,4),(938,4), (939,4), (940,4), (941,4), (942,4), (943,4),    (944,3), (945,3), (946,4), (947,3), (948,3), (949,3), (950,4),    (951,4), (952,4), (953,3),(954,4), (955,4), (956,3), (957,3),    (958,3), (959,4), (960,4), (961,4), (962,4), (963,4), (964,4),    (965,4), (966,4), (967,4), (968,4), (969,4),(970,3), (971,4),    (972,4), (973,3), (974,4), (975,3), (976,4), (977,3), (978,3),    (979,4), (980,4), (981,4), (982,4), (983,3), (984,3), (985,4), (986,4), (987,3), (988,3), (989,4), (990,3), (991,3), (992,4),    (993,4), (994,3), (995,3), (996,3), (997,4), (998,4), (999,4),    (1000,3), (1001,3),(1002,3), (1003,3), (1004,3), (1005,3), (1006,4),    (1007,2), (1008,4), (1009,2), (1010,2), (1011,2), (1012,3), (1013,3),    (1014,3), (1015,4),(1016,4), (1017,3), (1018,3), (1019,3), (1020,3),    (1021,4), (1022,3), (1023,3), (1024,3), (1025,4), (1026,4), (1027,4),    (1028,3), (1029,4),(1030,4), (1031,4), (1032,2), (1033,4), (1034,3),    (1035,3), (1036,3), (1037,3), (1038,3), (1039,4), (1040,3), (1041,4),    (1042,3), (1043,4),(1044,3), (1045,3), (1046,4), (1047,4), (1048,4),    (1049,3), (1050,4), (1051,4), (1052,4), (1053,4), (1054,4), (1055,4),    (1056,3), (1057,3),(1058,4), (1059,4), (1060,3), (1061,4), (1062,3),    (1063,3), (1064,3), (1065,4), (1066,3), (1067,3), (1068,3), (1069,4),    (1070,3), (1071,4),(1072,3), (1073,3), (1074,3), (1075,3), (1076,3),    (1077,3), (1078,4), (1079,3), (1080,4), (1081,3), (1082,4), (1083,4),    (1084,3), (1085,3),(1086,3), (1087,3), (1088,2), (1089,4), (1090,3),    (1091,4), (1092,3), (1093,4), (1094,3), (1095,3), (1096,3), (1097,4),    (1098,3), (1099,3),(1100,3), (1101,4), (1102,3), (1103,3), (1104,3),    (1105,3), (1106,3), (1107,2), (1108,3), (1109,3), (1110,3), (1111,3),    (1112,3), (1113,3),(1114,3), (1115,3), (1116,3), (1117,4), (1118,4),    (1119,3), (1120,3), (1121,4), (1122,3), (1123,3), (1124,3), (1125,3),    (1126,3), (1127,4),(1128,3), (1129,3), (1130,3), (1131,3), (1132,3),    (1133,3), (1134,3), (1135,3), (1136,3), (1137,3), (1138,3), (1139,3),    (1140,2), (1141,4),(1142,4), (1143,3), (1144,3), (1145,4), (1146,3),    (1147,3), (1148,3), (1149,3), (1150,4), (1151,3), (1152,3), (1153,3),    (1154,4), (1155,3),(1156,3), (1157,3), (1158,3), (1159,3), (1160,4),    (1161,3), (1162,3), (1163,3), (1164,2), (1165,3), (1166,3), (1167,3),    (1168,3), (1169,3), (1170,3), (1171,2), (1172,1), (1173,3), (1174,3),    (1175,3), (1176,3), (1177,3), (1178,3), (1179,3), (1180,3), (1181,3),    (1182,3), (1183,3), (1184,2), (1185,3), (1186,3), (1187,4), (1188,2),    (1189,3), (1190,4), (1191,3), (1192,3), (1193,3), (1194,3), (1195,3),    (1196,3), (1197,3), (1198,3), (1199,3), (1200,3), (1201,3), (1202,3),    (1203,3), (1204,3), (1205,3), (1206,2), (1207,3), (1208,2), (1209,3),    (1210,3), (1211,2), (1212,3), (1213,2), (1214,3), (1215,3), (1216,2),    (1217,3), (1218,3), (1219,3), (1220,3), (1221,3), (1222,3), (1223,3),    (1224,3), (1225,3), (1226,2), (1227,3), (1228,3), (1229,2), (1230,2),    (1231,3), (1232,3), (1233,2), (1234,2), (1235,3), (1236,3), (1237,2),    (1238,2), (1239,3), (1240,2), (1241,3), (1242,2), (1243,2), (1244,2),    (1245,2), (1246,3), (1247,2), (1248,2), (1249,2), (1250,2), (1251,2),    (1252,1)]
# import networkx.generators.atlas
# atlas_graphs = [Graph(i) for i in networkx.generators.atlas.graph_atlas_g()]


# from sage.combinat.subset import Subsets



# def Test():
#     looks_good = True
#     print "Testing to ensure that our new zero forcing function's results are indeed below all minimum ranks found before..."
#     atlas_graphs = [Graph(i) for i in networkx.generators.atlas.graph_atlas_g()]
#     for i,k in min_ranks:
#         if (len(atlas_graphs[i].vertices()) - zero_forcing_set(atlas_graphs[i])[0]) > k:
#             print "Utoh at atlas graph: " , i
#             print "Zero forcing set length should be less than or equal to ", k, " but is instead ", zero_forcing_set(atlas_graphs[i],False)
#             looks_good = False
#     if(looks_good):
#         print "Passed first test"
#     if not looks_good:
#         print "Failed first test"
#     looks_good = True
#     print "Testing to see if our our new zero forcing function's results are the exactly the same as the proven old zero forcing set function's results..."
#     for i in range(1,1252):
#         if zero_forcing_set(atlas_graphs[i])[0]!=len(find_zero_forcing_set(atlas_graphs[i])):
#             looks_good = False;
#             print " Something is wrong at ", i
#             print "     New zfs is ",zero_forcing_set(atlas_graphs[i],False)
#             print "     Old zfs is ",find_zero_forcing_set(atlas_graphs[i])
#     if(looks_good):
#         print "Passed second test"
#     if not looks_good:
#         print "Failed second test"


