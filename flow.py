#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 16:34:29 2023

@author: sofiat
"""


import networkx as nx

G = nx.DiGraph()
G.add_edges_from( [
    ('s', 's1', {"capacity": 1}),
    ('s', 's2', {"capacity": 1}),
    ('s', 's3', {"capacity": 1}),
    ('s', 's4', {"capacity": 1}),
    
    ('s1', 'p1', {"capacity": 1}),
    ('s2', 'p1', {"capacity": 1}),
    ('s3', 'p1', {"capacity": 1}),
    
    ('s4', 'p2', {"capacity": 1}),
    
    ('p1', 'l1', {"capacity": 1}),
    ('p2', 'l1', {"capacity": 2}),
    
    ('l1', 't', {"capacity": 2}),
    
    ])

max_flow = nx.max_flow_min_cost(G, 's', 't')
print(max_flow)
# print(G.nodes)
# print(G.edges)

