#!/usr/bin/python
#coding=utf-8
"""
SMARTS expressions for identifying sinks.
Less complex data structure than for pKa since we have to check all source and sink types
when we're looking at a molecule.
Logic to handle atomic tagging (see specs - object that maps the "pieces" of a source
to the "letters" in the source name) should be left in the segmentation module.
"""
#TODO: Figure out how to exclude pKa 40+ stuff from consideration; maybe prune after the source/sink list is built?
#TODO: Add carbon acids to the H-L sink list
SINKS = {
        "H-L":"[O!v1,S,N,F,Cl,Br,I]([H])",
    #anything we know can bond to hydrogen, bonded with it
    #NOTE: This ignores stuff like NaH (because H is "alone" and not bonded)
        "C+":"[#6+]", #carbocations - only take in the carbon! Need 6 in case of carbocations in a ring.
        "C-L":"[#6X4]([!#1&!#6])" #C-L is an sp3 carbon (i.e. 4 connections), and the L is not a C (6) or an H (1), since atomic number is the only way to get SMARTS to understand that.
    }
