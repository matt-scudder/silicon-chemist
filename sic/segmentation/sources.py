#!/usr/bin/python
#coding=utf-8
"""
SMARTS expressions for identifying sources.
Less complex data structure than for pKa since we have to check all source and sink types
when we're looking at a molecule.
Logic to handle atomic tagging (see specs - object that maps the "pieces" of a source
to the "letters" in the source name) should be left in the segmentation module.
"""

#TODO: Add support for C-
SOURCES = {
            "Y":"[O!v3,S!v3,N!v4,F-,Cl-,Br-,I-]",
            #finds anything with lone pairs, whether - charge or not. The awkward v is because you can't group
            #the v primitives as far as I've been able to find out.
            "C-":"[#6-]" #Need 6 to take in ring carbanions.
        }
            
