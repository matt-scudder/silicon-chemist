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
            "Y:":"[O!X3,S!X3,N!X4,FX0,ClX0,BrX0,IX0]"
            #finds anything with lone pairs, whether - charge or not. The awkward X is because you can't group
            #the X primitives as far as I've been able to find out.
        }
            
