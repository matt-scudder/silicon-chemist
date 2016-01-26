#!/usr/bin/python
#coding=utf-8
"""
Lists the table of source-sink interactions, i.e. which mechanisms a particular pair can perform.
By convention, we map sources to sinks, though it can be done either way.
We assume that all sources can map to all sinks, and vice versa.
"""

INTERACTIONS = {
        ("Y:","H-L"):["proton_transfer"]
        }
