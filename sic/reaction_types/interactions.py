#!/usr/bin/python
#coding=utf-8
"""
Lists the table of source-sink interactions, i.e. which mechanisms a particular pair can perform.
By convention, we map sources to sinks, though it can be done either way.
We assume that all sources can map to all sinks, and vice versa.
"""

INTERACTIONS = {
        ("Y","H-L"):["proton_transfer"],
        ("Y","C-L"):["SN2","E2"], #substitution vs. elimination!
        ("Y","C+"):["AN","E1"],
        ("DUM","C-L"):["DN"], #because DN is basically "the sink interacting with itself and falling apart", we want one DUM for every C-L
        ("C-","C-L"):["EB"]
        }
