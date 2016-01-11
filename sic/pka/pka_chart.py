#!/usr/bin/python
#coding=utf-8
"""
Format for the pKa chart:

{"SMARTS_STRING OF FRAGMENT": {"pKa_HA":val,"pKa_BH":val}}

If either pKa_HA or BH is unavailable for the particular fragment,
whether due to physical impossibility (pentavalent carbon)
or it not being a number we have, we put None for that slot.
The rest of the program should read these accordingly.
"""
PKA_CHART = [
    #Water
    {"[OH2X2]":{"pKa_HA":15.7,"pKa_BH":-1.7}}, #X2 required because it would match e.g. protonated carboxylics otherwise, which have diff. pKa
    {"[OH3]":{"pKa_HA":-1.7,"pKa_BH":None}}, #can't add another H to OH3+, so None for BH
    #Halogens - no pKa_BH for the acids, pretty sure stuff like H2I doesn't exist.
    {"[IH1]":{"pKa_HA":-10,"pKa_BH":None}},
    {"[I!H1]":{"pKa_HA":None,"pKa_BH":-10}}, #Since other stuff that makes halogens better Ls doesn't go into the numeric pKa estimate, just make all Is -10
    {"[BrH1]":{"pKa_HA":-9,"pKa_BH":None}},
    {"[Br!H1]":{"pKa_HA":None,"pKa_BH":-9}},#!H1 because we want halogens that AREN'T attached to H (and thus acidic...)
    {"[ClH1]":{"pKa_HA":-7,"pKa_BH":None}},
    {"[Cl!H1]":{"pKa_HA":None,"pKa_BH":-7}},
    {"[FH1]":{"pKa_HA":3.2,"pKa_BH":None}},
    {"[F!H1]":{"pKa_HA":None,"pKa_BH":3.2}}
]
