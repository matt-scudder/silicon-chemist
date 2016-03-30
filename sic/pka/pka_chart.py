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
#NOTE: The order here is of paramount importance such that species are matched correctly. Start with the least specific forms, and then go to the most specific.
#NOTE: The huge amount of ([H]) below are because otherwise the hydrogens we so desperately care about won't get matched in the expression, and we run into problems.
PKA_CHART = [
    #Water
    {"[OHv1]([H])":{"pKa_HA": 52, "pKa_BH": 15.7}}, #v1 because otherwise every alcohol ever would match it. Put in a separate case for alcohols...
    {"[OH2v2]([H])([H])":{"pKa_HA":15.7,"pKa_BH":-1.7}}, #v2 required because it would match e.g. protonated carboxylics otherwise, which have diff. pKa
    {"[OH3]([H])([H])([H])":{"pKa_HA":-1.7,"pKa_BH":None}}, #can't add another H to OH3+, so None for BH
    #Halogens - no pKa_BH for the acids, pretty sure stuff like H2I doesn't exist.
    {"[IH1]([H])":{"pKa_HA":-10,"pKa_BH":None}},
    {"[I!H1]":{"pKa_HA":None,"pKa_BH":-10}}, #Since other stuff that makes halogens better Ls doesn't go into the numeric pKa estimate, just make all Is -10
    {"[BrH1]([H])":{"pKa_HA":-9,"pKa_BH":None}},
    {"[Br!H1]":{"pKa_HA":None,"pKa_BH":-9}},#!H1 because we want halogens that AREN'T attached to H (and thus acidic...)
    {"[ClH1]([H])":{"pKa_HA":-7,"pKa_BH":None}},
    {"[Cl!H1]":{"pKa_HA":None,"pKa_BH":-7}},
    {"[FH1]([H])":{"pKa_HA":3.2,"pKa_BH":None}},
    {"[F!H1]":{"pKa_HA":None,"pKa_BH":3.2}},
    #Alcohols
    {"[CH3][O]([H])":{"pKa_HA": 15.5, "pKa_BH": -2.4}}, #methyl
    {"[CH2][O]([H])":{"pKa_HA": 16, "pKa_BH": -2.4}}, #primary
    {"[CH1][O]([H])":{"pKa_HA": 18, "pKa_BH": -3.5}}, #secondary
    {"[CH0][O]([H])":{"pKa_HA": 19, "pKa_BH": -4}}, #tertiary
    #Protonated alcohols
    {"[CH3][OH2]([H])([H])":{"pKa_HA": -2.4, "pKa_BH": None}}, #methyl
    {"[CH2][OH2]([H])([H])":{"pKa_HA": -2.4, "pKa_BH": None}}, #primary
    {"[CH1][OH2]([H])([H])":{"pKa_HA": -3.5, "pKa_BH": None}}, #secondary
    {"[CH0][OH2]([H])([H])":{"pKa_HA": -4, "pKa_BH": None}} #tertiary

]
