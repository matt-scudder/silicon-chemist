"""
Format for the pKa chart:

{"SMARTS_STRING OF FRAGMENT": {"pKa_HA":val,"pKa_BH":val}}

If either pKa_HA or BH is unavailable for the particular fragment,
whether due to physical impossibility (pentavalent carbon)
or it not being a number we have, we put None for that slot.
The rest of the program should read these accordingly.
"""
#NOTE: This file will get REALLY big really fast, but only needs to be run once per ReactionState
#NOTE: The order here is of paramount importance such that species are matched correctly. Start with the least specific forms, and then go to the most specific.
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
    #protonated
    {"[CH3][O]([H])":{"pKa_HA": 15.5, "pKa_BH": -2.4}}, #methyl
    {"[CH2][O]([H])":{"pKa_HA": 16, "pKa_BH": -2.4}}, #primary
    {"[CH1][O]([H])":{"pKa_HA": 18, "pKa_BH": -3.5}}, #secondary
    {"[CH0][O]([H])":{"pKa_HA": 19, "pKa_BH": -4}}, #tertiary
    #deprotonated
    {"[CH3][O-]":{"pKa_HA": None, "pKa_BH": 15.5}}, #methyl
    {"[CH2][O-]":{"pKa_HA": None, "pKa_BH": 16}}, #primary
    {"[CH1][O-]":{"pKa_HA": None, "pKa_BH": 18}}, #secondary
    {"[CH0][O-]":{"pKa_HA": None, "pKa_BH": 19}}, #tertiary
    #Protonated alcohols
    {"[CH3,CH2][OH2]([H])":{"pKa_HA": -2.4, "pKa_BH": None}}, #methyl + primary
    {"[CH1][OH2]([H])":{"pKa_HA": -3.5, "pKa_BH": None}}, #secondary
    {"[CH0][OH2]([H])":{"pKa_HA": -4, "pKa_BH": None}}, #tertiary
    #Nitriles (carbon acid) and Cyanides
    {"[NX1]#[CX2]([H])":{"pKa_HA": 9.2, "pKa_BH": None}}, #cyanide, with proton
    {"[NX1]#[#6-]":{"pKa_HA": None, "pKa_BH": 9.2}}, #cyanide, deprotonated
    {"[NX1]#[CX2][CH3,CH2,CH]([H])":{"pKa_HA":25,"pKa_BH":None}}, #nitrile CH, protonated, primary, secondary, tertiary
    {"[NX1]#[CX2][#6-]":{"pKa_HA":None,"pKa_BH":25}}, #nitrile CH, deprotonated, primary, secondary, tertiary
    {"[O-][N+X3](=O)[#6H3,#6H2,#6H1]([H])":{"pKa_HA":10.2,"pKa_BH":None}}, #Nitro, protonated
    {"[O-][N+X3](=O)[#6-]":{"pKa_HA":None,"pKa_BH":10.2}}, #Nitro, deprotonated
    #Aldehydes
    {"[#6X3H1](=O)[#6H3,#6H2,#6H1]([H])":{"pKa_HA":16.7,"pKa_BH":None}}, #primary, secondary, tertiary
    {"[#6X3H1](=O)[#6-]":{"pKa_HA":None,"pKa_BH":16.7}},
    # Oxygen Acids
    {"[#6][#6X3](C)=[OX2H1+]":{"pKa_HA":-7,"pKa_BH":None}},
    {"[#6][#6X3](C)=[OX1]":{"pKa_HA":None,"pKa_BH":-7}},
    {"[#6][#6X3H1]=[OX2H1+]":{"pKa_HA":-10,"pKa_BH":None}},
    {"[#6!-][#6X3H1]=[OX1]":{"pKa_HA":None,"pKa_BH":-10}},
    {"[CH3][CH2][O+H2]([H])([H])":{"pKa_HA":-2.4,"pKa_BH":None}},
    {"[CH3][CH2][O+H1]([H])":{"pKa_HA":None,"pKa_BH":-2.4}},
    #Sulfuric acid
    {"[$([SX4](=O)(=O)(O)O)][O]([H])":{"pKa_HA":-9, "pKa_BH": None}}, #protonated
    {"[$([SX4+1](O)([O-])(O)O)]":{"pKa_HA":None, "pKa_BH": -9}}, #deprotonated
    #Ketones
    {"[CX3C1](=O)[#6H3,#6H2,#6H1]([H])":{"pKa_HA":19.2,"pKa_BH":None}},#protonated
    {"[CX3C1](=O)[#6-]":{"pKa_HA":None,"pKa_BH":19.2}},#deprotonated
    #Esters
    {"[CX3O1](=O)[#6H3,#6H2,#6H1]([H])":{"pKa_HA":25.6,"pKa_BH":None}},#protonated
    {"[CX3O1](=O)[#6-]":{"pKa_HA":None,"pKa_BH":25.6}},#deprotonated
    #Amides
    {"[CX3N1](=O)[#6H3,#6H2,#6H1]([H])":{"pKa_HA":28,"pKa_BH":None}},#protonated
    {"[CX3N1](=O)[#6-]":{"pKa_HA":None,"pKa_BH":28}}#deprotonated
]
