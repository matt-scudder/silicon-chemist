"""
Format for the pKa chart:

"SMARTS_STRING OF FRAGMENT": {"pKa_HA": val,"pKa_BH": val}

If either pKa_HA or BH is unavailable for the particular fragment,
whether due to physical impossibility (pentavalent carbon)
or it not being a number we have, we put None for that slot.
The rest of the program should read these accordingly.
"""
#NOTE: This file will get REALLY big really fast, but only needs to be run once per ReactionState
#NOTE: The order here is of paramount importance such that species are matched correctly. 
# Start with the least specific forms, and then go to the most specific.
PKA_CHART = {
    #Water
    "[Ov1-][H]": {"pKa_HA": None, "pKa_BH": 14.0}, #v1 because otherwise every alcohol would match it. 
    "[Ov2]([H])[H]": {"pKa_HA": 14.0,"pKa_BH":0.0}, #v2 required because it would match e.g. protonated carboxylics otherwise, which have diff. pKa
    "[Ov3+]([H])([H])[H]": {"pKa_HA":0.0,"pKa_BH": None}, #v3 restricts to 3 bonds.  can't add another H to OH3+, so None for BH
    #Halogen acids - no pKa_BH for the acids, pretty sure stuff like H2I doesn't exist.
    "[Iv1][H]": {"pKa_HA":-10,"pKa_BH": None}, #for HI
    "[I!H1]": {"pKa_HA": None,"pKa_BH":-10}, #for I as an L Since other stuff that makes halogens better Ls doesn't go into the numeric pKa estimate, just make all Is -10
    "[Brv1][H]": {"pKa_HA":-9,"pKa_BH": None}, #for HBr
    "[Br!H1]": {"pKa_HA": None,"pKa_BH":-9},#!H1 for Br as an L because we want halogens that AREN'T attached to H (and thus acidic...)
    "[Clv1][H]": {"pKa_HA":-7,"pKa_BH": None}, #for HCl
    "[Cl!H1]": {"pKa_HA": None,"pKa_BH":-7}, # for Cl as an L
    "[Fv1][H]": {"pKa_HA": 3.2,"pKa_BH": None}, #for HF
    "[F!H1]": {"pKa_HA": None,"pKa_BH": 3.2}, # for F as an L
    #Sulfuric acid
    "[SX4](=O)(=O)([O][H])[O][H]": {"pKa_HA":-9, "pKa_BH": None}, #neutral H2SO4
    "[SX4](=O)(=O)([O][H])[O-]": {"pKa_HA": None, "pKa_BH": -9}, #bisulfate monoanion HSO4-
    "[#6][SX4](=O)(=O)[O-]": {"pKa_HA": None, "pKa_BH": -6.5}, #tosylate or mesylate
    "[#6][SX4](=O)(=O)[O][H]": {"pKa_HA":  -6.5, "pKa_BH": None}, #tosyl acid or mesyl acid
    #neutral alcohols
    "[CH3][O]([H])": {"pKa_HA": 15.5, "pKa_BH": -2.4}, #methyl all same PHS book pKa value for ETOH2 since all should be close
    "[CH2][O]([H])": {"pKa_HA": 16, "pKa_BH": -2.4}, #primary
    "[CH1][O]([H])": {"pKa_HA": 18, "pKa_BH": -2.4}, #secondary
    "[CH0][O]([H])": {"pKa_HA": 19, "pKa_BH": -2.4}, #tertiary
    "[c]:[cX3][OX2][H]": {"pKa_HA": 10, "pKa_BH": -6.4}, #phenol 
    #deprotonated alcohols
    "[CH3][O-]": {"pKa_HA": None, "pKa_BH": 15.2}, #methyl with +H (all from PHS book pKa values)
    "[CH2][O-]": {"pKa_HA": None, "pKa_BH": 16}, #primary alkoxide
    "[CH1][O-]": {"pKa_HA": None, "pKa_BH": 17}, #secondary alkoxide
    "[CH0][O-]": {"pKa_HA": None, "pKa_BH": 19}, #tertiary alkoxide
    "[c]:[cX3][O-]": {"pKa_HA": None, "pKa_BH": 10}, #phenolate
    #Protonated alcohols
    "[CH3][O+H]([H])": {"pKa_HA": -2.4, "pKa_BH": None}, #methyl +H
    "[CH2][O+H]([H])": {"pKa_HA": -2.4, "pKa_BH": None}, #primary +H
    "[CH1][O+H]([H])": {"pKa_HA": -2.4, "pKa_BH": None}, #secondary +H
    "[CH0][O+H]([H])": {"pKa_HA": -2.4, "pKa_BH": None}, #tertiary +H
    "[c]:[cX3][O+]([H])[H]": {"pKa_HA": -6.4, "pKa_BH": None}, #phenol +H
    #Tetrahedral intermediate
    "[#6][CX4]([OD2])[O][H]": {"pKa_HA": 13.3, "pKa_BH": -3}, #neutral ester tetrahedral intermediate
    "[#6][CX4]([OD2])[O-]": {"pKa_HA": None, "pKa_BH": 13.3}, #anionic ester tetrahedral intermediate
    "[#6][CX4]([OD2])[O+H2]": {"pKa_HA": -3, "pKa_BH": None}, #cationic ester tetrahedral intermediate
    #Nitriles (carbon acid) and Cyanides
    "[NX1]#[CX2]([H])": {"pKa_HA": 9.2, "pKa_BH": None}, #cyanide, with proton
    "[NX1]#[#6-]": {"pKa_HA": None, "pKa_BH": 9.2}, #cyanide, deprotonated
    "[NX1]#[CX2][CH2,CH1,CH0]([H])": {"pKa_HA": 25,"pKa_BH": None}, #nitrile CH, with proton, primary, secondary, tertiary
    "[NX1]#[CX2][#6-]": {"pKa_HA": None,"pKa_BH": 25}, #nitrile CH, deprotonated, primary, secondary, tertiary
    "[O-][N+X3](=O)[#6H2,#6H1,#6H0]([H])": {"pKa_HA": 10.2,"pKa_BH": None}, #Nitroalkane, with proton
    "[O-][N+X3](=O)[#6-]": {"pKa_HA": None,"pKa_BH": 10.2}, #Nitroalkane, deprotonated
    #Aldehydes
    "[#6X3H1](=O)[#6H2,#6H1,#6H0]([H])": {"pKa_HA": 16.7,"pKa_BH": -10}, #primary, secondary, tertiary
    "[#6X3H1](=O)[#6-]": {"pKa_HA": None,"pKa_BH": 16.7}, #aldehyde enolate pKa_BH
    "[#6][#6X3H1]=[OX2H1+]": {"pKa_HA":-10,"pKa_BH": None}, #protonated aldehyde pKa
    "[#6!-][#6X3H1]=[OX1]": {"pKa_HA": None,"pKa_BH":-10}, #neutral aldehyde pKa_BH 
    #Ketones
    "[#6][CX3](=O)[#6H2,#6H1,#6H0]([H])": {"pKa_HA": 19.2,"pKa_BH": None}, #neutral ketone pKa
    "[#6][CX3](=O)[#6-]": {"pKa_HA": None,"pKa_BH": 19.2}, #ketone enolate pKa_BH
    "[#6][CX3](C)=[OX2H1+]": {"pKa_HA":-7,"pKa_BH": None}, #protonated ketone pKa
    "[#6!-][CX3](C)=[OX1]": {"pKa_HA": None,"pKa_BH":-7}, #neutral ketone pKa_BH
    #Esters
    "[#6H1][CX3](=O)[OX2H0][#6H3,#6H2,#6H1,#6H0]": {"pKa_HA": 25.6,"pKa_BH": -6.5}, #neutral ester pKa and pKa_BH
    "[#6-][CX3](=O)[OX2H0][#6H3,#6H2,#6H1,#6H0]": {"pKa_HA": None,"pKa_BH": 25.6}, #ester enolate pKa_BH
    "[#6!-][CX3](=[OX2H1+])[OX2H0][#6H3,#6H2,#6H1,#6H0]": {"pKa_HA": -6.5,"pKa_BH": None}, #protonated ester pKa
    "[#6-]([CX3](=[OX1])[OX2H0][#6])[CX3](=[OX1])[OX2H0][#6]": {"pKa_HA": None,"pKa_BH": 13}, #Malonic ester pKa_BH
    "[C]([H])([CX3](=[OX1])[OX2H0][#6])[CX3](=[OX1])[OX2H0][#6]": {"pKa_HA": None,"pKa_BH": 13}, #Malonic ester pKa
    "[#6-]([CX3](=[OX1])[OX2H0][#6])[CX3](=[OX1])[#6]": {"pKa_HA": None,"pKa_BH": 10.7}, #Acetoacetic ester pKa_BH
    "[C]([H])([CX3](=[OX1])[OX2H0][#6])[CX3](=[OX1])[#6]": {"pKa_HA": 10.7,"pKa_BH": None}, #Acetoacetic ester pKa
    #Carboxylic acids
    "[#6][CX3](=O)[OX2H1]": {"pKa_HA": 4.8,"pKa_BH": -6}, #neutral carboxylic acid pKa and pKa_BH
    "[#6!-][CX3](=[OX2H1+])[OX2H1]": {"pKa_HA": -6,"pKa_BH": None}, #protonated carboxylic acid pKa
    "[#6][CX3](=O)[OX1-]": {"pKa_HA": None,"pKa_BH": 4.8}, #carboxylate anion pKa_BH
    #Amides
    "[NH0][CX3](=O)[#6H3,#6H2,#6H1]": {"pKa_HA": 28,"pKa_BH": -0.5}, #neutral tertiary amide CH pKa and pKa_BH
    "[NH0][CX3](=O)[#6-]": {"pKa_HA": None,"pKa_BH": 28}, #Tertiary Amide C enolate pKa_BH
    "[NH1,NH2][CX3](=O)[#6]": {"pKa_HA": 17,"pKa_BH": -0.5}, #primary or secondary Amide NH pKa and pKa_BH
    "[N-][CX3](=O)[#6]": {"pKa_HA": None,"pKa_BH": 17}, #Amidate anion pKa_BH
    #Thiols
    "[SHv1]([H])": {"pKa_HA": 12.90, "pKa_BH": 6.97}, #v1 means only one bond to sulfur 
    "[SH2v2]([H])([H])": {"pKa_HA": 6.97,"pKa_BH":-6}, #v2 means only two bonds to sulfur 
    "[CH3][S-]": {"pKa_HA": None, "pKa_BH": 10.6}, #methyl
    "[CH3][S]([H])": {"pKa_HA": 10.6, "pKa_BH": -7},
    "[c]:[cX3][SX2][H]": {"pKa_HA": 6.5, "pKa_BH": -7}, #Thiophenol 
    "[c]:[cX3][S-]": {"pKa_HA": None, "pKa_BH": 6.5}, #Thiophenolate 
    #protonated ethers
    "[#6][O][#6]": {"pKa_HA": None, "pKa_BH": -2.4}, #methyl
    "[#6][O+]([H])[#6]": {"pKa_HA": -3.5, "pKa_BH": None}, #methyl
    #Amines
    "[NX3]([H])([H])[H]": {"pKa_HA": 35,"pKa_BH": 9.2}, #Ammonia NH pKa and pKa_BH
    "[NX2-]([H])[H]": {"pKa_HA": None,"pKa_BH": 35}, #Ammonia anion NH pKa and pKa_BH
    "[NH2][CX4]": {"pKa_HA": 36,"pKa_BH": 10.7}, #primary amine NH pKa and pKa_BH
    "[NH-][CX4]": {"pKa_HA": None,"pKa_BH": 36}, #primary amine nion pKa_BH
    "[NH]([CX4])[CX4]": {"pKa_HA": 36,"pKa_BH": 10.7}, #secondary amine NH pKa and pKa_BH
    "[N-]([CX4])[CX4]": {"pKa_HA": None,"pKa_BH": 36}, #secondary amine anion pKa_BH
    "[N]([CX4])([CX4])[CX4]": {"pKa_HA": None,"pKa_BH": 10.7}, #tertiary amine pKa_BH
    "[c]:[cX3][NX3]([H])[H]": {"pKa_HA": 28, "pKa_BH": 4.6}, #aniline 
    "[c]:[cX3][NX2-][H]": {"pKa_HA": 28, "pKa_BH": 4.6}, #aniline anion
     #Carbon acids without EWG
    "[H,C][CX2]#[CX2]([H])": {"pKa_HA": 25, "pKa_BH": None}, #Acetylene with proton
    "[H,C][CX2]#[#6-]": {"pKa_HA": None, "pKa_BH": 25}, #Acetylene, deprotonated
    "[H][C]-c1=cc=cc=c1": {"pKa_HA": 40, "pKa_BH": None}, #Benzyl with proton
    "[C-]-c1=cc=cc=c1": {"pKa_HA": None, "pKa_BH": 40}, #Benzyl, deprotonated
    "[H]-c1=cc=cc=c1": {"pKa_HA": 43, "pKa_BH": None}, #Benzene with proton
    "c1=[c-]c=cc=c1": {"pKa_HA": None, "pKa_BH": 43}, #Benzene, deprotonated
    "C=C[C][H]": {"pKa_HA": 43, "pKa_BH": None}, #Allyl with proton
    "C=C[C-]": {"pKa_HA": None, "pKa_BH": 44}, #Allyl, deprotonated
    "[CH2]=[C,CH,CH2]": {"pKa_HA": 44, "pKa_BH": None}, #Ethylene with proton
    "[C-]=[C,CH,CH2]": {"pKa_HA": None, "pKa_BH": 44}, #Ethylene, deprotonated
    "[H][CH2]-[C,CH,CH2]": {"pKa_HA": 50, "pKa_BH": None}, #Ethyl with proton
    "[CH2-]-[C,CH,CH2]": {"pKa_HA": None, "pKa_BH": 50}, #Ethyl, deprotonated
    #Hydride
    "[H-]": {"pKa_HA": None, "pKa_BH": 35} #Hydride anion
}
 
