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
        ("C-","C-L"):["EB"],
        ("C=C", "H-L"):["AE","ADE3"], 
        ("C=C", "Y-L"):["AE", "NuL"],
        ("Z=C", "H-L"):["proton_transfer","ADE3"],
        # This is one square on page 248 of the book "Electron Flow of Organic Chemistry", but we have two lines because SMARTS doesn't have "OR" 
        ("Y","ewg(O)-C=C"):["ADN"], 
        ("Y","ewg(N)-C=C"):["ADN"],
        # This is one square on page 248 of the book "Electron Flow of Organic Chemistry", but we have two lines because SMARTS doesn't have "OR" 
        ("C-","ewg(O)-C=C"):["ADN"],
        ("C-","ewg(N)-C=C"):["ADN"],

        ("Y","Z=C"):["ADN"],
        ("C-","Z=C"):["ADN"]
        }
