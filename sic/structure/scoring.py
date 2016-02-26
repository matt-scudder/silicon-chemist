#!/usr/bin/python
#coding=utf-8
import math
def score_pka(limit,k,dpka):
    """
    Makes a score for pKa on a [0,1] scale.
    Limit is the cutoff point for your pKa scale (10 for proton transfers, 6 for nucleophilicity), at which point
    the result of the function should be 1. This can be a negative number if you want to flip which "side"
    of the number line the 1 is on.
    k is a constant, usually defining your midpoint. The function at 0 will not be exactly k, but will be
    pretty close.
    dpka is the difference in pKa that you want to test.
    """
    return 1.0 / (1.0 + k * (math.exp((6.0 / limit)*-dpka)))
