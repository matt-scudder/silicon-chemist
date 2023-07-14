"""
Base class that all reaction types extend.
Defines basic methods that all reaction types
should be able to execute, as well as data members.
"""

class Reaction(object):
    """
    A Reaction is defined as an interaction between a list of sources
    and a list of sinks. While usually only one source and one sink
    will be involved, this is left open for strange reaction types.

    A Reaction has a type, which is stored as a static string for
    all instances.

    A Reaction has a cross check, and a reversible rearrangement.

    A cross check is defined as a set of checks
    on the interaction between the source and the sink which produce
    a score between 0 and 1, usually involving the ΔpKa rule. A score
    of 0 means that the reaction cannot take place (and so should
    not be added to the list of possibilities), and a score of 1 means
    that the reaction is exceedingly likely to take place, and so should
    have high priority among all the possibilities.

    A rearrangement is defined as the shift in bonds between the source
    and the sink, resulting in a new molecular state. 

    Because of the fact that we need to keep discrete molecular state in order
    to list the final reaction pathway, as well as walk back and forth between it,
    a Reaction should take as input a #copy# of the current molecular state;
    the copy can be disposed of later on if "showing our work" is not required.

    All Reactions also include a method to convert themselves to a dictionary,
    which only really stores the cross check and a string denoting the Reaction's type,
    in order to make this class palatable by the browser and usable for teaching
    students.

    # "two_products" = True only for reactions that produce two equally products. This is whenthe Mark rule is equal for both carbons in Elimenations or Additions reactions
    """
    reaction_type = None
    def __init__(self, sources, sinks, second_product = False):  
        """
        Attaches sources and sinks to this object, and sets the initial cross check
        score to a sentinel value, so that cross_check() needs only to be run once.
        """
        self.sources = sources
        self.sinks = sinks
        self.cross_check_score = -2.0 
        self.second_product = second_product

    def __repr__(self):
        """
        Returns a representation of this object for easy debugging.
        """
        return "Reaction<Type:{},Cross-Check:{}>".format(self.reaction_type,self.cross_check_score)

    def cross_check(self):
        """
        Produces a score on a [0,1] scale based on how likely the source-sink
        interaction attached to this reaction is. Exceedingly unlikely (e.g. ΔpKa < -10)
        reactions should be given a score of 0, and exceedingly likely (e.g. ΔpKa > 10)
        reactions should be given a score of 1, though the latter is subject to change.
        """
        raise NotImplementedError("This is an abstract class and should not be used.")

    def rearrange(self):
        """
        Rearranges the atoms in the source and sink, in effect carrying out the reaction
        by making the connectivity changes necessary.
        """
        raise NotImplementedError("This is an abstract class and should not be used.")

    def to_json_dict(self):
        """
        Creates a dictionary that will be equivalent to a JSON representation of this object,
        for regeneration in the browser.
        """
        json_dict = {}
        json_dict["reaction_type"] = self.reaction_type
        json_dict["cross_check_score"] = self.cross_check_score
        return json_dict
