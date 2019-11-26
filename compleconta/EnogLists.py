#!/usr/bin/env python

class EnogList():
    """
    Object that stores the orthologous groups and their weights (if applicable)
    """

    def __init__(self, enog_list, enog_dict):
        """
        at initialization, the "EnogList" sorts the information that is required later. i.e. the dictionary of weights
        as used in completeness/contamination calculation. Additionally all used OGs (from the parameter enog_list) make
        up the total maximal weight (self.total)

        :param enog_list: a list of orthoglogous groups
        :param enog_dict: a dictionary containing all information from the weights file per orthologous group
        """

        self.weights={}
        self.enogs=enog_list[:]

        self.total=0

        for enog in self.enogs:
            if not enog_dict.get(enog):
                self.weights[enog] = 1
            else:
                percent_presence=enog_dict[enog].get("%present", 1)
                average_count=enog_dict[enog].get("av.count_if_present", 1)
                self.weights[enog]=float(percent_presence)/float(average_count)
            self.total=self.total+self.weights[enog]


    def get_weight(self, enog):
        """

        :param enog: id of the orthologous group
        :return: specific weight for the orthologous group id
        """
        weight=self.weights.get(enog, 0)

        return weight

    def get_total(self):
        """
        :return: total maximal score that can be reached (sum of weights)
        """

        return self.total 

    def get_dict(self):
        """
        :return: dictionary of weights as calculated
        """
        return self.weights
