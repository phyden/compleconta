#!/usr/bin/env python

import sys, os
#import numpy as np
#import numpy.random as random

class EnogList():

    def __init__(self, enog_list, enog_dict):
    
        self.weights={}
        self.enogs=enog_list[:]

        self.total=0

        for enog in self.enogs:
            percent_presence=enog_dict[enog]["%present"]
            average_count=enog_dict[enog]["av.count_if_present"]
            self.weights[enog]=float(percent_presence)/float(average_count)
            self.total=self.total+self.weights[enog]


    def get_weight(self, enog):
        
        try:
            weight=self.weights[enog]
        except KeyError:
            weight=0
        return weight

    def get_total(self):

        return self.total 

    def get_dict(self):

        return self.weights
