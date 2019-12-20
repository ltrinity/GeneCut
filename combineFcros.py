# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:28:19 2019
Script to add the output from an fcross script to an existing dictionary
@author: luket
"""
import os
import pickle
changesDict = pickle.load( open( "temp\\changesDict.p", "rb" ) )
