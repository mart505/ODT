#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 23:23:34 2020

Dataprocessing script for ODT.py
This script reads data and formats it into the right format for ODT and 
provides the tree topology variables

@author: martijnkorf
"""

# stores all the different data instances per set
import random
from collections import namedtuple
import math
import numpy
import itertools


def readData(flag,setting, maxRank,imbalanceCor):
    
    splitter = ","
    
    #readfile
    if(flag == "mush") :
        datafile = "agaricus_lepiota.txt"
        features = ["b","c","x","f","k","s",
                    "f","g","y","s",
                    "n","b","c","g","r","p","u","e","w","y",
                    "t","f",
                    "a","l","c","y","f","m","n","p","s",
                    "a","d","f","n",
                    "c","w","d",
                    "b","n",
                    "k","n","b","h","g","r","o","p","u","e","w","y",
                    "e","t",
                    "b","c","u","e","z","r","?",
                    "f","y","k","s",
                    "f","y","k","s",
                    "n","b","c","g","o","p","e","w","y",
                    "n","b","c","g","o","p","e","w","y",
                    "p","u",
                    "n","o","w","y",
                    "n","o","t",
                    "c","e","f","l","n","p","s","z",
                    "k","n","b","h","r","o","u","w","y",
                    "a","c","n","s","v","y",
                    "g","l","m","p","u","w","d"]
        G = 22
        J = len(features)
        placeYinData = 0
        yPositiveValue = "e"
        C = 1.083333333
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 6   : g = 2  ; numerical = False ; max1= []
            if j == 10  : g = 3  ; numerical = False ; max1= []
            if j == 20  : g = 4  ; numerical = False ; max1= []
            if j == 22  : g = 5  ; numerical = False ; max1= []
            if j == 31  : g = 6  ; numerical = False ; max1= []
            if j == 35  : g = 7  ; numerical = False ; max1= []
            if j == 38  : g = 8  ; numerical = False ; max1= []
            if j == 40  : g = 9  ; numerical = False ; max1= []
            if j == 52  : g = 10 ; numerical = False ; max1= []
            if j == 54  : g = 11 ; numerical = False ; max1= []
            if j == 61  : g = 12 ; numerical = False ; max1= []
            if j == 65  : g = 13 ; numerical = False ; max1= []
            if j == 69  : g = 14 ; numerical = False ; max1= []
            if j == 78  : g = 15 ; numerical = False ; max1= []
            if j == 87  : g = 16 ; numerical = False ; max1= []
            if j == 89  : g = 17 ; numerical = False ; max1= []
            if j == 93  : g = 18 ; numerical = False ; max1= []
            if j == 96  : g = 19 ; numerical = False ; max1= []
            if j == 104 : g = 20 ; numerical = False ; max1= []
            if j == 112 : g = 21 ; numerical = False ; max1= []
            if j == 119 : g = 22 ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if(flag == "mush2") :
        datafile = "agaricus_lepiota2"
        features = ["b","c","x","f","k","s",
                    "f","g","y","s",
                    "n","b","c","g","r","p","u","e","w","y",
                    "t","f",
                    "a","l","c","y","f","m","n","p","s",
                    
                    "a","f", # d, n
                    "c","w", # d
                    "b","n",
                    "k","n","b","h","g","r","o","p","u","e","w","y",
                    "e","t",
                    
                    "b","c","e","r","?", # u,z
                    "f","y","k","s",
                    "f","y","k","s",
                    "n","b","c","g","o","p","e","w","y",
                    "n","b","c","g","o","p","e","w","y",
                     #u, which makes this group redundant
                    "n","o","w","y",
                    "n","o","t",
                    "e","f","l","n","p", #c,s,z
                    "k","n","b","h","r","o","u","w","y",
                    "a","c","n","s","v","y",
                    "g","l","m","p","u","w","d"] # g wel g
        G = 21
        J = len(features)
        placeYinData = 0
        yPositiveValue = "e"
        splitter = ";"
        C = 1.083333333
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 6   : g = 2  ; numerical = False ; max1= []
            if j == 10  : g = 3  ; numerical = False ; max1= []
            if j == 20  : g = 4  ; numerical = False ; max1= []
            if j == 22  : g = 5  ; numerical = False ; max1= []
            if j == 31  : g = 6  ; numerical = False ; max1= []
            if j == 33  : g = 7  ; numerical = False ; max1= []
            if j == 35  : g = 8  ; numerical = False ; max1= []
            if j == 37  : g = 9  ; numerical = False ; max1= []
            if j == 49  : g = 10 ; numerical = False ; max1= []
            if j == 51  : g = 11 ; numerical = False ; max1= []
            if j == 56  : g = 12 ; numerical = False ; max1= []
            if j == 60  : g = 13 ; numerical = False ; max1= []
            if j == 64  : g = 14 ; numerical = False ; max1= []
            if j == 73  : g = 15 ; numerical = False ; max1= []
            if j == 82  : g = 16 ; numerical = False ; max1= []
            if j == 86  : g = 17 ; numerical = False ; max1= []
            if j == 89  : g = 18 ; numerical = False ; max1= []
            if j == 94  : g = 19 ; numerical = False ; max1= []
            if j == 103 : g = 20 ; numerical = False ; max1= []
            if j == 109 : g = 21 ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "ttt") :
        datafile = "tic-tac-toe.data.txt"
        features = ["x","o","b",
                    "x","o","b",
                    "x","o","b",
                    "x","o","b",
                    "x","o","b",
                    "x","o","b",
                    "x","o","b",
                    "x","o","b",
                    "x","o","b"]
        G = 9
        J = len(features)
        placeYinData = 9
        yPositiveValue = "positive"
        C = 1.85714
        
        features2 = []
        for j in range(J) :
            if j == 0  : g = 1 ; numerical = False ; max1= []
            if j == 3  : g = 2 ; numerical = False ; max1= []
            if j == 6  : g = 3 ; numerical = False ; max1= []
            if j == 9  : g = 4 ; numerical = False ; max1= []
            if j == 12 : g = 5 ; numerical = False ; max1= []
            if j == 15 : g = 6 ; numerical = False ; max1= []
            if j == 18 : g = 7 ; numerical = False ; max1= []
            if j == 21 : g = 8 ; numerical = False ; max1= []
            if j == 24 : g = 9 ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "heloc"):
        datafile = "Heloc_data"
        features = ["q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10"]
        J = len(features)
        G = 23
        placeYinData = 0
        yPositiveValue = "Good"
        splitter = ";"
        C = 0.612903
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = True ; max1= "89"
            if j == 10  : g = 2  ; numerical = True ; max1= "604"
            if j == 20  : g = 3  ; numerical = True ; max1= "65"
            if j == 30  : g = 4  ; numerical = True ; max1= "224"
            if j == 40  : g = 5  ; numerical = True ; max1= "78"
            if j == 50  : g = 6  ; numerical = True ; max1= "16"
            if j == 60  : g = 7  ; numerical = True ; max1= "16"
            if j == 70  : g = 8  ; numerical = True ; max1= "99"
            if j == 80  : g = 9  ; numerical = True ; max1= "83"
            if j == 90  : g = 10 ; numerical = True ; max1= "7"
            if j == 100 : g = 11 ; numerical = True ; max1= "6"
            if j == 110 : g = 12 ; numerical = True ; max1= "100"
            if j == 120 : g = 13 ; numerical = True ; max1= "19"
            if j == 130 : g = 14 ; numerical = True ; max1= "92"
            if j == 140 : g = 15 ; numerical = True ; max1= "24"
            if j == 150 : g = 16 ; numerical = True ; max1= "24"
            if j == 160 : g = 17 ; numerical = True ; max1= "24"
            if j == 170 : g = 18 ; numerical = True ; max1= "154"
            if j == 180 : g = 19 ; numerical = True ; max1= "153"
            if j == 190 : g = 20 ; numerical = True ; max1= "32"
            if j == 200 : g = 21 ; numerical = True ; max1= "23"
            if j == 210 : g = 22 ; numerical = True ; max1= "12"
            if j == 220 : g = 23 ; numerical = True ; max1= "100"
            features2.append([features[j], g, numerical, max1])
        
    if (flag == "votes"):  
        datafile = "house-votes-84.data"
        features = ["y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?",
                    "y","n","?"]
        J = len(features)
        G = 16
        placeYinData = 0
        yPositiveValue = "democrat"
        C = 1.564102
        
        features2 = []
        for j in range(J) :
            if j == 0  : g = 1  ; numerical = False ; max1= []
            if j == 3  : g = 2  ; numerical = False ; max1= []
            if j == 6  : g = 3  ; numerical = False ; max1= []
            if j == 9  : g = 4  ; numerical = False ; max1= []
            if j == 12 : g = 5  ; numerical = False ; max1= []
            if j == 15 : g = 6  ; numerical = False ; max1= []
            if j == 18 : g = 7  ; numerical = False ; max1= []
            if j == 21 : g = 8  ; numerical = False ; max1= []
            if j == 24 : g = 9  ; numerical = False ; max1= []
            if j == 27 : g = 10 ; numerical = False ; max1= []
            if j == 30 : g = 11 ; numerical = False ; max1= []
            if j == 33 : g = 12 ; numerical = False ; max1= []
            if j == 36 : g = 13 ; numerical = False ; max1= []
            if j == 39 : g = 14 ; numerical = False ; max1= []
            if j == 42 : g = 15 ; numerical = False ; max1= []
            if j == 45 : g = 16 ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
    
    if (flag == "student") :
        datafile = "student_mat.data"
        features = ["GP","MS",
                    "F","M",
                    "15","16","17","18","19","20","21","22",
                    "U","R",
                    "LE3","GT3",
                    "T","A",
                    "0","1","2","3","4",
                    "0","1","2","3","4",
                    "teacher","health","services","at_home","other",
                    "teacher","health","services","at_home","other",
                    "home","reputation","course","other",
                    "mother","father","other",
                    "1","2","3","4",
                    "1","2","3","4",
                    "0","1","2","3","4",
                    "yes","no",
                    "yes","no",
                    "yes","no",
                    "yes","no",
                    "yes","no",
                    "yes","no",
                    "yes","no",
                    "yes","no",
                    "1","2","3","4","5",
                    "1","2","3","4","5",
                    "1","2","3","4","5",
                    "1","2","3","4","5",
                    "1","2","3","4","5",
                    "1","2","3","4","5",
                    "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10"]
        J = len(features)
        G = 30
        placeYinData = 32
        yPositiveValue = 10 # bigger than 10
        splitter = ";"
        C = 2.030303
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 2   : g = 2  ; numerical = False ; max1= []
            if j == 4   : g = 3  ; numerical = True  ; max1= []
            if j == 12  : g = 4  ; numerical = False ; max1= []
            if j == 14  : g = 5  ; numerical = False ; max1= []
            if j == 16  : g = 6  ; numerical = False ; max1= []
            if j == 18  : g = 7  ; numerical = True  ; max1= []
            if j == 23  : g = 8  ; numerical = True  ; max1= []
            if j == 28  : g = 9  ; numerical = False ; max1= []
            if j == 33  : g = 10 ; numerical = False ; max1= []
            if j == 38  : g = 11 ; numerical = False ; max1= []
            if j == 42  : g = 12 ; numerical = False ; max1= []
            if j == 45  : g = 13 ; numerical = True  ; max1= []
            if j == 49  : g = 14 ; numerical = True  ; max1= []
            if j == 53  : g = 15 ; numerical = True  ; max1= []
            if j == 58  : g = 16 ; numerical = False ; max1= []
            if j == 60  : g = 17 ; numerical = False ; max1= []
            if j == 62  : g = 18 ; numerical = False ; max1= []
            if j == 64  : g = 19 ; numerical = False ; max1= []
            if j == 66  : g = 20 ; numerical = False ; max1= []
            if j == 68  : g = 21 ; numerical = False ; max1= []
            if j == 70  : g = 22 ; numerical = False ; max1= []
            if j == 72  : g = 23 ; numerical = False ; max1= []
            if j == 74  : g = 24 ; numerical = True  ; max1= [] 
            if j == 79  : g = 25 ; numerical = True  ; max1= []
            if j == 84  : g = 26 ; numerical = True  ; max1= []
            if j == 89  : g = 27 ; numerical = True  ; max1= []
            if j == 94  : g = 28 ; numerical = True  ; max1= []
            if j == 99  : g = 29 ; numerical = True  ; max1= []
            if j == 104 : g = 30 ; numerical = True  ; max1= "93" 
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "bc") :
        datafile = "bc.data"
        features = ["1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    ]
        J = len(features)
        G = 9
        placeYinData = 0
        yPositiveValue = "2"
        splitter = ";"
        C = 1.8571429
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = True ; max1= []
            if j == 10  : g = 2  ; numerical = True ; max1= []
            if j == 20  : g = 3  ; numerical = True ; max1= []
            if j == 30  : g = 4  ; numerical = True ; max1= []
            if j == 40  : g = 5  ; numerical = True ; max1= []
            if j == 50  : g = 6  ; numerical = True ; max1= []
            if j == 60  : g = 7  ; numerical = True ; max1= []
            if j == 70  : g = 8  ; numerical = True ; max1= []
            if j == 80  : g = 9  ; numerical = True ; max1= []
            features2.append([features[j], g, numerical, max1])
    
    if (flag == "a1a") :
        datafile = "adults.data"
        features = ["q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                    
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    "1","2","3","4","5","6","7","8","9","10",
                    ]
        J = len(features)
        G = 9
        placeYinData = 0
        yPositiveValue = 2
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = True ; max1= []
            if j == 10  : g = 2  ; numerical = True ; max1= []
            if j == 20  : g = 3  ; numerical = True ; max1= []
            if j == 30  : g = 4  ; numerical = True ; max1= []
            if j == 40  : g = 5  ; numerical = True ; max1= []
            if j == 50  : g = 6  ; numerical = True ; max1= []
            if j == 60  : g = 7  ; numerical = True ; max1= []
            if j == 70  : g = 8  ; numerical = True ; max1= []
            if j == 80  : g = 9  ; numerical = True ; max1= []
            features2.append([features[j], g, numerical, max1])
    
    if (flag == "heart") :
        datafile = "SPECT.data"
        features = ["1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0",
                    "1","0"]
        J = len(features)
        G = 22
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ","
        C = 3.7619048
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 2   : g = 2  ; numerical = False ; max1= []
            if j == 4   : g = 3  ; numerical = False ; max1= []
            if j == 6   : g = 4  ; numerical = False ; max1= []
            if j == 8   : g = 5  ; numerical = False ; max1= []
            if j == 10  : g = 6  ; numerical = False ; max1= []
            if j == 12  : g = 7  ; numerical = False ; max1= []
            if j == 14  : g = 8  ; numerical = False ; max1= []
            if j == 16  : g = 9  ; numerical = False ; max1= []
            if j == 18  : g = 10 ; numerical = False ; max1= []
            if j == 20  : g = 11 ; numerical = False ; max1= []
            if j == 22  : g = 12 ; numerical = False ; max1= []
            if j == 24  : g = 13 ; numerical = False ; max1= []
            if j == 26  : g = 14 ; numerical = False ; max1= []
            if j == 28  : g = 15 ; numerical = False ; max1= []
            if j == 30  : g = 16 ; numerical = False ; max1= []
            if j == 32  : g = 17 ; numerical = False ; max1= []
            if j == 34  : g = 18 ; numerical = False ; max1= []
            if j == 36  : g = 19 ; numerical = False ; max1= []
            if j == 38  : g = 20 ; numerical = False ; max1= []
            if j == 40  : g = 21 ; numerical = False ; max1= []
            if j == 42  : g = 22 ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "monks-1") :
        datafile = "monks-1.data"
        features = ["1","2","3",
                    "1","2","3",
                    "1","2",
                    "1","2","3",
                    "1","2","3","4",
                    "1","2"]
        J = len(features)
        G = 6
        placeYinData = 0
        yPositiveValue = "1"
        splitter = " "
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 3   : g = 2  ; numerical = False ; max1= []
            if j == 6   : g = 3  ; numerical = False ; max1= []
            if j == 8   : g = 4  ; numerical = False ; max1= []
            if j == 11  : g = 5  ; numerical = False ; max1= []
            if j == 15  : g = 6  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "monks-12") :
        datafile = "monks-1_2.data.csv"
        features = ["1","2","3","4"]
        J = len(features)
        G = 1
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "monks-13") :
        datafile = "monks-1_3.data.csv"
        features = ["1","2","3",
                    "1","2","3","4"]
        J = len(features)
        G = 2
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 3   : g = 2  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "monks-14") :
        datafile = "monks-1_4.data.csv"
        features = ["1","2","3",
                    "1","2","3",
                    "1","2","3",
                    "1","2","3","4"
                    ]
        J = len(features)
        G = 4
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 3   : g = 2  ; numerical = False ; max1= []
            if j == 6   : g = 3  ; numerical = False ; max1= []
            if j == 9   : g = 4  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
    
    if (flag == "monks-15") :
        datafile = "monks-1_5.data.csv"
        features = ["1","2","3",
                    "1","2","3",
                    "1","2","3",
                    "1","2","3","4",
                    "1","2"]
        J = len(features)
        G = 5
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 3   : g = 2  ; numerical = False ; max1= []
            if j == 6   : g = 3  ; numerical = False ; max1= []
            if j == 9   : g = 4  ; numerical = False ; max1= []
            if j == 13  : g = 5  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
    
    if (flag == "monks-16") :
        datafile = "monks-1_6.data.csv"
        features = ["1","2","3","4",
                    "1","2"]
        J = len(features)
        G = 2
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 4   : g = 2  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
        
    if (flag == "monks-17") :
        datafile = "monks-1_7.data.csv"
        features = ["1","2","3",
                    "1","2","3","4",
                    "1","2"]
        J = len(features)
        G = 3
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 3   : g = 2  ; numerical = False ; max1= []
            if j == 7   : g = 3  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
            
    if (flag == "monks-18") :
        datafile = "monks-1_8.data.csv"
        features = ["1","2","3",
                    "1","2","3",
                    "1","2",
                    "1","2","3",
                    "1","2"]
        J = len(features)
        G = 5
        placeYinData = 0
        yPositiveValue = "1"
        splitter = ";"
        C = 1
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 3   : g = 2  ; numerical = False ; max1= []
            if j == 6   : g = 3  ; numerical = False ; max1= []
            if j == 8   : g = 4  ; numerical = False ; max1= []
            if j == 11  : g = 5  ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
        
    if (flag == "krkp") :
        datafile = "kr-vs-kp.data"
        features = ["f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "g","l",
                    "f","t",
                    "b","n","w",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "f","t",
                    "n","t"]
        J = len(features)
        G = 36
        placeYinData = 36 # so place 37
        yPositiveValue = "won"
        splitter = ","
        C = 1.0833333
        
        features2 = []
        for j in range(J) :
            if j == 0   : g = 1  ; numerical = False ; max1= []
            if j == 2   : g = 2  ; numerical = False ; max1= []
            if j == 4   : g = 3  ; numerical = False ; max1= []
            if j == 6   : g = 4  ; numerical = False ; max1= []
            if j == 8   : g = 5  ; numerical = False ; max1= []
            if j == 10  : g = 6  ; numerical = False ; max1= []
            if j == 12  : g = 7  ; numerical = False ; max1= []
            if j == 14  : g = 8  ; numerical = False ; max1= []
            if j == 16  : g = 9  ; numerical = False ; max1= []
            if j == 18  : g = 10 ; numerical = False ; max1= []
            if j == 20  : g = 11 ; numerical = False ; max1= []
            if j == 22  : g = 12 ; numerical = False ; max1= []
            if j == 24  : g = 13 ; numerical = False ; max1= []
            if j == 26  : g = 14 ; numerical = False ; max1= []
            if j == 28  : g = 15 ; numerical = False ; max1= []
            if j == 31  : g = 16 ; numerical = False ; max1= []
            if j == 33  : g = 17 ; numerical = False ; max1= []
            if j == 35  : g = 18 ; numerical = False ; max1= []
            if j == 37  : g = 19 ; numerical = False ; max1= []
            if j == 39  : g = 20 ; numerical = False ; max1= []
            if j == 41  : g = 21 ; numerical = False ; max1= []
            if j == 43  : g = 22 ; numerical = False ; max1= []
            if j == 45  : g = 23 ; numerical = False ; max1= []
            if j == 47  : g = 24 ; numerical = False ; max1= []
            if j == 49  : g = 25 ; numerical = False ; max1= []
            if j == 51  : g = 26 ; numerical = False ; max1= []
            if j == 53  : g = 27 ; numerical = False ; max1= []
            if j == 55  : g = 28 ; numerical = False ; max1= []
            if j == 57  : g = 29 ; numerical = False ; max1= []
            if j == 59  : g = 30 ; numerical = False ; max1= []
            if j == 61  : g = 31 ; numerical = False ; max1= []
            if j == 63  : g = 32 ; numerical = False ; max1= []
            if j == 65  : g = 33 ; numerical = False ; max1= []
            if j == 67  : g = 34 ; numerical = False ; max1= []
            if j == 69  : g = 35 ; numerical = False ; max1= []
            if j == 71  : g = 36 ; numerical = False ; max1= []
            features2.append([features[j], g, numerical, max1])
    #-------------------------------------------------------------------------
    #general code for all data sets
    #-------------------------------------------------------------------------
    aName = namedtuple( "features", ["data","group","numerical","max1"])
    featuresNamed = [aName(*ar) for ar in features2]
    
    f = open(datafile,"r")
    yourResults = [line.split(splitter) for line in f.readlines()]
    N = len(yourResults)
    for i in range(N): # get rid of \n in last item in each sample
       yourResults[i][-1] = yourResults[i][G].replace("\n","")
       if (flag == "student") :
           for g in range(G+3) :
              yourResults[i][g] = yourResults[i][g].strip('"')
    random.shuffle(yourResults)
    
    #creation of ysamples (binary vector)
    if (flag == "student") :
        ysamples = [ int(int(yourResults[i][placeYinData]) >= yPositiveValue) for i in range(N)]
    else :
        ysamples = [ int(yourResults[i][placeYinData] == yPositiveValue) for i in range(N)]
        
    # creation of aMatrix (binary matrix)
    aMatrix = []
    for i in range(N) :
        
        # deals with different places of y in the dataset
        if (placeYinData == 0) : g = 1 ; lastGroup = 1 ; thisGroup = 0
        else :                   g = 0 ; lastGroup = 0 ; thisGroup = 1
          
        # converting categorical yourResults row into binary aMatrix row
        newline = []
        for j in range(J) :
            if (g == G+lastGroup) : newline.append(0) # last
            elif (len(featuresNamed[j].max1) > 0) & (featuresNamed[j].group == g+thisGroup) : # numerical variables /wout grouping
                step = int(featuresNamed[j].max1) / 10
                var = math.ceil(int(int(yourResults[i][g]))/ step)
                for l in range(10):
                    newline.append(1) if (l == var-1) or ((var == 0) & (l==0)) else newline.append(0)
                g += 1
                if (g == G+lastGroup) : break
            elif (yourResults[i][g] == featuresNamed[j].data) & (featuresNamed[j].group == g+thisGroup): #categorical variables or numerical with grouping
                newline.append(1)
                g += 1
            elif (len(featuresNamed[j].max1) == 0) : newline.append(0)
        aMatrix.append(newline)
        
    # max rank extension
    if (maxRank == True) :
        print("constructing max rank matrix")
        WIPMatrix = [] ; indexRank = [] ; rankMatrix = 0
        restMatrix = [] ; indexRest = []
        
        for i in range(len(aMatrix)):
            WIPMatrix.append(aMatrix[i])
            newRank = numpy.linalg.matrix_rank(WIPMatrix)
            if (newRank > rankMatrix) :
                rankMatrix = newRank
                indexRank.append(i)
            else :
                WIPMatrix.pop()
                restMatrix.append(aMatrix[i])
                indexRest.append(i)
        completeIndex = indexRank+indexRest
        aMatrix = WIPMatrix + restMatrix 
        #also reorder ysamples
        if (flag == "student") :
            ysamples = [ int(int(yourResults[i][placeYinData]) >= yPositiveValue) for i in completeIndex]
        else :
            ysamples = [ int(yourResults[i][placeYinData] == yPositiveValue) for i in completeIndex]

    #set of groups with numerical features
    GNumerical = set((g for g in range(G) if featuresNamed[j].numerical))
    # need lastjperg for numerical constraints
    lastjperg = [ j for j in range(J-1) if (featuresNamed[j].group != featuresNamed[j+1].group)]
    lastjperg.append(J-1)

    #anchoring array
    anchorFeatures = [j for j in range(J) if (featuresNamed[j].group != featuresNamed[j-1].group) ]
    
    # caculate A Matrix sparsity
    counter0 = sum(1 for i,j in itertools.product(range(N), range(J)) if aMatrix[i][j] == 0)
    datasetSparsity = counter0 / (J*N)
    
    # set training sample size and determine Iplus and Imin
    if (setting >= N) : print("error"); return
    Ntrain = min(int(N*0.9),600) if (setting == "90%/600") else setting
    Iplus = [i for i in range(Ntrain) if ysamples[i] == 1]
    Imin =  [i for i in range(Ntrain) if ysamples[i] == 0]
    
    if (imbalanceCor): C = 1

    return features, featuresNamed, anchorFeatures, ysamples, J, G, aMatrix, yourResults, GNumerical, lastjperg, datasetSparsity, C, Iplus, Imin, Ntrain

# function to generate variables for each tree topology
def Topology(treeTopology):
    
    topologies = ["depth2","depth2.5","depth3","imbalanced"]
    indexTop = topologies.index(treeTopology)
    
    KAll = [3,5,7,7]
    BAll = [4,6,8,8]
    KnoleavesALL = [{0},{0,1},{0,1,4},{0,1,2}]
    KleavesAll = [{1,2},{2,3,4},{2,3,5,6},{3,4,5,6}]
    KleftAll   = [[{0,1},{0},{2},{}],[{0,1,2},{0,1},{0,3},{0},{4},{}],
                  [{0,1,2},{0,1},{0,3},{0},{4,5},{4},{6},{}],
                  [{0,1,2,3},{0,1,2},{0,1,4},{0,1},{0,5},{0},{6},{}]]           
    KrightAll  = [[{},{1},{0},{0,2}],[{},{2},{1},{1,3},{0},{0,4}],
                  [{},{2},{1},{1,3},{0},{0,5},{0,4},{0,4,6}],
                  [{},{3},{2},{2,4},{1},{1,5},{0},{0,6}]]
    BplusAll   = [[0,2], [0,2,4], [0,2,4,6], [0,2,4,6]]
    BminAll    = [[1,3], [1,3,5], [1,3,5,7], [1,3,5,7]]
    BleftAll   = [[{0,1},{0},{2}],[{0,1,2,3},{0,1},{0},{2},{4}],
                  [{0,1,2,3},{0,1},{0},{2},{4,5},{4},{6}],
                  [{0,1,2,3,4,5},{0,1,2,3},{0,1},{0},{2},{4},{6}]]
    BrightAll  = [[{2,3},{1},{3}],[{4,5},{2,3},{1},{3},{5}],
                  [{4,5,6,7},{2,3},{1},{3},{5,7},{5},{7}],
                  [{6,7},{4,5},{2,3},{1},{3},{5},{7}]]
    
    K,B               = KAll[indexTop],BAll[indexTop]
    Knoleaves,Kleaves = KnoleavesALL[indexTop],KleavesAll[indexTop]
    Kleft,Kright      = KleftAll[indexTop],KrightAll[indexTop]
    Bplus,Bmin        = BplusAll[indexTop],BminAll[indexTop]
    Bleft,Bright      = BleftAll[indexTop],BrightAll[indexTop]
    
    return B, K, Bplus, Bmin, Kleft, Kright, Knoleaves, Kleaves, Bright, Bleft, treeTopology
