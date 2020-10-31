#! /usr/bin/env python3
# coding: UTF-8

import matplotlib.pyplot as plt
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-input", type = open, help = "Input file (fasta format)")
parser.add_argument("-interval", type = int, default = 20, help = "Interval for calculating moving average (default: 20)")
parser.add_argument("-output", type = str, help = "Name for output file (If not specified, file is not created)")
args = parser.parse_args()

interval = args.interval
input = args.input.read()
output = args.output

hydropathyIndex = {"F":2.8,
"A":1.8,
"M":1.9,
"I":4.5,
"L":3.8,
"P":-1.6,
"V":4.2,
"W":-0.9,
"G":-0.4,
"S":-0.8,
"Y":-1.3,
"N":-3.5,
"Q":-3.5,
"T":-0.7,
"C":2.5,
"K":-3.9,
"R":-4.5,
"H":-3.2,
"D":-3.5,
"E":-3.5}

residualMw = {"F":147.174,
"A":71.078,
"M":131.196,
"I":113.158,
"L":113.158,
"P":97.115,
"V":99.131,
"W":186.210,
"G":57.051,
"S":87.077,
"Y":163.173,
"N":114.103,
"Q":128.129,
"T":101.104,
"C":103.143,
"K":128.172,
"R":156.186,
"H":137.139,
"D":115.087,
"E":129.114}

pKa = {"K":10.5,
"R":12.5,
"H":6.0,
"D":3.9,
"E":4.3,
"Y":10.1,
"C":8.3}

seqName = re.findall(">.*\n", input)[0][1:-1]
# remove sequence name, convert to CAPITAL, and remove characters other than amino acid
sequence = re.sub("[^FAMILPVWGSYNQTCKRHDE]", "", input.replace(seqName, "").upper())

hydropathyArray = []
hydropathy = 0

# calcultate initial value of moving average
for i in range(interval):
    hydropathy += hydropathyIndex[sequence[i]]/interval
hydropathyArray.append(hydropathy)

for i in range(1, len(sequence) - interval + 1):
    hydropathy = hydropathy - hydropathyIndex[sequence[i]]/interval + hydropathyIndex[sequence[i + interval - 1]]/interval
    hydropathyArray.append(hydropathy)

residueArray = range(int(interval/2), len(sequence) - interval + int(interval/2) + 1)

def getAminoAcidComposition(sequence):
    aaComposition = {}
    for char in "FAMILPVWGSYNQTCKRHDE":
        aaComposition[char] = sequence.count(char)
    return aaComposition

aaComposition = getAminoAcidComposition(sequence)

def getCharge(pH):
    # N,C terminal
    charge = 10**(-pH)/(10**(-pH) + 10**(-7.5)) - 10**(-3.5)/(10**(-pH) + 10**(-3.5))
    # Basic amino acids
    for char in "KRH":
        charge += 10**(-pH)*aaComposition[char]/(10**(-pH) + 10**(-pKa[char]))
    # Acidic amino acids
    for char in "DEYC":
        charge -= 10**(-pKa[char])*aaComposition[char]/(10**(-pH) + 10**(-pKa[char]))
    return charge

# calculate pI
def getPI():
    pI = 7
    if getCharge(pI) < 0:
        while getCharge(pI) < 0:
            pI -= 1
        pI += 1
        while getCharge(pI) < 0:
            pI -= 0.1
        pI += 0.1
        while getCharge(pI) < 0:
            pI -= 0.01
    else:
        while getCharge(pI) > 0:
            pI += 1
        pI -= 1
        while getCharge(pI) > 0:
            pI += 0.1
        pI += 0.1
        while getCharge(pI) > 0:
            pI += 0.01
    return pI

def getMw():
    mw = 18
    for char in "FAMILPVWGSYNQTCKRHDE":
        mw += residualMw[char]*aaComposition[char]
    return mw

averageHydropathy = 0
for char in "FAMILPVWGSYNQTCKRHDE":
    averageHydropathy += hydropathyIndex[char]*aaComposition[char]/len(sequence)

print("Length:",len(sequence))
print("Molecular weight:", round(getMw(),1))
print("pI:",round(getPI(),2))
print("Average of hydropathy:", round(averageHydropathy,2))
print("Amino acid composition")
for char in "FAMILPVWGSYNQTCKRHDE":
    print(char + ":", aaComposition[char])

plt.plot(residueArray, hydropathyArray)
plt.xlabel("Residue")
plt.ylabel("Hydropathy")

if output != None:
    plt.savefig(output)

plt.show()
