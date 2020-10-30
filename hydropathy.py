#! /usr/bin/env python3
# coding: UTF-8

import matplotlib.pyplot as plt
import argparse
import re
import math

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

residualMw = {"F":147.1,
"A":71.0,
"M":131.0,
"I":113.1,
"L":113.1,
"P":97.1,
"V":99.1,
"W":186.1,
"G":57.0,
"S":87.0,
"Y":163.1,
"N":114.0,
"Q":128.1,
"T":101.0,
"C":103.0,
"K":128.1,
"R":156.1,
"H":137.1,
"D":115.0,
"E":129.0}

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

aaComposition = {}
for char in "FAMILPVWGSYNQTCKRHDE":
    aaComposition[char] = sequence.count(char)

hydropathy = []
movingAverage = 0

# calcultate initial value of moving average
for i in range(interval):
    movingAverage += hydropathyIndex[sequence[i]]/interval
hydropathy.append(movingAverage)

for i in range(1, len(sequence) - interval + 1):
    movingAverage = movingAverage - hydropathyIndex[sequence[i]]/interval + hydropathyIndex[sequence[i + interval - 1]]/interval
    hydropathy.append(movingAverage)

midResidue = range(int(interval/2), len(sequence) - interval + int(interval/2) + 1)

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
pI = 7
if getCharge(pI) < 0:
    while getCharge(pI) < 0:
        pI -= 0.01
else:
    while getCharge(pI) > 0:
        pI += 0.01

molecularWeight = 18
for char in "FAMILPVWGSYNQTCKRHDE":
    molecularWeight += residualMw[char]*aaComposition[char]

avgHydropathy = 0
for char in "FAMILPVWGSYNQTCKRHDE":
    avgHydropathy += hydropathyIndex[char]*aaComposition[char]/len(sequence)

print("Length:",len(sequence))
print("Molecular weight:", round(molecularWeight,1))
print("pI:",round(pI,2))
print("Average of hydropathy:", round(avgHydropathy,2))
print("Amino acid composition")
for char in "FAMILPVWGSYNQTCKRHDE":
    print(char + ":", aaComposition[char])

plt.plot(midResidue, hydropathy)
plt.xlabel("Residue")
plt.ylabel("Hydropathy")

if output != None:
    plt.savefig(output)

plt.show()
