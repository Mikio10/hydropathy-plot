#! /usr/bin/env python3
# coding: UTF-8

import matplotlib.pyplot as plt
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-input", type = open, help = "Input file (fasta format)")
parser.add_argument("-interval", type = int, default = 20, help = "Interval for calculating moving average (default: 20)")
parser.add_argument("-output", type = str, help = "Name for output file (File is not created if not specified)")
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

seqName = re.findall(">.*\n", input)[0][1:-1]
# remove sequence name, convert to CAPITAL, and remove characters other than amino acid
sequence = re.sub("[^FAMILPVWGSYNQTCKRHDE]", "", input.replace(seqName, "").upper())

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

plt.plot(midResidue, hydropathy)
plt.xlabel("Residue")
plt.ylabel("Hydropathy")

if output != None:
    plt.savefig(output)

plt.show()
