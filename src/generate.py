#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Marcus Ritt <marcus.ritt@inf.ufrgs.br>
"""

import os
import sys
import random
import numpy as np

random.seed(1)

## generate an instance set with operations in `os` and workers drawn from `ms`
def battarra(op,ms):
    for t in ['ap','mp']:
        for o in range(op[0],op[1]+1):
            for s in range(4):
                m = random.randint(ms[0],ms[1])
                print(f'./generate --n {o} --m {m} --seed {s} --type {t} > {t}-{m}-{o}-{s}.txt')
                os.system(f'./generate --n {o} --m {m} --seed {s} --type {t} > ../dat/setIIIa/{t}-{m}-{o}-{s}.txt')

## instance set I
#battarra([18,42],[9,15])

## instance set II
#battarra([43,67],[16,22])

## instance set III
battarra([68,92],[23,29])
