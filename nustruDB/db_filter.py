#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

df = pd.read_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/E_COLI_K12_02.csv', index_col=False)


df_reduced = df.drop_duplicates(subset=['primary_id','nucleotide_sequence'])

print(df_reduced.count())
