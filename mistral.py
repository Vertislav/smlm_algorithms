#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-

import json
import pandas as pd
import numpy as np

def pmtToPandas(pmt):
    dataFrames = {}
    for rank in pmt:
        columns = []
        for coefficient in pmt[rank]:
            if coefficient not in ['child', 'children']:
                columns.append((coefficient, pmt[rank][coefficient]))
#                print(columns)
        dataFrames[rank] = pd.DataFrame.from_items(columns)
    for rank in pmt:
        if 'child' in pmt[rank]:
            childRank = pmt[rank]['child']
            parent_ids =  []
            ids = []
            sizes = []
            for id, children in enumerate(pmt[rank]['children']):
                ids += children
                parent_ids += [id]*len(children)
                sizes.append(len(children))
            dataFrames[childRank][rank+'_id'] = pd.Series(parent_ids, index=ids)
            dataFrames[rank]['size'] = pd.Series(sizes)
                
    return dataFrames

def pandasToPmt(dataFrames):
    pmt = {}
    for key in dataFrames:
        rank = {}
        for coefficient in dataFrames[key]:
            if coefficient not in ['size', 'children'] and not coefficient.endswith('_id'):
                rank[coefficient] = list(dataFrames[key][coefficient].values)
        pmt[key] = rank
    for key in dataFrames:
        for coefficient in dataFrames[key]:
            if coefficient.endswith('_id'):
                children = [[]]*len(dataFrames[coefficient[:-3]])
                data = dataFrames[key].groupby([coefficient]).apply(lambda x: list(x.index)).to_dict()
                for id in data:
                    children[int(id)] = data[id]
                pmt[coefficient[:-3]]['children'] = children
                pmt[coefficient[:-3]]['child'] = key
    return pmt


def load(filename):
    with open(filename, 'r') as f:
        return pmtToPandas(json.load(f))

def dump(o, filename):

    class _NPEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                return super(_NPEncoder, self).default(obj)

    with open(filename, 'w') as f:
        json.dump(pandasToPmt(o), f, indent='\t', cls=_NPEncoder)
