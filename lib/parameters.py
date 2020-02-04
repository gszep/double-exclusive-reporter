from __future__ import print_function

def parse_mle(file_path):
    store = {}
    file = open(file_path, 'r')
    lines = file.readlines()
    for line in lines:
        kws = line.replace('e76','eps76').replace('e81','eps81').split(' ')
        if (kws[0] == 'parameterSummary') & kws[2].startswith('{mle='):
            v = float(kws[2].replace('{mle=','').replace(';',''))
            store[kws[1]] = v
    return store