from warnings import filterwarnings
from matplotlib import MatplotlibDeprecationWarning
from anndata import ImplicitModificationWarning

filterwarnings('ignore', category=ImplicitModificationWarning)
filterwarnings('ignore', category=MatplotlibDeprecationWarning)

from FlowCytometryTools import FCPlate
from numpy import arcsinh

from glob import glob
from re import search

markers = ['Y','C','R']
channel_map = {
    'FJComp-B488 530_30-A' :'Y',
    'FJComp-V405 525_50-A' :'C',
    'FJComp-YG561 610_20-A':'R'
}

row_layout = {
    'A':25000, 'B':5000, 'C':1000, 'D':200,
    'E':40,    'F':8,    'G':1.6,  'H':0
}

col_layout = {
    '1':25000, '2':8333, '3':2777,'4':925,
    '5':308,   '6':102,  '7':34,  '8':11,
    '9':3.8,  '10':1.3, '11':0.4, '12':0
}

def parser(file_path) :
    index = search('_[A-Z][0-9]+_',file_path).group()

    col,row = index[2:-1],index[1]
    return (row_layout[row],col_layout[col])

plates = {'C6':[], 'C12':[]}
data_dirs = [ dir for dir in glob('*/*') if '.py' not in dir ]

for data_dir in data_dirs :
    primed,date = data_dir.split('/')

    plate = FCPlate.from_dir( position_mapper= lambda x: x,
            ID=date, path=data_dir, pattern='*.fcs', parser=parser,
            col_labels=list(col_layout.values()), row_labels = list(row_layout.values()))

    for idx in plate :
        
        plate[idx].data = plate[idx].data.rename(columns=channel_map)
        channel_differences = set(plate[idx].data.columns)-set(markers)

        if len(channel_differences) > 0 :
            raise ValueError('''{}\nChannel names not in markers. Change marker list or use channel_map'''.format(channel_differences))

        # map channel names
        plate[idx].data = plate[idx].data.reindex(columns=markers)
        plate[idx].meta['_channel_names_'] = tuple(markers)
        plate[idx].data = plate[idx].data[~plate[idx].data.isna()]

        # normalise to rfp channel
        plate[idx].data.Y /= plate[idx].data.R.abs()
        plate[idx].data.C /= plate[idx].data.R.abs()
        plate[idx].data = arcsinh(plate[idx].data)

    for idx in plate :

        # centre around autofluorescence signal
        plate[idx].data.Y -= plate[(0,0)].data.Y.mean()
        plate[idx].data.C -= plate[(0,0)].data.C.mean()
        plate[idx].data.R -= plate[(0,0)].data.R.mean()

    plates[primed] += [ plate ]