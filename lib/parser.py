from re import sub,search,finditer
from yaml import load

from .model import Model
from .utils import isnumber


def fromcrn(file_path) :
    '''Creates model object from parameters in a given crn file

    --- parameters ---
    file_path : <str>
        path to crn file to read parameters from

    --- returns ---
    model : <Model>
        model object initialised with parameters
    '''

    # open and read file as string
    kwargs = { }
    with open(file_path, 'r') as file:

        # read and remove comments
        file_string = file.read()
        file_string = remove_comments(file_string)

        # parse rules separated by pipes
        kwargs['reactions'] = parse_delimited(file_string,contains='->')
        kwargs['inits'] = parse_delimited(file_string,contains='init')

        # convert to yaml compatible syntax
        file_string = file_string.replace(']','}')
        file_string = file_string.replace('[','{')
        file_string = file_string.replace(';',',')
        file_string = file_string.replace('=',':')

        # find assigned parameters
        pattern = r'([0-9|\w]+)[ ]*:[ ]*[\w0-9.eE+-]+'
        for match in finditer(pattern,file_string) :

                item = match.group()
                name,value = item.split(':')

                # format as yaml
                value = value.strip() if isnumber(value) else '"'+value.strip()+'"'
                yaml_format = '"'+name.strip()+'":'+value
                file_string = file_string.replace(item,yaml_format)

        # parse each directive to dictionary
        pattern = r'(?<=directive ).*?(?=(directive|init|\|))'
        for match in finditer(pattern,file_string.replace('\n',' ')) :
            item = match.group()

            # attempt to parse automagically
            if 'simulator' not in item :
                key,value = item.split(" ",1)
                try : kwargs[key] = load(value)

                # otherwise parse manually
                except :

                    kwargs[key] = {}
                    for s in value[1:-1].split(',') :
                        name,value = s.strip().split(':')

                        # convert to python compatible syntax
                        value = value.replace('^','**')

                        value = value.replace('"','').replace('{','').replace('}','')
                        name = name.replace('"','').replace('{','').replace('}','')

                        value = float(value) if isnumber(value) else value.strip()
                        kwargs[key][name.strip()] = value

        # parse initial conditions
        kwargs['init'] = {}
        pattern = r'init[ ]*\w+[ ]*[0-9.eE+-]+'
        for match in finditer(pattern,file_string) :

                item = match.group()
                _,name,value = item.split(' ')
                kwargs['init'][name.strip()] = float(value)

    return Model(**kwargs)


def remove_comments(file_string,comment_symbol='//'):
    '''removes commented lines from string'''
    return sub(r'({}).*'.format(comment_symbol),' ',file_string)


def parse_delimited(file_string,delimiter=r'\|',contains=''):
    '''find all strings containing symbols between delimiters'''

    # regex pattern for matching reaction syntax
    pattern = r'(?:[^{}\n]*{}[^{}\n]*)'.format(delimiter,contains,delimiter)
    reactions = []

    for match in finditer(pattern,file_string) :
            reaction = match.group()

            # parse with python compatible syntax
            reactions += [ reaction.replace('^','**') ]
    return reactions
