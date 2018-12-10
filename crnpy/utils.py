

def isnumber(value):
    '''retrun boolean depending on whether the
    current value is parsable as a number.'''

    try:
        float(value)
        return True

    except (ValueError,TypeError):
        return False
