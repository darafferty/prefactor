def plugin_main(args, **kwargs):
    """
    Makes a losoto parset file for given inputs
    """
    steps = kwargs['steps'].lstrip('[').rstrip(']').replace(' ', '').split(',')

    parset_dict = {}
    parset_dict['global'] = {}
    for step in steps:
        parset_dict[step.strip()] = {}

    for key in kwargs.keys():
        keyword = key.split('.')[0]
        option  = key.split(keyword + '.')[-1]
        value   = kwargs[key]
        if keyword in steps:
            parset_dict[keyword][option] = value
        elif keyword == 'global':
            parset_dict[keyword][option] = value

    parset_file = open(kwargs['filename'], 'w')
    for option in parset_dict['global']:
        parset_file.write(option + ' = ' + parset_dict['global'][option] + '\n')

    for step in steps:
        parset_file.write('[' + step + ']\n')
        if step in parset_dict.keys():
            for option in parset_dict[step]:
                parset_file.write(option.replace(';', '.') + ' = ' + parset_dict[step][option] + '\n')
    parset_file.close()

    return 0
