import shelve

# UNIVERSE #####################################################################################################################################################

# Set/update variable to universe
def add(varName, value):
    shelve.open('universe')[varName] = value

# Check whether universe contains a certain varName
def has(varName):
    return varName in shelve.open('universe')

# Retrieve variable from universe
def get(varName):
    if has(varName):
        return shelve.open('universe')[varName]

    data = eval(input("couldn't retrieve var \"{0}\" from  Enter manually: ".format(varName)))
    print("add {0} = {1} {2}".format(varName, data, type(data)))
    add(varName, data)
    return data

# Display all variables (name, data, type) stored in the universe
def inspect():
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in sorted(shelf):
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}] ({3}) {4}".format(item, shelf[item][0], shelf[item][-1], len(shelf[item]), type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1} {2}".format(item, shelf[item], type(shelf[item]), arg=longest))

# OUTPUT #######################################################################################################################################################

# Prints debug information for the user.
def pedantic(tool, message):
    if (get('d_verbosity') > 2):
        print("{:18s} : {:s}".format(tool, message))

# Prints update for the user.
def update(tool, message):
    if (get('d_verbosity') > 1):
        print("{:18s} : {:s}".format(tool, message))

# Prints warning for the user.
def warning(tool, message):
    if (get('d_verbosity') > 0):
        print("{:18s} : WARNING - {:s}".format(tool, message))

# Prints error for the users and quits program.
def error(tool, message):
    if (get('d_verbosity') > 0):
        print("{:18s} : ERROR - {:s} quiting...".format(tool, message))
    quit()
