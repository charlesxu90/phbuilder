import universe

# Prints debug information for the user.
def pedantic(tool, message):
    if (universe.get('d_verbosity') > 2):
        print("{:18s} : {:s}".format(tool, message))

# Prints update for the user.
def update(tool, message):
    if (universe.get('d_verbosity') > 1):
        print("{:18s} : {:s}".format(tool, message))

# Prints warning for the user.
def warning(tool, message):
    if (universe.get('d_verbosity') > 0):
        print("{:18s} : WARNING - {:s}".format(tool, message))

# Prints error for the users and quits program.
def error(tool, message):
    if (universe.get('d_verbosity') > 0):
        print("{:18s} : ERROR - {:s} quiting...".format(tool, message))
    quit()
