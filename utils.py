import universe

# Print user update.
def update(message, verbosity=2):
    
    if (universe.get('d_verbosity') >= verbosity):
    
        if (type(message) == type([])):
            
            print("phbuilder : ", end='')
            print(message)
        
        else:
            print("phbuilder : {}".format(message))

# Prints error for the users and quits program.
def error(message, verbosity=2):
    
    if (universe.get('d_verbosity') >= verbosity):
    
        print("phbuilder : ERROR - {} quiting...".format(message))
    
    quit()
