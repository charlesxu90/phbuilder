import universe

# Print user update.
def update(message, verbosity=2):

    # If specified verbosity is larger than verbosity set for the specific message, print.
    if (universe.get('d_verbosity') >= verbosity):

        # Formatting for if we want to print a list.
        if type(message) == type([]) or type(message) == type({}):

            print("phbuilder : ", end='')
            print(message)

        # Formatting for if we want a newline at the start.
        elif message[0] == '\n':
            print("\nphbuilder : {}".format(message[1:]))

        else:
            print("phbuilder : {}".format(message))

# Prints error for the users and quits program.
def error(message, verbosity=2):

    if (universe.get('d_verbosity') >= verbosity):

        print("phbuilder : ERROR - {} quiting...".format(message))

    quit()
