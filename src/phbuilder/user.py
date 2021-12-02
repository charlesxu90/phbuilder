# Provides the user with debug, updates, warnings, and errors.
class User:
    def __init__(self, verbosity):
        self.d_verbosity = verbosity

    # Used for messages helpful for debugging.
    def verbose(self, message):
        if self.d_verbosity == 3:
            # Format list or dict.
            if type(message) == type([]) or type(message) == type({}):
                print("phbuilder : ", end='')
                print(message)
            # Format newline at start.
            elif message[0] == '\n':
                print("\nphbuilder : {}".format(message[1:]))
            # Format normal message.
            else:
                print("phbuilder : {}".format(message))

    # Prints regular user update.
    def update(self, message):
        if self.d_verbosity >= 2:
            # Format list or dict.
            if type(message) == type([]) or type(message) == type({}):
                print("phbuilder : ", end='')
                print(message)
            # Format newline at start.
            elif message[0] == '\n':
                print("\nphbuilder : {}".format(message[1:]))
            # Format normal message.
            else:
                print("phbuilder : {}".format(message))

    # Prints warning message for the user.
    def warning(self, message):
        if self.d_verbosity >= 1:
            print("phbuilder : WARNING - {}".format(message))

    # Prints error message for the user and quits program.
    def error(self, message):
        if self.d_verbosity >= 1:
            print("phbuilder : ERROR - {} quiting...".format(message))
        quit()
