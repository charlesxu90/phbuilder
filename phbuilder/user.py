import sys

class User:
    """Provides and formats user updates, warnings, and errors.
    """

    def __init__(self, verbosity: int, logFileName='', maxLineLength=90):
        """Initialize User object.

        Args:
            verbosity (int): verbosity value. 0 = supress all, 1 = only warnings and errors, 2 = regular, 3 = verbose.
            maxLineLength (int): maximum length of sentence before a newline (does not consider length of preMessage etc).
        """

        self.verbosity = verbosity
        self.logFileName = logFileName
        self.maxLineLength = maxLineLength

    def doPrint(self, message: str):
        """Handle printing internals.

        Args:
            message (str): message to be printen/logged.
        """

        if self.logFileName:
            with open(self.logFileName, 'a+') as logfile:
                print(message, file=logfile)

        print(message)

    def base(self, message: str, preMessage: str = ''):
        """Base method for handeling user messages.

        Args:
            message (str): message.
            preMessage (str, optional): pre-message. Defaults to ''.
        """

        # If we input a dictionary of list as a message, print and return.
        if isinstance(message, list) or isinstance(message, dict):
            self.doPrint(message)
            return

        firstLineWritten = False
        charCount = 0
        currentLine = ''

        for word in message.split(' '):

            charCount += len(word) + 1

            if charCount < self.maxLineLength:
                currentLine += word + ' '
            else:
                self.doPrint(f"phbuilder : {preMessage}{currentLine}")
                firstLineWritten = True
                charCount = len(word) + 1
                currentLine = word + ' '

            # This we add so that not every subsequent line has preMessage,
            # but it is aligned nonetheless.
            if firstLineWritten and preMessage != '':
                preMessage = ' '.ljust(len(preMessage))

        self.doPrint(f"phbuilder : {preMessage}{currentLine.lstrip()}")  # Flush.

    def verbose(self, message: str):
        """Print message when verbosity is high.

        Args:
            message (str): message.
        """

        if self.verbosity > 2:
            self.base(message)

    def update(self, message: str):
        """Print default message.

        Args:
            message (str): message.
        """

        if self.verbosity > 1:
            self.base(message)

    def warning(self, message: str):
        """Print warning message.

        Args:
            message (str): message.
        """

        if self.verbosity > 0:
            self.doPrint('phbuilder : ')
            self.base(message, preMessage='WARNING - ')
            self.doPrint('phbuilder : ')

    def error(self, message: str):
        """Print error message and sys.exit() the program.

        Args:
            message (str): message.
        """

        if self.verbosity > 0:
            self.doPrint('phbuilder : ')
            self.base(message, preMessage='ERROR - ')
            self.doPrint('phbuilder : ')

        sys.exit()

    def inputOptionHandler(self, message: str, options: list) -> int:
        """Handles user input.

        Args:
            message (str): message.
            options (list): list of strings describing the options. Starts counting from 0.

        Returns:
            int: number corresponding to the option selected by the user.
        """

        valids = []
        msgstring = "phbuilder : {}:".format(message)

        # Loop through the options list and create string for display
        for idx in range(0, len(options)):
            msgstring += "\nphbuilder : {}. {}".format(idx, options[idx])
            valids.append(str(idx))

        while True:
            self.doPrint(msgstring)
            val = input("phbuilder : Type a number: ")

            if val in valids:
                self.doPrint(f"phbuilder : selected {val}")
                self.doPrint('')
                return int(val)

            self.doPrint("phbuilder : {} is not a valid option, please try again:\n".format(val))
