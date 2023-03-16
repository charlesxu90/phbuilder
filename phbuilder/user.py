import sys

class User:
    """Provides and formats user updates, warnings, and errors.
    """

    def __init__(self, verbosity: int, maxLineLength=90):
        """Initialize User object.

        Args:
            verbosity (int): verbosity value. 0 = supress all, 1 = only warnings and errors, 2 = regular, 3 = verbose.
            maxLineLength (int): maximum length of sentence before a newline (does not consider length of preMessage etc).
        """

        self.verbosity = verbosity
        self.maxLineLength = maxLineLength

    def base(self, message: str, preMessage: str = ''):
        """Base method for handeling user messages.

        Args:
            message (str): message.
            preMessage (str, optional): pre-message. Defaults to ''.
        """

        # If we input a dictionary of list as a message, print and return.
        if isinstance(message, list) or isinstance(message, dict):
            print(message)
            return

        firstLineWritten = False
        charCount = 0
        currentLine = ''

        for word in message.split(' '):

            charCount += len(word) + 1

            if charCount < self.maxLineLength:
                currentLine += word + ' '
            else:
                print(f"phbuilder : {preMessage}{currentLine}")
                firstLineWritten = True
                charCount = len(word) + 1
                currentLine = word + ' '

            # This we add so that not every subsequent line has preMessage,
            # but it is aligned nonetheless.
            if firstLineWritten and preMessage != '':
                preMessage = ' '.ljust(len(preMessage))                

        print(f"phbuilder : {preMessage}{currentLine.lstrip()}")  # Flush.

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
            print('phbuilder : ')
            self.base(message, preMessage='WARNING - ')
            print('phbuilder : ')

    def error(self, message: str):
        """Print error message and sys.exit() the program.

        Args:
            message (str): message.
        """

        if self.verbosity > 0:
            print('phbuilder : ')
            self.base(message, preMessage='ERROR - ')
            print('phbuilder : ')

        sys.exit()
