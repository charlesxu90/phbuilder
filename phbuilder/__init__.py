# PYTHON_ARGCOMPLETE_OK

def entryFunction():
    # PYTHON_ARGCOMPLETE_OK

    # Execute this as soon as possible because of the TAB autocomplete thing.
    from . import parsecmdline

    # Parse command line input.
    CLI = parsecmdline.parsecmdline()

    # Import everything else.
    from . import phbuilder

    # Run main program.
    phbuilder.phbuilder(CLI).runner()
