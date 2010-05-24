# my_logger.py
"""This is 'my_logger' module, which is imported into all
the other modules of my application."""
import logging
import logging.handlers
from functools import wraps
# Create a global logger
_vs_logger = None
# Set default DEBUG_ON True - in which case debug messages
# are saved to the log file. Or set it to False - in which
# case only INFO and ERROR messages are saved to the log file
_DEBUG_ON = True
def set_logger():
    "Set up the logger"
    global _vs_logger
    _vs_logger = logging.getLogger("my_logger")
    # Set the logger level
    if _DEBUG_ON:
        _vs_logger.setLevel(logging.DEBUG)
    else:
        _vs_logger.setLevel(logging.INFO)
    # Set the format
    form = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        # Add the log message handler to the logger
        # The location of the logfile is given here as 'logfile.txt'
        # in an actual application I would take a bit of care
        # where this is located
    handler = logging.handlers.RotatingFileHandler("logfile.txt",
                                                       maxBytes=20000,
                                                       backupCount=5)
    handler.setFormatter(form)
    _vs_logger.addHandler(handler)
        
def info_log(message):
    "Log message with level info"
    if _vs_logger:
        _vs_logger.info(str(message))
def debug_log(message):
    "Log message with level debug"
    if _DEBUG_ON and _vs_logger:
        _vs_logger.debug(str(message))
def logmethod(f):
    "Creates a decorator to log a method"
    @wraps(f)
    def wrapper(self, *args, **kwds):
        debug_log("%s in %s called" % (f.__name__, self.__class__.__name__))
        return f(self, *args, **kwds)
    return wrapper
def exception_log(message):
    "Log message with level error plus exception traceback"
    if _vs_logger:
        _vs_logger.exception(str(message))


from my_logger import *
# set up the logger
set_logger()
# now, in any part of my program I can record
# a 'debug' message - if the global _DEBUG_ON
# is True - otherwise nothing will be recorded
debug_log("my program has reached this point")
# I can also record an 'info' message - which is
# recorded whether _DEBUG_ON is True or False
# presumably because more important information
# is needed to be saved
info_log("This is important")
# I also like to record exceptions
try:
    # a block of code which may cause an error
    x = 1/0
    # That usually does it!
except Exception, e:
    # handle the error
    # and record it
    exception_log(e)

info_log("so - we are now past the exception")
# as a final example, I sometimes like
# to record if a method
# in a class has been called

class MyClass:
    @logmethod
    def mymethod(self):
        pass

x = MyClass()
# and calling its method should
# create a log if _DEBUG_ON is True
x.mymethod()
