from sys import exit

def ERROR(msg) :
    """
    Output Error message (msg) and stop the program
    """
    print("[Error] " + msg + "\n")
    exit(1)