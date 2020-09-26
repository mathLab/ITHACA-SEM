from NekPy.LibUtilities import SessionReader
import sys

session = SessionReader.CreateInstance(sys.argv)

print("Loaded session: %s" % session.GetSessionName())
