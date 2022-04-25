import sys
from System.MSStaff import CStaff
import multiprocessing

if __name__=="__main__":
    multiprocessing.freeze_support()
    staff=CStaff(sys.argv)
    staff.start()