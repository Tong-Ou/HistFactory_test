import glob, os, sys
import commands
from ROOT import *
from glob import *
import time

def help():
    print " "
    print "python runInspector.py <WS file>  <WS name>  <data name> "
    print "* <WS file> is MANDATORY"
    print "* <WS name> (combined)"
    print "* <data name> (obsData)"
    print " "

if len(sys.argv)==1:
    print " "
    print " NO filename specified " 
    help()
    sys.exit()

WSfile  =""
WSname  ="combined"
dataName="obsData"

WSfile=sys.argv[1]
if not os.path.isfile(WSfile):
    print " "
    print " file: '"+WSfile+"' does NOT exists ... please check"
    print " "
    sys.exit()

if len(sys.argv)>2:
    WSname=sys.argv[2]
Rfile=TFile(WSfile)
theWS=None
theWS=Rfile.Get(WSname)
if theWS==None:
    print " "
    print " WS '"+WSname+"' could not be found inside file '"+WSfile+"' ... please check"
    print " "
    sys.exit()
Rfile.Close()

#### and now the real command
gROOT.ProcessLine(".L WSinspector.C++g")
command="LimitCrossCheck::PlotFitCrossChecks(\""+WSfile+"\",\"./test\",\""+WSname+"\" ,\"ModelConfig\" ,\""+dataName+"\")"
print " "
print "=========================================================================================================================================================================================="
print command
gROOT.ProcessLine(command)
#gROOT->ProcessLine("LimitCrossCheck::PlotFitCrossChecks(\"WSs/125.root\" ,\"./test\"   ,\"combined\" ,\"ModelConfig\" ,\"obsData\")" );
