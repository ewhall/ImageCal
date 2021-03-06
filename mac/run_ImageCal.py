import os,sys

if len(sys.argv) !=4:
	print()
	print("CONFIG_FILE 	= str(sys.argv[1])")
	print("IMG_FILE 	= str(sys.argv[2])")
	print("OUTPUT_DIR	= str(sys.argv[3])")
	print()
	sys.exit(1)

import ROOT
from ROOT import std
from ROOT import larcv

ROOT.gROOT.SetBatch(False) #set to true before production

CONFIG_FILE = str(sys.argv[1])
IMG_FILE	= str(sys.argv[2])
OUTPUT_DIR 	= str(sys.argv[3])

#num = int(os.path.basename(IMG_FILE).split(".")[0].split("_")[-1])

BASE_PATH = os.path.realpath(__file__)
BASE_PATH = os.path.dirname(BASE_PATH)
sys.path.insert(0,BASE_PATH)


print("Process Driver")
proc = larcv.ProcessDriver('ProcessDriver')

print("Load Config File")
proc.configure(ROOT.std.string(CONFIG_FILE))
print("flist")
flist=ROOT.std.vector('std::string')()
flist.push_back(ROOT.std.string(IMG_FILE))
proc.override_input_file(flist)

proc.override_ana_file(ROOT.std.string(os.path.join(OUTPUT_DIR,'tracker_anaout.root')))

alg_id	= proc.process_id('ImageCal')
alg   	= proc.process_ptr(alg_id)
#print('GOT: ',alg,'@ id='alg_id)

alg.SetOutDir(OUTPUT_DIR)
#need to think of an actual name for the output file, but that's not super important right now
alg.SetLLOutName(ROOT.std.string(os.path.join(OUTPUT_DIR,'new_image_file.root')))

proc.initialize()
proc.batch_process(0,1)
proc.finalize()





