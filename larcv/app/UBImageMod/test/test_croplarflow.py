import ROOT as rt
from larcv import larcv
import sys
sys.argv.append("-b")


superafile = "/media/hdd2/taritree/larflow/xfer/larcv_5477923_0.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "baka.root" )
io.initialize()

# -------------------------------------
# UBSplitDetector

split_cfg="""Verbosity:0
InputProducer: \"wire\"
OutputBBox2DProducer: \"detsplit\"
CropInModule: true
OutputCroppedProducer: \"detsplit\"
BBoxPixelHeight: 512
BBoxPixelWidth: 832
CoveredZWidth: 310
FillCroppedYImageCompletely: true
DebugImage: false
MaxImages: 20
RandomizeCrops: true
MaxRandomAttempts: 50
MinFracPixelsInCrop: 0.0001
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,split_cfg
fcfg.close()
split_pset = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )

# -------------------------------------
# UBLArFlowCropDetector

lfcrop_cfg="""Verbosity:0
InputBBoxProducer: \"detsplit\"
InputADCProducer: \"wire\"
InputCroppedADCProducer: \"detsplit\"
InputVisiProducer: \"pixvisi\"
InputFlowProducer: \"pixflow\"
OutputCroppedADCProducer: \"adc\"
OutputCroppedVisiProducer: \"visi\"
OutputCroppedFlowProducer: \"flow\"
OutputFilename: \"baka_lf.root\"
CheckFlow: false
MakeCheckImage: false
"""

lfcfg = open("ublarflowcrop.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "ublarflowcrop.cfg", "UBCropLArFlow" )

# -------------------------------------
# ALGOS

split_algo = larcv.UBSplitDetector()
split_algo.configure(split_pset)
split_algo.initialize()

lfcrop_algo = larcv.UBCropLArFlow()
lfcrop_algo.configure(lfpset)
lfcrop_algo.initialize()

# -------------------------------------

nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 1

for n in range(nentries):
    io.read_entry(n)
    split_algo.process( io )
    lfcrop_algo.process( io );

lfcrop_algo.finalize()
io.finalize()


print "FIN"
