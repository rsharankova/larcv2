import ROOT
from larcv import larcv


superafile = "/media/hdd2/larflow/xfer/larcv_5477923_0.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "baka.root" )
io.initialize()

# -------------------------------------
# UBSplitDetector

scfg="""InputProducer: \"wire\"
OutputBBox2DProducer: \"detsplit\"
BBoxPixelHeight: 512
BBoxPixelWidth: 832
CoveredZWidth: 310
OutputCroppedProducer: \"detsplit\"
CropInModule: false
DebugImage: true
MaxImages: 2
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,scfg
fcfg.close()


cfg = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )
algo = larcv.UBSplitDetector()
algo.configure(cfg)
algo.initialize()

# -------------------------------------

cfgcrop="""
"""


nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 1

for n in range(nentries):
    io.read_entry(n)
    algo.process( io )
    io.save_entry()

io.finalize()

print "FIN"
