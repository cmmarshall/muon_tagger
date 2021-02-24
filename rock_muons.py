#!/usr/bin/env python
import ROOT
import random
from optparse import OptionParser

# These are the total POT in the input files
lar_pot = 3.51505E17
rock_pot = 8.89006E15 # this is only ~100 spills, not great

# hard-coded position of the LAr active volume relative to (0,0,0) at the center of the upstream face
offset = [ 0., 5., 410. ]

class Muon:
    def __init__(self, entry=None, exit=None):
        self.entry = entry
        self.exit = exit

    def setEntry(self, v):
        self.entry = ROOT.TVector3(v.x()/10. - offset[0], v.y()/10. - offset[1], v.z()/10. - offset[2])

    def setExit(self, v):
        self.exit = ROOT.TVector3(v.x()/10. - offset[0], v.y()/10. - offset[1], v.z()/10. - offset[2])

# Process a single event
def DoEvent( event ):

    # Find hits in ArgonCube active regions
    hits = []
    for key in event.SegmentDetectors:
       if key.first == "ArgonCube":
            hits += key.second

    # Loop over ArgonCube energy deposits
    muon = None
    for hit in hits: # use the hits to find muons that go into the detector
        hitpdg = event.Trajectories[hit.Contrib[0]].PDGCode # pdg code of the thing that made the hit
        if abs(hitpdg) == 13: # hit due to a muon
            if muon is None: # first muon hit
                muon = Muon( hit.Start.Vect(), hit.Stop.Vect() ) # the end position will be updated
            else:
                muon.setExit( hit.Stop.Vect() )

    # there is one or zero muons, will be none if no muon hits LAr
    return muon

# Main event loop
def loop( RockEvents, h, spill_pot ):

    event = ROOT.TG4Event()
    RockEvents.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = RockEvents.GetEntries() # number of rock events in rock_pot POT
    rock_per = N * spill_pot / rock_pot # mean number of rock events per spill at whatever intensity

    # RNG for poisson fluctuating the events per spill
    rando = ROOT.TRandom3(12345)

    # Build spills until we run out of events
    ispill = 0 # spill index
    i = 0 # event index
    while i < N:

        # Poisson fluctuate number of events in this spill
        this_spill = rando.Poisson(rock_per)
        last_event = i + this_spill # actually last event + 1
        print "Generating spill %d with %d rock events; %d remaining..." % (ispill, this_spill, N-i)

        # if we run out of events just stop, don't do a partial spill
        if last_event >= N: break

        # loop over rock muons
        muons = []
        while i < last_event:
            RockEvents.GetEntry(i)
            muon = DoEvent( event )
            if muon is not None:
                muons.append( muon )
            i += 1

        # Fill histograms
        h["mult"].Fill( len(muons) )
        for muon in muons:
            h["entry_xy"].Fill( muon.entry.x(), muon.entry.y() )
            if muon.entry.z() > 10. and muon.entry.y() > 140.:
                h["entry_ztop"].Fill( muon.entry.z() )
            elif muon.entry.z() > 10. and abs(muon.entry.y()) < 140.:
                h["entry_zside"].Fill( muon.entry.z() )

        ispill += 1

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)

    (args, dummy) = parser.parse_args()

    loaded = False
    RockEvents = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"

    # we have 246 "runs"
    for run in range(10,246):
        rockname = "/pnfs/dune/persistent/users/marshalc/full_hall/rock.%s.%d.edep.root" % (neutrino, run)

        tf = ROOT.TFile( rockname )
        if not tf.TestBit(ROOT.TFile.kRecovered): 
            RockEvents.Add( rockname )

        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        tf.Close()

    # output histograms
    h = {}
    h["mult"] = ROOT.TH1D( "muon_multiplicity", ";Number of muons", 50, 0., 50. )
    h["entry_xy"] = ROOT.TH2D( "entry_xy", ";Entry x (cm);Entry y (cm)", 180, -360., 360., 80, -160., 160. )
    h["entry_ztop"] = ROOT.TH1D( "entry_ztop", ";Entry z (cm)", 265, -10., 520. )
    h["entry_zside"] = ROOT.TH1D( "entry_zside", ";Entry z (cm)", 265, -10., 520. )

    loop( RockEvents, h, 7.5E13 )

    fout = ROOT.TFile( args.outfile, "RECREATE" )
    for key in h:
        h[key].Write()


