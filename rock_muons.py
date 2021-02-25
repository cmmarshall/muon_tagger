#!/usr/bin/env python
import ROOT
import random
from optparse import OptionParser

# These are the total POT in the input files
lar_pot = 3.51505E17
rock_pot = 8.89006E15 # this is only ~100 spills, not great

# hard-coded position of the LAr active volume relative to (0,0,0) at the center of the upstream face
offset = [ 0., 5., 410. ]

# where did the muon enter the detector from
entries = ["front", "side", "top", "bottom"]
colors = [ROOT.kGreen+2, ROOT.kRed, ROOT.kBlue, ROOT.kCyan+1]

class Muon:
    def __init__(self, entry=None, exit=None):
        self.entry = ROOT.TVector3(entry.x()/10. - offset[0], entry.y()/10. - offset[1], entry.z()/10. - offset[2])
        self.exit =  ROOT.TVector3(exit.x()/10. - offset[0], exit.y()/10. - offset[1], exit.z()/10. - offset[2])
        # Determine whether the muon entered from the front, sides, top/bottom
        dist_to_edge = [abs(self.entry.z()-0.), min(abs(self.entry.x()+357.), abs(self.entry.x()-357.)), abs(self.entry.y()-150.), abs(self.entry.y()+150.)]
        smallest = min(dist_to_edge)
        self.where = entries[dist_to_edge.index(smallest)]

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
    nonMuonEnergy = 0.
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
        n = {}
        for entry in entries: n[entry] = 0
        for muon in muons:
            h["entry_xy"].Fill( muon.entry.x(), muon.entry.y() )
            if muon.where == "top":
                h["entry_ztop"].Fill( muon.entry.z() )
            elif muon.where == "side":
                h["entry_zside"].Fill( muon.entry.z() )

            n[muon.where] += 1

        for entry in entries:
            h["mult_%s" % entry].Fill( n[entry] )

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
    h["mult"] = ROOT.TH1D( "muon_multiplicity", ";Number of muons;Fraction of 1.2MW spills", 30, 0., 30. )
    for entry in entries:
        h["mult_%s" % entry] = ROOT.TH1D( "mult_%s" % entry, ";Muons per spill", 30, 0., 30. )
    h["entry_xy"] = ROOT.TH2D( "entry_xy", ";Entry x (cm);Entry y (cm)", 72, -360., 360., 32, -160., 160. )
    h["entry_ztop"] = ROOT.TH1D( "entry_ztop", ";Entry z (cm)", 25, 0., 500. )
    h["entry_zside"] = ROOT.TH1D( "entry_zside", ";Entry z (cm)", 25, 0., 500. )

    loop( RockEvents, h, 7.5E13 )

    fout = ROOT.TFile( args.outfile, "RECREATE" )
    for key in h:
        h[key].Write()

    c = ROOT.TCanvas()
    leg = ROOT.TLegend(0.6, 0.5, 0.846, 0.846)
    leg.AddEntry( h["mult"], "All", "l" )
    h["mult"].Scale( 1. / h["mult"].Integral() )
    h["mult"].SetMaximum(0.4)
    h["mult"].Draw("hist")
    for i,entry in enumerate(entries):
        h["mult_%s" % entry].Scale( 1. / h["mult_%s" % entry].Integral() )
        h["mult_%s" % entry].SetLineColor(colors[i])
        leg.AddEntry( h["mult_%s" % entry], entry.capitalize(), "l" )
        h["mult_%s" % entry].Draw("hist same")
    leg.Draw()
    c.Print( "muon_multiplicity.png" )


