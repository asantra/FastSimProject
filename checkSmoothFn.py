#### check the ROOT smooth function
import os, sys
from ROOT import *
import argparse
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *


def main():
    gROOT.SetBatch()

    parser = argparse.ArgumentParser(description='Code to check 2D plots and smoothening function')
    parser.add_argument('-det', action="store", dest="detid", type=str, default="33")
    parser.add_argument('-nSmooth', action="store", dest="numSmooth", type=int, default=1)
    args = parser.parse_args()

    outDir = "TwoDimensionalPlots_Smoothed"+str(args.numSmooth)+"times"

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    inDir = "/Users/arkasantra/arka/Sasha_Work/OutputFile/"
    inFileName = inDir+"/LUXEDumpFiles_FullSim_0p06BX_DetId"+args.detid+"_NoECutNtrn_CoarseBinning.root"
    inFile = TFile(inFileName,"READ")

    rUp_theta_neutron = inFile.Get("dump_plane_bkg_track_rUp_track_theta_backward_neutron_cut")
    rDn_theta_neutron = inFile.Get("dump_plane_bkg_track_rDn_track_theta_backward_neutron_cut")
    rUp_theta_photon  = inFile.Get("dump_plane_bkg_track_rUp_track_theta_backward_photon_cut")
    rDn_theta_photon  = inFile.Get("dump_plane_bkg_track_rDn_track_theta_backward_photon_cut")
    
    clone_rUp_theta_neutron = rUp_theta_neutron.Clone('clone_rUp_theta_neutron')
    clone_rDn_theta_neutron = rDn_theta_neutron.Clone('clone_rDn_theta_neutron')
    clone_rUp_theta_photon  = rUp_theta_photon.Clone('clone_rUp_theta_photon')
    clone_rDn_theta_photon  = rDn_theta_photon.Clone('clone_rDn_theta_photon')
    
    clone_rUp_theta_neutron.Smooth(args.numSmooth) 
    clone_rDn_theta_neutron.Smooth(args.numSmooth) 
    clone_rUp_theta_photon.Smooth(args.numSmooth)  
    clone_rDn_theta_photon.Smooth(args.numSmooth)


    LegendName   = ['Unsmoothed', 'Smoothed ('+str(args.numSmooth)+' times)']
    PlotColor    = [kBlue, kRed]
    yAxisName    = "#theta [rad]"
    xrange1down  = 0
    xrange1up    = 1000
    yrange1down  = 2.8
    yrange1up    = 3.2


    drawline     = False
    logy         = False
    latexName2   = ''
    latexName3   = ''
    leftLegend   = False
    doAtlas      = False
    doLumi       = False
    noRatio      = False
    do80         = False
    do59         = False
    drawPattern  = "COLZ"
    logz         = True

    latexName    = 'neutron'
    xAxisName    = "r_{up} [mm]" 
    FirstTH1     = [rUp_theta_neutron, clone_rUp_theta_neutron]
    ### draw unsmooth and smoothed TH2D together
    DrawHists2Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, outDir+"/rUp_theta_"+latexName, 1.0, 1.0, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz)

    xAxisName    = "r_{dn} [mm]" 
    FirstTH1     = [rDn_theta_neutron, clone_rDn_theta_neutron]
    ### draw unsmooth and smoothed TH2D together
    DrawHists2Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, outDir+"/rDn_theta_"+latexName, 1.0, 1.0, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz)

    latexName    = 'photon'
    xAxisName    = "r_{up} [mm]" 
    FirstTH1     = [rUp_theta_photon, clone_rUp_theta_photon]
    ### draw unsmooth and smoothed TH2D together
    DrawHists2Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, outDir+"/rUp_theta_"+latexName, 1.0, 1.0, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz)

    xAxisName    = "r_{dn} [mm]" 
    FirstTH1     = [rDn_theta_photon, clone_rDn_theta_photon]
    ### draw unsmooth and smoothed TH2D together
    DrawHists2Canvas(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, outDir+"/rDn_theta_"+latexName, 1.0, 1.0, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz)







if __name__=="__main__":
    main()