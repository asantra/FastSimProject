#### run: python3 makeDumpPlotsFromText.py -l  <event file in text mode>  -b <bx number> -d <detId>
#### Here the histograms are made from FullSim, coming from Sasha
import os
import sys
import time
import pprint
import math
from array import array
from ROOT import *
from collections import OrderedDict
import argparse
    
def getTrackThetaPhi(x, y, z):
    phi   = math.atan2(y, x)
    dist  = math.sqrt(x**2+y**2+z**2) 
    theta = math.acos(z/dist)
    return theta, phi
    
    
def energyBins(nbins):
  #//// variable binned in X axis histograms
  logmin           = -7;
  logmax           = 1;
  logbinwidth      = (logmax-logmin)/float(nbins);
  xpoints          = []
  
  for i in range(0, nbins+1):
      #print((logmin + i*logbinwidth), pow( 10,(logmin + i*logbinwidth) ))
      xpoints.append(pow( 10,(logmin + i*logbinwidth) ))
                     
  xpoints.append(2*pow( 10,1))
  return xpoints



def timeBins(nbins):
  #//// variable binned in X axis histograms
  logmin           = -1;
  logmax           = 9;
  logbinwidth      = (logmax-logmin)/float(nbins);
  tpoints          = []
  
  for i in range(0, nbins+1):
      #print((logmin + i*logbinwidth), pow( 10,(logmin + i*logbinwidth) ))
      tpoints.append(pow( 10,(logmin + i*logbinwidth) ))
                     
  tpoints.append(2*pow( 10,9))
  return tpoints

### send 6300 bins 
def rBins(nbins):
    rPoints = []
    for i in range(0, nbins+1):
        if i < 6000:
            rPoints.append(i*1)
        elif i < 6200:
            rPoints.append(6000+(i-6000)*10)
        else:
            rPoints.append(8000+(i-6200)*20)
    return rPoints




def getR(x,y):
    return math.sqrt(x**2+y**2)


def main():
    
    parser = argparse.ArgumentParser(description='Code to get 2D plots')
    parser.add_argument('-l', action="store", dest="inFile", type=str, default="RestrictedDumpOnlyFiles_DetId33_trackInfo.txt")
    parser.add_argument('-b', action="store", dest="bx", type=float, default=1.0)
    parser.add_argument('-d', action="store", dest="det", type=int, default=33)
    args = parser.parse_args()
    
    inDir = '/Users/arkasantra/arka/Sasha_Work/OutputFile'
    #inDir = "/Volumes/Study/Weizmann_PostDoc/Sasha_work/OutputFile/ReprocessedBkgTracksAfterTDR"
    
    bkgFileName   = open(inDir+'/'+args.inFile)
    withoutText   = args.inFile.split('.txt')[0]
    #rootFile      = withoutText+"_DetId"+str(args.det)+"_timelt1000.root"
    rootFile      = withoutText+"_DetId"+str(args.det)+".root"
    # rootFile      = withoutText+".root"
    #rootFile      = withoutText+"_timelt1000.root"
    #rootFile      = withoutText+"_MainFullSimFile.root"
    nbx           = args.bx
    
    if args.det==33:
        zPos = -350
    elif args.det==32:
        zPos = -5000
    elif args.det==31:
        zPos = -10000
    elif args.det==30:
        zPos = -15000
    else:
        zPos = -350
        
    print("zPos: ", zPos)
    
    
    print('BX selected: ',nbx)
    
    outFile       = TFile(inDir+'/'+rootFile, "RECREATE")
    outFile.cd()
    #treeIn = TTree("treeIn", "treeIn")
    
    #bxNumberBranch      = array('i', [0])
    #pdgIdBranch         = array('i', [0])
    #trackIdBranch       = array('i', [0])
    ##staveIdBranch       = array('i', [0])
    #xPosBranch          = array('d', [0.0])
    #yPosBranch          = array('d', [0.0])
    #energyValBranch     = array('d', [0.0])
    #weightBranch        = array('d', [0.0])
    #vtx_xBranch         = array('d', [0.0])
    #vtx_yBranch         = array('d', [0.0])
    #vtx_zBranch         = array('d', [0.0])
    #parent_idBranch     = array('i', [0])
    #pxxBranch           = array('d', [0.0])
    #pyyBranch           = array('d', [0.0])
    #pzzBranch           = array('d', [0.0])
    #theta_{p}Branch         = array('d', [0.0])
    #phi_{p}Branch           = array('d', [0.0])
    #physprocessBranch   = array('i', [0])
    #timeBranch          = array('d', [0.0])
    #rValueBranch        = array('d', [0.0])
    #rValueWeightBranch  = array('d', [0.0])
    
    
    #treeIn.Branch('bxNumberBranch', bxNumberBranch, 'bxNumberBranch/I')
    #treeIn.Branch('pdgIdBranch', pdgIdBranch, 'pdgIdBranch/I')
    #treeIn.Branch('trackIdBranch', trackIdBranch, 'trackIdBranch/I')
    ##treeIn.Branch('staveIdBranch', staveIdBranch, 'staveIdBranch/I')
    #treeIn.Branch('xPosBranch', xPosBranch, 'xPosBranch/D')
    #treeIn.Branch('yPosBranch', yPosBranch, 'yPosBranch/D')
    #treeIn.Branch('energyValBranch', energyValBranch, 'energyValBranch/D')
    #treeIn.Branch('weightBranch', weightBranch, 'weightBranch/D')
    #treeIn.Branch('vtx_xBranch', vtx_xBranch, 'vtx_xBranch/D')
    #treeIn.Branch('vtx_yBranch', vtx_yBranch, 'vtx_yBranch/D')
    #treeIn.Branch('vtx_zBranch', vtx_zBranch, 'vtx_zBranch/D')
    #treeIn.Branch('parent_idBranch', parent_idBranch, 'parent_idBranch/I')
    #treeIn.Branch('pxxBranch', pxxBranch, 'pxxBranch/D')
    #treeIn.Branch('pyyBranch', pyyBranch, 'pyyBranch/D')
    #treeIn.Branch('pzzBranch', pzzBranch, 'pzzBranch/D')
    #treeIn.Branch('thetaBranch', thetaBranch, 'thetaBranch/D')
    #treeIn.Branch('phiBranch', phiBranch, 'phiBranch/D')
    #treeIn.Branch('physprocessBranch', physprocessBranch, 'physprocessBranch/I')
    #treeIn.Branch('timeBranch', timeBranch, 'timeBranch/D')
    #treeIn.Branch('rValueBranch', rValueBranch, 'rValueBranch/D')
    #treeIn.Branch('rValueWeightBranch', rValueWeightBranch, 'rValueWeightBranch/D')


    nbins  = 450
    xbins  = energyBins(nbins)
    xarray = array('d',xbins)
    
    tbins  = timeBins(nbins)
    tarray = array('d',tbins)
    

    nRBins = 6300
    rbins = rBins(nRBins)
    rarray = array('d',rbins)
        
    # print(rarray)
    allHistoDict  = {}
    
    
    #### 2D plot
    #allHistoDict.update({"dump_plane_bkg_track_x_track_y_photon":TH2D("dump_plane_bkg_track_x_track_y_photon","dump_plane_bkg_track_x_track_y_photon; track x [mm]; track y [mm]",500,  -10000.0, 10000.0, 500, -10000.0, 10000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_x_track_y_neutron":TH2D("dump_plane_bkg_track_x_track_y_neutron","dump_plane_bkg_track_x_track_y_neutron; track x [mm]; track y [mm]",500,  -10000.0, 10000.0, 500, -10000.0, 10000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_phi_photon":TH2D("dump_plane_bkg_track_theta_track_phi_photon","dump_plane_bkg_track_theta_track_phi_photon; #theta_{p} [rad]; #phi_{p} [rad]", 6400, 1.6, 3.2, 640, -3.2, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_phi_neutron":TH2D("dump_plane_bkg_track_theta_track_phi_neutron","dump_plane_bkg_track_theta_track_phi_neutron; #theta_{p} [rad]; #phi_{p} [rad]",6400, 1.6, 3.2, 640, -3.2, 3.2)})
    
    
    
    #### 1D plot
    #allHistoDict.update({"dump_plane_bkg_track_time_photon":TH1D("dump_plane_bkg_track_time_photon","dump_plane_bkg_track_time_photon; t [ns]; Events",nbins, tarray)})
    #allHistoDict.update({"dump_plane_bkg_track_time_neutron":TH1D("dump_plane_bkg_track_time_neutron","dump_plane_bkg_track_time_neutron; t [ns]; Events",nbins, tarray)})
    #allHistoDict.update({"dump_plane_bkg_track_energy_photon":TH1D("dump_plane_bkg_track_energy_photon","dump_plane_bkg_track_energy_photon; E [GeV]; Events",nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_energy_neutron":TH1D("dump_plane_bkg_track_energy_neutron","dump_plane_bkg_track_energy_neutron; E [GeV]; Events",nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_r_photon":TH1D("dump_plane_bkg_track_r_photon","dump_plane_bkg_track_r_photon; r [mm]; Events",nRBins, rarray)})
    #allHistoDict.update({"dump_plane_bkg_track_r_neutron":TH1D("dump_plane_bkg_track_r_neutron","dump_plane_bkg_track_r_neutron; r [mm]; Events",nRBins, rarray)})
    
    #allHistoDict.update({"dump_plane_bkg_track_theta_neutron":TH1D("dump_plane_bkg_track_theta_neutron","dump_plane_bkg_track_theta_neutron; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_phi_neutron":TH1D("dump_plane_bkg_track_phi_neutron","dump_plane_bkg_track_phi_neutron; #phi_{p} [rad]; Events",640, -3.2, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_photon":TH1D("dump_plane_bkg_track_theta_photon","dump_plane_bkg_track_theta_photon; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_phi_photon":TH1D("dump_plane_bkg_track_phi_photon","dump_plane_bkg_track_phi_photon; #phi_{p} [rad]; Events",640, -3.2, 3.2)})
    
    
    #### histograms after applying some cuts
    
    ### 2D plot
    # allHistoDict.update({"dump_plane_bkg_track_x_track_y_neutron_cut":TH2D("dump_plane_bkg_track_x_track_y_neutron_cut","dump_plane_bkg_track_x_track_y_neutron_cut; track x [mm]; track y [mm]",500,  -10000.0, 10000.0, 500,  -10000.0, 10000.0)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_neutron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_neutron_weighted_cut","dump_plane_bkg_track_r_track_theta_neutron_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_neutron_cut":TH2D("dump_plane_bkg_track_r_track_theta_neutron_cut","dump_plane_bkg_track_r_track_theta_neutron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    # allHistoDict.update({"dump_plane_bkg_track_r_track_vtxz_neutron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxz_neutron_weighted_cut","dump_plane_bkg_track_r_track_vtxz_neutron_weighted_cut; r [mm]; vtx_{z} [mm]", nRBins, rarray, 1000, -2000.0,  2000.0)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxz_neutron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxz_neutron_weighted_cut","dump_plane_bkg_track_theta_track_vtxz_neutron_weighted_cut; #theta_{p} [rad]; vtx_{z} [mm]", 1600, 1.6, 3.2, 1000, -2000.0,  2000.0)})
    # allHistoDict.update({"dump_plane_bkg_track_r_track_vtxx_neutron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxx_neutron_weighted_cut","dump_plane_bkg_track_r_track_vtxx_neutron_weighted_cut; r [mm]; vtx_{x} [mm]", nRBins, rarray, 260,  -260.0,  260.0)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxx_neutron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxx_neutron_weighted_cut","dump_plane_bkg_track_theta_track_vtxx_neutron_weighted_cut; #theta_{p} [rad]; vtx_{x} [mm]", 1600, 1.6, 3.2, 260,  -260.0,  260.0)})
    # allHistoDict.update({"dump_plane_bkg_track_vtxz_track_vtxx_neutron_weighted_cut":TH2D("dump_plane_bkg_track_vtxz_track_vtxx_neutron_weighted_cut","dump_plane_bkg_track_vtxz_track_vtxx_neutron_weighted_cut; vtx_{z} [mm]; vtx_{x} [mm]", 1000, -2000.0, 2000.0, 260,  -260.0,  260.0)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_neutron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_E_neutron_weighted_cut","dump_plane_bkg_track_r_track_E_neutron_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_E_neutron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_E_neutron_weighted_cut","dump_plane_bkg_track_theta_track_E_neutron_weighted_cut; #theta_{p} [rad]; E [GeV]", 1600, 1.6, 3.2, nbins, xarray)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_pos_theta_neutron_cut":TH2D("dump_plane_bkg_track_theta_pos_theta_neutron_cut","dump_plane_bkg_track_theta_pos_theta_neutron_cut; #theta_{p} [rad]; #theta_{pos} [rad]", 320, 0.0, 3.2, 320, 0.0, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_phi_pos_phi_neutron_cut":TH2D("dump_plane_bkg_track_phi_pos_phi_neutron_cut","dump_plane_bkg_track_phi_pos_phi_neutron_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_neutron_cut":TH2D("dump_plane_bkg_track_r_track_E_neutron_cut","dump_plane_bkg_track_r_track_E_neutron_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_E_neutron_cut":TH2D("dump_plane_bkg_time_track_E_neutron_cut","dump_plane_bkg_time_track_E_neutron_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_neutron_cut":TH2D("dump_plane_bkg_time_track_r_neutron_cut","dump_plane_bkg_time_track_r_neutron_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_neutron_cut":TH2D("dump_plane_bkg_time_track_theta_neutron_cut","dump_plane_bkg_time_track_theta_neutron_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_neutron_weighted_cut":TH2D("dump_plane_bkg_time_track_r_neutron_weighted_cut","dump_plane_bkg_time_track_r_neutron_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_neutron_weighted_cut":TH2D("dump_plane_bkg_time_track_theta_neutron_weighted_cut","dump_plane_bkg_time_track_theta_neutron_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    
    
    
    
    # allHistoDict.update({"dump_plane_bkg_track_x_track_y_photon_cut":TH2D("dump_plane_bkg_track_x_track_y_photon_cut","dump_plane_bkg_track_x_track_y_photon_cut; track x [mm]; track y [mm]",500,  -10000.0, 10000.0, 500,  -10000.0, 10000.0)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_phi_photon_cut":TH2D("dump_plane_bkg_track_theta_track_phi_photon_cut","dump_plane_bkg_track_theta_track_phi_photon_cut; #theta_{p} [rad]; #phi_{p} [rad]", 1600, 1.6, 3.2, 640, -3.2, 3.2)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_phi_neutron_cut":TH2D("dump_plane_bkg_track_theta_track_phi_neutron_cut","dump_plane_bkg_track_theta_track_phi_neutron_cut; #theta_{p} [rad]; #phi_{p} [rad]",1600, 1.6, 3.2, 640, -3.2, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_weighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_weighted_cut","dump_plane_bkg_track_r_track_theta_photon_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_rweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_rweighted_cut","dump_plane_bkg_track_r_track_theta_photon_rweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut","dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_cut","dump_plane_bkg_track_r_track_theta_photon_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    # allHistoDict.update({"dump_plane_bkg_track_r_track_vtxz_photon_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxz_photon_weighted_cut","dump_plane_bkg_track_r_track_vtxz_photon_weighted_cut; r [mm]; vtx_{z} [mm]", nRBins, rarray, 1000, -2000.0,  2000.0)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxz_photon_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxz_photon_weighted_cut","dump_plane_bkg_track_theta_track_vtxz_photon_weighted_cut; #theta_{p} [rad]; vtx_{z} [mm]", 1600, 1.6, 3.2, 1000, -2000.0,  2000.0)})
    # allHistoDict.update({"dump_plane_bkg_track_r_track_vtxx_photon_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxx_photon_weighted_cut","dump_plane_bkg_track_r_track_vtxx_photon_weighted_cut; r [mm]; vtx_{x} [mm]", nRBins, rarray, 260,  -260.0,  260.0)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxx_photon_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxx_photon_weighted_cut","dump_plane_bkg_track_theta_track_vtxx_photon_weighted_cut; #theta_{p} [rad]; vtx_{x} [mm]", 1600, 1.6, 3.2, 260,  -260.0,  260.0)})
    # allHistoDict.update({"dump_plane_bkg_track_vtxz_track_vtxx_photon_weighted_cut":TH2D("dump_plane_bkg_track_vtxz_track_vtxx_photon_weighted_cut","dump_plane_bkg_track_vtxz_track_vtxx_photon_weighted_cut; vtx_{z} [mm]; vtx_{x} [mm]", 1000, -2000.0, 2000.0, 260,  -260.0,  260.0)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_photon_weighted_cut":TH2D("dump_plane_bkg_track_r_track_E_photon_weighted_cut","dump_plane_bkg_track_r_track_E_photon_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_track_E_photon_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_E_photon_weighted_cut","dump_plane_bkg_track_theta_track_E_photon_weighted_cut; #theta_{p} [rad]; E [GeV]", 1600, 1.6, 3.2, nbins, xarray)})
    # allHistoDict.update({"dump_plane_bkg_track_theta_pos_theta_photon_cut":TH2D("dump_plane_bkg_track_theta_pos_theta_photon_cut","dump_plane_bkg_track_theta_pos_theta_photon_cut; #theta_{p} [rad]; #theta_{pos} [rad]", 320, 0.0, 3.2, 320, 0.0, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_phi_pos_phi_photon_cut":TH2D("dump_plane_bkg_track_phi_pos_phi_photon_cut","dump_plane_bkg_track_phi_pos_phi_photon_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_photon_cut":TH2D("dump_plane_bkg_track_r_track_E_photon_cut","dump_plane_bkg_track_r_track_E_photon_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_E_photon_cut":TH2D("dump_plane_bkg_time_track_E_photon_cut","dump_plane_bkg_time_track_E_photon_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_photon_cut":TH2D("dump_plane_bkg_time_track_r_photon_cut","dump_plane_bkg_time_track_r_photon_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_photon_cut":TH2D("dump_plane_bkg_time_track_theta_photon_cut","dump_plane_bkg_time_track_theta_photon_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_photon_weighted_cut":TH2D("dump_plane_bkg_time_track_r_photon_weighted_cut","dump_plane_bkg_time_track_r_photon_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_photon_weighted_cut":TH2D("dump_plane_bkg_time_track_theta_photon_weighted_cut","dump_plane_bkg_time_track_theta_photon_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    
    
    
    
    #allHistoDict.update({"dump_plane_bkg_track_x_track_y_electron_cut":TH2D("dump_plane_bkg_track_x_track_y_electron_cut","dump_plane_bkg_track_x_track_y_electron_cut; track x [mm]; track y [mm]",500,  -10000.0, 10000.0, 500,  -10000.0, 10000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_phi_electron_cut":TH2D("dump_plane_bkg_track_theta_track_phi_electron_cut","dump_plane_bkg_track_theta_track_phi_electron_cut; #theta_{p} [rad]; #phi_{p} [rad]", 6400, 1.6, 3.2, 640, -3.2, 3.2)})    
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_electron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_electron_weighted_cut","dump_plane_bkg_track_r_track_theta_electron_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_electron_rweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_electron_rweighted_cut","dump_plane_bkg_track_r_track_theta_electron_rweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_electron_thetaweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_electron_thetaweighted_cut","dump_plane_bkg_track_r_track_theta_electron_thetaweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_electron_cut":TH2D("dump_plane_bkg_track_r_track_theta_electron_cut","dump_plane_bkg_track_r_track_theta_electron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_vtxz_electron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxz_electron_weighted_cut","dump_plane_bkg_track_r_track_vtxz_electron_weighted_cut; r [mm]; vtx_{z} [mm]", nRBins, rarray, 1000, -2000.0,  2000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxz_electron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxz_electron_weighted_cut","dump_plane_bkg_track_theta_track_vtxz_electron_weighted_cut; #theta_{p} [rad]; vtx_{z} [mm]", 6400, 1.6, 3.2, 1000, -2000.0,  2000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_vtxx_electron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxx_electron_weighted_cut","dump_plane_bkg_track_r_track_vtxx_electron_weighted_cut; r [mm]; vtx_{x} [mm]", nRBins, rarray, 260,  -260.0,  260.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxx_electron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxx_electron_weighted_cut","dump_plane_bkg_track_theta_track_vtxx_electron_weighted_cut; #theta_{p} [rad]; vtx_{x} [mm]", 6400, 1.6, 3.2, 260,  -260.0,  260.0)})
    #allHistoDict.update({"dump_plane_bkg_track_vtxz_track_vtxx_electron_weighted_cut":TH2D("dump_plane_bkg_track_vtxz_track_vtxx_electron_weighted_cut","dump_plane_bkg_track_vtxz_track_vtxx_electron_weighted_cut; vtx_{z} [mm]; vtx_{x} [mm]", 1000, -2000.0, 2000.0, 260,  -260.0,  260.0)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_E_electron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_E_electron_weighted_cut","dump_plane_bkg_track_r_track_E_electron_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_E_electron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_E_electron_weighted_cut","dump_plane_bkg_track_theta_track_E_electron_weighted_cut; #theta_{p} [rad]; E [GeV]", 6400, 1.6, 3.2, nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_E_electron_cut":TH2D("dump_plane_bkg_track_r_track_E_electron_cut","dump_plane_bkg_track_r_track_E_electron_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    
    
    
    
    #allHistoDict.update({"dump_plane_bkg_track_x_track_y_positron_cut":TH2D("dump_plane_bkg_track_x_track_y_positron_cut","dump_plane_bkg_track_x_track_y_positron_cut; track x [mm]; track y [mm]",500,  -10000.0, 10000.0, 500,  -10000.0, 10000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_phi_positron_cut":TH2D("dump_plane_bkg_track_theta_track_phi_positron_cut","dump_plane_bkg_track_theta_track_phi_positron_cut; #theta_{p} [rad]; #phi_{p} [rad]", 6400, 1.6, 3.2, 640, -3.2, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_positron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_positron_weighted_cut","dump_plane_bkg_track_r_track_theta_positron_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_positron_rweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_positron_rweighted_cut","dump_plane_bkg_track_r_track_theta_positron_rweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_positron_thetaweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_positron_thetaweighted_cut","dump_plane_bkg_track_r_track_theta_positron_thetaweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_theta_positron_cut":TH2D("dump_plane_bkg_track_r_track_theta_positron_cut","dump_plane_bkg_track_r_track_theta_positron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_vtxz_positron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxz_positron_weighted_cut","dump_plane_bkg_track_r_track_vtxz_positron_weighted_cut; r [mm]; vtx_{z} [mm]", nRBins, rarray, 1000, -2000.0,  2000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxz_positron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxz_positron_weighted_cut","dump_plane_bkg_track_theta_track_vtxz_positron_weighted_cut; #theta_{p} [rad]; vtx_{z} [mm]", 6400, 1.6, 3.2, 1000, -2000.0,  2000.0)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_vtxx_positron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_vtxx_positron_weighted_cut","dump_plane_bkg_track_r_track_vtxx_positron_weighted_cut; r [mm]; vtx_{x} [mm]", nRBins, rarray, 260,  -260.0,  260.0)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_vtxx_positron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_vtxx_positron_weighted_cut","dump_plane_bkg_track_theta_track_vtxx_positron_weighted_cut; #theta_{p} [rad]; vtx_{x} [mm]", 6400, 1.6, 3.2, 260,  -260.0,  260.0)})
    #allHistoDict.update({"dump_plane_bkg_track_vtxz_track_vtxx_positron_weighted_cut":TH2D("dump_plane_bkg_track_vtxz_track_vtxx_positron_weighted_cut","dump_plane_bkg_track_vtxz_track_vtxx_positron_weighted_cut; vtx_{z} [mm]; vtx_{x} [mm]", 1000, -2000.0, 2000.0, 260,  -260.0,  260.0)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_E_positron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_E_positron_weighted_cut","dump_plane_bkg_track_r_track_E_positron_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_track_E_positron_weighted_cut":TH2D("dump_plane_bkg_track_theta_track_E_positron_weighted_cut","dump_plane_bkg_track_theta_track_E_positron_weighted_cut; #theta_{p} [rad]; E [GeV]", 6400, 1.6, 3.2, nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_r_track_E_positron_cut":TH2D("dump_plane_bkg_track_r_track_E_positron_cut","dump_plane_bkg_track_r_track_E_positron_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    
    
    
    
    
    ### 1D plot
    allHistoDict.update({"dump_plane_bkg_track_time_neutron_cut":TH1D("dump_plane_bkg_track_time_neutron_cut","dump_plane_bkg_track_time_neutron_cut; t [ns]; Events",nbins, tarray)})
    allHistoDict.update({"dump_plane_bkg_track_energy_neutron_cut":TH1D("dump_plane_bkg_track_energy_neutron_cut","dump_plane_bkg_track_energy_neutron_cut; E [GeV]; Events",nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_track_r_neutron_cut":TH1D("dump_plane_bkg_track_r_neutron_cut","dump_plane_bkg_track_r_neutron_cut; r [mm]; Events",nRBins, rarray)})
    allHistoDict.update({"dump_plane_bkg_track_theta_neutron_cut":TH1D("dump_plane_bkg_track_theta_neutron_cut","dump_plane_bkg_track_theta_neutron_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_theta_neutron_weighted_cut":TH1D("dump_plane_bkg_track_theta_neutron_weighted_cut","dump_plane_bkg_track_theta_neutron_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_phi_neutron_cut":TH1D("dump_plane_bkg_track_phi_neutron_cut","dump_plane_bkg_track_phi_neutron_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)})

    
    
    allHistoDict.update({"dump_plane_bkg_track_time_photon_cut":TH1D("dump_plane_bkg_track_time_photon_cut","dump_plane_bkg_track_time_photon_cut; t [ns]; Events",nbins, tarray)})
    allHistoDict.update({"dump_plane_bkg_track_energy_photon_cut":TH1D("dump_plane_bkg_track_energy_photon_cut","dump_plane_bkg_track_energy_photon_cut; E [GeV]; Events",nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_track_r_photon_cut":TH1D("dump_plane_bkg_track_r_photon_cut","dump_plane_bkg_track_r_photon_cut; r [mm]; Events",nRBins, rarray)})
    allHistoDict.update({"dump_plane_bkg_track_theta_photon_weighted_cut":TH1D("dump_plane_bkg_track_theta_photon_weighted_cut","dump_plane_bkg_track_theta_photon_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_theta_photon_cut":TH1D("dump_plane_bkg_track_theta_photon_cut","dump_plane_bkg_track_theta_photon_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_phi_photon_cut":TH1D("dump_plane_bkg_track_phi_photon_cut","dump_plane_bkg_track_phi_photon_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)})
    
    
    #allHistoDict.update({"dump_plane_bkg_track_time_electron_cut":TH1D("dump_plane_bkg_track_time_electron_cut","dump_plane_bkg_track_time_electron_cut; t [ns]; Events",nbins, tarray)})
    #allHistoDict.update({"dump_plane_bkg_track_energy_electron_cut":TH1D("dump_plane_bkg_track_energy_electron_cut","dump_plane_bkg_track_energy_electron_cut; E [GeV]; Events",nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_r_electron_cut":TH1D("dump_plane_bkg_track_r_electron_cut","dump_plane_bkg_track_r_electron_cut; r [mm]; Events",nRBins, rarray)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_electron_weighted_cut":TH1D("dump_plane_bkg_track_theta_electron_weighted_cut","dump_plane_bkg_track_theta_electron_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_electron_cut":TH1D("dump_plane_bkg_track_theta_electron_cut","dump_plane_bkg_track_theta_electron_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_phi_electron_cut":TH1D("dump_plane_bkg_track_phi_electron_cut","dump_plane_bkg_track_phi_electron_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)})
    
    
    #allHistoDict.update({"dump_plane_bkg_track_time_positron_cut":TH1D("dump_plane_bkg_track_time_positron_cut","dump_plane_bkg_track_time_positron_cut; t [ns]; Events",nbins, tarray)})
    #allHistoDict.update({"dump_plane_bkg_track_energy_positron_cut":TH1D("dump_plane_bkg_track_energy_positron_cut","dump_plane_bkg_track_energy_positron_cut; E [GeV]; Events",nbins, xarray)})
    #allHistoDict.update({"dump_plane_bkg_track_r_positron_cut":TH1D("dump_plane_bkg_track_r_positron_cut","dump_plane_bkg_track_r_positron_cut; r [mm]; Events",nRBins, rarray)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_positron_weighted_cut":TH1D("dump_plane_bkg_track_theta_positron_weighted_cut","dump_plane_bkg_track_theta_positron_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_theta_positron_cut":TH1D("dump_plane_bkg_track_theta_positron_cut","dump_plane_bkg_track_theta_positron_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)})
    #allHistoDict.update({"dump_plane_bkg_track_phi_positron_cut":TH1D("dump_plane_bkg_track_phi_positron_cut","dump_plane_bkg_track_phi_positron_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)})
    
    
    
    
    
    lineCounter   = 0
    uniqueTrackList = []

    ### write the bkg as it is
    for lines in bkgFileName.readlines():
        if '#' in lines:
            continue
        
        lineCounter += 1
        if(lineCounter%100000==0): 
            print("processed: ", lineCounter)

        # if(lineCounter > 1000000):
        #     break
        
        ### this is to reject bkg particles from calorimeter
        
        ###### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z << parentid << pxx << pyy << pzz << physicsprocess << time
        
        ##if(lineCounter > 200000):break
        try:
            lines        = lines.rstrip()
            eachWord     = lines.split()
            bxNumber     = int(eachWord[0])
            pdgId        = int(eachWord[1])
            trackId      = int(eachWord[2])
            
            detId        = int(eachWord[3])
            
            if detId != args.det:
                continue
            
            #print("detId: ", detId)
            
            xPos         = float(eachWord[4])
            yPos         = float(eachWord[5])
            energyVal    = float(eachWord[6])
            weight       = float(eachWord[7])
            vtx_x        = float(eachWord[8])
            vtx_y        = float(eachWord[9])
            vtx_z        = float(eachWord[10])
            parent_id    = int(eachWord[11])
            pxx          = float(eachWord[12])
            pyy          = float(eachWord[13])
            pzz          = float(eachWord[14])
            
            #### applying selection cuts to remove unwanted particles
            if(pzz > 0): continue
            
            physprocess  = int(eachWord[15])
            time         = float(eachWord[16])
        except:
            print("Some problem in line: ", lineCounter)
            continue
        
        
        #if(trackId!=1):
            #continue
        #if((pdgId != 2112 or pdgId != 22)):
            #continue
            
        theta, phi   = getTrackThetaPhi(pxx, pyy, pzz)
        thetaPos, phiPos = getTrackThetaPhi(xPos, yPos, zPos)
        rValue       = getR(xPos, yPos)

        # if rValue < 2.01:
        #     rWeight      = 1./(2*math.pi*1.0)
        # else:
        #    rWeight = 1./(2*math.pi*rValue)

        rWeight      = 1./(2*math.pi*rValue)
        
        #### sin theta should not be 0, otherwise the weight will blow up
        # if theta < 0.01:
        #     thetaWeight  = 1./(2*math.pi*math.sin(0.005))
        # elif theta > 3.13:
        #     thetaWeight  = 1./(2*math.pi*math.sin(3.135))
        # else:
        #     thetaWeight  = 1./(2*math.pi*math.sin(theta))

        thetaWeight  = 1./(2*math.pi*math.sin(theta))
            
        
        #bxNumberBranch[0]     = bxNumber
        #pdgIdBranch[0]        = pdgId
        #trackIdBranch[0]      = trackId
        #xPosBranch[0]         = xPos
        #yPosBranch[0]         = yPos
        #energyValBranch[0]    = energyVal
        #weightBranch[0]       = weight
        #vtx_xBranch[0]        = vtx_x
        #vtx_yBranch[0]        = vtx_y
        #vtx_zBranch[0]        = vtx_z
        #parent_idBranch[0]    = parent_id
        #pxxBranch[0]          = pxx
        #pyyBranch[0]          = pyy
        #pzzBranch[0]          = pzz
        #theta_{p}Branch[0]        = theta
        #phi_{p}Branch[0]          = phi
        #physprocessBranch[0]  = physprocess
        #timeBranch[0]         = time
        #rValueBranch[0]       = rValue
        #rValueWeightBranch[0] = rWeight
        
        
        #treeIn.Fill()
        
        
        #### neutrons
        #if(pdgId == 2112):
            #allHistoDict["dump_plane_bkg_track_x_track_y_neutron"].Fill(xPos, yPos, weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_phi_neutron"].Fill(theta, phi, weight)
            #allHistoDict["dump_plane_bkg_track_theta_neutron"].Fill(theta, weight)
            #allHistoDict["dump_plane_bkg_track_phi_neutron"].Fill(phi, weight)
            #allHistoDict["dump_plane_bkg_track_time_neutron"].Fill(time, weight)
            #allHistoDict["dump_plane_bkg_track_r_neutron"].Fill(rValue, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_energy_neutron"].Fill(energyVal, weight)
            
        #### photon
        #if(pdgId == 22):
            #allHistoDict["dump_plane_bkg_track_x_track_y_photon"].Fill(xPos, yPos, weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_phi_photon"].Fill(theta, phi, weight)
            #allHistoDict["dump_plane_bkg_track_theta_photon"].Fill(theta, weight)
            #allHistoDict["dump_plane_bkg_track_phi_photon"].Fill(phi, weight)
            #allHistoDict["dump_plane_bkg_track_time_photon"].Fill(time, weight)
            #allHistoDict["dump_plane_bkg_track_r_photon"].Fill(rValue, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_energy_photon"].Fill(energyVal, weight)
            
            
        
        #### select unique tracks only
        #uniqueTrack = str(bxNumber)+"_"+str(trackId)
        ##if(uniqueTrack in uniqueTrackList): continue
        ##uniqueTrackList.append(uniqueTrack)
        ##if(time > 1000): continue
    
    
    
        ### neutrons
        if(pdgId == 2112):
            # allHistoDict["dump_plane_bkg_track_x_track_y_neutron_cut"].Fill(xPos, yPos, weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_phi_neutron_cut"].Fill(theta, phi, weight)
            allHistoDict["dump_plane_bkg_track_theta_neutron_cut"].Fill(theta, weight)
            allHistoDict["dump_plane_bkg_track_theta_neutron_weighted_cut"].Fill(theta, thetaWeight*weight)
            allHistoDict["dump_plane_bkg_track_phi_neutron_cut"].Fill(phi, weight)
            allHistoDict["dump_plane_bkg_track_time_neutron_cut"].Fill(time, weight)
            allHistoDict["dump_plane_bkg_track_r_neutron_cut"].Fill(rValue, rWeight*weight)
            allHistoDict["dump_plane_bkg_track_energy_neutron_cut"].Fill(energyVal, weight)
            allHistoDict["dump_plane_bkg_track_r_track_theta_neutron_cut"].Fill(rValue, theta, weight)
            allHistoDict["dump_plane_bkg_track_r_track_theta_neutron_weighted_cut"].Fill(rValue, theta, rWeight*thetaWeight*weight)
            allHistoDict["dump_plane_bkg_track_r_track_E_neutron_cut"].Fill(rValue, energyVal, weight)
            allHistoDict["dump_plane_bkg_track_r_track_E_neutron_weighted_cut"].Fill(rValue, energyVal, rWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_E_neutron_weighted_cut"].Fill(theta, energyVal, thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_r_track_vtxz_neutron_weighted_cut"].Fill(rValue, vtx_z, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_vtxz_neutron_weighted_cut"].Fill(theta, vtx_z, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_r_track_vtxx_neutron_weighted_cut"].Fill(rValue, vtx_x, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_vtxx_neutron_weighted_cut"].Fill(theta, vtx_x, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_vtxz_track_vtxx_neutron_weighted_cut"].Fill(vtx_z, vtx_x, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_pos_theta_neutron_cut"].Fill(theta, thetaPos, weight)
            allHistoDict["dump_plane_bkg_track_phi_pos_phi_neutron_cut"].Fill(phi, phiPos, weight)
            allHistoDict["dump_plane_bkg_time_track_E_neutron_cut"].Fill(time,energyVal, weight)
            allHistoDict["dump_plane_bkg_time_track_r_neutron_cut"].Fill(time,rValue, weight)
            allHistoDict["dump_plane_bkg_time_track_theta_neutron_cut"].Fill(time,theta, weight)
            allHistoDict["dump_plane_bkg_time_track_r_neutron_weighted_cut"].Fill(time,rValue, rWeight*weight)
            allHistoDict["dump_plane_bkg_time_track_theta_neutron_weighted_cut"].Fill(time,theta, thetaWeight*weight)
            
            
            
        ### photon
        if(pdgId == 22):
            # allHistoDict["dump_plane_bkg_track_x_track_y_photon_cut"].Fill(xPos, yPos, weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_phi_photon_cut"].Fill(theta, phi, weight)
            allHistoDict["dump_plane_bkg_track_theta_photon_cut"].Fill(theta, weight)
            allHistoDict["dump_plane_bkg_track_theta_photon_weighted_cut"].Fill(theta, thetaWeight*weight)
            allHistoDict["dump_plane_bkg_track_phi_photon_cut"].Fill(phi, weight)
            allHistoDict["dump_plane_bkg_track_time_photon_cut"].Fill(time, weight)
            allHistoDict["dump_plane_bkg_track_r_photon_cut"].Fill(rValue, rWeight*weight)
            allHistoDict["dump_plane_bkg_track_energy_photon_cut"].Fill(energyVal, weight)
            allHistoDict["dump_plane_bkg_track_r_track_theta_photon_cut"].Fill(rValue, theta, weight)
            allHistoDict["dump_plane_bkg_track_r_track_theta_photon_weighted_cut"].Fill(rValue, theta, rWeight*thetaWeight*weight)
            allHistoDict["dump_plane_bkg_track_r_track_E_photon_cut"].Fill(rValue, energyVal, weight)
            allHistoDict["dump_plane_bkg_track_r_track_E_photon_weighted_cut"].Fill(rValue, energyVal, rWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_E_photon_weighted_cut"].Fill(theta, energyVal, thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_r_track_vtxz_photon_weighted_cut"].Fill(rValue, vtx_z, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_vtxz_photon_weighted_cut"].Fill(theta, vtx_z, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_r_track_vtxx_photon_weighted_cut"].Fill(rValue, vtx_x, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_track_vtxx_photon_weighted_cut"].Fill(theta, vtx_x, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_vtxz_track_vtxx_photon_weighted_cut"].Fill(vtx_z, vtx_x, rWeight*thetaWeight*weight)
            # allHistoDict["dump_plane_bkg_track_theta_pos_theta_photon_cut"].Fill(theta, thetaPos, weight)
            allHistoDict["dump_plane_bkg_track_phi_pos_phi_photon_cut"].Fill(phi, phiPos, weight)
            allHistoDict["dump_plane_bkg_time_track_E_photon_cut"].Fill(time,energyVal, weight)
            allHistoDict["dump_plane_bkg_time_track_r_photon_cut"].Fill(time,rValue, weight)
            allHistoDict["dump_plane_bkg_time_track_theta_photon_cut"].Fill(time,theta, weight)
            allHistoDict["dump_plane_bkg_time_track_r_photon_weighted_cut"].Fill(time,rValue, rWeight*weight)
            allHistoDict["dump_plane_bkg_time_track_theta_photon_weighted_cut"].Fill(time,theta, thetaWeight*weight)
            
            
            
        #### positron
        #if(pdgId == -11):
            #allHistoDict["dump_plane_bkg_track_x_track_y_positron_cut"].Fill(xPos, yPos, weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_phi_positron_cut"].Fill(theta, phi, weight)
            #allHistoDict["dump_plane_bkg_track_theta_positron_cut"].Fill(theta, weight)
            #allHistoDict["dump_plane_bkg_track_theta_positron_weighted_cut"].Fill(theta, thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_phi_positron_cut"].Fill(phi, weight)
            #allHistoDict["dump_plane_bkg_track_time_positron_cut"].Fill(time, weight)
            #allHistoDict["dump_plane_bkg_track_r_positron_cut"].Fill(rValue, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_energy_positron_cut"].Fill(energyVal, weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_positron_cut"].Fill(rValue, theta, weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_positron_weighted_cut"].Fill(rValue, theta, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_positron_rweighted_cut"].Fill(rValue, theta, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_positron_thetaweighted_cut"].Fill(rValue, theta, thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_E_positron_cut"].Fill(rValue, energyVal, weight)
            #allHistoDict["dump_plane_bkg_track_r_track_E_positron_weighted_cut"].Fill(rValue, energyVal, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_E_positron_weighted_cut"].Fill(theta, energyVal, thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_vtxz_positron_weighted_cut"].Fill(rValue, vtx_z, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_vtxz_positron_weighted_cut"].Fill(theta, vtx_z, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_vtxx_positron_weighted_cut"].Fill(rValue, vtx_x, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_vtxx_positron_weighted_cut"].Fill(theta, vtx_x, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_vtxz_track_vtxx_positron_weighted_cut"].Fill(vtx_z, vtx_x, rWeight*thetaWeight*weight)
            
            
            
        #### electron
        #if(pdgId == 11):
            #allHistoDict["dump_plane_bkg_track_x_track_y_electron_cut"].Fill(xPos, yPos, weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_phi_electron_cut"].Fill(theta, phi, weight)
            #allHistoDict["dump_plane_bkg_track_theta_electron_cut"].Fill(theta, weight)
            #allHistoDict["dump_plane_bkg_track_theta_electron_weighted_cut"].Fill(theta, thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_phi_electron_cut"].Fill(phi, weight)
            #allHistoDict["dump_plane_bkg_track_time_electron_cut"].Fill(time, weight)
            #allHistoDict["dump_plane_bkg_track_r_electron_cut"].Fill(rValue, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_energy_electron_cut"].Fill(energyVal, weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_electron_cut"].Fill(rValue, theta, weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_electron_weighted_cut"].Fill(rValue, theta, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_electron_rweighted_cut"].Fill(rValue, theta, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_theta_electron_thetaweighted_cut"].Fill(rValue, theta, thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_E_electron_cut"].Fill(rValue, energyVal, weight)
            #allHistoDict["dump_plane_bkg_track_r_track_E_electron_weighted_cut"].Fill(rValue, energyVal, rWeight*weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_E_electron_weighted_cut"].Fill(theta, energyVal, thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_vtxz_electron_weighted_cut"].Fill(rValue, vtx_z, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_vtxz_electron_weighted_cut"].Fill(theta, vtx_z, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_r_track_vtxx_electron_weighted_cut"].Fill(rValue, vtx_x, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_theta_track_vtxx_electron_weighted_cut"].Fill(theta, vtx_x, rWeight*thetaWeight*weight)
            #allHistoDict["dump_plane_bkg_track_vtxz_track_vtxx_electron_weighted_cut"].Fill(vtx_z, vtx_x, rWeight*thetaWeight*weight)
            
            
    
        
        
    for keys in allHistoDict:
        allHistoDict[keys].Scale(1./nbx)
        allHistoDict[keys].Write()
        
    #treeIn.Write()
    outFile.Close()
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
