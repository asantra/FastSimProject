#### run: python3 makeDumpParticlesFromHistogram.py -l  <event file in text mode>  -ev <bx number> -nphot <mean number of photons> -nntrn <mean number of neutrons>
#### this code prepares the fastsim particles from histogram
#### this is using dump geometry
import os
import sys
import time
import pprint
import math
import random
from array import array
import ctypes
from ROOT import *
from collections import OrderedDict
import argparse


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

def main():
    rnd = TRandom()
    rnd.SetSeed()

    parser = argparse.ArgumentParser(description='Code to get 2D plots')
    parser.add_argument('-l', action="store", dest="inFile", type=str, default="RestrictedDumpOnlyFiles_DetId33_trackInfo.root")
    parser.add_argument('-nphot', action="store", dest="nphot", type=int, default=7696824) ### average number of photons per BX
    parser.add_argument('-nntrn', action="store", dest="nntrn", type=int, default=1076188) ### average number of neutrons per BX
    parser.add_argument('-p', action="store", dest="part", type=int, default=1) ### average number of neutrons per BX
    parser.add_argument('-b', action="store", dest="bx", type=int, default=1)
    args = parser.parse_args()
    
    n_neutron   = int(random.gauss(args.nntrn, 1))    ### randomly smearing the number of neutrons per BX
    n_photon    = int(random.gauss(args.nphot, 1))   ### randomly smearing the number of photons per BX

    print("n_neutron: ", n_neutron, " and n_photon: ", n_photon)

    inDir = '/Users/arkasantra/arka/Sasha_Work/OutputFile'
    
    
    withoutText   = args.inFile.split('.root')[0]
    #rootFile      = withoutText+"_MainFullSimFile.root"
    rootFile      = withoutText+".root"
    inFile        = TFile(inDir+"/"+rootFile, "READ")

    #### load histograms from the inRoot file
    dump_plane_bkg_track_r_track_theta_neutron_cut          = inFile.Get("dump_plane_bkg_track_r_track_theta_neutron_cut")
    dump_plane_bkg_track_r_track_theta_neutron_cut.Smooth()
    dump_plane_bkg_track_energy_neutron_cut                 = inFile.Get("dump_plane_bkg_track_energy_neutron_cut")
    dump_plane_bkg_track_phi_pos_phi_neutron_cut            = inFile.Get("dump_plane_bkg_track_phi_pos_phi_neutron_cut")
    dump_plane_bkg_track_phi_pos_phi_neutron_cut.Smooth()
    dump_plane_bkg_time_track_E_neutron_cut                 = inFile.Get("dump_plane_bkg_time_track_E_neutron_cut")
    dump_plane_bkg_time_track_E_neutron_cut.Smooth()

    dump_plane_bkg_track_r_track_theta_photon_cut           = inFile.Get("dump_plane_bkg_track_r_track_theta_photon_cut")
    dump_plane_bkg_track_r_track_theta_photon_cut.Smooth()
    dump_plane_bkg_track_energy_photon_cut                  = inFile.Get("dump_plane_bkg_track_energy_photon_cut")
    dump_plane_bkg_track_phi_pos_phi_photon_cut             = inFile.Get("dump_plane_bkg_track_phi_pos_phi_photon_cut")
    dump_plane_bkg_track_phi_pos_phi_photon_cut.Smooth()
    dump_plane_bkg_track_time_photon_cut                    = inFile.Get("dump_plane_bkg_track_time_photon_cut")
    
    allHistoDict = {}
    nbins  = 450
    xbins  = energyBins(nbins)
    xarray = array('d',xbins)


    tbins  = timeBins(nbins)
    tarray = array('d',tbins)

    nRBins = 6300
    rbins = rBins(nRBins)
    rarray = array('d',rbins)
    
    
    ### output root file containing the photon and neutron information
    # outFile       = TFile(inDir+'/'+withoutText+"_RandomGeneration_v7.root", "RECREATE")
    outFile       = TFile(inDir+'/'+withoutText+"_RandomGeneration_v7_Part"+str(args.part)+".root", "RECREATE")
    outFile.cd()
    Tracks        = TTree("Tracks", "Tracks")


    eventid         = array('i', [0])
    trackid         = std.vector('int')();
    ptrackid        = std.vector('int')();
    detid           = std.vector('int')();
    pdg             = std.vector('int')();
    physproc        = std.vector('int')();
    #r               = std.vector('double')();
    #phi_position    = std.vector('double')();
    E               = std.vector('double')();
    x               = std.vector('double')();
    y               = std.vector('double')();
    z               = std.vector('double')();
    t               = std.vector('double')();
    vtxx            = std.vector('double')();
    vtxy            = std.vector('double')();
    vtxz            = std.vector('double')();
    px              = std.vector('double')();
    py              = std.vector('double')();
    pz              = std.vector('double')();
    theta           = std.vector('double')();
    phi             = std.vector('double')();
    xlocal          = std.vector('double')();
    ylocal          = std.vector('double')();
    zlocal          = std.vector('double')();
    weight          = array('d', [0])
    nsecondary      = std.vector('int')();
    esecondary      = std.vector('double')();
    
    
    Tracks.Branch('eventid', eventid, 'eventid/I');
    Tracks.Branch('trackid', trackid);
    Tracks.Branch('ptrackid', ptrackid);
    Tracks.Branch('detid', detid);
    Tracks.Branch('pdg', pdg);
    Tracks.Branch('physproc', physproc);
    #Tracks.Branch('r', r);
    #Tracks.Branch('phi_position', phi_position);
    Tracks.Branch('E', E);
    Tracks.Branch('x', x);
    Tracks.Branch('y', y);
    Tracks.Branch('z', z);
    Tracks.Branch('t', t);
    Tracks.Branch('vtxx', vtxx);
    Tracks.Branch('vtxy', vtxy);
    Tracks.Branch('vtxz', vtxz);
    Tracks.Branch('px', px);
    Tracks.Branch('py', py);
    Tracks.Branch('pz', pz);
    Tracks.Branch('theta', theta);
    Tracks.Branch('phi', phi);
    Tracks.Branch('xlocal', xlocal);
    Tracks.Branch('ylocal', ylocal);
    Tracks.Branch('zlocal', zlocal);
    Tracks.Branch('weight', weight, 'weight/D');
    Tracks.Branch('nsecondary', nsecondary);
    Tracks.Branch('esecondary', esecondary);

    allHistoDict.update({"dump_plane_bkg_track_E_neutron":TH1D("dump_plane_bkg_track_E_neutron","dump_plane_bkg_track_E_neutron; E [GeV]; Events", nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_track_E_photon":TH1D("dump_plane_bkg_track_E_photon","dump_plane_bkg_track_E_photon; E [GeV]; Events", nbins, xarray)})
    
    allHistoDict.update({"dump_plane_bkg_track_E_time_neutron":TH2D("dump_plane_bkg_track_E_time_neutron","dump_plane_bkg_track_E_time_neutron; t [ns]; E [GeV]", nbins, tarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_track_time_neutron":TH1D("dump_plane_bkg_track_time_neutron","dump_plane_bkg_track_time_neutron; t [ns]; Events", nbins, tarray)})


    allHistoDict.update({"dump_plane_bkg_track_phi_pos_phi_neutron":TH2D("dump_plane_bkg_track_phi_pos_phi_neutron","dump_plane_bkg_track_phi_pos_phi_neutron; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_phi_pos_phi_photon":TH2D("dump_plane_bkg_track_phi_pos_phi_photon","dump_plane_bkg_track_phi_pos_phi_photon; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)})
    
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_neutron_weighted":TH2D("dump_plane_bkg_track_r_track_theta_neutron_weighted","dump_plane_bkg_track_r_track_theta_neutron_weighted; r [mm]; #theta_{p} [rad]",nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_weighted":TH2D("dump_plane_bkg_track_r_track_theta_photon_weighted","dump_plane_bkg_track_r_track_theta_photon_weighted; r [mm]; #theta_{p} [rad]",nRBins, rarray, 6400, 1.6, 3.2)})
    
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_neutron":TH2D("dump_plane_bkg_track_r_track_theta_neutron","dump_plane_bkg_track_r_track_theta_neutron; r [mm]; #theta_{p} [rad]",nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon":TH2D("dump_plane_bkg_track_r_track_theta_photon","dump_plane_bkg_track_r_track_theta_photon; r [mm]; #theta_{p} [rad]",nRBins, rarray, 6400, 1.6, 3.2)})
    
    allHistoDict.update({"dump_plane_bkg_track_theta_neutron_weighted":TH1D("dump_plane_bkg_track_theta_neutron_weighted","dump_plane_bkg_track_theta_neutron_weighted; #theta_{p} [rad]; Events", 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_theta_photon_weighted":TH1D("dump_plane_bkg_track_theta_photon_weighted","dump_plane_bkg_track_theta_photon_weighted; #theta_{p} [rad]; Events", 6400, 1.6, 3.2)})
    
    allHistoDict.update({"dump_plane_bkg_track_r_neutron_weighted":TH1D("dump_plane_bkg_track_r_neutron_weighted","dump_plane_bkg_track_r_neutron_weighted; r [mm]; Events", nRBins, rarray)})
    allHistoDict.update({"dump_plane_bkg_track_r_photon_weighted":TH1D("dump_plane_bkg_track_r_photon_weighted","dump_plane_bkg_track_r_photon_weighted; r [mm]; Events", nRBins, rarray)})
    
    allHistoDict.update({"dump_plane_bkg_track_theta_neutron":TH1D("dump_plane_bkg_track_theta_neutron","dump_plane_bkg_track_theta_neutron; #theta_{p} [rad]; Events", 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_theta_photon":TH1D("dump_plane_bkg_track_theta_photon","dump_plane_bkg_track_theta_photon; #theta_{p} [rad]; Events", 6400, 1.6, 3.2)})
    
    allHistoDict.update({"dump_plane_bkg_track_r_neutron":TH1D("dump_plane_bkg_track_r_neutron","dump_plane_bkg_track_r_neutron; r [mm]; Events", nRBins, rarray)})
    allHistoDict.update({"dump_plane_bkg_track_r_photon":TH1D("dump_plane_bkg_track_r_photon","dump_plane_bkg_track_r_photon; r [mm]; Events", nRBins, rarray)})
    
    allHistoDict.update({"dump_plane_bkg_track_time_photon":TH1D("dump_plane_bkg_track_time_photon","dump_plane_bkg_track_time_photon; t [ns]; Events", nbins, tarray)})




    for nevents in range(args.bx):
        print("BX: ", nevents)
        ### clearing the vector
        
        trackid.clear()         
        ptrackid.clear()        
        detid.clear()           
        pdg.clear()             
        physproc.clear()        
        #r.clear()               
        #phi_position.clear()    
        E.clear()
        x.clear()
        y.clear()
        z.clear()             
        t.clear()               
        vtxx.clear()            
        vtxy.clear()            
        vtxz.clear()            
        px.clear()              
        py.clear()              
        pz.clear()              
        theta.clear()           
        phi.clear()
        xlocal.clear()            
        ylocal.clear()            
        zlocal.clear()
        nsecondary.clear()
        esecondary.clear()            
        
        
        
        #### filling the branches
        eventid[0] = nevents
        weight[0]  = 1.0      
        for nneutron in range(n_neutron):
            if(nneutron%10000==0): print("neutron: working on ", nneutron)
            #### work with neutrons
            p_neutron       = dump_plane_bkg_track_energy_neutron_cut.GetRandom()
            m_neutron       = 0.939565 # GeV
            #E_neutron       = p_neutron #math.sqrt(p_neutron**2+m_neutron**2) ##3 we only want kinematic energy
            r_neutron       = ctypes.c_double(0)
            thetap_neutron  = ctypes.c_double(0)
            phi_neutron     = ctypes.c_double(0)
            phip_neutron    = ctypes.c_double(0)
            #dump_plane_bkg_track_r_track_theta_neutron_weighted_cut.GetRandom2(r_neutron, thetap_neutron) ### from weighted distribution
            
            dump_plane_bkg_track_r_track_theta_neutron_cut.GetRandom2(r_neutron, thetap_neutron) ### from unweighted distribution
            dump_plane_bkg_track_phi_pos_phi_neutron_cut.GetRandom2(phip_neutron, phi_neutron) 
            # p_neutron = math.sqrt(E_neutron**2 - m_neutron**2)
            pxx            = p_neutron*math.sin(thetap_neutron.value)*math.cos(phip_neutron.value)
            pyy            = p_neutron*math.sin(thetap_neutron.value)*math.sin(phip_neutron.value)
            pzz            = p_neutron*math.cos(thetap_neutron.value)
            vtxx_pos       = r_neutron.value*math.cos(phi_neutron.value)
            vtxy_pos       = r_neutron.value*math.sin(phi_neutron.value)
            time_neutron   = ctypes.c_double(0)
            energy_neutron = ctypes.c_double(0)
            dump_plane_bkg_time_track_E_neutron_cut.GetRandom2(time_neutron, energy_neutron)

            E_neutron      = energy_neutron.value
            trackid.push_back(nneutron)
            ptrackid.push_back(-10)
            detid.push_back(-10)
            pdg.push_back(2112)
            physproc.push_back(7000)
            #r.push_back(r_neutron.value)
            #phi_position.push_back(phi_neutron.value);
            E.push_back(E_neutron);
            x.push_back(vtxx_pos)
            y.push_back(vtxy_pos)
            z.push_back(-350.0)
            # t.push_back(random.uniform(0, 1000));  ### uniform timing upto 1000 ns
            t.push_back(time_neutron.value)
            vtxx.push_back(vtxx_pos);
            vtxy.push_back(vtxy_pos);
            vtxz.push_back(-350.0);
            px.push_back(pxx);
            py.push_back(pyy);
            pz.push_back(pzz);
            theta.push_back(thetap_neutron.value);
            phi.push_back(phip_neutron.value);
            xlocal.push_back(vtxx_pos)
            ylocal.push_back(vtxy_pos)
            zlocal.push_back(-350.0)
            nsecondary.push_back(0)
            esecondary.push_back(0.0)
            
            if r_neutron.value < 2.01:
                rWeight      = 1./(2*math.pi*1.0)
            else:
                rWeight      = 1./(2*math.pi*r_neutron.value)
        
        
            #### sin theta should not be 0, otherwise the weight will blow up
            if thetap_neutron.value < 0.01:
                thetaWeight  = 1./(2*math.pi*math.sin(0.005))
            elif thetap_neutron.value > 3.13:
                thetaWeight  = 1./(2*math.pi*math.sin(3.135))
            else:
                thetaWeight  = 1./(2*math.pi*math.sin(thetap_neutron.value))

            allHistoDict["dump_plane_bkg_track_E_neutron"].Fill(E_neutron)
            allHistoDict["dump_plane_bkg_track_time_neutron"].Fill(time_neutron.value)
            allHistoDict["dump_plane_bkg_track_E_time_neutron"].Fill(time_neutron.value, E_neutron)
            allHistoDict["dump_plane_bkg_track_theta_neutron_weighted"].Fill(thetap_neutron.value, thetaWeight)
            allHistoDict["dump_plane_bkg_track_r_neutron_weighted"].Fill(r_neutron.value, rWeight)
            allHistoDict["dump_plane_bkg_track_theta_neutron"].Fill(thetap_neutron.value)
            allHistoDict["dump_plane_bkg_track_r_neutron"].Fill(r_neutron.value)
            allHistoDict["dump_plane_bkg_track_phi_pos_phi_neutron"].Fill(phip_neutron.value, phi_neutron.value)
            allHistoDict["dump_plane_bkg_track_r_track_theta_neutron_weighted"].Fill(r_neutron.value, thetap_neutron.value, rWeight*thetaWeight)
            allHistoDict["dump_plane_bkg_track_r_track_theta_neutron"].Fill(r_neutron.value, thetap_neutron.value)


        for nphoton in range(n_photon):
            if(nphoton%10000==0): print("photon: working on ", nphoton)
            #### work with photons
            p_photon       = dump_plane_bkg_track_energy_photon_cut.GetRandom()
            m_photon       = 0.0 # GeV
            r_photon       = ctypes.c_double(0)
            thetap_photon  = ctypes.c_double(0)
            phi_photon     = ctypes.c_double(0)
            phip_photon    = ctypes.c_double(0)
            #dump_plane_bkg_track_r_track_theta_photon_weighted_cut.GetRandom2(r_photon, thetap_photon)
            dump_plane_bkg_track_r_track_theta_photon_cut.GetRandom2(r_photon, thetap_photon)
            dump_plane_bkg_track_phi_pos_phi_photon_cut.GetRandom2(phip_photon, phi_photon)
            E_photon       = math.sqrt(p_photon**2 + m_photon**2)
            pxx            = p_photon*math.sin(thetap_photon.value)*math.cos(phip_photon.value)
            pyy            = p_photon*math.sin(thetap_photon.value)*math.sin(phip_photon.value)
            pzz            = p_photon*math.cos(thetap_photon.value)
            vtxx_pos       = r_photon.value*math.cos(phi_photon.value)
            vtxy_pos       = r_photon.value*math.sin(phi_photon.value)

            trackid.push_back(nphoton)
            ptrackid.push_back(-10)
            detid.push_back(-10)
            pdg.push_back(22)
            physproc.push_back(7000)
            #r.push_back(r_photon.value)
            #phi_position.push_back(phi_photon.value);
            #r.push_back(r_neutron.value)
            #phi_position.push_back(phi_neutron.value);
            time_photon = dump_plane_bkg_track_time_photon_cut.GetRandom()
            E.push_back(E_photon);
            x.push_back(vtxx_pos)
            y.push_back(vtxy_pos)
            z.push_back(-350.0)
            # t.push_back(random.uniform(0, 1000));  ### uniform timing upto 1000 ns
            t.push_back(time_photon);  ### uniform timing upto 1000 ns
            # vtxx.push_back(vtxx_pos);
            vtxx.push_back(vtxx_pos);
            vtxy.push_back(vtxy_pos);
            vtxz.push_back(-350.0);
            px.push_back(pxx);
            py.push_back(pyy);
            pz.push_back(pzz);
            theta.push_back(thetap_photon.value);
            phi.push_back(phip_photon.value);
            xlocal.push_back(vtxx_pos)
            ylocal.push_back(vtxy_pos)
            zlocal.push_back(-350.0)
            nsecondary.push_back(0)
            esecondary.push_back(0.0)
            
            if r_photon.value < 2.01:
                rWeight      = 1./(2*math.pi*1.0)
            else:
                rWeight      = 1./(2*math.pi*r_photon.value)
        
            #### sin theta should not be 0, otherwise the weight will blow up
            if thetap_photon.value < 0.01:
                thetaWeight  = 1./(2*math.pi*math.sin(0.005))
            elif thetap_photon.value > 3.13:
                thetaWeight  = 1./(2*math.pi*math.sin(3.135))
            else:
                thetaWeight  = 1./(2*math.pi*math.sin(thetap_photon.value))


            allHistoDict["dump_plane_bkg_track_E_photon"].Fill(p_photon)
            allHistoDict["dump_plane_bkg_track_time_photon"].Fill(time_photon)
            allHistoDict["dump_plane_bkg_track_theta_photon_weighted"].Fill(thetap_photon.value, thetaWeight)
            allHistoDict["dump_plane_bkg_track_r_photon_weighted"].Fill(r_photon.value, rWeight)
            allHistoDict["dump_plane_bkg_track_theta_photon"].Fill(thetap_photon.value)
            allHistoDict["dump_plane_bkg_track_r_photon"].Fill(r_photon.value)
            allHistoDict["dump_plane_bkg_track_phi_pos_phi_photon"].Fill(phip_photon.value, phi_photon.value)
            allHistoDict["dump_plane_bkg_track_r_track_theta_photon_weighted"].Fill(r_photon.value, thetap_photon.value, rWeight*thetaWeight)
            allHistoDict["dump_plane_bkg_track_r_track_theta_photon"].Fill(r_photon.value, thetap_photon.value)
            
        Tracks.Fill()
        
    Tracks.Write()
    
    for keys in allHistoDict:
        allHistoDict[keys].Scale(1./args.bx)
        # allHistoDict[keys].Scale(1./10.0)
        allHistoDict[keys].Write()
    outFile.Close()
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
