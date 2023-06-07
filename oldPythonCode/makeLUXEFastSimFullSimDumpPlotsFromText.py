#### this code mixes the background tracks and signal tracks per bunch crossing
#### run: python3 plotBackgroundFilesFromText.py -l  <event file in text mode>  -b <bx number>
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
    
    nNt = 0
    nPh = 0
    nNtLim = 999999999999
    nPhLim = 999999999999

    if (args.det==32):
       nNtLim = 345114
       nPhLim = 105606
    elif (args.det==31):
       nNtLim = 268118
       nPhLim = 97020
    else:
        pass

    inDir = '/Users/arkasantra/arka/Sasha_Work/OutputFile'
    #inDir = "/Volumes/Study/Weizmann_PostDoc/Sasha_work/OutputFile/ReprocessedBkgTracksAfterTDR"
    
    bkgFileName   = open(inDir+'/'+args.inFile)
    withoutText   = args.inFile.split('.txt')[0]
    #rootFile      = withoutText+"_DetId"+str(args.det)+"_timelt1000.root"
    # rootFile      = withoutText+"_DetId"+str(args.det)+"_Egt0p1GeV_rlt300mm.root"
    # rootFile      = withoutText+"_DetId"+str(args.det)+"_rlt300mm.root"
    rootFile      = withoutText+"_DetId"+str(args.det)+".root"
    # rootFile      = withoutText+".root"
    #rootFile      = withoutText+"_timelt1000.root"
    #rootFile      = withoutText+"_MainFullSimFile.root"
    nbx           = args.bx
    
    ### this is Dump only geometry
    # if args.det==33:
    #     zPos = -350
    # elif args.det==32:
    #     zPos = -5000
    # elif args.det==31:
    #     zPos = -10000
    # elif args.det==30:
    #     zPos = -15000
    # else:
    #     zPos = -350

    ### this is LUXE geometry
    if args.det==33:
        zPos = 6621.91;
    elif args.det==32:
        zPos = 5450.25;
    elif args.det==31:
        zPos = 4125;
    elif args.det==30:
        zPos = 0
    else:
        zPos = 6621.91;
        
    print("zPos: ", zPos)
    print('BX selected: ',nbx)
    
    outFile       = TFile(inDir+'/'+rootFile, "RECREATE")
    outFile.cd()
    
    nbins  = 450
    xbins  = energyBins(nbins)
    xarray = array('d',xbins)
    
    tbins  = timeBins(nbins)
    tarray = array('d',tbins)
    

    nRBins = 6300
    rbins = rBins(nRBins)
    rarray = array('d',rbins)
        
    allHistoDict  = {}
    
    
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_neutron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_neutron_weighted_cut","dump_plane_bkg_track_r_track_theta_neutron_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_neutron_cut":TH2D("dump_plane_bkg_track_r_track_theta_neutron_cut","dump_plane_bkg_track_r_track_theta_neutron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_neutron_weighted_cut":TH2D("dump_plane_bkg_track_r_track_E_neutron_weighted_cut","dump_plane_bkg_track_r_track_E_neutron_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_track_phi_pos_phi_neutron_cut":TH2D("dump_plane_bkg_track_phi_pos_phi_neutron_cut","dump_plane_bkg_track_phi_pos_phi_neutron_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_neutron_cut":TH2D("dump_plane_bkg_track_r_track_E_neutron_cut","dump_plane_bkg_track_r_track_E_neutron_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_E_neutron_cut":TH2D("dump_plane_bkg_time_track_E_neutron_cut","dump_plane_bkg_time_track_E_neutron_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_neutron_cut":TH2D("dump_plane_bkg_time_track_r_neutron_cut","dump_plane_bkg_time_track_r_neutron_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_neutron_cut":TH2D("dump_plane_bkg_time_track_theta_neutron_cut","dump_plane_bkg_time_track_theta_neutron_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_neutron_weighted_cut":TH2D("dump_plane_bkg_time_track_r_neutron_weighted_cut","dump_plane_bkg_time_track_r_neutron_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_neutron_weighted_cut":TH2D("dump_plane_bkg_time_track_theta_neutron_weighted_cut","dump_plane_bkg_time_track_theta_neutron_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    
    
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_weighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_weighted_cut","dump_plane_bkg_track_r_track_theta_photon_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_rweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_rweighted_cut","dump_plane_bkg_track_r_track_theta_photon_rweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut","dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_theta_photon_cut":TH2D("dump_plane_bkg_track_r_track_theta_photon_cut","dump_plane_bkg_track_r_track_theta_photon_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_photon_weighted_cut":TH2D("dump_plane_bkg_track_r_track_E_photon_weighted_cut","dump_plane_bkg_track_r_track_E_photon_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_track_phi_pos_phi_photon_cut":TH2D("dump_plane_bkg_track_phi_pos_phi_photon_cut","dump_plane_bkg_track_phi_pos_phi_photon_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)})
    allHistoDict.update({"dump_plane_bkg_track_r_track_E_photon_cut":TH2D("dump_plane_bkg_track_r_track_E_photon_cut","dump_plane_bkg_track_r_track_E_photon_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_E_photon_cut":TH2D("dump_plane_bkg_time_track_E_photon_cut","dump_plane_bkg_time_track_E_photon_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_photon_cut":TH2D("dump_plane_bkg_time_track_r_photon_cut","dump_plane_bkg_time_track_r_photon_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_photon_cut":TH2D("dump_plane_bkg_time_track_theta_photon_cut","dump_plane_bkg_time_track_theta_photon_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    allHistoDict.update({"dump_plane_bkg_time_track_r_photon_weighted_cut":TH2D("dump_plane_bkg_time_track_r_photon_weighted_cut","dump_plane_bkg_time_track_r_photon_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)})
    allHistoDict.update({"dump_plane_bkg_time_track_theta_photon_weighted_cut":TH2D("dump_plane_bkg_time_track_theta_photon_weighted_cut","dump_plane_bkg_time_track_theta_photon_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)})
    
    
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

            # if energyVal < 0.1:
            #     continue

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

        # if rValue > 300:
        #     continue

        if rValue < 1.01:
            rWeight      = 1./(2*math.pi*0.5)
        else:
           rWeight = 1./(2*math.pi*rValue)

        # rWeight      = 1./(2*math.pi*rValue)
        
        #### sin theta should not be 0, otherwise the weight will blow up
        if theta < 0.01:
            thetaWeight  = 1./(2*math.pi*math.sin(0.005))
        elif theta > 3.13:
            thetaWeight  = 1./(2*math.pi*math.sin(3.135))
        else:
            thetaWeight  = 1./(2*math.pi*math.sin(theta))

        # thetaWeight  = 1./(2*math.pi*math.sin(theta))
    
    
        ### neutrons
        if(pdgId == 2112):
            nNt += 1
            if (nNt < nNtLim):
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
                allHistoDict["dump_plane_bkg_track_phi_pos_phi_neutron_cut"].Fill(phi, phiPos, weight)
                allHistoDict["dump_plane_bkg_time_track_E_neutron_cut"].Fill(time,energyVal, weight)
                allHistoDict["dump_plane_bkg_time_track_r_neutron_cut"].Fill(time,rValue, weight)
                allHistoDict["dump_plane_bkg_time_track_theta_neutron_cut"].Fill(time,theta, weight)
                allHistoDict["dump_plane_bkg_time_track_r_neutron_weighted_cut"].Fill(time,rValue, rWeight*weight)
                allHistoDict["dump_plane_bkg_time_track_theta_neutron_weighted_cut"].Fill(time,theta, thetaWeight*weight)
            
            
            
        ### photon
        if(pdgId == 22):
            nPh += 1
            if (nPh < nPhLim):
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
                allHistoDict["dump_plane_bkg_track_phi_pos_phi_photon_cut"].Fill(phi, phiPos, weight)
                allHistoDict["dump_plane_bkg_time_track_E_photon_cut"].Fill(time,energyVal, weight)
                allHistoDict["dump_plane_bkg_time_track_r_photon_cut"].Fill(time,rValue, weight)
                allHistoDict["dump_plane_bkg_time_track_theta_photon_cut"].Fill(time,theta, weight)
                allHistoDict["dump_plane_bkg_time_track_r_photon_weighted_cut"].Fill(time,rValue, rWeight*weight)
                allHistoDict["dump_plane_bkg_time_track_theta_photon_weighted_cut"].Fill(time,theta, thetaWeight*weight)
            

        if (nPh > nPhLim and nNt > nNtLim):
            break
        
    for keys in allHistoDict:
        allHistoDict[keys].Scale(1./nbx)
        allHistoDict[keys].Write()
        
    print("total Neutron: ", nNt)
    print("total Photon: ", nPh)
    #treeIn.Write()
    outFile.Close()
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
