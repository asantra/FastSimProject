//#### run: inside root: .L makeLUXEFastSimFullSimDumpPlotsFromText.C++ && makeLUXEFastSimFullSimDumpPlotsFromText()
//#### Here the histograms are made from FullSim, coming from Sasha
// This is used to run FastSim and FullSim samples at test surfaces. 
// The plots will be later compared, so normalization (number of particles) are important.
// This uses LUXE dump geometry.

#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <map>
#include <chrono>


using namespace std;
using namespace std::chrono;


/// energy and time bins
const int nbins=650;
/// r bins
const int nRBins = 6300;
/// x/y bins
const int nXBins = 6600;


vector<double> getTrackThetaPhi(double x, double y, double z){
    double phi   = atan2(y, x);
    double dist  = sqrt(x*x+y*y+z*z);
    double theta = acos(z/dist);
    vector<double> angles;
    angles.push_back(theta);
    angles.push_back(phi);
    return angles;
}

double getR(double x,double y){
    return sqrt(x*x+y*y);
};


void makeLUXEFastSimFullSimDumpPlotsFromText(string bkgFileName="/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", int det=33, float bx=1.0, bool isFullSim=false){
    
    auto start = high_resolution_clock::now();
    int nNt = 0;
    int nPh = 0;
    long nNtLim = 99999999999999;
    long nPhLim = 99999999999999;

    /// need less number of events from FullSim
    if(isFullSim){
        if(det==33){
            /// this is limited to save time. They should not be used when comparing with the full stat
            /// these numbers should match the numbers used in makeDumpParticlesFromHistogramLUXE
            nPhLim = 1353100;
            nNtLim = 5497735;
            /// if we use full 0.06BX, then these are the numbers:
            // nPhLim=13531005;
            // nNtLim=54977353;
        }
        else if (det==32){
            /// fullsim entries
               nNtLim = 44919187;
               nPhLim = 12559752;
            /// fastsim entries, unweighted
            // nNtLim = 9035488;
            // nPhLim = 1265180;
            /// fastsim entries, weighted
            // nNtLim = 483526;
            // nPhLim = 126613;
        }
        else if (det==31){
            /// fullsim entries
               nNtLim = 34844633;
               nPhLim = 11513244;
            // /// fastsim entries, unweighted
            // nNtLim = 7020570;
            // nPhLim = 1160168;
            /// fastsim entries, unweighted
            // nNtLim = 388300;
            // nPhLim = 116090;
        }
        else
            ;
    }

    /// for local files
    // string inDir = "/Users/arkasantra/arka/Sasha_Work/OutputFile/";
    /// for DESY files
    // string inDir = "/nfs/dust/luxe/user/santraar/TextSamples_October2021/DumpFromFastSim/FastSimFiles_LUXEDump/";
    // string inDir = "/Volumes/Study/Weizmann_PostDoc/Sasha_work/OutputFile/ReprocessedBkgTracksAfterTDR"
    string inDir = "/srv01/agrp/arkas/GANFastSim/";
    

    ifstream inFile;
    
    inFile.open(bkgFileName);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    std::string foutname    = bkgFileName.substr(bkgFileName.find_last_of("/")+1);
    foutname                = foutname.substr(0, foutname.find_last_of("."));
    std::string rootoutname = inDir+foutname;
    if(isFullSim)
        rootoutname             += std::string("_NoECutNtrn_DetId"+to_string(det)+".root");
    else 
        rootoutname             += std::string("_DetId"+to_string(det)+".root");
    
    cout << "The output file: " << rootoutname << endl;

    //### this is LUXE geometry
    double zPos = 0.0;
    if (det==33)
        zPos = 6621.91;
    else if (det==32)
        zPos = 5450.25;
    else if (det==31)
        zPos = 4125;
    else if (det==30)
        zPos = 0;
    else
        zPos = 6621.91;
        
    cout << "zPos: " << zPos << endl;
    cout << "BX selected: " << bx << endl;
    
    TFile *outFile       = new TFile(rootoutname.c_str(), "RECREATE");
    outFile->cd();
    
    //#//// variable binned in E;
    double xarray[nbins + 1]  = {0.0};    
    int logmin                = -12;
    int logmax                = 0;
    double logbinwidth        = (logmax-logmin)/float(nbins);
  
    for(int i=0; i < nbins+1; ++i){
      xarray[i] = pow( 10,(logmin + i*logbinwidth) );
    }
    
    
    //#//// variable binned in t;
    double tarray[nbins + 1]  = {0.0};    
    logmin                = 1;
    logmax                = 12;
    logbinwidth        = (logmax-logmin)/float(nbins);
  
    for(int i=0; i < nbins+1; ++i){
      tarray[i] = pow( 10,(logmin + i*logbinwidth) );
    }
    
    
    double rarray[nRBins+1] = {0.0};
    for(int i=0; i < nRBins+1; ++i){
        if (i < 6000)
            rarray[i] = (i*1);
        else if (i < 6200)
            rarray[i] = (6000+(i-6000)*10);
        else
            rarray[i] = (8000+(i-6200)*20);
    }
    
    double xyarray[nXBins+1] = {0.0};
    for(int i=0; i < nXBins+1; ++i){
        if(i<100)
            xyarray[i] = (-10000+i*20);
        else if(i<300)
            xyarray[i] = (-8000+(i-100)*10);
        else if (i >= 300 && i < 6300)
            xyarray[i] = -6000 + ((i-300)*2);
        else if (i < 6500)
            xyarray[i] = (6000+(i-6300)*10);
        else
            xyarray[i] = (8000+(i-6500)*20);
    }
    
    map<string, TH1D*> allHisto1Dict;
    map<string, TH2D*> allHisto2Dict;
    
    
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_neutron_weighted_cut", new TH2D("dump_plane_bkg_track_r_track_theta_neutron_weighted_cut","dump_plane_bkg_track_r_track_theta_neutron_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_neutron_cut", new TH2D("dump_plane_bkg_track_r_track_theta_neutron_cut","dump_plane_bkg_track_r_track_theta_neutron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_neutron_weighted_cut", new TH2D("dump_plane_bkg_track_r_track_E_neutron_weighted_cut","dump_plane_bkg_track_r_track_E_neutron_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_phi_pos_phi_neutron_cut", new TH2D("dump_plane_bkg_track_phi_pos_phi_neutron_cut","dump_plane_bkg_track_phi_pos_phi_neutron_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_neutron_cut", new TH2D("dump_plane_bkg_track_r_track_E_neutron_cut","dump_plane_bkg_track_r_track_E_neutron_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_E_neutron_cut", new TH2D("dump_plane_bkg_time_track_E_neutron_cut","dump_plane_bkg_time_track_E_neutron_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_neutron_cut", new TH2D("dump_plane_bkg_time_track_r_neutron_cut","dump_plane_bkg_time_track_r_neutron_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_neutron_cut", new TH2D("dump_plane_bkg_time_track_theta_neutron_cut","dump_plane_bkg_time_track_theta_neutron_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_neutron_weighted_cut", new TH2D("dump_plane_bkg_time_track_r_neutron_weighted_cut","dump_plane_bkg_time_track_r_neutron_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_neutron_weighted_cut", new TH2D("dump_plane_bkg_time_track_theta_neutron_weighted_cut","dump_plane_bkg_time_track_theta_neutron_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    
    
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_photon_weighted_cut", new TH2D("dump_plane_bkg_track_r_track_theta_photon_weighted_cut","dump_plane_bkg_track_r_track_theta_photon_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_photon_rweighted_cut", new TH2D("dump_plane_bkg_track_r_track_theta_photon_rweighted_cut","dump_plane_bkg_track_r_track_theta_photon_rweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut", new TH2D("dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut","dump_plane_bkg_track_r_track_theta_photon_thetaweighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_photon_cut", new TH2D("dump_plane_bkg_track_r_track_theta_photon_cut","dump_plane_bkg_track_r_track_theta_photon_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 6400, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_photon_weighted_cut", new TH2D("dump_plane_bkg_track_r_track_E_photon_weighted_cut","dump_plane_bkg_track_r_track_E_photon_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_phi_pos_phi_photon_cut", new TH2D("dump_plane_bkg_track_phi_pos_phi_photon_cut","dump_plane_bkg_track_phi_pos_phi_photon_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_photon_cut", new TH2D("dump_plane_bkg_track_r_track_E_photon_cut","dump_plane_bkg_track_r_track_E_photon_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_E_photon_cut", new TH2D("dump_plane_bkg_time_track_E_photon_cut","dump_plane_bkg_time_track_E_photon_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_photon_cut", new TH2D("dump_plane_bkg_time_track_r_photon_cut","dump_plane_bkg_time_track_r_photon_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_photon_cut", new TH2D("dump_plane_bkg_time_track_theta_photon_cut","dump_plane_bkg_time_track_theta_photon_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_photon_weighted_cut", new TH2D("dump_plane_bkg_time_track_r_photon_weighted_cut","dump_plane_bkg_time_track_r_photon_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_photon_weighted_cut", new TH2D("dump_plane_bkg_time_track_theta_photon_weighted_cut","dump_plane_bkg_time_track_theta_photon_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    
    
    
    
    
    
    
    ///### 1D plot
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_neutron_cut", new TH1D("dump_plane_bkg_track_time_neutron_cut","dump_plane_bkg_track_time_neutron_cut; t [ns]; Events",nbins, tarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_energy_neutron_cut", new TH1D("dump_plane_bkg_track_energy_neutron_cut","dump_plane_bkg_track_energy_neutron_cut; E [GeV]; Events",nbins, xarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_neutron_cut", new TH1D("dump_plane_bkg_track_r_neutron_cut","dump_plane_bkg_track_r_neutron_cut; r [mm]; Events",nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_small_neutron_cut", new TH1D("dump_plane_bkg_track_r_small_neutron_cut","dump_plane_bkg_track_r_small_neutron_cut; r [mm]; Events",500,0, 500)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_neutron_cut", new TH1D("dump_plane_bkg_track_theta_neutron_cut","dump_plane_bkg_track_theta_neutron_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_neutron_weighted_cut", new TH1D("dump_plane_bkg_track_theta_neutron_weighted_cut","dump_plane_bkg_track_theta_neutron_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_neutron_cut", new TH1D("dump_plane_bkg_track_phi_neutron_cut","dump_plane_bkg_track_phi_neutron_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)));

    
    
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_photon_cut", new TH1D("dump_plane_bkg_track_time_photon_cut","dump_plane_bkg_track_time_photon_cut; t [ns]; Events",nbins, tarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_energy_photon_cut", new TH1D("dump_plane_bkg_track_energy_photon_cut","dump_plane_bkg_track_energy_photon_cut; E [GeV]; Events",nbins, xarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_photon_cut", new TH1D("dump_plane_bkg_track_r_photon_cut","dump_plane_bkg_track_r_photon_cut; r [mm]; Events",nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_r_small_photon_cut", new TH1D("dump_plane_bkg_track_r_small_photon_cut","dump_plane_bkg_track_r_small_photon_cut; r [mm]; Events",500, 0, 500)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_photon_weighted_cut", new TH1D("dump_plane_bkg_track_theta_photon_weighted_cut","dump_plane_bkg_track_theta_photon_weighted_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_photon_cut", new TH1D("dump_plane_bkg_track_theta_photon_cut","dump_plane_bkg_track_theta_photon_cut; #theta_{p} [rad]; Events",6400, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_photon_cut", new TH1D("dump_plane_bkg_track_phi_photon_cut","dump_plane_bkg_track_phi_photon_cut; #phi_{p} [rad]; Events",640, -3.2, 3.2)));
    
    long lineCounter   = 0;

    int bxNumber(-99999), pdgId(-99999), trackId(-99999), detId(-99999);
    double xPos(-99999.0), yPos(-99999.0), energyVal(-99999.0), weight(-99999.0);
    double vtx_x(-99999.0), vtx_y(-99999.0), vtx_z(-99999.0);
    int parent_id(-99999), physprocess(-99999);
    double pxx(-99999.0), pyy(-99999.0), pzz(-99999.0), time(-99999.0);
    
    // ###### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z << parentid << pxx << pyy << pzz << physicsprocess << time
        
    while(inFile >> bxNumber >> pdgId >> trackId >> detId >> xPos >> yPos >> energyVal >> weight >> vtx_x >> vtx_y >> vtx_z >> parent_id >> pxx >> pyy >> pzz >> physprocess >> time){
        
        lineCounter += 1;
        if(lineCounter%100000==0)
            cout << "processed: " << lineCounter << endl;
        
        if(pzz > 0) continue;
        if(detId!=det) continue;

        // # if(lineCounter > 1000000)
        // #     break
        
        // ### this is to reject bkg particles from calorimeter
            
        vector<double> angles;
        angles.clear();
        angles          = getTrackThetaPhi(pxx, pyy, pzz);
        double theta           = angles.at(0);
        double phi             = angles.at(1);
        
        vector<double> anglesPos;
        anglesPos.clear();
        anglesPos              = getTrackThetaPhi(xPos, yPos, zPos); 
        double thetaPos        = anglesPos.at(0);
        double phiPos          = anglesPos.at(1); 
        double rValue          = getR(xPos, yPos);
        double rWeight         = 1./(2*TMath::Pi()*rValue);
        double thetaWeight     = 1./(2*TMath::Pi()*TMath::Sin(theta));


        if (rValue < 1.01)
            rWeight      = 1./(2*TMath::Pi()*0.5);
        else
            rWeight      = 1./(2*TMath::Pi()*rValue);

        //# rWeight      = 1./(2*math.pi*rValue)
        
        // #### sin theta should not be 0, otherwise the weight will blow up
        if (theta < 0.01)
            thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(0.005));
        else if (theta > 3.13)
            thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(3.135));
        else
            thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(theta));

        // # thetaWeight  = 1./(2*math.pi*math.sin(theta))


        
    
        //### neutrons
        if(pdgId == 2112){
            nNt += 1;
            if (nNt < nNtLim){
            // if (true){
                allHisto1Dict["dump_plane_bkg_track_theta_neutron_cut"]->Fill(theta, weight);
                allHisto1Dict["dump_plane_bkg_track_theta_neutron_weighted_cut"]->Fill(theta, thetaWeight*weight);
                allHisto1Dict["dump_plane_bkg_track_phi_neutron_cut"]->Fill(phi, weight);
                allHisto1Dict["dump_plane_bkg_track_time_neutron_cut"]->Fill(time, weight);
                allHisto1Dict["dump_plane_bkg_track_r_neutron_cut"]->Fill(rValue, rWeight*weight);
                allHisto1Dict["dump_plane_bkg_track_r_small_neutron_cut"]->Fill(rValue, rWeight*weight);
                allHisto1Dict["dump_plane_bkg_track_energy_neutron_cut"]->Fill(energyVal, weight);


                allHisto2Dict["dump_plane_bkg_track_r_track_theta_neutron_cut"]->Fill(rValue, theta, weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_theta_neutron_weighted_cut"]->Fill(rValue, theta, rWeight*thetaWeight*weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_E_neutron_cut"]->Fill(rValue, energyVal, weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_E_neutron_weighted_cut"]->Fill(rValue, energyVal, rWeight*weight);
                allHisto2Dict["dump_plane_bkg_track_phi_pos_phi_neutron_cut"]->Fill(phi, phiPos, weight);
                allHisto2Dict["dump_plane_bkg_time_track_E_neutron_cut"]->Fill(time,energyVal, weight);
                allHisto2Dict["dump_plane_bkg_time_track_r_neutron_cut"]->Fill(time,rValue, weight);
                allHisto2Dict["dump_plane_bkg_time_track_theta_neutron_cut"]->Fill(time,theta, weight);
                allHisto2Dict["dump_plane_bkg_time_track_r_neutron_weighted_cut"]->Fill(time,rValue, rWeight*weight);
                allHisto2Dict["dump_plane_bkg_time_track_theta_neutron_weighted_cut"]->Fill(time,theta, thetaWeight*weight);
            }
            
        }
            
        //### photon
        if(pdgId == 22){
            nPh += 1;
            if (nPh < nPhLim){
            // if(true){
                allHisto1Dict["dump_plane_bkg_track_theta_photon_cut"]->Fill(theta, weight);
                allHisto1Dict["dump_plane_bkg_track_theta_photon_weighted_cut"]->Fill(theta, thetaWeight*weight);
                allHisto1Dict["dump_plane_bkg_track_phi_photon_cut"]->Fill(phi, weight);
                allHisto1Dict["dump_plane_bkg_track_time_photon_cut"]->Fill(time, weight);
                allHisto1Dict["dump_plane_bkg_track_r_photon_cut"]->Fill(rValue, rWeight*weight);
                allHisto1Dict["dump_plane_bkg_track_r_small_photon_cut"]->Fill(rValue, rWeight*weight);
                allHisto1Dict["dump_plane_bkg_track_energy_photon_cut"]->Fill(energyVal, weight);

                
                allHisto2Dict["dump_plane_bkg_track_r_track_theta_photon_cut"]->Fill(rValue, theta, weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_theta_photon_weighted_cut"]->Fill(rValue, theta, rWeight*thetaWeight*weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_E_photon_cut"]->Fill(rValue, energyVal, weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_E_photon_weighted_cut"]->Fill(rValue, energyVal, rWeight*weight);
                allHisto2Dict["dump_plane_bkg_track_phi_pos_phi_photon_cut"]->Fill(phi, phiPos, weight);
                allHisto2Dict["dump_plane_bkg_time_track_E_photon_cut"]->Fill(time,energyVal, weight);
                allHisto2Dict["dump_plane_bkg_time_track_r_photon_cut"]->Fill(time,rValue, weight);
                allHisto2Dict["dump_plane_bkg_time_track_theta_photon_cut"]->Fill(time,theta, weight);
                allHisto2Dict["dump_plane_bkg_time_track_r_photon_weighted_cut"]->Fill(time,rValue, rWeight*weight);
                allHisto2Dict["dump_plane_bkg_time_track_theta_photon_weighted_cut"]->Fill(time,theta, thetaWeight*weight);
            }
        }

        if (nPh > nPhLim && nNt > nNtLim)
            break;
    }
        
    for (map<string, TH2D*>::iterator it = allHisto2Dict.begin(); it != allHisto2Dict.end(); ++it){
        it->second->Scale(1./bx);
        it->second->Write();
    }
    for (map<string, TH1D*>::iterator it = allHisto1Dict.begin(); it != allHisto1Dict.end(); ++it){
        it->second->Scale(1./bx);
        it->second->Write();
    }
    
    outFile->Close();
    delete outFile;
    std::cout << "total neutron: " << nNt << std::endl;
    std::cout << "total photon: " << nPh << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "time taken for the running: " << duration.count() << " s" << endl;
    
}
