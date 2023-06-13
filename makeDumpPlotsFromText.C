//#### run: inside root: .L makeDumpPlotsFromText.C++ && makeDumpPlotsFromText()
//#### Here the histograms are made from FullSim text files at the sampling surface, coming from Sasha.
// This code prepares the 2D and 1D plots from where we sample.
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
const int nbins=450;
/// r bins
// const int nRBins = 6300;
// coarse binning
const int nRBins = 1800;
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
    

void makeDumpPlotsFromText(string bkgFileName="/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33.txt", float bx=1.0, int det=33, bool cmpPlt=false){
    
    auto start = high_resolution_clock::now();
    // "/nfs/dust/luxe/user/santraar/TextSamples_October2021/DumpFromFastSim/LUXEDumpFiles_FullSim_0p06BX_DetId33.txt"
    string inDir = "/Users/arkasantra/arka/Sasha_Work/OutputFile/";
    
    const double x0 = -92.65;

    
    ifstream inFile;
    
    inFile.open(bkgFileName);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
    
    std::string foutname    = bkgFileName.substr(bkgFileName.find_last_of("/")+1);
    foutname                = foutname.substr(0, foutname.find_last_of("."));
    std::string rootoutname = inDir+foutname;
    if (cmpPlt)
        rootoutname             += std::string("_NoECutNtrn_CoarseBinning_1DComparePlot.root");
    else
        rootoutname             += std::string("_NoECutNtrn_CoarseBinning.root");

    cout << "The output file: " << rootoutname.c_str() << endl;

    long nntrnLim = 549773;
    long nphoLim = 135310;
    
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
    int logmin                = -7;
    int logmax                = 0;
    double logbinwidth        = (logmax-logmin)/float(nbins);
  
    for(int i=0; i < nbins+1; ++i){
      xarray[i] = pow( 10,(logmin + i*logbinwidth) );
    }
    
    
    //#//// variable binned in t;
    double tarray[nbins + 1]  = {0.0};    
    logmin                = 1;
    logmax                = 8;
    logbinwidth        = (logmax-logmin)/float(nbins);
  
    for(int i=0; i < nbins+1; ++i){
      tarray[i] = pow( 10,(logmin + i*logbinwidth) );
    }
    
    /// fine binning in r, hence bin width of 1 upto 6000
    // double rarray[nRBins+1] = {0.0};
    // for(int i=0; i < nRBins+1; ++i){
    //     if (i*4 < 6000)
    //         rarray[i] = (i*4);
    //     else if (i < 6200)
    //         rarray[i] = (6000+(i-6000)*10);
    //     else
    //         rarray[i] = (8000+(i-6200)*20);
    // }
    /// coarse binning in r, hence bin width of 1 upto 6000
    double rarray[nRBins+1] = {0.0};
    for(int i=0; i < nRBins+1; ++i){
        if (i*4 < 6000)
            rarray[i] = (i*4);
        else if (i >= 1500 && i < 1700)
            rarray[i] = (6000+(i-1500)*10);
        else
            rarray[i] = (8000+(i-1700)*20);
    }
    // cout << "The r bins: " << endl;
    // for(int nRit=0; nRit < nRBins+1; nRit++)
    //    cout << "rarray["<< nRit << "]: " << rarray[nRit] << endl;
    
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
    
    //////// histograms after applying some cuts;
    
    ////// 2D plot;
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_x_track_y_neutron_cut",  new TH2D("dump_plane_bkg_track_x_track_y_neutron_cut","dump_plane_bkg_track_x_track_y_neutron_cut; track x [mm]; track y [mm]",nXBins, xyarray, nXBins, xyarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_neutron_weighted_cut",  new TH2D("dump_plane_bkg_track_r_track_theta_neutron_weighted_cut","dump_plane_bkg_track_r_track_theta_neutron_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_neutron_cut",  new TH2D("dump_plane_bkg_track_rUp_track_theta_neutron_cut","dump_plane_bkg_track_rUp_track_theta_neutron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_neutron_cut",  new TH2D("dump_plane_bkg_track_rDn_track_theta_neutron_cut","dump_plane_bkg_track_rDn_track_theta_neutron_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_backward_neutron_cut",  new TH2D("dump_plane_bkg_track_rUp_track_theta_backward_neutron_cut","dump_plane_bkg_track_rUp_track_theta_backward_neutron_cut; r [mm]; #theta_{p} [rad]", 250, 0, 1000, 200, 2.8, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_backward_neutron_cut",  new TH2D("dump_plane_bkg_track_rDn_track_theta_backward_neutron_cut","dump_plane_bkg_track_rDn_track_theta_backward_neutron_cut; r [mm]; #theta_{p} [rad]", 250, 0, 1000, 200, 2.8, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_neutron_weighted_cut",  new TH2D("dump_plane_bkg_track_r_track_E_neutron_weighted_cut","dump_plane_bkg_track_r_track_E_neutron_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_phi_pos_phi_neutron_cut",  new TH2D("dump_plane_bkg_track_phi_pos_phi_neutron_cut","dump_plane_bkg_track_phi_pos_phi_neutron_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_neutron_cut",  new TH2D("dump_plane_bkg_track_r_track_E_neutron_cut","dump_plane_bkg_track_r_track_E_neutron_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_E_neutron_cut",  new TH2D("dump_plane_bkg_time_track_E_neutron_cut","dump_plane_bkg_time_track_E_neutron_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_neutron_cut",  new TH2D("dump_plane_bkg_time_track_r_neutron_cut","dump_plane_bkg_time_track_r_neutron_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_neutron_cut",  new TH2D("dump_plane_bkg_time_track_theta_neutron_cut","dump_plane_bkg_time_track_theta_neutron_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_neutron_weighted_cut",  new TH2D("dump_plane_bkg_time_track_r_neutron_weighted_cut","dump_plane_bkg_time_track_r_neutron_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_neutron_weighted_cut",  new TH2D("dump_plane_bkg_time_track_theta_neutron_weighted_cut","dump_plane_bkg_time_track_theta_neutron_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_x_track_y_photon_cut",  new TH2D("dump_plane_bkg_track_x_track_y_photon_cut","dump_plane_bkg_track_x_track_y_photon_cut; track x [mm]; track y [mm]",nXBins, xyarray, nXBins, xyarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_theta_photon_weighted_cut",  new TH2D("dump_plane_bkg_track_r_track_theta_photon_weighted_cut","dump_plane_bkg_track_r_track_theta_photon_weighted_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_photon_cut",  new TH2D("dump_plane_bkg_track_rUp_track_theta_photon_cut","dump_plane_bkg_track_rUp_track_theta_photon_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_photon_cut",  new TH2D("dump_plane_bkg_track_rDn_track_theta_photon_cut","dump_plane_bkg_track_rDn_track_theta_photon_cut; r [mm]; #theta_{p} [rad]", nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_backward_photon_cut",  new TH2D("dump_plane_bkg_track_rUp_track_theta_backward_photon_cut","dump_plane_bkg_track_rUp_track_theta_backward_photon_cut; r [mm]; #theta_{p} [rad]", 250, 0, 1000, 200, 2.8, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_backward_photon_cut",  new TH2D("dump_plane_bkg_track_rDn_track_theta_backward_photon_cut","dump_plane_bkg_track_rDn_track_theta_backward_photon_cut; r [mm]; #theta_{p} [rad]", 250, 0, 1000, 200, 2.8, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_photon_weighted_cut",  new TH2D("dump_plane_bkg_track_r_track_E_photon_weighted_cut","dump_plane_bkg_track_r_track_E_photon_weighted_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_phi_pos_phi_photon_cut",  new TH2D("dump_plane_bkg_track_phi_pos_phi_photon_cut","dump_plane_bkg_track_phi_pos_phi_photon_cut; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_r_track_E_photon_cut",  new TH2D("dump_plane_bkg_track_r_track_E_photon_cut","dump_plane_bkg_track_r_track_E_photon_cut; r [mm]; E [GeV]", nRBins, rarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_E_photon_cut",  new TH2D("dump_plane_bkg_time_track_E_photon_cut","dump_plane_bkg_time_track_E_photon_cut; time [ns]; E [GeV]", nbins, tarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_photon_cut",  new TH2D("dump_plane_bkg_time_track_r_photon_cut","dump_plane_bkg_time_track_r_photon_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_photon_cut",  new TH2D("dump_plane_bkg_time_track_theta_photon_cut","dump_plane_bkg_time_track_theta_photon_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_r_photon_weighted_cut",  new TH2D("dump_plane_bkg_time_track_r_photon_weighted_cut","dump_plane_bkg_time_track_r_photon_weighted_cut; time [ns]; r [mm]", nbins, tarray, 1000, 0.0, 2000.0)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_time_track_theta_photon_weighted_cut",  new TH2D("dump_plane_bkg_time_track_theta_photon_weighted_cut","dump_plane_bkg_time_track_theta_photon_weighted_cut; time [ns]; #theta_{p} [rad]", nbins, tarray, 1600, 1.6, 3.2)));
    
    ////// 1D plot;
    
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_x_neutron_cut",  new TH1D("dump_plane_bkg_track_x_neutron_cut","dump_plane_bkg_track_x_neutron_cut; track x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_y_neutron_cut",  new TH1D("dump_plane_bkg_track_y_neutron_cut","dump_plane_bkg_track_y_neutron_cut; track y [mm]; Events",nXBins, xyarray)));
    
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_neutron_cut",  new TH1D("dump_plane_bkg_track_time_neutron_cut","dump_plane_bkg_track_time_neutron_cut; t [ns]; Events",nbins, tarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_energy_neutron_cut",  new TH1D("dump_plane_bkg_track_energy_neutron_cut","dump_plane_bkg_track_energy_neutron_cut; E [GeV]; Events",nbins, xarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rUp_neutron_cut",  new TH1D("dump_plane_bkg_track_rUp_neutron_cut","dump_plane_bkg_track_rUp_neutron_cut; r [mm]; Events",nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rDn_neutron_cut",  new TH1D("dump_plane_bkg_track_rDn_neutron_cut","dump_plane_bkg_track_rDn_neutron_cut; r [mm]; Events",nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xUp_neutron_cut",  new TH1D("dump_plane_bkg_track_xUp_neutron_cut","dump_plane_bkg_track_xUp_neutron_cut; x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xDn_neutron_cut",  new TH1D("dump_plane_bkg_track_xDn_neutron_cut","dump_plane_bkg_track_xDn_neutron_cut; x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_neutron_cut",  new TH1D("dump_plane_bkg_track_theta_neutron_cut","dump_plane_bkg_track_theta_neutron_cut; #theta_{p} [rad]; Events",800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaUp_neutron_weighted_cut",  new TH1D("dump_plane_bkg_track_thetaUp_neutron_weighted_cut","dump_plane_bkg_track_thetaUp_neutron_weighted_cut; #theta_{p} [rad]; Events",800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaDn_neutron_weighted_cut",  new TH1D("dump_plane_bkg_track_thetaDn_neutron_weighted_cut","dump_plane_bkg_track_thetaDn_neutron_weighted_cut; #theta_{p} [rad]; Events",800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_neutron_cut",  new TH1D("dump_plane_bkg_track_phi_neutron_cut","dump_plane_bkg_track_phi_neutron_cut; #phi_{p} [rad]; Events",6400, -3.2, 3.2)));

    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_x_photon_cut",  new TH1D("dump_plane_bkg_track_x_photon_cut","dump_plane_bkg_track_x_photon_cut; track x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_y_photon_cut",  new TH1D("dump_plane_bkg_track_y_photon_cut","dump_plane_bkg_track_y_photon_cut; track y [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_photon_cut",  new TH1D("dump_plane_bkg_track_time_photon_cut","dump_plane_bkg_track_time_photon_cut; t [ns]; Events",nbins, tarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_energy_photon_cut",  new TH1D("dump_plane_bkg_track_energy_photon_cut","dump_plane_bkg_track_energy_photon_cut; E [GeV]; Events",nbins, xarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rUp_photon_cut",  new TH1D("dump_plane_bkg_track_rUp_photon_cut","dump_plane_bkg_track_rUp_photon_cut; r [mm]; Events",nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rDn_photon_cut",  new TH1D("dump_plane_bkg_track_rDn_photon_cut","dump_plane_bkg_track_rDn_photon_cut; r [mm]; Events",nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xUp_photon_cut",  new TH1D("dump_plane_bkg_track_xUp_photon_cut","dump_plane_bkg_track_xUp_photon_cut; x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xDn_photon_cut",  new TH1D("dump_plane_bkg_track_xDn_photon_cut","dump_plane_bkg_track_xDn_photon_cut; x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaUp_photon_weighted_cut",  new TH1D("dump_plane_bkg_track_thetaUp_photon_weighted_cut","dump_plane_bkg_track_thetaUp_photon_weighted_cut; #theta_{p} [rad]; Events",800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaDn_photon_weighted_cut",  new TH1D("dump_plane_bkg_track_thetaDn_photon_weighted_cut","dump_plane_bkg_track_thetaDn_photon_weighted_cut; #theta_{p} [rad]; Events",800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_photon_cut",  new TH1D("dump_plane_bkg_track_theta_photon_cut","dump_plane_bkg_track_theta_photon_cut; #theta_{p} [rad]; Events",800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_photon_cut",  new TH1D("dump_plane_bkg_track_phi_photon_cut","dump_plane_bkg_track_phi_photon_cut; #phi_{p} [rad]; Events",6400, -3.2, 3.2)));
    
    long lineCounter   = 0;
    
    
    int bxNumber(-99999), pdgId(-99999), trackId(-99999), detId(-99999);
    double xPos(-99999.0), yPos(-99999.0), energyVal(-99999.0), weight(-99999.0);
    double vtx_x(-99999.0), vtx_y(-99999.0), vtx_z(-99999.0);
    int parent_id(-99999), physprocess(-99999);
    double pxx(-99999.0), pyy(-99999.0), pzz(-99999.0), time(-99999.0);
    
    int nntrn(0), npho(0);
    int niteration(1); /// the number of times we want to fill the high theta bins
    while(inFile >> bxNumber >> pdgId >> trackId >> detId >> xPos >> yPos >> energyVal >> weight >> vtx_x >> vtx_y >> vtx_z >> parent_id >> pxx >> pyy >> pzz >> physprocess >> time){
        
        
        lineCounter += 1;
        if(lineCounter%100000==0)
            cout << "processed: " << lineCounter << endl;

        // if(lineCounter > 1000000)
        //     break;
                
        
        if(pzz>0)
            continue;
            
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
        
        
        
        ////// if(time > 1000): continue;
        ////// neutrons;
        if(pdgId == 2112){
            // if(energyVal<1e-7)continue;
            nntrn++;
            /// limit the number of particles in case of comaprison with small stat
            bool needNeutron = false;
            if(cmpPlt){
               needNeutron = (nntrn < nntrnLim);
            }
            else
               needNeutron = true;

            if(needNeutron){
                allHisto2Dict["dump_plane_bkg_track_x_track_y_neutron_cut"]->Fill(xPos, yPos, weight);
                allHisto1Dict["dump_plane_bkg_track_x_neutron_cut"]->Fill(xPos, weight);
                allHisto1Dict["dump_plane_bkg_track_y_neutron_cut"]->Fill(yPos, weight);
                // allHistoDict["dump_plane_bkg_track_theta_track_phi_neutron_cut"]->Fill(theta, phi, weight);
                allHisto1Dict["dump_plane_bkg_track_theta_neutron_cut"]->Fill(theta, weight);
                
                allHisto1Dict["dump_plane_bkg_track_phi_neutron_cut"]->Fill(phi, weight);
                allHisto1Dict["dump_plane_bkg_track_time_neutron_cut"]->Fill(time, weight);
                allHisto1Dict["dump_plane_bkg_track_energy_neutron_cut"]->Fill(energyVal, weight);
                if(vtx_x > x0){
                    allHisto1Dict["dump_plane_bkg_track_rUp_neutron_cut"]->Fill(rValue, rWeight*weight);
                    allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_neutron_cut"]->Fill(rValue, theta, weight);
                    allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_backward_neutron_cut"]->Fill(rValue, theta, weight);
                    allHisto1Dict["dump_plane_bkg_track_thetaUp_neutron_weighted_cut"]->Fill(theta, thetaWeight*weight);
                    allHisto1Dict["dump_plane_bkg_track_xUp_neutron_cut"]->Fill(xPos, weight);
                }
                else{
                    allHisto1Dict["dump_plane_bkg_track_rDn_neutron_cut"]->Fill(rValue, rWeight*weight);
                    allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_neutron_cut"]->Fill(rValue, theta, weight);
                    allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_backward_neutron_cut"]->Fill(rValue, theta, weight);
                    allHisto1Dict["dump_plane_bkg_track_thetaDn_neutron_weighted_cut"]->Fill(theta, thetaWeight*weight);
                    allHisto1Dict["dump_plane_bkg_track_xDn_neutron_cut"]->Fill(xPos, weight);
                }
                allHisto2Dict["dump_plane_bkg_track_r_track_theta_neutron_weighted_cut"]->Fill(rValue, theta, rWeight*thetaWeight*weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_E_neutron_cut"]->Fill(rValue, energyVal, weight);
                allHisto2Dict["dump_plane_bkg_track_r_track_E_neutron_weighted_cut"]->Fill(rValue, energyVal, rWeight*weight);
                allHisto2Dict["dump_plane_bkg_time_track_E_neutron_cut"]->Fill(time,energyVal, weight);
                allHisto2Dict["dump_plane_bkg_time_track_r_neutron_cut"]->Fill(time,rValue, weight);
                allHisto2Dict["dump_plane_bkg_time_track_theta_neutron_cut"]->Fill(time,theta, weight);
                allHisto2Dict["dump_plane_bkg_time_track_r_neutron_weighted_cut"]->Fill(time,rValue, rWeight*weight);
                allHisto2Dict["dump_plane_bkg_time_track_theta_neutron_weighted_cut"]->Fill(time,theta, thetaWeight*weight);
            }
        }
        ////// photon;
        if(pdgId == 22){
            npho++;
            /// limit the number of particles in case of comaprison with small stat
            bool needPhoton = false;
            if(cmpPlt){
               needPhoton = (npho < nphoLim);
            }
            else
               needPhoton = true;

            if(needPhoton){
                allHisto2Dict["dump_plane_bkg_track_x_track_y_photon_cut"]->Fill(xPos, yPos, weight);
                allHisto1Dict["dump_plane_bkg_track_x_photon_cut"]->Fill(xPos, weight);
                allHisto1Dict["dump_plane_bkg_track_y_photon_cut"]->Fill(yPos, weight);
                // allHistoDict["dump_plane_bkg_track_theta_track_phi_photon_cut"]->Fill(theta, phi, weight);
                allHisto1Dict["dump_plane_bkg_track_theta_photon_cut"]->Fill(theta, weight);
                
                allHisto1Dict["dump_plane_bkg_track_phi_photon_cut"]->Fill(phi, weight);
                allHisto1Dict["dump_plane_bkg_track_time_photon_cut"]->Fill(time, weight);
                allHisto1Dict["dump_plane_bkg_track_energy_photon_cut"]->Fill(energyVal, weight);
                
                if(vtx_x > x0){
                    allHisto1Dict["dump_plane_bkg_track_rUp_photon_cut"]->Fill(rValue, rWeight*weight);
                    allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_photon_cut"]->Fill(rValue, theta, weight);
                    allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_backward_photon_cut"]->Fill(rValue, theta, weight);
                    allHisto1Dict["dump_plane_bkg_track_thetaUp_photon_weighted_cut"]->Fill(theta, thetaWeight*weight);
                    allHisto1Dict["dump_plane_bkg_track_xUp_photon_cut"]->Fill(xPos, weight);
                }
                else{
                    allHisto1Dict["dump_plane_bkg_track_rDn_photon_cut"]->Fill(rValue, rWeight*weight);
                    allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_photon_cut"]->Fill(rValue, theta, weight);
                    allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_backward_photon_cut"]->Fill(rValue, theta, weight);
                    allHisto1Dict["dump_plane_bkg_track_thetaDn_photon_weighted_cut"]->Fill(theta, thetaWeight*weight);
                    allHisto1Dict["dump_plane_bkg_track_xDn_photon_cut"]->Fill(xPos, weight);
                }

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
        if(cmpPlt && (nntrn > nntrnLim) && (npho > nphoLim)) break;
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
    std::cout << "total neutron: " << nntrn << std::endl;
    std::cout << "total photon: " << npho << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "time taken for the running: " << duration.count() << " s" << endl;
    
    
}
