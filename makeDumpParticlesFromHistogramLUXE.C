//////// run: inside root, .L makeDumpParticlesFromHistogramLUXE.C++ && makeDumpParticlesFromHistogramLUXE()
// #### this code prepares the fastsampling particles from histogram
// #### this is using LUXE dump geometry
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <map>
#include <chrono>
#include <TRandom.h>


using namespace std;
using namespace std::chrono;


/// energy and time bins
const int nbins=650;
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

/// for no energy cut on neutron (all stat)
/// total neutron: 54977353
/// total photon: 13531005 
/// for energy cut > 1e-7 on neutron (small stat)
/// total neutron: 210410
/// total photon: 56555
/// for comparing plot 
/// nneutron = 549773;
/// nphoton = 135310;

void makeDumpParticlesFromHistogramLUXE(string bkgFileName="/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn.root", int nntrn=5497735, int nphot=1353100, int bx=1, string part="v1"){
// void makeDumpParticlesFromHistogramLUXE(string bkgFileName, int bx, int nphot, int nntrn, string part){
    
    
    auto start = high_resolution_clock::now();
    ////// the position which determines r_up or r_down, the distribution we need to take
    const double x0 = -92.65;

    
    long n_neutron   = nntrn;//long(random->gauss(nntrn, 1));   ////// randomly smearing the number of neutrons per BX
    long n_photon    = nphot;//long(random->gauss(nphot, 1));   ////// randomly smearing the number of photons per BX

    cout << "n_neutron: " << n_neutron << " and n_photon: " << n_photon << endl;
    cout << "From the code: working on " << part << endl;

    string inDir = "/Users/arkasantra/arka/Sasha_Work/OutputFile/";

    /// get the 2D plots required for sampling
    std::string foutname    = bkgFileName.substr(bkgFileName.find_last_of("/")+1);
    foutname                = foutname.substr(0, foutname.find_last_of("."));
    std::string inrootname  = inDir+foutname+".root";
    
    TFile *inFile           = new TFile(inrootname.c_str(), "READ");

    cout << "Reading " << inrootname.c_str() << endl;

    //////// load histograms from the inRoot file;
    TH2D* dump_plane_bkg_track_x_track_y_neutron_cut                = (TH2D*)inFile->Get("dump_plane_bkg_track_x_track_y_neutron_cut");
    TH2D* dump_plane_bkg_track_rUp_track_theta_neutron_cut          = (TH2D*)inFile->Get("dump_plane_bkg_track_rUp_track_theta_neutron_cut");
    TH2D* dump_plane_bkg_track_rDn_track_theta_neutron_cut          = (TH2D*)inFile->Get("dump_plane_bkg_track_rDn_track_theta_neutron_cut");
    TH1D* dump_plane_bkg_track_energy_neutron_cut                   = (TH1D*)inFile->Get("dump_plane_bkg_track_energy_neutron_cut");
    TH2D* dump_plane_bkg_track_phi_pos_phi_neutron_cut              = (TH2D*)inFile->Get("dump_plane_bkg_track_phi_pos_phi_neutron_cut");
    TH2D*  dump_plane_bkg_time_track_E_neutron_cut                  = (TH2D*)inFile->Get("dump_plane_bkg_time_track_E_neutron_cut");

    TH2D* dump_plane_bkg_track_x_track_y_photon_cut                 = (TH2D*)inFile->Get("dump_plane_bkg_track_x_track_y_photon_cut");
    TH2D* dump_plane_bkg_track_rUp_track_theta_photon_cut           = (TH2D*)inFile->Get("dump_plane_bkg_track_rUp_track_theta_photon_cut");
    TH2D* dump_plane_bkg_track_rDn_track_theta_photon_cut           = (TH2D*)inFile->Get("dump_plane_bkg_track_rDn_track_theta_photon_cut");
    TH1D* dump_plane_bkg_track_energy_photon_cut                    = (TH1D*)inFile->Get("dump_plane_bkg_track_energy_photon_cut");
    TH2D* dump_plane_bkg_track_phi_pos_phi_photon_cut               = (TH2D*)inFile->Get("dump_plane_bkg_track_phi_pos_phi_photon_cut");
    TH1D* dump_plane_bkg_track_time_photon_cut                      = (TH1D*)inFile->Get("dump_plane_bkg_track_time_photon_cut");



    //// first smoothen the 2D histograms before sampling
    cout << "Smothening the 2D histograms " << endl;
    // dump_plane_bkg_track_x_track_y_neutron_cut->Smooth();
    // dump_plane_bkg_track_rUp_track_theta_neutron_cut->Smooth();
    // dump_plane_bkg_track_rDn_track_theta_neutron_cut->Smooth();
    // dump_plane_bkg_track_phi_pos_phi_neutron_cut->Smooth();
    // dump_plane_bkg_time_track_E_neutron_cut->Smooth();

    for(int l=0; l < 2; l++){

        dump_plane_bkg_track_x_track_y_photon_cut->Smooth();
        dump_plane_bkg_track_rUp_track_theta_photon_cut->Smooth();
        dump_plane_bkg_track_rDn_track_theta_photon_cut->Smooth();
        dump_plane_bkg_track_phi_pos_phi_photon_cut->Smooth();
    }
    
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
    
//     for(int j=0; j < nbins+1; ++j){
//         cout << "xarray [ " << j << " ] " << xarray[j] << endl;
//         cout << "tarray [ " << j << " ] " << tarray[j] << endl;
//     }
//     
//     for(int j=0; j < nRBins+1; ++j){
//         cout << "rarray [ " << j << " ] " << rarray[j] << endl;
//     }
//     for(int j=0; j < nXBins+1; ++j){
//         cout << "xyarray [ " << j << " ] " << xyarray[j] << endl;
//     }
//     exit(-1);
    map<string, TH1D*> allHisto1Dict;
    map<string, TH2D*> allHisto2Dict;
    
    
    ////// output root file containing the photon and neutron information    
    
    std::string foutname2    = bkgFileName.substr(bkgFileName.find_last_of("/")+1);
    foutname2                = foutname2.substr(0, foutname2.find_last_of("."));
    std::string rootoutname2 = inDir+foutname2;
    rootoutname2             += std::string("_RandomGeneration_")+part+std::string(".root");
    
    TFile *outFile       = new TFile(rootoutname2.c_str(), "RECREATE");
    outFile->cd();
    TTree *Tracks = new TTree("Tracks", "Tracks");
    std::cout << "The output file is " << rootoutname2.c_str() << std::endl;


    int eventid;
    vector<int> trackid;
    trackid.push_back(30);
    vector<int> ptrackid;
    vector<int> detid;
    vector<int> pdg;
    vector<int> physproc;
    vector<double> E;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> t;
    vector<double> vtxx;
    vector<double> vtxy;
    vector<double> vtxz;
    vector<double> px;
    vector<double> py;
    vector<double> pz;
    vector<double> theta;
    vector<double> phi;
    vector<double> xlocal;
    vector<double> ylocal;
    vector<double> zlocal;
    double weight;
    vector<int> nsecondary;
    vector<double> esecondary;
    
    
    Tracks->Branch("eventid", &eventid, "eventid/I");
    Tracks->Branch("trackid", &trackid);
    Tracks->Branch("ptrackid", &ptrackid);
    Tracks->Branch("detid", &detid);
    Tracks->Branch("pdg", &pdg);
    Tracks->Branch("physproc", &physproc);
    Tracks->Branch("E", &E);
    Tracks->Branch("x", &x);
    Tracks->Branch("y", &y);
    Tracks->Branch("z", &z);
    Tracks->Branch("t", &t);
    Tracks->Branch("vtxx", &vtxx);
    Tracks->Branch("vtxy", &vtxy);
    Tracks->Branch("vtxz", &vtxz);
    Tracks->Branch("px", &px);
    Tracks->Branch("py", &py);
    Tracks->Branch("pz", &pz);
    Tracks->Branch("theta", &theta);
    Tracks->Branch("phi", &phi);
    Tracks->Branch("xlocal", &xlocal);
    Tracks->Branch("ylocal", &ylocal);
    Tracks->Branch("zlocal", &zlocal);
    Tracks->Branch("weight", &weight, "weight/D");
    Tracks->Branch("nsecondary", &nsecondary);
    Tracks->Branch("esecondary", &esecondary);


    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_x_neutron_cut", new TH1D("dump_plane_bkg_track_x_neutron_cut","dump_plane_bkg_track_x_neutron_cut; track x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_y_neutron_cut", new TH1D("dump_plane_bkg_track_y_neutron_cut","dump_plane_bkg_track_y_neutron_cut; track y [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_x_photon_cut", new TH1D("dump_plane_bkg_track_x_photon_cut","dump_plane_bkg_track_x_photon_cut; track x [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_y_photon_cut", new TH1D("dump_plane_bkg_track_y_photon_cut","dump_plane_bkg_track_y_photon_cut; track y [mm]; Events",nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_E_neutron", new TH1D("dump_plane_bkg_track_E_neutron","dump_plane_bkg_track_E_neutron; E [GeV]; Events", nbins, xarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_E_photon", new TH1D("dump_plane_bkg_track_E_photon","dump_plane_bkg_track_E_photon; E [GeV]; Events", nbins, xarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_neutron", new TH1D("dump_plane_bkg_track_time_neutron","dump_plane_bkg_track_time_neutron; t [ns]; Events", nbins, tarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_neutron_weighted", new TH1D("dump_plane_bkg_track_theta_neutron_weighted","dump_plane_bkg_track_theta_neutron_weighted; #theta_{p} [rad]; Events", 800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_theta_photon_weighted", new TH1D("dump_plane_bkg_track_theta_photon_weighted","dump_plane_bkg_track_theta_photon_weighted; #theta_{p} [rad]; Events", 800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rUp_neutron_weighted", new TH1D("dump_plane_bkg_track_rUp_neutron_weighted","dump_plane_bkg_track_rUp_neutron_weighted; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rUp_photon_weighted", new TH1D("dump_plane_bkg_track_rUp_photon_weighted","dump_plane_bkg_track_rUp_photon_weighted; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rDn_neutron_weighted", new TH1D("dump_plane_bkg_track_rDn_neutron_weighted","dump_plane_bkg_track_rDn_neutron_weighted; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rDn_photon_weighted", new TH1D("dump_plane_bkg_track_rDn_photon_weighted","dump_plane_bkg_track_rDn_photon_weighted; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xUp_neutron", new TH1D("dump_plane_bkg_track_xUp_neutron","dump_plane_bkg_track_xUp_neutron; x [mm]; Events", nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xUp_photon", new TH1D("dump_plane_bkg_track_xUp_photon","dump_plane_bkg_track_xUp_photon; x [mm]; Events", nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xDn_neutron", new TH1D("dump_plane_bkg_track_xDn_neutron","dump_plane_bkg_track_xDn_neutron; x [mm]; Events", nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_xDn_photon", new TH1D("dump_plane_bkg_track_xDn_photon","dump_plane_bkg_track_xDn_photon; x [mm]; Events", nXBins, xyarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaUp_neutron", new TH1D("dump_plane_bkg_track_thetaUp_neutron","dump_plane_bkg_track_thetaUp_neutron; #theta_{p} [rad]; Events", 800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaUp_neutron_morebins", new TH1D("dump_plane_bkg_track_thetaUp_neutron_morebins","dump_plane_bkg_track_thetaUp_neutron_morebins; #theta_{p} [rad]; Events", 1600, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaDn_neutron", new TH1D("dump_plane_bkg_track_thetaDn_neutron","dump_plane_bkg_track_thetaDn_neutron; #theta_{p} [rad]; Events", 800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaDn_neutron_morebins", new TH1D("dump_plane_bkg_track_thetaDn_neutron_morebins","dump_plane_bkg_track_thetaDn_neutron_morebins; #theta_{p} [rad]; Events", 1600, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaUp_photon", new TH1D("dump_plane_bkg_track_thetaUp_photon","dump_plane_bkg_track_thetaUp_photon; #theta_{p} [rad]; Events", 800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaUp_photon_morebins", new TH1D("dump_plane_bkg_track_thetaUp_photon_morebins","dump_plane_bkg_track_thetaUp_photon_morebins; #theta_{p} [rad]; Events", 1600, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaDn_photon", new TH1D("dump_plane_bkg_track_thetaDn_photon","dump_plane_bkg_track_thetaDn_photon; #theta_{p} [rad]; Events", 800, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_thetaDn_photon_morebins", new TH1D("dump_plane_bkg_track_thetaDn_photon_morebins","dump_plane_bkg_track_thetaDn_photon_morebins; #theta_{p} [rad]; Events", 1600, 1.6, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_neutron", new TH1D("dump_plane_bkg_track_phi_neutron","dump_plane_bkg_track_phi_neutron; #phi_{p} [rad]; Events", 6400, -3.2, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_phi_photon", new TH1D("dump_plane_bkg_track_phi_photon","dump_plane_bkg_track_phi_photon; #phi_{p} [rad]; Events", 6400, -3.2, 3.2)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rUp_neutron", new TH1D("dump_plane_bkg_track_rUp_neutron","dump_plane_bkg_track_rUp_neutron; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rUp_photon", new TH1D("dump_plane_bkg_track_rUp_photon","dump_plane_bkg_track_rUp_photon; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rDn_neutron", new TH1D("dump_plane_bkg_track_rDn_neutron","dump_plane_bkg_track_rDn_neutron; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_rDn_photon", new TH1D("dump_plane_bkg_track_rDn_photon","dump_plane_bkg_track_rDn_photon; r [mm]; Events", nRBins, rarray)));
    allHisto1Dict.insert(make_pair("dump_plane_bkg_track_time_photon", new TH1D("dump_plane_bkg_track_time_photon","dump_plane_bkg_track_time_photon; t [ns]; Events", nbins, tarray)));
    
    
    
    
    
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_E_time_neutron", new TH2D("dump_plane_bkg_track_E_time_neutron","dump_plane_bkg_track_E_time_neutron; t [ns]; E [GeV]", nbins, tarray, nbins, xarray)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_phi_pos_phi_neutron", new TH2D("dump_plane_bkg_track_phi_pos_phi_neutron","dump_plane_bkg_track_phi_pos_phi_neutron; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_phi_pos_phi_photon", new TH2D("dump_plane_bkg_track_phi_pos_phi_photon","dump_plane_bkg_track_phi_pos_phi_photon; #phi_{p} [rad]; #phi_{pos} [rad]", 6400, -3.2, 3.2, 6400, -3.2, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_neutron_weighted", new TH2D("dump_plane_bkg_track_rUp_track_theta_neutron_weighted","dump_plane_bkg_track_rUp_track_theta_neutron_weighted; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_photon_weighted", new TH2D("dump_plane_bkg_track_rUp_track_theta_photon_weighted","dump_plane_bkg_track_rUp_track_theta_photon_weighted; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_neutron", new TH2D("dump_plane_bkg_track_rUp_track_theta_neutron","dump_plane_bkg_track_rUp_track_theta_neutron; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rUp_track_theta_photon", new TH2D("dump_plane_bkg_track_rUp_track_theta_photon","dump_plane_bkg_track_rUp_track_theta_photon; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_neutron_weighted", new TH2D("dump_plane_bkg_track_rDn_track_theta_neutron_weighted","dump_plane_bkg_track_rDn_track_theta_neutron_weighted; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_photon_weighted", new TH2D("dump_plane_bkg_track_rDn_track_theta_photon_weighted","dump_plane_bkg_track_rDn_track_theta_photon_weighted; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_neutron", new TH2D("dump_plane_bkg_track_rDn_track_theta_neutron","dump_plane_bkg_track_rDn_track_theta_neutron; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    allHisto2Dict.insert(make_pair("dump_plane_bkg_track_rDn_track_theta_photon", new TH2D("dump_plane_bkg_track_rDn_track_theta_photon","dump_plane_bkg_track_rDn_track_theta_photon; r [mm]; #theta_{p} [rad]",nRBins, rarray, 800, 1.6, 3.2)));
    


    for(int nevents=0; nevents < bx; nevents++){
        cout << "BX: " << nevents << endl;
        ////// clearing the vector
        trackid.clear();
        trackid.push_back(20);
        ptrackid.clear();
        detid.clear();
        pdg.clear();
        physproc.clear();
        E.clear();
        x.clear();
        y.clear();
        z.clear();
        t.clear();
        vtxx.clear();
        vtxy.clear();
        vtxz.clear();
        px.clear();
        py.clear();
        pz.clear();
        theta.clear();
        phi.clear();
        xlocal.clear();
        ylocal.clear();
        zlocal.clear();
        nsecondary.clear();
        esecondary.clear();          
        
        
        //////// filling the branches
        eventid = nevents;
        weight  = 1.0;  
        for(int nneutron=0; nneutron < n_neutron; ++nneutron){
            if(nneutron%10000==0) cout << "neutron: working on " << nneutron << endl;
            // if(nneutron>10000)break;
            //////// work with neutrons
            //p_neutron       = dump_plane_bkg_track_energy_neutron_cut.GetRandom()
            float m_neutron       = 0.939565; // GeV
            //E_neutron       = p_neutron //math.sqrt(p_neutron**2+m_neutron**2) ////3 we only want kinematic energy
            //r_neutron       = ctypes.c_double(0)
            
            double x_neutron       = 0;
            double y_neutron       = 0;
            
            
            //dump_plane_bkg_track_r_track_theta_neutron_weighted_cut.GetRandom2(r_neutron, thetap_neutron) ////// from weighted distribution
            
            ////// get r from x and y plot
            dump_plane_bkg_track_x_track_y_neutron_cut->GetRandom2(x_neutron, y_neutron);
            
            allHisto1Dict["dump_plane_bkg_track_x_neutron_cut"]->Fill(x_neutron);
            allHisto1Dict["dump_plane_bkg_track_y_neutron_cut"]->Fill(y_neutron);
            double r_neutron   = sqrt(x_neutron*x_neutron+y_neutron*y_neutron);
            double phi_neutron = atan2(y_neutron, x_neutron);
            double rWeight     = -999;
            double thetaWeight = -999;
            double thetap_neutron = -999;
            double phip_neutron   = -999;


            ////// project from 2D plots
            if(x_neutron > x0){
                try{
                    int rBinValue =  dump_plane_bkg_track_rUp_track_theta_neutron_cut->GetXaxis()->FindBin(r_neutron);
                    TH1D *dump_plane_bkg_rUp_track_theta_neutron_cut_1D = dump_plane_bkg_track_rUp_track_theta_neutron_cut->ProjectionY("dump_plane_bkg_rUp_track_theta_neutron_cut_1D", rBinValue, rBinValue);
                    thetap_neutron = dump_plane_bkg_rUp_track_theta_neutron_cut_1D->GetRandom();
                    throw(rBinValue);
                }
                catch(int rBin){
                    //cout << "neutron problem at rUp: " << r_neutron << " and " << rBin << endl;
                    ;
                }

                /// now modify the theta_p according to the target distribution
                
                // if(thetap_neutron > 3.0){
                //     float randomThetaCrn = -999.0;
                //     // draw from 1D distribution, keep it only when it is > 3.0
                //     // while(randomThetaCrn < 3.0){
                //     //     randomThetaCrn = thetapUpNtrnTrgt->GetRandom();
                //     // }
                //     /// just draw from 1D distribution
                //     randomThetaCrn = thetapUpNtrnTrgt->GetRandom();
                //     // cout << "modified theta " << randomThetaCrn << endl;
                //     thetap_neutron = randomThetaCrn;
                // }

                if (r_neutron < 2.01)
                    rWeight      = 1./(2*TMath::Pi()*1.0);
                else
                    rWeight      = 1./(2*TMath::Pi()*r_neutron);
                    
                //////// sin theta should not be 0, otherwise the weight will blow up
                if (thetap_neutron < 0.01)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(0.005));
                else if (thetap_neutron > 3.13)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(3.135));
                else
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(thetap_neutron));
                
                allHisto1Dict["dump_plane_bkg_track_rUp_neutron_weighted"]->Fill(r_neutron, rWeight);
                allHisto1Dict["dump_plane_bkg_track_rUp_neutron"]->Fill(r_neutron);
                allHisto1Dict["dump_plane_bkg_track_xUp_neutron"]->Fill(x_neutron);
                allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_neutron"]->Fill(r_neutron, thetap_neutron);
                allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_neutron_weighted"]->Fill(r_neutron, thetap_neutron, rWeight*thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_thetaUp_neutron"]->Fill(thetap_neutron, thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_thetaUp_neutron_morebins"]->Fill(thetap_neutron, thetaWeight);
            }
            else{
                try{
                    int rBinValue =  dump_plane_bkg_track_rDn_track_theta_neutron_cut->FindBin(r_neutron);
                    TH1D* dump_plane_bkg_rDn_track_theta_neutron_cut_1D = dump_plane_bkg_track_rDn_track_theta_neutron_cut->ProjectionY("dump_plane_bkg_rDn_track_theta_neutron_cut_1D", rBinValue, rBinValue);
                    thetap_neutron = dump_plane_bkg_rDn_track_theta_neutron_cut_1D->GetRandom();
                    throw(rBinValue);
                }
                catch(int rBin){
                    // cout << "neutron problem at rDn: " << r_neutron << " and " << rBin << endl;
                    ;
                }

                /// now modify the theta_p according to the target distribution
                
                // if(thetap_neutron > 3.0){
                //     float randomThetaCrn = -999.0;
                //     // draw from 1D distribution, keep it only when it is > 3.0
                //     // while(randomThetaCrn < 3.0){
                //     //     randomThetaCrn = thetapDnNtrnTrgt->GetRandom();
                //     // }
                //     /// just draw from 1D distribution
                //     randomThetaCrn = thetapDnNtrnTrgt->GetRandom();
                //     // cout << "modified theta " << randomThetaCrn << endl;
                //     thetap_neutron = randomThetaCrn;
                // }

                if (r_neutron < 2.01)
                    rWeight      = 1./(2*TMath::Pi()*1.0);
                else
                    rWeight      = 1./(2*TMath::Pi()*r_neutron);
                    
                //////// sin theta should not be 0, otherwise the weight will blow up
                if (thetap_neutron < 0.01)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(0.005));
                else if (thetap_neutron > 3.13)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(3.135));
                else
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(thetap_neutron));
                    
                allHisto1Dict["dump_plane_bkg_track_rDn_neutron_weighted"]->Fill(r_neutron, rWeight);
                allHisto1Dict["dump_plane_bkg_track_rDn_neutron"]->Fill(r_neutron);
                allHisto1Dict["dump_plane_bkg_track_xDn_neutron"]->Fill(x_neutron);
                allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_neutron"]->Fill(r_neutron, thetap_neutron);
                allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_neutron_weighted"]->Fill(r_neutron, thetap_neutron, rWeight*thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_thetaDn_neutron"]->Fill(thetap_neutron, thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_thetaDn_neutron_morebins"]->Fill(thetap_neutron, thetaWeight);
            }
            //dump_plane_bkg_track_phi_pos_phi_neutron_cut.GetRandom2(phip_neutron, phi_neutron)
            try{
                int phiBin = dump_plane_bkg_track_phi_pos_phi_neutron_cut->GetYaxis()->FindBin(phi_neutron);
                TH1D *dump_plane_bkg_track_phi_neutron_cut_1D = dump_plane_bkg_track_phi_pos_phi_neutron_cut->ProjectionX("dump_plane_bkg_track_phi_neutron_cut_1D", phiBin, phiBin);
                phip_neutron = dump_plane_bkg_track_phi_neutron_cut_1D->GetRandom();
                throw(phiBin);
            }
            catch(int phBin){
                // cout << "neutron problem at phi: " << phi_neutron << " and " << phBin << endl; 
                ;
            }
            allHisto1Dict["dump_plane_bkg_track_phi_neutron"]->Fill(phip_neutron);
            
            // p_neutron = math.sqrt(E_neutron**2 - m_neutron**2)
            double time_neutron   = 0;
            double energy_neutron = 0;
            
            dump_plane_bkg_time_track_E_neutron_cut->GetRandom2(time_neutron, energy_neutron);

            double E_neutron      = energy_neutron;
            
            double p_neutron      = E_neutron;
            double pxx            = p_neutron*TMath::Sin(thetap_neutron)*TMath::Cos(phip_neutron);
            double pyy            = p_neutron*TMath::Sin(thetap_neutron)*TMath::Sin(phip_neutron);
            double pzz            = p_neutron*TMath::Cos(thetap_neutron);
            double vtxx_pos       = r_neutron*TMath::Cos(phi_neutron);
            double vtxy_pos       = r_neutron*TMath::Sin(phi_neutron);


            // vector<double> angles;
            // angles.clear();
            // angles                 = getTrackThetaPhi(pxx, pyy, pzz);
            // double theta           = angles.at(0);
            // double phi             = angles.at(1);
            
            // vector<double> anglesPos;
            // anglesPos.clear();
            // anglesPos              = getTrackThetaPhi(xPos, yPos, zPos); 
            // double thetaPos        = anglesPos.at(0);
            // double phiPos          = anglesPos.at(1); 
            // double rValue          = getR(xPos, yPos);
            // double rWeight         = 1./(2*TMath::Pi()*rValue);
            // double thetaWeight     = 1./(2*TMath::Pi()*TMath::Sin(theta));
            
            trackid.push_back(nneutron);
            ptrackid.push_back(-10);
            detid.push_back(-10);
            pdg.push_back(2112);
            physproc.push_back(7000);
            E.push_back(E_neutron);
            x.push_back(vtxx_pos);
            y.push_back(vtxy_pos);
            z.push_back(6621.91);
            t.push_back(time_neutron);
            vtxx.push_back(vtxx_pos);
            vtxy.push_back(vtxy_pos);
            vtxz.push_back(6621.91);
            px.push_back(pxx);
            py.push_back(pyy);
            pz.push_back(pzz);
            theta.push_back(thetap_neutron);
            phi.push_back(phip_neutron);
            xlocal.push_back(vtxx_pos);
            ylocal.push_back(vtxy_pos);
            zlocal.push_back(6621.91);
            nsecondary.push_back(0);
            esecondary.push_back(0.0);

            allHisto1Dict["dump_plane_bkg_track_E_neutron"]->Fill(E_neutron);
            allHisto1Dict["dump_plane_bkg_track_time_neutron"]->Fill(time_neutron);
            allHisto2Dict["dump_plane_bkg_track_E_time_neutron"]->Fill(time_neutron, E_neutron);
            allHisto1Dict["dump_plane_bkg_track_theta_neutron_weighted"]->Fill(thetap_neutron, thetaWeight);
            allHisto2Dict["dump_plane_bkg_track_phi_pos_phi_neutron"]->Fill(phip_neutron, phi_neutron);
        }

        for(long nphoton=0; nphoton < n_photon; ++nphoton){
            if(nphoton%10000==0) cout << "photon: working on " << nphoton << endl;
            // if(nphoton>10000)break;
            //////// work with photons
            double p_photon       = dump_plane_bkg_track_energy_photon_cut->GetRandom();
            double m_photon       = 0.0; // GeV
            
            double x_photon       = 0;
            double y_photon       = 0;
            
            ////// get r from x and y plot
            dump_plane_bkg_track_x_track_y_photon_cut->GetRandom2(x_photon, y_photon);
            allHisto1Dict["dump_plane_bkg_track_x_photon_cut"]->Fill(x_photon);
            allHisto1Dict["dump_plane_bkg_track_y_photon_cut"]->Fill(y_photon);
            double r_photon         = sqrt(x_photon*x_photon+y_photon*y_photon);
            double phi_photon       = atan2(y_photon, x_photon);
            double rWeight          = -999.0;
            double thetaWeight      = -999.0;
            double thetap_photon    = -999.0;
            double phip_photon      = -999.0;
            ////// project from 2D plots
            if(x_photon > x0){
                try{
                    int rBinValue =  dump_plane_bkg_track_rUp_track_theta_photon_cut->GetXaxis()->FindBin(r_photon);
                    TH1D *dump_plane_bkg_rUp_track_theta_photon_cut_1D = dump_plane_bkg_track_rUp_track_theta_photon_cut->ProjectionY("dump_plane_bkg_rUp_track_theta_photon_cut_1D", rBinValue, rBinValue);
                    thetap_photon = dump_plane_bkg_rUp_track_theta_photon_cut_1D->GetRandom();
                    throw(rBinValue);
                }
                catch(int rBin){
                    // cout << "photon problem at rUp: " << r_photon << " and " << rBin << endl;
                    ;
                }
                /// now modify the theta_p according to the target distribution
                
                // if(thetap_photon > 3.0){
                //     float randomThetaCrn = -999.0;
                //     // draw from 1D distribution, keep it only when it is > 3.0
                //     // while(randomThetaCrn < 3.0){
                //     //     randomThetaCrn = thetapUpPhotTrgt->GetRandom();
                //     // }
                //     /// just draw from 1D distribution
                //     randomThetaCrn = thetapUpPhotTrgt->GetRandom();
                //     // cout << "modified theta in photon " << randomThetaCrn << endl;
                //     thetap_photon = randomThetaCrn;
                // }

                if (r_photon < 2.01)
                    rWeight      = 1./(2*TMath::Pi()*1.0);
                else
                    rWeight      = 1./(2*TMath::Pi()*r_photon);
        
                //////// sin theta should not be 0, otherwise the weight will blow up
                
                if (thetap_photon < 0.01)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(0.005));
                else if (thetap_photon > 3.13)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(3.135));
                else
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(thetap_photon));
                
                allHisto1Dict["dump_plane_bkg_track_rUp_photon_weighted"]->Fill(r_photon, rWeight);
                allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_photon_weighted"]->Fill(r_photon, thetap_photon, rWeight*thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_rUp_photon"]->Fill(r_photon);
                allHisto1Dict["dump_plane_bkg_track_xUp_photon"]->Fill(x_photon);
                allHisto2Dict["dump_plane_bkg_track_rUp_track_theta_photon"]->Fill(r_photon, thetap_photon);
                allHisto1Dict["dump_plane_bkg_track_thetaUp_photon"]->Fill(thetap_photon, thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_thetaUp_photon_morebins"]->Fill(thetap_photon, thetaWeight);
            }

            else{
                try{
                    int rBinValue =  dump_plane_bkg_track_rDn_track_theta_photon_cut->FindBin(r_photon);
                    TH1D* dump_plane_bkg_rDn_track_theta_photon_cut_1D = dump_plane_bkg_track_rDn_track_theta_photon_cut->ProjectionY("dump_plane_bkg_rDn_track_theta_photon_cut_1D", rBinValue, rBinValue);
                    thetap_photon = dump_plane_bkg_rDn_track_theta_photon_cut_1D->GetRandom();
                    throw(rBinValue);
                }
                catch(int rBin){
                    // cout << "photon problem at rDn: " << r_photon << " and " << rBin << endl;
                    ;
                }
                
                // if(thetap_photon > 3.0){
                //     float randomThetaCrn = -999.0;
                //     // draw from 1D distribution, keep it only when it is > 3.0
                //     // while(randomThetaCrn < 3.0){
                //     //     randomThetaCrn = thetapDnPhotTrgt->GetRandom();
                //     // }
                //     /// just draw from 1D distribution
                //     randomThetaCrn = thetapDnPhotTrgt->GetRandom();
                //     // cout << "modified theta in photon " << randomThetaCrn << endl;
                //     thetap_photon = randomThetaCrn;
                // }

                if (r_photon < 2.01)
                    rWeight      = 1./(2*TMath::Pi()*1.0);
                else
                    rWeight      = 1./(2*TMath::Pi()*r_photon);
        
                //////// sin theta should not be 0, otherwise the weight will blow up
                if (thetap_photon < 0.01)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(0.005));
                else if (thetap_photon > 3.13)
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(3.135));
                else
                    thetaWeight  = 1./(2*TMath::Pi()*TMath::Sin(thetap_photon));
            
                allHisto1Dict["dump_plane_bkg_track_rDn_photon_weighted"]->Fill(r_photon, rWeight);
                allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_photon_weighted"]->Fill(r_photon, thetap_photon, rWeight*thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_rDn_photon"]->Fill(r_photon);
                allHisto1Dict["dump_plane_bkg_track_xDn_photon"]->Fill(x_photon);
                allHisto2Dict["dump_plane_bkg_track_rDn_track_theta_photon"]->Fill(r_photon, thetap_photon);
                allHisto1Dict["dump_plane_bkg_track_thetaDn_photon"]->Fill(thetap_photon, thetaWeight);
                allHisto1Dict["dump_plane_bkg_track_thetaDn_photon_morebins"]->Fill(thetap_photon, thetaWeight);
            }
            //dump_plane_bkg_track_phi_pos_phi_photon_cut.GetRandom2(phip_photon, phi_photon) 
            try{
                int phiBin = dump_plane_bkg_track_phi_pos_phi_photon_cut->GetYaxis()->FindBin(phi_photon);
                TH1D* dump_plane_bkg_track_phi_photon_cut_1D = dump_plane_bkg_track_phi_pos_phi_photon_cut->ProjectionX("dump_plane_bkg_track_phi_photon_cut_1D", phiBin, phiBin);
                phip_photon = dump_plane_bkg_track_phi_photon_cut_1D->GetRandom();
                throw(phiBin);
            }
            catch(int phBin){
                // cout << "photon problem at phi: " << phi_photon << " and " << phBin << endl;
                ;
            }
            
            allHisto1Dict["dump_plane_bkg_track_phi_photon"]->Fill(phip_photon);
            
            
            //dump_plane_bkg_track_r_track_theta_photon_weighted_cut.GetRandom2(r_photon, thetap_photon)
            double E_photon       = sqrt(p_photon*p_photon + m_photon*m_photon);
            double pxx            = p_photon*TMath::Sin(thetap_photon)*TMath::Cos(phip_photon);
            double pyy            = p_photon*TMath::Sin(thetap_photon)*TMath::Sin(phip_photon);
            double pzz            = p_photon*TMath::Cos(thetap_photon);
            double vtxx_pos       = r_photon*TMath::Cos(phi_photon);
            double vtxy_pos       = r_photon*TMath::Sin(phi_photon);

            trackid.push_back(nphoton);
            ptrackid.push_back(-10);
            detid.push_back(-10);
            pdg.push_back(22);
            physproc.push_back(7000);
            double time_photon = dump_plane_bkg_track_time_photon_cut->GetRandom();
            E.push_back(E_photon);
            x.push_back(vtxx_pos);
            y.push_back(vtxy_pos);
            z.push_back(6621.91);
            t.push_back(time_photon);  ////// uniform timing upto 1000 ns
            vtxx.push_back(vtxx_pos);
            vtxy.push_back(vtxy_pos);
            vtxz.push_back(6621.91);
            px.push_back(pxx);
            py.push_back(pyy);
            pz.push_back(pzz);
            theta.push_back(thetap_photon);
            phi.push_back(phip_photon);
            xlocal.push_back(vtxx_pos);
            ylocal.push_back(vtxy_pos);
            zlocal.push_back(6621.91);
            nsecondary.push_back(0);
            esecondary.push_back(0.0);


            allHisto1Dict["dump_plane_bkg_track_E_photon"]->Fill(p_photon);
            allHisto1Dict["dump_plane_bkg_track_time_photon"]->Fill(time_photon);
            allHisto1Dict["dump_plane_bkg_track_theta_photon_weighted"]->Fill(thetap_photon, thetaWeight);
            allHisto2Dict["dump_plane_bkg_track_phi_pos_phi_photon"]->Fill(phip_photon, phi_photon);
        }   
        
        Tracks->Fill();
    }
        
    Tracks->Write();
    /// write the histograms
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
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "time taken for the running: " << duration.count() << " s" << endl;
    
}
