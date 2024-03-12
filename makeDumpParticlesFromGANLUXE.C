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
#include <sstream>
#include <map>
#include <string>
#include <chrono>
#include <TRandom.h>


using namespace std;
using namespace std::chrono;

/// see the answer in https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
class CSVRow{
    public:
        std::string_view operator[](std::size_t index) const
        {
            return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data){
    data.readNextRow(str);
    return str;
}   


class CSVIterator{   
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:nullptr) { ++(*this); }
        CSVIterator()                   :m_str(nullptr) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = nullptr;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == nullptr) && (rhs.m_str == nullptr)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};

class CSVRange{
    std::istream&   stream;
    public:
        CSVRange(std::istream& str)
            : stream(str)
        {}
        CSVIterator begin() const {return CSVIterator{stream};}
        CSVIterator end()   const {return CSVIterator{};}
};

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

void makeDumpParticlesFromGANLUXE(string bkgFileName="/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn.root", int nntrn=5497735, int nphot=1353100, int bx=1, string part="v1"){
    
    //  int nntrn=5497735, int nphot=1353100
    auto start = high_resolution_clock::now();
    ////// the position which determines r_up or r_down, the distribution we need to take
    const double x0 = -92.65;

    
    long n_neutron   = nntrn;//long(random->gauss(nntrn, 1));   ////// randomly smearing the number of neutrons per BX
    long n_photon    = nphot;//long(random->gauss(nphot, 1));   ////// randomly smearing the number of photons per BX

    // cout << "n_neutron: " << n_neutron << " and n_photon: " << n_photon << endl;
    cout << "From the code: working on " << part << endl;

    // string inDir = "/Users/arkasantra/arka/Sasha_Work/OutputFile/";
    string inDir = "/srv01/agrp/arkas/GANFastSim/";

    /// get the 2D plots required for sampling
    std::string foutname    = bkgFileName.substr(bkgFileName.find_last_of("/")+1);
    foutname                = foutname.substr(0, foutname.find_last_of("."));
    std::string inrootname  = inDir+foutname+".root";

    fstream fin;
  
    // Open an existing file
    fin.open(bkgFileName, ios::in);
    string line;
    
    
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


    vector<string> row;
    int nevents = 0;
    int neutronNum = 0;
    int lineNum = 0;
    /// the csv file should have the following columns
    /// , rx, xx, yy, rp, phi_p, pzz, eneg, time, phi_x, pxx, pyy,theta
    /// new column from Alon
    /// ' xx', ' yy', ' rp', ' phi_p', ' pzz', ' eneg', ' time', ' phi_x', ' rx', ' pxx', ' pyy', 'theta'
    ///' xx', ' yy', ' rp', ' phi_p', ' pzz', ' eneg', ' time', ' phi_x', ' rx', ' pxx', ' pyy', 'theta'
    std::ifstream       file(bkgFileName);

    for(auto& row: CSVRange(file)){
        // if(lineNum%10000==0)cout << "Processed: " << lineNum << endl;
        // cout << "BX: " << nevents << endl;
        string row_str = string(row[1]);
        if (row_str.find(" rx")==0)continue;
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
        neutronNum++;
        if(neutronNum%10000==0) cout << "neutron: working on " << neutronNum << endl;
        if(neutronNum>n_neutron)break;
            
        // float m_neutron        = 0.939565; // GeV
        // double x_neutron       = stod(string(row[2]));
        // double y_neutron       = stod(string(row[3]));
        
        // double r_neutron   = stod(string(row[1]));
        // double phi_neutron = atan2(y_neutron, x_neutron);
        // double rWeight     = -999;
        // double thetaWeight = -999;
        // double thetap_neutron = stod(string(row[12]));
        // double phip_neutron   = stod(string(row[5]));
            
        // // p_neutron = math.sqrt(E_neutron**2 - m_neutron**2)
        // double time_neutron   = stod(string(row[8]));
        // double energy_neutron = stod(string(row[7]));
        
        
        // double p_neutron      = stod(string(row[4]));
        // double pxx            = stod(string(row[10]));
        // double pyy            = stod(string(row[11]));
        // double pzz            = stod(string(row[6]));
        // double vtxx_pos       = r_neutron*TMath::Cos(phi_neutron);
        // double vtxy_pos       = r_neutron*TMath::Sin(phi_neutron);

        /// new columns
        float m_neutron        = 0.939565; // GeV
        double x_neutron       = stod(string(row[1]));
        double y_neutron       = stod(string(row[2]));
        
        double r_neutron   = stod(string(row[9]));
        double phi_neutron = atan2(y_neutron, x_neutron);
        double rWeight     = -999;
        double thetaWeight = -999;
        double thetap_neutron = stod(string(row[12]));
        double phip_neutron   = stod(string(row[4]));
            
        // p_neutron = math.sqrt(E_neutron**2 - m_neutron**2)
        double time_neutron   = stod(string(row[7]));
        
        
        
        double p_neutron      = stod(string(row[3]));
        double pxx            = stod(string(row[10]));
        double pyy            = stod(string(row[11]));
        double pzz            = stod(string(row[5]));
        double energy_neutron = stod(string(row[6]));
        double vtxx_pos       = r_neutron*TMath::Cos(phi_neutron);
        double vtxy_pos       = r_neutron*TMath::Sin(phi_neutron);


            
        trackid.push_back(neutronNum);
        ptrackid.push_back(-10);
        detid.push_back(-10);
        pdg.push_back(2112);
        physproc.push_back(7000);
        E.push_back(energy_neutron);
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

        /// fill the tree
        Tracks->Fill();
        lineNum++;
    }
        
    Tracks->Write(); 
    outFile->Close();
    delete outFile;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "time taken for the running: " << duration.count() << " s" << endl;
    
}
