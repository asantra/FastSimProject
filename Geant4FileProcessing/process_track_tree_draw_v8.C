//
// This is for electron+laser samples, only print the numbers in a text file, no root file, input is a text file which is a list of text files containing root files, distinction in bunch crossing
// this can be used for dump files as well
// Used for fastsim project



#ifndef __RUN_PROC_TREE__

void process_track_tree_draw_v8(const char *fnlist = 0,  std::string bxNumber=0, std::string inFileSuffix=0)
{
   gROOT->ProcessLineSync("#define __RUN_PROC_TREE__ 1");
   gROOT->ProcessLineSync(".L MHists.C+");
   gROOT->ProcessLine("#include \"process_track_tree_draw_v8.C\"");
   gROOT->ProcessLine(Form("run_process_track_tree_draw(\"%s\", \"%s\", \"%s\")", fnlist, bxNumber.c_str(), inFileSuffix.c_str()));
   gROOT->ProcessLine("#undef __RUN_PROC_TREE__");
}

#else


int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);


int run_process_track_tree_draw(const char *fnlist = 0,  std::string bxNumber=0, std::string inFileSuffix=0)
{
  int debugl = 0; //1;

  if (!fnlist) {
    std::cout << "Usage: root -l process_track_tree_draw_v8.C\'(\"file_with_list_of_mc_files\")\'\n";
    return -1;
  }
  
  
  
  std::string fnamelist(fnlist);
  std::vector<std::string>  flist;
  ProcessList(fnamelist, flist);  
  if (debugl) {
    std::cout << "The following files will be processed:\n";
    std::for_each(flist.begin(), flist.end(), [](const std::string ss) {std::cout << ss << std::endl;});
  }
  
  /// put your outDir name
  std::string outDir = "/nfs/dust/luxe/user/santraar/TextSamples_October2021/DumpFromFastSim/FastSimFiles_LUXEDump/";
  std::string suffix("");
  std::string foutname = inFileSuffix.substr(inFileSuffix.find_last_of("/")+1);
  foutname = foutname.substr(0, foutname.find_last_of("."));
  std::string textoutname = outDir+foutname;
  
  /// write to a text file 
  std::cout << "output text name: " << textoutname << std::endl;
  ofstream trackFile;
  textoutname += std::string("_trackInfo") + std::string(".txt");
  trackFile.open(textoutname, fstream::in | fstream::out | fstream::app);

  if(debugl)std::cout << "After creating the histograms" << std::endl; 
  
  TChain *track_tree = new TChain("Tracks");
  std::for_each(flist.begin(), flist.end(), [track_tree](const std::string ss) {track_tree->Add(ss.c_str());} );

  double ev_weight;
  int ev_id;

  std::vector<int> det_idv, pdgv, track_idv, n_secondv, ptidv, physprocv;
  std::vector<int> *det_id_ptr = &det_idv;
  std::vector<int> *pdg_ptr = &pdgv;
  std::vector<int> *track_id_ptr = &track_idv;
  std::vector<int> *n_second_ptr = &n_secondv;
  std::vector<int> *ptid_ptr = &ptidv;
  std::vector<int> *physproc_ptr = &physprocv;

  std::vector<double> enegv, xxv, yyv, zzv, ttv, pxxv, pyyv, pzzv, vtx_xv, vtx_yv, vtx_zv;
  std::vector<double> *xx_ptr = &xxv;
  std::vector<double> *yy_ptr = &yyv;
  std::vector<double> *zz_ptr = &zzv;
  std::vector<double> *tt_ptr = &ttv;
  std::vector<double> *eneg_ptr = &enegv;
  std::vector<double> *pxx_ptr = &pxxv;
  std::vector<double> *pyy_ptr = &pyyv;
  std::vector<double> *pzz_ptr = &pzzv;
  std::vector<double> *vtx_x_ptr = &vtx_xv;
  std::vector<double> *vtx_y_ptr = &vtx_yv;
  std::vector<double> *vtx_z_ptr = &vtx_zv;
  
  track_tree->SetBranchAddress("eventid", &ev_id);
  track_tree->SetBranchAddress("weight", &ev_weight);
  track_tree->SetBranchAddress("detid", &det_id_ptr);
  track_tree->SetBranchAddress("pdg", &pdg_ptr);
  track_tree->SetBranchAddress("trackid", &track_id_ptr);
  track_tree->SetBranchAddress("x", &xx_ptr);
  track_tree->SetBranchAddress("y", &yy_ptr);
  track_tree->SetBranchAddress("z", &zz_ptr);
  track_tree->SetBranchAddress("t", &tt_ptr);
  track_tree->SetBranchAddress("E",  &eneg_ptr);
  track_tree->SetBranchAddress("px", &pxx_ptr);
  track_tree->SetBranchAddress("py", &pyy_ptr);
  track_tree->SetBranchAddress("pz", &pzz_ptr);
  track_tree->SetBranchAddress("vtxx", &vtx_x_ptr);
  track_tree->SetBranchAddress("vtxy", &vtx_y_ptr);
  track_tree->SetBranchAddress("vtxz", &vtx_z_ptr);
  track_tree->SetBranchAddress("physproc", &physproc_ptr);
  track_tree->SetBranchAddress("ptrackid", &ptid_ptr);
  track_tree->SetBranchAddress("nsecondary", &n_second_ptr);


  int pevent = -1;
  int nevproc = track_tree->GetEntries();
  if(true)std::cout << "Total number of events: " << nevproc << std::endl;
  
  for (Long64_t ii = 0; ii < nevproc; ++ii) {
    track_tree->GetEntry(ii);
    if(debugl)std::cout << "Inside the loop. Processing " << ii << std::endl;
    if (!(ii % 100)) { std::cout << "Event " << ii << std::endl;}
    
    size_t ntrck = det_idv.size();
    int primary_pdg = -99;

    if(debugl)std::cout << "After getting primary particle energies" << std::endl;
    for (size_t jj = 0; jj < ntrck; ++jj){
        int det_id, pdg, track_id, parent_id;
        det_id = det_idv[jj];
        pdg = pdgv[jj];
        track_id = track_idv[jj];
        
        double eneg, xx, yy, zz, tt, pxx, pyy, pzz, vtx_x, vtx_y, vtx_z, physprocess;
        eneg = enegv[jj];
        xx = xxv[jj];
        yy = yyv[jj];
        zz = zzv[jj];
        tt = ttv[jj];
        pxx = pxxv[jj];
        pyy = pyyv[jj];
        pzz = pzzv[jj];
        vtx_x = vtx_xv[jj];
        vtx_y = vtx_yv[jj];
        vtx_z = vtx_zv[jj];
        physprocess = physprocv[jj];
        parent_id  = ptidv[jj];
        
        if(debugl)std::cout << "Getting information from the test surfaces and sampling plane" << std::endl;
        if ((det_id >=30 && det_id <= 33)) {
            trackFile << bxNumber << " " << pdg << " " << track_id << " " << det_id << " " << xx << " " << yy << " " << eneg << " " << ev_weight << " " << vtx_x << " " << vtx_y << " " << vtx_z << " " << parent_id << " " << pxx << " " << pyy << " " << pzz << " " << physprocess << " " << tt << std::endl;
        }
    } //tracks loop
  } // tree loop
  
  //closing the text file
  trackFile.close();
  
  delete track_tree;
  return 0;  
}






int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist)
{
  std::fstream  fdata;
  fdata.open(fnamelist, std::ios::in);
  if (!fdata.is_open()) {
    throw std::runtime_error(std::string("Error reding data from the file ") + fnamelist);
  }
  
  unsigned long lid = 0;
  while (!fdata.eof()) {
    std::string  ffname;
    double fweight;
    fdata >> ffname;
    if (!fdata.fail()) { 
      flist.push_back(ffname);
    }
    else if (fdata.eof()) { break; }
    else {
      std::cout << "ProcessList(..)  :  Error reading data from the file " << fnamelist 
                << ",  line: " << lid << ". Exit." << std::endl;
      fdata.close();          
      return -2;
    }
    ++lid;
  }
  
  fdata.close();

  return 0;
}  

#endif
  

