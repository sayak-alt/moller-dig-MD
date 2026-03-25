#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <chrono>

#include <TROOT.h>
#include <TString.h>
#include <TObjString.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>


#include "remoll_tree.h"
#include "remollDigAuxi.h"
#include "GEMSimDigitizer.h"


using namespace std;
int main(int argc, char **argv)
{
  string db_file, inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  //UShort_t Nbkgd = 0;//number of background files to add to each event
  double BkgdTimeWindow = 0, LumiFrac = 0;
  bool pmtbkgddig = false;
      
  if(argc<3 || argc>4){
    std::cout << "*** Inadequate number of arguments! ***" << std::endl
	 << " Arguments: database (mandatory); " << std::endl
	 << "           list_of_sig_input_files (str, mandatory); " << std::endl
	 << "          nb_of_sig_evts_to_process (int, def=-1); " << std::endl;
	 std::cerr << "Usgae: ./mollerdig db_file file_list.txt #Events"<<std::endl;
      //<< "         bkgd_histo_input_file (str, def=''); " << endl
      // << "        bkgd_lumi_frac (double, def=0); " << endl;
    return(-1);
  }
  
  db_file = argv[1];
  cout << " database file " << db_file << endl;
  inputsigfile = argv[2];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>3)Nentries = atoi(argv[3]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  /*
  if(argc>5){
    inputbkgdfile = argv[4];
    cout << " Background histgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[5]));
    cout << " Fraction of background to superimpose to signal = " << LumiFrac << endl;
   }*/
/*
   if (argc < 2) {
       std::cerr << "Usage: ./mollerdig file.root\n";
       return 1;
   }

   TString infile = argv[1];
   std::cout << "Input file = " << infile << std::endl;
*/
  std::vector<GEMDet*> GEMdetectors;
  std::vector<GEMSimDigitizer*> GEMSimDig;
  std::vector<int> gemdetmap;

  // Variable parameters. 
  // Can be configured with the database, but are provided with defaults.
  Int_t Rseed = 0;
  Double_t TriggerJitter = 3.0;
  
  std::vector<TString> detectors_list;
 //Number of parameters in the database 
  const int nparam_pmtdet_adc = 12;
  const int nparam_pmtdet_fadc = 11;
  const int nparam_gemdet = 12;

  int nparam_mollergem_read = 0;
  Int_t NPlanes_mollergem = 32;// number of planes/modules/readout
  Double_t gatewidth_mollergem = 400.;
  Double_t ZsupThr_mollergem = 240.;
  Int_t Nlayers_mollergem = 4;
  std::vector<Double_t> mollergem_layer_z;
  Int_t* layer_mollergem;
  Int_t* nstrips_mollergem;
  Double_t* offset_mollergem;
  Double_t* RO_angle_mollergem;
  Double_t* triggeroffset_mollergem;
  Double_t* gain_mollergem;//one gain per module
  Double_t* commonmode_array_mollergem;
  UShort_t nAPV_mollergem = 0;

  //-----------------------------
  //  Read database
  //-----------------------------
  cout << "read database: " << db_file.c_str() << endl;
  ifstream in_db(db_file.c_str());
  if(!in_db.is_open()){
    cout << "database " << db_file.c_str() << " does not exist!!!" << endl;
    exit(-1);
  }
  
  TString currentline;
  while( currentline.ReadLine(in_db) && !currentline.BeginsWith("endconfig")){
    if( !currentline.BeginsWith("#") ){
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
      }
      //TObjArray *tokens = currentline.Tokenize(" ");//vg: def lost => versions prior to 6.06; should be fixed! ??? 
      //int ntokens = tokens->GetEntries();
      
      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	
	if(skey=="Rseed"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Rseed = stemp.Atoi();
	}
	
	if(skey=="TriggerJitter"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TriggerJitter = stemp.Atof();
	}
	
	if(skey=="detectors_list"){
	  for(int k = 1; k<ntokens; k++){
	    TString sdet = ( (TObjString*) (*tokens)[k] )->GetString();
	    detectors_list.push_back(sdet);
	  }
	}

	//GEMs
	if(skey=="NPlanes_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_mollergem = stemp.Atoi();
	  
	  layer_mollergem = new Int_t[NPlanes_mollergem];
	  nstrips_mollergem = new Int_t[NPlanes_mollergem];
	  offset_mollergem = new Double_t[NPlanes_mollergem];
	  RO_angle_mollergem = new Double_t[NPlanes_mollergem];
	  triggeroffset_mollergem = new Double_t[NPlanes_mollergem/2];
	  gain_mollergem = new Double_t[NPlanes_mollergem/2];
	  nparam_mollergem_read++;
	}
	
	if(skey=="gatewidth_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_mollergem = stemp.Atof();
	  nparam_mollergem_read++;
	}
		
	if(skey=="ZsupThr_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_mollergem = stemp.Atof();
	  nparam_mollergem_read++;
	}

	if(skey=="nlayers_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Nlayers_mollergem = stemp.Atof();
	  nparam_mollergem_read++;
	}
	
	if(skey=="mollergem_layer_z"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  if(ntokens==Nlayers_mollergem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      mollergem_layer_z.push_back(stemp.Atof());
	    }
	  }else{
	    cout << "number of entries for mollergem_layer_z = " << ntokens-1 
		 << " don't match nlayers = " << Nlayers_mollergem << endl;
	    cout << "fix your db " << endl;
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="layer_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_mollergem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      layer_mollergem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for layer_mollergem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_mollergem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="nstrips_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_mollergem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_mollergem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_mollergem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_mollergem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="offset_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_mollergem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_mollergem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_mollergem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_mollergem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="RO_angle_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_mollergem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_mollergem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_mollergem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_mollergem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="triggeroffset_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_mollergem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_mollergem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for triggeroffset_mollergem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_mollergem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="gain_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_mollergem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      gain_mollergem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for gain_mollergem = " << ntokens-1 
		 << " don't match Nplanes/2 = " << NPlanes_mollergem/2 << endl;
	    if(ntokens>=2){
	      cout << "applying first value on all planes " << endl;
	      TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	      for(int k = 0; k<NPlanes_mollergem/2; k++){
		gain_mollergem[k] = stemp.Atof();
	      }
	    }else{
	      cout << "fix your db " << endl;
	      exit(-1);
	    }
	  }
	  nparam_mollergem_read++;
	}
	
	if(skey=="commonmode_array_mollergem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_mollergem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV_mollergem++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_mollergem[k-1] = stemp.Atof();
	  }
	  nparam_mollergem_read++;
	}
	
      }//end if( ntokens >= 2 )
      tokens->~TObjArray();// ineffective... :(
    }//end if( !currentline.BeginsWith("#"))
  }//end while


  //-----------------------------
  //  Declare detectors
  //-----------------------------
  cout << " declaring detectors " << endl;
  for(int k = 0; k<detectors_list.size(); k++){
    cout << "detector: " << detectors_list[k].Data() << "... " << endl;
    if(detectors_list[k] == "mollergem"){
      if(nparam_mollergem_read!=nparam_gemdet){
	cout << detectors_list[k] <<  " does not have the right number of parameters!!! " << endl << " fix database and retry! " << endl;
	exit(-1);
      }

      GEMDet* mollergem = new GEMDet(0, NPlanes_mollergem, layer_mollergem, nstrips_mollergem, offset_mollergem, RO_angle_mollergem, 6, ZsupThr_mollergem);
      GEMSimDigitizer* gemdig = new GEMSimDigitizer(NPlanes_mollergem/2, triggeroffset_mollergem, gain_mollergem, ZsupThr_mollergem, nAPV_mollergem, commonmode_array_mollergem);
      for(int m = 0; m<Nlayers_mollergem; m++){
	mollergem->fZLayer.push_back(mollergem_layer_z[m]);
      }
      mollergem->fGateWidth = gatewidth_mollergem;
      
      GEMdetectors.push_back(mollergem);
      gemdetmap.push_back(0);
      GEMSimDig.push_back(gemdig);

      cout << " set up! " << endl;
    }
  }
  
  TRandom3* R = new TRandom3(Rseed);
  
  // Step 1: read input files build the input chains
  // build signal chain
  ifstream sig_inputfile(inputsigfile);
  TChain *C_s = new TChain("T");
  while( currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C_s->Add(currentline.Data());
    }
  }
  TObjArray *fileElements_s=C_s->GetListOfFiles();
  TIter next_s(fileElements_s);
  TChainElement *chEl_s=0;

  ULong64_t Nev_fs;
  ULong64_t ev_s;
  
  ULong64_t NEventsTotal = 0;
  
  int i_fs = 0;
  bool has_data;
  
  double timeZero;
  
  std::chrono::time_point<std::chrono::steady_clock> start = 
    std::chrono::steady_clock::now();
  
  while (( chEl_s=(TChainElement*)next_s() )) {
    if(NEventsTotal>=Nentries){
      break;
    }
    //open input file in READ mode
    TFile f_s(chEl_s->GetTitle(), "READ");
    if(f_s.IsZombie())cout << "File " << chEl_s->GetTitle() << " cannot be found. Please check the path of your file." << endl; 

    C_s = (TChain*)f_s.Get("T");

    TString outname = TString(chEl_s->GetTitle());
    outname.ReplaceAll(".root", "_dig.root");
    TFile f_out(outname, "RECREATE");
    TTree *T_out = new TTree("T", "Digitized GEM output");


    remoll_tree *T_s = new remoll_tree();
    T_s->Init(C_s);
    T_s->SetupOutputBranches(T_out);
    has_data = false;                            

    Nev_fs = C_s->GetEntries();

	

    //remollDigAuxi auxi;

    for(ev_s = 0; ev_s<Nev_fs; ev_s++, NEventsTotal++){
      if(NEventsTotal>=Nentries)break;
      if(NEventsTotal%100==0)
	cout << NEventsTotal << "/" << Nentries << endl;
      
    
      T_s->ClearDigBranches();
      for(int k = 0; k<GEMdetectors.size(); k++){
		GEMdetectors[k]->Clear();
		}
      timeZero = R->Gaus(0.0, TriggerJitter);

      T_s->GetEntry(ev_s);
      has_data = UnfoldEvent(T_s->hit, R, GEMdetectors,gemdetmap,timeZero,0);

      if(has_data){
        std::cout << "Event = " << ev_s << std::endl;

        for(int k = 0; k<GEMdetectors.size();k++){
	 GEMSimDig[k]->Digitize(GEMdetectors[k],R, 0);
	 GEMSimDig[k]->CheckOut(GEMdetectors[k],gemdetmap[k], R, T_s);
	}
      }
/*
      for (int detID : auxi.GetActiveDetectors()) {

          const auto& hits = auxi.GetHits(detID);
          int counter = 0;
          for (size_t ihit=0; ihit<hits.size();ihit++ ) {
      	  const remollSimHit& h = hits[ihit];
              std::cout << "detID = " << detID
                        << " x = " << h.x
                        << " edep = " << h.edep
                        << " px = " << h.px
                        << " py = " << h.py
                        << " pz = " << h.pz
                        << std::endl;
          }
	}*/
      T_s->FillDigBranches();
      T_out->Fill();
      }
/*
    for(int k = 0; k<GEMdetectors.size(); k++){
      cout << "GEM det ID: " << GEMdetectors[k]->fUniqueID << endl;
      GEMSimDig[k]->write_histos();
      GEMSimDig[k]->print_time_execution();
    }    
*/
    f_out.cd();
    T_out->Write();
    f_out.Close();
    //f_s.Write();
    f_s.Close();
    std::cout<<"Output written to: "<<outname.Data() <<std::endl;
    i_fs++;
  }// end loop on signal files
   
  std::chrono::time_point<std::chrono::steady_clock> end = 
    std::chrono::steady_clock::now();
  
  std::chrono::duration<double> diff = end-start;
  cout << " Total time " << std::setprecision(9) << diff.count() << " s "<< endl;
  
  exit(0);
/*
   TFile *f = TFile::Open(infile, "UPDATE");
   if (!f || f->IsZombie()) {
       std::cerr << "Error opening file\n";
       return 1;
   }

   TTree *rawTree = (TTree*) f->Get("T");
   if (!rawTree) {
       std::cerr << "Tree 'T' not found!\n";
       return 1;
   }

   remoll_tree T;
   T.Init(rawTree);

   remollDigAuxi auxi;

   Long64_t N = T.GetEntries();
   std::cout << "Processing " << N << " events\n";

   //for (Long64_t iev = 0; iev < N; iev++) {
   for (Long64_t iev = 0; iev < 10; iev++) {

       T.GetEntry(iev);

       auxi.ProcessEvent(T.hit);

       std::cout << "Event = " << iev << std::endl;

       for (int detID : auxi.GetActiveDetectors()) {

           const auto& hits = auxi.GetHits(detID);
           int counter = 0;
           for (size_t ihit=0; ihit<hits.size();ihit++ ) {
       	const remollSimHit& h = hits[ihit];
               std::cout << "detID = " << detID
                         << " x = " << h.x
                         << " edep = " << h.edep
                         << " px = " << h.px
                         << " py = " << h.py
                         << " pz = " << h.pz
                         << std::endl;
       	counter++;
           }
           //for (const auto& h : hits ) {
           //  std::cout << "detID = " << detID
               //          << " x = " << h.x
                 //        << " edep = " << h.edep
                   //      << std::endl;
          // }
       }
   T.FillDigBranches();
   }
   f->Write();
   f->Close();
   return 0;
*/
}

