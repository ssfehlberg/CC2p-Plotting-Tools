#define twoproton_pelee_overlay_cxx
#include "twoproton_pelee_overlay.h"
#include "histogram_funcs.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void twoproton_pelee_overlay::Loop()
{

  //Making a new Root File that will contain all the histograms that we will want to plot:
  ///////////////////////////////////////////////////////////////////////////////////////
  TFile *tfile = new TFile("root_files/histograms_filtered_wgt.root","RECREATE");

  //Files with RSE's in them                                                                            
  ofstream myfile;//File that will contain RSE of good events                                          
  ofstream cc2p; //File that will contain good cc2p events                                                                          
  //ofstream ccNp0pi_file;
                     
  myfile.open("lists/files_filtered_wgt.list");
  cc2p.open("lists/files_filtered_wgt_cc2p.list");
  //ccNp0pi_file.open("lists/files_filtered_wgt_ccNp0pi.list"); 
  myfile<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  cc2p<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;
  //ccNp0pi_file<<"Run"<<" "<<"Subrun"<<" "<<"Event"<<endl;  

  /*
  //Define all the histograms I am going to fill                                
  ////////////////////////////////////////////
  Define_Histograms();

  //Defining all the constans we will use later
  //////////////////////////////
  bool _debug = false; //debug statements

  //Counters
  int fvcntr = 0; //Number of events with reconstructed vertex within the FV
  int isfromnucntr = 0; //how many pfp's are from the neutrino slice
  int has3pfp = 0; //how many events have exactly 3 pfps
  int has0shower = 0;//how many events has 0 showers (i.e. three tracks)
  int threetrkcntr = 0; //Number of events with three tracks    
  int vectorsize3 = 0; //Number of events with 3 tracks whose start is < 5cm from the reco vertex
  int secondtrkgood = 0; //Number of events where the second shortest/longest track is contained
  int shortesttrkgood=0; //Number of events where the shortest track is contained
  int events_remaining = 0; //sanity check for number of events remaining
  int pid_cut0 = 0; //sanity pid cut
  int pid_cut1 = 0; //sannity pid cut
  int proton_cut = 0; //number of events with 2 protons in final state
  int muon_cut = 0; //number of events with 1 muon in final state
  int events_chi2p = 0; //sanity checks on chi2
  int events_chi2mu = 0; //sanity cheks on chi2
  int events_mc = 0; //sanity checks o chi2
  int events_2mu = 0; //events with two muons
  int events_2other = 0; //events with 2 others
  int n_mom_mu = 0;
  int n_mom_p1 = 0;
  int n_mom_p2 = 0;

  //neutrino counters
  int neutrinos_0 = 0;
  int neutrinos_1 = 0;
  int neutrinos_else = 0;

  //stupid counters
  int a = 0;
  int b = 0;

  //FV Stuff
  float_t FV_edge = 10.0;
  float_t xmin = 0.0 + FV_edge;
  float_t xmax = 256.35 - FV_edge;
  float_t ymin = -116.5 + FV_edge;
  float_t ymax = 116.5 - FV_edge;
  float_t zmin = 0.0 + FV_edge;
  float_t zmax = 1036.8 - FV_edge;

*/
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout<<"Total Number of Entries: "<<nentries<<std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry); nbytes += nb;

      std::cout<<"-----------------------------------"<<std::endl;
      std::cout<<"BEGINNING TO PROCESS RUN: "<<run<<"  SUBRUN: "<<sub<<"  EVENT: "<<evt<<std::endl;
      std::cout<<"-----------------------------------"<<std::endl;
      


   } //end of loop through entries
} //end of program
