#ifndef analysis_h
#define analysis_h

#include "helper_funcs.h"
#include "paul_tol_colors.hpp"
#include <string>

//////////////////////////////
// Functions to perform plotting are defined here.
// Some helper functions can be found in helper_funcs.h
// Histograms are also defined here to create global variables
//////////////////////////////

class analysis
{

public:
  virtual void main();
  virtual void Grab_Histograms( TFile* f1,  TFile* f2,  TFile* f3, TFile* f4, TFile* f_dirt);
  virtual void Plot_Histograms(Color_t colors[], std::vector<TH1D*> h_overlay, TH1D* h_overlay00,  TH1D* h_overlay1,  TH1D* h_overlay2,  TH1D* h_ext,  TH1D* h_ext1,  TH1D* h_ext2,  TH1D* h_dirt,  TH1D* h_dirt1,  TH1D* h_dirt2,  TH1D* h_bnb,  TH1D* h_bnb1, TCanvas* canv, THStack* h, TPad* pad, TPad* pad0, TLegend* legend, std::vector<const char*> channel_legend, int ymax, int ymin, int num_channels, const char* titles, string path, const char* titles2 = "", const char* plots="", const char* cut = "", bool plot_total = false, bool plot_ccnue = false, double pad_lim = 0.0, double pad0_lim = 0.19,double bnb_min = 0.0, double bnb_max = 2.7);

  

  private:

  //////////////////
  //GENERAL VARIABLES
  //////////////////
  
  //POT number and In Progress
  char const * pot_num="#scale[0.6]{Accumulated POT: 4.566e+19}";//pot number printed on the plots
  char const * sample_name="#scale[0.6]{MicroBooNE In-Progress}";//sample name printed on the plots

  //latex stuff
  TLatex* t = new TLatex();
  
  //Stuff for time/date
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;
  
  
  ////////////////////////////////////////////////////////////////////
  //PLOTS THAT DON'T HAVE EITHER RAQUEL'S OR MY MONTE-CARLO BREAKDOWNS
  ////////////////////////////////////////////////////////////////////

 //Plots of variables concerning the pfp particles:
 /////////////////////////////////////////////////
 static const int num_group = 4;
 const char* group[num_group] = {"npfp","vtx_npfp","ntrack","nshower"};
 const char* titles_pfp[num_group] = {"Number of PFP","Number of PFP Attached to the Vertex","Number of Tracks","Number of Showers"};
 int ylim_pfp[num_group] = {6000,6000,6000,6000}; 
 TH1D* h_overlay_pfp[num_group];
 TH1D* h_bnb_pfp[num_group];
 TH1D* h_ext_pfp[num_group];
 THStack* h_pfp[num_group];
 TCanvas* canv_pfp[num_group];
 TLegend* legend_pfp[num_group];

 //Plots of 2D Variables:
 //////////////////////////
 static const int num_group2d = 3;
 const char* group2d[num_group2d] = {"reco","truth","truth_sce"};
 TH2D* h_overlay2D[num_group2d];
 TCanvas* canv_2d[num_group2d];

 //Plots of Efficiency
 ///////////////////////
 static const int num_eff = 6;
 const char* eff[num_eff] = {"nue","q2","X","Y","muon_mom","proton_mom"};//,"recoil_mom","leading_mom"};
 const char* titles_eff[num_eff] = {"True Neutrino Energy (GeV)","True Q^{2} (GeV^{2})","True X","True Y","True Muon Momentum (GeV/c)","True Proton Momentum (GeV/c)"};//"True Recoil Proton Momentum (GeV/c)","True Leading Proton Momentum (GeV/c)"};
 int ylim_eff[num_eff] = {8000,8000,8000,8000,8000,8000};
 double xlim_eff[num_eff] = {-9999,-9999,-9999,-9999,0.1,0.25};
 TLine* a[num_eff];
 TH1D* h_num[num_eff]; //numerator
 TH1D* h_num0[num_eff]; //numerator clone
 TH1D* h_num1[num_eff]; //numerator clone
 TH1D* h_denom[num_eff]; //denominator
 TH1D* h_denom0[num_eff]; //denominator clone
 TH1D* h_denom1[num_eff]; //denominator clone
 TCanvas* canv_eff[num_eff];
 TCanvas* canv_both[num_eff];
 TLegend* legend_eff[num_eff];

 //Chi2 Plots
 ///////////////////////////////////
 int z0=0; //dumb indices                                                                                                                                                     
 int f0 = 0; //dumb indices
 static const int num_planes = 3;
 static const int num_particles = 6;
 static const int num_hypothesis = 2;
 const char* plane[num_planes] = {"Plane_0","Plane_1","Plane_2"};
 const char* particles[num_particles] = {"total","proton","muon","pion","electron","other"};
 const char* hypothesis[num_hypothesis] = {"chi2p","chi2mu"};
  //const char* channel_legend_chi2[num_particles] = {"Total Overlay","Proton","Muon","Pion","Electron","Other"};
  std::vector<const char*> channel_legend_chi2 = {"Total Overlay","Proton","Muon","Pion","Electron","Other"};
  const char* titles_chi2[num_hypothesis] = {"#chi^{2}_{p}","#chi^{2}_{#mu}"};
 const char* titles_planes[num_planes] ={"U Plane","V Plane","Y Plane"};
 int ylim_chi2[num_planes][num_hypothesis] = {{250,400},
					      {250,250},
					      {600,300}};
  int ymin_chi2[num_planes][num_hypothesis] = {{-2,-2},
					       {-2,-2},
					       {-2,-2}};

  Color_t colors_chi2[6] = {kBlack,kRed,kBlue,kGreen,kYellow,204};
 bool plot_total_chi2 = false;
 TH1D* h_overlay_chi2[num_planes][num_hypothesis][num_particles];
  std::vector<TH1D*> h_overlay_chi2_vec;
  TH1D* h_overlay0_chi2[num_planes][num_hypothesis][3]; //copies for statistics
 TH1D* h_dirt_chi2[num_planes][num_hypothesis][3];
 TH1D* h_ext_chi2[num_planes][num_hypothesis][3]; 
 TH1D* h_bnb_chi2[num_planes][num_hypothesis][2];
 THStack* h_chi2[num_planes][num_hypothesis];
 TCanvas* canv_chi2[num_planes][num_hypothesis];
 TLegend* legend_chi2[num_planes][num_hypothesis];
 TPad* pad_chi2[num_planes][num_hypothesis];
 TPad* pad0_chi2[num_planes][num_hypothesis];

 //3D Chi2 plots
 ///////////////////////////////////////////
 static const int num_cuts_3D = 2;
 const char* cuts_3D[num_cuts_3D] = {"_before_selection","_after_selection"};
 const char* cuts_3D_titles[num_cuts_3D] = {"Before Selection","After Selection"};
 static const int num_hypothesis_3D = 3;
 const char* hypothesis_3D[num_hypothesis_3D] = {"chi2p","chi2mu","chi2pi"};
 const char* titles_chi2_3D[num_hypothesis_3D] = {"3D #chi^{2}_{p}","3D #chi^{2}_{#mu}", "3D #chi^{2}_{#pi}"};
 int ylim_chi2_3D[num_cuts_3D][num_hypothesis_3D] = {{600,300,300},
						     {100,100,100}};
  int ymin_chi2_3D[num_cuts_3D][num_hypothesis_3D] = {{-2,-2,-2},
						      {-2,-2,-2}};
 double xlim_chi2_3D[num_cuts_3D][num_hypothesis_3D] = {{70,-9999,-9999},
							{-9999,-9999,-9999}};
  TLine* a_3D[num_cuts_3D][num_hypothesis_3D];
  bool plot_total_chi2_3D = false;
  TH1D* h_overlay_chi2_3D[num_cuts_3D][num_hypothesis_3D][num_particles];
  std::vector<TH1D*> h_overlay_chi2_3D_vec;
  TH1D* h_overlay0_chi2_3D[num_cuts_3D][num_hypothesis_3D][3]; //copies for statistics
 TH1D* h_dirt_chi2_3D[num_cuts_3D][num_hypothesis_3D][3]; 
 TH1D* h_ext_chi2_3D[num_cuts_3D][num_hypothesis_3D][3]; //copies for statistics
 TH1D* h_bnb_chi2_3D[num_cuts_3D][num_hypothesis_3D][2]; //copies for statistics
 THStack* h_chi2_3D[num_cuts_3D][num_hypothesis_3D];
 TCanvas* canv_chi2_3D[num_cuts_3D][num_hypothesis_3D];
 TLegend* legend_chi2_3D[num_cuts_3D][num_hypothesis_3D];
 TPad* pad_chi2_3D[num_cuts_3D][num_hypothesis_3D];
 TPad* pad0_chi2_3D[num_cuts_3D][num_hypothesis_3D];

 ///////////////////////////////////////
 //PLOTS WITH MY MONTEC-CARLO DEFINITIONS
 ////////////////////////////////////////

 //General
 /////////
 static const int num_cuts = 4; //number of applied cuts
 static const int num_channels = 11; //number of various overlay channels
 const char* cut[num_cuts] = {"_before_selection","_after_3goodtrks","_after_containment","_after_PID"}; 
 const char* channel[num_channels] = {"_total", "_cc2p0pi","_ccNp1pi","_ccNp0pi","_cc1p0pi","_nc","_ccNpNpi","_cc0p0pi","_other","_outfv","_ccnue"};
  //const char* channel_legend[num_channels] = {"Total Overlay","CC2p0#pi (Signal)","CC(N>=0)p1#pi","CC(N>2)p0#pi","CC1p0#pi","NC","CC(N>=0)p(N>1)#pi","CC0p#pi","Other","OOFV","CC#nu_{e}"};
  std::vector<const char*> channel_legend = {"Total Overlay","CC2p0#pi (Signal)","CC(N>=0)p1#pi","CC(N>2)p0#pi","CC1p0#pi","NC","CC(N>=0)p(N>1)#pi","CC0p#pi","Other","OOFV","CC#nu_{e}"};
  const char* cut_titles[num_cuts] = {"Before Selection", "After 3 Track Quality Cuts","After Containment Cut","After PID Cuts"};
  int z = 0; //dumb indices
  int f = 0; //dumb indices

  //Plots of the reconstructed vertex: MC, Dirt, BNB, & EXT
  /////////////////////////////////////////////////////////
  static const int num_variables = 3; //number of various variables to plot
 const char* plots[num_variables] = {"_vtx_x","_vtx_y","_vtx_z"};
 const char* titles[num_variables] = {"Reconstructed X of Vertex (cm)", "Reconstructed Y of Vertex (cm)", "Reconstructed Z of Vertex (cm)"};
 int ylim[num_cuts][num_variables] = {{80,80,80},
				      {80,80,80},
				      {80,80,80},
				      {40,40,40}};
  int ymin[num_cuts][num_variables] = {{-3,-3,-4},
				      {-4,-3,-3},
		 		      {-3,-3,-3},
				      {-2,-2,-2}};
  bool plot_total = false; //sanity check on total overlay
  bool plot_ccnue = false; //do you want to plot ccnue?
  std::vector<TH1D*> h_overlay_vec; //overlay plots
  TH1D* h_overlay[num_cuts][num_variables][num_channels];
  TH1D* h_overlay0[num_cuts][num_variables][3]; //copies of the overlay plots for statistics
  TH1D* h_dirt[num_cuts][num_variables][3];
  TH1D* h_bnb[num_cuts][num_variables][2]; //two copies
  TH1D* h_ext[num_cuts][num_variables][3]; //two copies
  THStack* h[num_cuts][num_variables];
  TCanvas* canv[num_cuts][num_variables];
  TLegend* legend[num_cuts][num_variables];
  TPad* pad[num_cuts][num_variables];
  TPad* pad0[num_cuts][num_variables];

 //Plots of reconstructed vertex: Truth Only
 /////////////////////////////////////////////
 static const int num_truth = 9;//10;
 const char* truth[num_truth] = {"_vtx_x_mc","_vtx_y_mc","_vtx_z_mc","_vtx_x_mc_sce","_vtx_y_mc_sce","_vtx_z_mc_sce","_q2","_X","_Y"};
 const char* truth_titles[num_truth] = {"True Vertex X (cm)","True Vertex Y (cm)","True Vertex Z (cm)",
					"True Vertex w/ SCE X (cm)","True Vertex w/ SCE Y (cm)","True Vertex w/ SCE Z (cm)",
					"True Q^{2} (GeV^{2})","True Bjorken X","True Bjorken Y"};
 int ylim_truth[num_cuts][num_truth] = {{80,80,80,80,80,80,200,200,200},
					{60,60,60,60,60,60,200,200,200},
					{60,60,60,60,60,60,200,200,200},
					{60,60,60,60,60,60,200,200,200}};
 bool plot_total_truth = true;
 bool plot_ccnue_truth = true;
 TH1D* h_mc[num_cuts][num_truth][num_channels];
  TCanvas* canv_truth[num_cuts][num_truth];
 THStack* h_truth[num_cuts][num_truth];
 TLegend* legend_truth[num_cuts][num_truth];

 //Plots of the individual particle quantities
 ////////////////////////////////////////////
 static const int num_var = 4;
 const char* var[num_var] = {"_mom","_E","_theta","_phi"};
 const char* titles_var[num_var] = {"Momentum (GeV/c)","Energy (GeV)","cos(#theta)","#phi (Rad)"};

 //Muon
 TH1D* h_muon_overlay[num_var][num_channels]; //actual overlay
  std::vector<TH1D*> h_muon_overlay_vec;
  TH1D* h_muon_overlay0[num_var][3]; //copies for statistics
 TH1D* h_muon_bnb[num_var][2];
 TH1D* h_muon_ext[num_var][3];
 TH1D* h_muon_dirt[num_var][3];
 THStack* h_muon[num_var];
 TCanvas* canv_muon[num_var];
 TLegend* legend_muon[num_var];
 TPad* pad_muon[num_var];
 TPad* pad0_muon[num_var];
 int muon_ylim[num_var] = {40,40,160,100};
  int muon_ymin[num_var] = {-2,-2,-2,-2};
  
 //Recoil Proton
 TH1D* h_recoil_overlay[num_var][num_channels];
  std::vector<TH1D*> h_recoil_overlay_vec;
 TH1D* h_recoil_overlay0[num_var][3];
 TH1D* h_recoil_bnb[num_var][2];
 TH1D* h_recoil_ext[num_var][3];
 TH1D* h_recoil_dirt[num_var][3];
 THStack* h_recoil[num_var];
 TCanvas* canv_recoil[num_var];
 TLegend* legend_recoil[num_var];
 TPad* pad_recoil[num_var];
 TPad* pad0_recoil[num_var];
 double recoil_ylim[num_var] = {60,60,80,150};
  int recoil_ymin[num_var] = {-2,-2,-2,-2};

 //Leading Proton
 TH1D* h_leading_overlay[num_var][num_channels];
  std::vector<TH1D*> h_leading_overlay_vec;
 TH1D* h_leading_overlay0[num_var][3];
 TH1D* h_leading_bnb[num_var][2];
 TH1D* h_leading_ext[num_var][3];
 TH1D* h_leading_dirt[num_var][3];
 THStack* h_leading[num_var];
 TCanvas* canv_leading[num_var];
 TLegend* legend_leading[num_var];
 TPad* pad_leading[num_var];
 TPad* pad0_leading[num_var];
 double leading_ylim[num_var] = {50,50,160,100};
  int leading_ymin[num_var] = {-2,-2,-2,-2};

 //Random Physics Variables
 ///////////////////////////
  static const int num_phys =5; //NOTE: These values are the same for Raquel's plots. If this changes in the future, be sure to chage this!
  const char* physics[num_phys] = {"_cos_gamma_cm","_opening_angle_protons","_opening_angle_mu_leading","_mom_struck_nuc","_tot_pz"};
  const char* physics_titles[num_phys] = {"cos(#gamma_{cm})","cos(#gamma_{Lab})","Opening Angle Between the Muon and Leading Proton (Rad.)","Momentum of Struck Nucleon (GeV)","Total P_{z} of the Two Protons"};
  int phys_ylim[num_phys] =  {130,60,60,60,80};
  int phys_ymin[num_phys] =  {-2,-2,-2,-2};
  TH1D* h_phys_overlay[num_phys][num_channels];
  std::vector<TH1D*> h_phys_overlay_vec;
  TH1D* h_phys_overlay0[num_phys][3];
 TH1D* h_phys_bnb[num_phys][2];
 TH1D* h_phys_ext[num_phys][3];
 TH1D* h_phys_dirt[num_phys][3];
 THStack* h_phys[num_phys];
 TCanvas* canv_phys[num_phys];
 TLegend* legend_phys[num_phys];
 TPad* pad_phys[num_phys];
 TPad* pad0_phys[num_phys];

  //STVs
  //////////////////////////
  static const int num_stv =3;
  const char* stv[num_stv] = {"_delta_PT","_delta_phiT","_delta_alphaT"};
  const char* stv_titles[num_stv] = {"#delta P_{T} (GeV/c)","#delta #phi_{T} (Deg.)","#delta #alpha_{T} (Deg.)"};
  int stv_ylim[num_stv] = {130,150,120};
  int stv_ymin[num_stv] = {-2,-2,-2};  
  TH1D* h_stv_overlay[num_stv][num_channels];
  std::vector<TH1D*> h_stv_overlay_vec;
  TH1D* h_stv_overlay0[num_stv][3];
  TH1D* h_stv_bnb[num_stv][2];
  TH1D* h_stv_ext[num_stv][3];
  TH1D* h_stv_dirt[num_stv][3];
  THStack* h_stv[num_stv];
  TCanvas* canv_stv[num_stv];
  TLegend* legend_stv[num_stv];
  TPad* pad_stv[num_stv];
  TPad* pad0_stv[num_stv];

 //////////////////////////////////////////////////
 //PLOTS WITH RAQUEL'S MOTE-CARLO DEFINITIONS
 //////////////////////////////////////////////////

 //General
 //////////////////////////
 bool plot_total_raquel = false; //sanity check on total overlay                                                                                                                                                                       
 bool plot_ccnue_raquel = false; //do you want to plot ccnue? 
 int z_raquel=0; //dumb indices                                                                                                                                                                                                        
 int f_raquel = 0; //dumb indices 
 static const int num_channels_raquel = 10;
 const char* channel_raquel[num_channels_raquel] = {"_total", "_ccRES","_ccQE","_ccMEC","_nc","_ccDIS","_ccCOH","_outfv","_other","_ccNue"};
  //const char* channel_legend_raquel[num_channels_raquel] = {"Total Overlay","CCRES","CCQE","CCMEC","NC","CCDIS","CCCOH","OOFV","Other","CC$\nu_{e}$"};
  std::vector<const char*> channel_legend_raquel = {"Total Overlay","CCRES","CCQE","CCMEC","NC","CCDIS","CCCOH","OOFV","Other","CC$\nu_{e}$"};
  
 //Plots of the individual particle quantities
 //Note: we are taking all the variables from my definintions since they are the same for Raquel's plots                                                                                                                              
 //We are just defining new cavases, legends, pads, and overlay histograms 
 ////////////////////////////////////////////

 //Muon
 TH1D* h_muon_overlay_raquel[num_var][num_channels_raquel]; //actual overlay
  std::vector<TH1D*> h_muon_overlay_raquel_vec;
 TH1D* h_muon_overlay0_raquel[num_var][3]; //copies for statistics
 THStack* h_muon_raquel[num_var];
 TCanvas* canv_muon_raquel[num_var];
 TLegend* legend_muon_raquel[num_var];
 TPad* pad_muon_raquel[num_var];
 TPad* pad0_muon_raquel[num_var];

 //Recoil Proton
 TH1D* h_recoil_overlay_raquel[num_var][num_channels_raquel];
  std::vector<TH1D*> h_recoil_overlay_raquel_vec;
 TH1D* h_recoil_overlay0_raquel[num_var][3];
 THStack* h_recoil_raquel[num_var];
 TCanvas* canv_recoil_raquel[num_var];
 TLegend* legend_recoil_raquel[num_var];
 TPad* pad_recoil_raquel[num_var];
 TPad* pad0_recoil_raquel[num_var];

 //Leading Proton
 TH1D* h_leading_overlay_raquel[num_var][num_channels_raquel];
  std::vector<TH1D*> h_leading_overlay_raquel_vec;
 TH1D* h_leading_overlay0_raquel[num_var][3];
 THStack* h_leading_raquel[num_var];
 TCanvas* canv_leading_raquel[num_var];
 TLegend* legend_leading_raquel[num_var];
 TPad* pad_leading_raquel[num_var];
 TPad* pad0_leading_raquel[num_var];

 //Rando plots
 //Note: we are taking all the variables from my definintions since they are the same for Raquel's plots
 //We are just defining new cavases, legends, pads, and overlay histograms 
 ////////////////////////////////////////////////////////
 TH1D* h_phys_overlay_raquel[num_phys][num_channels_raquel];
  std::vector<TH1D*> h_phys_overlay_raquel_vec;
 TH1D* h_phys_overlay0_raquel[num_phys][3];
 THStack* h_phys_raquel[num_phys];
 TCanvas* canv_phys_raquel[num_phys];
 TLegend* legend_phys_raquel[num_phys];
 TPad* pad_phys_raquel[num_phys];
 TPad* pad0_phys_raquel[num_phys];

 //STVs
 //Note: we are taking all the variables from my definintions since they are the same for Raquel's plots
 //We are just defining new cavases, legends, pads, and overlay histograms
 ////////////////////////////   
 TH1D* h_stv_overlay_raquel[num_stv][num_channels_raquel];
  std::vector<TH1D*> h_stv_overlay_raquel_vec;
 TH1D* h_stv_overlay0_raquel[num_stv][3];
 THStack* h_stv_raquel[num_stv];
 TCanvas* canv_stv_raquel[num_stv];
 TLegend* legend_stv_raquel[num_stv];
 TPad* pad_stv_raquel[num_stv];
 TPad* pad0_stv_raquel[num_stv];		

}; //end of class definition

#endif

#ifdef analysis_cxx

//////////////////////////
//FUNCTIONS
////////////////////////
void analysis::Grab_Histograms( TFile* f1,  TFile* f2,  TFile* f3, TFile* f4, TFile* f_dirt){
  for(int i = 0; i < num_cuts; i++){
    for(int j = 0; j < num_variables; j++){
      h_bnb[i][j][0] = (TH1D*)f2->Get(Form("h%s%s_bnb",plots[j],cut[i]));
      h_ext[i][j][0] = (TH1D*)f3->Get(Form("h%s%s_ext",plots[j],cut[i]));
      h_dirt[i][j][0] = (TH1D*)f_dirt->Get(Form("h%s%s_dirt",plots[j],cut[i]));
      for(int k = 0; k < num_channels; k++){
	h_overlay[i][j][k] = (TH1D*)f1->Get(Form("h%s%s%s",plots[j],cut[i],channel[k]));
      }
    }
    
    //grabbing truth stuff
    for(int j=0; j < num_truth; j++){
      for(int k=0; k < num_channels; k++){
	h_mc[i][j][k] = (TH1D*)f1->Get(Form("h%s%s%s",truth[j],cut[i],channel[k]));
      }
    } 
  }

  //grabbing things related the PFP's
  for(int i=0; i < num_group; i++){
    h_overlay_pfp[i] = (TH1D*)f1->Get(Form("h_%s_overlay",group[i]));
    h_bnb_pfp[i] = (TH1D*)f2->Get(Form("h_%s_overlay",group[i]));
    h_ext_pfp[i] =(TH1D*)f3->Get(Form("h_%s_overlay",group[i]));
  }
  
  //grabbing the 2D histograms
  for(int i=0; i < num_group2d; i++){
    h_overlay2D[i] = (TH2D*)f1->Get(Form("h_correlation_overlay_%s",group2d[i])); 
  }

  //grabbing efficiency plots:
  for(int i=0; i < num_eff; i++){
    h_num[i] =  (TH1D*)f4->Get(Form("h_num_overlay_%s",eff[i]));
    h_denom[i] =  (TH1D*)f4->Get(Form("h_denom_overlay_%s",eff[i]));
    
  }
  
  //grabbing the chi2 plots
  for(int i = 0; i < num_planes; i ++){
    for(int j = 0; j < num_hypothesis; j++){
      h_ext_chi2[i][j][0] = (TH1D*)f3->Get(Form("h_%s_%s_ext",hypothesis[j],plane[i]));
      h_bnb_chi2[i][j][0] = (TH1D*)f2->Get(Form("h_%s_%s_bnb",hypothesis[j],plane[i]));
      h_dirt_chi2[i][j][0] = (TH1D*)f_dirt->Get(Form("h_%s_%s_dirt",hypothesis[j],plane[i]));    
      for(int k =0; k < num_particles; k++){
	h_overlay_chi2[i][j][k] = (TH1D*)f1->Get(Form("h_%s_%s_%s",hypothesis[j],plane[i],particles[k]));
      }
    }
  }
  
  //grabbing 3D chi2 plots
  for(int i = 0; i < num_cuts_3D; i++){
    for(int j = 0; j < num_hypothesis_3D; j++){
      h_ext_chi2_3D[i][j][0] = (TH1D*)f3->Get(Form("h_%s_3D%s_ext",hypothesis_3D[j],cuts_3D[i]));
      h_bnb_chi2_3D[i][j][0] = (TH1D*)f2->Get(Form("h_%s_3D%s_bnb",hypothesis_3D[j],cuts_3D[i]));
      h_dirt_chi2_3D[i][j][0] = (TH1D*)f_dirt->Get(Form("h_%s_3D%s_dirt",hypothesis_3D[j],cuts_3D[i]));    
      for(int k =0; k < num_particles; k++){
	h_overlay_chi2_3D[i][j][k] = (TH1D*)f1->Get(Form("h_%s_3D%s_%s",hypothesis_3D[j],cuts_3D[i],particles[k]));
      }
    }
  }
  
  //grabbing particle plots
  for(int i = 0; i < num_var; i++){
    h_muon_ext[i][0] = (TH1D*)f3->Get(Form("h_muon%s",var[i]));
    h_recoil_ext[i][0] = (TH1D*)f3->Get(Form("h_recoil%s",var[i]));
    h_leading_ext[i][0] = (TH1D*)f3->Get(Form("h_leading%s",var[i]));
    h_muon_bnb[i][0] = (TH1D*)f2->Get(Form("h_muon%s",var[i]));
    h_recoil_bnb[i][0] = (TH1D*)f2->Get(Form("h_recoil%s",var[i]));
    h_leading_bnb[i][0] = (TH1D*)f2->Get(Form("h_leading%s",var[i]));
    h_muon_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h_muon%s",var[i]));
    h_recoil_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h_recoil%s",var[i]));
    h_leading_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h_leading%s",var[i]));
    for(int j =0; j < num_channels; j++){
      h_muon_overlay[i][j] = (TH1D*)f1->Get(Form("h_muon%s%s",var[i],channel[j]));
      h_recoil_overlay[i][j] = (TH1D*)f1->Get(Form("h_recoil%s%s",var[i],channel[j]));
      h_leading_overlay[i][j] = (TH1D*)f1->Get(Form("h_leading%s%s",var[i],channel[j]));
    }
  }

  //grabbing physics plots
  for(int i=0; i < num_phys; i++){
    h_phys_ext[i][0] = (TH1D*)f3->Get(Form("h%s",physics[i]));
    h_phys_bnb[i][0] = (TH1D*)f2->Get(Form("h%s",physics[i]));
    h_phys_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h%s",physics[i]));
    for(int j=0; j < num_channels; j++){
      h_phys_overlay[i][j] = (TH1D*)f1->Get(Form("h%s%s",physics[i],channel[j])); 
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_phys_overlay_raquel[i][j] = (TH1D*)f1->Get(Form("h%s_raquel%s",physics[i],channel_raquel[j])); 
    }  
  }

  //grabbing stvs plots
  for(int i=0; i < num_stv; i++){
    h_stv_ext[i][0] = (TH1D*)f3->Get(Form("h%s",stv[i]));
    h_stv_bnb[i][0] = (TH1D*)f2->Get(Form("h%s",stv[i]));
    h_stv_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h%s",stv[i]));
    for(int j=0; j < num_channels; j++){
      h_stv_overlay[i][j] = (TH1D*)f1->Get(Form("h%s%s",stv[i],channel[j])); 
    }
   for(int j=0; j < num_channels_raquel; j++){
      h_stv_overlay_raquel[i][j] = (TH1D*)f1->Get(Form("h%s_raquel%s",stv[i],channel_raquel[j]));
    }    
  }

}

void analysis::Plot_Histograms(Color_t colors[],std::vector<TH1D*> h_overlay, TH1D* h_overlay00,  TH1D* h_overlay1,  TH1D* h_overlay2,  TH1D* h_ext,  TH1D* h_ext1,  TH1D* h_ext2,  TH1D* h_dirt,  TH1D* h_dirt1,  TH1D* h_dirt2,  TH1D* h_bnb,  TH1D* h_bnb1,TCanvas* canv, THStack* h, TPad* pad, TPad* pad0, TLegend* legend, std::vector<const char*> channel_legend,int ymax, int ymin, int num_channels, const char* titles, string path, const char* titles2 ="", const char* plots="", const char* cut = "", bool plot_total = false, bool plot_ccnue = false, double pad_lim = 0.0, double pad0_lim = 0.19, double bnb_min = 0.0, double bnb_max = 2.7)
{      

  h_ext1 = (TH1D*)h_ext->Clone();
  h_ext2 = (TH1D*)h_ext->Clone();
  h_dirt1 = (TH1D*)h_dirt->Clone();
  h_dirt2 = (TH1D*)h_dirt->Clone();
  h_bnb1 = (TH1D*)h_bnb->Clone();
  h_overlay00 = (TH1D*)h_overlay[0]->Clone();
  h_overlay1 = (TH1D*)h_overlay[0]->Clone();
  h_overlay2 = (TH1D*)h_overlay[0]->Clone();
      
  canv = new TCanvas(Form("C%s%s",plots,cut),Form("C%s%s",plots,cut),2000,1500);
  canv->cd(1);
  h = new THStack(Form("h%s%s",plots,cut),Form("h%s%s",plots,cut));
  pad = new TPad(Form("pad%s%s",plots,cut),Form("pad%s%s",plots,cut),0,0.35,1,1);
  pad->SetBottomMargin(pad_lim);
  pad->SetGridx();
  pad->SetBorderMode(0);
  pad->Draw();
  pad->cd();
  
  //Stacked Histrogram parameters
  h->Draw("HIST");
  h->SetTitle("");
  h->SetMaximum(ymax); 
  h->SetMinimum(ymin);
  
  //MC
  int z;
  if(plot_ccnue){
    z = num_channels;
  }else{ z = num_channels - 1;}
  int f = num_channels-1;
      
  for(int k=1; k < z ; k++){
    h_overlay[k]->SetLineColor(colors[k]);
    h_overlay[k]->SetFillColor(colors[k]);
    h_overlay[k]->SetLineWidth(1);
    h->Add(h_overlay[k]);
    }
  
  //Dirt
  h->Add(h_dirt);
  h_dirt->SetFillColor(kOrange-8);
  h_dirt->SetLineColor(kOrange-8);
  h_dirt->SetLineWidth(1);
  
  //EXT
  h->Add(h_ext);
  h_ext->SetFillColor(kViolet-7);
  h_ext->SetFillStyle(3005);
  h_ext->SetLineWidth(1);
  
  //BNB
  h_bnb->Draw("e1SAME");
  h_bnb->SetLineColor(kBlack);
  h_bnb->SetLineWidth(1);

  h->GetYaxis()->SetTitle("No. Events");
  
  //if you want to plot the total for sanity sake:
  if(plot_total){
    h_overlay[0]->Draw("SAME");
  }
      
  //Make sure to do overlay statistical uncertainty
  OverlayStatistics(h_overlay00,h_ext2, h_dirt2);
  h_overlay00->Draw("e2SAME");
  h_overlay00->SetLineColor(kBlack);
  h_overlay00->SetFillColor(kBlack);
  h_overlay00->SetFillStyle(3004);
  h_overlay00->SetMarkerSize(0);
  h_overlay00->SetLineWidth(1);
  
  legend = new TLegend(0.71, 0.54, 0.899, 0.89);
  legend->AddEntry(h_bnb,"Data (Beam-On)","lepf");
  legend->AddEntry(h_overlay00,"Stat. Unc.","f");
  legend->AddEntry(h_ext,"Data (Beam-Off)","f");
  legend->AddEntry(h_dirt,"Dirt","f");
  for(int k =1; k < z; k++){	
  legend->AddEntry(h_overlay[f-k],Form("%s",channel_legend[f-k]),"f");
  }
  legend->SetLineWidth(0);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.03);
  legend->Draw("same");
  t->SetNDC();
  t->SetTextAlign(22);
  t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s: %s}",titles,titles2));
  t->DrawLatex(0.195,0.92,Form("%s",pot_num));
  t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
  canv->cd();
  pad0 = new TPad(Form("pad0%s%s",plots,cut),Form("pad0%s%s",plots,cut),0,0,1,0.35);
  pad0->SetTopMargin(0);                
  pad0->SetBottomMargin(pad0_lim);
  pad0->SetGridx();
  pad0->Draw();
  pad0->cd();

  h_ext1->Add(h_overlay1); //causing problems
  h_ext1->Add(h_dirt1);
  h_ext1->Sumw2();
  h_bnb1->Divide(h_bnb1,h_ext1,1,1);
  h_bnb1->Sumw2();
  h_bnb1->Draw("e1p");
  h_bnb1->SetStats(kFALSE);
  h_bnb1->SetTitle("");
      
  TF1 *a = new TF1("a","1", -150000 , 150000);
  a->SetLineColor(kRed);
  a->Draw("SAME");

  //Calculate the Chi2
  double chisqv1 = calculatePearsonChiSq(h_bnb, h_overlay[0]);
  int nBins1 = h_overlay[0]->GetXaxis()->GetNbins();
  //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
  t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1/nBins1));
  t->SetNDC();
  t->SetTextAlign(22);
 
  h_bnb1->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
  h_bnb1->GetYaxis()->CenterTitle();
  h_bnb1->GetYaxis()->SetTitleSize(28);
  h_bnb1->GetYaxis()->SetTitleFont(43);
  h_bnb1->GetYaxis()->SetTitleOffset(1.5);
  h_bnb1->GetYaxis()->SetLabelFont(43);
  h_bnb1->GetYaxis()->SetLabelSize(30);
  h_bnb1->GetXaxis()->SetTitle(Form("%s",titles));
  h_bnb1->GetXaxis()->SetTitleSize(35);
  h_bnb1->GetXaxis()->SetTitleFont(43);
  h_bnb1->GetXaxis()->SetTitleOffset(3);
  h_bnb1->GetXaxis()->SetLabelFont(43);
  h_bnb1->GetXaxis()->SetLabelSize(35);
  h_bnb1->SetMaximum(bnb_max);
  h_bnb1->SetMinimum(bnb_min);
      
  canv->Print(Form("%s%s%s.png",path.c_str(),plots,cut));
  canv->Print(Form("%s%s%s.pdf",path.c_str(),plots,cut));
  
} //end of plot histograms

#endif
