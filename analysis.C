#include "analysis.h"
#include "paul_tol_colors.hpp"
#include <iostream>
#include <ctime>
#include <string>

void analysis(){

  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);

  //Define All the Histogram Files
  /////////////////////////////////////
  TFile *f1=new TFile("histograms_filtered.root");//overlay histograms
  TFile *f2=new TFile("histograms_filtered_bnb.root");//bnb histograms  
  TFile *f3=new TFile("histograms_filtered_ext.root");//extbnb histograms 
  TFile *f4=new TFile("histograms_efficiency.root"); //efficiency histograms

  //Things Common to all the Plots
  //////////////////////////////////////////
  tolcols::init();
  Color_t colors[] = {0,9031, 9030, 9029, 9028, 9026, 9025, 9024, 9032, kGray+2, 9027};  //Black, light pink, light orange, light yellow, olive, mint, light cyan, light blue, light grey, darker grey, olive
  
  //POT number and In Progress
  char const * pot_num="#scale[0.7]{Accumulated POT: 3.575e+19}";//pot number printed on the plots
  char const * sample_name="#scale[0.7]{MicroBooNE In-Progress}";//sample name printed on the plots

  //Stuff for time/date
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;

  //Make directory with Today's Date:
  const char* pathname = Form("images/%d%d%d/",Month,Day,Year);
  string path(pathname);
  int dir_exists = dirExists(pathname);
  if(dir_exists == 0){
    mkdir(pathname,0777);
    std::cout<<"New Directory Succesfully Created"<<std::endl;
  } else if (dir_exists < 0){
    std::cout<<"An Error has Occured. Please Check the Input Path Name."<<std::endl;
  }else if(dir_exists > 0){
    std::cout<<"Directory Already Exists. Continuing with Analysis"<<std::endl;
    }
  
  //latex stuff
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  
  //Plots of all the reconstructed varaibles with MC, EXT, and BNB
  ////////////////////////////////////////////////////////////////
  int z=0; //dumb indices
  int f = 0; //dumb indices

  const int num_cuts = 2; //number of applied cuts
  const int num_variables = 3; //number of various variables to plot 
  const int num_channels = 11; //number of various overlay channels
  const char* cut[num_cuts] = {"_before_selection","_after_selection"};
  const char* plots[num_variables] = {"_vtx_x","_vtx_y","_vtx_z"};
  const char* channel[num_channels] = {"_total", "_cc2p0pi","_ccNp1pi","_ccNp0pi","_cc1p0pi","_nc","_ccNpNpi","_cc0p0pi","_other","_outfv","_ccnue"};
  const char* channel_legend[num_channels] = {"Total Overlay","CC2p0#pi (Signal)","CC(N>=0)p1#pi","CC(N>2p)0#pi","CC1p0#pi","NC","CC(N>=0)p(N>2)#pi","CC0p#pi","Other","OOFV","CC#nu_{e}"};
  const char* cut_titles[num_cuts] = {"Before Selection", "After Selection"};
  const char* titles[num_variables] = {"Reconstructed X of Vertex (cm)", "Reconstructed Y of Vertex (cm)", "Reconstructed Z of Vertex (cm)"};
  int ylim[num_cuts][num_variables] = {{80,80,80},{80,80,80}};
  bool plot_total = false; //sanity check on total overlay
  bool plot_ccnue = false; //do you want to plot ccnue?
  TH1D* h_overlay[num_cuts][num_variables][num_channels]; //overlay plots
  TH1D* h_overlay0[num_cuts][num_variables][3]; //copies of the overlay plots for statistics
  TH1D* h_bnb[num_cuts][num_variables][2]; //two copies  
  TH1D* h_ext[num_cuts][num_variables][3]; //two copies
  THStack* h[num_cuts][num_variables];
  TCanvas* canv[num_cuts][num_variables];
  TLegend* legend[num_cuts][num_variables];
  TPad* pad[num_cuts][num_variables];
  TPad* pad0[num_cuts][num_variables];

  //Plots of Truth Only
  //////////////////////////////
  const int num_truth = 9;//10;
  const char* truth[num_truth] = {"_vtx_x_mc","_vtx_y_mc","_vtx_z_mc","_vtx_x_mc_sce","_vtx_y_mc_sce","_vtx_z_mc_sce","_q2","_X","_Y"};
  const char* truth_titles[num_truth] = {"True Vertex X (cm)","True Vertex Y (cm)","True Vertex Z (cm)",
					 "True Vertex w/ SCE X (cm)","True Vertex w/ SCE Y (cm)","True Vertex w/ SCE Z (cm)",
                                         "True Q^{2} (GeV^{2})","True Bjorken X","True Bjorken Y"};
  int ylim_truth[num_cuts][num_truth] = {{80,80,80,80,80,80,200,200,200},{60,60,60,60,60,60,200,200,200}};
  bool plot_total_truth = true;
  bool plot_ccnue_truth = true;
  TH1D* h_mc[num_cuts][num_truth][num_channels];
  TCanvas* canv_truth[num_cuts][num_truth];
  THStack* h_truth[num_cuts][num_truth];
  TLegend* legend_truth[num_cuts][num_truth];
  
  //Plots of variables concerning the pfp particles:
  /////////////////////////////////////////////////
  const int num_group = 4;
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
  const int num_group2d = 3;
  const char* group2d[num_group2d] = {"reco","truth","truth_sce"};
  TH2D* h_overlay2D[num_group2d];
  TCanvas* canv_2d[num_group2d];

  //Plots of Efficiency
  ///////////////////////
  const int num_eff = 4;
  const char* eff[num_eff] = {"nue","q2","X","Y"};
  const char* titles_eff[num_eff] = {"True Neutrino Energy","True Q^{2} (GeV^{2})","True X","True Y"};
  int ylim_eff[num_eff] = {8000,8000,8000,8000};
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
  ///////////////////////////
  int z0=0; //dumb indices                                                                                                                                                     
  int f0 = 0; //dumb indices
  const int num_planes = 3;
  const int num_particles = 6;
  const int num_hypothesis = 2;
  const char* plane[num_planes] = {"Plane_0","Plane_1","Plane_2"};
  const char* particles[num_particles] = {"total","proton","muon","pion","electron","other"};
  const char* hypothesis[num_hypothesis] = {"chi2p","chi2mu"};
  const char* channel_legend_chi2[num_particles] = {"Total Overlay","Proton","Muon","Pion","Electron","Other"};
  const char* titles_chi2[num_hypothesis] = {"#chi^{2}_{p}","#chi^{2}_{#mu}"};
  const char* titles_planes[num_planes] ={"U Plane","V Plane","Y Plane"};
  int ylim_chi2[num_planes][num_hypothesis] = {{250,400},{250,250},{150,250}};
  Color_t colors_chi2[] = {kBlack,kRed,kBlue,kGreen,kYellow,204};
  bool plot_total_chi2 = false;
  TH1D* h_overlay_chi2[num_planes][num_hypothesis][num_particles];
  TH1D* h_overlay0_chi2[num_planes][num_hypothesis][3]; //copies for statistics
  TH1D* h_ext_chi2[num_planes][num_hypothesis][3]; 
  TH1D* h_bnb_chi2[num_planes][num_hypothesis][2];
  THStack* h_chi2[num_planes][num_hypothesis];
  TCanvas* canv_chi2[num_planes][num_hypothesis];
  TLegend* legend_chi2[num_planes][num_hypothesis];
  TPad* pad_chi2[num_planes][num_hypothesis];
  TPad* pad0_chi2[num_planes][num_hypothesis];
  
  
  //Plots from Raquels Files:
  //////////////////////////
  bool plot_raquel = false;
  
  //Grab all the histograms from the files
  //////////////////////////////////////
  for(int i = 0; i < num_cuts; i++){
    for(int j = 0; j < num_variables; j++){
      h_bnb[i][j][0] = (TH1D*)f2->Get(Form("h%s%s_bnb",plots[j],cut[i]));
      h_ext[i][j][0] = (TH1D*)f3->Get(Form("h%s%s_ext",plots[j],cut[i]));

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
      
      for(int k =0; k < num_particles; k++){
	h_overlay_chi2[i][j][k] = (TH1D*)f1->Get(Form("h_%s_%s_%s",hypothesis[j],plane[i],particles[k]));
	
      }

    }
  }
  
  //Now to make some plots: Reconstructed MC,EXT,and BNB
  ////////////////////////////////////////
  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_variables; j++){
      
      h_ext[i][j][1] = (TH1D*)h_ext[i][j][0]->Clone();
      h_ext[i][j][2] = (TH1D*)h_ext[i][j][0]->Clone();
      h_bnb[i][j][1] = (TH1D*)h_bnb[i][j][0]->Clone();
      h_overlay0[i][j][0] = (TH1D*)h_overlay[i][j][0]->Clone();
      h_overlay0[i][j][1] = (TH1D*)h_overlay[i][j][0]->Clone();
      h_overlay0[i][j][2] = (TH1D*)h_overlay[i][j][0]->Clone();
      
      canv[i][j] = new TCanvas(Form("C%s%s",plots[j],cut[i]),Form("C%s%s",plots[j],cut[i]),2000,1500);
      canv[i][j]->cd(1);
      h[i][j] = new THStack(Form("h%s%s",plots[j],cut[i]),Form("h%s%s",plots[j],cut[i]));
      pad[i][j] = new TPad(Form("pad%s%s",plots[j],cut[i]),Form("pad%s%s",plots[j],cut[i]),0,0.28,1.0,1.0);
      pad[i][j]->SetBottomMargin(0.0);
      pad[i][j]->SetGridx();
      pad[i][j]->SetBorderMode(0);
      pad[i][j]->Draw();
      pad[i][j]->cd();

      //Stacked Histrogram parameters
      h[i][j]->Draw("HIST");
      h[i][j]->SetTitle("");
      h[i][j]->SetMaximum(ylim[i][j]); 

      //MC
      if(plot_ccnue){
	z = num_channels;
	f = 11;
      }else{
	z = num_channels -1;
	f = 10;
      }
      
      for(int k=1; k < z ; k++){
	h_overlay[i][j][k]->SetLineColor(colors[k]);
	h_overlay[i][j][k]->SetFillColor(colors[k]);
	h_overlay[i][j][k]->SetLineWidth(1);
	h[i][j]->Add(h_overlay[i][j][k]);
      }

      //EXT
      h[i][j]->Add(h_ext[i][j][0]);
      h_ext[i][j][0]->SetFillColor(kViolet-7);
      h_ext[i][j][0]->SetFillStyle(3005);
      h_ext[i][j][0]->SetLineWidth(1);
      
      //BNB
      h_bnb[i][j][0]->Draw("e1SAME");
      h_bnb[i][j][0]->SetLineColor(kBlack);
      h_bnb[i][j][0]->SetLineWidth(1);
      
      //if you want to plot the total for sanity sake:
      if(plot_total){
      h_overlay[i][j][0]->Draw("SAME");
      }
      
      //Make sure to do overlay statistical uncertainty
      OverlayStatistics(h_overlay0[i][j][0],h_ext[i][j][2]);
      h_overlay0[i][j][0]->Draw("e2SAME");
      h_overlay0[i][j][0]->SetLineColor(kBlack);
      h_overlay0[i][j][0]->SetFillColor(kBlack);
      h_overlay0[i][j][0]->SetFillStyle(3004);
      h_overlay0[i][j][0]->SetMarkerSize(0);
      h_overlay0[i][j][0]->SetLineWidth(1);

      legend[i][j] = new TLegend(0.71, 0.54, 0.899, 0.89);
      legend[i][j]->AddEntry(h_bnb[i][j][0],"Data (Beam-On)","lepf");
      legend[i][j]->AddEntry(h_overlay0[i][j][0],"Stat. Unc.","f");
      legend[i][j]->AddEntry(h_ext[i][j][0],"Data (Beam-Off)","f");
      for(int k =1; k < z; k++){	
	legend[i][j]->AddEntry(h_overlay[i][j][f-k],Form("%s",channel_legend[f-k]),"f");

      }
      legend[i][j]->SetLineWidth(0);
      //legend[i][j]->SetFillStyle(1);
      legend[i][j]->SetFillColor(kWhite);
      legend[i][j]->SetTextSize(0.03);
      legend[i][j]->Draw("same");
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s: %s}",cut_titles[i],titles[j]));
      t->DrawLatex(0.195,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
      canv[i][j]->cd();
      pad0[i][j] = new TPad(Form("pad0%s%s",plots[j],cut[i]),Form("pad0%s%s",plots[j],cut[i]),0,0.0,1.0,0.294);
      pad0[i][j]->SetTopMargin(0);                
      pad0[i][j]->SetBottomMargin(0.19);
      pad0[i][j]->SetGridx();
      pad0[i][j]->Draw();
      pad0[i][j]->cd();

      h_ext[i][j][1]->Add(h_overlay0[i][j][1]);
      h_bnb[i][j][1]->Divide(h_bnb[i][j][1],h_ext[i][j][1],1,1,"B");
      h_bnb[i][j][1]->Draw("e1p");
      h_bnb[i][j][1]->SetStats(kFALSE);
      h_bnb[i][j][1]->SetTitle("");
      
      TF1 *a = new TF1("a","1", -150000 , 150000);
      a->SetLineColor(kRed);
      a->Draw("SAME");

      //Calculate the Chi2
      double chisqv1 = calculatePearsonChiSq(h_bnb[i][j][0], h_overlay0[i][j][2]);
      int nBins1 = h_overlay0[i][j][2]->GetXaxis()->GetNbins();
      //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
      t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1/nBins1));
      
      h_bnb[i][j][1]->GetYaxis()->SetTitle("Ratio");
      h_bnb[i][j][1]->GetYaxis()->SetTitleSize(25);
      h_bnb[i][j][1]->GetYaxis()->SetTitleFont(43);
      h_bnb[i][j][1]->GetYaxis()->SetTitleOffset(1.5);
      h_bnb[i][j][1]->GetYaxis()->SetLabelFont(43);
      h_bnb[i][j][1]->GetYaxis()->SetLabelSize(25);
      h_bnb[i][j][1]->GetXaxis()->SetTitle(Form("%s",titles[j]));
      h_bnb[i][j][1]->GetXaxis()->SetTitleSize(25);
      h_bnb[i][j][1]->GetXaxis()->SetTitleFont(43);
      h_bnb[i][j][1]->GetXaxis()->SetTitleOffset(3);
      h_bnb[i][j][1]->GetXaxis()->SetLabelFont(43);
      h_bnb[i][j][1]->GetXaxis()->SetLabelSize(25);
      
      canv[i][j]->Print(Form("%s%s%s.png",path.c_str(),plots[j],cut[i]));
      canv[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),plots[j],cut[i]));
      
    }
  }

  //Plots of the Truth Variables
  ////////////////////////////////////
  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_truth; j++){
      
      canv_truth[i][j] = new TCanvas(Form("C_truth%s%s",truth[j],cut[i]),Form("C_truth%s%s",truth[j],cut[i]),2000,1500);
      h_truth[i][j] = new THStack(Form("h_truth%s%s",truth[j],cut[i]),Form("h_truth%s%s",truth[j],cut[i]));
      h_truth[i][j]->Draw("HIST");
      h_truth[i][j]->SetTitle(Form(" ; %s ; Number of Events",truth_titles[j]));
      h_truth[i][j]->SetMaximum(ylim_truth[i][j]);
      //h_truth[i][j]->GetXaxis()->SetTitle(Form("%s",truth_titles[j]));

      if(plot_ccnue_truth){
	z = num_channels;
	f = num_channels;
      }else{
	z = num_channels - 1;
	f = num_channels - 1;
      }
      
      for(int k=1; k < z ; k++){
	h_mc[i][j][k]->SetLineColor(colors[k]);
	h_mc[i][j][k]->SetFillColor(colors[k]);
	h_mc[i][j][k]->SetLineWidth(1);
	h_truth[i][j]->Add(h_mc[i][j][k]);
      }
      
      if(plot_total_truth){
      h_mc[i][j][0]->Draw("e1SAME");
      }

      legend_truth[i][j] = new TLegend(0.71, 0.54, 0.899, 0.89);
      for(int k =1; k < z; k++){	
	legend_truth[i][j]->AddEntry(h_mc[i][j][f-k],Form("%s",channel_legend[f-k]),"f");

      }
      legend_truth[i][j]->SetLineWidth(0);
      legend_truth[i][j]->SetFillColor(kWhite);
      legend_truth[i][j]->SetTextSize(0.03);
      legend_truth[i][j]->Draw("same");
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s: %s}",cut_titles[i],truth_titles[j]));
      t->DrawLatex(0.195,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
      canv_truth[i][j]->Print(Form("%s%s%s.png",path.c_str(),truth[j],cut[i]));
      canv_truth[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),truth[j],cut[i]));
      
    }
  }
  
  //Plots of the PFP Stuff:
  //////////////////////////
  for(int i = 0; i < num_group; i ++){
    
      canv_pfp[i] = new TCanvas(Form("C_pfp_%s",group[i]),Form("C_pfp_%s",group[i]),2000,1500);
      h_pfp[i] = new THStack(Form("h_pfp%s",group[i]),Form("h_pfp%s",group[i]));
      h_pfp[i]->Add(h_overlay_pfp[i]);
      h_pfp[i]->Add(h_ext_pfp[i]);
      h_pfp[i]->Draw("hist");
      h_pfp[i]->GetYaxis()->SetTitle("Number of Events");
      h_pfp[i]->GetXaxis()->SetTitle(Form("%s",titles_pfp[i]));
      h_pfp[i]->SetTitle(Form("%s",titles_pfp[i]));
      h_pfp[i]->SetMaximum(ylim_pfp[i]);

      h_ext_pfp[i]->SetFillColor(kRed);
      h_overlay_pfp[i]->SetFillColor(kBlue);

      h_bnb_pfp[i]->Draw("e1SAME");
      h_bnb_pfp[i]->SetLineColor(kBlack);
      h_bnb_pfp[i]->SetLineWidth(1);

      legend_pfp[i] = new TLegend(0.63,0.58,1.0,0.89);
      legend_pfp[i]->AddEntry(h_bnb_pfp[i], "BNB", "lepf");
      legend_pfp[i]->AddEntry(h_ext_pfp[i], "EXTBNB", "f");
      legend_pfp[i]->AddEntry(h_overlay_pfp[i], "Overlay MC", "f");
      legend_pfp[i]->Draw("same");
      
      canv_pfp[i]->Print(Form("%s_%s.png",path.c_str(),group[i]));
      canv_pfp[i]->Print(Form("%s_%s.pdf",path.c_str(),group[i]));
  }

  //Plot the 2D Histograms:
  //////////////////////////////////
  for(int i =0; i < num_group2d; i++){
    canv_2d[i] = new TCanvas(Form("C_2D_%s",group2d[i]),Form("C_2D_%s",group2d[i]),2000,1500);
    h_overlay2D[i]->Draw("colz"); 
    canv_2d[i]->Print(Form("%s_%s.png",path.c_str(),group2d[i]));
    canv_2d[i]->Print(Form("%s_%s.pdf",path.c_str(),group2d[i]));
  }


  //Plot the Efficiency Stuff:
  //////////////////////////
  for(int i=0; i < num_eff; i++){

    h_num0[i] = (TH1D*)h_num[i]->Clone();
    h_denom0[i] = (TH1D*)h_denom[i]->Clone();
    canv_both[i] = new TCanvas(Form("C_both_%s",eff[i]),Form("C_both_%s",eff[i]),2000,1500);
    h_denom0[i]->Draw("hist");
    h_denom0[i]->SetFillColor(kBlue);
    h_denom0[i]->SetTitle(Form("Efficiency(%s); %s ; Number of Entries",titles_eff[i],titles_eff[i]));
    h_denom0[i]->SetMaximum(ylim_eff[i]);
    h_num0[i]->Draw("SAME");
    h_num0[i]->SetFillColor(kRed);

    legend_eff[i] = new TLegend(0.63,0.58,1.0,0.89);
    legend_eff[i]->AddEntry(h_num0[i], "Numerator", "f");
    legend_eff[i]->AddEntry(h_denom0[i], "Denomiator", "f");
    legend_eff[i]->Draw("same");

    std::cout<<Form("Number of Entries in %s Denominator: ",titles_eff[i])<<h_denom0[i]->GetEntries()<<std::endl;
    std::cout<<Form("Number of Entries in %s  Numerator: ",titles_eff[i])<<h_num0[i]->GetEntries()<<std::endl;
    
    canv_both[i]->Print(Form("%s_%s_both.png",path.c_str(),eff[i]));
    canv_both[i]->Print(Form("%s_%s_both.pdf",path.c_str(),eff[i]));

    h_num1[i] = (TH1D*)h_num[i]->Clone();
    h_denom1[i] = (TH1D*)h_denom[i]->Clone();
    canv_eff[i] = new TCanvas(Form("C_eff_%s",eff[i]),Form("C_eff_%s",eff[i]),2000,1500);
    h_num1[i]->Divide(h_num1[i],h_denom1[i],1.0,1.0,"B");
    h_num1[i]->Draw("1e1p");
    h_num1[i]->SetTitle(Form("Efficiency(%s); %s ; Efficiency",titles_eff[i],titles_eff[i]));
    h_num1[i]->SetLineColor(kViolet);
    h_num1[i]->SetMaximum(1);
    h_num1[i]->SetMinimum(0);
    
    canv_eff[i]->Print(Form("%s_%s.png",path.c_str(),eff[i]));
    canv_eff[i]->Print(Form("%s_%s.pdf",path.c_str(),eff[i]));
  }

   //Chi2 Plots
  ////////////////////////////////////////
  for(int i = 0; i< num_planes; i++){
    for(int j = 0; j< num_hypothesis; j++){
      
      h_ext_chi2[i][j][1] = (TH1D*)h_ext_chi2[i][j][0]->Clone();
      h_ext_chi2[i][j][2] = (TH1D*)h_ext_chi2[i][j][0]->Clone();
      h_bnb_chi2[i][j][1] = (TH1D*)h_bnb_chi2[i][j][0]->Clone();
      h_overlay0_chi2[i][j][0] = (TH1D*)h_overlay_chi2[i][j][0]->Clone();
      h_overlay0_chi2[i][j][1] = (TH1D*)h_overlay_chi2[i][j][0]->Clone();
      h_overlay0_chi2[i][j][2] = (TH1D*)h_overlay_chi2[i][j][0]->Clone();
      
      canv_chi2[i][j] = new TCanvas(Form("C_chi2_%s_%s",plane[i],hypothesis[j]),Form("C_chi2_%s_%s",plane[i],hypothesis[j]),2000,1500);
      canv_chi2[i][j]->cd(1);
      h_chi2[i][j] = new THStack(Form("h_chi2_%s_%s",plane[i],hypothesis[j]),Form("h_chi2_%s_%s",plane[i],hypothesis[j]));
      pad_chi2[i][j] = new TPad(Form("pad_chi2_%s_%s",plane[i],hypothesis[j]),Form("pad_chi2_%s_%s",plane[i],hypothesis[j]),0,0.28,1.0,1.0);
      pad_chi2[i][j]->SetBottomMargin(0.0);
      pad_chi2[i][j]->SetGridx();
      pad_chi2[i][j]->SetBorderMode(0);
      pad_chi2[i][j]->Draw();
      pad_chi2[i][j]->cd();

      //Stacked Histrogram parameters
      h_chi2[i][j]->Draw("HIST");
      h_chi2[i][j]->SetTitle("");
      h_chi2[i][j]->SetMaximum(ylim_chi2[i][j]); 
      
      for(int k=1; k < num_particles ; k++){
	h_overlay_chi2[i][j][k]->SetLineColor(colors_chi2[k]);
	h_overlay_chi2[i][j][k]->SetFillColor(colors_chi2[k]);
	h_overlay_chi2[i][j][k]->SetLineWidth(1);
	h_chi2[i][j]->Add(h_overlay_chi2[i][j][k]);
      }

      //EXT
      h_chi2[i][j]->Add(h_ext_chi2[i][j][0]);
      h_ext_chi2[i][j][0]->SetFillColor(kViolet-7);
      h_ext_chi2[i][j][0]->SetFillStyle(3005);
      h_ext_chi2[i][j][0]->SetLineWidth(1);
      
      //BNB
      h_bnb_chi2[i][j][0]->Draw("e1SAME");
      h_bnb_chi2[i][j][0]->SetLineColor(kBlack);
      h_bnb_chi2[i][j][0]->SetLineWidth(1);
      
      //if you want to plot the total for sanity sake:
      if(plot_total_chi2){
      h_overlay_chi2[i][j][0]->Draw("SAME");
      }
      
      //Make sure to do overlay statistical uncertainty
      OverlayStatistics(h_overlay0_chi2[i][j][0],h_ext_chi2[i][j][2]);
      h_overlay0_chi2[i][j][0]->Draw("e2SAME");
      h_overlay0_chi2[i][j][0]->SetLineColor(kBlack);
      h_overlay0_chi2[i][j][0]->SetFillColor(kBlack);
      h_overlay0_chi2[i][j][0]->SetFillStyle(3004);
      h_overlay0_chi2[i][j][0]->SetMarkerSize(0);
      h_overlay0_chi2[i][j][0]->SetLineWidth(1);

      legend_chi2[i][j] = new TLegend(0.71, 0.54, 0.899, 0.89);
      legend_chi2[i][j]->AddEntry(h_bnb_chi2[i][j][0],"Data (Beam-On)","lepf");
      legend_chi2[i][j]->AddEntry(h_overlay0_chi2[i][j][0],"Stat. Unc.","f");
      legend_chi2[i][j]->AddEntry(h_ext_chi2[i][j][0],"Data (Beam-Off)","f");
      for(int k =1; k < num_particles; k++){	
	legend_chi2[i][j]->AddEntry(h_overlay_chi2[i][j][num_particles-k],Form("%s",channel_legend_chi2[num_particles-k]),"f");

      }
      legend_chi2[i][j]->SetLineWidth(0);
      //legend_chi2[i][j]->SetFillStyle(1);
      legend_chi2[i][j]->SetFillColor(kWhite);
      legend_chi2[i][j]->SetTextSize(0.03);
      legend_chi2[i][j]->Draw("same");
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s: %s}",titles_planes[i],titles_chi2[j]));
      t->DrawLatex(0.195,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
      canv_chi2[i][j]->cd();
      pad0_chi2[i][j] = new TPad(Form("pad0_chi2_%s_%s",plane[i],hypothesis[j]),Form("pad0_%s_%s",plane[i],hypothesis[j]),0,0.0,1.0,0.294);
      pad0_chi2[i][j]->SetTopMargin(0);                
      pad0_chi2[i][j]->SetBottomMargin(0.19);
      pad0_chi2[i][j]->SetGridx();
      pad0_chi2[i][j]->Draw();
      pad0_chi2[i][j]->cd();

      h_ext_chi2[i][j][1]->Add(h_overlay0_chi2[i][j][1]);
      h_bnb_chi2[i][j][1]->Divide(h_bnb_chi2[i][j][1],h_ext_chi2[i][j][1],1,1,"B");
      h_bnb_chi2[i][j][1]->Draw("e1p");
      h_bnb_chi2[i][j][1]->SetStats(kFALSE);
      h_bnb_chi2[i][j][1]->SetTitle("");
      
      TF1 *a_chi2 = new TF1("a_chi2","1", -150000 , 150000);
      a_chi2->SetLineColor(kRed);
      a_chi2->Draw("SAME");

      //Calculate the Chi2
      double chisqv1_chi2 = calculatePearsonChiSq(h_bnb_chi2[i][j][0], h_overlay0_chi2[i][j][2]);
      int nBins1_chi2 = h_overlay0_chi2[i][j][2]->GetXaxis()->GetNbins();
      //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv_chi2, nBins_chi2, chisqv/nBins_chi2);
      t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_chi2/nBins1_chi2));
      
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitle("Ratio");
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitleSize(25);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitleFont(43);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitleOffset(1.5);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetLabelFont(43);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetLabelSize(25);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitle(Form("%s",titles_chi2[j]));
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitleSize(25);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitleFont(43);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitleOffset(3);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetLabelFont(43);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetLabelSize(25);
      
      canv_chi2[i][j]->Print(Form("%s_%s_%s.png",path.c_str(),plane[i],hypothesis[j]));
      canv_chi2[i][j]->Print(Form("%s_%s_%s.pdf",path.c_str(),plane[i],hypothesis[j]));
      
    }
  }
  
} //end of program
