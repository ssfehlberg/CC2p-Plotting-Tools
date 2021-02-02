#define analysis_cxx
#include "analysis.h"

#include "paul_tol_colors.hpp"
#include <iostream>
#include <ctime>
#include <string>


void analysis::main(){

  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);

  ///////////////////////////////////////////////////////////
  //FIRST: CHECK WHICH SAMPLE THIS IS: FILTERED OR UNFILTERED
  //////////////////////////////////////////////////////////
  const char* sample= which_sample();

  ///////////////////////////////////////////////////////
  //SECOND: DEFINE HISTOGRAMS FILES AND GENERAL VARIABLES
  //////////////////////////////////////////////////////
  TFile *f1=new TFile(Form("../pelee/histograms_%s_wgt.root",sample));//overlay histograms. Note: be sure to use weighted for both dirt and overlay
  TFile *f_dirt=new TFile(Form("../pelee/histograms_%s_dirt_wgt.root",sample));//dirt histograms 
  TFile *f2=new TFile(Form("../pelee/histograms_%s_bnb.root",sample));//bnb histograms  
  TFile *f3=new TFile(Form("../pelee/histograms_%s_ext.root",sample));//extbnb histograms 
  TFile *f4=new TFile("../histograms_efficiency.root"); //efficiency histograms
  
  //Color Scheme
  tolcols::init();
  Color_t colors[] = {0,9031, 9030, 9029, 9028, 9026, 9025, 9024, 9032, kGray+2, 9027};  //Black, light pink, light orange, light yellow, olive, mint, light cyan, light blue, light grey, darker grey, olive
  Color_t colors_raquel[] = {0,9012,9011,9010,9009,9008,9007,kGray+2,9032,9027}; //black, magenta, red, orange, mint, cyan, blue, dark gray, light gray, ccnue
  
  //POT number and In Progress
  char const * pot_num="#scale[0.6]{Accumulated POT: 4.566e+19}";//pot number printed on the plots
  char const * sample_name="#scale[0.6]{MicroBooNE In-Progress}";//sample name printed on the plots

  //latex stuff
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  
  //Stuff for time/date
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;

  /////////////////////////////////////////////////////////////
  //THIRD: MAKE DIRECTORY WITH TODAY'S DATE TO STORE ALL IMAGES
  //////////////////////////////////////////////////////////////
  const char* pathname = Form("images/%d%d%d_%s/",Month,Day,Year,sample);
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

  /////////////////////////////////////////
  //GRAB ALL THE HISTOGRAMS FROM THE FILES
  ////////////////////////////////////////
  Grab_Histograms(f1, f2, f3, f4, f_dirt);
  
  /////////////////////////////////////////////
  //PLOTTING TIME
  ///////////////////////////////////////////
  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_variables; j++){

      for(int k=0; k < num_channels; k++){
	h_overlay_vec.push_back(h_overlay[i][j][k]);
      }
      
      Plot_Histograms(colors, h_overlay_vec, h_overlay0[i][j][0],h_overlay0[i][j][1],h_overlay0[i][j][2], h_ext[i][j][0], h_ext[i][j][1], h_ext[i][j][2], h_dirt[i][j][0],h_dirt[i][j][1],h_dirt[i][j][2],h_bnb[i][j][0],h_bnb[i][j][1],
		      canv[i][j], h[i][j], pad[i][j], pad0[i][j], legend[i][j],ylim[i][j],ymin[i][j], num_channels, titles[j], path, plots[j], cut[i], false, false);

      h_overlay_vec.clear();
      
    }
  }
  
  /* 
 //Now to make some plots: Reconstructed MC,EXT,and BNB
  ////////////////////////////////////////
  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_variables; j++){

      h_ext[i][j][1] = (TH1D*)h_ext[i][j][0]->Clone();
      h_ext[i][j][2] = (TH1D*)h_ext[i][j][0]->Clone();
      h_dirt[i][j][1] = (TH1D*)h_dirt[i][j][0]->Clone();
      h_dirt[i][j][2] = (TH1D*)h_dirt[i][j][0]->Clone();
      h_bnb[i][j][1] = (TH1D*)h_bnb[i][j][0]->Clone();
      h_overlay0[i][j][0] = (TH1D*)h_overlay[i][j][0]->Clone();
      h_overlay0[i][j][1] = (TH1D*)h_overlay[i][j][0]->Clone();
      h_overlay0[i][j][2] = (TH1D*)h_overlay[i][j][0]->Clone();
      
      canv[i][j] = new TCanvas(Form("C%s%s",plots[j],cut[i]),Form("C%s%s",plots[j],cut[i]),2000,1500);
      canv[i][j]->cd(1);
      h[i][j] = new THStack(Form("h%s%s",plots[j],cut[i]),Form("h%s%s",plots[j],cut[i]));
      pad[i][j] = new TPad(Form("pad%s%s",plots[j],cut[i]),Form("pad%s%s",plots[j],cut[i]),0,0.35,1.0,1.0);
      pad[i][j]->SetBottomMargin(0.0);
      pad[i][j]->SetGridx();
      pad[i][j]->SetBorderMode(0);
      pad[i][j]->Draw();
      pad[i][j]->cd();

      //Stacked Histrogram parameters
      h[i][j]->Draw("HIST");
      h[i][j]->SetTitle("");
      h[i][j]->SetMaximum(ylim[i][j]); 
      h[i][j]->SetMinimum(ymin[i][j]);
   
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

      //Dirt
      h[i][j]->Add(h_dirt[i][j][0]);
      h_dirt[i][j][0]->SetFillColor(kOrange-8);
      h_dirt[i][j][0]->SetLineColor(kOrange-8);
      h_dirt[i][j][0]->SetLineWidth(1);
      
      //EXT
      h[i][j]->Add(h_ext[i][j][0]);
      h_ext[i][j][0]->SetFillColor(kViolet-7);
      h_ext[i][j][0]->SetFillStyle(3005);
      h_ext[i][j][0]->SetLineWidth(1);
      
      //BNB
      h_bnb[i][j][0]->Draw("e1SAME");
      h_bnb[i][j][0]->SetLineColor(kBlack);
      h_bnb[i][j][0]->SetLineWidth(1);

      h[i][j]->GetYaxis()->SetTitle("No. Events");
   
      //if you want to plot the total for sanity sake:
      if(plot_total){
      h_overlay[i][j][0]->Draw("SAME");
      }
      
      //Make sure to do overlay statistical uncertainty
      OverlayStatistics(h_overlay0[i][j][0],h_ext[i][j][2], h_dirt[i][j][2]);
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
      legend[i][j]->AddEntry(h_dirt[i][j][0],"Dirt","f");
      for(int k =1; k < z; k++){	
	legend[i][j]->AddEntry(h_overlay[i][j][f-k],Form("%s",channel_legend[f-k]),"f");

      }
      legend[i][j]->SetLineWidth(0);
      //legend[i][j]->SetFillStyle(1);
      legend[i][j]->SetFillColor(kWhite);
      legend[i][j]->SetTextSize(0.03);
      legend[i][j]->Draw("same");
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s}",titles[j]));
      t->DrawLatex(0.195,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
      canv[i][j]->cd();
      pad0[i][j] = new TPad(Form("pad0%s%s",plots[j],cut[i]),Form("pad0%s%s",plots[j],cut[i]),0,0.0,1.0,0.35);
      pad0[i][j]->SetTopMargin(0);                
      pad0[i][j]->SetBottomMargin(0.19);
      pad0[i][j]->SetGridx();
      pad0[i][j]->Draw();
      pad0[i][j]->cd();

      h_ext[i][j][1]->Add(h_overlay0[i][j][1]);
      h_ext[i][j][1]->Add(h_dirt[i][j][1]);
      h_ext[i][j][1]->Sumw2();
      h_bnb[i][j][1]->Divide(h_bnb[i][j][1],h_ext[i][j][1],1,1);

      // for(int p=0; p < h_bnb[i][j][1]->GetNBins(); p++){a
      //	h_ext = h_ext[i][j][1]->GetBinContent(p);
      //	h_bnb = h_bnb[i][j][1]->GetBinContent(p);
      //	h_bnb[i][j][1]->SetBinError(p, std::sqrt(h_bnb/std::pow(h_ext,2) + (std::pow(h_bnb,2)/std::pow(h_ext,3))));
	
      //      }

      h_bnb[i][j][1]->Sumw2();
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
      t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1/nBins1));
      
      h_bnb[i][j][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
      h_bnb[i][j][1]->GetYaxis()->CenterTitle();
      h_bnb[i][j][1]->GetYaxis()->SetTitleSize(28);
      h_bnb[i][j][1]->GetYaxis()->SetTitleFont(43);
      h_bnb[i][j][1]->GetYaxis()->SetTitleOffset(1.5);
      h_bnb[i][j][1]->GetYaxis()->SetLabelFont(43);
      h_bnb[i][j][1]->GetYaxis()->SetLabelSize(30);
      h_bnb[i][j][1]->GetXaxis()->SetTitle(Form("%s",titles[j]));
      h_bnb[i][j][1]->GetXaxis()->SetTitleSize(35);
      h_bnb[i][j][1]->GetXaxis()->SetTitleFont(43);
      h_bnb[i][j][1]->GetXaxis()->SetTitleOffset(3);
      h_bnb[i][j][1]->GetXaxis()->SetLabelFont(43);
      h_bnb[i][j][1]->GetXaxis()->SetLabelSize(35);
      h_bnb[i][j][1]->SetMaximum(2.7);
      h_bnb[i][j][1]->SetMinimum(0);
      
      canv[i][j]->Print(Form("%s%s%s.png",path.c_str(),plots[j],cut[i]));
      canv[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),plots[j],cut[i]));
      
    }
  }
  */
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
    h_num1[i]->Divide(h_num1[i],h_denom1[i],1.0,1.0, "B");
    h_num1[i]->Draw("1e1p");
    h_num1[i]->SetTitle(Form(" ; %s ; Efficiency",titles_eff[i]));
    h_num1[i]->SetLineColor(kViolet);
    h_num1[i]->SetMaximum(1);
    h_num1[i]->SetMinimum(0);
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s}",titles_eff[i]));
    t->DrawLatex(0.23,0.92,"#scale[0.5]{Accumulated POT: 4.566e+19}");
    t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
    a[i] = new TLine(xlim_eff[i],0,xlim_eff[i],1);
    a[i]->Draw("same");
    a[i]->SetLineColor(kBlack);
    a[i]->SetLineWidth(4);
    canv_eff[i]->Print(Form("%s_%s_eff.png",path.c_str(),eff[i]));
    canv_eff[i]->Print(Form("%s_%s_eff.pdf",path.c_str(),eff[i]));
  }

  //Chi2 Plots
  ////////////////////////////////////////
  for(int i = 0; i< num_planes; i++){
    for(int j = 0; j< num_hypothesis; j++){
      
      h_ext_chi2[i][j][1] = (TH1D*)h_ext_chi2[i][j][0]->Clone();
      h_ext_chi2[i][j][2] = (TH1D*)h_ext_chi2[i][j][0]->Clone();
      h_dirt_chi2[i][j][1] = (TH1D*)h_dirt_chi2[i][j][0]->Clone();
      h_dirt_chi2[i][j][2] = (TH1D*)h_dirt_chi2[i][j][0]->Clone();
      h_bnb_chi2[i][j][1] = (TH1D*)h_bnb_chi2[i][j][0]->Clone();
      h_overlay0_chi2[i][j][0] = (TH1D*)h_overlay_chi2[i][j][0]->Clone();
      h_overlay0_chi2[i][j][1] = (TH1D*)h_overlay_chi2[i][j][0]->Clone();
      h_overlay0_chi2[i][j][2] = (TH1D*)h_overlay_chi2[i][j][0]->Clone();
      
      canv_chi2[i][j] = new TCanvas(Form("C_chi2_%s_%s",plane[i],hypothesis[j]),Form("C_chi2_%s_%s",plane[i],hypothesis[j]),2000,1500);
      canv_chi2[i][j]->cd(1);
      h_chi2[i][j] = new THStack(Form("h_chi2_%s_%s",plane[i],hypothesis[j]),Form("h_chi2_%s_%s",plane[i],hypothesis[j]));
      pad_chi2[i][j] = new TPad(Form("pad_chi2_%s_%s",plane[i],hypothesis[j]),Form("pad_chi2_%s_%s",plane[i],hypothesis[j]),0,0.35,1.0,1.0);
      pad_chi2[i][j]->SetBottomMargin(0.0);
      pad_chi2[i][j]->SetGridx();
      pad_chi2[i][j]->SetBorderMode(0);
      pad_chi2[i][j]->Draw();
      pad_chi2[i][j]->cd();

      //Stacked Histrogram parameters
      h_chi2[i][j]->Draw("HIST");
      h_chi2[i][j]->SetTitle("");
      h_chi2[i][j]->SetMaximum(ylim_chi2[i][j]);
      h_chi2[i][j]->SetMinimum(-2);
      
      for(int k=1; k < num_particles ; k++){
	h_overlay_chi2[i][j][k]->SetLineColor(colors_chi2[k]);
	h_overlay_chi2[i][j][k]->SetFillColor(colors_chi2[k]);
	h_overlay_chi2[i][j][k]->SetLineWidth(1);
	h_chi2[i][j]->Add(h_overlay_chi2[i][j][k]);
      }
      
      //Dirt
      h_chi2[i][j]->Add(h_dirt_chi2[i][j][0]);
      h_dirt_chi2[i][j][0]->SetFillColor(kOrange-8);
      h_dirt_chi2[i][j][0]->SetLineColor(kOrange-8);
      h_dirt_chi2[i][j][0]->SetLineWidth(1);

	
      //EXT
      h_chi2[i][j]->Add(h_ext_chi2[i][j][0]);
      h_ext_chi2[i][j][0]->SetFillColor(kViolet-7);
      h_ext_chi2[i][j][0]->SetFillStyle(3005);
      h_ext_chi2[i][j][0]->SetLineWidth(1);
	
      //BNB
      h_bnb_chi2[i][j][0]->Draw("e1SAME");
      h_bnb_chi2[i][j][0]->SetLineColor(kBlack);
      h_bnb_chi2[i][j][0]->SetLineWidth(1);

      h_chi2[i][j]->GetYaxis()->SetTitle("No. Tracks");
   
      
      //if you want to plot the total for sanity sake:
      if(plot_total_chi2){
	h_overlay_chi2[i][j][0]->Draw("SAME");
      }
      
      //Make sure to do overlay statistical uncertainty
      OverlayStatistics(h_overlay0_chi2[i][j][0],h_ext_chi2[i][j][2], h_dirt_chi2[i][j][2]);
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
      legend_chi2[i][j]->AddEntry(h_dirt_chi2[i][j][0],"Dirt","f");
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
      pad0_chi2[i][j] = new TPad(Form("pad0_chi2_%s_%s",plane[i],hypothesis[j]),Form("pad0_%s_%s",plane[i],hypothesis[j]),0,0.0,1.0,0.35);
      pad0_chi2[i][j]->SetTopMargin(0);                
      pad0_chi2[i][j]->SetBottomMargin(0.19);
      pad0_chi2[i][j]->SetGridx();
      pad0_chi2[i][j]->Draw();
      pad0_chi2[i][j]->cd();
      
      h_ext_chi2[i][j][1]->Add(h_overlay0_chi2[i][j][1]);
      h_ext_chi2[i][j][1]->Add(h_dirt_chi2[i][j][1]);
      h_ext_chi2[i][j][1]->Sumw2();
      h_bnb_chi2[i][j][1]->Divide(h_bnb_chi2[i][j][1],h_ext_chi2[i][j][1],1,1);
      h_bnb_chi2[i][j][1]->Sumw2();
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
      t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_chi2/nBins1_chi2));
	
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
      h_bnb_chi2[i][j][1]->GetYaxis()->CenterTitle();
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitleSize(28);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitleFont(43);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetTitleOffset(1.5);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetLabelFont(43);
      h_bnb_chi2[i][j][1]->GetYaxis()->SetLabelSize(30);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitle(Form("%s",titles_chi2[j]));
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitleSize(35);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitleFont(43);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetTitleOffset(3);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetLabelFont(43);
      h_bnb_chi2[i][j][1]->GetXaxis()->SetLabelSize(35);
      
      canv_chi2[i][j]->Print(Form("%s_%s_%s.png",path.c_str(),plane[i],hypothesis[j]));
      canv_chi2[i][j]->Print(Form("%s_%s_%s.pdf",path.c_str(),plane[i],hypothesis[j]));
      
    }
  }
  
  //3D Chi2 Plots
  for(int i = 0; i < num_cuts_3D; i++){
    for(int j = 0; j< num_hypothesis_3D; j++){
    
      h_ext_chi2_3D[i][j][1] = (TH1D*)h_ext_chi2_3D[i][j][0]->Clone();
      h_ext_chi2_3D[i][j][2] = (TH1D*)h_ext_chi2_3D[i][j][0]->Clone();
      h_dirt_chi2_3D[i][j][1] = (TH1D*)h_dirt_chi2_3D[i][j][0]->Clone();
      h_dirt_chi2_3D[i][j][2] = (TH1D*)h_dirt_chi2_3D[i][j][0]->Clone();
      h_bnb_chi2_3D[i][j][1] = (TH1D*)h_bnb_chi2_3D[i][j][0]->Clone();
      h_overlay0_chi2_3D[i][j][0] = (TH1D*)h_overlay_chi2_3D[i][j][0]->Clone();
      h_overlay0_chi2_3D[i][j][1] = (TH1D*)h_overlay_chi2_3D[i][j][0]->Clone();
      h_overlay0_chi2_3D[i][j][2] = (TH1D*)h_overlay_chi2_3D[i][j][0]->Clone();
      
      canv_chi2_3D[i][j] = new TCanvas(Form("C_chi2_3D_%s%s",hypothesis_3D[j],cuts_3D[i]),Form("C_chi2_3D_%s%s",hypothesis_3D[j],cuts_3D[i]),2000,1500);
      canv_chi2_3D[i][j]->cd(1);
      h_chi2_3D[i][j] = new THStack(Form("h_chi2_3D_%s%s",hypothesis_3D[j],cuts_3D[i]),Form("h_chi2__3D_%s%s",hypothesis_3D[j],cuts_3D[i]));
      pad_chi2_3D[i][j] = new TPad(Form("pad_chi2_3D_%s%s",hypothesis_3D[j],cuts_3D[i]),Form("pad_chi2__3D_%s%s",hypothesis_3D[j],cuts_3D[i]),0,0.35,1.0,1.0);
      pad_chi2_3D[i][j]->SetBottomMargin(0.01);
      pad_chi2_3D[i][j]->SetGridx();
      pad_chi2_3D[i][j]->SetBorderMode(0);
      pad_chi2_3D[i][j]->Draw();
      pad_chi2_3D[i][j]->cd();
  
      //Stacked Histrogram parameters
      h_chi2_3D[i][j]->Draw("HIST");
      h_chi2_3D[i][j]->SetTitle("");
      h_chi2_3D[i][j]->SetMaximum(ylim_chi2_3D[i][j]);
      h_chi2_3D[i][j]->SetMinimum(-2);
      
      for(int k=1; k < num_particles ; k++){
	h_overlay_chi2_3D[i][j][k]->SetLineColor(colors_chi2[k]);
	h_overlay_chi2_3D[i][j][k]->SetFillColor(colors_chi2[k]);
	h_overlay_chi2_3D[i][j][k]->SetLineWidth(1);
	h_chi2_3D[i][j]->Add(h_overlay_chi2_3D[i][j][k]);
      }
      
      //Dirt
      h_chi2_3D[i][j]->Add(h_dirt_chi2_3D[i][j][0]);
      h_dirt_chi2_3D[i][j][0]->SetFillColor(kOrange-8);
      h_dirt_chi2_3D[i][j][0]->SetLineColor(kOrange-8);
      h_dirt_chi2_3D[i][j][0]->SetLineWidth(1);

      //EXT
      h_chi2_3D[i][j]->Add(h_ext_chi2_3D[i][j][0]);
      h_ext_chi2_3D[i][j][0]->SetFillColor(kViolet-7);
      h_ext_chi2_3D[i][j][0]->SetFillStyle(3005);
      h_ext_chi2_3D[i][j][0]->SetLineWidth(1);
      
      //BNB
      h_bnb_chi2_3D[i][j][0]->Draw("e1SAME");
      h_bnb_chi2_3D[i][j][0]->SetLineColor(kBlack);
      h_bnb_chi2_3D[i][j][0]->SetLineWidth(1);

      h_chi2_3D[i][j]->GetYaxis()->SetTitle("No. Tracks");
      
      //if you want to plot the total for sanity sake:
      if(plot_total_chi2_3D){
	h_overlay_chi2_3D[i][j][0]->Draw("SAME");
      }

      a_3D[i][j] = new TLine(xlim_chi2_3D[i][j],0,xlim_chi2_3D[i][j],ylim_chi2_3D[i][j]+30);
      a_3D[i][j]->Draw("same");
      a_3D[i][j]->SetLineColor(kBlack);
      a_3D[i][j]->SetLineWidth(4);
      
      //Make sure to do overlay statistical uncertainty
      OverlayStatistics(h_overlay0_chi2_3D[i][j][0],h_ext_chi2_3D[i][j][2], h_dirt_chi2_3D[i][j][2]);
      h_overlay0_chi2_3D[i][j][0]->Draw("e2SAME");
      h_overlay0_chi2_3D[i][j][0]->SetLineColor(kBlack);
      h_overlay0_chi2_3D[i][j][0]->SetFillColor(kBlack);
      h_overlay0_chi2_3D[i][j][0]->SetFillStyle(3004);
      h_overlay0_chi2_3D[i][j][0]->SetMarkerSize(0);
      h_overlay0_chi2_3D[i][j][0]->SetLineWidth(1);

      legend_chi2_3D[i][j] = new TLegend(0.71, 0.54, 0.899, 0.89);
      legend_chi2_3D[i][j]->AddEntry(h_bnb_chi2_3D[i][j][0],"Data (Beam-On)","lepf");
      legend_chi2_3D[i][j]->AddEntry(h_overlay0_chi2_3D[i][j][0],"Stat. Unc.","f");
      legend_chi2_3D[i][j]->AddEntry(h_ext_chi2_3D[i][j][0],"Data (Beam-Off)","f");
      legend_chi2_3D[i][j]->AddEntry(h_dirt_chi2_3D[i][j][0],"Dirt","f");
      for(int k =1; k < num_particles; k++){
	legend_chi2_3D[i][j]->AddEntry(h_overlay_chi2_3D[i][j][num_particles-k],Form("%s",channel_legend_chi2[num_particles-k]),"f");
      }
      legend_chi2_3D[i][j]->SetLineWidth(0);
      //legend_chi2_3D[i][j]->SetFillStyle(1);
      legend_chi2_3D[i][j]->SetFillColor(kWhite);
      legend_chi2_3D[i][j]->SetTextSize(0.03);
      legend_chi2_3D[i][j]->Draw("same");
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s}",titles_chi2_3D[j])); //,cuts_3D_titles[i]
      t->DrawLatex(0.195,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));

      canv_chi2_3D[i][j]->cd();
      pad0_chi2_3D[i][j] = new TPad(Form("pad0_3D_%s%s",hypothesis_3D[j],cuts_3D[i]),Form("pad0_3D_%s%s",hypothesis_3D[j],cuts_3D[i]),0,0.0,1.0,0.35);
      pad0_chi2_3D[i][j]->SetTopMargin(0);                
      pad0_chi2_3D[i][j]->SetBottomMargin(0.15);
      pad0_chi2_3D[i][j]->SetGridx();
      pad0_chi2_3D[i][j]->Draw();
      pad0_chi2_3D[i][j]->cd();

      h_ext_chi2_3D[i][j][1]->Add(h_overlay0_chi2_3D[i][j][1]);
      h_ext_chi2_3D[i][j][1]->Add(h_dirt_chi2_3D[i][j][1]);
      h_ext_chi2_3D[i][j][1]->Sumw2();
      h_bnb_chi2_3D[i][j][1]->Divide(h_bnb_chi2_3D[i][j][1],h_ext_chi2_3D[i][j][1],1,1);
      h_bnb_chi2_3D[i][j][1]->Sumw2();
      h_bnb_chi2_3D[i][j][1]->Draw("e1p");
      h_bnb_chi2_3D[i][j][1]->SetStats(kFALSE);
      h_bnb_chi2_3D[i][j][1]->SetTitle("");
      h_bnb_chi2_3D[i][j][1]->SetMaximum(2.7); //just to make 3D chi2p look nice
      
      TF1 *a_chi2_3D = new TF1("a_chi2_3D","1", -150000 , 150000);
      a_chi2_3D->SetLineColor(kRed);
      a_chi2_3D->Draw("SAME");

      //Calculate the Chi2
      double chisqv1_chi2 = calculatePearsonChiSq(h_bnb_chi2_3D[i][j][0], h_overlay0_chi2_3D[i][j][2]);
      int nBins1_chi2 = h_overlay0_chi2_3D[i][j][2]->GetXaxis()->GetNbins();
      //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv_chi2, nBins_chi2, chisqv/nBins_chi2);
      t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_chi2/nBins1_chi2));
      
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->CenterTitle();
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->SetTitleSize(28);
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->SetTitleFont(43);
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->SetTitleOffset(1.5);
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->SetLabelFont(43);
      h_bnb_chi2_3D[i][j][1]->GetYaxis()->SetLabelSize(30);
      h_bnb_chi2_3D[i][j][1]->GetXaxis()->SetTitle(Form("%s",titles_chi2_3D[j]));
      h_bnb_chi2_3D[i][j][1]->GetXaxis()->SetTitleSize(35);
      h_bnb_chi2_3D[i][j][1]->GetXaxis()->SetTitleFont(43);
      h_bnb_chi2_3D[i][j][1]->GetXaxis()->SetTitleOffset(3);
      h_bnb_chi2_3D[i][j][1]->GetXaxis()->SetLabelFont(43);
      h_bnb_chi2_3D[i][j][1]->GetXaxis()->SetLabelSize(35);
      
      canv_chi2_3D[i][j]->Print(Form("%s_3D_%s%s.png",path.c_str(),hypothesis_3D[j],cuts_3D[i]));
      canv_chi2_3D[i][j]->Print(Form("%s_3D_%s%s.pdf",path.c_str(),hypothesis_3D[j],cuts_3D[i]));
      
    } //end of number of hypothesis
  }
  

  ////////////////
  //Particle specific plots
  ////////////////
  for(int i = 0; i < num_var; i++){

    //cloning histograms for statistics
    h_muon_ext[i][1] = (TH1D*)h_muon_ext[i][0]->Clone();
    h_muon_ext[i][2] = (TH1D*)h_muon_ext[i][0]->Clone();
    h_muon_dirt[i][1] = (TH1D*)h_muon_dirt[i][0]->Clone();
    h_muon_dirt[i][2] = (TH1D*)h_muon_dirt[i][0]->Clone();   
    h_muon_bnb[i][1] = (TH1D*)h_muon_bnb[i][0]->Clone();
    h_muon_overlay0[i][0] = (TH1D*)h_muon_overlay[i][0]->Clone();
    h_muon_overlay0[i][1] = (TH1D*)h_muon_overlay[i][0]->Clone();
    h_muon_overlay0[i][2] = (TH1D*)h_muon_overlay[i][0]->Clone();
      
    canv_muon[i] = new TCanvas(Form("C_muon%s",var[i]),Form("C_muon%s",var[i]),2000,1500);
    canv_muon[i]->cd(1);
    h_muon[i] = new THStack(Form("h_muon%s",var[i]),Form("h_muon%s",var[i]));
    pad_muon[i] = new TPad(Form("pad_muon%s",var[i]),Form("pad_muon%s",var[i]),0,0.35,1.0,1.0);
    pad_muon[i]->SetBottomMargin(0.0);
    pad_muon[i]->SetGridx();
    pad_muon[i]->SetBorderMode(0);
    pad_muon[i]->Draw();
    pad_muon[i]->cd();

    //Stacked Histrogram parameters
    h_muon[i]->Draw("HIST");
    h_muon[i]->SetTitle("");
    h_muon[i]->SetMaximum(muon_ylim[i]);
    h_muon[i]->SetMinimum(-2);

    //MC
    if(plot_ccnue){
      z = num_channels;
      f = 11;
    }else{
      z = num_channels -1;
      f = 10;
    }
      
    for(int k=1; k < z ; k++){
      h_muon_overlay[i][k]->SetLineColor(colors[k]);
      h_muon_overlay[i][k]->SetFillColor(colors[k]);
      h_muon_overlay[i][k]->SetLineWidth(1);
      h_muon[i]->Add(h_muon_overlay[i][k]);
    }

  
    //Dirt
    h_muon[i]->Add(h_muon_dirt[i][0]);
    h_muon_dirt[i][0]->SetFillColor(kOrange-8);
    h_muon_dirt[i][0]->SetLineColor(kOrange-8);
    h_muon_dirt[i][0]->SetLineWidth(1);
    
    //EXT
    h_muon[i]->Add(h_muon_ext[i][0]);
    h_muon_ext[i][0]->SetFillColor(kViolet-7);
    h_muon_ext[i][0]->SetFillStyle(3005);
    h_muon_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_muon_bnb[i][0]->Draw("e1SAME");
    h_muon_bnb[i][0]->SetLineColor(kBlack);
    h_muon_bnb[i][0]->SetLineWidth(1);

    h_muon[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total){
      h_muon_overlay[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_muon_overlay0[i][0],h_muon_ext[i][2], h_muon_dirt[i][2]);
    h_muon_overlay0[i][0]->Draw("e2SAME");
    h_muon_overlay0[i][0]->SetLineColor(kBlack);
    h_muon_overlay0[i][0]->SetFillColor(kBlack);
    h_muon_overlay0[i][0]->SetFillStyle(3004);
    h_muon_overlay0[i][0]->SetMarkerSize(0);
    h_muon_overlay0[i][0]->SetLineWidth(1);

    legend_muon[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_muon[i]->AddEntry(h_muon_bnb[i][0],"Data (Beam-On)","lepf");
    legend_muon[i]->AddEntry(h_muon_overlay0[i][0],"Stat. Unc.","f");
    legend_muon[i]->AddEntry(h_muon_ext[i][0],"Data (Beam-Off)","f");
    legend_muon[i]->AddEntry(h_muon_dirt[i][0],"Dirt","f");
    for(int k =1; k < z; k++){	
      legend_muon[i]->AddEntry(h_muon_overlay[i][f-k],Form("%s",channel_legend[f-k]),"f");  
    }
    legend_muon[i]->SetLineWidth(0);
    //legend_muon[i]->SetFillStyle(1);
    legend_muon[i]->SetFillColor(kWhite);
    legend_muon[i]->SetTextSize(0.03);
    legend_muon[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Muon: %s}",titles_var[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_muon[i]->cd();
    pad0_muon[i] = new TPad(Form("pad0_muon%s",var[i]),Form("pad0_muon%s",var[i]),0,0.0,1.0,0.35);
    pad0_muon[i]->SetTopMargin(0);                
    pad0_muon[i]->SetBottomMargin(0.19);
    pad0_muon[i]->SetGridx();
    pad0_muon[i]->Draw();
    pad0_muon[i]->cd();

    h_muon_ext[i][1]->Add(h_muon_overlay0[i][1]);
    h_muon_ext[i][1]->Add(h_muon_dirt[i][1]);
    h_muon_ext[i][1]->Sumw2();
    h_muon_bnb[i][1]->Divide(h_muon_bnb[i][1],h_muon_ext[i][1],1,1);
    h_muon_bnb[i][1]->Sumw2();
    h_muon_bnb[i][1]->Draw("e1p");
    h_muon_bnb[i][1]->SetStats(kFALSE);
    h_muon_bnb[i][1]->SetTitle("");
      
    TF1 *a = new TF1("a","1", -150000 , 150000);
    a->SetLineColor(kRed);
    a->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_muon = calculatePearsonChiSq(h_muon_bnb[i][0], h_muon_overlay0[i][2]);
    int nBins1_muon = h_muon_overlay0[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_muon/nBins1_muon));
      
    h_muon_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_muon_bnb[i][1]->GetYaxis()->CenterTitle();
    h_muon_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_muon_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_muon_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_muon_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_muon_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_muon_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",titles_var[i]));
    h_muon_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_muon_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_muon_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_muon_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_muon_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    
    canv_muon[i]->Print(Form("%s_muon%s.png",path.c_str(),var[i]));
    canv_muon[i]->Print(Form("%s_muon%s.pdf",path.c_str(),var[i]));

    //Proton 1
    //cloning histograms for statistics
    h_recoil_ext[i][1] = (TH1D*)h_recoil_ext[i][0]->Clone();
    h_recoil_ext[i][2] = (TH1D*)h_recoil_ext[i][0]->Clone();
    h_recoil_dirt[i][1] = (TH1D*)h_recoil_dirt[i][0]->Clone();
    h_recoil_dirt[i][2] = (TH1D*)h_recoil_dirt[i][0]->Clone(); 
    h_recoil_bnb[i][1] = (TH1D*)h_recoil_bnb[i][0]->Clone();
    h_recoil_overlay0[i][0] = (TH1D*)h_recoil_overlay[i][0]->Clone();
    h_recoil_overlay0[i][1] = (TH1D*)h_recoil_overlay[i][0]->Clone();
    h_recoil_overlay0[i][2] = (TH1D*)h_recoil_overlay[i][0]->Clone();
      
    canv_recoil[i] = new TCanvas(Form("C_recoil%s",var[i]),Form("C_recoil%s",var[i]),2000,1500);
    canv_recoil[i]->cd(1);
    h_recoil[i] = new THStack(Form("h_recoil%s",var[i]),Form("h_recoil%s",var[i]));
    pad_recoil[i] = new TPad(Form("pad_recoil%s",var[i]),Form("pad_recoil%s",var[i]),0,0.35,1.0,1.0);
    pad_recoil[i]->SetBottomMargin(0.0105);
    pad_recoil[i]->SetGridx();
    pad_recoil[i]->SetBorderMode(0);
    pad_recoil[i]->Draw();
    pad_recoil[i]->cd();

    //Stacked Histrogram parameters
    h_recoil[i]->Draw("HIST");
    h_recoil[i]->SetTitle("");
    h_recoil[i]->SetMaximum(recoil_ylim[i]);
    h_recoil[i]->SetMinimum(-2);

    //MC
    if(plot_ccnue){
      z = num_channels;
      f = 11;
    }else{
      z = num_channels -1;
      f = 10;
    }
      
    for(int k=1; k < z ; k++){
      h_recoil_overlay[i][k]->SetLineColor(colors[k]);
      h_recoil_overlay[i][k]->SetFillColor(colors[k]);
      h_recoil_overlay[i][k]->SetLineWidth(1);
      h_recoil[i]->Add(h_recoil_overlay[i][k]);
    }
  
    //Dirt
    h_recoil[i]->Add(h_recoil_dirt[i][0]);
    h_recoil_dirt[i][0]->SetFillColor(kOrange-8);
    h_recoil_dirt[i][0]->SetLineColor(kOrange-8);
    h_recoil_dirt[i][0]->SetLineWidth(1);
   
    //EXT
    h_recoil[i]->Add(h_recoil_ext[i][0]);
    h_recoil_ext[i][0]->SetFillColor(kViolet-7);
    h_recoil_ext[i][0]->SetFillStyle(3005);
    h_recoil_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_recoil_bnb[i][0]->Draw("e1SAME");
    h_recoil_bnb[i][0]->SetLineColor(kBlack);
    h_recoil_bnb[i][0]->SetLineWidth(1);

    h_recoil[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total){
      h_recoil_overlay[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_recoil_overlay0[i][0],h_recoil_ext[i][2], h_recoil_dirt[i][2]);
    h_recoil_overlay0[i][0]->Draw("e2SAME");
    h_recoil_overlay0[i][0]->SetLineColor(kBlack);
    h_recoil_overlay0[i][0]->SetFillColor(kBlack);
    h_recoil_overlay0[i][0]->SetFillStyle(3004);
    h_recoil_overlay0[i][0]->SetMarkerSize(0);
    h_recoil_overlay0[i][0]->SetLineWidth(1);

    legend_recoil[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_recoil[i]->AddEntry(h_recoil_bnb[i][0],"Data (Beam-On)","lepf");
    legend_recoil[i]->AddEntry(h_recoil_overlay0[i][0],"Stat. Unc.","f");
    legend_recoil[i]->AddEntry(h_recoil_ext[i][0],"Data (Beam-Off)","f");
    legend_recoil[i]->AddEntry(h_recoil_dirt[i][0],"Dirt","f");
    for(int k =1; k < z; k++){	
      legend_recoil[i]->AddEntry(h_recoil_overlay[i][f-k],Form("%s",channel_legend[f-k]),"f");  
    }
    legend_recoil[i]->SetLineWidth(0);
    //legend_recoil[i]->SetFillStyle(1);
    legend_recoil[i]->SetFillColor(kWhite);
    legend_recoil[i]->SetTextSize(0.03);
    legend_recoil[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Recoil: %s}",titles_var[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_recoil[i]->cd();
    pad0_recoil[i] = new TPad(Form("pad0_recoil%s",var[i]),Form("pad0_recoil%s",var[i]),0,0.0,1.0,0.35);
    pad0_recoil[i]->SetTopMargin(0);                
    pad0_recoil[i]->SetBottomMargin(0.19);
    pad0_recoil[i]->SetGridx();
    pad0_recoil[i]->Draw();
    pad0_recoil[i]->cd();

    h_recoil_ext[i][1]->Add(h_recoil_overlay0[i][1]);
    h_recoil_ext[i][1]->Add(h_recoil_dirt[i][1]);
    h_recoil_ext[i][1]->Sumw2();
    h_recoil_bnb[i][1]->Divide(h_recoil_bnb[i][1],h_recoil_ext[i][1],1,1);
    h_recoil_bnb[i][1]->Sumw2();
    h_recoil_bnb[i][1]->Draw("e1p");
    h_recoil_bnb[i][1]->SetStats(kFALSE);
    h_recoil_bnb[i][1]->SetTitle("");
      
    TF1 *a2 = new TF1("a2","1", -150000 , 150000);
    a2->SetLineColor(kRed);
    a2->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_recoil = calculatePearsonChiSq(h_recoil_bnb[i][0], h_recoil_overlay0[i][2]);
    int nBins1_recoil = h_recoil_overlay0[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_recoil/nBins1_recoil));
      
    h_recoil_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_recoil_bnb[i][1]->GetYaxis()->CenterTitle();
    h_recoil_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_recoil_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_recoil_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_recoil_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_recoil_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_recoil_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",titles_var[i]));
    h_recoil_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_recoil_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_recoil_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_recoil_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_recoil_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    h_recoil_bnb[i][1]->SetMaximum(3.7);
    h_recoil_bnb[i][1]->SetMinimum(-1.0);
    
    canv_recoil[i]->Print(Form("%s_recoil%s.png",path.c_str(),var[i]));
    canv_recoil[i]->Print(Form("%s_recoil%s.pdf",path.c_str(),var[i]));
  
    //Proton 2
    //cloning histograms for statistics
    h_leading_ext[i][1] = (TH1D*)h_leading_ext[i][0]->Clone();
    h_leading_ext[i][2] = (TH1D*)h_leading_ext[i][0]->Clone();
    h_leading_dirt[i][1] = (TH1D*)h_leading_dirt[i][0]->Clone();
    h_leading_dirt[i][2] = (TH1D*)h_leading_dirt[i][0]->Clone();  
    h_leading_bnb[i][1] = (TH1D*)h_leading_bnb[i][0]->Clone();
    h_leading_overlay0[i][0] = (TH1D*)h_leading_overlay[i][0]->Clone();
    h_leading_overlay0[i][1] = (TH1D*)h_leading_overlay[i][0]->Clone();
    h_leading_overlay0[i][2] = (TH1D*)h_leading_overlay[i][0]->Clone();
      
    canv_leading[i] = new TCanvas(Form("C_leading%s",var[i]),Form("C_leading%s",var[i]),2000,1500);
    canv_leading[i]->cd(1);
    h_leading[i] = new THStack(Form("h_leading%s",var[i]),Form("h_leading%s",var[i]));
    pad_leading[i] = new TPad(Form("pad_leading%s",var[i]),Form("pad_leading%s",var[i]),0,0.35,1.0,1.0);
    pad_leading[i]->SetBottomMargin(0.0);
    pad_leading[i]->SetGridx();
    pad_leading[i]->SetBorderMode(0);
    pad_leading[i]->Draw();
    pad_leading[i]->cd();

    //Stacked Histrogram parameters
    h_leading[i]->Draw("HIST");
    h_leading[i]->SetTitle("");
    h_leading[i]->SetMaximum(leading_ylim[i]);
    h_leading[i]->SetMinimum(-2);

    //MC
    if(plot_ccnue){
      z = num_channels;
      f = 11;
    }else{
      z = num_channels -1;
      f = 10;
    }
      
    for(int k=1; k < z ; k++){
      h_leading_overlay[i][k]->SetLineColor(colors[k]);
      h_leading_overlay[i][k]->SetFillColor(colors[k]);
      h_leading_overlay[i][k]->SetLineWidth(1);
      h_leading[i]->Add(h_leading_overlay[i][k]);
    }

      
    //Dirt
    h_leading[i]->Add(h_leading_dirt[i][0]);
    h_leading_dirt[i][0]->SetFillColor(kOrange-8);
    h_leading_dirt[i][0]->SetLineColor(kOrange-8);
    h_leading_dirt[i][0]->SetLineWidth(1);
   
    //EXT
    h_leading[i]->Add(h_leading_ext[i][0]);
    h_leading_ext[i][0]->SetFillColor(kViolet-7);
    h_leading_ext[i][0]->SetFillStyle(3005);
    h_leading_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_leading_bnb[i][0]->Draw("e1SAME");
    h_leading_bnb[i][0]->SetLineColor(kBlack);
    h_leading_bnb[i][0]->SetLineWidth(1);

    h_leading[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total){
      h_leading_overlay[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_leading_overlay0[i][0],h_leading_ext[i][2],h_leading_dirt[i][2]);
    h_leading_overlay0[i][0]->Draw("e2SAME");
    h_leading_overlay0[i][0]->SetLineColor(kBlack);
    h_leading_overlay0[i][0]->SetFillColor(kBlack);
    h_leading_overlay0[i][0]->SetFillStyle(3004);
    h_leading_overlay0[i][0]->SetMarkerSize(0);
    h_leading_overlay0[i][0]->SetLineWidth(1);

    legend_leading[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_leading[i]->AddEntry(h_leading_bnb[i][0],"Data (Beam-On)","lepf");
    legend_leading[i]->AddEntry(h_leading_overlay0[i][0],"Stat. Unc.","f");
    legend_leading[i]->AddEntry(h_leading_ext[i][0],"Data (Beam-Off)","f");
    legend_leading[i]->AddEntry(h_leading_dirt[i][0],"Dirt","f"); 
    for(int k =1; k < z; k++){	
      legend_leading[i]->AddEntry(h_leading_overlay[i][f-k],Form("%s",channel_legend[f-k]),"f");  
    }
    legend_leading[i]->SetLineWidth(0);
    //legend_leading[i]->SetFillStyle(1);
    legend_leading[i]->SetFillColor(kWhite);
    legend_leading[i]->SetTextSize(0.03);
    legend_leading[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Leading: %s}",titles_var[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_leading[i]->cd();
    pad0_leading[i] = new TPad(Form("pad0_leading%s",var[i]),Form("pad0_leading%s",var[i]),0,0.0,1.0,0.35);
    pad0_leading[i]->SetTopMargin(0);                
    pad0_leading[i]->SetBottomMargin(0.19);
    pad0_leading[i]->SetGridx();
    pad0_leading[i]->Draw();
    pad0_leading[i]->cd();

    h_leading_ext[i][1]->Add(h_leading_overlay0[i][1]);
    h_leading_ext[i][1]->Add(h_leading_dirt[i][1]);
    h_leading_ext[i][1]->Sumw2();
    h_leading_bnb[i][1]->Divide(h_leading_bnb[i][1],h_leading_ext[i][1],1,1);
    h_leading_bnb[i][1]->Sumw2();
    h_leading_bnb[i][1]->Draw("e1p");
    h_leading_bnb[i][1]->SetStats(kFALSE);
    h_leading_bnb[i][1]->SetTitle("");
      
    TF1 *a4 = new TF1("a4","1", -150000 , 150000);
    a4->SetLineColor(kRed);
    a4->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_leading = calculatePearsonChiSq(h_leading_bnb[i][0], h_leading_overlay0[i][2]);
    int nBins1_leading = h_leading_overlay0[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_leading/nBins1_leading));
      
    h_leading_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_leading_bnb[i][1]->GetYaxis()->CenterTitle();
    h_leading_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_leading_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_leading_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_leading_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_leading_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_leading_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",titles_var[i]));
    h_leading_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_leading_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_leading_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_leading_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_leading_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    h_leading_bnb[i][1]->SetMaximum(4.5);
    h_leading_bnb[i][1]->SetMinimum(-2);
    
    canv_leading[i]->Print(Form("%s_leading%s.png",path.c_str(),var[i]));
    canv_leading[i]->Print(Form("%s_leading%s.pdf",path.c_str(),var[i]));
      
   } //end of loop over the particle plots

  
  //PHYSICS PLOTS!!!
  ////////////////////////////////
  for(int i=0; i < num_phys; i++){
    h_phys_ext[i][1] = (TH1D*)h_phys_ext[i][0]->Clone();
    h_phys_ext[i][2] = (TH1D*)h_phys_ext[i][0]->Clone();
    h_phys_dirt[i][1] = (TH1D*)h_phys_dirt[i][0]->Clone();
    h_phys_dirt[i][2] = (TH1D*)h_phys_dirt[i][0]->Clone();  
    h_phys_bnb[i][1] = (TH1D*)h_phys_bnb[i][0]->Clone();
    h_phys_overlay0[i][0] = (TH1D*)h_phys_overlay[i][0]->Clone();
    h_phys_overlay0[i][1] = (TH1D*)h_phys_overlay[i][0]->Clone();
    h_phys_overlay0[i][2] = (TH1D*)h_phys_overlay[i][0]->Clone();
      
    canv_phys[i] = new TCanvas(Form("C_phys%s",var[i]),Form("C_phys%s",var[i]),2000,1500);
    canv_phys[i]->cd(1);
    h_phys[i] = new THStack(Form("h_phys%s",var[i]),Form("h_phys%s",var[i]));
    pad_phys[i] = new TPad(Form("pad_phys%s",var[i]),Form("pad_phys%s",var[i]),0,0.35,1.0,1.0);
    pad_phys[i]->SetBottomMargin(0.0);
    pad_phys[i]->SetGridx();
    pad_phys[i]->SetBorderMode(0);
    pad_phys[i]->Draw();
    pad_phys[i]->cd();

    //Stacked Histrogram parameters
    h_phys[i]->Draw("HIST");
    h_phys[i]->SetTitle("");
    h_phys[i]->SetMaximum(physics_ylim[i]);
    h_phys[i]->SetMinimum(-2);

    //MC
    if(plot_ccnue){
      z = num_channels;
      f = 11;
    }else{
      z = num_channels -1;
      f = 10;
    }
      
    for(int k=1; k < z ; k++){
      h_phys_overlay[i][k]->SetLineColor(colors[k]);
      h_phys_overlay[i][k]->SetFillColor(colors[k]);
      h_phys_overlay[i][k]->SetLineWidth(1);
      h_phys[i]->Add(h_phys_overlay[i][k]);
    }

      
    //Dirt
    h_phys[i]->Add(h_phys_dirt[i][0]);
    h_phys_dirt[i][0]->SetFillColor(kOrange-8);
    h_phys_dirt[i][0]->SetLineColor(kOrange-8);
    h_phys_dirt[i][0]->SetLineWidth(1);
   
    //EXT
    h_phys[i]->Add(h_phys_ext[i][0]);
    h_phys_ext[i][0]->SetFillColor(kViolet-7);
    h_phys_ext[i][0]->SetFillStyle(3005);
    h_phys_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_phys_bnb[i][0]->Draw("e1SAME");
    h_phys_bnb[i][0]->SetLineColor(kBlack);
    h_phys_bnb[i][0]->SetLineWidth(1);

    h_phys[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total){
      h_phys_overlay[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_phys_overlay0[i][0],h_phys_ext[i][2],h_phys_dirt[i][2]);
    h_phys_overlay0[i][0]->Draw("e2SAME");
    h_phys_overlay0[i][0]->SetLineColor(kBlack);
    h_phys_overlay0[i][0]->SetFillColor(kBlack);
    h_phys_overlay0[i][0]->SetFillStyle(3004);
    h_phys_overlay0[i][0]->SetMarkerSize(0);
    h_phys_overlay0[i][0]->SetLineWidth(1);

    legend_phys[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_phys[i]->AddEntry(h_phys_bnb[i][0],"Data (Beam-On)","lepf");
    legend_phys[i]->AddEntry(h_phys_overlay0[i][0],"Stat. Unc.","f");
    legend_phys[i]->AddEntry(h_phys_ext[i][0],"Data (Beam-Off)","f");
    legend_phys[i]->AddEntry(h_phys_dirt[i][0],"Dirt","f"); 
    for(int k =1; k < z; k++){	
      legend_phys[i]->AddEntry(h_phys_overlay[i][f-k],Form("%s",channel_legend[f-k]),"f");  
    }
    legend_phys[i]->SetLineWidth(0);
    //legend_phys[i]->SetFillStyle(1);
    legend_phys[i]->SetFillColor(kWhite);
    legend_phys[i]->SetTextSize(0.03);
    legend_phys[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s}",physics_titles[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_phys[i]->cd();
    pad0_phys[i] = new TPad(Form("pad0_phys%s",var[i]),Form("pad0_phys%s",var[i]),0,0.0,1.0,0.35);
    pad0_phys[i]->SetTopMargin(0);                
    pad0_phys[i]->SetBottomMargin(0.19);
    pad0_phys[i]->SetGridx();
    pad0_phys[i]->Draw();
    pad0_phys[i]->cd();

    h_phys_ext[i][1]->Add(h_phys_overlay0[i][1]);
    h_phys_ext[i][1]->Add(h_phys_dirt[i][1]);
    h_phys_ext[i][1]->Sumw2();
    h_phys_bnb[i][1]->Divide(h_phys_bnb[i][1],h_phys_ext[i][1],1,1);
    h_phys_bnb[i][1]->Sumw2();
    h_phys_bnb[i][1]->Draw("e1p");
    h_phys_bnb[i][1]->SetStats(kFALSE);
    h_phys_bnb[i][1]->SetTitle("");
      
    TF1 *a4 = new TF1("a4","1", -150000 , 150000);
    a4->SetLineColor(kRed);
    a4->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_phys = calculatePearsonChiSq(h_phys_bnb[i][0], h_phys_overlay0[i][2]);
    int nBins1_phys = h_phys_overlay0[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_phys/nBins1_phys));
      
    h_phys_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_phys_bnb[i][1]->GetYaxis()->CenterTitle();
    h_phys_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_phys_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_phys_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_phys_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_phys_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_phys_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",physics_titles[i]));
    h_phys_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_phys_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_phys_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_phys_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_phys_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    
    canv_phys[i]->Print(Form("%s%s.png",path.c_str(),physics[i]));
    canv_phys[i]->Print(Form("%s%s.pdf",path.c_str(),physics[i]));

  }
  
  //PHYSICS PLOTS RAQUEL
  //////////////////////////////
  for(int i=0; i < num_phys; i++){
    h_phys_ext[i][1] = (TH1D*)h_phys_ext[i][0]->Clone();
    h_phys_ext[i][2] = (TH1D*)h_phys_ext[i][0]->Clone();
    h_phys_dirt[i][1] = (TH1D*)h_phys_dirt[i][0]->Clone();
    h_phys_dirt[i][2] = (TH1D*)h_phys_dirt[i][0]->Clone();  
    h_phys_bnb[i][1] = (TH1D*)h_phys_bnb[i][0]->Clone();
    h_phys_overlay0_raquel[i][0] = (TH1D*)h_phys_overlay_raquel[i][0]->Clone();
    h_phys_overlay0_raquel[i][1] = (TH1D*)h_phys_overlay_raquel[i][0]->Clone();
    h_phys_overlay0_raquel[i][2] = (TH1D*)h_phys_overlay_raquel[i][0]->Clone();
      
    canv_phys_raquel[i] = new TCanvas(Form("C_phys_raquel%s",var[i]),Form("C_phys_raquel%s",var[i]),2000,1500);
    canv_phys_raquel[i]->cd(1);
    h_phys_raquel[i] = new THStack(Form("h_phys%s",var[i]),Form("h_phys%s",var[i]));
    pad_phys_raquel[i] = new TPad(Form("pad_phys_raquel%s",var[i]),Form("pad_phys_raquel%s",var[i]),0,0.35,1.0,1.0);
    pad_phys_raquel[i]->SetBottomMargin(0.0);
    pad_phys_raquel[i]->SetGridx();
    pad_phys_raquel[i]->SetBorderMode(0);
    pad_phys_raquel[i]->Draw();
    pad_phys_raquel[i]->cd();

    //Stacked Histrogram parameters
    h_phys_raquel[i]->Draw("HIST");
    h_phys_raquel[i]->SetTitle("");
    h_phys_raquel[i]->SetMaximum(physics_ylim[i]);
    h_phys_raquel[i]->SetMinimum(-2);

    //MC
    if(plot_ccnue_raquel){
      z_raquel = num_channels_raquel;
      f_raquel = 10;
    }else{
      z_raquel = num_channels_raquel -1;
      f_raquel = 9;
    }
    
    for(int k=1; k < z_raquel ; k++){
      h_phys_overlay_raquel[i][k]->SetLineColor(colors_raquel[k]);
      h_phys_overlay_raquel[i][k]->SetFillColor(colors_raquel[k]);
      h_phys_overlay_raquel[i][k]->SetLineWidth(1);
      h_phys_raquel[i]->Add(h_phys_overlay_raquel[i][k]);
    }

    //Dirt
    h_phys_raquel[i]->Add(h_phys_dirt[i][0]);
    h_phys_dirt[i][0]->SetFillColor(kOrange-8);
    h_phys_dirt[i][0]->SetLineColor(kOrange-8);
    h_phys_dirt[i][0]->SetLineWidth(1);
   
    //EXT
    h_phys_raquel[i]->Add(h_phys_ext[i][0]);
    h_phys_ext[i][0]->SetFillColor(kViolet-7);
    h_phys_ext[i][0]->SetFillStyle(3005);
    h_phys_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_phys_bnb[i][0]->Draw("e1SAME");
    h_phys_bnb[i][0]->SetLineColor(kBlack);
    h_phys_bnb[i][0]->SetLineWidth(1);

    h_phys_raquel[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total_raquel){
      h_phys_overlay_raquel[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_phys_overlay0_raquel[i][0],h_phys_ext[i][2],h_phys_dirt[i][2]);
    h_phys_overlay0_raquel[i][0]->Draw("e2SAME");
    h_phys_overlay0_raquel[i][0]->SetLineColor(kBlack);
    h_phys_overlay0_raquel[i][0]->SetFillColor(kBlack);
    h_phys_overlay0_raquel[i][0]->SetFillStyle(3004);
    h_phys_overlay0_raquel[i][0]->SetMarkerSize(0);
    h_phys_overlay0_raquel[i][0]->SetLineWidth(1);

    legend_phys_raquel[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_phys_raquel[i]->AddEntry(h_phys_bnb[i][0],"Data (Beam-On)","lepf");
    legend_phys_raquel[i]->AddEntry(h_phys_overlay0[i][0],"Stat. Unc.","f");
    legend_phys_raquel[i]->AddEntry(h_phys_ext[i][0],"Data (Beam-Off)","f");
    legend_phys_raquel[i]->AddEntry(h_phys_dirt[i][0],"Dirt","f"); 
    for(int k =1; k < z_raquel; k++){	
      legend_phys_raquel[i]->AddEntry(h_phys_overlay_raquel[i][f_raquel-k],Form("%s",channel_legend_raquel[f_raquel-k]),"f");  
    }
    legend_phys_raquel[i]->SetLineWidth(0);
    //legend_phys[i]->SetFillStyle(1);
    legend_phys_raquel[i]->SetFillColor(kWhite);
    legend_phys_raquel[i]->SetTextSize(0.03);
    legend_phys_raquel[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s}",physics_titles[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_phys_raquel[i]->cd();
    pad0_phys_raquel[i] = new TPad(Form("pad0_phys_raquel%s",var[i]),Form("pad0_phys_raquel%s",var[i]),0,0.0,1.0,0.35);
    pad0_phys_raquel[i]->SetTopMargin(0);                
    pad0_phys_raquel[i]->SetBottomMargin(0.19);
    pad0_phys_raquel[i]->SetGridx();
    pad0_phys_raquel[i]->Draw();
    pad0_phys_raquel[i]->cd();

    h_phys_ext[i][1]->Add(h_phys_overlay0_raquel[i][1]);
    h_phys_ext[i][1]->Add(h_phys_dirt[i][1]);
    h_phys_ext[i][1]->Sumw2();
    h_phys_bnb[i][1]->Divide(h_phys_bnb[i][1],h_phys_ext[i][1],1,1);
    h_phys_bnb[i][1]->Sumw2();
    h_phys_bnb[i][1]->Draw("e1p");
    h_phys_bnb[i][1]->SetStats(kFALSE);
    h_phys_bnb[i][1]->SetTitle("");
      
    TF1 *a4 = new TF1("a4","1", -150000 , 150000);
    a4->SetLineColor(kRed);
    a4->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_phys = calculatePearsonChiSq(h_phys_bnb[i][0], h_phys_overlay0_raquel[i][2]);
    int nBins1_phys = h_phys_overlay0_raquel[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_phys/nBins1_phys));
      
    h_phys_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_phys_bnb[i][1]->GetYaxis()->CenterTitle();
    h_phys_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_phys_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_phys_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_phys_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_phys_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_phys_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",physics_titles[i]));
    h_phys_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_phys_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_phys_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_phys_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_phys_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    
    canv_phys_raquel[i]->Print(Form("%s%s_raquel.png",path.c_str(),physics[i]));
    canv_phys_raquel[i]->Print(Form("%s%s_raquel.pdf",path.c_str(),physics[i]));

  }
  

  //STV PLOTS
  ///////////////////////////////////////////
  for(int i=0; i < num_stv; i++){
    h_stv_ext[i][1] = (TH1D*)h_stv_ext[i][0]->Clone();
    h_stv_ext[i][2] = (TH1D*)h_stv_ext[i][0]->Clone();
    h_stv_dirt[i][1] = (TH1D*)h_stv_dirt[i][0]->Clone();
    h_stv_dirt[i][2] = (TH1D*)h_stv_dirt[i][0]->Clone();  
    h_stv_bnb[i][1] = (TH1D*)h_stv_bnb[i][0]->Clone();
    h_stv_overlay0[i][0] = (TH1D*)h_stv_overlay[i][0]->Clone();
    h_stv_overlay0[i][1] = (TH1D*)h_stv_overlay[i][0]->Clone();
    h_stv_overlay0[i][2] = (TH1D*)h_stv_overlay[i][0]->Clone();
      
    canv_stv[i] = new TCanvas(Form("C_stv%s",var[i]),Form("C_stv%s",var[i]),2000,1500);
    canv_stv[i]->cd(1);
    h_stv[i] = new THStack(Form("h_stv%s",var[i]),Form("h_stv%s",var[i]));
    pad_stv[i] = new TPad(Form("pad_stv%s",var[i]),Form("pad_stv%s",var[i]),0,0.35,1.0,1.0);
    pad_stv[i]->SetBottomMargin(0.01);
    pad_stv[i]->SetGridx();
    pad_stv[i]->SetBorderMode(0);
    pad_stv[i]->Draw();
    pad_stv[i]->cd();

    //Stacked Histrogram parameters
    h_stv[i]->Draw("HIST");
    h_stv[i]->SetTitle("");
    h_stv[i]->SetMaximum(stv_ylim[i]);
    h_stv[i]->SetMinimum(-2);

    //MC
    if(plot_ccnue){
      z = num_channels;
      f = 11;
    }else{
      z = num_channels -1;
      f = 10;
    }
      
    for(int k=1; k < z ; k++){
      h_stv_overlay[i][k]->SetLineColor(colors[k]);
      h_stv_overlay[i][k]->SetFillColor(colors[k]);
      h_stv_overlay[i][k]->SetLineWidth(1);
      h_stv[i]->Add(h_stv_overlay[i][k]);
    }

      
    //Dirt
    h_stv[i]->Add(h_stv_dirt[i][0]);
    h_stv_dirt[i][0]->SetFillColor(kOrange-8);
    h_stv_dirt[i][0]->SetLineColor(kOrange-8);
    h_stv_dirt[i][0]->SetLineWidth(1);
   
    //EXT
    h_stv[i]->Add(h_stv_ext[i][0]);
    h_stv_ext[i][0]->SetFillColor(kViolet-7);
    h_stv_ext[i][0]->SetFillStyle(3005);
    h_stv_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_stv_bnb[i][0]->Draw("e1SAME");
    h_stv_bnb[i][0]->SetLineColor(kBlack);
    h_stv_bnb[i][0]->SetLineWidth(1);

    h_stv[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total){
      h_stv_overlay[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_stv_overlay0[i][0],h_stv_ext[i][2],h_stv_dirt[i][2]);
    h_stv_overlay0[i][0]->Draw("e2SAME");
    h_stv_overlay0[i][0]->SetLineColor(kBlack);
    h_stv_overlay0[i][0]->SetFillColor(kBlack);
    h_stv_overlay0[i][0]->SetFillStyle(3004);
    h_stv_overlay0[i][0]->SetMarkerSize(0);
    h_stv_overlay0[i][0]->SetLineWidth(1);

    legend_stv[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_stv[i]->AddEntry(h_stv_bnb[i][0],"Data (Beam-On)","lepf");
    legend_stv[i]->AddEntry(h_stv_overlay0[i][0],"Stat. Unc.","f");
    legend_stv[i]->AddEntry(h_stv_ext[i][0],"Data (Beam-Off)","f");
    legend_stv[i]->AddEntry(h_stv_dirt[i][0],"Dirt","f"); 
    for(int k =1; k < z; k++){	
      legend_stv[i]->AddEntry(h_stv_overlay[i][f-k],Form("%s",channel_legend[f-k]),"f");  
    }
    legend_stv[i]->SetLineWidth(0);
    //legend_stv[i]->SetFillStyle(1);
    legend_stv[i]->SetFillColor(kWhite);
    legend_stv[i]->SetTextSize(0.03);
    legend_stv[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s}",stv_titles[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_stv[i]->cd();
    pad0_stv[i] = new TPad(Form("pad0_stv%s",var[i]),Form("pad0_stv%s",var[i]),0,0.0,1.0,0.35);
    pad0_stv[i]->SetTopMargin(0);                
    pad0_stv[i]->SetBottomMargin(0.19);
    pad0_stv[i]->SetGridx();
    pad0_stv[i]->Draw();
    pad0_stv[i]->cd();

    h_stv_ext[i][1]->Add(h_stv_overlay0[i][1]);
    h_stv_ext[i][1]->Add(h_stv_dirt[i][1]);
    h_stv_ext[i][1]->Sumw2();
    h_stv_bnb[i][1]->Divide(h_stv_bnb[i][1],h_stv_ext[i][1],1,1);
    h_stv_bnb[i][1]->Sumw2();
    h_stv_bnb[i][1]->Draw("e1p");
    h_stv_bnb[i][1]->SetStats(kFALSE);
    h_stv_bnb[i][1]->SetTitle("");
      
    TF1 *a4 = new TF1("a4","1", -150000 , 150000);
    a4->SetLineColor(kRed);
    a4->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_stv = calculatePearsonChiSq(h_stv_bnb[i][0], h_stv_overlay0[i][2]);
    int nBins1_stv = h_stv_overlay0[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_stv/nBins1_stv));
      
    h_stv_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_stv_bnb[i][1]->GetYaxis()->CenterTitle();
    h_stv_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_stv_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_stv_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_stv_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_stv_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_stv_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",stv_titles[i]));
    h_stv_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_stv_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_stv_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_stv_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_stv_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    
    canv_stv[i]->Print(Form("%s%s.png",path.c_str(),stv[i]));
    canv_stv[i]->Print(Form("%s%s.pdf",path.c_str(),stv[i]));

  }
  
  //STV PLOTS: RAQUEL
  ///////////////////////////////////////////
  for(int i=0; i < num_stv; i++){
    h_stv_ext[i][1] = (TH1D*)h_stv_ext[i][0]->Clone();
    h_stv_ext[i][2] = (TH1D*)h_stv_ext[i][0]->Clone();
    h_stv_dirt[i][1] = (TH1D*)h_stv_dirt[i][0]->Clone();
    h_stv_dirt[i][2] = (TH1D*)h_stv_dirt[i][0]->Clone();  
    h_stv_bnb[i][1] = (TH1D*)h_stv_bnb[i][0]->Clone();
    h_stv_overlay0_raquel[i][0] = (TH1D*)h_stv_overlay_raquel[i][0]->Clone();
    h_stv_overlay0_raquel[i][1] = (TH1D*)h_stv_overlay_raquel[i][0]->Clone();
    h_stv_overlay0_raquel[i][2] = (TH1D*)h_stv_overlay_raquel[i][0]->Clone();
      
    canv_stv_raquel[i] = new TCanvas(Form("C_stv_raquel%s",var[i]),Form("C_stv_raquel%s",var[i]),2000,1500);
    canv_stv_raquel[i]->cd(1);
    h_stv_raquel[i] = new THStack(Form("h_stv_raquel%s",var[i]),Form("h_stv_raquel%s",var[i]));
    pad_stv_raquel[i] = new TPad(Form("pad_stv_raquel%s",var[i]),Form("pad_stv_raquel%s",var[i]),0,0.35,1.0,1.0);
    pad_stv_raquel[i]->SetBottomMargin(0.01);
    pad_stv_raquel[i]->SetGridx();
    pad_stv_raquel[i]->SetBorderMode(0);
    pad_stv_raquel[i]->Draw();
    pad_stv_raquel[i]->cd();

    //Stacked Histrogram parameters
    h_stv_raquel[i]->Draw("HIST");
    h_stv_raquel[i]->SetTitle("");
    h_stv_raquel[i]->SetMaximum(stv_ylim[i]);
    h_stv_raquel[i]->SetMinimum(-2);

    //MC                                                                                                                                                                                                                                
    if(plot_ccnue_raquel){
      z_raquel = num_channels_raquel;
      f_raquel = 10;
    }else{
      z_raquel = num_channels_raquel -1;
      f_raquel = 9;
    }
    
    for(int k=1; k < z_raquel ; k++){
      h_stv_overlay_raquel[i][k]->SetLineColor(colors_raquel[k]);
      h_stv_overlay_raquel[i][k]->SetFillColor(colors_raquel[k]);
      h_stv_overlay_raquel[i][k]->SetLineWidth(1);
      h_stv_raquel[i]->Add(h_stv_overlay_raquel[i][k]);
    }

      
    //Dirt
    h_stv_raquel[i]->Add(h_stv_dirt[i][0]);
    h_stv_dirt[i][0]->SetFillColor(kOrange-8);
    h_stv_dirt[i][0]->SetLineColor(kOrange-8);
    h_stv_dirt[i][0]->SetLineWidth(1);
   
    //EXT
    h_stv_raquel[i]->Add(h_stv_ext[i][0]);
    h_stv_ext[i][0]->SetFillColor(kViolet-7);
    h_stv_ext[i][0]->SetFillStyle(3005);
    h_stv_ext[i][0]->SetLineWidth(1);
      
    //BNB
    h_stv_bnb[i][0]->Draw("e1SAME");
    h_stv_bnb[i][0]->SetLineColor(kBlack);
    h_stv_bnb[i][0]->SetLineWidth(1);

    h_stv_raquel[i]->GetYaxis()->SetTitle("No. Events");
      
    //if you want to plot the total for sanity sake:
    if(plot_total){
      h_stv_overlay_raquel[i][0]->Draw("SAME");
    }
      
    //Make sure to do overlay statistical uncertainty
    OverlayStatistics(h_stv_overlay0_raquel[i][0],h_stv_ext[i][2],h_stv_dirt[i][2]);
    h_stv_overlay0_raquel[i][0]->Draw("e2SAME");
    h_stv_overlay0_raquel[i][0]->SetLineColor(kBlack);
    h_stv_overlay0_raquel[i][0]->SetFillColor(kBlack);
    h_stv_overlay0_raquel[i][0]->SetFillStyle(3004);
    h_stv_overlay0_raquel[i][0]->SetMarkerSize(0);
    h_stv_overlay0_raquel[i][0]->SetLineWidth(1);

    legend_stv_raquel[i] = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_stv_raquel[i]->AddEntry(h_stv_bnb[i][0],"Data (Beam-On)","lepf");
    legend_stv_raquel[i]->AddEntry(h_stv_overlay0[i][0],"Stat. Unc.","f");
    legend_stv_raquel[i]->AddEntry(h_stv_ext[i][0],"Data (Beam-Off)","f");
    legend_stv_raquel[i]->AddEntry(h_stv_dirt[i][0],"Dirt","f"); 
    for(int k =1; k < z_raquel; k++){	
      legend_stv_raquel[i]->AddEntry(h_stv_overlay_raquel[i][f_raquel-k],Form("%s",channel_legend_raquel[f_raquel-k]),"f");  
    }
    legend_stv_raquel[i]->SetLineWidth(0);
    //legend_stv_raquel[i]->SetFillStyle(1);
    legend_stv_raquel[i]->SetFillColor(kWhite);
    legend_stv_raquel[i]->SetTextSize(0.03);
    legend_stv_raquel[i]->Draw("same");
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s}",stv_titles[i]));
    t->DrawLatex(0.195,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
    canv_stv_raquel[i]->cd();
    pad0_stv_raquel[i] = new TPad(Form("pad0_stv_raquel%s",var[i]),Form("pad0_stv_raquel%s",var[i]),0,0.0,1.0,0.35);
    pad0_stv_raquel[i]->SetTopMargin(0);                
    pad0_stv_raquel[i]->SetBottomMargin(0.19);
    pad0_stv_raquel[i]->SetGridx();
    pad0_stv_raquel[i]->Draw();
    pad0_stv_raquel[i]->cd();

    h_stv_ext[i][1]->Add(h_stv_overlay0_raquel[i][1]);
    h_stv_ext[i][1]->Add(h_stv_dirt[i][1]);
    h_stv_ext[i][1]->Sumw2();
    h_stv_bnb[i][1]->Divide(h_stv_bnb[i][1],h_stv_ext[i][1],1,1);
    h_stv_bnb[i][1]->Sumw2();
    h_stv_bnb[i][1]->Draw("e1p");
    h_stv_bnb[i][1]->SetStats(kFALSE);
    h_stv_bnb[i][1]->SetTitle("");
      
    TF1 *a4 = new TF1("a4","1", -150000 , 150000);
    a4->SetLineColor(kRed);
    a4->Draw("SAME");

    //Calculate the Chi2
    double chisqv1_stv = calculatePearsonChiSq(h_stv_bnb[i][0], h_stv_overlay0_raquel[i][2]);
    int nBins1_stv = h_stv_overlay0_raquel[i][2]->GetXaxis()->GetNbins();
    //t->DrawLatex(0.83,0.92,Form("#scale[1.5]{#chi^{2}_{stat}/No. Bins: %g / %i = %g"}",chisqv, nBins, chisqv/nBins);
    t->DrawLatex(0.83,0.92,Form("#scale[1.3]{#chi^{2}_{stat}/No. Bins: %.2f}", chisqv1_stv/nBins1_stv));
      
    h_stv_bnb[i][1]->GetYaxis()->SetTitle("Beam-On/(Simulation + Beam-Off)");
    h_stv_bnb[i][1]->GetYaxis()->CenterTitle();
    h_stv_bnb[i][1]->GetYaxis()->SetTitleSize(28);
    h_stv_bnb[i][1]->GetYaxis()->SetTitleFont(43);
    h_stv_bnb[i][1]->GetYaxis()->SetTitleOffset(1.5);
    h_stv_bnb[i][1]->GetYaxis()->SetLabelFont(43);
    h_stv_bnb[i][1]->GetYaxis()->SetLabelSize(30);
    h_stv_bnb[i][1]->GetXaxis()->SetTitle(Form("%s",stv_titles[i]));
    h_stv_bnb[i][1]->GetXaxis()->SetTitleSize(35);
    h_stv_bnb[i][1]->GetXaxis()->SetTitleFont(43);
    h_stv_bnb[i][1]->GetXaxis()->SetTitleOffset(3);
    h_stv_bnb[i][1]->GetXaxis()->SetLabelFont(43);
    h_stv_bnb[i][1]->GetXaxis()->SetLabelSize(35);
    
    canv_stv_raquel[i]->Print(Form("%s%s_raquel.png",path.c_str(),stv[i]));
    canv_stv_raquel[i]->Print(Form("%s%s_raquel.pdf",path.c_str(),stv[i]));

  }
  


  
  
} //end of program
