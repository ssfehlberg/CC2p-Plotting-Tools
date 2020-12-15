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
  TFile *f1=new TFile(Form("histograms_%s_wgt.root",sample));//overlay histograms. Note: be sure to use weighted for both dirt and overlay
  TFile *f_dirt=new TFile(Form("histograms_%s_dirt_wgt.root",sample));//dirt histograms 
  TFile *f2=new TFile(Form("histograms_%s_bnb.root",sample));//bnb histograms  
  TFile *f3=new TFile(Form("histograms_%s_ext.root",sample));//extbnb histograms 
  TFile *f4=new TFile("histograms_efficiency.root"); //efficiency histograms
  
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

  ///////////////////////////////////
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

  //////////////////////////
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

  //////////////////////////////////
  //Plot the 2D Histograms:
  //////////////////////////////////
  for(int i =0; i < num_group2d; i++){
    canv_2d[i] = new TCanvas(Form("C_2D_%s",group2d[i]),Form("C_2D_%s",group2d[i]),2000,1500);
    h_overlay2D[i]->Draw("colz"); 
    canv_2d[i]->Print(Form("%s_%s.png",path.c_str(),group2d[i]));
    canv_2d[i]->Print(Form("%s_%s.pdf",path.c_str(),group2d[i]));
  }

  /////////////////////////////
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

  ///////////////////////////////////
  //Plots of the Reconstructed Vertex
  ///////////////////////////////////
  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_variables; j++){
      for(int k=0; k < num_channels; k++){
	h_overlay_vec.push_back(h_overlay[i][j][k]);
      }
      Plot_Histograms(colors, h_overlay_vec, h_overlay0[i][j][0],h_overlay0[i][j][1],h_overlay0[i][j][2], h_ext[i][j][0], h_ext[i][j][1], h_ext[i][j][2], h_dirt[i][j][0],h_dirt[i][j][1],h_dirt[i][j][2],h_bnb[i][j][0],h_bnb[i][j][1],canv[i][j], h[i][j], pad[i][j], pad0[i][j], legend[i][j],ylim[i][j],ymin[i][j], num_channels, titles[j], path, "", plots[j], cut[i], false, false);
      h_overlay_vec.clear();
    }
  }

  ////////////////////////////////////////
  //Chi2 Plots
  ////////////////////////////////////////
  for(int i = 0; i< num_planes; i++){
    for(int j = 0; j< num_hypothesis; j++){
      for(int k=0; k < num_particles; k++){
	h_overlay_chi2_vec.push_back(h_overlay_chi2[i][j][k]);
      } 
      Plot_Histograms(colors_chi2, h_overlay_chi2_vec, h_overlay0_chi2[i][j][0],h_overlay0_chi2[i][j][1],h_overlay0_chi2[i][j][2], h_ext_chi2[i][j][0], h_ext_chi2[i][j][1], h_ext_chi2[i][j][2], h_dirt_chi2[i][j][0],h_dirt_chi2[i][j][1],h_dirt_chi2[i][j][2],h_bnb_chi2[i][j][0],h_bnb_chi2[i][j][1],canv_chi2[i][j], h_chi2[i][j], pad_chi2[i][j], pad0_chi2[i][j], legend_chi2[i][j],ylim_chi2[i][j], ymin_chi2[i][j], num_particles, titles_chi2[j], path,titles_planes[i], plane[i], hypothesis[j], false, false);
      h_overlay_chi2_vec.clear();
    }
  }

  ///////////////////////////////
  //3D Chi2 Plots
  ///////////////////////////////   
  for(int i = 0; i < num_cuts_3D; i++){
    for(int j = 0; j< num_hypothesis_3D; j++){
      for(int k=0; k < num_particles; k++){
	h_overlay_chi2_3D_vec.push_back(h_overlay_chi2_3D[i][j][k]);
      } 
      Plot_Histograms(colors_chi2, h_overlay_chi2_3D_vec, h_overlay0_chi2_3D[i][j][0],h_overlay0_chi2_3D[i][j][1],h_overlay0_chi2_3D[i][j][2],
		      h_ext_chi2_3D[i][j][0], h_ext_chi2_3D[i][j][1], h_ext_chi2_3D[i][j][2],
		      h_dirt_chi2_3D[i][j][0],h_dirt_chi2_3D[i][j][1],h_dirt_chi2_3D[i][j][2],
		      h_bnb_chi2_3D[i][j][0],h_bnb_chi2_3D[i][j][1],canv_chi2_3D[i][j], h_chi2_3D[i][j],
		      pad_chi2_3D[i][j], pad0_chi2_3D[i][j], legend_chi2_3D[i][j],
		      ylim_chi2_3D[i][j], ymin_chi2_3D[i][j], num_particles, titles_chi2_3D[j], path, cuts_3D_titles[i], Form("3D_%s",hypothesis_3D[j]), cuts_3D[i],false, false);
      h_overlay_chi2_3D_vec.clear();
    } 
  }
  
  //////////////////////////////
  //Particle specific plots
  /////////////////////////////
  for(int i = 0; i < num_var; i++){
    for(int k=0; k < num_channels; k++){
      h_muon_overlay_vec.push_back(h_muon_overlay[i][k]);
      h_recoil_overlay_vec.push_back(h_recoil_overlay[i][k]);
      h_leading_overlay_vec.push_back(h_leading_overlay[i][k]);
    }
    for(int k=0; k < num_channels_raquel; k++){
      h_muon_overlay_raquel_vec.push_back(h_muon_overlay_raquel[i][k]);
      h_recoil_overlay_raquel_vec.push_back(h_recoil_overlay_raquel[i][k]);
      h_leading_overlay_raquel_vec.push_back(h_leading_overlay_raquel[i][k]);
      
    }
    
    //muon:mine
    Plot_Histograms(colors, h_muon_overlay_vec, h_muon_overlay0[i][0],h_muon_overlay0[i][1],h_muon_overlay0[i][2],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon[i], h_muon[i],
		    pad_muon[i], pad0_muon[i], legend_muon[i],
		    muon_ylim[i], muon_ymin[i], num_channels, Form("Muon: %s",titles_var[i]), path,"", Form("_muon%s",var[i]), "",false, false);
    h_muon_overlay_vec.clear();
    
    //muon:raquel
    /* Plot_Histograms(colors_raquel, h_muon_overlay_raquel_vec, h_muon_overlay0_raquel[i][0],h_muon_overlay0_raquel[i][1],h_muon_overlay0_raquel[i][2],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon_raquel[i], h_muon_raquel[i],
		    pad_muon_raquel[i], pad0_muon_raquel[i], legend_muon_raquel[i],
		    muon_ylim[i], muon_ymin[i], num_channels_raquel, Form("Muon: %s",titles_var[i]), path,"", Form("_muon_raquel%s",var[i]), "",false, false);
    h_muon_overlay_raquel_vec.clear();
    */

    //recoil proton:mine
    Plot_Histograms(colors, h_recoil_overlay_vec, h_recoil_overlay0[i][0],h_recoil_overlay0[i][1],h_recoil_overlay0[i][2],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil[i], h_recoil[i],
		    pad_recoil[i], pad0_recoil[i], legend_recoil[i],
		    recoil_ylim[i], recoil_ymin[i], num_channels, Form("Recoil Proton: %s",titles_var[i]), path,"", Form("_recoil%s",var[i]), "",false, false);
    h_recoil_overlay_vec.clear();

    //recoil proton:raquel
    /*Plot_Histograms(colors_raquel, h_recoil_overlay_raquel_vec, h_recoil_overlay0_raquel[i][0],h_recoil_overlay0_raquel[i][1],h_recoil_overlay0_raquel[i][2],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil_raquel[i], h_recoil_raquel[i],
		    pad_recoil_raquel[i], pad0_recoil_raquel[i], legend_recoil_raquel[i],
		    recoil_ylim[i], recoil_ymin[i], num_channels_raquel, Form("Recoil Proton: %s",titles_var[i]), path,"", Form("_recoil_raquel%s",var[i]), "",false, false);
    h_recoil_overlay_raquel_vec.clear();
    */
    //leading proton:mine
    Plot_Histograms(colors, h_leading_overlay_vec, h_leading_overlay0[i][0],h_leading_overlay0[i][1],h_leading_overlay0[i][2],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading[i], h_leading[i],
		    pad_leading[i], pad0_leading[i], legend_leading[i],
		    leading_ylim[i], leading_ymin[i], num_channels, Form("Leading Proton: %s",titles_var[i]), path,"", Form("_leading%s",var[i]), "",false, false);
    h_leading_overlay_vec.clear();

    //leading proton:raquel
    /*Plot_Histograms(colors_raquel, h_leading_overlay_raquel_vec, h_leading_overlay0_raquel[i][0],h_leading_overlay0_raquel[i][1],h_leading_overlay0_raquel[i][2],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading_raquel[i], h_leading_raquel[i],
		    pad_leading_raquel[i], pad0_leading_raquel[i], legend_leading_raquel[i],
		    leading_ylim[i], leading_ymin[i], num_channels_raquel, Form("Leading Proton: %s",titles_var[i]), path,"", Form("_leading_raquel%s",var[i]), "",false, false);
    h_leading_overlay_raquel_vec.clear();
    */
  }
    
  //PHYSICS PLOTS!!!
  ////////////////////////////////
  for(int i=0; i < num_phys; i++){
    for(int j=0; j < num_channels; j++){
      h_phys_overlay_vec.push_back(h_phys_overlay[i][j]);
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_phys_overlay_raquel_vec.push_back(h_phys_overlay_raquel[i][j]);
    }
    
    //mine 
    Plot_Histograms(colors, h_phys_overlay_vec, h_phys_overlay0[i][0],h_phys_overlay0[i][1],h_phys_overlay0[i][2],
                    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys[i], h_phys[i],
                    pad_phys[i], pad0_phys[i], legend_phys[i],
                    phys_ylim[i], phys_ymin[i], num_channels, physics_titles[i], path,"", physics[i], "",false, false);
    h_phys_overlay_vec.clear();

    //raquel
    Plot_Histograms(colors_raquel, h_phys_overlay_raquel_vec, h_phys_overlay0_raquel[i][0],h_phys_overlay0_raquel[i][1],h_phys_overlay0_raquel[i][2],
                    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys_raquel[i], h_phys_raquel[i],
                    pad_phys_raquel[i], pad0_phys_raquel[i], legend_phys_raquel[i],
                    phys_ylim[i], phys_ymin[i], num_channels_raquel, physics_titles[i], path,"", Form("%s_raquel",physics[i]), "",false, false);
    h_phys_overlay_raquel_vec.clear();
  }


  //STV PLOTS
  ///////////////////////////////////////////
  for(int i=0; i < num_stv; i++){
    for(int j=0; j < num_channels; j++){
      h_stv_overlay_vec.push_back(h_stv_overlay[i][j]);
      h_stv_overlay[i][j]->Draw("hist");
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_stv_overlay_raquel_vec.push_back(h_stv_overlay_raquel[i][j]);
    }

    //mine 
    Plot_Histograms(colors, h_stv_overlay_vec, h_stv_overlay0[i][0],h_stv_overlay0[i][1],h_stv_overlay0[i][2],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv[i], h_stv[i],
                    pad_stv[i], pad0_stv[i], legend_stv[i],
                    stv_ylim[i], stv_ymin[i], num_channels, stv_titles[i], path,"", stv[i], "",false, false);
    h_stv_overlay_vec.clear();
    
    //raquel
    Plot_Histograms(colors_raquel, h_stv_overlay_raquel_vec, h_stv_overlay0_raquel[i][0],h_stv_overlay0_raquel[i][1],h_stv_overlay0_raquel[i][2],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv_raquel[i], h_stv_raquel[i],
                    pad_stv_raquel[i], pad0_stv_raquel[i], legend_stv_raquel[i],
                    stv_ylim[i], stv_ymin[i], num_channels_raquel, stv_titles[i], path,"", Form("%s_raquel",stv[i]), "",false, false);
    h_stv_overlay_raquel_vec.clear();
    

  }
    
  /*
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
  

    */
  
  
} //end of program
