#define analysis_cxx
#include "analysis.h"
//#include "paul_tol_colors.hpp"
//#include <iostream>
//#include <ctime>
//#include <string>


void analysis::main(){

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //FIRST: CHECK WHICH SAMPLE THIS IS: PELEE,FILTERED,OR UNFILTERED & WHAT RUN: JAN, RUN 1,RUN 2,RUN or 3
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  const char* sample = which_sample();
  std::pair<const char*, const char*> r = which_run();
  const char* run = r.first;
  const char* run_title = r.second;

  ///////////////////////////////////////////////////////
  //SECOND: DEFINE HISTOGRAMS FILES AND GENERAL VARIABLES
  //////////////////////////////////////////////////////
  TFile *f1=new TFile(Form("root_files/%s/%s/histograms_%s_overlay_wgt.root",sample,run,sample));//overlay histograms.
  TFile *f_dirt=new TFile(Form("root_files/%s/%s/histograms_%s_dirt_wgt.root",sample,run,sample));//dirt histograms 
  TFile *f2=new TFile(Form("root_files/%s/%s/histograms_%s_bnb.root",sample,run,sample));//bnb histograms  
  TFile *f3=new TFile(Form("root_files/%s/%s/histograms_%s_ext.root",sample,run,sample));//extbnb histograms 
  TFile *f4=new TFile("root_files/unfiltered/histograms_efficiency.root"); //This is only for unfiltered. efficinecy is now baked into the overlay pelee file

  //Define Parameters and color scheme
  ////////////////////////////////////
  Define_Parameters(run);
  
  //Color Scheme
  tolcols::init();
  Color_t colors[] = {0,9031, 9030, 9029, 9028, 9026, 9025, 9024, 9032, kGray+2, 9027};  //Black, light pink, light orange, light yellow, olive, mint, light cyan, light blue, light grey, darker grey, olive
  Color_t colors_raquel[] = {0,9012,9011,9010,9009,9008,9007,kGray+2,9032,9027}; //black, magenta, red, orange, mint, cyan, blue, dark gray, light gray, ccnue
  Color_t colors_chi2[10] = {kBlack,kRed,kRed+2,kBlue,kGreen,kGreen+1,kYellow,kYellow-3,kOrange+8,204};

  /////////////////////////////////////////////////////////////
  //THIRD: MAKE DIRECTORY WITH TODAY'S DATE TO STORE ALL IMAGES
  //////////////////////////////////////////////////////////////
  const char* pathname = Form("images/%d%d%d_%s_%s/",Month,Day,Year,sample,run);
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
      h_truth[i][j]->SetMaximum(ylim_truth[run_num][i][j]);
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
  /*
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
  */
  
  //////////////////////////////////
  //Plot the 2D Histograms:
  //////////////////////////////////
  /*  for(int i =0; i < num_group2d; i++){
    canv_2d[i] = new TCanvas(Form("C_2D_%s",group2d[i]),Form("C_2D_%s",group2d[i]),2000,1500);
    h_overlay2D[i]->Draw("colz"); 
    canv_2d[i]->Print(Form("%s_%s.png",path.c_str(),group2d[i]));
    canv_2d[i]->Print(Form("%s_%s.pdf",path.c_str(),group2d[i]));
  }
  */
  
  /////////////////////////////
  //Plot the Efficiency Stuff:
  //////////////////////////
  /*for(int i=0; i < num_eff; i++){

    h_num0[i] = (TH1D*)h_num[i]->Clone();
    h_denom0[i] = (TH1D*)h_denom[i]->Clone();
    canv_both[i] = new TCanvas(Form("C_both_%s",eff[i]),Form("C_both_%s",eff[i]),2000,1500);
    h_denom0[i]->Draw("hist");
    h_denom0[i]->SetFillColor(kBlue);
    h_denom0[i]->SetTitle(Form("Efficiency(%s); %s ; Number of Entries",titles_eff[i],titles_eff[i]));
    h_denom0[i]->SetMaximum(ylim_eff[run_num][i]);
    h_num0[i]->Draw("histSAME");
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
    t->DrawLatex(0.23,0.92,Form("%s",pot_num));
    t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
    a[i] = new TLine(xlim_eff[run_num][i],0,xlim_eff[run_num][i],1);
    a[i]->Draw("same");
    a[i]->SetLineColor(kBlack);
    a[i]->SetLineWidth(4);
    canv_eff[i]->Print(Form("%s_%s_eff.png",path.c_str(),eff[i]));
    canv_eff[i]->Print(Form("%s_%s_eff.pdf",path.c_str(),eff[i]));
    }
  */
  //efficiency and purity as function of cuts
  TCanvas* canv_eff_pur = new TCanvas("canv_eff_pur","canv_eff_pur",2000,1500);
  eff_graph->Draw("alp");
  eff_graph->SetLineColor(kBlue);
  eff_graph->SetMarkerColor(kBlue);
  eff_graph->SetLineWidth(3);
  eff_graph->SetMarkerStyle(20);
  //eff_graph->GetXaxis()->SetLabel("Cut Number");
  eff_graph->SetTitle("Efficiency and Purity as a Function of Cuts");
  eff_graph->GetYaxis()->SetRangeUser(0.,1.);
  eff_graph->GetYaxis()->SetTitle("Efficiency or Purity");
  const char* eff_cut_label[6] = {"No Cuts","FV","3 PFP","Trk Scr","Vtx","PID Cut"};
  TAxis *xax = eff_graph->GetXaxis();
  for(int i = 1; i < xax->GetXmax()-1; i++){
    int bin_index_up = xax->FindBin(i+1);
    int bin_index_low = xax->FindBin(i);
    xax->SetBinLabel(bin_index_low,Form("%s", eff_cut_label[i-1]));
    xax->SetBinLabel(bin_index_up,Form("%s", eff_cut_label[i]));
  }
  xax->Draw("hist");
  
  pur_graph->Draw("SAME lp");
  pur_graph->SetLineColor(kRed);
  pur_graph->SetMarkerColor(kRed);
  pur_graph->SetLineWidth(3);
  pur_graph->SetMarkerStyle(20);

  TLegend* lg = new TLegend(0.65,0.65,0.85,0.85);
  lg->AddEntry(eff_graph,"Efficiency","lp");
  lg->AddEntry(pur_graph,"Purity","lp");
  lg->Draw("same");
 
  canv_eff_pur->Print(Form("%s_efficiency_and_purity.png",path.c_str()));
  canv_eff_pur->Print(Form("%s_efficiency_and_purity.pdf",path.c_str()));

  //Effieincy of other variables
  for(int i =0; i < num_particles_eff; i++){
    for (int j=0; j <num_particles_eff_plots; j++){
      canv_particle_eff[i][j] = new TCanvas(Form("C_eff%s%s",particles_eff[j],particles_eff_var[i]),Form("C_eff%s%s",particles_eff[j],particles_eff_var[i]),2000,1500);

      h_particle_num[i][j]->Divide(h_particle_num[i][j],h_particle_denom[i][j],1.0,1.0, "B");
      h_particle_num[i][j]->Draw("1e1p");
      h_particle_num[i][j]->SetTitle(Form(" ; %s ; Efficiency",particles_eff_var_titles[j]));
      h_particle_num[i][j]->SetLineColor(kViolet);
      h_particle_num[i][j]->SetMaximum(1);
      h_particle_num[i][j]->SetMinimum(0);
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s of %s}",particles_eff_var_titles[j],particles_eff_titles[i]));
      t->DrawLatex(0.23,0.92,Form("%s",pot_num));
      t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
      //a[i] = new TLine(xlim_eff[run_num][i],0,xlim_eff[run_num][i],1);
      //a[i]->Draw("same");
      //a[i]->SetLineColor(kBlack);
      //a[i]->SetLineWidth(4);

      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff.png",path.c_str(),particles_eff[i],particles_eff_var[j]));
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff.pdf",path.c_str(),particles_eff[i],particles_eff_var[j]));
 
    }
  }
  

  for(int i =0; i < num_other_eff; i++){
    canv_other_eff[i] = new TCanvas(Form("C_other_eff%s",other_eff[i]),Form("C_other_eff%s",other_eff[i]),2000,1500);
    h_other_eff_num[i]->Divide(h_other_eff_num[i],h_other_eff_denom[i],1.0,1.0, "B");
    h_other_eff_num[i]->Draw("1e1p");
    h_other_eff_num[i]->SetTitle(Form(" ; %s ; Efficiency",other_eff_titles[i]));
    h_other_eff_num[i]->SetLineColor(kViolet);
    h_other_eff_num[i]->SetMaximum(1);
    h_other_eff_num[i]->SetMinimum(0);
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s}",other_eff_titles[i]));
    t->DrawLatex(0.23,0.92,Form("%s",pot_num));
    t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
    //a[i] = new TLine(xlim_eff[run_num][i],0,xlim_eff[run_num][i],1);
    //a[i]->Draw("same");
    //a[i]->SetLineColor(kBlack);
    //a[i]->SetLineWidth(4);
    canv_other_eff[i]->Print(Form("%s%s_eff.png",path.c_str(),other_eff[i]));
    canv_other_eff[i]->Print(Form("%s%s_eff.pdf",path.c_str(),other_eff[i]));

  }

  

  ///////////////////////////////////
  //Plots of the Reconstructed Vertex
  ///////////////////////////////////
  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_variables; j++){
      for(int k=0; k < num_channels; k++){
	h_overlay_vec.push_back(h_overlay[i][j][k]);
      }
      for(int k=0; k < num_channels_raquel; k++){
	h_overlay_raquel_vec.push_back(h_overlay_raquel[i][j][k]);
      }
      
      
      Plot_Histograms(pot_num,sample_name,colors, h_overlay_vec, h_overlay0[i][j][0],h_overlay0[i][j][1],h_overlay0[i][j][2], h_ext[i][j][0], h_ext[i][j][1], h_ext[i][j][2], h_dirt[i][j][0],h_dirt[i][j][1],h_dirt[i][j][2],h_bnb[i][j][0],h_bnb[i][j][1],canv[i][j], h[i][j], pad[i][j], pad0[i][j], legend[i][j], channel_legend, ylim[run_num][i][j],ymin[run_num][i][j], num_channels, titles[j], path, a1 ,"", plots[j], cut[i], false, false, false, 0.0, 0.19, 0.5, 1.5);
      h_overlay_vec.clear();

      Plot_Histograms(pot_num,sample_name,colors_raquel, h_overlay_raquel_vec, h_overlay0_raquel[i][j][0],h_overlay0_raquel[i][j][1],h_overlay0_raquel[i][j][2], h_ext[i][j][0], h_ext[i][j][1], h_ext[i][j][2], h_dirt[i][j][0],h_dirt[i][j][1],h_dirt[i][j][2],h_bnb[i][j][0],h_bnb[i][j][1],canv_raquel[i][j], h_raquel[i][j], pad_raquel[i][j], pad0_raquel[i][j], legend_raquel[i][j], channel_legend_raquel, ylim[run_num][i][j],ymin[run_num][i][j], num_channels_raquel, titles[j], path,   a1,"", plots[j], Form("%s_raquel",cut[i]), false, false, false, 0.0, 0.19, 0.5, 1.5);
      h_overlay_raquel_vec.clear();
     
    }
  }
    
  /*
  ////////////////////////////////////////
  //Chi2 Plots
  ////////////////////////////////////////
  for(int i = 0; i< num_planes; i++){
    for(int j = 0; j< num_hypothesis; j++){
      for(int k=0; k < num_particles; k++){
	h_overlay_chi2_vec.push_back(h_overlay_chi2[i][j][k]);
      } 
      Plot_Histograms(pot_num,sample_name,colors_chi2, h_overlay_chi2_vec, h_overlay0_chi2[i][j][0],h_overlay0_chi2[i][j][1],h_overlay0_chi2[i][j][2], h_ext_chi2[i][j][0], h_ext_chi2[i][j][1], h_ext_chi2[i][j][2], h_dirt_chi2[i][j][0],h_dirt_chi2[i][j][1],h_dirt_chi2[i][j][2],h_bnb_chi2[i][j][0],h_bnb_chi2[i][j][1],canv_chi2[i][j], h_chi2[i][j], pad_chi2[i][j], pad0_chi2[i][j], legend_chi2[i][j],channel_legend_chi2,ylim_chi2[i][j], ymin_chi2[i][j], num_particles, titles_chi2[j], path,titles_planes[i], plane[i], hypothesis[j], false, false);
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
      Plot_Histograms(pot_num,sample_name,colors_chi2, h_overlay_chi2_3D_vec, h_overlay0_chi2_3D[i][j][0],h_overlay0_chi2_3D[i][j][1],h_overlay0_chi2_3D[i][j][2],
		      h_ext_chi2_3D[i][j][0], h_ext_chi2_3D[i][j][1], h_ext_chi2_3D[i][j][2],
		      h_dirt_chi2_3D[i][j][0],h_dirt_chi2_3D[i][j][1],h_dirt_chi2_3D[i][j][2],
		      h_bnb_chi2_3D[i][j][0],h_bnb_chi2_3D[i][j][1],canv_chi2_3D[i][j], h_chi2_3D[i][j],
		      pad_chi2_3D[i][j], pad0_chi2_3D[i][j], legend_chi2_3D[i][j],channel_legend_chi2,
		      ylim_chi2_3D[i][j], ymin_chi2_3D[i][j], num_particles, titles_chi2_3D[j], path, cuts_3D_titles[i], Form("3D_%s",hypothesis_3D[j]), cuts_3D[i],false, false, 0.01, 0.15);
      h_overlay_chi2_3D_vec.clear();
    } 
  }
 
  */

  //////////////////////////
  //Track plots, such as PID
  //////////////////////////
  for(int i=0; i < num_track; i++){
    for(int k=0; k < track_cut; k++){
      for(int j=0; j < num_particles; j ++){
	h_track_overlay_vec.push_back(h_track_overlay[i][k][j]);
      }
    
      Plot_Histograms(pot_num,sample_name,colors_chi2, h_track_overlay_vec, h_track_overlay0[i][k][0],h_track_overlay0[i][k][1],h_track_overlay0[i][k][2],
		      h_track_ext[i][k][0], h_track_ext[i][k][1], h_track_ext[i][k][2],
		      h_track_dirt[i][k][0],h_track_dirt[i][k][1],h_track_dirt[i][k][2],
		      h_track_bnb[i][k][0],h_track_bnb[i][k][1],canv_track[i][k], h_track[i][k],
		      pad_track[i][k], pad0_track[i][k], legend_track[i][k], channel_legend_chi2,
		      ymax_track[i], ymin_track[i], num_particles, Form("%s",titles_track[i]), path, a_track[i],Form(""), Form("%s%s",variable[i],which_track_cut[k]), Form(""),false, true, true, 0.01, 0.15, xlim_track[i]);
      h_track_overlay_vec.clear();
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
    Plot_Histograms(pot_num,sample_name,colors, h_muon_overlay_vec, h_muon_overlay0[i][0],h_muon_overlay0[i][1],h_muon_overlay0[i][2],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon[i], h_muon[i],
		    pad_muon[i], pad0_muon[i], legend_muon[i],channel_legend,
		    muon_ylim[run_num][i], muon_ymin[run_num][i], num_channels, Form("Muon: %s",titles_var[i]), path,  a1,"", Form("_muon%s",var[i]), "",false, false,flip_muon[i],0.0,0.19,0.0,2.7);
    h_muon_overlay_vec.clear();
       
    //muon:raquel
     Plot_Histograms(pot_num,sample_name,colors_raquel, h_muon_overlay_raquel_vec, h_muon_overlay0_raquel[i][0],h_muon_overlay0_raquel[i][1],h_muon_overlay0_raquel[i][2],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon_raquel[i], h_muon_raquel[i],
		     pad_muon_raquel[i], pad0_muon_raquel[i], legend_muon_raquel[i], channel_legend_raquel,
		     muon_ylim[run_num][i], muon_ymin[run_num][i], num_channels_raquel, Form("Muon: %s",titles_var[i]), path,  a1,"", Form("_muon_raquel%s",var[i]), "",false, false,flip_muon[i],0.0,0.19,0.0,2.7);
    h_muon_overlay_raquel_vec.clear();
    
      
    //recoil proton:mine
    Plot_Histograms(pot_num,sample_name,colors, h_recoil_overlay_vec, h_recoil_overlay0[i][0],h_recoil_overlay0[i][1],h_recoil_overlay0[i][2],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil[i], h_recoil[i],
		    pad_recoil[i], pad0_recoil[i], legend_recoil[i],channel_legend,
		    recoil_ylim[run_num][i], recoil_ymin[run_num][i], num_channels, Form("Recoil Proton: %s",titles_var[i]), path,  a1,"", Form("_recoil%s",var[i]), "",false, false, flip_recoil[i],0.0105, 0.19, 0.0, 2.7);
    h_recoil_overlay_vec.clear();

    //recoil proton:raquel
    Plot_Histograms(pot_num,sample_name,colors_raquel, h_recoil_overlay_raquel_vec, h_recoil_overlay0_raquel[i][0],h_recoil_overlay0_raquel[i][1],h_recoil_overlay0_raquel[i][2],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil_raquel[i], h_recoil_raquel[i],
		    pad_recoil_raquel[i], pad0_recoil_raquel[i], legend_recoil_raquel[i], channel_legend_raquel,
		    recoil_ylim[run_num][i], recoil_ymin[run_num][i], num_channels_raquel, Form("Recoil Proton: %s",titles_var[i]), path,  a1,"", Form("_recoil_raquel%s",var[i]), "",false, false,flip_recoil[i],0.0,0.19,0.0,2.7);
    h_recoil_overlay_raquel_vec.clear();
    
    //leading proton:mine
    Plot_Histograms(pot_num,sample_name,colors, h_leading_overlay_vec, h_leading_overlay0[i][0],h_leading_overlay0[i][1],h_leading_overlay0[i][2],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading[i], h_leading[i],
		    pad_leading[i], pad0_leading[i], legend_leading[i],channel_legend,
		    leading_ylim[run_num][i], leading_ymin[run_num][i], num_channels, Form("Leading Proton: %s",titles_var[i]), path,  a1,"", Form("_leading%s",var[i]), "",false, false, flip_lead[i],0.0,0.19,0.0,2.7);
    h_leading_overlay_vec.clear();

    //leading proton:raquel
    Plot_Histograms(pot_num,sample_name,colors_raquel, h_leading_overlay_raquel_vec, h_leading_overlay0_raquel[i][0],h_leading_overlay0_raquel[i][1],h_leading_overlay0_raquel[i][2],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading_raquel[i], h_leading_raquel[i],
		    pad_leading_raquel[i], pad0_leading_raquel[i], legend_leading_raquel[i],channel_legend_raquel,
		    leading_ylim[run_num][i], leading_ymin[run_num][i], num_channels_raquel, Form("Leading Proton: %s",titles_var[i]), path,  a1,"", Form("_leading_raquel%s",var[i]), "",false, false, flip_lead[i],0.0,0.19,0.0,2.7);
    h_leading_overlay_raquel_vec.clear();
    
  }

  
  //PHYSICS PLOTS!!!
  ////////////////////////////////
  for(int i=0; i < num_phys; i++){
    for(int j=0; j < num_channels; j++){
      h_phys_overlay_vec.push_back(h_phys_overlay[i][j]);
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_phys_overlay_raquel_vec.push_back(h_phys_overlay_raquel[i][j]);
      h_phys_overlay_raquel[i][j]->Draw("hist");
    }

    //mine 
    Plot_Histograms(pot_num,sample_name,colors, h_phys_overlay_vec, h_phys_overlay0[i][0],h_phys_overlay0[i][1],h_phys_overlay0[i][2],
                    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys[i], h_phys[i],
                    pad_phys[i], pad0_phys[i], legend_phys[i],channel_legend,
                    phys_ylim[run_num][i], phys_ymin[run_num][i], num_channels, physics_titles[i], path,  a1,"", physics[i], "",false, false, false,0.0, 0.19,0.0,2.7);
    h_phys_overlay_vec.clear();
     
    //raquel
    Plot_Histograms(pot_num,sample_name,colors_raquel, h_phys_overlay_raquel_vec, h_phys_overlay0_raquel[i][0],h_phys_overlay0_raquel[i][1],h_phys_overlay0_raquel[i][2],
		    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys_raquel[i], h_phys_raquel[i],
                    pad_phys_raquel[i], pad0_phys_raquel[i], legend_phys_raquel[i],channel_legend_raquel,
                    phys_ylim[run_num][i], phys_ymin[run_num][i], num_channels_raquel, physics_titles[i], path,  a1,"", Form("%s_raquel",physics[i]), "",false, false, false,0.0,0.19,0.0,2.7);
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
    Plot_Histograms(pot_num,sample_name,colors, h_stv_overlay_vec, h_stv_overlay0[i][0],h_stv_overlay0[i][1],h_stv_overlay0[i][2],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv[i], h_stv[i],
                    pad_stv[i], pad0_stv[i], legend_stv[i],channel_legend,
                    stv_ylim[run_num][i], stv_ymin[run_num][i], num_channels, stv_titles[i], path, a1,"", stv[i], "", a1,false, false,false,0.0,0.19);
    h_stv_overlay_vec.clear();
    
    //raquel
    Plot_Histograms(pot_num,sample_name,colors_raquel, h_stv_overlay_raquel_vec, h_stv_overlay0_raquel[i][0],h_stv_overlay0_raquel[i][1],h_stv_overlay0_raquel[i][2],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv_raquel[i], h_stv_raquel[i],
                    pad_stv_raquel[i], pad0_stv_raquel[i], legend_stv_raquel[i],channel_legend_raquel,
                    stv_ylim[run_num][i], stv_ymin[run_num][i], num_channels_raquel, stv_titles[i], path, a1,"", Form("%s_raquel",stv[i]), "",false, false,false,0.0,0.19);
    h_stv_overlay_raquel_vec.clear();
  }



} //end of program
