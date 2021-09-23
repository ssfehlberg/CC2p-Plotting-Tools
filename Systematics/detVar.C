#define detVar_cxx
#include "../paul_tol_colors.hpp"
#include <iostream>
#include <ctime>
#include <string>

void area_normalize(TH1D* h_CV, TH1D* hist){

  //This sets the integral to 1
  int n_bins = hist->GetNbinsX();
  double hist_integral = hist->Integral();

  for(int i=1; i < n_bins+1; i++){
    double hist_bin_content = hist->GetBinContent(i);
    double hist_value = hist_bin_content/hist_integral;
    hist->SetBinContent(i,hist_value);
  }
  
  /*
  //This Scales to CV
  double CV_integral = h_CV->GetEntries();
  double hist_integral = hist->GetEntries();
  double scale_factor = CV_integral / hist_integral;
  hist->Scale(scale_factor);

  std::cout<<"CV_Integral: "<<CV_integral<<std::endl;
  std::cout<<"hist_Integral: "<<hist_integral<<std::endl;
  std::cout<<"Scale Factor: "<<scale_factor<<std::endl;
  */
  
}


void detVar(){

  //Stuff for Drawing
  //////////////////////
  TLatex* t = new TLatex(); //latex
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  t->SetNDC();
  t->SetTextAlign(22);
  char const* pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";
  char const* sample_name="#scale[0.6]{MicroBooNE In-Progress}";
  tolcols::init();
  Color_t colors[] = {0,9031, 9025, 9030, 9024, 9029, kViolet-8, 9028, 9032, 9026, kGray+2,kOrange+7};

  //Which binning do we want
  ///////////////////////////
  bool use_xsec_binning = false;
  const char* which_binning;
  
  if(use_xsec_binning){
    which_binning = "pelee_xsec";
  } else{
    which_binning = "pelee";
  }
  
  //Detector Variation Samples
  /////////////////////////////
  static const int number_of_files = 12;
  const char* samples[number_of_files] = {"detVar_CV", "detVar_dEdx","detVar_LY_Attenuation",
					  "detVar_LY_Down","detVar_LY_Rayleigh","detVar_Overlay",
					  "detVar_Recombination","detVar_SCE","detVar_ThetaXZ", 
					  "detVar_ThetaYZ","detVar_X","detVar_YZ"};
  const char* sample_titles[number_of_files] = {"Central Value","dE/dx","LY Attenuation",
						"LY Down","LY Rayleigh","Overlay",
						"Recombination","SCE","ThetaXZ",
						"ThetaYZ","X","YZ"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};
  
  //Muon
  TH1D* h_muon[number_of_files][num_var];
  TH1D* h_muon0[num_var];
  double ymax_muon[num_var] = {0.05,0.08,0.05};
  double ymin_muon[num_var] = {-0.05,-0.08,-0.05};
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
 
  //Lead
  TH1D* h_leading[number_of_files][num_var];
  TH1D* h_leading0[num_var];
  double ymax_leading[num_var] = {0.05,0.08,0.05};
  double ymin_leading[num_var] = {-0.05,-0.08,-0.05};
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  
  //Recoil
  TH1D* h_recoil[number_of_files][num_var];
  TH1D* h_recoil0[num_var];
  double ymax_recoil[num_var] = {0.05,0.08,0.05};
  double ymin_recoil[num_var] = {-0.05,-0.08,-0.05};
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  
  //Other Variables
  /////////////////////////
  static const int num_other_var = 9;
  const	char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading",
                                          "_opening_angle_mu_both","_delta_PT","_delta_alphaT",
					  "_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})","cos(#gamma_{#mu,p_{L}})",
						 "cos(#gamma_{#mu,p_{L} + p_{R}})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)",
						 "#delta #phi_{T} (Deg.)","Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"}; 
  TH1D* h_other[number_of_files][num_other_var];
  TH1D* h_other0[num_other_var];
  double ymax_other[num_other_var] = {0.05,0.08,0.05,0.05,0.05,0.05,0.05,0.05,0.05};
  double ymin_other[num_other_var] = {-0.05,-0.08,-0.05,-0.05,-0.05,-0.05,-0.05,-0.05,-0.05};
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];

  
  //Empty file to grab everything
  TFile* file;
  
  //Grab the Files
  ///////////////////////////
  for(int i=0; i < number_of_files; i++){
    file = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/Systematics/root_files/%s/detVar/histograms_%s_%s.root",which_binning,which_binning,samples[i]));

    for(int j=0; j < num_var; j++){
      h_muon[i][j] = (TH1D*)file->Get(Form("h_muon%s_total",var0[j]));
      h_leading[i][j] = (TH1D*)file->Get(Form("h_leading%s_total",var0[j]));
      h_recoil[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total",var0[j]));
    } //end of loop over particle variables

    for(int j=0; j < num_other_var; j++){
      h_other[i][j] = (TH1D*)file->Get(Form("h%s_total",other_var[j]));
    } //end of loop over other variables

  }//end of loop over number of files


  ///////////////////////////////
  //Now to draw all of this garbage
  ////////////////////////////////
  

  for(int j=0; j < num_var; j++){

    //Muon
    ////////////////////
    canv_muon[j] = new TCanvas(Form("canv_muon%s",var0[j]),Form("canv_muon%s",var0[j]),2000,1500);

    for(int k=0; k < number_of_files; k++){
      area_normalize(h_muon[0][j],h_muon[k][j]);
    }
    
    h_muon0[j] = (TH1D*)h_muon[0][j]->Clone();
    h_muon[0][j]->Add(h_muon0[j],-1);
    h_muon[0][j]->Draw("hist");
    h_muon[0][j]->SetLineColor(kBlack);
    h_muon[0][j]->SetLineWidth(3);
    h_muon[0][j]->SetTitle(Form("Muon %s",var_titles[j]));
    h_muon[0][j]->SetXTitle(Form("Muon %s",var_titles[j]));
    h_muon[0][j]->GetXaxis()->SetTitleSize(40);
    h_muon[0][j]->GetXaxis()->SetTitleFont(43);
    h_muon[0][j]->GetXaxis()->SetTitleOffset(1.5);
    h_muon[0][j]->GetXaxis()->SetLabelFont(43);
    h_muon[0][j]->GetXaxis()->SetLabelSize(30);
      
    h_muon[0][j]->SetYTitle("Fractional Uncertainty (%)");
    h_muon[0][j]->GetYaxis()->SetTitleSize(40);
    h_muon[0][j]->GetYaxis()->SetTitleFont(43);
    h_muon[0][j]->GetYaxis()->SetTitleOffset(1.5);
    h_muon[0][j]->GetYaxis()->SetLabelFont(43);
    h_muon[0][j]->GetYaxis()->SetLabelSize(30);
  
    h_muon[0][j]->SetTitle(Form("Muon %s",var_titles[j]));
    h_muon[0][j]->SetMaximum(1.0);//ymax_muon[j]);
    h_muon[0][j]->SetMinimum(-1.0);//ymin_muon[j]);

    legend_muon[j] = new TLegend(0.48, 0.65, 0.89, 0.89);
    legend_muon[j]->SetNColumns(3);
    legend_muon[j]->AddEntry(h_muon[0][j],Form("%s",sample_titles[0]),"L");
   
    for(int k=1; k < number_of_files; k++){
      h_muon[k][j]->Add(h_muon0[j],-1);
      h_muon[k][j]->Draw("hist same");
      h_muon[k][j]->SetLineColor(colors[k]);
      h_muon[k][j]->SetLineWidth(3);
      legend_muon[j]->AddEntry(h_muon[k][j],Form("%s",sample_titles[k]),"L");
    }

    legend_muon[j]->Draw("same");
    canv_muon[j]->Print(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/Systematics/images/%s/detVar/_muon%s.png",which_binning,var0[j]));
    
    //Leading Proton
    /////////////////

    canv_leading[j] = new TCanvas(Form("canv_leading%s",var0[j]),Form("canv_leading%s",var0[j]),2000,1500);

    for(int k=0; k < number_of_files; k++){
      area_normalize(h_leading[0][j],h_leading[k][j]);
    }
    
    h_leading0[j] = (TH1D*)h_leading[0][j]->Clone();
    h_leading[0][j]->Add(h_leading0[j],-1);
    h_leading[0][j]->Draw("hist");
    h_leading[0][j]->SetLineColor(kBlack);
    h_leading[0][j]->SetLineWidth(3);
    h_leading[0][j]->SetTitle(Form("Leading %s",var_titles[j]));
    h_leading[0][j]->SetXTitle(Form("Leading %s",var_titles[j]));
    h_leading[0][j]->GetXaxis()->SetTitleSize(40);
    h_leading[0][j]->GetXaxis()->SetTitleFont(43);
    h_leading[0][j]->GetXaxis()->SetTitleOffset(1.5);
    h_leading[0][j]->GetXaxis()->SetLabelFont(43);
    h_leading[0][j]->GetXaxis()->SetLabelSize(30);
      
    h_leading[0][j]->SetYTitle("Fractional Uncertainty (%)");
    h_leading[0][j]->GetYaxis()->SetTitleSize(40);
    h_leading[0][j]->GetYaxis()->SetTitleFont(43);
    h_leading[0][j]->GetYaxis()->SetTitleOffset(1.5);
    h_leading[0][j]->GetYaxis()->SetLabelFont(43);
    h_leading[0][j]->GetYaxis()->SetLabelSize(30);
  
    h_leading[0][j]->SetTitle(Form("Leading %s",var_titles[j]));
    h_leading[0][j]->SetMaximum(ymax_leading[j]);
    h_leading[0][j]->SetMinimum(ymin_leading[j]);

    legend_leading[j] = new TLegend(0.48, 0.65, 0.89, 0.89);
    legend_leading[j]->SetNColumns(3);
    legend_leading[j]->AddEntry(h_leading[0][j],Form("%s",sample_titles[0]),"L");
   
    for(int k=1; k < number_of_files; k++){
      h_leading[k][j]->Add(h_leading0[j],-1);
      h_leading[k][j]->Draw("hist same");
      h_leading[k][j]->SetLineColor(colors[k]);
      h_leading[k][j]->SetLineWidth(3);
      legend_leading[j]->AddEntry(h_leading[k][j],Form("%s",sample_titles[k]),"L");
    }

    legend_leading[j]->Draw("same");
    canv_leading[j]->Print(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/Systematics/images/%s/detVar/_leading%s.png",which_binning,var0[j]));

    //Recoil Proton
    ///////////////

    canv_recoil[j] = new TCanvas(Form("canv_recoil%s",var0[j]),Form("canv_recoil%s",var0[j]),2000,1500);

    for(int k=0; k < number_of_files; k++){
      area_normalize(h_recoil[0][j],h_recoil[k][j]);
    }
    
    h_recoil0[j] = (TH1D*)h_recoil[0][j]->Clone();
    h_recoil[0][j]->Add(h_recoil0[j],-1);
    h_recoil[0][j]->Draw("hist");
    h_recoil[0][j]->SetLineColor(kBlack);
    h_recoil[0][j]->SetLineWidth(3);
    h_recoil[0][j]->SetTitle(Form("Recoil %s",var_titles[j]));
    h_recoil[0][j]->SetXTitle(Form("Recoil %s",var_titles[j]));
    h_recoil[0][j]->GetXaxis()->SetTitleSize(40);
    h_recoil[0][j]->GetXaxis()->SetTitleFont(43);
    h_recoil[0][j]->GetXaxis()->SetTitleOffset(1.5);
    h_recoil[0][j]->GetXaxis()->SetLabelFont(43);
    h_recoil[0][j]->GetXaxis()->SetLabelSize(30);
      
    h_recoil[0][j]->SetYTitle("Fractional Uncertainty (%)");
    h_recoil[0][j]->GetYaxis()->SetTitleSize(40);
    h_recoil[0][j]->GetYaxis()->SetTitleFont(43);
    h_recoil[0][j]->GetYaxis()->SetTitleOffset(1.5);
    h_recoil[0][j]->GetYaxis()->SetLabelFont(43);
    h_recoil[0][j]->GetYaxis()->SetLabelSize(30);
  
    h_recoil[0][j]->SetTitle(Form("Recoil %s",var_titles[j]));
    h_recoil[0][j]->SetMaximum(ymax_recoil[j]);
    h_recoil[0][j]->SetMinimum(ymin_recoil[j]);

    legend_recoil[j] = new TLegend(0.48, 0.65, 0.89, 0.89);
    legend_recoil[j]->SetNColumns(3);
    legend_recoil[j]->AddEntry(h_recoil[0][j],Form("%s",sample_titles[0]),"L");
   
    for(int k=1; k < number_of_files; k++){
      h_recoil[k][j]->Add(h_recoil0[j],-1);
      h_recoil[k][j]->Draw("hist same");
      h_recoil[k][j]->SetLineColor(colors[k]);
      h_recoil[k][j]->SetLineWidth(3);
      legend_recoil[j]->AddEntry(h_recoil[k][j],Form("%s",sample_titles[k]),"L");
    }

    legend_recoil[j]->Draw("same");
    canv_recoil[j]->Print(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/Systematics/images/%s/detVar/_recoil%s.png",which_binning,var0[j]));

  } //end of loop over num_var
  

  //Now for the Other Variables
  ///////////////////////////
  for(int j=0; j < num_other_var; j++){
    canv_other[j] = new TCanvas(Form("canv_%s",other_var[j]),Form("canv_%s",other_var[j]),2000,1500);

    for(int k=0; k < number_of_files; k++){
      area_normalize(h_other[0][j],h_other[k][j]);
    }
    
    h_other0[j] = (TH1D*)h_other[0][j]->Clone();
    h_other[0][j]->Add(h_other0[j],-1);
    h_other[0][j]->Draw("hist");
    h_other[0][j]->SetLineColor(kBlack);
    h_other[0][j]->SetLineWidth(3);
    h_other[0][j]->SetTitle(Form("%s",other_var_titles[j]));
    h_other[0][j]->SetXTitle(Form("%s",other_var_titles[j]));
    h_other[0][j]->GetXaxis()->SetTitleSize(40);
    h_other[0][j]->GetXaxis()->SetTitleFont(43);
    h_other[0][j]->GetXaxis()->SetTitleOffset(1.5);
    h_other[0][j]->GetXaxis()->SetLabelFont(43);
    h_other[0][j]->GetXaxis()->SetLabelSize(30);
      
    h_other[0][j]->SetYTitle("Fractional Uncertainty (%)");
    h_other[0][j]->GetYaxis()->SetTitleSize(40);
    h_other[0][j]->GetYaxis()->SetTitleFont(43);
    h_other[0][j]->GetYaxis()->SetTitleOffset(1.5);
    h_other[0][j]->GetYaxis()->SetLabelFont(43);
    h_other[0][j]->GetYaxis()->SetLabelSize(30);
  
    h_other[0][j]->SetTitle(Form("%s",other_var_titles[j])); 
    h_other[0][j]->SetMaximum(ymax_other[j]);
    h_other[0][j]->SetMinimum(ymin_other[j]);

    legend_other[j] = new TLegend(0.48, 0.65, 0.89, 0.89);
    legend_other[j]->SetNColumns(3);
    legend_other[j]->AddEntry(h_other[0][j],Form("%s",sample_titles[0]),"L");
   
    for(int k=1; k < number_of_files; k++){
      h_other[k][j]->Add(h_other0[j],-1);
      h_other[k][j]->Draw("hist same");
      h_other[k][j]->SetLineColor(colors[k]);
      h_other[k][j]->SetLineWidth(3);
      legend_other[j]->AddEntry(h_other[k][j],Form("%s",sample_titles[k]),"L");
    }

    legend_other[j]->Draw("same");
    canv_other[j]->Print(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/Systematics/images/%s/detVar/%s.png",which_binning,other_var[j]));
}
  
  
} //end of program
