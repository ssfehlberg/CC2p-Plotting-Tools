#define xsec_cxx
#include "xsec.h"
#include "helper_funcs.h"

void xsec::main(){

  TFile *f_overlay=new TFile("root_files/pelee/Run_all/histograms_pelee_overlay_wgt.root");//overlay histograms.
  TFile *f_dirt=new TFile("root_files/pelee/Run_all/histograms_pelee_dirt_wgt.root");//dirt histograms 
  TFile *f_bnb=new TFile("root_files/pelee/Run_all/histograms_pelee_bnb.root");//bnb histograms  
  TFile *f_ext=new TFile("root_files/pelee/Run_all/histograms_pelee_ext.root");//extbnb histograms 

   TFile *tfile = new TFile("root_files/pelee/Run_all/xsec_extraction.root","RECREATE"); //output root file

  //Grabbing Correct Date and Time and Directory
  ///////////////////////////////////////////////////
  const char* pathname = Form("images/%d%d%d_pelee_Run_all/",Month,Day,Year);
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

  //Stuff for Drawing
  //////////////////////
  TLatex* t = new TLatex(); //latex
  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  t->SetNDC();
  t->SetTextAlign(22);
  char const* pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";
  char const* sample_name="#scale[0.6]{MicroBooNE In-Progress}";
  
  //Calculate the Constants for the XSec
  ///////////////////////////////////////////
  Define_Constants();

  //Get the Histograms and Define the XSec ones to Save
  ///////////////////////////////////////////////////
  Grab_Histograms(f_bnb,f_ext,f_overlay,f_dirt);
   
  //Plotting the Migration Matrices
  //////////////////////////////////
  for(int i=0; i < num_particles_matrices; i++){
    for(int j=0; j < num_particles_matrices_plots; j++){
      canv_particle_matrices[i][j] = new TCanvas(Form("canv_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]),Form("canv_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]),2000,1500);
      gStyle->SetPalette(1);
      h_particle_matrices[i][j]->Draw("colz text"); //will add text later
      h_particle_matrices[i][j]->SetTitle(Form("%s: %s ; True %s; Reco. %s",titles_particles_matrices[i],titles_particle_matrices[j],titles_particle_matrices[j],titles_particle_matrices[j]));
      t->DrawLatex(0.3,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      h_particle_matrices[i][j]->Write();
      canv_particle_matrices[i][j]->Update();
      canv_particle_matrices[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));
      canv_particle_matrices[i][j]->Print(Form("%s%s%s.png",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));   
    }
  }

  for(int i=0; i < num_other_matrices; i++){
    canv_other_matrices[i] = new TCanvas(Form("canv_other_matrices%s",other_matrices[i]),Form("canv_other_matrices%s",other_matrices[i]),2000,1500);
    gStyle->SetPalette(1);
    h_other_matrices[i]->Draw("colz text"); //will add text later
    h_other_matrices[i]->SetTitle(Form("%s ; True %s; Reco. %s",titles_other_matrices[i],titles_other_matrices[i],titles_other_matrices[i]));
    t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
    canv_other_matrices[i]->Update();
    h_other_matrices[i]->Write();
    canv_other_matrices[i]->Print(Form("%s%s.pdf",path.c_str(),other_matrices[i]));
    canv_other_matrices[i]->Print(Form("%s%s.png",path.c_str(),other_matrices[i]));
  } 
  

  //Have to make the efficiency histograms and save them to the file for smearing matrix
  ////////////////////////////////////////////////////////////////////////////////////

  //Momentum theta and phi
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
      /* a[i] = new TLine(xlim_eff[i],0,xlim_eff[i],1);
      a[i]->Draw("same");
      a[i]->SetLineColor(kBlack);
      a[i]->SetLineWidth(4);
      a1[i] = new TLine(xlim_high_eff[i],0,xlim_high_eff[i],1);
      a1[i]->Draw("same");
      a1[i]->SetLineColor(kBlack);
      a1[i]->SetLineWidth(4);
      */
      h_particle_num[i][j]->Write();
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff_4_smearing.png",path.c_str(),particles_eff[i],particles_eff_var[j]));
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff_4_smearing.pdf",path.c_str(),particles_eff[i],particles_eff_var[j])); 
    }
  }

  //all other variables
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
    h_other_eff_num[i]->Write();
    canv_other_eff[i]->Print(Form("%s%s_eff_4_smearing.png",path.c_str(),other_eff[i]));
    canv_other_eff[i]->Print(Form("%s%s_eff_4_smearing.pdf",path.c_str(),other_eff[i]));

  }

  //Creating the Smearing Matrices
  //////////////////////////////////
  for(int i=0; i < num_particles_matrices; i++){
    for(int j=0; j < num_particles_matrices_plots; j++){

      canv_particle_smearing[i][j] =  new TCanvas(Form("canv_particle_smearing%s%s",particles_matrices[i],particles_matrices_plots[j]),Form("canv_particle_smearing%s%s",particles_matrices[i],particles_matrices_plots[j]),2000,1500);
      h_particle_smearing[i][j] = (TH2D*)h_particle_matrices[i][j]->Clone();
      double num_bins_true = double(h_particle_smearing[i][j]->GetXaxis()->GetNbins());
      double num_bins_reco = double(h_particle_smearing[i][j]->GetYaxis()->GetNbins()); 
      
      for(int jj=1; jj < num_bins_reco+1; jj++){
	double sum = 0;
	for(int k=1; k < num_bins_true+1; k++){
	  double bin_content = h_particle_smearing[i][j]->GetBinContent(k,jj);
	  sum += 0 + bin_content;
	}
	for(int k=1; k < num_bins_true+1; k++){
	  double bin_content = h_particle_smearing[i][j]->GetBinContent(k,jj);
	  double eff = h_particle_num[i][j]->GetBinContent(k);
	  double eff_correction;
	  if(eff == 0.0 | eff < 0.025){
	    eff_correction = 0.0;
	  }else{
	    eff_correction = 1.0/eff;
	  }
	  h_particle_smearing[i][j]->SetBinContent(k,jj,(bin_content/sum)*eff_correction); //remember to apply eff_correction
	}
      } //end of loop over true bins

      h_particle_smearing[i][j]->Draw("colz"); //will add text later
      h_particle_smearing[i][j]->SetTitle(Form("%s: %s ; True %s; Reco. %s",titles_particles_matrices[i],titles_particle_matrices[j],titles_particle_matrices[j],titles_particle_matrices[j]));
      t->DrawLatex(0.3,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      canv_particle_smearing[i][j]->Update();
      canv_particle_smearing[i][j]->Print(Form("%s%s%s_smearing.pdf",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));
      canv_particle_smearing[i][j]->Print(Form("%s%s%s_smearing.png",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));        
    }
  }

  for(int i=0; i <num_other_matrices; i++){

    canv_other_smearing[i] = new TCanvas(Form("canv_other_smearing%s",other_matrices[i]),Form("canv_other_smearing%s",other_matrices[i]),2000,1500);
    h_other_smearing[i] = (TH2D*)h_other_matrices[i]->Clone();
    double num_bins_true = double(h_other_smearing[i]->GetXaxis()->GetNbins());
    double num_bins_reco = double(h_other_smearing[i]->GetYaxis()->GetNbins());
  
    for(int jj=1; jj < num_bins_reco+1; jj++){
      double sum = 0;
      for(int k=1; k < num_bins_true+1; k++){
	double bin_content = h_other_smearing[i]->GetBinContent(k,jj);
	sum += 0 + bin_content;
      }
      for(int k=1; k < num_bins_true+1; k++){
	double bin_content = h_other_smearing[i]->GetBinContent(k,jj);
	double eff = h_other_eff_num[i]->GetBinContent(k);
	double eff_correction;
	if(eff == 0.0 | eff < 0.025){
	  eff_correction = 0.0;
	}else{
	  eff_correction = 1.0/eff;
	}
	h_other_smearing[i]->SetBinContent(k,jj,(bin_content/sum)*eff_correction); //remember to apply eff_correction
      }
    } //end of loop over reco bins

    h_other_smearing[i]->Draw("colz text"); //add text later
    h_other_smearing[i]->SetTitle(Form("%s ; True %s; Reco. %s",titles_other_matrices[i],titles_other_matrices[i],titles_other_matrices[i]));
    t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
    canv_other_smearing[i]->Update();
    canv_other_smearing[i]->Print(Form("%s%s_smearing.pdf",path.c_str(),other_matrices[i]));
    canv_other_smearing[i]->Print(Form("%s%s_smearing.png",path.c_str(),other_matrices[i]));
  }
   

  /*  //Trying to Determine Optimum Bin Size
  //////////////////////////////////////
  for(int i=0; i <25; i++){
    const char* x_low = linspaced_char[i];
    double xlow = linspaced[i];
 
    for(int j=0; j < 24; j++){
      const char* x_high = linspaced_char[j+1];
      double xhigh = linspaced[j+1];
      
      if(xlow < xhigh){
	double num_entries = h_test[i][j]->GetEntries();
	std::cout<<"Number of Entries: "<<num_entries<<std::endl;

	if(num_entries != 0){
	  canv_bin[i][j] = new TCanvas(Form("canv_bin_%s-%s",x_low,x_high),Form("canv_bin_%s-%s",x_low,x_high),2000,1500);
	  h_test[i][j]->Draw("hist");
	  h_test[i][j]->Fit("gaus");
	  g[i][j] = (TF1*)h_test[i][j]->GetListOfFunctions()->FindObject("gaus");
	  double sigma = g[i][j]->GetParameter(2);
	  std::cout <<"Value of Sigma: "<< sigma<<std::endl;
	  std::cout<<"Value of 2 Sigma: "<<2*sigma<<std::endl;

	  if((2.0* sigma) <= 0.2 && (2.0*sigma) >= 0.0){
	    std::cout<<Form("Value of x_low: %s",x_low)<<Form(" Value of x_high: %s",x_high)<<std::endl;
	  }//end of if sigma
	
	} //end of entry loop
      } //end of high loop
    } //end of j loop
  } //end of i loop
  
  */

  
  //Close File Before Finishing:
  /////////////////////////////////////
  tfile->Close();
  
} //end of program
