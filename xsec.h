#ifndef xsec_h
#define xsec_h

class xsec{

public:
  virtual void main();
  virtual void Grab_Histograms(TFile* f_bnb,TFile* f_ext,TFile* f_overlay,TFile* f_dirt);
  virtual void Define_Constants();
  virtual void Plot_XSec();
  
private:

  //Stuff for date and time
  /////////////////////////
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;

  //Migration Matrices
  //////////////////////////////////////////////
  static const int num_particles_matrices_plots = 3;
  const char* particles_matrices_plots[num_particles_matrices_plots] ={"_mom","_costheta","_phi"};
  static const int num_particles_matrices = 5;
  const char* particles_matrices[num_particles_matrices] = {"_muon_all","_muon_contained","_muon_uncontainied","_lead_proton","_recoil_proton"};
  TH2D* h_particle_matrices[num_particles_matrices][num_particles_matrices_plots];
  TCanvas* canv_particle_matrices[num_particles_matrices][num_particles_matrices_plots];
  const char* titles_particle_matrices[num_particles_matrices_plots] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};
  const char* titles_particles_matrices[num_particles_matrices] = {"All Muons","Contained Muons","Uncontained Muons","Leading Proton","Recoil Proton"};

  static const int num_other_matrices =8;
  const char* other_matrices[num_other_matrices] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading","_opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_nu_E"};
  TH2D* h_other_matrices[num_other_matrices];
  TCanvas* canv_other_matrices[num_other_matrices];
  const char* titles_other_matrices[num_other_matrices] = {"cos(#gamma_{Lab})","cos(#gamma_{COM})","cos(Opening Angle Muon and Leading)","cos(Opening Angle Muon and Both Protons)","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)", "#delta #phi_{T} (Deg.)","Neutrino Energy"};

  //Efficiency Plots. Needed to make the Smearing Matrices
  ///////////////////////////////////////////////////////

  //Momentum, theta and phi of particles
  static const int num_particles_eff_plots = 3;
  const char* particles_eff_var[num_particles_eff_plots] = {"_mom","_costheta","_phi"};
  const char* particles_eff_var_titles[num_particles_eff_plots] = {"True Momentum (GeV/c)","True cos(#theta)","True #phi (Rad.)"};
  static const int num_particles_eff = 5;
  const char* particles_eff[num_particles_eff] = {"_muon_all","_muon_contained","_muon_uncontainied","_lead_proton","_recoil_proton"};
  const char* particles_eff_titles[num_particles_eff] = {"All Muons","Contained Muons","Uncontained Muons","Leading Proton","Recoil Proton"};
  TH1D* h_particle_num[num_particles_eff][num_particles_eff_plots];
  TH1D* h_particle_denom[num_particles_eff][num_particles_eff_plots];     
  TCanvas* canv_particle_eff[num_particles_eff][num_particles_eff_plots];
  /*double xlim_eff[num_eff] = {0.1,0.1,0.1,0.25,0.25};
  double xlim_high_eff[num_eff] = {-9999.,-9999.,-9999.,1.2,1.2};
  TLine* a[num_eff];
  TLine* a1[num_eff];
  */
  
  //efficiency plots of other variables to determine potential xsec candidates
  static const int num_other_eff = 8;
  const char* other_eff[num_other_eff] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading","_opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_nu_E"};
  const char* other_eff_titles[num_other_eff] = {"cos(Opening Angle Protons_{Lab})","cos(Opening Angle Protons_{COM}","cos(Opening Angle Muon and Leading_{Lab})","cos(Opening Angle Muon and Both Protons_{Lab})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)","True Neutrino Energy (GeV)"};
  TH1D* h_other_eff_num[num_other_eff];
  TH1D* h_other_eff_denom[num_other_eff];
  TCanvas* canv_other_eff[num_other_eff];

  //Smearing Matrices
  /////////////////////
  TH2D* h_particle_smearing[num_particles_matrices][num_particles_matrices_plots];
  TCanvas* canv_particle_smearing[num_particles_matrices][num_particles_matrices_plots];
  TH2D* h_other_smearing[num_other_matrices];
  TCanvas* canv_other_smearing[num_other_matrices];

  //Grabbin Histograms to Determine Bin Size
  /////////////////////////////////////////
  TH1D* h_test[25][25];
  std::vector<double> linspaced{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  const char* linspaced_char[25] = {"0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9","1_0","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9","2_0","2_1","2_2","2_3","2_4","2_5"};

  static const int bleh = 25;
  TCanvas* canv_bin[bleh][bleh];
  TF1* g[bleh][bleh];

  //Things for the XSec
  //////////////////////////
  static const int num_particles = 5;
  static const int num_particle_plots = 3;
  const char* particles[num_particles] = {"_muon_all","_muon_contained","_muon_uncontained","leading_proton","_recoil_proton"};
  const char* particle_plots[num_particle_plots] = {"_mom","_cos_theta","phi"};
  TH1D* h_particle[num_particles][num_particle_plots];
  TH1D* h_particle_bnb[num_particles][num_particle_plots];
  TH1D* h_particle_ext[num_particles][num_particle_plots];
  TH1D* h_particle_overlay[num_particles][num_particle_plots];
  TH1D* h_particle_dirt[num_particles][num_particle_plots];

  static const int num_plots = 8;
  const char* plots[num_plots] = {"_opening_angle_lab","opening_angle_com","opening_angle_mu_leading","opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_neutrino_energy"};
  TH1D* h[num_plots];
  TH1D* h_bnb[num_plots];
  TH1D* h_ext[num_plots];
  TH1D* h_overlay[num_plots];
  TH1D* h_dirt[num_plots];

  //List to Write Histograms to
  ////////////////////////////
  std::vector<TH2D*> h_2D_list;
  
};
#endif
#ifdef xsec_cxx

void xsec::Grab_Histograms(TFile* f_bnb,TFile* f_ext,TFile* f_overlay,TFile* f_dirt){

  //Migration Matrices and creating new Smearing Matrices
  ////////////////////////////
  for(int i=0; i < num_particles_matrices; i++){
    for(int j=0; j < num_particles_matrices_plots; j++){
      h_particle_matrices[i][j] = (TH2D*)f_overlay->Get(Form("h_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]));
      //h_particle_smearing[i][j] = new (TH2D*)Form();
      //h_2D_list.push_back(h_particle_smearing[i][j]);
    }
  }
  
  for(int i=0; i < num_other_matrices; i++){
    h_other_matrices[i] = (TH2D*)f_overlay->Get(Form("h_other_matrices%s",other_matrices[i]));
    //h_other_smearing[i] = new (TH2D*);
    //h_2D_list.push_back(h_other_smearing[i]);
  }


  //Efficiency Plots
  ///////////////////////////
  for(int i=0; i < num_particles_eff; i++){
    for(int j=0; j < num_particles_eff_plots; j++){
      h_particle_num[i][j] = (TH1D*)f_overlay->Get(Form("h_particle_num%s%s",particles_eff[i],particles_eff_var[j]));
      h_particle_denom[i][j] = (TH1D*)f_overlay->Get(Form("h_particle_denom%s%s",particles_eff[i],particles_eff_var[j]));
    }
  }

  for(int i = 0; i < num_other_eff; i++){
    h_other_eff_num[i] = (TH1D*)f_overlay->Get(Form("h_other_eff_num%s",other_eff[i]));
    h_other_eff_denom[i] = (TH1D*)f_overlay->Get(Form("h_other_eff_denom%s",other_eff[i]));
  }


  //Determine Bin Size Plots
  ///////////////////////////
  for(int i=0; i <25; i++){
    const char* x_low = linspaced_char[i];
    for(int j=0; j < 24; j++){
      const char* x_high = linspaced_char[j+1];
      if(x_low < x_high){
	h_test[i][j] = (TH1D*)f_overlay->Get(Form("h_test_%s-%s",x_low,x_high));
      }
    }
  }


  
  /* //XSec Histograms
     /////////////////////////
  for(int i=0; i < num_particles; i++){
    for(int j=0; j < num_particle_plots; j++){
      h_particle[i][j] = new TH1D*(Form(),Form(),);
      h_particle_bnb[i][j] = f_bnb->Get(Form(""));
      h_particle_ext[i][j] = f_ext->Get(Form(""));
      h_particle_overlay[i][j] = f_overlay->Get(Form("_total"));
      h_particle_dirt[i][j] = f_dirt->Get(Form(""));
    }
  }

  for(int i=0; i < num_plots; i++){
    h[i] = new TH1D*(Form(),Form(),);
    h_bnb[i] = f_bnb->Get(Form(""));
    h_ext[i] = f_ext->Get(Form(""));
    h_overlay[i] = f_overlay->Get(Form("_total"));
    h_dirt[i] = f_dirt->Get(Form(""));
  }
  */ 
} //end of grab histograms

void xsec::Define_Constants(){

  //Calculate the Constants for the XSec
  ///////////////////////////////////////////
  double rho_Ar = 1.3836; //density of argon in g*cm^-3
  double V = (256.35 - 20.0) * (233.0 - 20.0) * (1036.8 - 20.0);//volume of detector cm: x*y*z
  double N_A = 6.023e+23; //avogadro's number
  double m_mol = 39.95; //mass of argon in g*mol^-1
  double N_targets = (rho_Ar*V*N_A)/m_mol; //number of target nuclei
  std::cout<<"Value of Volume: "<<V<<std::endl;
  std::cout<<"Value of N_targets assuming 10cm FV border: "<<N_targets<<std::endl;
  
  double pot = 6.79e+20;
  double flux = pot*(1.16e+11/1.6e+20); //fraction is MCC8 CCInclusive flux/POT
  std::cout<<"Value of Flux: "<<flux<<""<<std::endl;

  double bin_wdith;

} //end of define constants


/*void xsec::Get_Bin_Sizes(TH2D* hist){

  std::cout<<"Starting Process to Find Bin Sizes"<<std::endl;
  
  const int num_bins = int(hist->GetXaxis()->GetNbins());
  std::cout<<"Value of num_bins:"<<num_bins<<std::endl;
  std::vector<double> xlim{double(hist->GetXaxis()->GetXmin()),double(hist->GetXaxis()->GetXmax())};
  std::cout<<"Min X Value: "<<xlim[0]<<" Max X Value: "<<xlim[1]<<std::endl;


  const int num_iterations = 50; //how many times are we going to make the plot
  double delta = (xlim[0]-xlim[1])/(num_iterations-1);
  std::vector<double> line[num_iterations];
  for(int i=0; i < num_iterations; i++){
    line[i] = xlim[0] + (i*delta);
  }
  std::cout<<"Number of Iterations: "<<num_iterations<<std::endl;
  std::cout<<"Length of Vector: "<<line.size()<<std::endl;
 
  TCanvas* canv[num_bins];
  TH1D* h_reco[num_bins];

  std::cout<<"Starting to Loop:"<<std::endl;  
  for(int i=0; i < num_bins; i++){
    h_reco[i] = new TH1D(Form("h_reco_%d",i),Form("h_reco_%d",i),num_bins,xlim[0],xlim[1]);

    if(){
      h_reco[i]->Fill(hist->GetBinContent(i,j));
    }

    canv[i] = new TCanvas(Form("canv_%d",i),Form("canv_%d",i),2000,1500);
    h_reco[i]->Draw("hist");
    h_reco[i]->Fit("gaus");
    h_reco[i]->Draw();   
  } 
  
}//end of GetBinSizes()
*/
void xsec::Plot_XSec(){


  
}
#endif
