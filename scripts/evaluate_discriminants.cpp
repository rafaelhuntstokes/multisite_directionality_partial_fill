#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DB.hh>
#include <RAT/GeoUtils.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>
#include <string>
#include <RAT/SunUtil.hh>

double multisite_discriminant(std::vector<double> time_residuals, TH1D* pdf_B8_multi, TH1D* pdf_Tl208_multi){

    /*
    Loop through each time residual and find the loglikelihood ratio for B8 and Tl208.
    */

    double loglikelihood  = 0.0;
    int contributing_hits = 0;  // normalise loglikelihood by number of contributing hits 
    for (int ires = 0; ires < time_residuals.size(); ires++){
        
        double residual = time_residuals.at(ires);

        // find probability of being signal and background
        double prob_signal     = pdf_B8_multi->GetBinContent(pdf_B8_multi->FindBin(residual));
        double prob_background = pdf_Tl208_multi->GetBinContent(pdf_Tl208_multi->FindBin(residual));

        // check if we will have divide by zero error or log(0) error
        if (prob_background <= 0 or prob_signal <= 0){
            continue;
        }
        contributing_hits++;

        // likelihood ratio
        double ratio = prob_signal / prob_background;
        ratio = log(ratio);

        loglikelihood += ratio;
    }
    
    // normalise by number of hits
    loglikelihood /= contributing_hits;

    return loglikelihood;
}

double directionality_fitter(TVector3 solar_dir, RAT::DU::Point3D event_position, std::vector<double> pmt_x, std::vector<double> pmt_y, std::vector<double> pmt_z, std::vector<double> time_residuals, TH2D* pdf_dir){
    /*
    Function performs a grid search over two spherical angles, theta and phi to propose an event direction.

    Proposed event direction is evaluated using 2D eneergy dependant directionality PDFs. 

    Highest likelihood event direction used to find cos_theta_sun, which is returned.

    If solar direction fails (no time residuals in -5 --> 15 ns window) -999.9 is returned.
    */

    std::cout << "Fitting direction." << std::endl;
    double pi = M_PI;
    
    // define the angles to grid search over
    std::vector<double> theta;
    std::vector<double> phi;

    // make event position into a TVector3
    TVector3 event_pos(event_position[0], event_position[1], event_position[2]);

    double step_theta = (pi - 0) / 50;
    double step_phi    = (2*pi - 0) / 90;
    for (int i = 0; i < 90; i++){
        double angle = 0 + i * step_phi;
        phi.push_back(angle);
    }
    for (int i = 0; i < 50; i++){
        double angle = 0 + i * step_theta;
        theta.push_back(angle);
    }

    // initial values for best direction and likelihood
    double best_likelihood = 10e30;
    TVector3 best_dir(1,1,1);

    // flag to check we actually have some time residuals inside the directionality window
    bool fit_succeed = false;

    // grid search over these angles
    for (int itheta = 0; itheta < theta.size(); itheta++){
        for (int iphi = 0; iphi < phi.size(); iphi++){
            // log likelihood starting angle set to zero for each proposed direction
            double loglikelihood = 0;

            // use spherical coordinate transform to propose an event direction
            double x = sin(theta.at(itheta)) * cos(phi.at(iphi));
            double y = sin(theta.at(itheta)) * sin(phi.at(iphi));
            double z = cos(theta.at(itheta));

            TVector3 proposed_dir(x, y, z);

            // calculate the likelihood from the pdfs and time residuals for this direction
            for (int ires = 0; ires < time_residuals.size(); ires++){
                double residual = time_residuals.at(ires);

                // check if we fall into the directionality time window
                if (residual < -5 or residual > 15){
                    continue;
                }

                fit_succeed = true;

                // else we get the pmt positions too and calculate the likelihood given by this hit
                TVector3 pmt_position(pmt_x.at(ires), pmt_y.at(ires), pmt_z.at(ires));

                // use pdf to get density
                TVector3 photon_angle = pmt_position - event_pos;
                photon_angle        = photon_angle.Unit();

                double cos_photon_event = photon_angle.Dot(proposed_dir);
                double density          = pdf_dir->GetBinContent(pdf_dir->FindBin(residual, cos_photon_event));
                loglikelihood           -= log(density);
                
            }

            // check if any residuals inside dir time window
            if (fit_succeed == false){
                // return out the function here DO NOT PROCEED
                return -999.9;
            }

            // finished this event direction evaluation --> see if -2log(L) is < best direction LL value
            loglikelihood = 2*loglikelihood;

            if (loglikelihood < best_likelihood){
                best_likelihood = loglikelihood;
                best_dir        = proposed_dir;
            }
        }
    }

    // we've now finished looping through all the proposed event directions, so find cos(theta_sun) with this direction
    double cos_theta = -solar_dir.Dot(best_dir);

    return cos_theta;
}

void evaluate_discriminants(std::string isotope, int run_number, float fv_cut, float z_cut, std::string input_fpath, std::string output_fpath){
    /*
    Function loads in a run_by_run test .root ratds file and loops over all the events. 
    
    Every event that passes the FV, Z, Recon and energy cuts is filtered by energy into a pre-defined energy bin.
    
    The multisite (dLog(L), Fisher, IQR) and directionality discriminants are calculted using previously created PDFs a given 
    energy.

    A run_by_run ntuple file is returned containing the [energy, position, discriminants] of each event and for each energy.
    */


    // define ratds input file location and output TTree file locations
    std::string input_fname  = input_fpath + "/PartialScint" + isotope + "_ScintRun_r" + std::to_string(run_number) + "*.root";
    std::string output_fname = output_fpath + "/" + isotope + "/" + std::to_string(run_number) + ".root"; 
    
    // load the PDF files for each isotope
    std::string working_directory = "/data/snoplus3/hunt-stokes/multisite_directionality_0p6gL_refactored";
    std::string pdf_B8_fname      = working_directory + "/run_by_run_pdf/B8_Solar_Nue/total.root";
    std::string pdf_Tl208_fname   = working_directory + "/run_by_run_pdf/Tl208/total.root";
    TFile *PDF_B8 = new TFile(pdf_B8_fname.c_str());
    TFile *PDF_Tl208 = new TFile(pdf_Tl208_fname.c_str());
    // get the directionality histograms from the files (use B8 PDF only!)

    /*
    BE CAREFUL ! if you put in the wrong name of the histograms, dynamic_cast returns a null_ptr.
    No error is reported until you try to access the histograms...
    */                                             
    TH2D *dir_pdf_2p5_5p0    = dynamic_cast<TH2D*>(PDF_B8->Get("directionality_2.5_5.0"));
    TH2D *dir_pdf_2p5_3p125  = dynamic_cast<TH2D*>(PDF_B8->Get("directionality_2.5_3.125"));
    TH2D *dir_pdf_3p125_3p75 = dynamic_cast<TH2D*>(PDF_B8->Get("directionality_3.125_3.75"));
    TH2D *dir_pdf_3p75_4p375 = dynamic_cast<TH2D*>(PDF_B8->Get("directionality_3.75_4.375"));
    TH2D *dir_pdf_4p375_5p0  = dynamic_cast<TH2D*>(PDF_B8->Get("directionality_4.375_5.0"));
    if (dir_pdf_2p5_5p0 == false){
        std::cout << "Broken!" << std::endl;
    }

    // get the multisite histograms from the files (for both isotopes!)
    TH1D *multi_pdf_B8_2p5_5p0       = dynamic_cast<TH1D*>(PDF_B8->Get("multi_2.5_5.0"));
    TH1D *multi_pdf_Tl208_2p5_5p0    = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_2.5_5.0"));
    TH1D *multi_pdf_B8_2p5_3p125     = dynamic_cast<TH1D*>(PDF_B8->Get("multi_2.5_3.125"));
    TH1D *multi_pdf_Tl208_2p5_3p125  = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_2.5_3.125"));
    TH1D *multi_pdf_B8_3p125_3p75    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_3.125_3.75"));
    TH1D *multi_pdf_Tl208_3p125_3p75 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_3.125_3.75"));
    TH1D *multi_pdf_B8_3p75_4p375    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_3.75_4.375"));
    TH1D *multi_pdf_Tl208_3p75_4p375 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_3.75_4.375"));
    TH1D *multi_pdf_B8_4p375_5p0     = dynamic_cast<TH1D*>(PDF_B8->Get("multi_4.375_5.0"));
    TH1D *multi_pdf_Tl208_4p375_5p0  = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_4.375_5.0"));

    // create the output ntuple with an output TTree containing discriminants for each energy PDF
    TFile *output_file    = new TFile(output_fname.c_str(), "RECREATE");

    TTree *info_2p5_5p0    = new TTree("2p5_5p0", "2p5_5p0");
    TTree *info_2p5_3p125  = new TTree("2p5_3p125", "2p5_3p125");
    TTree *info_3p125_3p75 = new TTree("3p125_3p75", "3p125_3p75");
    TTree *info_3p75_4p375 = new TTree("3p75_4p375", "3p75_4p375");
    TTree *info_4p375_5p0  = new TTree("4p375_5p0", "4p375_5p0");

    // define the variables to fill branches with
    double energy, x, y, z, dlogL, fisher, IQR, cos_theta_sun;
    
    // add the relevant branches to each TTRee
    info_2p5_5p0->Branch("energy", &energy);
    info_2p5_5p0->Branch("x", &x);
    info_2p5_5p0->Branch("y", &y);
    info_2p5_5p0->Branch("z", &z);
    info_2p5_5p0->Branch("dlogL", &dlogL);
    info_2p5_5p0->Branch("fisher", &fisher);
    info_2p5_5p0->Branch("IQR", &IQR);
    info_2p5_5p0->Branch("cos_theta_sun", &cos_theta_sun);

    info_2p5_3p125->Branch("energy", &energy);
    info_2p5_3p125->Branch("x", &x);
    info_2p5_3p125->Branch("y", &y);
    info_2p5_3p125->Branch("z", &z);
    info_2p5_3p125->Branch("dlogL", &dlogL);
    info_2p5_3p125->Branch("fisher", &fisher);
    info_2p5_3p125->Branch("IQR", &IQR);
    info_2p5_3p125->Branch("cos_theta_sun", &cos_theta_sun);

    info_3p125_3p75->Branch("energy", &energy);
    info_3p125_3p75->Branch("x", &x);
    info_3p125_3p75->Branch("y", &y);
    info_3p125_3p75->Branch("z", &z);
    info_3p125_3p75->Branch("dlogL", &dlogL);
    info_3p125_3p75->Branch("fisher", &fisher);
    info_3p125_3p75->Branch("IQR", &IQR);
    info_3p125_3p75->Branch("cos_theta_sun", &cos_theta_sun);

    info_3p75_4p375->Branch("energy", &energy);
    info_3p75_4p375->Branch("x", &x);
    info_3p75_4p375->Branch("y", &y);
    info_3p75_4p375->Branch("z", &z);
    info_3p75_4p375->Branch("dlogL", &dlogL);
    info_3p75_4p375->Branch("fisher", &fisher);
    info_3p75_4p375->Branch("IQR", &IQR);
    info_3p75_4p375->Branch("cos_theta_sun", &cos_theta_sun);

    info_4p375_5p0->Branch("energy", &energy);
    info_4p375_5p0->Branch("x", &x);
    info_4p375_5p0->Branch("y", &y);
    info_4p375_5p0->Branch("z", &z);
    info_4p375_5p0->Branch("dlogL", &dlogL);
    info_4p375_5p0->Branch("fisher", &fisher);
    info_4p375_5p0->Branch("IQR", &IQR);
    info_4p375_5p0->Branch("cos_theta_sun", &cos_theta_sun);

    /*
    From here we load the RATDS file for this run and evaluate the event by event discriminants!
    */
    std::cout << "Loading File: " << input_fname << std::endl;
    RAT::DU::DSReader dsReader(input_fname.c_str());
    std::cout << "Loaded File." << std::endl;

    int num_entries = dsReader.GetEntryCount();
    std::cout << "Found " << num_entries << " events." << std::endl;

    // setup time residual calculator, Point3D and get PMT info
    RAT::DU::TimeResidualCalculator timeResCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
    size_t fPSUPSystemId = RAT::DU::Point3D::GetSystemId("innerPMT");
    size_t av_system_id  = RAT::DU::Point3D::GetSystemId("av");
    const RAT::DU::PMTInfo& pmt_info = RAT::DU::Utility::Get()->GetPMTInfo();

    // begin loop over every entry
    for (int ientry = 0; ientry < num_entries; ientry++){

        const RAT::DS::Entry& rDS = dsReader.GetEntry(ientry);

        // check event triggered detector
        if (rDS.GetEVCount() == 0){
            // did not trigger!
            continue;
        }

        // get primary trigger
        const RAT::DS::EV& rEV = rDS.GetEV(0);

        // check the reconstruction succeeded
        RAT::DS::FitVertex fVertex = rEV.GetDefaultFitVertex();
      
        // check that a position/time fit exists                 
        if( !fVertex.ContainsPosition() ) continue;
        if( !fVertex.ValidPosition() ) continue;
        if( !fVertex.ContainsTime() ) continue;
        if( !fVertex.ValidTime() ) continue;

        // get the reconstructed variables
        RAT::DU::Point3D event_position_recon(fPSUPSystemId);
        event_position_recon.SetXYZ(fPSUPSystemId, fVertex.GetPosition());

        energy            = fVertex.GetEnergy();
        double recon_time = fVertex.GetTime();

        // convert to AV coordinates and apply selection cuts
        event_position_recon.SetCoordinateSystem(av_system_id);
        x = event_position_recon.X();
        y = event_position_recon.Y();
        z = event_position_recon.Z();
        if (event_position_recon.Mag() > fv_cut){
            // failed FV cut!
            continue;
        }

        if (z < z_cut){
            // not in the scintillator!
            continue;
        }

        // convert back to PSUP Coords for time residual calculation
        event_position_recon.SetCoordinateSystem(fPSUPSystemId);

        // check that event energy fits into the large domain
        if (energy < 2.5 or energy > 5.0){
            // outside energy ROI!
            continue;
        }

        // passed all event selection cuts! Now we evaluate the discriminants

        // get the solar direction for this run
        RAT::DS::UniversalTime event_time = rEV.GetUniversalTime();
        TVector3 solar_dir                = RAT::SunDirection(event_time.GetDays(), event_time.GetSeconds(), event_time.GetNanoSeconds()).Unit();

        // check result exists for the directionality
        bool dir_flag = false;

        // work out which PDF to use for this energy
        TH2D *pdf_dir;
        TH1D *pdf_B8_multi;
        TH1D *pdf_Tl208_multi;
        if (energy >= 2.5 and energy < 3.125){
            pdf_dir         = dir_pdf_2p5_3p125;
            pdf_B8_multi    = multi_pdf_B8_2p5_3p125;
            pdf_Tl208_multi = multi_pdf_Tl208_2p5_3p125;  
        }
        if (energy >= 3.125 and energy < 3.75){
            pdf_dir         = dir_pdf_3p125_3p75;
            pdf_B8_multi    = multi_pdf_B8_3p125_3p75;
            pdf_Tl208_multi = multi_pdf_Tl208_3p125_3p75;  
        }
        if (energy >= 3.75 and energy < 4.375){
            pdf_dir         = dir_pdf_3p75_4p375;
            pdf_B8_multi    = multi_pdf_B8_3p75_4p375;
            pdf_Tl208_multi = multi_pdf_Tl208_3p75_4p375;  
        }
        if (energy >= 4.375 and energy < 5.0){
            pdf_dir         = dir_pdf_4p375_5p0;
            pdf_B8_multi    = multi_pdf_B8_4p375_5p0;
            pdf_Tl208_multi = multi_pdf_Tl208_4p375_5p0;  
        }

        // each discriminant involves the time residuals, so evaluate and store (temporarily) the time residuals and PMT positions
        std::vector<double> time_residuals;
        std::vector<double> pmt_x;
        std::vector<double> pmt_y;
        std::vector<double> pmt_z;
        double tres_recon; 

        const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
        for (size_t ipmt = 0; ipmt < calibratedPMTs.GetCount(); ipmt++){
            const RAT::DS::PMTCal& pmtCal             = calibratedPMTs.GetPMT(ipmt);
            RAT::DU::Point3D pmt_point(fPSUPSystemId, pmt_info.GetPosition(pmtCal.GetID()));
            const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
            // check the hit is 'good'
            if (PMTCalStatus.GetHitStatus(pmtCal) != 0){
                continue;
            }

            tres_recon = timeResCalc.CalcTimeResidual(pmtCal, event_position_recon, recon_time, true);
            
            time_residuals.push_back(tres_recon);
            pmt_x.push_back(pmt_point[0]);
            pmt_y.push_back(pmt_point[1]);
            pmt_z.push_back(pmt_point[2]);
        }
    
        // now have all the time residuals and PMT positions saved in vector --> can fit the direction and find multisite dLog(l)
        std::cout << "Fitting direction for event " << ientry << std::endl;
        cos_theta_sun = directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, pdf_dir); // return false if no in time residuals
        std::cout << "Cos(theta_sun) = " << cos_theta_sun << std::endl;
        // only calculate multisite discriminants if cos_theta_sun fit succeeded
        if (cos_theta_sun != -999.9){
            std::cout << "Evaluating multisite discriminant." << std::endl;
            dlogL         = multisite_discriminant(time_residuals, pdf_B8_multi, pdf_Tl208_multi);
            IQR           = 0; // not implemented yet
            fisher        = 0; // not implemented yet

            // fill the ntuples

            if (energy >= 2.5 and energy < 3.125){
                info_2p5_3p125->Fill();
            }
            if (energy >= 3.125 and energy < 3.75){
                info_3p125_3p75->Fill();
            }
            if (energy >= 3.75 and energy < 4.375){
                info_3p75_4p375->Fill();
            }
            if (energy >= 4.375 and energy < 5.0){
                info_4p375_5p0->Fill();
            }
        }

        // now repeat it using the full 2.5 --> 5.0 PDF
        cos_theta_sun = directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, dir_pdf_2p5_5p0);
        // only calculate multisite discriminants if cos_theta_sun fit succeeded
        if (cos_theta_sun != -999.9){
            dlogL         = multisite_discriminant(time_residuals, multi_pdf_B8_2p5_5p0, multi_pdf_Tl208_2p5_5p0);
            IQR           = 0; // not implemented yet
            fisher        = 0; // not implemented yet

            info_2p5_5p0->Fill();
        }
    }

    std::cout << "Writing results to file." << std::endl;
    // finished loop over all events and TTrees are filled --> write them to the output file
    info_2p5_5p0->Write();
    info_2p5_3p125->Write();
    info_3p125_3p75->Write();
    info_3p75_4p375->Write();
    info_4p375_5p0->Write();
    output_file->Close();
}

int main(int argc, char* argv[]){

    std::string isotope      = argv[1];
    int run_number           = std::stoi(argv[2]);
    double fv_cut            = std::stod(argv[3]);
    double z_cut             = std::stod(argv[4]);
    std::string input_fpath  = argv[5];
    std::string output_fpath = argv[6];

    evaluate_discriminants(isotope, run_number, fv_cut, z_cut, input_fpath, output_fpath);
}