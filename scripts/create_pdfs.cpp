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
#include <TFile.h>
#include <string>

void create_pdf(std::string isotope, int run_number, float fv_cut, float z_cut, std::string input_fpath, std::string output_fpath){
    /*
    Create an output .root file containing PDF information given an input run_by_run MC .root file.

    INPUTS: isotope          - for naming and saving the file + loading correct input file, str
            run_number       - for loading correct MC file, int
            fv_cut           - fiducial volume cut in mm to apply to events, float
            z_cut            - z cut in mm, float
            output_fpath     - where to save the output to, str 
            input_fpath      - where to load the input from, str

    OUTPUT: .root file containing PDF level information for a given isotope and event cuts
    */

    // define the run_by_run input MC file
    std::string input_fname = input_fpath + "/PartialScint" + isotope + "_ScintRun_r" + std::to_string(run_number) + "*.root";

    // create the output file
    std::string output_path = output_fpath + "/" + isotope + "/" + std::to_string(run_number) + ".root";
    TFile *output           = new TFile(output_path.c_str(), "RECREATE");

    // populate the output file with a TTree and branches to save relevant information
    TTree *hit_level             = new TTree("PDF_hit_by_hit_Information", "PDF hit_by_hit Information");
    TTree *event_level           = new TTree("PDF_event_by_event_Information", "PDF event_by_event Information");
    // define type of variables
    double posX, posY, posZ, energy, tres_recon, tres_mc, cerenkov_angle;
    
    // pass by reference to each branch
    hit_level->Branch("tRes_recon_full", &tres_recon);
    hit_level->Branch("tRes_mc_full", &tres_mc);
    hit_level->Branch("cerenkov_angle", &cerenkov_angle);

    event_level->Branch("x", &posX);
    event_level->Branch("y", &posY);
    event_level->Branch("z", &posZ);
    event_level->Branch("energy", &energy);

    // create PDF histograms for multisite and directionality for each energy bin
    TH2D* dir_pdf_2p5_5p0    = new TH2D("directionality_2.5_5.0", "", 401, -100.5, 300.5, 20, -1.0, 1.0);
    TH2D* dir_pdf_2p5_3p125  = new TH2D("directionality_2.5_3.125", "", 401, -100.5, 300.5, 20, -1.0, 1.0);
    TH2D* dir_pdf_3p125_3p75 = new TH2D("directionality_3.125_3.75", "", 401, -100.5, 300.5, 20, -1.0, 1.0);
    TH2D* dir_pdf_3p75_4p375 = new TH2D("directionality_3.75_4.375", "", 401, -100.5, 300.5, 20, -1.0, 1.0);
    TH2D* dir_pdf_4p375_5p0  = new TH2D("directionality_4.375_5.0", "", 401, -100.5, 300.5, 20, -1.0, 1.0);

    TH1D* multi_pdf_2p5_5p0    = new TH1D("multi_2.5_5.0", "", 401, -100.5, 300.5);
    TH1D* multi_pdf_2p5_3p125  = new TH1D("multi_2.5_3.125", "", 401, -100.5, 300.5);
    TH1D* multi_pdf_3p125_3p75 = new TH1D("multi_3.125_3.75", "", 401, -100.5, 300.5);
    TH1D* multi_pdf_3p75_4p375 = new TH1D("multi_3.75_4.375", "", 401, -100.5, 300.5);
    TH1D* multi_pdf_4p375_5p0  = new TH1D("multi_4.375_5.0", "", 401, -100.5, 300.5);

    // push back these histograms into vectors --> this is important when deciding which energy histo to fill
    // efficienctly inside the time residual loop
    std::vector<TH1D*> multi_pdf_vec;
    std::vector<TH2D*> dir_pdf_vec;

    multi_pdf_vec.push_back(multi_pdf_2p5_5p0);
    multi_pdf_vec.push_back(multi_pdf_2p5_3p125);
    multi_pdf_vec.push_back(multi_pdf_3p125_3p75);
    multi_pdf_vec.push_back(multi_pdf_3p75_4p375);
    multi_pdf_vec.push_back(multi_pdf_4p375_5p0);

    dir_pdf_vec.push_back(dir_pdf_2p5_5p0);
    dir_pdf_vec.push_back(dir_pdf_2p5_3p125);
    dir_pdf_vec.push_back(dir_pdf_3p125_3p75);
    dir_pdf_vec.push_back(dir_pdf_3p75_4p375);
    dir_pdf_vec.push_back(dir_pdf_4p375_5p0);

    // load the input MC file
    std::cout << "Reading file: " << input_fname << std::endl;
    RAT::DU::DSReader dsReader(input_fname.c_str());
    std:: cout << "Read file." <<std::endl;
    // Get AV offset for this run                   
    std::vector<double> av_pos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
    double av_offset = av_pos[2];
    std::cout << "AV offset = " << av_offset << std::endl;

    int num_entries = dsReader.GetEntryCount();
    std::cout << "Extracting PDF information from "  << num_entries << " events." << std::endl;

    // setup time residual calculator
    RAT::DU::TimeResidualCalculator timeResCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();                                                             
    
    // point3D use now ENFORCED
    size_t fPSUPSystemId = RAT::DU::Point3D::GetSystemId("innerPMT");
    size_t av_system_id  = RAT::DU::Point3D::GetSystemId("av");
    
    // Get PMT info (positions etc)
    const RAT::DU::PMTInfo& pmt_info = RAT::DU::Utility::Get()->GetPMTInfo();

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

        energy = fVertex.GetEnergy();
        double event_time   = fVertex.GetTime();

        // convert to AV coordinates and apply selection cuts
        event_position_recon.SetCoordinateSystem(av_system_id);
        posX       = event_position_recon.X();
        posY       = event_position_recon.Y();
        posZ       = event_position_recon.Z();
        std::cout << "posZ before: " << posZ;
        if (event_position_recon.Mag() > fv_cut){
            // failed FV cut!
            continue;
        }

        if (posZ < z_cut){
            // not in the scintillator!
            continue;
        }

        // convert back to PSUP Coords for time residual calculation
        event_position_recon.SetCoordinateSystem(fPSUPSystemId);
        std::cout << "posZ after: " << posZ;
        // check that event energy fits into the large domain
        if (energy < 2.5 or energy > 5.0){
            // outside energy ROI!
            continue;
        }

        // find what energy bin pdf to fill
        int idx_to_fill_pdf;
        if (energy >= 2.5 and energy < 3.125){
            idx_to_fill_pdf = 1;
        }
        else if (energy >= 3.125 and energy < 3.75){
            idx_to_fill_pdf = 2;
        }
        else if (energy >= 3.75 and energy < 4.375){
            idx_to_fill_pdf = 3;
        }
        else{
            idx_to_fill_pdf = 4;
        }

        /* 
        we now have an event that falls inside the energy ROI and within the FV!
        Now calculate the time residual and cerenkov angle for the multisite and directionality PDF.
        */

        // obtain the event position and direction from MC
        RAT::DU::Point3D event_position_mc(fPSUPSystemId);
        event_position_mc.SetXYZ(fPSUPSystemId, rDS.GetMC().GetMCParticle(0).GetPosition());
        TVector3 event_dir         = rDS.GetMC().GetMCParticle(0).GetMomentum().Unit();

        // fill event level info
        event_level->Fill();
        const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
        for (size_t ipmt = 0; ipmt < calibratedPMTs.GetCount(); ipmt++){
            const RAT::DS::PMTCal& pmtCal             = calibratedPMTs.GetPMT(ipmt);
            RAT::DU::Point3D pmt_point(fPSUPSystemId, pmt_info.GetPosition(pmtCal.GetID()));
            const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
            // check the hit is 'good'
            if (PMTCalStatus.GetHitStatus(pmtCal) != 0){
                continue;
            }

            // calculate time residuals using MC and reconstructed position (but both using recon time (?!))
            tres_recon = timeResCalc.CalcTimeResidual(pmtCal, event_position_recon, event_time, true);
            tres_mc    = timeResCalc.CalcTimeResidual(pmtCal, event_position_mc, event_time, true);
            
            // work out the cerenkov angle for this hit
            // TVector3 event_pos(event_position_mc[0], event_position_mc[1], event_position_mc[2]);
            TVector3 photon_vec = pmt_point - event_position_mc;
            TVector3 photon_dir = photon_vec.Unit();
            cerenkov_angle      = photon_dir.Dot(event_dir);

            // write to the full branch tracking information for the PDF (x, y, z, energy, tres, cerenkov angle)
            hit_level->Fill();

            // fill the histograms
            dir_pdf_2p5_5p0->Fill(tres_mc, cerenkov_angle);
            dir_pdf_vec.at(idx_to_fill_pdf)->Fill(tres_mc, cerenkov_angle);
            
            multi_pdf_2p5_5p0->Fill(tres_recon);
            multi_pdf_vec.at(idx_to_fill_pdf)->Fill(tres_recon);
        }
    }

    // now write the histograms and TTree to the file
    dir_pdf_vec.at(0)->Write();
    dir_pdf_2p5_3p125->Write();
    dir_pdf_3p125_3p75->Write();
    dir_pdf_3p75_4p375->Write();
    dir_pdf_4p375_5p0->Write();

    multi_pdf_2p5_5p0->Write();
    multi_pdf_2p5_3p125->Write();
    multi_pdf_3p125_3p75->Write();
    multi_pdf_3p75_4p375->Write();
    multi_pdf_4p375_5p0->Write();

    hit_level->Write();
    event_level->Write();
}

int main(int argc, char *argv[]){

    std::string isotope = argv[1];;
    int run_number      = std::stoi(argv[2]);
    float fv_cut        = std::stof(argv[3]);
    float z_cut         = std::stof(argv[4]);
    std::string input   = argv[5];
    std::string output  = argv[6];

    create_pdf(isotope, run_number, fv_cut, z_cut, input, output);

    return 0;
}