#ifndef __IMAGECAL_CXX__
#define __IMAGECAL_CXX__

#include "ImageCal.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "../Reco3D/AStarUtils.h"
#include <string>

namespace larcv {

	static ImageCalProcessFactory __global_ImageCalProcessFactory__;
	ImageCal::ImageCal(const std::string name)
	: ProcessBase(name)
	{}
	void ImageCal::configure(const PSet& cfg) {

		_img2d_producer		= cfg.get<std::string>("Image2DProducer");
		//	_par_pix_producer	= cfg.get<std::string>("ParPixelProducer");
	}
	void ImageCal::initialize() {
		std::string file_loc = "/uboone/app/users/ewhall/dllee_unified/LArCV/app/CalibMaps/CalibrationMaps_MCC9.root";
		//LoadCalibrationMaps(file_loc);
	
//		TCanvas *halp = new TCanvas("halp","halp",400,750);
//		halp->Divide(1,3);
//		std::cout << "Draw Map 1" << std::endl;
//		halp->cd(1); hCalibrationMap[0]->Draw("BOX2");
//		std::cout << "Draw Map 2" << std::endl;
//		halp->cd(2); hCalibrationMap[1]->Draw("BOX2");
//		std::cout << "Draw Map 3" << std::endl;
//		halp->cd(3); hCalibrationMap[2]->Draw("BOX2");	
		LARCV_INFO() << "[ImageCal]" << std::endl;

		//assert(!_spline_file.empty());
		std::string filename;
		//_storage {nullptr};
		track_file = "/uboone/app/users/ewhall/test_files/tracker_reco_63120148.root";

		std::cout << filename << std::endl;
		
		if (_foutll.empty()) throw larbys("specify larlite file output name");

      		_storage.set_io_mode(larlite::storage_manager::kWRITE);
       		_storage.set_out_filename(_foutll);
		xmin = .5; 	xmax = 255.;
		ymin = -116.;	ymax = 116.;
		zmin = 0.9; 	zmax = 1036.;
		
		point_vector = CreatePoints(.25);
		std::cout << "Problem Here?" << std::endl;
		TrackInformation();
		std::cout << "Number of Events in Track_Info: " << track_info.size();
		EventCount = 0;
		std::cout << "Nope!" << std::endl;
		std::cout << track_info.size() << std::endl;
		double wireRange[3] = {3456,3456,3456};
		double tickRange = 8448 - 2400;
		
		std::cout << "Reserving meta size" << std::endl;
		Standard_meta_v.reserve(3);
		//Calibration_meta_v.reserve(3);
		for (size_t iPlane{0}; iPlane <3; iPlane++) {
			std::cout << "Wirerange: " << wireRange[iPlane] << std::endl;
			std::cout << "Tickrange: " << tickRange << std::endl;
			//Standard_meta_v[iPlane] = larcv::ImageMeta(,
			Standard_meta_v[iPlane]	= larcv::ImageMeta(wireRange[iPlane],tickRange,(int)(tickRange)/6,(int)(wireRange[iPlane]),0,8448,iPlane);
			std::cout << "standard meta " << iPlane << std::endl;
			Calibration_meta_v.push_back(larcv::Image2D(Standard_meta_v[iPlane]));
			Track_image_v.push_back(larcv::Image2D(Standard_meta_v[iPlane]));
			std::cout << "calib meta " << iPlane << std::endl;
			//Calibration_meta_v[iPlane].paint(1);
		}
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			Calibration_meta_v[iPlane].paint(1);
		}
		std::cout << "minx: " << Calibration_meta_v[0].meta().min_x() << std::endl;
		std::cout << "maxx: " << Calibration_meta_v[0].meta().max_x() << std::endl;
		std::cout << "miny: " << Calibration_meta_v[0].meta().min_y() << std::endl;
		std::cout << "maxy: " << Calibration_meta_v[0].meta().max_y() << std::endl;
		std::cout << "pixel width: " << Calibration_meta_v[0].meta().pixel_width() << std::endl;
		std::cout << "pixel height: " << Calibration_meta_v[0].meta().pixel_height() << std::endl;
		
		LoadCalibrationMaps(file_loc);
		CreateCalibration2D();

		// for (size_t x_pos{0}; x_pos < point_vector[0].size(); x_pos++) {
		// 	for (size_t y_pos{0}; y_pos < point_vector[1].size(); y_pos++) {
		// 		for (size_t z_pos{0}; z_pos < point_vector[2].size(); z_pos++) {
		// 			for (size_t iPlane{0}; iPlane < 3; iPlane++) {
		// 				TVector3 point(point_vector[0][x_pos], point_vector[1][y_pos], point_vector[2][z_pos]);
		// 				int bin = hCalibrationMap[iPlane]->FindBin(point.X(),point.Y(),point.Z());		
		// 				//std::cout << "Point at : (" << point_vector[0][x_pos] << ", " << point_vector[1][y_pos] << ", " << point_vector[2][z_pos] << ") is in bin number: " << bin << std::endl;	
		// 				double calFac = hCalibrationMap[iPlane]->GetBinContent(bin);
		// 				ProjectTo3D(Standard_meta_v[iPlane],
		// 						point_vector[0][x_pos],
		// 						point_vector[1][y_pos],
		// 						point_vector[2][z_pos],
		// 						0, iPlane, xpixel, ypixel);
		// 				ypixel = (Calibration_meta_v[iPlane].meta().rows() - ypixel);
		// 				//std:: cout << "(ypix, xpix): (" << ypixel << ", " << xpixel << "), CalFac: " << calFac << std::endl;
		// 				if (calFac !=0) { Calibration_meta_v[iPlane].set_pixel(ypixel,xpixel,calFac); }
		// 				//if (calFac == 0) { Calibration_meta_v[iPlane].set_pixel(ypixel,xpixel,1); }
						
												
		// 				//if (point_vector[0][x_pos] < 15) { Calibration_meta_v[iPlane].set_pixel(ypixel,xpixel,1); }
		// 				//else if (point_vector[0][x_pos] > 240) { Calibration_meta_v[iPlane].set_pixel(ypixel,xpixel,1); }
		// 				//else { Calibration_meta_v[iPlane].set_pixel(ypixel,xpixel,1); }
		// 			}
		// 		}
		// 	}
		// }
//		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
//			TGraph *g = new TGraph();
//			track_points.push_back(*g); 
//		} 
		std::vector<std::vector<double>> pixel_x (3);
		std::vector<std::vector<double>> pixel_y (3);
		for (size_t idx{0}; idx < track_info[0].size(); idx++) {
			double x_p = track_info[0][idx][0];
			double y_p = track_info[0][idx][1];
			double z_p = track_info[0][idx][2];
			for (size_t iPlane{0}; iPlane < 3; iPlane++) {
				ProjectTo3D(Standard_meta_v[iPlane],
						x_p, y_p, z_p, 0, iPlane, xpixel, ypixel);
				//track_points[iPlane].Set(track_points[iPlane].GetN()+1);
				//track_points[iPlane].SetPoint(track_points[iPlane].GetN()-1,ypixel,xpixel);
				//Track_image_v[iPlane].set_pixel(ypixel, xpixel, 1);
				pixel_x[iPlane].push_back(xpixel);
				pixel_y[iPlane].push_back((Calibration_meta_v[iPlane].meta().rows()*6)-(ypixel*6)+2400);
			}
		}
		std::cout << "Number of Rows: " << Calibration_meta_v[0].meta().rows();
		//TGraph *gTrack[3]; 
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			gTracks[iPlane] = new TGraph(pixel_x[iPlane].size(), &pixel_x[iPlane][0], &pixel_y[iPlane][0]);
			TGraph *g = new TGraph(pixel_x[iPlane].size(), &pixel_y[iPlane][0], &pixel_x[iPlane][0]);
			track_points.push_back(*g);
		}
		TCanvas *help = new TCanvas("help","help",1000,750);
		help->Divide(3,1);
		help->cd(1); gTracks[0]->Draw("AL");
		help->cd(2); gTracks[1]->Draw("AL");
		help->cd(3); gTracks[2]->Draw("AL");
		
		std::cout << "Draw ROI" << std::endl;
		_run = 1;
		_subrun = 1;
		_event = 1;
		_entry = 1;
		DrawROI(Calibration_meta_v, Track_image_v);
		
        if(!_storage.open()) {
            LARCV_CRITICAL() << "ERROR, larlite output file could not open" << std::endl;
            throw larbys("die");
		}
	}

	bool ImageCal::process(IOManager &mgr) {

		//ClearEvent();
		std::cout << std::endl;
		std::cout << "============================================" << std::endl;
		std::cout << "Entry " << mgr.current_entry() << " / " << mgr.get_n_entries() << std::endl;
		std::cout << "============================================" << std::endl;
		gStyle->SetOptStat();
	
		std::cout << "Event Count is at: " << EventCount << std::endl;
		//EventCount++;

		auto ev_img_v = (EventImage2D*)mgr.get_data(kProductImage2D,_img2d_producer);
		std::cout << ev_img_v->Image2DArray().size() << std::endl;
		//std::cout << "above" << std::endl;
		_run	= (int) ev_img_v->run();
		_subrun = (int) ev_img_v->subrun();
		_event 	= (int) ev_img_v->event();
		_entry 	= (int) mgr.current_entry();

		auto full_adc_img_v = &(ev_img_v->Image2DArray());
	//	std::cout << full_adc_img_v->size() << std::endl;

		// std::vector<larcv::ImageMeta> Full_meta_v(3);
		// std::vector<larcv::Image2D> Tagged_Image(3);
		std::vector<larcv::Image2D> Full_image_v(3);
		std::vector<larcv::Image2D> Edited_image_v(3);
		//double wireRange[3] = {3456,3456,3456};
		//double tickRange = 8448 - 2400 ;
	//	std::cout << "Entering loop 1!" << std::endl;
		int filled {0};

		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			//Full_meta_v[iPlane]		= larcv::ImageMeta(wireRange[iPlane],tickRange,(int)(tickRange)/6,(int)(wireRange[iPlane]),0,tickRange);
			Full_image_v[iPlane] 	= larcv::Image2D(Standard_meta_v[iPlane]);
			Edited_image_v[iPlane]	= larcv::Image2D(Standard_meta_v[iPlane]);
			if(full_adc_img_v->size() == 3)Full_image_v[iPlane].overlay( (*full_adc_img_v)[iPlane]);
		//	Edited_image_v[iPlane] = Full_image_v[iPlane]*Calibration_meta_v[iPlane];
			//Full_image_v[iPlane].paint(0);
		}
		
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			for (size_t xpix {0}; xpix < Full_image_v[0].meta().cols(); xpix++) { 
				for (size_t ypix {0}; ypix < Full_image_v[1].meta().rows(); ypix++) {
					double new_pix_val = Full_image_v[iPlane].pixel(ypix,xpix)*Calibration_meta_v[iPlane].pixel(ypix,xpix);
					Edited_image_v[iPlane].set_pixel(ypix,xpix,new_pix_val);
				}
			}
		}
		std::cout << Standard_meta_v[0].dump() << std::endl;
		std::cout << (*full_adc_img_v)[0].meta().dump() << std::endl;
		//std::cout: << Full_meta_v[0].dump() << std::endl;
	//	std::cout << Full_image_v[1].size() << std::endl;
		//std::cout << Full_meta_v[1].size() << std::endl;	
//		for (size_t x_pos {0}; x_pos < point_vector[0].size(); x_pos++) {
//			for (size_t y_pos {0}; y_pos < point_vector[1].size(); y_pos++) {
//				for (size_t z_pos {0}; z_pos < point_vector[2].size(); z_pos++) {
//					int bin = hCalibrationMap[0]->FindBin(point_vector[0][x_pos],point_vector[1][y_pos],point_vector[2][z_pos]);
//					for (size_t iPlane {0}; iPlane < 3; iPlane++) {
//						ProjectTo3D(Full_meta_v[iPlane],														point_vector[0][x_pos],														point_vector[1][y_pos],
//								point_vector[2][z_pos],
//								0, iPlane, xpixel, ypixel);
//						double calFac = hCalibrationMap[iPlane]->GetBinContent(bin);
						
						//std::cout << xpixel << " " << ypixel << std::endl;
						//std::cout << point_vector[0][x_pos] << ", " << point_vector[1][y_pos] << ", " << point_vector[2][z_pos] << std::endl;
					//	std::cout << Full_image_v[iPlane].pixel(ypixel,xpixel) << std::endl;
//						Edited_image_v[iPlane].set_pixel(ypixel, xpixel, Full_image_v[iPlane].pixel(ypixel,xpixel)*calFac);
//						if (Edited_image_v[iPlane].pixel(ypixel,xpixel) != 0) { filled++; }
//					}	
//				}
//			}
//		}
//		std::cout << "filled bins " << filled << std::endl;
//		std::cout << "zero bins" << std::endl;	
//		std::vector<int> zero_bins {0,0,0};
//		for (size_t x {500}; x < 1000; x++) {
//			for (size_t y {2500}; y < 3000; y++) {
//				for (size_t iPlane {0}; iPlane < 3; iPlane++) {
//					if (Full_image_v[iPlane].pixel(x,y) == 0) { zero_bins[iPlane]++; }
//				}
//			}
//		}
//		std::cout << "Number of zero bins in Plane 0 = " << zero_bins[0] << std::endl;
//		std::cout << "Number of zero bins in Plane 1 = " << zero_bins[1] << std::endl;
//		std::cout << "Number of zero bins in Plane 2 = " << zero_bins[2] << std::endl;

		std::cout << "Hello!" << std::endl;
		DrawROI(Full_image_v,Edited_image_v);
		std::cout << "Goodbye!" << std::endl;
		EventCount++;

		return true;
	}

	void ImageCal::finalize() {

	}

	void ImageCal::CreateCalibration2D(){

		point_vector = CreatePoints(.25); //step size for iterating through 3D space in detector 

		TH2D* hCalibrationMeta[3];
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			Calibration_meta_v[iPlane].paint(1);
			hCalibrationMeta[iPlane] = new TH2D(Form("hCalibrationMeta_%zu,iPlane"),Form("hCalibrationMeta_%zu",iPlane),1008,3456);
		}


		for (size_t x_pos{0}; x_pos < point_vector[0].size(); x_pos++) {
			for (size_t y_pos{0}; y_pos < point_vector[1].size(); y_pos++) {
				for (size_t z_pos{0}; z_pos < point_vector[2].size(); z_pos++) {
					for (size_t iPlane{0}; iPlane < 3; iPlane++) {
						TVector3 point(point_vector[0][x_pos], point_vector[1][y_pos], point_vector[2][z_pos]);
						int bin = hCalibrationMap[iPlane]->FindBin(point.X(),point.Y(),point.Z());
						double calFac = hCalibrationMap[iPlane]->GetBinContent(bin);
						ProjectTo3D(Standard_meta_v[iPlane],
							point.X(), point.Y(), point.Z(),
							0, iPlane, xpixel, ypixel);
						ypixel = Calibration_meta_v[iPlane].meta().rows() - ypixel;
						if (calFac != 0) { 
							Calibration_meta_v[iPlane].set_pixel(ypixel,xpixel,calFac); 
							hCalibrationMeta[iPlane]->SetBinContent(ypixel, xpixel, calFac);
						}

					}
				}
			}
		}

		hCalibrationMeta[0]->GetBinContent

		TFile f("CalibrationImage2D.root","RECREATE");
		hCalibrationMeta[0]->Write();
		hCalibrationMeta[1]->Write();
		hCalibrationMeta[2]->Write();

		TH2D

		
	}

	void ImageCal::LoadCalibrationMaps(std::string file_loc) {

		TFile *fCalibration = TFile::Open(Form("%s",file_loc.c_str()),"READ");

		if (!(fCalibration->IsOpen())) {std::cout << "Error: Could not open Calibration File." << std::endl;}
		else {std::cout << "Calibrtion File opened." << std::endl;}

		hCalibrationMap[0] = (TH3D*)fCalibration->Get("hImageCalibrationMap_00");
		hCalibrationMap[1] = (TH3D*)fCalibration->Get("hImageCalibrationMap_01");
		hCalibrationMap[2] = (TH3D*)fCalibration->Get("hImageCalibrationMap_02");

		for (int iPlane {0}; iPlane < 3; iPlane++) {
			if (hCalibrationMap[iPlane]==nullptr) {std::cout << "Error: Calibration Map for Plane " << iPlane << " was not found!" << std::endl;}
			else { std::cout << "Calibration Map for Plane " << iPlane << " was loaded!" << std::endl;}		}

//		TCanvas *cMaps = new TCanvas("cMaps","cMaps",400,750);
//		cMaps->Divide(1,3);
//		std::cout << "Draw Map 1" << std::endl;
//		cMaps->cd(1); hCalibrationMap[0]->Draw("BOX2");
//		std::cout << "Draw Map 2" << std::endl;
//		cMaps->cd(2); hCalibrationMap[1]->Draw("BOX2");
//		std::cout << "Draw Map 3" << std::endl;
//		cMaps->cd(3); hCalibrationMap[2]->Draw("BOX2");
	}

  void ImageCal::DrawROI(std::vector<larcv::Image2D> Full_image_v, std::vector<larcv::Image2D> Edit_image_v) {
		std::cout << "DrawROI()!!" << std::endl;
		TH2D *hImage[3];
		TH2D *hEdit[3];
		// TGraph *gStartNend[3];
		// TGraph *gStart[3];
		std::cout << "Start Loop!" << std::endl;
		std::cout << Full_image_v.size() <<std::endl;
		
		std::vector<int> zero_bins {0,0,0};
		std::vector<int> one_bins {0,0,0};
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
	//	  std::cout << "In Loop!" << iPlane << std::endl;
	//	  std::cout << "Number of cols" << Full_image_v[iPlane].meta().cols() << std::endl;
	//	  std::cout << "Number of rows" << Full_image_v[iPlane].meta().rows() << std::endl;
			hImage[iPlane] = new TH2D(Form("hImage_%d_%d_%d_%zu",_run,_subrun,_event,iPlane),
				";wire;time (sample)",
				Full_image_v[iPlane].meta().cols(),
				Full_image_v[iPlane].meta().tl().x,
				Full_image_v[iPlane].meta().tl().x+Full_image_v[iPlane].meta().width(),
				Full_image_v[iPlane].meta().rows(),
				Full_image_v[iPlane].meta().br().y,
				Full_image_v[iPlane].meta().br().y+Full_image_v[iPlane].meta().height());
			hEdit[iPlane] = new TH2D(Form("hEdit_%d_%d_%d_%zu",_run,_subrun,_event,iPlane),
				";wire;time (sample)",
				Edit_image_v[iPlane].meta().cols(),
				Edit_image_v[iPlane].meta().tl().x,
				Edit_image_v[iPlane].meta().tl().x+Edit_image_v[iPlane].meta().width(),
				Edit_image_v[iPlane].meta().rows(),
				Edit_image_v[iPlane].meta().br().y,
				Edit_image_v[iPlane].meta().br().y+Edit_image_v[iPlane].meta().height());
	//		std::cout << "Loop rows and cols!" << std::endl;
			for (size_t icol{0}; icol < Full_image_v[iPlane].meta().cols(); icol++) {
			  //std::cout << "In cols" << std::endl;
			  for (size_t irow{0}; irow < Full_image_v[iPlane].meta().rows(); irow++) {
				  //std::cout << "In rows" << std::endl;
					hImage[iPlane]->SetBinContent(icol+1, irow+1, Full_image_v[iPlane].pixel(irow,icol));
					hEdit[iPlane]->SetBinContent(icol+1, irow+1, Edit_image_v[iPlane].pixel(irow,icol));
					//if (icol > 500 && icol < 1000) {
					//	if (irow > 4000/6 && irow < 5000/6) {
					//		if ( Full_image_v[iPlane].pixel(irow,icol) == 0 ) {zero_bins[iPlane]++;}
					//		if ( Full_image_v[iPlane].pixel(irow,icol) == 1 ) {one_bins[iPlane]++;}
			//			}
			//		}
				}
			}
			// gStartNend[iPlane] = new TGraph();
			// gStart[iPlane] = new TGraph();
			// double x_pixel_st, y_pixel_st, x_pixel_end, y_pixel_end;
			// double Xstart, Ystart, Xend, Yend;
		}
		//std::cout << zero_bins[0] << " " << zero_bins[1] << " " << zero_bins[2] << std::endl;
		//std::cout << one_bins[0] << " " << one_bins[1] << " " << one_bins[2] << std::endl;
		//std::cout << "Make Canvas!" << std::endl;
		std::vector<std::vector<double>> pixel_x (3);
		std::vector<std::vector<double>> pixel_y (3);
		for (size_t idx{0}; idx < track_info[EventCount].size(); idx++) {
			double x_p = track_info[EventCount][idx][0];
			double y_p = track_info[EventCount][idx][1];
			double z_p = track_info[EventCount][idx][2];
			for (size_t iPlane{0}; iPlane < 3; iPlane++) {
				ProjectTo3D(Standard_meta_v[iPlane],
						x_p, y_p, z_p, 0, iPlane, xpixel, ypixel);
				//track_points[iPlane].Set(track_points[iPlane].GetN()+1);
				//track_points[iPlane].SetPoint(track_points[iPlane].GetN()-1,ypixel,xpixel);
				//Track_image_v[iPlane].set_pixel(ypixel, xpixel, 1);
				pixel_x[iPlane].push_back(xpixel);
				pixel_y[iPlane].push_back((Calibration_meta_v[iPlane].meta().rows()*6)-(ypixel*6)+2400);
			}
		}
		std::cout << "Number of Rows: " << Calibration_meta_v[0].meta().rows();
		//TGraph *gTrack[3]; 
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			gTracks[iPlane] = new TGraph(pixel_x[iPlane].size(), &pixel_x[iPlane][0], &pixel_y[iPlane][0]);
			TGraph *g = new TGraph(pixel_x[iPlane].size(), &pixel_y[iPlane][0], &pixel_x[iPlane][0]);
			track_points.push_back(*g);
		}
		TCanvas *c = new TCanvas(Form("ROI_%d_%d_%d",_run,_subrun,_event),Form("ROI_%d_%d_%d",_run,_subrun,_event),1000,850);
		c->Divide(3,5);
		std::cout << "Paint the canvas" << std::endl;
		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
		  std::cout << "switching to " << iPlane << std::endl;
			c->cd(iPlane+1);
			std::cout << "drawing canvas " << iPlane << std::endl;
			hImage[iPlane]->Draw("colz");
			gTracks[iPlane]->SetMarkerStyle(1);
			gTracks[iPlane]->SetMarkerColor(6);
			//gTracks[iPlane]->SetLineWidth(5);
			//gTracks[iPlane]->SetLineColor(2);
			gTracks[iPlane]->Draw("P SAME");
			c->cd(iPlane+4);
			//gTracks[iPlane]->Draw("AP");
			hEdit[iPlane]->Draw("colz");

		//	track_points[iPlane].SetMarkerStyle(7);
		//	track_points[iPlane].SetLineColor(2);
			c->cd(iPlane+7);
			gTracks[iPlane]->Draw("AP");
			c->cd(iPlane+10);
			gTracks[iPlane]->Draw("AP");
			//hImage[iPlane]->GetZaxis()->SetRangeUser(0,200);
			hImage[iPlane]->Draw("colz same");
			c->cd(iPlane+13);
			gTracks[iPlane]->Draw("AP");
			//hEdit[iPlane]->GetZaxis()->SetRangeUser(0,200);
			hEdit[iPlane]->Draw("colz same");
		}
		std::cout << "Save Canvas" << std::endl;
		c->SaveAs(Form("%s.pdf",c->GetName()));
		//std::cout << "Saved Canvas" << std::endl;
		// for (size_t iPlane{0}; iPlane < 3; iPlane++) {
		// 	hImage[iPlane]->Delete();
		// }
	}
	std::vector<std::vector<double>> ImageCal::CreatePoints(double spacing) {
		std::vector<std::vector<double>> all_points;
		std::vector<double> x_points, y_points, z_points;
		double x_point {xmin+spacing}, y_point {ymin+sppacing}, z_point {zmin+spacing};
		while (	x_point < xmax ) {
			x_points.push_back(x_point);
			x_point+=spacing;
		}
		all_points.push_back(x_points);
		while ( y_point < ymax) {
			y_points.push_back(y_point);
			y_point+=spacing;
		}
		all_points.push_back(y_points);
		while ( z_point < zmax) {
			z_points.push_back(z_point);
			z_point+=spacing;
		}
		all_points.push_back(z_points);
		return all_points;
			
	}
	void ImageCal::TrackInformation() {
//		std::cout << "In Track Info" << std::endl;
		larlite::storage_manager mgr;
		std::cout << "storage manager" << std::endl;
		mgr.set_io_mode(larlite::storage_manager::kREAD);
		std::cout << "set_io_mode" << std::endl;
		mgr.add_in_filename(track_file);
		std::cout << "add_in_filename" << std::endl;
		mgr.open();
		if (!mgr.is_open()) {
			std::cerr << "Failed to open ROOT file. Aborting" << std::endl;
		}
		else std::cout << "Yay!" << std::endl;
		while (mgr.next_event()) {
			
			std::cout << "Next Event!" << std::endl;
			std::cout << "run_id: = " << mgr.run_id() << std::endl;
			std::cout << "subrun_id = " << mgr.subrun_id() << std::endl;
			std::cout << "event_id = " << mgr.event_id() << std::endl;
			std::vector<std::vector<double>> track_event;
			const larlite::event_track *track_v = (larlite::event_track *)mgr.get_data<larlite::event_track>("trackReco");
			std::cout << "vertex info: " << track_v->size() << std::endl;
			if (track_v->size() == 0) {std::cout << "Moving to next event!" << std::endl; std::vector<std::vector<double>> empty; track_info.push_back(empty); continue;}
			std::cout << "Filling 3D Points" << std::endl;
			//larlite::event_track* ev_trk = nullptr;
			//auto const& vtx_to_trk = mgr.find_one_ass(track_v->id(), ev_trk, track_v->id().second);
			std::cout << "track info" << std::endl;
		//	if(!ev_trk) throw larlite::DataFormatException("Could not find associated track data product!");		
			//std::cout << "ev_trk size: " << ev_trk->size() << std::endl;

			for (size_t idx{0}; idx < track_v->size(); idx++) {
				std::cout << idx << std::endl;
				std::cout << "Number of Trajectory Points: " << (*track_v)[idx].NumberTrajectoryPoints() << std::endl;
				for (size_t idx2{0}; idx2 < (*track_v)[idx].NumberTrajectoryPoints(); idx2++) {
//					std::cout << idx2 << std::endl;
					TVector3 point = (*track_v)[idx].LocationAtPoint(idx2);
					std::vector<double> pos;
					pos.push_back(point.X());
					pos.push_back(point.Y());
					pos.push_back(point.Z());
					//std::vector<double> pos = (point.X(), point.Y(), point.Z());
					track_event.push_back(pos);
				}
				//track_info.push_back(track_event);
			}
			track_info.push_back(track_event);
		}
	}
}

#endif
