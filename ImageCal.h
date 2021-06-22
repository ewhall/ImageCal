#ifndef __IMAGECAL_H__
#define __IMAGECAL_H__

#include "TH3D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"

#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/storage_manager.h"
#include "DataFormat/ImageMeta.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "UBWireTool/UBWireTool.h"

#include "Processor/ProcessBase.h"
#include "Processor/ProcessFactory.h"

#include <string>

namespace larcv {

	class ImageCal : public ProcessBase {

	public:

		ImageCal(const std::string name="ImageCal");

		~ImageCal(){}

		void initialize();
		void finalize();
		void configure(const PSet &cfg);
		bool process(IOManager &mgr);
		void LoadCalibrationMaps(std::string file_loc);
		void DrawROI(std::vector<larcv::Image2D> Full_image_v, std::vector<larcv::Image2D> Edit_image_v);
		void SetOutDir(std::string s){out_dir = s;}
		void SetLLOutName(const std::string& foutll){_foutll = foutll;}
		void TrackInformation();
		void CreateCalibration2D();
		std::vector<std::vector<double>> CreatePoints(double spacing);


	protected: 

		int _run;
		int _subrun;
		int _event;
		int _entry;
		int EventCount;
		double xmin, xmax;
		double ymin, ymax;
		double zmin, zmax;
		double xpixel, ypixel;

		std::vector<std::vector<double>> point_vector;

		std::string _img2d_producer;
		std::string _par_pix_producer;

		std::vector<larcv::ImageMeta> Standard_meta_v;
		std::vector<larcv::Image2D> Calibration_meta_v;
		std::vector<larcv::Image2D> Track_image_v;

		std::string out_dir;

		std::vector<larcv::Image2D> Full_image_v;
		TH3D *hCalibrationMap[3];

		std::string _foutll;
		larlite::storage_manager _storage;
		std::vector<std::vector<std::vector<double>>> track_info;
		std::string track_file;
		std::vector<TGraph> track_points;
		TGraph *gTracks[3];
	};	

	class ImageCalProcessFactory : public ProcessFactoryBase {
	public:
		ImageCalProcessFactory() { ProcessFactory::get().add_factory("ImageCal",this); }
		~ImageCalProcessFactory() {}
		ProcessBase* create(const std::string instance_name) { return new ImageCal(instance_name); }
	};

}

#endif
