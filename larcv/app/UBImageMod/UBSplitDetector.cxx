#ifndef __UBSPLITDETECTOR_CXX__
#define __UBSPLITDETECTOR_CXX__

#include "UBSplitDetector.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"

//larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larcv {

  static UBSplitDetectorProcessFactory __global_UBSplitDetectorProcessFactory__;

  UBSplitDetector::UBSplitDetector(const std::string name)
    : ProcessBase(name)
  {}

  void UBSplitDetector::configure(const PSet& cfg)
  {
    _input_producer        = cfg.get<std::string>("InputProducer");    
    _output_bbox_producer  = cfg.get<std::string>("OutputBBox2DProducer");
    _output_img_producer   = cfg.get<std::string>("OutputCroppedProducer");
    _enable_img_crop       = cfg.get<bool>("CropInModule");
    // _ub_trig_tick          = cfg.get<std::vector<float> >("UBTrigTick",3200.0);
    // _ub_driftv_cmperusec   = cfg.get<float>("UBDriftVel_cmperusec", 0.111);
    // _ub_usec_per_tick      = cfg.get<float>("UBUsecPerTick", 0.5 );
    // _ub_dims_xyz           = cfg.get<std::vector<float> >("UBDimsXYZ_cm");
    _box_pixel_height      = cfg.get<int>("BBoxPixelHeight");
    _box_pixel_width       = cfg.get<int>("BBoxPixelWidth");
    // _box_depth_x_cm        = cfg.get<float>("BoxXDepthcm");    
    // _box_height_y_cm       = cfg.get<float>("BoxYHeightcm");
    // _box_width_z_cm        = cfg.get<float>("BoxZWidthcm");
  }

  void UBSplitDetector::initialize()
  {}

  bool UBSplitDetector::process(IOManager& mgr)
  {
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // get data
    auto input_image  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_producer));
    if (!input_image) {
      LARCV_CRITICAL() << "No Image2D found with a name: " << _input_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_v = input_image->image2d_array();

    auto output_bbox  = mgr.get_data<larcv::EventBBox2D>(_output_bbox_producer);

    auto output_imgs  = mgr.get_data<larcv::EventImage2D>(_output_img_producer);

    // ----------------------------------------------------------------

    // get classes UB info
    const larutil::Geometry* geo       = larutil::Geometry::GetME();
    
    // first define the lattice of 3D points
    // set lattice (y,z) pitch using width of image

    // --- (y,z) ----
    float dy = _box_pixel_width*0.3;
    float dz = dy;

    int nz = (1036.0-dz)/dz;
    if ( fabs(nz*dz - (1036.0-dz) )>0.5 ) nz++;
    int ny = (117.0*2 - dy)/dy;
    if ( fabs(ny*dy - (117.0*2-dy) )>0.5 ) ny++;

    float zstep = (1036.0-dz)/nz;
    float ystep = (117.0*2-dy)/ny;
    
    float startz = 0.5*dz;
    float starty = -117.0+0.5*dy;

    // --- x/tick ----
    const larcv::ImageMeta& meta = img_v.front().meta();
    float dtick = _box_pixel_height*meta.pixel_height();

    float dtickimg = (meta.max_y()-meta.min_y() - dtick);
    int nt = dtickimg/dtick;
    if ( fabs(nt*dtick-dtickimg)>1.0 ) nt++;
    float tstep  = dtickimg/nt;
    float startt = meta.min_y() + 0.5*dtick;


    std::cout << "nx,ny,nz: " << nt  << " " << ny << " " << nz << std::endl;
    std::cout << "start: (" << startt << ", " << starty << ", " << startz << ")" << std::endl;
    
    std::vector< std::vector<float> > lattice;
    lattice.reserve( nt*nz*ny );

    
    for (int it=0; it<=nt; it++) {
      for (int iz=0; iz<=nz; iz++) {
	for (int iy=0; iy<=ny; iy++) {

	  float pz = startz + iz*zstep;
	  float py = starty + iy*ystep;
	  float pt = startt + it*tstep;

	  std::vector<float> latpt(3);
	  latpt[0] = pt;
	  latpt[1] = py;
	  latpt[2] = pz;

	  lattice.emplace_back( std::move(latpt) );
	  
	}
      }
    }
    std::cout << "Num lattice points: " << lattice.size() << std::endl;
    
    // create bounding boxes around lattice points
    for ( auto const& latpt : lattice ) {
      float center_t = latpt[0];
      Double_t pos[3];
      pos[0] = latpt[0];
      pos[1] = latpt[1];
      pos[2] = latpt[2];
      
      for ( int p=0; p<(int)img_v.size(); p++ ) {

	float center_w = geo->NearestWire( pos, p );
	
	int center_r = img_v[p].meta().row( center_t );
	int center_c = img_v[p].meta().col( center_w );

	int minr = center_r-_box_pixel_height/2;
	int maxr = center_r+_box_pixel_height/2;
	int minc = center_c-_box_pixel_width/2;
	int maxc = center_c+_box_pixel_width/2;
	float minx = img_v[p].meta().pos_x(minc);
	float maxx = img_v[p].meta().pos_x(maxc);
	float miny = img_v[p].meta().pos_y(minr);
	float maxy = img_v[p].meta().pos_y(maxr);

	if ( miny>maxy ) {
	  float tmp = miny;
	  maxy = miny;
	  miny = tmp;
	}

	larcv::BBox2D bbox( minx, miny, maxx, maxy, img_v[p].meta().id() );

	larcv::Image2D cropped = img_v[p].crop( bbox );
	
	output_bbox.emplace_back( std::move(bbox) );
	output_imgs.emplace( std::move(cropped) );
      }
    }
    

    std::cout << "Number of cropped images: " << output_imgs.image2d_array().size() << std::endl;
    
    return true;
  }

  void UBSplitDetector::finalize()
  {}

}
#endif
