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
    _enable_img_crop       = cfg.get<bool>("CropInModule",true);
    _box_pixel_height      = cfg.get<int>("BBoxPixelHeight");
    _box_pixel_width       = cfg.get<int>("BBoxPixelWidth");
    _debug_img             = cfg.get<bool>("DebugImage",false);
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
    // for u,v, distance between wires as we move in (Y,Z) coordinates
    float dudz = 2.0*0.3/sqrt(3.0);
    float dudy = 2.0*0.3;

    float startz = 0.25*(_box_pixel_width-1.0)*dudy;
    float endz   = 1037.0-startz;
    float starty = 117.0-0.25*(_box_pixel_height-1.0)*dudz;
    float endy   = -117.0+0.25*(_box_pixel_height-1.0)*dudz;
    
    float zspan = endz-startz;
    float yspan = endy-starty;
    float dz = 0.9*_box_pixel_width*dudz;
    float dy = 1.0*_box_pixel_height*dudy;
  
    int nz = zspan/dz;
    if ( fabs(nz*dz - zspan)>0.3 )  nz++;
    int ny = yspan/dy;
    if ( fabs(ny*dy - yspan )>0.3 ) ny++;

    float zstep = zspan/nz;
    float ystep = yspan/ny;


    // --- x/tick ----
    const larcv::ImageMeta& meta = img_v.front().meta();
    float dtick = _box_pixel_height*meta.pixel_height();

    float dtickimg = (meta.max_y()-meta.min_y() - dtick);
    int nt = dtickimg/dtick;
    if ( fabs(nt*dtick-dtickimg)>0.5 ) nt++;
    float tstep  = dtickimg/nt;
    float startt = meta.min_y() + 0.5*dtick;


    LARCV_DEBUG() << "nx,ny,nz: " << nt  << " " << ny << " " << nz << std::endl;
    LARCV_DEBUG() << "start: (" << startt << ", " << starty << ", " << startz << ")" << std::endl;
    
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
    LARCV_DEBUG() << "Num lattice points: " << lattice.size() << std::endl;

    // debug
    std::vector<larcv::Image2D> coverage_v;
    if ( _debug_img ) {
      for ( int p=0; p<3; p++) {
	larcv::Image2D cov( img_v[p].meta() );
	cov.paint(0.0);
	coverage_v.emplace_back( std::move(cov) );
      }
    }
    
    // create bounding boxes around lattice points
    for ( auto const& latpt : lattice ) {
      float center_t = latpt[0];
      Double_t pos[3];
      pos[0] = latpt[0];
      pos[1] = latpt[1];
      pos[2] = latpt[2];

      LARCV_DEBUG() << "======= CROP ===================" << std::endl;
      
      for ( int p=0; p<(int)img_v.size(); p++ ) {

	float center_w = geo->NearestWire( pos, p );
	
	int center_r = img_v[p].meta().row( center_t );
	int center_c = img_v[p].meta().col( center_w );

	int minr = center_r-_box_pixel_height/2;
	int maxr = center_r+_box_pixel_height/2;
	int minc = center_c-_box_pixel_width/2;
	int maxc = center_c+_box_pixel_width/2;

	if ( maxr-minr!=_box_pixel_height )
	  maxr = minr+_box_pixel_height;
	if ( maxc-minc!=_box_pixel_width )
	  maxc = minc+_box_pixel_width;
	
	if ( minr<0 ) {
	  minr = 0;
	  maxr = _box_pixel_height;
	}
	if ( maxr>=(int)img_v[p].meta().rows() ) {
	  maxr = img_v[p].meta().rows()-1;
	  minr = maxr-_box_pixel_height;
	}
	if ( minc<0 ) {
	  minc = 0;
	  maxc = _box_pixel_width;
	}
	if ( maxc>=(int)img_v[p].meta().cols() ) {
	  maxc = img_v[p].meta().cols()-1;
	  minc = maxc - _box_pixel_width;
	}
	
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

	if ( _enable_img_crop ) {
	  larcv::Image2D cropped = img_v[p].crop( bbox );
	  output_imgs.emplace( std::move(cropped) );
	}
	output_bbox.emplace_back( std::move(bbox) );

	if ( _debug_img ) {
	  // fill in coverage map
	  for (int r=minr; r<maxr; r++) {
	    for (int c=minc; c<maxc; c++) {
	      coverage_v[p].set_pixel( r,c, coverage_v[p].pixel(r,c)+1.0 );
	    }
	  }
	}
	
      }
    }

    LARCV_DEBUG() << "Number of cropped images: " << output_imgs.image2d_array().size() << std::endl;
    LARCV_DEBUG() << "Number of cropped images per plane: " << output_imgs.image2d_array().size()/3 << std::endl;

    if ( _debug_img ) {
      auto outev_coverage = (larcv::EventImage2D*)(mgr.get_data("image2d", "coverage"));
      int nuncovered[3] = {0};
      float meancoverage[3] = {0};
      for (int p=0; p<3; p++) {
	int maxc = 3456;
	int maxr = img_v[p].meta().rows();
	if ( p<2 )
	  maxc = 2400;
	for (int r=0; r<(int)img_v[p].meta().rows(); r++) {
	  for (int c=0; c<maxc; c++) {
	    if ( coverage_v[p].pixel(r,c)<0.5 )
	      nuncovered[p]++;
	    meancoverage[p] += coverage_v[p].pixel(r,c)/float(maxc*maxr);
	  }
	}
	LARCV_INFO() << "plane " << p << ": uncovered=" << nuncovered[p] << "  meancoverage=" << meancoverage[p] << std::endl;
	outev_coverage->emplace( std::move(coverage_v[p]) );
      }
    }
    
    return true;
  }

  void UBSplitDetector::finalize()
  {}

}
#endif
