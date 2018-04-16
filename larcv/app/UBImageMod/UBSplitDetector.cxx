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
    _box_pixel_height      = cfg.get<int>("BBoxPixelHeight",512);
    _box_pixel_width       = cfg.get<int>("BBoxPixelWidth",832);
    _covered_z_width       = cfg.get<int>("CoveredZWidth",310);
    _debug_img             = cfg.get<bool>("DebugImage",false);
    _max_images            = cfg.get<int>("MaxImages",-1);
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
    
    larcv::EventBBox2D*  output_bbox  = (larcv::EventBBox2D*)mgr.get_data( "bbox2d",_output_bbox_producer);
    larcv::EventImage2D* output_imgs  = (larcv::EventImage2D*)mgr.get_data("image2d",_output_img_producer);

    // ----------------------------------------------------------------

    // get classes UB info
    const larutil::Geometry* geo       = larutil::Geometry::GetME();
    
    // first define the lattice of 3D points
    // set lattice (y,z) pitch using width of image

    // --- image parameters ---
    // we aim to make an image where all y-charge has a partner to match against
    //const float dudz = 0.3*2;
    //const float dudy = 0.3*2/sqrt(3);
    //const float detheight = 117.0*2.0;
    int zwidth     = _covered_z_width;
    
    // --- x/tick divisions ----

    const larcv::ImageMeta& meta = img_v.front().meta();
    float dtick = _box_pixel_height*meta.pixel_height();

    float dtickimg = (meta.max_y()-meta.min_y() - dtick);
    int nt = dtickimg/dtick;
    if ( fabs(nt*dtick-dtickimg)>0.5 ) nt++;
    float tstep  = dtickimg/nt;
    float startt = meta.min_y() + 0.5*dtick;

    // --- z divisions ---------

    int zcols      = img_v.front().meta().cols();

    //int zwidth     = _box_pixel_width;
    int zend       = zcols-zwidth/2;
    int zstart     = zwidth/2;
    int zspan      = zend-zstart;
    int nz         = zspan/(zwidth/2);
    if ( abs( (zwidth/2)*nz - zspan )!=0 )
      nz++;
    float zstep = float(zspan)/float(nz);
    
    LARCV_DEBUG() << "nt,nz: " << nt  << " " << nz << std::endl;
    LARCV_DEBUG() << "start (z,t): (" << zstart << ", " << startt << ")" << std::endl;
    
    std::vector< std::vector<int> > lattice;
    lattice.reserve( nt*nz );

    Double_t xyzStart[3];
    Double_t xyzEnd[3];

    for (int it=0; it<=nt; it++) {

      float tmid = startt + it*tstep;
      float t1 = tmid-0.5*dtick;
      float t2 = tmid+0.5*dtick;
      int r1 = meta.row( t1 );
      int r2 = meta.row( t2 );

      if ( r2-r1!=_box_pixel_height ) {
	r1 = r2-_box_pixel_height;
      }
      
      if ( r1<0 ) {
	r1 = 0;
	r2 = r1 + _box_pixel_height;
      }
      if ( r2>= (int)meta.rows() ) {
	r2 = (int)meta.rows()-1;
	r1 = r2-meta.rows();
      }
      
      for (int iz=0; iz<=nz; iz++) {

	int zwire = zstart + iz*zstep;
	int zcol0 = zwire - zwidth/2;
	int zcol1 = zwire + zwidth/2;

	if ( zcol1>3455 )
	  zcol1 = 3455;

	// determine range for u-plane
	geo->WireEndPoints( 2, zcol0, xyzStart, xyzEnd );
	float z0 = xyzStart[2];
	Double_t zupt0[3] = { 0,+117.5, z0 };
	int ucol0 = geo->NearestWire( zupt0, 0 );

	geo->WireEndPoints( 2, zcol1, xyzStart, xyzEnd );
	float z1 = xyzStart[2];
	Double_t zupt1[3] = { 0,-117.5, z1-0.1 };
	int ucol1 = 0;
	if ( iz!=nz )
	  ucol1 = geo->NearestWire( zupt1, 0 );
	else
	  ucol1 = 2399;

	if ( ucol0>ucol1 ) {
	  // this happens on the detector edge
	  ucol0 = 0;
	}
	
	// must fit in _box_pixe_width
	int ddu = ucol1-ucol0;
	int rdu = _box_pixel_width%ddu;
	int ndu = ddu/_box_pixel_width;
	if ( rdu!= 0) {
	  if ( ndu==0 ) {
	    // short, extend the end (or lower the start if near end)
	    if ( ucol1+rdu<2400 )
	      ucol1+=rdu;
	    else
	      ucol0-=rdu;
	  }
	  else {
	    rdu = ddu%_box_pixel_width;
	    // long, reduce the end
	    ucol1 -= rdu;
	  }
	}

	// determine v-plane
	geo->WireEndPoints( 2, zcol0, xyzStart, xyzEnd );
	z0 = xyzStart[2];
	Double_t zvpt0[3] = { 0,-115.5, z0 };
	int vcol0 = geo->NearestWire( zvpt0, 1 );
	
	geo->WireEndPoints( 2, zcol1, xyzStart, xyzEnd );
	z1 = xyzStart[2];
	Double_t zvpt1[3] = { 0,+117.5, z1-0.1 };
	int vcol1 = 0;
	if ( iz!=nz )
	  vcol1 = geo->NearestWire( zvpt1, 1 );
	else
	  vcol1 = 2399;
	int ddv = vcol1-vcol0;
	int rdv = _box_pixel_width%ddv;
	int ndv = ddv/_box_pixel_width;
	if ( rdv!= 0) {
	  if ( ndv==0 ) {
	    // short, extend the end (or lower the start if near end)
	    if ( vcol1+rdv<2400 )
	      vcol1+=rdv;
	    else
	      vcol0-=rdv;
	  }
	  else {
	    // long, redvce the end
	    rdv = ddv%_box_pixel_width;
	    vcol1 -= rdv;
	  }
	}
	
	
	std::cout << "Crop: z=[" << z0 << "," << z1 << "] zcol=[" << zcol0 << "," << zcol1 << "] "
		  << "u=[" << ucol0 << "," << ucol1 << "] du=" << ucol1-ucol0 << " "
		  << "v=[" << vcol0 << "," << vcol1 << "] dv=" << vcol1-vcol0 << " "
		  << "t=[" << r1 << "," << r2 << "]"
		  << std::endl;

	std::vector<int> crop_coords(8);
	crop_coords[0] = zcol0;
	crop_coords[1] = zcol1;
	crop_coords[2] = ucol0;
	crop_coords[3] = ucol1;
	crop_coords[4] = vcol0;
	crop_coords[5] = vcol1;
	crop_coords[6] = r1;
	crop_coords[7] = r2;
	
	lattice.emplace_back( std::move(crop_coords) );

      }
    }
    std::cout << "Num lattice points: " << lattice.size() << std::endl;

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
    for ( auto const& cropcoords : lattice ) {

      if ( _max_images>0 && _max_images<=(int)output_imgs->image2d_array().size() )
	break;
      
      
      int y1 = cropcoords[0];
      int y2 = cropcoords[1];
      int u1 = cropcoords[2];
      int u2 = cropcoords[3];
      int v1 = cropcoords[4];
      int v2 = cropcoords[5];
      int t1 = cropcoords[6];
      int t2 = cropcoords[7];

      int nrows = _box_pixel_height;

      float mint = meta.pos_y(t1);
      float maxt = meta.pos_y(t2);

      // we crop an image with W x H = maxdu x _box_pixel_height
      // we embed in the center, the Y-plane source image with zwidth across
      // we crop the entire range for the U or V plane, the target images

      LARCV_DEBUG() << "======= CROP ===================" << std::endl;
      
      // prepare the u-plane
      const larcv::ImageMeta& umeta = img_v[0].meta();
      float minu = umeta.pos_x( u1 );
      float maxu = umeta.pos_x( u2 );
      larcv::BBox2D bbox_u( minu, mint, maxu, maxt, img_v[0].meta().id() );
      larcv::ImageMeta metacropu( minu, mint, maxu, maxt, nrows, _box_pixel_width, img_v[0].meta().id() );

      if ( _enable_img_crop ) {
	// crop if we are asked to save the image, and not just the bounding box
	larcv::Image2D crop_up = img_v[0].crop( bbox_u );
	output_imgs->emplace( std::move(crop_up) );
      }
	

      // prepare the v-plane
      const larcv::ImageMeta& vmeta = img_v[1].meta();
      float minv = vmeta.pos_x( v1 );
      float maxv = vmeta.pos_x( v2 );
      larcv::BBox2D bbox_v( minv, mint, maxv, maxt, img_v[1].meta().id() );
      if ( _enable_img_crop ) {
	larcv::Image2D crop_vp = img_v[1].crop( bbox_v );
	output_imgs->emplace( std::move(crop_vp) );
      }
      
      // prepare the y-plane
      // we take the narrow range and try to put it in the center of the y-plane image
      const larcv::ImageMeta& ymeta = img_v[2].meta();
      int ycenter = (y1+y2)/2;
      int ycmin   = ycenter - (int)metacropu.cols()/2;
      int ycmax   = ycmin + (int)metacropu.cols();
      float miny = 0;
      float maxy = 0;
      if ( ycmin>=0 && ycmax<(int)ymeta.cols() ) {
	miny = ymeta.pos_x( ycmin );
	maxy = ymeta.pos_x( ycmax );
      }
      if ( ycmin<0 ) {
	float pw = ymeta.pixel_width();
	int diffy = ycmax-ycmin;
	maxy = ymeta.pos_x( ycmax );
	miny = maxy - float(diffy)*pw;
      }
      if ( ycmax>(int)ymeta.cols() ) {
	miny = ymeta.pos_x( ycmin );
	maxy = miny + (ycmax-ycmin)*ymeta.pixel_width();
      }
      larcv::ImageMeta crop_yp( miny, mint, maxy, maxt,
				(maxt-mint)/ymeta.pixel_height(),
				ycmax-ycmin,
				ymeta.id() );
      larcv::BBox2D bbox_y( miny, miny, maxy, maxt, ymeta.id() );

      if ( _enable_img_crop ) {
	larcv::Image2D ytarget( crop_yp );
	ytarget.paint(0.0);
	for (int c=0; c<(int)crop_yp.cols(); c++) {
	  float cropx = crop_yp.pos_x(c);
	  if ( cropx<ymeta.min_x() || cropx>=ymeta.max_x() )
	    continue;
	  int cropc = ymeta.col(cropx);
	  for (int r=0; r<(int)crop_yp.rows(); r++) {
	    ytarget.set_pixel( r, c, img_v[2].pixel( t1+r, cropc ) );
	  }
	}
	std::cout << "Cropped Target-Y: " << ytarget.meta().dump() << std::endl;
      }
      
      output_bbox->emplace_back( std::move(bbox_u) );
      output_bbox->emplace_back( std::move(bbox_v) );
      output_bbox->emplace_back( std::move(bbox_y) );

      // if ( _debug_img ) {
      // 	// fill in coverage map
      // 	for (int r=minr; r<maxr; r++) {
      // 	  for (int c=minc; c<maxc; c++) {
      // 	    coverage_v[p].set_pixel( r,c, coverage_v[p].pixel(r,c)+1.0 );
      // 	  }
      // 	}
      // }
      
    }///end of loop over lattice

    LARCV_DEBUG() << "Number of cropped images: " << output_imgs->image2d_array().size() << std::endl;
    LARCV_DEBUG() << "Number of cropped images per plane: " << output_imgs->image2d_array().size()/3 << std::endl;

    // if ( _debug_img ) {
    //   auto outev_coverage = (larcv::EventImage2D*)(mgr.get_data("image2d", "coverage"));
    //   int nuncovered[3] = {0};
    //   float meancoverage[3] = {0};
    //   for (int p=0; p<3; p++) {
    // 	int maxc = 3456;
    // 	int maxr = img_v[p].meta().rows();
    // 	if ( p<2 )
    // 	  maxc = 2400;
    // 	for (int r=0; r<(int)img_v[p].meta().rows(); r++) {
    // 	  for (int c=0; c<maxc; c++) {
    // 	    if ( coverage_v[p].pixel(r,c)<0.5 )
    // 	      nuncovered[p]++;
    // 	    meancoverage[p] += coverage_v[p].pixel(r,c)/float(maxc*maxr);
    // 	  }
    // 	}
    // 	LARCV_INFO() << "plane " << p << ": uncovered=" << nuncovered[p] << "  meancoverage=" << meancoverage[p] << std::endl;
    // 	outev_coverage->emplace( std::move(coverage_v[p]) );
    //   }
    // }
    
    return true;
  }

  void UBSplitDetector::finalize()
  {}

}
#endif
