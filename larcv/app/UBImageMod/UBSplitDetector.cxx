#ifndef __UBSPLITDETECTOR_CXX__
#define __UBSPLITDETECTOR_CXX__

#include <ctime>

#include "UBSplitDetector.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"

//larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

// ROOT TRandom3
#include "TRandom3.h"

namespace larcv {

  static UBSplitDetectorProcessFactory __global_UBSplitDetectorProcessFactory__;

  UBSplitDetector::UBSplitDetector(const std::string name)
    : ProcessBase(name)
  {}

  void UBSplitDetector::configure(const PSet& cfg)
  {
    // operating parameters
    // name of tree from which to get ADC images and their meta
    _input_producer        = cfg.get<std::string>("InputProducer");

    // name of producer to store output bounding boxes for each crop
    // all bboxes stored in unrolled vector in (u,v,y) plane order
    // i.e. (bbox0_u, bbox0_v, bbox0_y, bbox1_u, bbox1_v, bbox1_y, ...)
    _output_bbox_producer  = cfg.get<std::string>("OutputBBox2DProducer");

    // we can ask this module to do the crop for us. intended to work
    // with ADC images. other truth-label images, you are on your own
    _enable_img_crop       = cfg.get<bool>("CropInModule",true);

    // if we do the crop, this is used to name the output tree
    _output_img_producer   = cfg.get<std::string>("OutputCroppedProducer");

    // set dimensions of bounding box
    // pixel height will set the tick bounds of the image
    _box_pixel_height      = cfg.get<int>("BBoxPixelHeight",512);
    
    // this parameter sets the width of the image
    // it is the width in the Z dimension that defines visible Y-wires
    // the U and V wires saved in the image are chosen to completely cover
    // this range of Y-wires. This means, the range of filled pixels in U,V
    // will be larger than Y. Y-wires are centered and the unused pixel columns
    // are blanked out
    _covered_z_width       = cfg.get<int>("CoveredZWidth",310);
    
    // enforces a maximum picture width. we clip the edges of the
    // Y-overlapping U,V wires
    _box_pixel_width       = cfg.get<int>("BBoxPixelWidth",832);

    // by default, we leave ends of y-image blank, but we have option to fill completely
    _complete_y_crop      = cfg.get<bool>("FillCroppedYImageCompletely",false);
    
    // dump a png for debuggin
    _debug_img             = cfg.get<bool>("DebugImage",false);

    // we will split the detector completely (if not in random mode)
    // this caps the number of images. useful for debug, as this is
    // a fairly slow module
    _max_images            = cfg.get<int>("MaxImages",-1);

    // we can also choose to randomly crop within the detector
    // random (t,z) coordinate will be chosen and that location
    // will be cropped
    _randomize_crops       = cfg.get<bool>("RandomizeCrops",false);

    // max number of bounding boxes to generate
    // but we only save the number of images according to MaxImages
    _randomize_attempts    = cfg.get<int>("MaxRandomAttempts",10);

    // max fraction of pixels with above threshold values
    // else we don't keep it
    // if value is 0 or less, then we do not enforce this cut
    _randomize_minfracpix  = cfg.get<float>("MinFracPixelsInCrop",-1.0);    
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
    output_bbox->clear();
    output_imgs->clear();
    
    // ----------------------------------------------------------------

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

    int zend       = zcols-zwidth/2;
    int zstart     = zwidth/2;
    int zspan      = zend-zstart;
    int nz         = zspan/(zwidth/2);
    if ( abs( (zwidth/2)*nz - zspan )!=0 )
      nz++;
    float zstep = float(zspan)/float(nz);

    // store cropped coordinates here
    std::vector< std::vector<int> > lattice;
      
    if ( !_randomize_crops ) {
      // crop in lattice pattern through detector
      LARCV_DEBUG() << "Lattice Crop" << std::endl;
      LARCV_DEBUG() << "nt,nz: " << nt  << " " << nz << std::endl;
      LARCV_DEBUG() << "start (z,t): (" << zstart << ", " << startt << ")" << std::endl;
      
      lattice.reserve( nt*nz );
      
      for (int it=0; it<=nt; it++) {
	
	float tmid = startt + it*tstep;
	
	for (int iz=0; iz<=nz; iz++) {

	  float zwire = zstart + zstep*iz;
	  std::vector<int> crop_coords = defineImageBoundsFromPosZT( zwire, tmid, zwidth, dtick,
								     _box_pixel_width, _box_pixel_height,
								     img_v );
	  lattice.emplace_back( std::move(crop_coords) );
	  
	}
      }
      std::cout << "Num lattice points: " << lattice.size() << std::endl;
    }
    else {
      // random cropping
      TRandom3 rand(time(NULL));
      for (int iatt=0; iatt<_randomize_attempts; iatt++) {
	
	float z;
	if ( !_complete_y_crop )
	  z = zstart + zspan*rand.Uniform();
	else
	  z = _box_pixel_width + (zcols-2*_box_pixel_width)*rand.Uniform();
	
	float t = startt + dtickimg*rand.Uniform();

	std::vector<int> crop_coords = defineImageBoundsFromPosZT( z, t, zwidth, dtick,
								   _box_pixel_width, _box_pixel_height,
								   img_v );
	lattice.emplace_back( std::move(crop_coords) );
      }
      //std::cout << "Num of randomized cropping points generated: " << lattice.size() << std::endl;
    }

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
    int nfilled = 0;
    int nrejected = 0;
    for ( auto const& cropcoords : lattice ) {

      if ( _max_images>0 && _max_images<=(int)(output_imgs->image2d_array().size()/3) )
	break;
      
      int y1 = cropcoords[0];
      int y2 = cropcoords[1];
      int u1 = cropcoords[2];
      int u2 = cropcoords[3];
      int v1 = cropcoords[4];
      int v2 = cropcoords[5];
      int t1 = cropcoords[6];
      int t2 = cropcoords[7];

      std::vector<larcv::BBox2D> bbox_vec = defineBoundingBoxFromCropCoords( img_v, _box_pixel_width, _box_pixel_height,
									     t1, t2, u1, u2, v1, v2, y1, y2 );

      bool filledimg = false;
      if ( _enable_img_crop ) {
	filledimg = cropUsingBBox2D( bbox_vec, img_v, y1, y2, _complete_y_crop, _randomize_minfracpix, *output_imgs );
      }

      if ( filledimg ) {
	nfilled ++;
	for ( auto &bbox : bbox_vec ) {
	  output_bbox->emplace_back( std::move(bbox) );
	}
      }
      else {
	nrejected ++;
      }
      
    }///end of loop over lattice

    LARCV_DEBUG() << "Number of cropped images: " << output_imgs->image2d_array().size() << std::endl;
    LARCV_DEBUG() << "Number of cropped images per plane: " << output_imgs->image2d_array().size()/3 << std::endl;
    //LARCV_INFO()  << "BBoxes gen'ed=" << lattice.size() << " filled=" << nfilled << " rejected=" << nrejected << std::endl;
    std::cout  << "BBoxes gen'ed=" << lattice.size() << " filled=" << nfilled << " rejected=" << nrejected << std::endl;    

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

  std::vector<larcv::BBox2D> UBSplitDetector::defineBoundingBoxFromCropCoords( const std::vector<larcv::Image2D>& img_v,
									       const int box_pixel_width, const int box_pixel_height, 
									       const int t1, const int t2,
									       const int u1, const int u2,
									       const int v1, const int v2,
									       const int y1, const int y2) {

    // takes pre-defined image bounds on all 3 planes (given in min/max row/col)
    // note, box_pixel_width and box_pixel_height are meant to be the same
    // variable as member variables _box_pixel_width and _box_pixel_height.
    // we pass them here in order to make this function static, so it can be used as a stand-alone function.
    
    // input
    // ------
    // img_v: source ADC images from which we are cropping
    // box_pixel_width: corresponds to _box_pixel_width
    // box_pixel_height: corresponds to _box_pixel_height
    // (x)1, (x)2, where x=t,u,v,y
    // row bounds (t) and col bounds (u,v,y) max-exclusive [x1,x2)
    //
    // output
    // -------
    // (return) vector of bounding boxes defined for (u,v,y) plane

    std::vector< larcv::BBox2D > bbox_vec; // we create one for each plane
    bbox_vec.reserve(img_v.size());
    
    const larcv::ImageMeta& meta = img_v.front().meta(); // it is assumed that time-coordinate meta same between planes

    // define tick and row bounds
    int nrows = box_pixel_height;
    float mint = meta.pos_y(t1);
    float maxt = meta.pos_y(t2);

    // we crop an image with W x H = maxdu x _box_pixel_height
    // we embed in the center, the Y-plane source image with zwidth across
    // we crop the entire range for the U or V plane, the target images
    
    //LARCV_DEBUG() << "Defining bounding boxes" << std::endl;
    
    // prepare the u-plane
    const larcv::ImageMeta& umeta = img_v[0].meta();
    float minu = umeta.pos_x( u1 );
    float maxu = umeta.pos_x( u2 );
    larcv::BBox2D bbox_u( minu, mint, maxu, maxt, img_v[0].meta().id() );
    larcv::ImageMeta metacropu( minu, mint, maxu, maxt, nrows, box_pixel_width, img_v[0].meta().id() );

    // prepare the v-plane
    const larcv::ImageMeta& vmeta = img_v[1].meta();
    float minv = vmeta.pos_x( v1 );
    float maxv = vmeta.pos_x( v2 );
    larcv::BBox2D bbox_v( minv, mint, maxv, maxt, img_v[1].meta().id() );

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
    larcv::BBox2D bbox_y( miny, mint, maxy, maxt, ymeta.id() );

    bbox_vec.emplace_back( std::move(bbox_u) );
    bbox_vec.emplace_back( std::move(bbox_v) );
    bbox_vec.emplace_back( std::move(bbox_y) );    

    return bbox_vec;
    
  }

  bool UBSplitDetector::cropUsingBBox2D( const std::vector<larcv::BBox2D>& bbox_vec,
					 const std::vector<larcv::Image2D>& img_v,
					 const int y1, const int y2, bool fill_y_image,
					 const float minpixfrac,
					 larcv::EventImage2D& output_imgs ) {
    // inputs
    // ------
    // bbox_v, vector of bounding boxes for (u,v,y)
    // img_v, source adc images
    // y1, y2: range of y-wires fully covered by U,V
    // fill_y_image: if true, we fill the entire cropped y-image.
    //               else we only fill the region of y-wires that
    //               are entirely covered by the u,v cropped images
    // minpixfrac: if value is >0, we enforce a minimum value on the
    //               number of pixels occupied in the Y-image
    //
    // outputs
    // --------
    // output_imgs, cropped output image2d instances filled into eventimage2d container


    // get bounding boxes
    const larcv::BBox2D& bbox_u = bbox_vec[0];
    const larcv::BBox2D& bbox_v = bbox_vec[1];
    const larcv::BBox2D& bbox_y = bbox_vec[2];

    // y-plane meta
    const larcv::ImageMeta& ymeta = img_v[2].meta();
    
    // Y copies range of y-wires into center of output crop
    
    larcv::ImageMeta crop_yp( bbox_y.min_x(), bbox_y.min_y(), bbox_y.max_x(), bbox_y.max_y(),
			      (int)(bbox_y.height()/ymeta.pixel_height()),
			      (int)(bbox_y.width()/ymeta.pixel_width()),
			      ymeta.id() );

    std::vector<larcv::Image2D> y_img_holder;
    
    if ( fill_y_image ) {
      // we fill y-columns will all values
      //std::cout << "Y-crop: " << bbox_y.dump() << std::endl;
      larcv::Image2D crop_yimg = img_v[2].crop( bbox_y );
      y_img_holder.emplace_back( std::move(crop_yimg) );
    }
    else {
      // we only fill y-wires that are fully covered by U,V wires
      std::cout << "center-filled y plane crop" << std::endl;
      int t1 = ymeta.row( bbox_y.min_y() );
      
      larcv::Image2D ytarget( crop_yp );
      ytarget.paint(0.0);
      for (int c=0; c<(int)crop_yp.cols(); c++) {
	float cropx = crop_yp.pos_x(c);
	if ( cropx<y1 || cropx>=y2 )
	  continue;
	int cropc = ymeta.col(cropx);
	for (int r=0; r<(int)crop_yp.rows(); r++) {
	  ytarget.set_pixel( r, c, img_v[2].pixel( t1+r, cropc ) );
	}
      }
      y_img_holder.emplace_back( std::move( ytarget ) );    
    }

    float frac_occupied = 0.;
    bool saveimg = true;

    if ( minpixfrac>0 ) {
      const larcv::Image2D& yimg = y_img_holder[0];
      for (int row=0; row<(int)yimg.meta().rows(); row++) {
	for (int col=0; col<(int)yimg.meta().cols(); col++) {
	  if ( yimg.pixel(row,col)>10.0 )
	    frac_occupied+=1.0;
	}
      }
      frac_occupied /= float(yimg.meta().rows()*yimg.meta().cols());

      if ( frac_occupied<minpixfrac )
	saveimg = false;
    }

    
    if ( !saveimg )
      return false;
    

    larcv::Image2D crop_up = img_v[0].crop( bbox_u );
    output_imgs.emplace( std::move(crop_up) );
    
    larcv::Image2D crop_vp = img_v[1].crop( bbox_v );
    output_imgs.emplace( std::move(crop_vp) );

    output_imgs.emplace( std::move(y_img_holder[0]) );
    
    return true;
  }

  std::vector<int> UBSplitDetector::defineImageBoundsFromPosZT( const float zwire, const float tmid, const float zwidth, const float dtick,
								const int box_pixel_width, const int box_pixel_height,
								const std::vector<larcv::Image2D>& img_v ) {
    const larcv::ImageMeta& meta = img_v.front().meta();
    const larutil::Geometry* geo       = larutil::Geometry::GetME();
    
    float t1 = tmid-0.5*dtick;
    float t2 = tmid+0.5*dtick;
    int r1 = meta.row( t1 );
    int r2 = meta.row( t2 );

    // fix tick bounds
    if ( r2-r1!=box_pixel_height ) {
      r1 = r2-box_pixel_height;
    }
    
    if ( r1<0 ) {
      r1 = 0;
      r2 = r1 + box_pixel_height;
    }
    if ( r2>= (int)meta.rows() ) {
      r2 = (int)meta.rows()-1;
      r1 = r2-meta.rows();
    }

    // set z range
    int zcol0 = zwire - zwidth/2;
    int zcol1 = zwire + zwidth/2;
    
    if ( zcol1>3455 )
      zcol1 = 3455;

    
    // determine range for u-plane
    Double_t xyzStart[3];
    Double_t xyzEnd[3];
    geo->WireEndPoints( 2, zcol0, xyzStart, xyzEnd );
    float z0 = xyzStart[2];
    Double_t zupt0[3] = { 0,+117.5, z0 };
    int ucol0 = geo->NearestWire( zupt0, 0 );

    geo->WireEndPoints( 2, zcol1, xyzStart, xyzEnd );
    float z1 = xyzStart[2];
    Double_t zupt1[3] = { 0,-117.5, z1-0.1 };
    int ucol1 = 0;
    try {
      ucol1 = geo->NearestWire( zupt1, 0 );
    }
    catch (...) {
      ucol1 = 2399;
    }

    if ( ucol0>ucol1 ) {
      // this happens on the detector edge
      ucol0 = 0;
    }
	
    // must fit in _box_pixe_width
    int ddu = ucol1-ucol0;
    int rdu = box_pixel_width%ddu;
    int ndu = ddu/box_pixel_width;
    if ( rdu!= 0) {
      if ( ndu==0 ) {
	// short, extend the end (or lower the start if near end)
	if ( ucol1+rdu<2400 )
	  ucol1+=rdu;
	else
	  ucol0-=rdu;
      }
      else {
	rdu = ddu%box_pixel_width;
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
    try {
      vcol1 = geo->NearestWire( zvpt1, 1 );
    }
    catch  (...) {
      vcol1 = 2399;
    }
    
    int ddv = vcol1-vcol0;
    int rdv = box_pixel_width%ddv;
    int ndv = ddv/box_pixel_width;
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
	rdv = ddv%box_pixel_width;
	vcol1 -= rdv;
      }
    }
    
	
    // LARCV_DEBUG() << "Crop: z=[" << z0 << "," << z1 << "] zcol=[" << zcol0 << "," << zcol1 << "] "
    // 		   << "u=[" << ucol0 << "," << ucol1 << "] du=" << ucol1-ucol0 << " "
    // 		   << "v=[" << vcol0 << "," << vcol1 << "] dv=" << vcol1-vcol0 << " "
    // 		   << "t=[" << r1 << "," << r2 << "]"
    // 		   << std::endl;

    std::vector<int> crop_coords(8);
    crop_coords[0] = zcol0;
    crop_coords[1] = zcol1;
    crop_coords[2] = ucol0;
    crop_coords[3] = ucol1;
    crop_coords[4] = vcol0;
    crop_coords[5] = vcol1;
    crop_coords[6] = r1;
    crop_coords[7] = r2;
	
    return crop_coords;
    
  }

  
  void UBSplitDetector::finalize()
  {}

}
#endif
