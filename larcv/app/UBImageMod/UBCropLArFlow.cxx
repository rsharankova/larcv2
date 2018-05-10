#ifndef __UBCROPLARFLOW_CXX__
#define __UBCROPLARFLOW_CXX__

#include "UBCropLArFlow.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/IOManager.h"

//larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

namespace larcv {

  static UBCropLArFlowProcessFactory __global_UBCropLArFlowProcessFactory__;

  UBCropLArFlow::UBCropLArFlow(const std::string name)
    : ProcessBase(name)
  {}

  void UBCropLArFlow::configure(const PSet& cfg)
  {

    _verbosity_             = cfg.get<int>("Verbosity");
    _input_bbox_producer    = cfg.get<std::string>("InputBBoxProducer");
    _input_adc_producer     = cfg.get<std::string>("InputADCProducer");
    _input_cropped_producer = cfg.get<std::string>("InputCroppedADCProducer");
    _input_vis_producer     = cfg.get<std::string>("InputVisiProducer");
    _input_flo_producer     = cfg.get<std::string>("InputFlowProducer");
    _output_adc_producer    = cfg.get<std::string>("OutputCroppedADCProducer");
    _output_vis_producer    = cfg.get<std::string>("OutputCroppedVisiProducer");
    _output_flo_producer    = cfg.get<std::string>("OutputCroppedFlowProducer");
    _output_filename        = cfg.get<std::string>("OutputFilename");

    _max_images             = cfg.get<int>("MaxImages",-1);
    _thresholds_v           = cfg.get< std::vector<float> >("Thresholds",std::vector<float>(3,10.0) );

    // debug
    _check_flow             = cfg.get<bool>("CheckFlow",false);
    _make_check_image       = cfg.get<bool>("MakeCheckImage",false);

    // output file
    foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
    foutIO->set_out_file( _output_filename );
    foutIO->initialize();
    
  }

  void UBCropLArFlow::initialize()
  {}

  bool UBCropLArFlow::process(IOManager& mgr)
  {
    // we split the full detector image into 3D subpieces

    // ---------------------------------------------------------------
    // get data

    // input ADC
    auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_adc_producer));
    if (!ev_in_adc) {
      LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& img_v = ev_in_adc->image2d_array();

    // input visibility/matchability
    auto ev_in_vis  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_vis_producer));
    if (!ev_in_vis) {
      LARCV_CRITICAL() << "No Input VIS Image2D found with a name: " << _input_vis_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& vis_v = ev_in_vis->image2d_array();

    // input flo
    auto ev_in_flo  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_flo_producer));
    if (!ev_in_flo) {
      LARCV_CRITICAL() << "No Input flo Image2D found with a name: " << _input_flo_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& flo_v = ev_in_flo->image2d_array();

    // input BBox
    auto ev_in_bbox  = (larcv::EventImage2D*)(mgr.get_data("bbox2d", _input_bbox_producer));
    if (!ev_in_bbox) {
      LARCV_CRITICAL() << "No Input BBox2D found with a name: " << _input_bbox_producer << std::endl;
      throw larbys();
    }

    // cropped input ADC
    auto ev_in_cropped  = (larcv::EventImage2D*)(mgr.get_data("image2d", _input_cropped_producer));
    if (!ev_in_cropped) {
      LARCV_CRITICAL() << "No Input Cropped ADC Image2D found with a name: " << _input_cropped_producer << std::endl;
      throw larbys();
    }
    const std::vector< larcv::Image2D >& cropped_v = ev_in_cropped->image2d_array();
    

    // ----------------------------------------------------------------

    // Output ADC containers
    larcv::EventImage2D* ev_out_adc  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_adc_producer);
    larcv::EventImage2D* ev_vis_adc  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_vis_producer);
    larcv::EventImage2D* ev_flo_adc  = (larcv::EventImage2D*)foutIO->get_data("image2d",_output_flo_producer);
    
    // ----------------------------------------------------------------

    // Run, subrun, event
    int run    = ev_in_adc->run();
    int subrun = ev_in_adc->subrun();
    int event  = ev_in_adc->event();
    
    // crop corresponding flow and visibility images from cropped images
    const int src_plane = 2;
    const larcv::ImageMeta& src_meta = img_v[2].meta();
    int ncrops = cropped_v.size()/3;
    for (int icrop=0; icrop<ncrops; icrop++) {
      // this is a copy. not greate. could swap if needed ...
      std::vector<larcv::Image2D> crop_v;
      for (int i=0; i<3; i++) {
	crop_v.push_back( cropped_v.at( 3*icrop+i ) );
      }
      
      LARCV_DEBUG() << "Start crop of Flow and Visibility images of image #" << icrop << std::endl;
      std::vector<larcv::Image2D> cropped_flow;
      std::vector<larcv::Image2D> cropped_visi;
      make_cropped_flow_images( src_plane, src_meta,
				crop_v, flo_v, vis_v,
				_thresholds_v,
				cropped_flow, cropped_visi,
				&logger() );

      // check the quality of the crop
      if ( _check_flow ) {
	check_cropped_images( src_plane, crop_v, _thresholds_v, cropped_flow, cropped_visi, &logger(), 0 );
      }
      
      ev_out_adc->emplace( std::move(crop_v) );
      ev_vis_adc->emplace( std::move(cropped_visi) );
      ev_flo_adc->emplace( std::move(cropped_flow) );
      
      foutIO->set_id( run, subrun, 100*event+icrop );
      foutIO->save_entry();
    }
    
    return true;
  }
  
  void UBCropLArFlow::make_cropped_flow_images( const int src_plane,
						const larcv::ImageMeta& srcmeta,
						const std::vector<larcv::Image2D>& croppedadc_v,
						const std::vector<larcv::Image2D>& srcflow,
						const std::vector<larcv::Image2D>& srcvisi,
						const std::vector<float>& thresholds,
						std::vector<larcv::Image2D>& cropped_flow,
						std::vector<larcv::Image2D>& cropped_visi,
						const larcv::logger* log ) {

    // make the cropped flow images
    // we do not use any member variables of UBCropFlow instance in order to
    //   let this function be static, so other code can use it
    //
    // inputs
    // ------
    // src_plane: source image plane ID
    // srcmeta: meta of source image
    // croppedadc_v: vector of cropped adc image, provides meta for target images
    // srcflow: uncropped flow images
    // srcvisi: uncropped visibility images
    // thresholds: threshold below which pixel not considered matchable
    //
    // output
    // ------
    // cropped_flow: cropped lar flow images
    // cropped_visi: cropped visibility images
    //
    // misc
    // ----
    // log: pointer to logger of calling processor

    if ( cropped_flow.size()!=croppedadc_v.size() ) {
      cropped_flow.clear();
      for ( auto const& adc : croppedadc_v ) {
	larcv::Image2D flo( adc.meta() );
	flo.paint(0.0);
	cropped_flow.emplace_back( std::move(flo) );
      }
    }
    if ( cropped_visi.size()!=croppedadc_v.size() ) {
      cropped_visi.clear();
      for ( auto const& adc : croppedadc_v ) {
	larcv::Image2D visi( adc.meta() );
	visi.paint(0.0);
	cropped_visi.emplace_back( std::move(visi) );
      }
    }
    

    const int targetplanes[3][2] = { {1,2},
				     {0,2},
				     {0,1} };

    const larcv::Image2D& adcimg = croppedadc_v[src_plane];
    const larcv::ImageMeta& meta = adcimg.meta();

    const larcv::ImageMeta* target_meta[2];
    target_meta[0] = &croppedadc_v[ targetplanes[src_plane][0] ].meta();
    target_meta[1] = &croppedadc_v[ targetplanes[src_plane][1] ].meta();
  
    // we scan across the adc image + flow images.
    // we check and correct the un-cropped flow value and visibility
    //  in order to fill the cropped flow and visi images
    
    for (int r=0; r<(int)meta.rows(); r++) {
    
      float tick  = meta.pos_y(r);
      int src_row = srcmeta.row(tick);
    
      for (int c=0; c<(int)meta.cols(); c++) {
      
	float wire  = meta.pos_x(c);
	int src_col = srcmeta.col(wire);
	
	float adc   = adcimg.pixel(r,c);
      
	// two target images
	for (int i=0; i<2; i++) {
	  
	  float visi = srcvisi[i].pixel( src_row, src_col );
	  float flow = srcflow[i].pixel( src_row, src_col );

	  // is the target pixel in the cropped image?
	  float target_wire = (wire + flow);
	  if ( target_wire<target_meta[i]->min_x() || target_wire>=target_meta[i]->max_x() ) {
	    // outside the target cropped image
	    cropped_flow[i].set_pixel( r, c, -4000.0 );
	    cropped_visi[i].set_pixel( r, c, 0.0 );
	  }
	  else {

	    // inside. need to adjust the flow.
	    int target_col = target_meta[i]->col( target_wire );
	    float target_adc = croppedadc_v[ targetplanes[src_plane][i] ].pixel( r, target_col );
									
	    float target_flow = target_col - c;
	    cropped_flow[i].set_pixel( r, c, target_flow );
	    cropped_visi[i].set_pixel( r, c, visi );

	    if ( target_adc<thresholds[ targetplanes[src_plane][i] ] )
	      cropped_visi[i].set_pixel( r, c, 0.0 );

	  }
	
	}// end of loop over target images
      
      }//end of loop over cols
    }//end of loop over rows

    return;
  }

  void UBCropLArFlow::check_cropped_images( const int src_plane,
					    const std::vector<larcv::Image2D>& cropped_adc_v,
					    const std::vector<float>& thresholds,
					    const std::vector<larcv::Image2D>& cropped_flow,
					    const std::vector<larcv::Image2D>& cropped_visi,
					    const larcv::logger* log, const int verbosity ) {

    // we follow the flow to the target image.
    // correct if
    //  source adc above threshold has flow value
    //  and ( source visi=1+target image pixel above threshold OR target image pixel below and source visi is 0 )
    
    const int targetplanes[3][2] = { {1,2},
				     {0,2},
				     {0,1} };
    
    int ncorrect[2] = {0};
    int nwrong_flow2nothing[2] = {0};
    int nwrong_badvisi[2] = {0};
    int nwrong_nolabel[2] = {0};
    int nabove = 0;
    
    const larcv::Image2D& src_adc = cropped_adc_v[src_plane];
    const larcv::ImageMeta& meta = src_adc.meta();
    
    for (int r=0; r<(int)meta.rows(); r++) {
      for (int c=0; c<(int)meta.cols(); c++) {
	
	if ( src_adc.pixel(r,c)<thresholds[src_plane] )
	  continue;

	nabove++;
	
      for (int i=0; i<2; i++) {
	int flow = cropped_flow[i].pixel(r,c);
	float visi = cropped_visi[i].pixel(r,c);
	
	if ( flow<=-4000 ) {
	  if ( visi<0.5)
	    ncorrect[i]++;
	  else
	    nwrong_nolabel[i]++;
	  continue;
	}
	
	int targetc = c+flow;
	
	//std::cout << "(" << r << "," << c << ") flow=" << flow << " targetc=" << targetc << std::endl;
	
	float targetadc = cropped_adc_v[ targetplanes[src_plane][i] ].pixel( r, targetc );

	if ( visi>0.5 ) {
	  if ( targetadc>=thresholds[ targetplanes[src_plane][1] ] )
	    ncorrect[i]++;
	  else
	    nwrong_flow2nothing[i]++;
	}
	else {
	  if ( thresholds[ targetplanes[src_plane][i] ]>targetadc )
	    ncorrect[i]++;
	  else
	    nwrong_badvisi[i]++;
	}
      }
      
      }
    }

    if ( verbosity==0 ) {
      (*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	<< "[source plane " << src_plane << "-> target planes (" << targetplanes[src_plane][0] << "," << targetplanes[src_plane][1] << ")] "
	<< "(ncorrect=" << float(ncorrect[0])/float(nabove) << "," << float(ncorrect[1])/float(nabove) << ")" << std::endl;
      (*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	<< "  badvisi: ("      << float(nwrong_badvisi[0])/float(nabove) << "," << float(nwrong_badvisi[1])/float(nabove) << ")" << std::endl;
      (*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)
	<< "  flow2nothing: (" << float(nwrong_flow2nothing[0])/float(nabove) << "," << float(nwrong_flow2nothing[1])/float(nabove) << ")" << std::endl;
      (*log).send(::larcv::msg::kDEBUG,    __FUNCTION__, __LINE__, __FILE__)      
	<< "  nolabel: (" << float(nwrong_nolabel[0])/float(nabove) << "," << float(nwrong_nolabel[1])/float(nabove) << ")" << std::endl;
    }
  }
  
  void UBCropLArFlow::finalize()
  {
    foutIO->finalize();
  }
  
}
#endif
