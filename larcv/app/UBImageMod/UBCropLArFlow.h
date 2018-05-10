/**
 * \file UBCropLArFlow.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBCropLArFlow
 *
 * @author twongjirad
 *
 * We use crops from UBSplitDetector to crop out
 * LArFlow images. This does the work we have to do
 * to change the pixel flow variables from one crop into another.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBCROPLARFLOW_H__
#define __UBCROPLARFLOW_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBCropLArFlow ... these comments are used to generate
     doxygen documentation!
  */
  class UBCropLArFlow : public ProcessBase {

  public:

    /// Default constructor
    UBCropLArFlow(const std::string name = "UBCropLArFlow");

    /// Default destructor
    ~UBCropLArFlow() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions
    static void make_cropped_flow_images( const int src_plane,
					  const larcv::ImageMeta& srcmeta,
					  const std::vector<larcv::Image2D>& croppedadc_v,
					  const std::vector<larcv::Image2D>& srcflow,
					  const std::vector<larcv::Image2D>& srcvisi,
					  const std::vector<float>& thresholds,
					  std::vector<larcv::Image2D>& cropped_flow,
					  std::vector<larcv::Image2D>& cropped_visi,
					  const larcv::logger* log );
    
    // ----------------------------------------------------------------------
    // we save data ourselves
    larcv::IOManager* foutIO;
    
  private:

    // config parameters
    std::string _input_bbox_producer;    
    std::string _input_adc_producer;
    std::string _input_vis_producer;
    std::string _input_flo_producer;
    std::string _input_cropped_producer;
    std::string _output_adc_producer;
    std::string _output_vis_producer;
    std::string _output_flo_producer;
    std::string _output_filename;
    int _max_images;
    std::vector<float> _thresholds_v;
    bool _check_flow;
    bool _make_check_image;
    
  };

  /**
     \class larcv::UBCropLArFlowFactory
     \brief A concrete factory class for larcv::UBCropLArFlow
  */
  class UBCropLArFlowProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBCropLArFlowProcessFactory() { ProcessFactory::get().add_factory("UBCropLArFlow", this); }
    /// dtor
    ~UBCropLArFlowProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBCropLArFlow(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group

