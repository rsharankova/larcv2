/**
 * \file UBSplitDetector.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBSplitDetector
 *
 * @author twongjirad
 *
 * We carve the detector into 3D regions and define vectors of bounding boxes to
 * then produce cropped images.
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBSPLITDETECTOR_H__
#define __UBSPLITDETECTOR_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "larcv/core/DataFormat/BBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class UBSplitDetector ... these comments are used to generate
     doxygen documentation!
  */
  class UBSplitDetector : public ProcessBase {

  public:

    /// Default constructor
    UBSplitDetector(const std::string name = "UBSplitDetector");

    /// Default destructor
    ~UBSplitDetector() {}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

    // algo functions
    // static functions are defined to allow them to be reused in a stand-alone manner
    static std::vector<larcv::BBox2D> defineBoundingBoxFromCropCoords( const std::vector<larcv::Image2D>& img_v,
								       const int box_pixel_width, const int box_pixel_height, 
								       const int t1, const int t2,
								       const int u1, const int u2,
								       const int v1, const int v2,
								       const int y1, const int y2);
    
    static void cropUsingBBox2D( const std::vector<larcv::BBox2D>& bbox_vec,
				 const std::vector<larcv::Image2D>& img_v,
				 const int y1, const int y2, bool fill_y_image,				 
				 larcv::EventImage2D& output_imgs );

    

  private:

    // config parameters
    std::string _input_producer;
    std::string _output_bbox_producer;
    std::string _output_img_producer;
    bool _enable_img_crop;
    int _box_pixel_height;
    int _box_pixel_width;
    int _covered_z_width;
    bool _complete_y_crop;
    bool _debug_img;
    int _max_images;
    // randomize crops
    bool _randomize_crops;
    int  _randomize_maxcrops;
    int  _randomize_attempts;
    float _randomize_minfracpix;
  };

  /**
     \class larcv::UBSplitDetectorFactory
     \brief A concrete factory class for larcv::UBSplitDetector
  */
  class UBSplitDetectorProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    UBSplitDetectorProcessFactory() { ProcessFactory::get().add_factory("UBSplitDetector", this); }
    /// dtor
    ~UBSplitDetectorProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new UBSplitDetector(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group

