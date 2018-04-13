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

  private:

    // config parameters
    std::string _input_producer;
    std::string _output_bbox_producer;
    float _box_depth_x_cm;
    float _box_height_y_cm;
    float _box_width_z_cm;
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

