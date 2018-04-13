#ifndef __UBSPLITDETECTOR_CXX__
#define __UBSPLITDETECTOR_CXX__

#include "UBSplitDetector.h"
#include "larcv/core/DataFormat/EventBBox.h"
#include "larcv/core/DataFormat/EventImage2D.h"

//larlite
#include "LArUtil/Geometry.h"

namespace larcv {

  static UBSplitDetectorProcessFactory __global_UBSplitDetectorProcessFactory__;

  UBSplitDetector::UBSplitDetector(const std::string name)
    : ProcessBase(name)
  {}

  void UBSplitDetector::configure(const PSet& cfg)
  {
    _input_producer        = cfg.get<std::string>("InputProducer");    
    _output_bbox_producer  = cfg.get<std::string>("OutputBBox2DProducer");
    // _ub_trig_tick          = cfg.get<std::vector<float> >("UBTrigTick",3200.0);
    // _ub_driftv_cmperusec   = cfg.get<float>("UBDriftVel_cmperusec", 0.111);
    // _ub_usec_per_tick      = cfg.get<float>("UBUsecPerTick", 0.5 );
    // _ub_dims_xyz           = cfg.get<std::vector<float> >("UBDimsXYZ_cm");
    _box_depth_x_cm        = cfg.get<float>("BoxXDepthcm");    
    _box_height_y_cm       = cfg.get<float>("BoxYHeightcm");
    _box_width_z_cm        = cfg.get<float>("BoxZWidthcm");
  }

  void UBSplitDetector::initialize()
  {}

  bool UBSplitDetector::process(IOManager& mgr)
  {
    /*
    auto input_image  = (EventImage2D*)(mgr.get_data("image2d", _input_producer));
    if (!input_image) {
      LARCV_CRITICAL() << "No Image2D found with a name: " << _input_producer << std::endl;
      throw larbys();
    }

    auto output_image = (EventImage2D*)(mgr.get_data("image2d", _output_producer));
    if (!output_image) {
      LARCV_CRITICAL() << "No Image2D found with a name: " << _output_producer << std::endl;
      throw larbys();
    }

    for (auto const& idx : _image_idx) {
      if (idx >= input_image->image2d_array().size()) {
        LARCV_CRITICAL() << "ImageIndex array contains index " << idx
                         << " not available in Image2DArray (" << input_image->image2d_array().size()
                         << ")!" << std::endl;
        throw larbys();
      }
    }

    auto const& bbox_v = mgr.get_data<larcv::EventBBox2D>(_bbox_producer);

    if (bbox_v.size() <= _bbox_idx) {
      LARCV_CRITICAL() << "BBOX index " << _bbox_idx << " not found!" << std::endl;
      throw larbys();
    }

    std::vector<larcv::Image2D> image_v;
    if (_input_producer == _output_producer) {
      std::vector<larcv::Image2D> tmp_v;
      input_image->move(tmp_v);
      for (auto const& idx : _image_idx)
        image_v.emplace_back(std::move(tmp_v[idx]));
    } else {
      auto const& tmp_v = input_image->image2d_array();
      for (auto const& idx : _image_idx)
        image_v.push_back(tmp_v[idx]);
    }

    // Now process
    for (size_t idx = 0; idx < _image_idx.size(); ++idx) {
      auto const& bbox = bbox_v[idx];
      image_v[idx] = image_v[idx].crop(bbox, bbox_v.unit);
    }
    output_image->emplace(std::move(image_v));

    */
    
    return true;
  }

  void UBSplitDetector::finalize()
  {}

}
#endif
