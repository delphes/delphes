#include "fastjet/ClusterSequenceArea.hh"

FASTJET_BEGIN_NAMESPACE

LimitedWarning ClusterSequenceArea::_range_warnings;
LimitedWarning ClusterSequenceArea::_explicit_ghosts_repeats_warnings;

/// print a warning if the range is unsuitable for the current
/// calculation of the area (e.g. because ghosts do not extend
/// far enough).
void ClusterSequenceArea::_warn_if_range_unsuitable(const Selector & selector) const {
  _check_selector_good_for_median(selector);

  bool no_ghosts = (_area_def.area_type() == voronoi_area)
    || (_area_def.area_type() == passive_area
        && jet_def().jet_algorithm() == kt_algorithm);
  if (! no_ghosts) {
    double rapmin, rapmax;
    selector.get_rapidity_extent(rapmin, rapmax);
    if (rapmin < -_area_def.ghost_spec().ghost_maxrap()+0.95*jet_def().R() ||
        rapmax >  _area_def.ghost_spec().ghost_maxrap()-0.95*jet_def().R()) {
      _range_warnings.warn("rapidity range for median (rho) extends beyond +-(ghost_maxrap - 0.95*R); this is likely to cause the results to be unreliable; safest option is to increase ghost_maxrap in the area definition");
    }
  }
}


FASTJET_END_NAMESPACE
