#pragma once

#include "ThreadedGeneralUserObject.h"

/**
 * Class that provides a bin index given a spatial coordinate
 */
class SpatialBinUserObject : public ThreadedGeneralUserObject
{
public:
  static InputParameters validParams();

  SpatialBinUserObject(const InputParameters & parameters);

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  virtual void threadJoin(const UserObject &) final {}
  virtual void subdomainSetup() final {}

  virtual Real spatialValue(const Point & p) const override;

  /**
   * Whether this distribution has the distanceFromGap method defined; this
   * check is needed for some side-volume bin combinations
   * @return whether distanceFromGap is defined
   */
  virtual bool isSideBinning() const { return false; }

  /**
   * Distance between a point and a gap
   * @param[in] point point
   * @param[in] gap_index gap index
   */
  virtual const Real distanceFromGap(const Point & point, const unsigned int & gap_index) const
  {
    mooseError("'" + name() + "' does not support the 'distanceFromGap' interface!");
  }

  virtual unsigned int gapIndex(const Point & point) const
  {
    mooseError("'" + name() + "' does not support the 'gapIndex' interface!");
  }

  virtual void gapIndexAndDistance(const Point & point, unsigned int & index, Real & distance) const
  {
    mooseError("'" + name() + "' does not support the gapIndexAndDistance interface!");
  }

  /**
   * Get the bin index from a spatial point
   * @param[in] p point
   * @return bin index
   */
  virtual const unsigned int bin(const Point & p) const = 0;

  /**
   * Get the total number of bins
   * @return total number of bins
   */
  virtual const unsigned int num_bins() const = 0;

  /**
   * Get the bin given a point in an array of bounding points between layers
   * @param[in] pt point along axis of the bounding points
   * @param[in] bounds vector of bounding points
   * @return layer
   */
  unsigned int binFromBounds(const Real & pt, const std::vector<Real> & bounds) const;

  /**
   * Get the bin centers
   * @return bin centers
   */
  virtual const std::vector<Point> & getBinCenters() const { return _bin_centers; }

  /**
   * Get the coordinate directions (x, y, z) along which the bin distribution specifies
   * the bins. For 1-D distributions, this will be just one of x, y, and z. For 2-D
   * distributions, this will be a combination of x-y, y-z, or x-z.
   */
  virtual const std::vector<unsigned int> directions() const = 0;

protected:
  /// Center coordinates of the bins
  std::vector<Point> _bin_centers;
};
