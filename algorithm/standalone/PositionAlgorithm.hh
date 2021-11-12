#ifndef POSITIONALGORITHM_HH
#define POSITIONALGORITHM_HH

/*! @file PositionAlgorithm.hh

    Class that stores the position reconstruction algorithms for the
    global fit and energy partitioning procedures.

    @author  David A. Sweigart
    @date    06/01/2018
    @version 06/02/2018

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

class PositionAlgorithm {
  public:
    /*!
        Estimates the horizontal and vertical position for a given
        matrix of weight values with the requested algorithm(s).

        @param  list of crystal indices
        @param  list of associated weight values
        @param  list of shower constants
        @param  inner and outer algorithm
        @param  number of rows and columns
        @return horizontal and vertical position
    */
    static std::pair<double, double> EstimatePosition(
        const std::vector<int>    &p_crystals,
        const std::vector<double> &p_values,
        const std::vector<double> &p_constants,
        const std::pair<int, int> &p_algorithm,
        const std::pair<int, int> &p_dimensions = {6, 9}) {
      // verify that the algorithms requested exist
      assert(p_algorithm.first  < 3);
      assert(p_algorithm.second < 3);

      // if there's only one crystal
      if (p_crystals.size() == 1) {
        // location of the only crystal
        const double y = std::floor(p_crystals.front() / p_dimensions.second);
        const double x = p_dimensions.second - 1 -
          (p_crystals.front() % p_dimensions.second);

        // return the crystal's center position
        return {x + 0.5, y + 0.5};
      } // end if there's only one crystal

      // initialize position pair
      std::pair<double, double> position(0.5, 0.5);

      // x and y indices of crystals
      std::vector<int> xIndex(p_crystals.size());
      std::vector<int> yIndex(p_crystals.size());

      // x and y totals
      std::vector<double> xTotal(p_dimensions.second, 0.0);
      std::vector<double> yTotal(p_dimensions.first , 0.0);

      // loop over the crystals
      for (std::size_t i_crystal = 0; i_crystal < p_crystals.size(); ++i_crystal) {
        // this crystal's location
        const double yLocation = std::floor(p_crystals.at(i_crystal) / p_dimensions.second);
        const double xLocation = p_dimensions.second - 1 -
          (p_crystals.at(i_crystal) % p_dimensions.second);

        // update row and column indices
        xIndex.at(i_crystal) = xLocation;
        yIndex.at(i_crystal) = yLocation;

        // update row and column totals
        xTotal.at(xLocation) += p_values.at(i_crystal);
        yTotal.at(yLocation) += p_values.at(i_crystal);
      } // end loop over the crystals

      // largest row and column total indices
      const int xLargest = std::distance(xTotal.begin(),
        std::max_element(xTotal.begin(), xTotal.end()));
      const int yLargest = std::distance(yTotal.begin(),
        std::max_element(yTotal.begin(), yTotal.end()));

      // x index offset for adjacent algorithm
      const int xOffset =
        xLargest == 0                                     ?  1 :
        xLargest == p_dimensions.second - 1               ? -1 :
        xTotal.at(xLargest + 1) < xTotal.at(xLargest - 1) ?  1 : -1;

      // y index offset for adjacent algorithm
      const int yOffset =
        yLargest == 0                                     ?  1 :
        yLargest == p_dimensions.first - 1                ? -1 :
        yTotal.at(yLargest + 1) < yTotal.at(yLargest - 1) ?  1 : -1;

      // select between outer and inner algorithm in x-direction
      const auto xAlgorithm =
          xLargest == 0 || xLargest == p_dimensions.second - 1
        ? p_algorithm.first
        : p_algorithm.second;

      // select between outer and inner algorithm in y-direction
      const auto yAlgorithm =
          yLargest == 0 || yLargest == p_dimensions.first - 1
        ? p_algorithm.first
        : p_algorithm.second;

      // perform the prescribed algorithm in x-direction
      switch (xAlgorithm) {
        case 0: position.first += xLargest;
                break;
        case 1: position.first += Algorithm(
                  xLargest, xTotal,
                  p_constants.at(0));
                break;
        case 2: position.first += Algorithm(
                  xLargest, xOffset, xTotal.at(xLargest), xTotal.at(xLargest + xOffset),
                  p_constants.at(1));
                break;
      } // end perform the prescribed algorithm in x-direction

      // perform the prescribed algorithm in y-direction
      switch (yAlgorithm) {
        case 0: position.second += yLargest;
                break;
        case 1: position.second += Algorithm(
                  yLargest, yTotal,
                  p_constants.at(0));
                break;
        case 2: position.second += Algorithm(
                  yLargest, yOffset, yTotal.at(yLargest), yTotal.at(yLargest + yOffset),
                  p_constants.at(1));
                break;
      } // end perform the prescribed algorithm in y-direction

      // return the estimated position
      return position;
    } // end function EstimatePosition

  private:
    /*!
        Estimates the one-dimentional position using the center-of-
        gravity method with logarithmic weights.  For details, refer
        to T. Awes et al., Nucl. Instr. Meth. A 311, 130 (1992).

        @param  default location
        @param  list of values
        @param  shower tail inclusion constant
        @return reconstructed position
    */
    static double Algorithm(
      const int                 &p_default,
      const std::vector<double> &p_values,
      const double              &p_constant) {
      // initialize position and weight and values total
      double position    = 0.0;
      double weightTotal = 0.0;
      double valuesTotal = std::accumulate(p_values.begin(), p_values.end(), 0.0);

      // loop over the entries
      for (std::size_t i = 0; i < p_values.size(); ++i) {
        // calculate logarithmic weight
        const double weight = (valuesTotal > 0)
          ? std::max(0.0, p_constant + std::log(p_values.at(i) / valuesTotal))
          : 0.0;

        // include the weighted position
        position += weight * i;

        // increment the weight total
        weightTotal += weight;
      } // end loop over the entries

      // return the reconstructed position
      return (weightTotal > 0) ? position / weightTotal : p_default;
    } // end function Algorithm

    /*!
        Estimates the one-dimentional position using the adjacent-
        crystal logarithimic weighting method.  For details, refer
        to T. Awes et al., Nucl. Instr. Meth. A 311, 130 (1992).

        @param  center position index
        @param  adjacent position offset
        @param  larger weight value
        @param  smaller weight value
        @param  shower profile constant
        @return reconstructed position
    */
    static double Algorithm(
        const int    &p_location,
        const int    &p_offset,
        const double &p_large,
        const double &p_small,
        const double &p_constant) {
      // return the reconstructed position
      return (p_small <= 0)
        ? p_location
        : p_location + p_offset / 2.0 -
          p_offset * p_constant * std::log(0.5 * p_large / p_small + 0.5);
    } // end function Algorithm

};

#endif
