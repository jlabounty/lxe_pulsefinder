#ifndef NORMALORDER_HH
#define NORMALORDER_HH

/*! @file NormalOrder.hh

    Class that stores the expected values of normal order statistics,
    as calculated in Table 1 of H. Harter, Biometrika 48, 151 (1961).
    A Monte Carlo program with Ranlux64Engine and 1x10^12 statistics
    was used to estimate the corresponding standard deviations.

    @author  David A. Sweigart
    @date    05/25/2018
    @version 06/05/2018

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

class NormalOrder {
  public:
    /*!
      Provides the expected value of the largest observation in the sample
      size provided, scaled by the normal population's standard deviation.

      @param  sample size
      @param  normal population's standard deviation
      @return expected value of the largest observation
    */
    static double ExpectedValue(
        const int    &p_size,
        const double &p_width = 1.0) {
      // check that sample size is non-zero
      if (p_size < 1) {
        // default to sample size of one
        return p_width * values[0];
      }

      // check that sample size is within bounds
      if (p_size > 100) {
        // default to a sample size of one hundred
        return p_width * values[99];
      }

      // return pre-calculated value
      return p_width * values[p_size - 1];
    } // end function ExpectedValue

    /*!
      Provides the expected standard deviation of the expected value of the
      largest observation in the sample size provided, scaled by the normal
      population's standard deviation.

      @param  sample size
      @param  normal population's standard deviation
      @return expected width of the expected value
    */
    static double ExpectedWidth(
        const int    &p_size,
        const double &p_width = 1.0) {
      // check that sample size is non-zero
      if (p_size < 1) {
        // default to sample size of one
        return p_width * widths[0];
      }

      // check that sample size is within bounds
      if (p_size > 100) {
        // default to a sample size of one hundred
        return p_width * widths[99];
      }

      // return pre-calculated value
      return p_width * widths[p_size - 1];
    } // end function ExpectedValue

  private:
    // the expected values of the largest observation in a sample size of
    // n from a standard normal population to five decimal places
    static double values[100];

    // the expected standard deviations of the expected value of the lagrest
    // observation in a sample size of n from a standard normal population
    // to five decimal places
    static double widths[100];
};

double NormalOrder::values[100] = {
  0.00000, 0.56419, 0.84628, 1.02938, 1.16296, 1.26721, 1.35218, 1.42360, 1.48501, 1.53875,
  1.58644, 1.62923, 1.66799, 1.70338, 1.73591, 1.76599, 1.79394, 1.82003, 1.84448, 1.86748,
  1.88917, 1.90969, 1.92916, 1.94767, 1.96531, 1.98216, 1.99827, 2.01371, 2.02852, 2.04276,
  2.05646, 2.06967, 2.08241, 2.09471, 2.10661, 2.11812, 2.12928, 2.14009, 2.15059, 2.16078,
  2.17068, 2.18032, 2.18969, 2.19882, 2.20772, 2.21639, 2.22486, 2.23312, 2.24119, 2.24907,
  2.25678, 2.26432, 2.27169, 2.27891, 2.28598, 2.29291, 2.29970, 2.30635, 2.31288, 2.31928,
  2.32556, 2.33173, 2.33778, 2.34373, 2.34958, 2.35532, 2.36097, 2.36652, 2.37199, 2.37736,
  2.38265, 2.38785, 2.39298, 2.39802, 2.40299, 2.40789, 2.41271, 2.41747, 2.42215, 2.42677,
  2.43133, 2.43582, 2.44026, 2.44463, 2.44894, 2.45320, 2.45741, 2.46156, 2.46565, 2.46970,
  2.47370, 2.47764, 2.48154, 2.48540, 2.48920, 2.49297, 2.49669, 2.50036, 2.50400, 2.50759};

double NormalOrder::widths[100] = {
  1.00002, 0.82567, 0.74799, 0.70124, 0.66900, 0.64494, 0.62605, 0.61069, 0.59782, 0.58683,
  0.57731, 0.56892, 0.56146, 0.55476, 0.54870, 0.54318, 0.53811, 0.53345, 0.52912, 0.52510,
  0.52134, 0.51782, 0.51452, 0.51140, 0.50847, 0.50568, 0.50304, 0.50053, 0.49813, 0.49585,
  0.49366, 0.49156, 0.48955, 0.48763, 0.48577, 0.48399, 0.48228, 0.48063, 0.47903, 0.47749,
  0.47600, 0.47456, 0.47315, 0.47180, 0.47049, 0.46921, 0.46797, 0.46677, 0.46559, 0.46445,
  0.46333, 0.46225, 0.46120, 0.46017, 0.45917, 0.45819, 0.45723, 0.45629, 0.45538, 0.45449,
  0.45361, 0.45276, 0.45192, 0.45110, 0.45030, 0.44951, 0.44874, 0.44798, 0.44724, 0.44651,
  0.44580, 0.44509, 0.44440, 0.44372, 0.44305, 0.44240, 0.44175, 0.44112, 0.44050, 0.43989,
  0.43929, 0.43870, 0.43812, 0.43754, 0.43698, 0.43642, 0.43587, 0.43533, 0.43480, 0.43428,
  0.43376, 0.43326, 0.43276, 0.43226, 0.43177, 0.43129, 0.43082, 0.43035, 0.42989, 0.42943};

#endif
