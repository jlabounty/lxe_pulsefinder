#ifndef PRINT_HH
#define PRINT_HH

/*! @file Print.hh

    Class that stores print functions to display status information
    in the global fit algorithm to avoid code congestion.

    @author  David A. Sweigart
    @date    05/29/2018
    @version 06/02/2018

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

// external includes
#include "TString.h"

// standard library includes
#include <iomanip>
#include <iostream>

class Print {
  public:
    /*!
      Displays a matrix of values with the dimensions provided.

      @param  number of rows
      @param  number of columns
      @param  list of values
      @param  indices of values to be colored
      @return void
    */
    static void Matrix(
        const int                 &p_rows,
        const int                 &p_columns,
        const std::vector<double> &p_values,
        const std::vector<int>    &p_colors,
        const std::vector<int>    &p_masked) {
      // loop over the rows
      for (int i = p_rows - 1; i >= 0; --i) {
        // loop over the columns
        for (int j = p_columns - 1; j >= 0; --j) {

          // whether the value should be in color
          const bool colors = (std::find_if(p_colors.begin(), p_colors.end(),
            [&](const auto &l_index) { return std::abs(l_index) == p_columns * i + j; })
            != p_colors.end());

          // whether the value should be masked off
          const bool masked = (std::find_if(p_masked.begin(), p_masked.end(),
            [&](const auto &l_index) { return std::abs(l_index) == p_columns * i + j; })
            != p_masked.end());

          // setup display
          std::cout << (colors ? "\033[0;31m" : "") << std::setw(5);

          // display the value
          if (masked) {
            std::cout << "x";
          } else {
            std::cout << std::round(p_values.at(p_columns * i + j));
          }

          // cleanup display
          std::cout << (colors ? "\033[0m" : "") << (j == 0 ? "\n" : "");

        } // end loop over the columns
      } // end loop over the rows
    } // end function Matrix

    /*!
      Displays a list of values using the number of lines requested.

      @param  list of values
      @param  sorted list of indices
      @param  number of vertical lines
      @return void
    */
    static void Plot(
        const std::vector<double> &p_values,
        const std::vector<double> &p_indices,
        const int                 &p_lines) {
      // abort if no data is provided
      if (p_values .empty() ||
          p_indices.empty()) {
        std::cout << std::endl;
        return;
      }

      // find extrema
      const double minElement = *std::min_element(p_values.begin(), p_values.end());
      const double maxElement = *std::max_element(p_values.begin(), p_values.end());

      // length to display
      const int length = p_indices.back() - p_indices.front() + 1;

      // display border with limit
      std::cout << (maxElement < 0 ? "-" : "+") << std::setw(7)
        << Form("%.1f ", std::abs(maxElement)) << std::string(length + 2, '-') << std::endl;

      // loop over the list
      std::vector<double> units(length);
      std::iota(units.begin(), units.end(), p_indices.front());
      std::pair<int, double> peak(0, 0.0);
      std::vector<bool> plotted(p_values.size(), false);
      for (int i = p_lines; i >= 0; --i) {
        // find value range
        const double yMin = minElement + (maxElement - minElement) / p_lines * i;
        const double yMax = minElement + (maxElement - minElement) / p_lines * (i + 1);

        // display the point
        std::cout << (yMin <= 0.0 && 0.0 < yMax ? "    0.0 " : std::string(8, ' '));
        std::cout << (yMin <= 0.0 && 0.0 < yMax ? "-" : "|");
        for (int j = 0; j < length; ++j) {
          // find list index
          const std::size_t k = std::distance(p_indices.begin(),
            std::find(p_indices.begin(), p_indices.end(), units.at(j)));

          // if this index's value is absent
          if (k == p_values.size()) {
            std::cout << (yMin <= 0.0 && 0.0 < yMax ? "-" : " ");
          } else {
            std::cout <<
              (yMin <= p_values.at(k) && p_values.at(k) < yMax ? "\033[0;31m+\033[0m" :
              (yMin <= 0.0            && 0.0 < yMax            ? "-" :
              (!plotted.at(k)         && yMin <= 0.0           ? "\033[0;31m|\033[0m" :
              (plotted.at(k)          && 0.0 < yMax            ? "\033[0;31m|\033[0m" :  " "))));
            if (k == 0 || p_values.at(k) > peak.second) {
              peak.first  = j;
              peak.second = p_values.at(k);
            }

            // record its been plotted
            plotted.at(k) = (yMin <= p_values.at(k) && p_values.at(k) < yMax
              ? true : plotted.at(k));
          } // end if this index's value is absent
        } // end display the point

        std::cout << (yMin <= 0.0 && 0.0 < yMax ? "-" : "|") << std::endl;
      } // end loop over the list

      // display border with limit
      std::cout << (minElement < 0 ? "-" : "+") << std::setw(7)
        << Form("%.1f ", std::abs(minElement)) << std::string(length + 2, '-') << std::endl;

      // display peak marker
      std::cout << std::setw(peak.first + 11) << "^ " << Form("%.1f ", units.at(peak.first)) << std::endl;
    } // end function Plot

    /*!
      Displays a table of values with the headers and data provided.

      @param  list of columns widths
      @param  list of headers
      @param  list of output formats
      @param  list of row data
      @return void
    */
    static void Table(
        const std::vector<int>                 &p_buffers,
        const std::vector<std::string>         &p_headers,
        const std::vector<std::string>         &p_formats,
        const std::vector<std::vector<double>> &p_values) {
      // calculate length
      const int length = std::accumulate(p_buffers.begin(), p_buffers.end(), 0) +
        3 * p_buffers.size() - 1;

      // display border
      std::cout << std::string(length, '-') << std::endl;

      // display headers
      for (std::size_t i = 0; i < p_headers.size(); ++i) {
        std::cout << (i == 0 ? "" : (p_headers.at(i) == "" ? " " : "|"))
          << std::setw(p_buffers.at(i) + 1) << p_headers.at(i) << " ";
      }

      // display border
      std::cout << std::endl << std::string(length, '-');

      // display values with format
      for (std::size_t i = 0; i < p_values.size(); ++i) {
        for (std::size_t j = 0; j < p_values.at(i).size(); ++j) {
          std::cout << (j == 0 ? "\n" : "|") << std::setw(p_buffers.at(j) + 1)
            << Form(p_formats.at(j).c_str(), p_values.at(i).at(j)) << " ";
        }
      }

      // display border
      std::cout << std::endl << std::string(length, '-') << std::endl;
    } // end function Table

    /*!
      Displays a list of values with the data provided.

      @param  column width
      @param  output format
      @param  column count
      @param  list of raw data
      @param  list of raw indices
      @return void
    */
    static void List(
        const int                 &p_buffer,
        const int                 &p_columns,
        const std::string         &p_format,
        const std::vector<double> &p_values,
        const std::vector<double> &p_indices) {
      // skip if no indices provided
      if (p_indices.size() == 0) {
        return;
      }

      // the last index processed
      double lastIndex = p_indices.front() - 1;

      // display header
      std::cout << "[";

      // display rows of values with format
      for (std::size_t i = 0; i < p_indices.size(); ++i) {
        std::cout << (i % p_columns == 0 ? "\n  " : " ")
          << std::setw(p_buffer + 1)
          << (p_indices.at(i) != lastIndex + 1 ? "" : Form(p_format.c_str(), p_values.at(i)))
          << (i < p_values.size() - 1 ? "," : "\n");

        // update the last index processed
        lastIndex = p_indices.at(i);
      }

      // display trailer
      std::cout << "]" << std::endl;
    } // end function List
};

#endif
