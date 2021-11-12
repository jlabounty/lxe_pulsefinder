#ifndef MODELFUNCTION_HH
#define MODELFUNCTION_HH

/*! @file ModelFunction.hh

    Class that implements the base class, ROOT::Minuit2::FCNBase,
    which calculates the function value to be minimized.

    @author  David A. Sweigart
    @date    01/25/2018
    @version 10/01/2018

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

// external includes
#include "Minuit2/FCNBase.h"
#include "TSpline.h"
#include <Eigen/Dense>

// standard library includes
#include <iostream>
#include <numeric>

class ModelFunction : public ROOT::Minuit2::FCNBase {
  public:
    ModelFunction(
        const int                                          &p_pedestalCount,
        const std::vector<std::vector<double>>             &p_traces,
        const std::vector<std::vector<double>>             &p_indices,
        const std::vector<std::shared_ptr<const TSpline3>> &p_templates,
        const std::vector<double>                          &p_noiseLevels,
        const double                                       &p_artificialDeadTime);

    /*!
        Calculates the chi-squared for all the crystals and clusters using
        the model with the provided parameters.

        @param  parameter model values
        @return chi-squared function value
    */
    double operator()(const std::vector<double> &p_parameters) const;

    /*!
        Calculates the chi-squared for all the clusters in the provided
        crystal and sample time window using the model with the provided
        parameters.

        @param  parameter model values
        @param  crystal index
        @param  list of sample indices
        @return chi-squared function value
    */
    double CrystalChiSquared(
        const std::vector<double> &p_parameters,
        const int                 &p_crystal,
        const std::vector<int>    &p_samples) const;

    /*!
        Provides the parameter error definiton.

        @return error definition value
    */
    double Up() const {
      return 1.0;
    } // end function Up

    /*!
        Appends a new cluster to the model for the provided crystals and
        in the provided sample time window.

        @param  pedestal index
        @param  list of crystals
        @param  list of sample indices
        @return void
    */
    void AddCluster(
        const int                 &p_pedestal,
        const std::vector<int>    &p_crystals,
        const std::vector<double> &p_samples) {
      // include the crystals to the new cluster
      ExtendCluster(p_pedestal, m_clusterCount++, p_crystals, p_samples);
    } // end function AddCluster

    /*!
        Extends a present cluster in the model for the provided crystals
        and in the provided sample time window.

        @param  pedestal index
        @param  cluster index
        @param  list of crystals
        @param  list of sample indices
        @return void
    */
    void ExtendCluster(
        const int                 &p_pedestal,
        const int                 &p_cluster,
        const std::vector<int>    &p_crystals,
        const std::vector<double> &p_samples);

    /*!
        Removes a cluster from the model for all the crystals.

        @param  cluster index
        @return void
    */
    void RemoveCluster(const int &p_cluster);

    /*!
        Merges a list of value goups into one with no repeated values.

        @param  lists of value groups
        @return void
    */
    std::vector<int> MergeGroups(
        const std::vector<std::vector<int>> &p_groups) const {
      // initialize flattened list
      std::vector<int> values;

      // loop over the groups
      for (const auto &group : p_groups) {
        // loop over the values
        for (const auto &value : group) {

          // check if repeated value
          if (std::find(values.begin(), values.end(), value) == values.end()) {
            // append new value
            values.push_back(value);
          }

        } // end loop over the values
      } // end loop over the groups

      // return flattened list
      return values;
    } // end function MergeGroups

    /*!
        Determines whether the provided crystal's sample index is too
        close to a present cluster and hence is masked off.

        @param  parameter model values
        @param  crystal index
        @param  sample index
        @return true if a masked index
    */
    bool MaskedIndex(
        const std::vector<double> &p_parameters,
        const int                 &p_crystal,
        const int                 &p_sample) const {
      // whether this crystal's sample index is masked
      bool masked = false;

      // loop over the clusters
      for (const auto &cluster : m_crystals.at(p_crystal).clusters) {
        // if this residual is within the artifical dead time of a cluster
        if (std::round(std::abs(m_indices.at(p_crystal).at(p_sample) -
            p_parameters.at(cluster))) <= m_artificialDeadTime) {

          // ignore this residual
          masked = true;

        } // end if this residual is within the artifical dead time of a cluster
      } // end loop over the clusters

      // return masked information
      return masked;
    } // end function MaskedIndex

    /*!
        Calculates the residual trace for each crystal using the model
        with the provided parameters.

        @param  parameter model values
        @return residual trace list
    */
    std::vector<std::vector<double>> Residuals(
        const std::vector<double> &p_parameters) const;

    /*!
        Provides the number of fit samples used in the model for the
        provided crystal index for the provided cluster.

        @param  cluster index
        @param  crystal index
        @return sample count
    */
    int FitSampleCount(
        const int &p_cluster,
        const int &p_crystal) const {
      // return this crystal's sample count for this cluster
      return m_crystals.at(p_crystal).sampleGroups.at(
        std::distance(m_crystals.at(p_crystal).clusters.begin(), std::find(
        m_crystals.at(p_crystal).clusters.begin(),
        m_crystals.at(p_crystal).clusters.end(), p_cluster))).size();
    } // end function FitSampleCount

    /*!
        Provides the number of samples in the model's signal trace for
        the provided crystal index, excluding any samples that are within
        the artifical dead time of a cluster.

        @param  parameter model values
        @param  crystal index
        @return sample count
    */
    int TraceSampleCount(
        const std::vector<double> &p_parameters,
        const int                 &p_crystal) const;

    /*!
        Calculates the number of degrees of freedom for all the crystals
        using the number of the model parameters provided.

        @param  parameter model values
        @return number of degrees of freedom
    */
    int DegreesOfFreedom(const std::vector<double> &p_parameters) const {
      // the number of samples minus the number of parameters
      return std::accumulate(m_crystals.begin(), m_crystals.end(), 0,
        [](const auto &l_value, const auto &l_crystal) {
        return l_value + l_crystal.samples.size(); }) - p_parameters.size();
    } // end function DegreesOfFreedom

    /*!
        Calculates and provides the pedestal and template scale for all
        crystals in each cluster in the model.

        @param  parameter model values
        @return pedestal lists
    */
    std::vector<std::pair<std::vector<double>, std::vector<double>>> LinearParameters(
        const std::vector<double> &p_parameters) const;

    /*!
        Provides the present number of clusters in the model.

        @return cluster count
    */
    int ClusterCount() const {
      // return number of clusters
      return m_clusterCount;
    } // end function ClusterCount

  private:
    // crystal parameter information
    typedef struct {
      // cluster, pedestal, and sample indices
      std::vector<int>              clusters;     // per cluster
      std::vector<int>              pedestals;    // per cluster
      std::vector<std::vector<int>> sampleGroups; // per cluster
      std::vector<int>              samples;      // merged between clusters
    } Crystal;

    /*!
        Analytically minimizes the linear parameters, i.e., template scales
        and pedestal, for the provided crystal.

        @param  parameter model values
        @param  crystal index
        @return result matrix
    */
    Eigen::MatrixXd AnalyticalMinimization(
        const std::vector<double> &p_parameters,
        const int                 &p_crystal) const;

    /*!
        Evaluates the model function value using the provided parameters.

        @param  crystal index
        @param  trace sample index
        @param  peak time value
        @return model function value
    */
    double Evaluate(
        const int    &p_crystal,
        const double &p_sample,
        const double &p_peakTime) const;

    /*!
        Provides the weight value for the provided crystal.

        @param  crystal index
        @return weight value
    */
    double Weight(const int &p_crystal) const;

    // sample trace lists
    const std::vector<std::vector<double>> m_traces;
    const std::vector<std::vector<double>> m_indices;

    // template models for the calorimeter
    const std::vector<std::shared_ptr<const TSpline3>> m_templates;

    // noise level list
    const std::vector<double> m_noiseLevels;

    // artificial dead time
    const double m_artificialDeadTime;

    // number of pedestals in model
    const int m_pedestalCount;

    // number of clusters in model
    int m_clusterCount;

    // crystal information list
    std::vector<Crystal> m_crystals;

    // pedestal index list
    std::vector<std::vector<int>> m_pedestals;
};

#endif
