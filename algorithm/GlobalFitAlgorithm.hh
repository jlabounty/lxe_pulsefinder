#ifndef GLOBALFITALGORITHM_HH
#define GLOBALFITALGORITHM_HH

/*! @file GlobalFitAlgorithm.hh

    Class that implements the global fit algorithm.  For details, refer
    to the dissertation of D. Sweigart, https://doi.org/10.2172/1581403.

    @author  David A. Sweigart
    @date    01/22/2018
    @version 04/10/2020

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

// external includes
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "ModelFunction.hh"
#include "NormalOrder.hh"
#include "PositionAlgorithm.hh"
#include "Print.hh"
#include "TSpline.h"

// standard library includes
#include <iostream>

class GlobalFitAlgorithm {
  public:
    // cluster parameter information
    typedef struct {
      std::size_t               index;
      double                    residualSeed;
      double                    fitSeed;
      double                    tweak;
      double                    time;
      std::pair<double, double> position;
      std::vector<double>       scales;
      std::vector<double>       pedestals;
      std::vector<double>       samples;
      std::vector<int>          crystals;
      std::vector<int>          centers;
    } Cluster;

    // cluster collection parameter information
    typedef struct {
      bool                 valid;
      double               chiSquared;
      int                  degreesOfFreedom;
      int                  functionCalls;
      double               edm;
      bool                 reachedCallLimit;
      bool                 aboveMaxEdm;
      int                  covarianceStatus;
      std::vector<Cluster> clusters;
    } ClusterCollection;

    GlobalFitAlgorithm(
        const std::vector<std::shared_ptr<const TSpline3>> &p_templates,
        const std::vector<std::shared_ptr<const TSpline3>> &p_pseudotimes,
        const std::vector<double>                          &p_noiseLevels,
        const int                                          &p_crystalRows,
        const int                                          &p_crystalColumns,
        const std::pair<int, int>                          &p_windowRows,
        const std::pair<int, int>                          &p_windowColumns,
        const std::vector<double>                          &p_showerConstants,
        const double                                       &p_clusterThreshold,
        const double                                       &p_artificialDeadTime,
        const std::vector<int>                             &p_maxFunctionCalls,
        const double                                       &p_minimizationTolerance,
        const double                                       &p_minimizationPrecision,
        const int                                          &p_strategyLevel,
        const int                                          &p_minimizerPrintLevel,
        const int                                          &p_algorithmPrintLevel,
        const bool                                         &p_refineTimeEstimate,
        const bool                                         &p_allowNegativeMetrics,
        const std::pair<int, int>                          &p_fitSamples,
        const std::pair<int, int>                          &p_positionAlgorithm,
        const std::vector<double>                          &p_seedTweaks,
        const double                                       &p_chiSquaredTolerance,
        const double                                       &p_limitTolerance);

    /*!
        Reconstructs the matrix list of signal traces provided.

        @param  lists of signal traces
        @param  lists of time indices
        @param  if time indices are non-sequential
        @return cluster collection
    */
    ClusterCollection Reconstruct(
        const std::vector<std::vector<double>> &p_traces,
        const std::vector<std::vector<double>> &p_indices,
        const bool                             &p_isSaturated) const;

    /*!
        Fits the discontiguous matrix list of signal traces provided.

        @param  lists of signal traces
        @param  lists of time indices
        @param  list of cluster indices not to fit
        @param  if a single pedestal in model
        @param  previously found cluster collections
        @return list of cluster collections
    */
    std::vector<ClusterCollection> Fit(
        const std::vector<std::vector<double>> &p_traces,
        const std::vector<std::vector<double>> &p_indices,
        const std::vector<int>                 &p_fixed,
        const bool                             &p_singlePedestal,
        const std::vector<ClusterCollection>   &p_clusters) const;

  private:
    // model parameter information
    typedef struct {
      const std::vector<std::shared_ptr<const TSpline3>> templates;
      const std::vector<double>                          noiseLevels;
    } ModelParameters;

    /*!
        Finds the crystal window with the largest metric sum, which is
        re-centered about the crystal in the largest metric.

        @param  model function
        @param  lists of signal traces
        @param  list of parameter model values
        @param  list of masked crystals
        @return crystal window
    */
    std::vector<std::pair<int, double>> FindCrystalWindow(
        const ModelFunction                    &p_model,
        const std::vector<std::vector<double>> &p_traces,
        const std::vector<double>              &p_parameters,
        const std::vector<int>                 &p_masked) const;

    /*!
        Estimates the peak time of the central crystal's signal trace.

        @param  lists of signal traces
        @param  lists of time indices
        @param  central crystal
        @return estimated peak time
    */
    double EstimateTime(
        const std::vector<std::vector<double>> &p_traces,
        const std::vector<std::vector<double>> &p_indices,
        const int                              &p_crystal) const;

    /*!
        Provides the pseudo-time correction to apply based on the peak
        neighboring samples of the central crystal's signal trace.

        @param  signal trace
        @param  central crystal
        @return pseudo-time correction
    */
    double PseudotimeCorrection(
        const std::vector<double> &p_trace,
        const int                 &p_crystal) const;

    /*!
        Determines the time indices to include in the fit based on
        the provided number of pre- and post-samples from the peak.

        @param  signal trace
        @param  time indices
        @return list of fit times
    */
    std::vector<double> TimeIndexList(
        const std::vector<double> &p_trace,
        const std::vector<double> &p_index) const;

    /*!
        Determines if two clusters are too close in time with a
        shared subset of window crystals and, if so, merges them
        together, modifying the argument objects appropriately.

        @param  pedestal index
        @param  collection offset
        @param  list of clusters
        @param  model function
        @return merging information
    */
    std::tuple<bool, bool, int> MergeClusters(
        const int                  &p_pedestal,
        const int                  &p_offset,
              std::vector<Cluster> &p_clusters,
              ModelFunction        &p_model) const;

    /*!
        Determines if any cluster has a template scale sum below the
        noise limit and, if so, provides the noise cluster's index.

        @param  collection offset
        @param  list of clusters
        @param  model function
        @return noise cluster index
    */
    int NoiseIndex(
        const int                  &p_offset,
        const std::vector<Cluster> &p_clusters,
        const ModelFunction        &p_model) const;

    // pseudo-time correction map
    const std::vector<std::shared_ptr<const TSpline3>> m_pseudotimes;

    // algorithm parameters
    const int                 m_crystalRows;
    const int                 m_crystalColumns;
    const std::pair<int, int> m_windowRows;
    const std::pair<int, int> m_windowColumns;
    const std::vector<double> m_showerConstants;
    const double              m_clusterThreshold;
    const double              m_artificialDeadTime;
    const std::vector<int>    m_maxFunctionCalls;
    const double              m_minimizationTolerance;
    const double              m_minimizationPrecision;
    const int                 m_strategyLevel;
    const int                 m_algorithmPrintLevel;
    const bool                m_refineTimeEstimate;
    const bool                m_allowNegativeMetrics;
    const std::pair<int, int> m_fitSamples;
    const std::pair<int, int> m_positionAlgorithm;
    const std::vector<double> m_seedTweaks;
    const double              m_chiSquaredTolerance;
    const double              m_limitTolerance;

    // model parameters
    ModelParameters m_modelParameters;
};

#endif
