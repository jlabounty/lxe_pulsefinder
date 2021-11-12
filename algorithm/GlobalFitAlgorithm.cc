/*! @file GlobalFitAlgorithm.cc

    @author  David A. Sweigart
    @date    01/22/2018
    @version 04/10/2020

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

#include "standalone/GlobalFitAlgorithm.hh"

GlobalFitAlgorithm::GlobalFitAlgorithm(
    const std::vector<std::shared_ptr<const TSpline3>> &p_templates,
    const std::vector<std::shared_ptr<const TSpline3>> &p_pseudotimes,
    const std::vector<double> &p_noiseLevels,
    const int &p_crystalRows,
    const int &p_crystalColumns,
    const std::pair<int, int> &p_windowRows,
    const std::pair<int, int> &p_windowColumns,
    const std::vector<double> &p_showerConstants,
    const double &p_clusterThreshold,
    const double &p_artificialDeadTime,
    const std::vector<int> &p_maxFunctionCalls,
    const double &p_minimizationTolerance,
    const double &p_minimizationPrecision,
    const int &p_strategyLevel,
    const int &p_minimizerPrintLevel,
    const int &p_algorithmPrintLevel,
    const bool &p_refineTimeEstimate,
    const bool &p_allowNegativeMetrics,
    const std::pair<int, int> &p_fitSamples,
    const std::pair<int, int> &p_positionAlgorithm,
    const std::vector<double> &p_seedTweaks,
    const double &p_chiSquaredTolerance,
    const double &p_limitTolerance)
    : m_pseudotimes(p_pseudotimes),
      m_crystalRows(p_crystalRows),
      m_crystalColumns(p_crystalColumns),
      m_windowRows(p_windowRows),
      m_windowColumns(p_windowColumns),
      m_showerConstants(p_showerConstants),
      m_clusterThreshold(p_clusterThreshold),
      m_artificialDeadTime(p_artificialDeadTime),
      m_maxFunctionCalls(p_maxFunctionCalls),
      m_minimizationTolerance(p_minimizationTolerance),
      m_minimizationPrecision(p_minimizationPrecision),
      m_strategyLevel(p_strategyLevel),
      m_algorithmPrintLevel(p_algorithmPrintLevel),
      m_refineTimeEstimate(p_refineTimeEstimate),
      m_allowNegativeMetrics(p_allowNegativeMetrics),
      m_fitSamples(p_fitSamples),
      m_positionAlgorithm(p_positionAlgorithm),
      m_seedTweaks(p_seedTweaks),
      m_chiSquaredTolerance(p_chiSquaredTolerance),
      m_limitTolerance(p_limitTolerance),
      m_modelParameters{
          p_templates,
          p_noiseLevels}
{
  // set minimizer's print level
  // ROOT::Minuit2::MnPrint::SetLevel(p_minimizerPrintLevel);
} // end function GlobalFitAlgorithm::GlobalFitAlgorithm

GlobalFitAlgorithm::ClusterCollection GlobalFitAlgorithm::Reconstruct(
    const std::vector<std::vector<double>> &p_traces,
    const std::vector<std::vector<double>> &p_indices,
    const bool &p_isSaturated) const
{
  // status printout
  if (m_algorithmPrintLevel > 1)
  {
    std::cout << "Reconstruction has started" << std::endl;
  }

  // setup cluster collections
  GlobalFitAlgorithm::ClusterCollection clusterCollection;
  GlobalFitAlgorithm::ClusterCollection previousClusterCollection;

  // setup model function
  ModelFunction modelFunction(
      1, p_traces, p_indices,
      m_modelParameters.templates,
      m_modelParameters.noiseLevels,
      m_artificialDeadTime);

  // inter-cluster parameters
  std::vector<int> maskedCrystals;

  // keep looking for another cluster while over threshold
  bool keepLooking = true;

  // while looking for another cluster
  while (keepLooking)
  {

    // extract the time parameters from the clusters
    std::vector<double> parameters(clusterCollection.clusters.size());
    std::transform(
        clusterCollection.clusters.begin(), clusterCollection.clusters.end(),
        parameters.begin(),
        [](const auto &l_cluster)
        { return l_cluster.residualSeed; });

    // calculate residuals for every trace
    const auto residuals = modelFunction.Residuals(parameters);

    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      std::cout << "Searching for the next crystal window" << std::endl;
    }

    // find crystal window
    auto window = FindCrystalWindow(modelFunction, residuals, parameters, maskedCrystals);

    // identify the negative crystal index
    std::size_t index = std::distance(window.begin(),
                                      std::find_if(window.begin(), window.end(),
                                                   [&](const auto &l_crystal)
                                                   { return l_crystal.first < 0; }));

    // if no negative crystal index is found
    if (index == window.size())
    {
      // the crystal index must be zero
      index = 0;
    }

    // reset the negative crystal index to be positive
    if (!window.empty())
    {
      window.at(index).first = std::abs(window.at(index).first);
    }

    // the central crystal in the window
    const int centralCrystal = (window.empty() ? -1 : window.at(index).first);

    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      std::cout
          << (window.empty() ? "No central crystal identified" : "Central crystal is ")
          << (window.empty() ? "" : std::to_string(centralCrystal)) << std::endl;
    }

    // plot residuals in terminal
    if (m_algorithmPrintLevel > 4 &&
        !window.empty())
    {
      std::cout << "Crystal " << centralCrystal << " residuals:" << std::endl;
      Print::Plot(residuals.at(centralCrystal), p_indices.at(centralCrystal), 15);
    }

    // list residuals in terminal
    if (m_algorithmPrintLevel > 5 &&
        !window.empty())
    {
      Print::List(8, 5, "%.3f", residuals.at(centralCrystal), p_indices.at(centralCrystal));
    }

    // the metric sum
    const double metricSum = std::accumulate(window.begin(), window.end(), 0.0,
                                             [](const auto &l_value, const auto &l_crystal)
                                             {
                                               return l_value + l_crystal.second;
                                             });

    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      std::cout << "Cluster's metric total is " << metricSum
                << " with a threshold of " << m_clusterThreshold << std::endl;
    }

    // if the metric total is above threshold
    if (metricSum < m_clusterThreshold)
    {
      // status printout
      if (m_algorithmPrintLevel > 1)
      {
        std::cout << "Next cluster is under threshold" << std::endl;
      }

      // no new cluster identified
      keepLooking = false;
    }
    else
    {
      // status printout
      if (m_algorithmPrintLevel > 1)
      {
        std::cout << "Found new cluster over threshold" << std::endl;
      }

      // extract the crystals in the window
      std::vector<int> crystals(window.size());
      std::transform(window.begin(), window.end(), crystals.begin(),
                     [](const auto &l_crystal)
                     { return l_crystal.first; });

      // sort the crystals in the window
      std::sort(crystals.begin(), crystals.end());

      // calculate sample indices to include in fit
      const auto samples = TimeIndexList(
          residuals.at(centralCrystal),
          p_indices.at(centralCrystal));

      // status printout
      if (m_algorithmPrintLevel > 1)
      {
        std::cout << "Adding cluster # "
                  << modelFunction.ClusterCount() << " to model" << std::endl;
      }

      // include this cluster in the model
      modelFunction.AddCluster(0, crystals, samples);

      // status printout
      if (m_algorithmPrintLevel > 1)
      {
        std::cout << Form("Estimated peak time is %.3f c.t.",
                          EstimateTime(residuals, p_indices, centralCrystal))
                  << std::endl;
      }

      // create new cluster object with estimated parameters
      GlobalFitAlgorithm::Cluster cluster;
      cluster.fitSeed = EstimateTime(residuals, p_indices, centralCrystal);
      cluster.position = {-1.0, -1.0};
      cluster.samples = samples;
      cluster.crystals = crystals;
      cluster.centers = {centralCrystal};

      // store this cluster in the collection
      clusterCollection.clusters.push_back(std::move(cluster));

      // two clusters have been merged
      bool haveMerged;

      do
      {

        // this cluster is overlapping with another cluster
        bool overlappingWindow = false;

        // loop over the previous clusters
        for (std::size_t i = 0; i < clusterCollection.clusters.size() - 1; ++i)
        {
          // loop over this cluster's crystals
          for (const auto &crystal : clusterCollection.clusters.at(i).crystals)
          {

            // if this crystal is shared
            if (std::find(
                    clusterCollection.clusters.back().crystals.begin(),
                    clusterCollection.clusters.back().crystals.end(), crystal) !=
                clusterCollection.clusters.back().crystals.end())
            {
              // this cluster's window overlaps
              overlappingWindow = true;
            } // end if this crystal is shared

          } // end loop over this cluster's crystals
        }   // end loop over the previous clusters

        // setup migrad
        ROOT::Minuit2::MnMigrad migrad(
            modelFunction,
            std::vector<double>(),
            std::vector<double>(),
            m_strategyLevel);

        // if minimizer precision is provided
        if (m_minimizationPrecision > 0)
        {
          // status printout
          if (m_algorithmPrintLevel > 1)
          {
            std::cout << "Minimization precision lowered to "
                      << m_minimizationPrecision << std::endl;
          }

          // lower precision
          migrad.SetPrecision(m_minimizationPrecision);
        } // end if minimizer precision is provided

        // loop over the clusters
        for (const auto &cluster : clusterCollection.clusters)
        {
          // status printout
          if (m_algorithmPrintLevel > 2)
          {
            std::cout << Form("Adding seed %.3f c.t. with limits [%.3f, %.3f]",
                              cluster.fitSeed, cluster.samples.front(), cluster.samples.back())
                      << std::endl;
          }

          // include this cluster's seed to migrad
          migrad.Add(Form("p%u", migrad.VariableParameters()), cluster.fitSeed, 0.1,
                     cluster.samples.front(), cluster.samples.back());
        } // end loop over the clusters

        // function call limit
        int maxFunctionCalls = 0;

        // loop over the orders
        for (std::size_t i = 0; i < m_maxFunctionCalls.size(); ++i)
        {
          maxFunctionCalls += m_maxFunctionCalls.at(i) *
                              std::pow(migrad.VariableParameters(), i);
        }

        // time seed tweaks
        std::vector<double> seedTweaks = {0.0};

        // include the provided seed tweaks
        seedTweaks.insert(seedTweaks.end(), m_seedTweaks.begin(), m_seedTweaks.end());

        // loop over the seed tweaks
        for (std::size_t i_tweak = 0; i_tweak < seedTweaks.size(); ++i_tweak)
        {
          // tweak the latest cluster seed
          migrad.SetValue(
              clusterCollection.clusters.size() - 1,
              clusterCollection.clusters.back().fitSeed + seedTweaks.at(i_tweak));

          // status printout
          if (m_algorithmPrintLevel > 1)
          {
            std::cout << "MIGRAD called with tweak "
                      << std::showpos << seedTweaks.at(i_tweak) << std::noshowpos << std::endl;
          }

          // simultaneous fit
          auto localMinimum = migrad(maxFunctionCalls, m_minimizationTolerance);

          // status printout
          if (m_algorithmPrintLevel > 6)
          {
            std::cout << "Function Minimum:\n"
                      << localMinimum << std::endl;
          }

          // extract the parameters and calculate the chi-squared
          const auto localParameters = localMinimum.UserState().Params();
          const auto localChiSquared = modelFunction(localParameters);

          // whether this fit has a lower chi-squared
          const bool lowerChiSquared =
              (clusterCollection.chiSquared - localChiSquared > m_chiSquaredTolerance);

          // if result has a lower chi-squared
          if ((i_tweak == 0) ||
              (lowerChiSquared && clusterCollection.valid && localMinimum.IsValid()) ||
              (lowerChiSquared && !clusterCollection.valid))
          {
            // retrieve linear parameters
            const auto linearParameters = modelFunction.LinearParameters(localParameters);

            // update cluster information
            for (std::size_t i = 0; i < clusterCollection.clusters.size(); ++i)
            {
              clusterCollection.clusters.at(i).residualSeed = localParameters.at(i);
              clusterCollection.clusters.at(i).time = localParameters.at(i);
              clusterCollection.clusters.at(i).scales = linearParameters.at(i).second;
              clusterCollection.clusters.at(i).pedestals = linearParameters.at(i).first;
            } // end update cluster information

            // update the seed tweak
            clusterCollection.clusters.back().tweak = seedTweaks.at(i_tweak);

            // update collection information
            clusterCollection.chiSquared = localChiSquared;
            clusterCollection.degreesOfFreedom = modelFunction.DegreesOfFreedom(localParameters);
            clusterCollection.valid = localMinimum.IsValid();
            clusterCollection.functionCalls = localMinimum.NFcn();
            clusterCollection.reachedCallLimit = localMinimum.HasReachedCallLimit();
            clusterCollection.edm = localMinimum.Edm();
            clusterCollection.aboveMaxEdm = localMinimum.IsAboveMaxEdm();
            clusterCollection.covarianceStatus = localMinimum.UserState().CovarianceStatus();
          } // end if result has a lower chi-squared

          // if converged for an isolated cluster fit
          if (localMinimum.IsValid() && !overlappingWindow)
          {
            // no need for further seed tweaks
            break;
          }
        } // end loop over the seed tweaks

        // status printout
        if (m_algorithmPrintLevel > 1)
        {
          std::cout << Form("Fit chi-squared/ndf is %.3f/%i", clusterCollection.chiSquared,
                            clusterCollection.degreesOfFreedom)
                    << " at tweak " << std::showpos
                    << clusterCollection.clusters.back().tweak << std::noshowpos << std::endl;

          std::cout << "Fit " << (clusterCollection.valid ? "converged" : "didn't converge") << std::endl;
        } // end status printout

        // window overlap between clusters to be merged
        bool shareWindow;

        // merge clusters, if necessary
        std::tie(haveMerged, shareWindow, std::ignore) = MergeClusters(
            0, 0, clusterCollection.clusters, modelFunction);

        // if two clusters have been merged
        if (haveMerged)
        {
          // skip ahead to fit again
          continue;
        }

        // the cluster index for any noise fit
        int noiseIndex = -1;

        // find the cluster index for any noise fit
        if (!shareWindow)
        {
          noiseIndex = NoiseIndex(0, clusterCollection.clusters, modelFunction);
        }

        // cluster is at the parameter limit
        bool parameterAtLimit = false;

        // find if the cluster at its parameter limit
        if (!shareWindow &&
            noiseIndex < 0 &&
            clusterCollection.clusters.size() > 0)
        {
          // this cluster
          const auto &cluster = clusterCollection.clusters.back();

          // the parameter limit tolerance
          const double limitTolerance = (m_limitTolerance > 0)
                                            ? m_limitTolerance
                                            : migrad.Precision().Eps2();

          // time is at limit within the provided precision
          parameterAtLimit =
              std::abs(cluster.time - cluster.samples.front()) < limitTolerance ||
              std::abs(cluster.time - cluster.samples.back()) < limitTolerance;
        } // end find if the cluster at its parameter limit

        // if to reject this cluster
        if (shareWindow ||
            noiseIndex >= 0 ||
            parameterAtLimit)
        {
          // mask this window's central crystal
          maskedCrystals.push_back(centralCrystal);

          // status printout
          if (m_algorithmPrintLevel > 1)
          {
            std::cout << (shareWindow ? "Clusters found close in time in the identical window" : (noiseIndex >= 0 ? "Cluster found with a template scale sum below threshold" : "Cluster found with a time at the parameter limit")) << std::endl;
            std::cout << "Masked off crystal # " << centralCrystal << std::endl;
          } // end status printout

          // check if no clusters were merged
          if (static_cast<int>(previousClusterCollection.clusters.size()) <
              modelFunction.ClusterCount())
          {
            // remove the latest cluster in the model
            modelFunction.RemoveCluster(clusterCollection.clusters.size() - 1);
          }
          else
          {
            // loop over the current clusters
            // note: they must be removed in decending order
            for (std::size_t i = 0; i < clusterCollection.clusters.size(); ++i)
            {
              // remove the current cluster in the model
              modelFunction.RemoveCluster(clusterCollection.clusters.size() - 1 - i);
            }

            // loop over the previous clusters
            for (std::size_t i = 0; i < previousClusterCollection.clusters.size(); ++i)
            {
              // include the previous cluster in the model
              modelFunction.AddCluster(0,
                                       previousClusterCollection.clusters.at(i).crystals,
                                       previousClusterCollection.clusters.at(i).samples);
            } // loop over the previous clusters
          }   // end check if no clusters were merged

          // revert back to the previous cluster collection
          clusterCollection = previousClusterCollection;
        } // end if to reject this cluster

        // status printout
        if (m_algorithmPrintLevel > 2 &&
            clusterCollection.clusters.size() > 0)
        {
          std::cout << (clusterCollection.valid ? "Valid" : "Invalid")
                    << " current state with chi-squared/ndf " << Form("%.3f/%i", clusterCollection.chiSquared, clusterCollection.degreesOfFreedom) << std::endl;

          std::vector<std::vector<double>> data(clusterCollection.clusters.size());
          for (std::size_t i = 0; i < data.size(); ++i)
          {
            data.at(i).push_back(i);
            data.at(i).push_back(
                clusterCollection.clusters.at(i).fitSeed +
                clusterCollection.clusters.at(i).tweak);
            data.at(i).push_back(clusterCollection.clusters.at(i).time);
            data.at(i).push_back(std::accumulate(
                clusterCollection.clusters.at(i).scales.begin(),
                clusterCollection.clusters.at(i).scales.end(), 0.0));
          }

          Print::Table(
              {2, 10, 10, 8},
              {"#", "Seed", "Time", "Energy"},
              {"%.0f", "%.3f", "%.3f", "%.3f"}, data);
        } // end status printout

        // if the time wandered away from the seed
        if (!p_isSaturated &&
            !overlappingWindow &&
            clusterCollection.clusters.size() > 0 &&
            std::abs(
                clusterCollection.clusters.back().time -
                clusterCollection.clusters.back().fitSeed -
                clusterCollection.clusters.back().tweak) > m_artificialDeadTime)
        {
          // status printout
          if (m_algorithmPrintLevel > 1)
          {
            std::cout << "Cluster # " << clusterCollection.clusters.size()
                      << "'s time has been reset back to its seed "
                      << clusterCollection.clusters.back().fitSeed
                      << " due to wandering away" << std::endl;
          } // end status printout

          // reset time back to the seed
          clusterCollection.clusters.back().residualSeed =
              clusterCollection.clusters.back().fitSeed;
        } // end if the time wandered away from the seed

      } while (haveMerged);

      // record current cluster collection
      previousClusterCollection = clusterCollection;
    } // end if the metric total is above threshold

  } // end while looking for another cluster

  // loop over the clusters
  for (auto &cluster : clusterCollection.clusters)
  {
    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      const auto position = PositionAlgorithm::EstimatePosition(
          cluster.crystals, cluster.scales,
          m_showerConstants, m_positionAlgorithm,
          {m_crystalRows, m_crystalColumns});
      std::cout << Form("Estimated position at [%.3f, %.3f] c.w.",
                        position.first, position.second)
                << std::endl;
    } // end status printout

    // estimate the position
    cluster.position = PositionAlgorithm::EstimatePosition(
        cluster.crystals, cluster.scales,
        m_showerConstants, m_positionAlgorithm,
        {m_crystalRows, m_crystalColumns});
  } // end loop over the clusters

  // status printout
  if (m_algorithmPrintLevel > 1)
  {
    std::cout << "Reconstruction has completed" << std::endl;
  }

  // status printout
  if (m_algorithmPrintLevel > 2 &&
      clusterCollection.clusters.size() > 0)
  {
    std::vector<std::vector<double>> data(clusterCollection.clusters.size());
    for (std::size_t i = 0; i < data.size(); ++i)
    {
      data.at(i).push_back(i);
      data.at(i).push_back(
          clusterCollection.clusters.at(i).fitSeed +
          clusterCollection.clusters.at(i).tweak);
      data.at(i).push_back(clusterCollection.clusters.at(i).time);
      data.at(i).push_back(std::accumulate(
          clusterCollection.clusters.at(i).scales.begin(),
          clusterCollection.clusters.at(i).scales.end(), 0.0));
    }

    Print::Table(
        {2, 10, 10, 8},
        {"#", "Seed", "Time", "Energy"},
        {"%.0f", "%.3f", "%.3f", "%.3f"}, data);
  } // end status printout

  // return the cluster collection
  return clusterCollection;
} // end function GlobalFitAlgorithm::Reconstruct

std::vector<GlobalFitAlgorithm::ClusterCollection> GlobalFitAlgorithm::Fit(
    const std::vector<std::vector<double>> &p_traces,
    const std::vector<std::vector<double>> &p_indices,
    const std::vector<int> &p_fixed,
    const bool &p_singlePedestal,
    const std::vector<GlobalFitAlgorithm::ClusterCollection> &p_clusters) const
{
  // status printout
  if (m_algorithmPrintLevel > 1)
  {
    std::cout << "Fit has started" << std::endl;
  }

  // setup cluster collection
  std::vector<GlobalFitAlgorithm::ClusterCollection> clusterCollections(p_clusters.size());

  // input cluster collection indices
  std::vector<double> clusterIndex;

  // store the clusters in the collections
  for (std::size_t i = 0; i < p_clusters.size(); ++i)
  {
    // include the clusters
    clusterCollections.at(i).clusters = p_clusters.at(i).clusters;

    // change the seed to be the time
    for (std::size_t j = 0; j < clusterCollections.at(i).clusters.size(); ++j)
    {
      // this cluster
      auto &cluster = clusterCollections.at(i).clusters.at(j);

      // change the seed to be the time
      cluster.fitSeed = cluster.time;

      // include the cluster index
      cluster.index = clusterIndex.size();

      // include the cluster index in collection
      clusterIndex.push_back(j);
    }
  } // end store the clusters in the collections

  // setup model function
  ModelFunction modelFunction(
      (p_singlePedestal ? 1 : p_clusters.size()), p_traces, p_indices,
      m_modelParameters.templates,
      m_modelParameters.noiseLevels,
      m_artificialDeadTime);

  // status printout
  if (m_algorithmPrintLevel > 1)
  {
    std::cout << "Adding clusters to model" << std::endl;
  }

  // loop over the collections
  for (std::size_t i = 0; i < clusterCollections.size(); ++i)
  {
    // loop over the clusters
    for (const auto &cluster : clusterCollections.at(i).clusters)
    {

      // include this cluster in the model
      modelFunction.AddCluster((p_singlePedestal ? 0 : i),
                               cluster.crystals, cluster.samples);

    } // end loop over the clusters
  }   // end loop over the collections

  // keep a copy of the clusters to not refine
  std::vector<int> fixedClusters = p_fixed;

  // two clusters have been merged or a cluster was rejected
  bool keepFitting;

  do
  {

    // setup migrad
    ROOT::Minuit2::MnMigrad migrad(
        modelFunction,
        std::vector<double>(),
        std::vector<double>(),
        m_strategyLevel);

    // if minimizer precision is provided
    if (m_minimizationPrecision > 0)
    {
      // status printout
      if (m_algorithmPrintLevel > 1)
      {
        std::cout << "Minimization precision lowered to "
                  << m_minimizationPrecision << std::endl;
      }

      // lower precision
      migrad.SetPrecision(m_minimizationPrecision);
    } // end if minimizer precision is provided

    // loop over the collections
    for (const auto &collection : clusterCollections)
    {
      // loop over the clusters
      for (const auto &cluster : collection.clusters)
      {

        // status printout
        if (m_algorithmPrintLevel > 2)
        {
          std::cout << Form("Adding seed %.3f c.t. with limits [%.3f, %.3f]",
                            cluster.fitSeed, cluster.samples.front(), cluster.samples.back())
                    << std::endl;
        }

        // include this cluster's seed to migrad
        migrad.Add(Form("p%u", migrad.VariableParameters()), cluster.fitSeed, 0.1,
                   cluster.samples.front(), cluster.samples.back());

      } // end loop over the clusters
    }   // end loop over the collections

    // loop over the clusters to not refine
    for (const auto &index : fixedClusters)
    {
      // fix this cluster's time in migrad
      migrad.Fix(index);
    }

    // function call limit
    int maxFunctionCalls = 0;

    // loop over the orders
    for (std::size_t i = 0; i < m_maxFunctionCalls.size(); ++i)
    {
      maxFunctionCalls += m_maxFunctionCalls.at(i) *
                          std::pow(migrad.VariableParameters(), i);
    }

    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      std::cout << "MIGRAD called" << std::endl;
    }

    // simultaneous fit
    auto functionMinimum = migrad(maxFunctionCalls, m_minimizationTolerance);

    // status printout
    if (m_algorithmPrintLevel > 6)
    {
      std::cout << "Function Minimum:\n"
                << functionMinimum << std::endl;
    }

    // extract the parameters and calculate the chi-squared
    const auto parameters = functionMinimum.UserState().Params();
    const auto chiSquared = modelFunction(parameters);
    const auto degreesOfFreedom = modelFunction.DegreesOfFreedom(parameters);

    // retrieve linear parameters
    const auto linearParameters = modelFunction.LinearParameters(parameters);

    // clusters' collection index
    std::vector<int> catalog;

    // identify cluster with collection
    for (std::size_t i = 0; i < clusterCollections.size(); ++i)
    {
      catalog.insert(catalog.end(), clusterCollections.at(i).clusters.size(), i);
    }

    // update cluster information
    for (std::size_t i = 0; i < parameters.size(); ++i)
    {
      // this cluster in the collection
      auto &cluster = clusterCollections.at(catalog.at(i)).clusters.at(std::distance(std::find(catalog.begin(), catalog.end(), catalog.at(i)), catalog.begin() + i));

      // update this cluster
      cluster.time = parameters.at(i);
      cluster.scales = linearParameters.at(i).second;
      cluster.pedestals = linearParameters.at(i).first;
      cluster.position = PositionAlgorithm::EstimatePosition(cluster.crystals, cluster.scales,
                                                             m_showerConstants, m_positionAlgorithm, {m_crystalRows, m_crystalColumns});
    } // end update cluster information

    // update collection information
    for (auto &collection : clusterCollections)
    {
      collection.chiSquared = chiSquared;
      collection.degreesOfFreedom = degreesOfFreedom;
      collection.valid = functionMinimum.IsValid();
      collection.functionCalls = functionMinimum.NFcn();
      collection.reachedCallLimit = functionMinimum.HasReachedCallLimit();
      collection.edm = functionMinimum.Edm();
      collection.aboveMaxEdm = functionMinimum.IsAboveMaxEdm();
      collection.covarianceStatus = functionMinimum.UserState().CovarianceStatus();
    } // end update collection information

    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      std::cout << Form("Fit chi-squared/ndf is %.3f/%i", chiSquared,
                        degreesOfFreedom)
                << "\nFit " << (functionMinimum.IsValid() ? "converged" : "didn't converge") << std::endl;
    }

    // reset inter-collection information
    keepFitting = false;

    // loop over the collections
    for (std::size_t i = 0; i < clusterCollections.size(); ++i)
    {
      // this cluster collection
      auto &clusterCollection = clusterCollections.at(i);

      // two clusters have been merged
      bool haveMerged;

      // window overlap between clusters to be merged
      bool shareWindow;

      // index of merged cluster
      int mergedIndex;

      // merge clusters, if necessary
      std::tie(haveMerged, shareWindow, mergedIndex) = MergeClusters(
          i, std::distance(catalog.begin(), std::find(catalog.begin(), catalog.end(), i)),
          clusterCollection.clusters, modelFunction);

      // if two clusters have been merged
      if (haveMerged)
      {
        // index of merged cluster offset by the collection start
        const int index = std::distance(catalog.begin(), std::find(
                                                             catalog.begin(), catalog.end(), i)) +
                          mergedIndex;

        // remove the collection index in the collection
        catalog.erase(catalog.begin() + index);

        // remove the cluster index in the collection
        clusterIndex.erase(clusterIndex.begin() + index);

        // find this cluster index in the fixed clusters list
        const auto clusterIterator = std::find(fixedClusters.begin(), fixedClusters.end(), index);

        // if this cluster index was found
        if (clusterIterator != fixedClusters.end())
        {
          // remove the cluster index in the fixed clusters list
          fixedClusters.erase(fixedClusters.begin() +
                              std::distance(fixedClusters.begin(), clusterIterator));
        }

        // lower any higher-indexed clusters by one
        std::for_each(fixedClusters.begin(), fixedClusters.end(),
                      [&](auto &l_cluster)
                      { l_cluster -= (l_cluster > index ? 1 : 0); });

        // update inter-collection information
        keepFitting = true;

        // skip ahead to fit again
        continue;
      } // end if two clusters have been merged

      // the cluster index for any noise fit
      int noiseIndex = -1;

      // find the cluster index for any noise fit
      if (!shareWindow)
      {
        noiseIndex = NoiseIndex(std::distance(catalog.begin(), std::find(
                                                                   catalog.begin(), catalog.end(), i)),
                                clusterCollection.clusters, modelFunction);
      }

      // if to reject this cluster
      if (shareWindow ||
          noiseIndex >= 0)
      {
        // status printout
        if (m_algorithmPrintLevel > 1)
        {
          std::cout << (shareWindow
                            ? "Clusters found close in time in the identical window"
                            : "Cluster found with a template scale sum below threshold")
                    << std::endl;
        } // end status printout

        // index of the offending cluster
        int index = (shareWindow) ? mergedIndex : noiseIndex;

        // remove the offending cluster in the model
        modelFunction.RemoveCluster(std::distance(catalog.begin(), std::find(
                                                                       catalog.begin(), catalog.end(), i)) +
                                    index);

        // remove the offending cluster in the collection
        clusterCollection.clusters.erase(clusterCollection.clusters.begin() + index);

        // offset index of offending cluster by the collection start
        index += std::distance(catalog.begin(),
                               std::find(catalog.begin(), catalog.end(), i));

        // remove the collection index in the collection
        catalog.erase(catalog.begin() + index);

        // remove the cluster index in the collection
        clusterIndex.erase(clusterIndex.begin() + index);

        // find this cluster index in the fixed clusters list
        const auto clusterIterator = std::find(fixedClusters.begin(), fixedClusters.end(), index);

        // if this cluster index was found
        if (clusterIterator != fixedClusters.end())
        {
          // remove the cluster index in the fixed clusters list
          fixedClusters.erase(fixedClusters.begin() +
                              std::distance(fixedClusters.begin(), clusterIterator));
        }

        // lower any higher-indexed clusters by one
        std::for_each(fixedClusters.begin(), fixedClusters.end(),
                      [&](auto &l_cluster)
                      { l_cluster -= (l_cluster > index ? 1 : 0); });

        // update inter-collection information
        keepFitting = true;
      } // end if to reject this cluster
    }   // end loop over the collections

    // status printout
    if (!keepFitting &&
        m_algorithmPrintLevel > 1)
    {
      std::cout << "Fit has completed" << std::endl;
    }

    // status printout
    if (!keepFitting &&
        m_algorithmPrintLevel > 2)
    {
      std::vector<std::vector<double>> data(parameters.size());
      for (std::size_t i = 0; i < data.size(); ++i)
      {
        const auto &scales = linearParameters.at(i).second;
        const auto &cluster = p_clusters.at(catalog.at(i)).clusters.at(clusterIndex.at(i));

        data.at(i).push_back(i);
        data.at(i).push_back(parameters.at(i));
        data.at(i).push_back(parameters.at(i) - cluster.time);
        data.at(i).push_back(std::accumulate(scales.begin(), scales.end(), 0.0));
        data.at(i).push_back(std::accumulate(scales.begin(), scales.end(), 0.0) -
                             std::accumulate(cluster.scales.begin(), cluster.scales.end(), 0.0));
      }

      Print::Table(
          {2, 10, 6, 8, 7},
          {"#", "Time", "", "Energy", ""},
          {"%.0f", "%.3f", "%.3f", "%.3f", "%.3f"}, data);
    } // end status printout

  } while (keepFitting);

  // return the cluster collections
  return clusterCollections;
} // end function GlobalFitAlgorithm::Fit

std::vector<std::pair<int, double>> GlobalFitAlgorithm::FindCrystalWindow(
    const ModelFunction &p_model,
    const std::vector<std::vector<double>> &p_traces,
    const std::vector<double> &p_parameters,
    const std::vector<int> &p_masked) const
{
  // the metric for each crystal
  std::vector<double> metrics(p_traces.size(), 0.0);

  // loop over the metrics
  for (std::size_t i = 0; i < metrics.size(); ++i)
  {
    // if a masked crystal
    if (std::find(p_masked.begin(), p_masked.end(), i) == p_masked.end() &&
        p_traces.at(i).size() > 0)
    {

      // calculate this crystal's metric
      metrics.at(i) =
          *std::max_element(p_traces.at(i).begin(), p_traces.at(i).end()) -
          NormalOrder::ExpectedValue(
              p_model.TraceSampleCount(p_parameters, i), m_modelParameters.noiseLevels.at(i));

      // if the metric is negative and not allowed
      if (metrics.at(i) < 0.0 &&
          !m_allowNegativeMetrics)
      {
        // reset metric to zero
        metrics.at(i) = 0.0;
      } // end if the metric is negative and not allowed

    } // end if a masked crystal
  }   // end loop over the metrics

  // window grid size
  const int windowRows = m_windowRows.first + m_windowRows.second + 1;
  const int windowColumns = m_windowColumns.first + m_windowColumns.second + 1;

  // largest window's metric sum
  std::pair<int, double> largestMetricSum;

  // the window's central crystal
  int centralCrystal = 0;

  // loop over the crystal rows
  for (int i_row = 0; i_row <= m_crystalRows - windowRows; ++i_row)
  {
    // loop over the crystal columns
    for (int i_column = 0; i_column <= m_crystalColumns - windowColumns; ++i_column)
    {

      // this crystal index
      const int crystalIndex = m_crystalColumns * i_row + i_column;

      // initialize the window's metric sum
      double metricSum = 0.0;

      // initialize the window's largest metric
      std::pair<int, double> largestMetric(0, 0.0);

      // loop over the window rows
      for (int i = 0; i < windowRows; ++i)
      {
        // loop over the window columns
        for (int j = 0; j < windowColumns; ++j)
        {

          // this metric value
          const double metric = metrics.at(crystalIndex + m_crystalColumns * i + j);

          // add this crystal metric to sum
          metricSum += metric;

          // if window's largest metric yet
          if ((i == 0 && j == 0) ||
              (largestMetric.second < metric))
          {
            // update the largest metric
            largestMetric.first = crystalIndex + m_crystalColumns * i + j;
            largestMetric.second = metric;
          } // end if window's largest metric yet

        } // end loop over the window columns
      }   // end loop over the window rows

      // if larger than largest metric sum
      if (largestMetricSum.second < metricSum)
      {
        // update the largest metric sum
        largestMetricSum.first = crystalIndex;
        largestMetricSum.second = metricSum;

        // update the central crystal
        centralCrystal = largestMetric.first;
      } // end if larger than largest metric sum

    } // end loop over the crystal columns
  }   // end loop over the crystal rows

  // initialize the crystal window
  std::vector<std::pair<int, double>> window;

  // location of the crystal with the largest metric
  const std::pair<int, int> metricLocation = {std::floor(centralCrystal / m_crystalColumns),
                                              m_crystalColumns - 1 - (centralCrystal % m_crystalColumns)};

  // loop over the crystals
  for (int i_crystal = 0; i_crystal < m_crystalRows * m_crystalColumns; ++i_crystal)
  {
    // skip this crystal if missing
    if (p_traces.at(i_crystal).empty())
    {
      continue;
    }

    // this crystal's location
    const std::pair<int, int> location = {std::floor(i_crystal / m_crystalColumns),
                                          m_crystalColumns - 1 - (i_crystal % m_crystalColumns)};

    // if crystal is within the window bounds
    if (location.first - metricLocation.first >= -1.0 * m_windowRows.first &&
        location.first - metricLocation.first <= 1.0 * m_windowRows.second &&
        location.second - metricLocation.second >= -1.0 * m_windowColumns.first &&
        location.second - metricLocation.second <= 1.0 * m_windowColumns.second)
    {
      // update the crystal window
      // note: central crystal is flagged with a negative index
      window.emplace_back(
          (i_crystal == centralCrystal ? -1.0 * i_crystal : i_crystal),
          metrics.at(i_crystal));
    } // end if crystal is within the window bounds
  }   // end loop over the crystals

  // status printout
  if (m_algorithmPrintLevel > 3)
  {
    std::cout << "Crystal metrics:" << std::endl;

    // extract the crystals in the window
    std::vector<int> crystals(window.size());
    std::transform(window.begin(), window.end(), crystals.begin(),
                   [](const auto &l_crystal)
                   { return l_crystal.first; });

    // display the matrix
    Print::Matrix(m_crystalRows, m_crystalColumns, metrics, crystals, p_masked);
  } // end status printout

  // return the found crystal window
  return window;
} // end function GlobalFitAlgorithm::FindCrystalWindow

double GlobalFitAlgorithm::EstimateTime(
    const std::vector<std::vector<double>> &p_traces,
    const std::vector<std::vector<double>> &p_indices,
    const int &p_crystal) const
{
  // this crystal's trace
  const auto &trace = p_traces.at(p_crystal);

  // peak time estimate
  double estimate = p_indices.at(p_crystal).at(std::distance(trace.begin(),
                                                             std::max_element(trace.begin(), trace.end())));

  // if to apply the pseudo-time correction
  if (m_refineTimeEstimate)
  {
    estimate += PseudotimeCorrection(trace, p_crystal);
  }

  // return estimated peak time
  return estimate;
} // end function GlobalFitAlgorithm::EstimateTime

double GlobalFitAlgorithm::PseudotimeCorrection(
    const std::vector<double> &p_trace,
    const int &p_crystal) const
{
  // index of the peak value
  const std::size_t peakIndex = std::distance(p_trace.begin(),
                                              std::max_element(p_trace.begin(), p_trace.end()));

  // if the first or last sample is the peak
  if (peakIndex == 0 ||
      peakIndex == p_trace.size() - 1)
  {
    // do not apply a psuedotime correction
    return 0.0;
  }

  // calculate the psuedotime
  const double psuedotime = (p_trace.at(peakIndex) == p_trace.at(peakIndex + 1))
                                ? 1.0
                                : 2.0 / M_PI * std::atan((p_trace.at(peakIndex) - p_trace.at(peakIndex - 1)) / (p_trace.at(peakIndex) - p_trace.at(peakIndex + 1)));

  // return the mapped real time
  return m_pseudotimes.at(p_crystal)->Eval(psuedotime) - 0.5;
  // return 0;
} // end function GlobalFitAlgorithm::PseudotimeCorrection

std::vector<double> GlobalFitAlgorithm::TimeIndexList(
    const std::vector<double> &p_trace,
    const std::vector<double> &p_index) const
{
  // index of the peak value
  const double peakIndex = p_index.at(std::distance(p_trace.begin(),
                                                    std::max_element(p_trace.begin(), p_trace.end())));

  // first and last sample index
  double beginIndex = peakIndex - m_fitSamples.first;
  double endIndex = peakIndex + m_fitSamples.second;

  // if all available pre-/post-samples requested
  beginIndex = (m_fitSamples.first < 0) ? p_index.front() : beginIndex;
  endIndex = (m_fitSamples.second < 0) ? p_index.back() : endIndex;

  // consider the traces' lower and upper limits
  beginIndex = (beginIndex < p_index.front()) ? p_index.front() : beginIndex;
  endIndex = (endIndex > p_index.back()) ? p_index.back() : endIndex;

  // sequential indices from first to last value
  std::vector<double> samples(endIndex - beginIndex + 1);
  std::iota(samples.begin(), samples.end(), beginIndex);

  // return the list of samples
  return samples;
} // end function GlobalFitAlgorithm::TimeIndexList

std::tuple<bool, bool, int> GlobalFitAlgorithm::MergeClusters(
    const int &p_pedestal,
    const int &p_offset,
    std::vector<Cluster> &p_clusters,
    ModelFunction &p_model) const
{
  // pair of clusters to be merged
  std::pair<int, int> closeClusters = {0, 0};

  // window overlap between clusters to be merged
  bool shareWindow = false;
  bool shareCrystals = false;

  // loop over the clusters
  for (std::size_t i = 0; i < p_clusters.size(); ++i)
  {
    // the first cluster in the pair
    const auto &firstCluster = p_clusters.at(i);

    // loop over the remaining clusters
    for (std::size_t j = i + 1; j < p_clusters.size(); ++j)
    {
      // the second cluster in the pair
      const auto &secondCluster = p_clusters.at(j);

      // cluster separation in time
      const double separation = std::abs(secondCluster.time - firstCluster.time);

      // if clusters are within the dead time
      if (separation <= m_artificialDeadTime)
      {
        // crystal count shared between the clusters
        std::size_t crystalCount = 0;

        // loop over the first cluster in the pair's crystals
        for (const auto &crystal : firstCluster.crystals)
        {
          // if this crystal is shared
          if (std::find(secondCluster.crystals.begin(), secondCluster.crystals.end(),
                        crystal) != secondCluster.crystals.end())
          {

            // increment shared crystal count
            ++crystalCount;

          } // end if this crystal is shared
        }   // end loop over the first cluster in the pair's crystals

        // if their windows don't overlap
        if (crystalCount == 0)
        {
          // nothing to merge
          continue;
        }

        // status printout
        if (m_algorithmPrintLevel > 1)
        {
          std::cout << Form("Clusters # %lu and %lu are separated by %.5f c.t.",
                            p_offset + i, p_offset + j, separation)
                    << std::endl;
        }

        // previous cluster separation in time
        const double previousSeparation = std::abs(
            p_clusters.at(closeClusters.second).time -
            p_clusters.at(closeClusters.first).time);

        // complete window overlap
        const bool completeOverlap =
            (firstCluster.crystals.size() == crystalCount) ||
            (secondCluster.crystals.size() == crystalCount);

        // if a closer cluster pair
        if (!shareCrystals ||
            (separation < previousSeparation && !shareWindow) ||
            (separation < previousSeparation && shareWindow && completeOverlap))
        {
          // update cluster information
          closeClusters.first = i;
          closeClusters.second = j;

          // update window information
          shareWindow = completeOverlap;
          shareCrystals = true;
        } // end if a closer cluster pair
      }   // end if clusters are within the dead time

    } // end loop over the remaining clusters
  }   // end loop over the clusters

  // if the clusters should be merged
  if (shareCrystals && !shareWindow)
  {
    // status printout
    if (m_algorithmPrintLevel > 1)
    {
      std::cout << "Merging cluster # " << p_offset + closeClusters.second
                << " into # " << p_offset + closeClusters.first << std::endl;
    }

    // extend the first cluster's model
    p_model.ExtendCluster(p_pedestal, p_offset + closeClusters.first,
                          p_clusters.at(closeClusters.second).crystals,
                          p_clusters.at(closeClusters.first).samples);

    // remove the second cluster's model
    p_model.RemoveCluster(p_offset + closeClusters.second);

    // extend the first cluster's crystals
    p_clusters.at(closeClusters.first).crystals =
        p_model.MergeGroups({p_clusters.at(closeClusters.first).crystals,
                             p_clusters.at(closeClusters.second).crystals});

    // sort the updated crystals list
    std::sort(
        p_clusters.at(closeClusters.first).crystals.begin(),
        p_clusters.at(closeClusters.first).crystals.end());

    // extend the first cluster's central crystals
    p_clusters.at(closeClusters.first).centers.insert(p_clusters.at(closeClusters.first).centers.end(), p_clusters.at(closeClusters.second).centers.begin(), p_clusters.at(closeClusters.second).centers.end());

    // remove the second cluster's object
    p_clusters.erase(p_clusters.begin() + closeClusters.second);
  } // end if the clusters should be merged

  // return merging information
  return {shareCrystals && !shareWindow, shareWindow, closeClusters.second};
} // end function GlobalFitAlgorithm::MergeClusters

int GlobalFitAlgorithm::NoiseIndex(
    const int &p_offset,
    const std::vector<Cluster> &p_clusters,
    const ModelFunction &p_model) const
{
  // noise cluster index and scale sum
  int noiseIndex = -1;
  double noiseTotal = 0.0;

  // loop over the clusters
  for (std::size_t i = 0; i < p_clusters.size(); ++i)
  {
    // this cluster
    const auto &cluster = p_clusters.at(i);

    // template scale sum threshold
    const double sumThreshold =
        std::accumulate(cluster.crystals.begin(), cluster.crystals.end(), 0.0,
                        [&](const auto &l_value, const auto &l_crystal)
                        { return l_value + NormalOrder::ExpectedValue(
                                               p_model.FitSampleCount(i + p_offset, l_crystal),
                                               m_modelParameters.noiseLevels.at(l_crystal)); });

    // the template scale sum
    const double scaleSum = std::accumulate(
        cluster.scales.begin(), cluster.scales.end(), 0.0);

    // if the template scale sum is too small
    if (scaleSum < sumThreshold && (noiseIndex == -1 || scaleSum < noiseTotal))
    {
      // update template scale sum information
      noiseIndex = i;
      noiseTotal = scaleSum;
    }
  } // end loop over the clusters

  // status printout
  if (noiseIndex >= 0 &&
      m_algorithmPrintLevel > 1)
  {
    std::cout << Form("Cluster # %i has a template scale sum of %.3f ADU",
                      p_offset + noiseIndex, std::accumulate(p_clusters.at(noiseIndex).scales.begin(), p_clusters.at(noiseIndex).scales.end(), 0.0))
              << std::endl;
  }

  // return noise fit information
  return noiseIndex;
} // end function GlobalFitAlgorithm::NoiseIndex
