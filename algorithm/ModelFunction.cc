/*! @file ModelFunction.cc

    @author  David A. Sweigart
    @date    01/25/2018
    @version 10/01/2018

    Copyright (C) 2018 David A. Sweigart.  All rights reserved.
 */

#include "ModelFunction.hh"

ModelFunction::ModelFunction(
  const int                                          &p_pedestalCount,
  const std::vector<std::vector<double>>             &p_traces,
  const std::vector<std::vector<double>>             &p_indices,
  const std::vector<std::shared_ptr<const TSpline3>> &p_templates,
  const std::vector<double>                          &p_noiseLevels,
  const double                                       &p_artificialDeadTime)
    : m_traces            (p_traces                  ),
      m_indices           (p_indices                 ),
      m_templates         (p_templates               ),
      m_noiseLevels       (p_noiseLevels             ),
      m_artificialDeadTime(p_artificialDeadTime      ),
      m_pedestalCount     (p_pedestalCount           ),
      m_clusterCount      (0                         ),
      m_crystals          (m_traces.size(), Crystal()) {
  // initialize pedestal indices to zero
  for (const auto &indices : m_indices) {
    m_pedestals.emplace_back(indices.size(), 0);
  }
} // end function ModelFunction::ModelFunction

double ModelFunction::operator()(const std::vector<double> &p_parameters) const {
  // initialize the chi-squared to zero
  double chiSquared = 0.0;

  // loop over the crystals
  for (std::size_t i = 0; i < m_crystals.size(); ++i) {
    // add this crystal's chi-squared term
    chiSquared += CrystalChiSquared(p_parameters, i, m_crystals.at(i).samples);
  }

  // return the chi-squared
  return chiSquared;
} // end function ModelFunction::operator

double ModelFunction::CrystalChiSquared(
    const std::vector<double> &p_parameters,
    const int                 &p_crystal,
    const std::vector<int>    &p_samples) const {
  // initialize the chi-squared to zero
  double chiSquared = 0.0;

  // if no samples for this crystal
  if (p_samples.size() == 0) {
    // return zero
    return chiSquared;
  }

  // solve for the template scales and pedestal
  const auto X = AnalyticalMinimization(p_parameters, p_crystal);

  // loop over the samples
  for (const auto &sample : p_samples) {
    // initialize model value to the pedestal
    double model = X(m_pedestals.at(p_crystal).at(sample));

    // loop over the clusters
    for (std::size_t i = 0; i < m_crystals.at(p_crystal).clusters.size(); ++i) {
      model += X(i + m_pedestalCount) * Evaluate(p_crystal,
        m_indices.at(p_crystal).at(sample),
        p_parameters.at(m_crystals.at(p_crystal).clusters.at(i)));
    } // end loop over the clusters

    // add the chi-squared term
    chiSquared += std::pow((m_traces.at(p_crystal).at(sample) - model) / Weight(p_crystal), 2);
  } // end loop over the samples

  // return the chi-squared
  return chiSquared;
} // end function ModelFunction::CrystalChiSquared

void ModelFunction::ExtendCluster(
    const int                 &p_pedestal,
    const int                 &p_cluster,
    const std::vector<int>    &p_crystals,
    const std::vector<double> &p_samples) {
  // loop over the crystals
  for (std::size_t i_crystal = 0; i_crystal < p_crystals.size(); ++i_crystal) {
    // this crystal
    auto &crystal = m_crystals.at(p_crystals.at(i_crystal));

    // if this is a new cluster index
    if (std::find(crystal.clusters.begin(), crystal.clusters.end(),
        p_cluster) == crystal.clusters.end()) {
      // include the new cluster index for this crystal
      crystal.clusters.push_back(p_cluster);

      // include the new pedestal index for this crystal
      crystal.pedestals.push_back(p_pedestal);

      // vector indices of new samples
      std::vector<int> sampleGroup;

      // the current time indices
      const auto &indices = m_indices.at(p_crystals.at(i_crystal));

      // loop over new time indices
      for (const auto &sample : p_samples) {
        // find this time index in the samples list
        const auto iterator = std::find_if(indices.begin(), indices.end(),
          [&](const auto &l_index) {
            return std::round(l_index) == std::round(sample);
          });

        // if this time index was found
        if (iterator != indices.end()) {
          const int index = std::distance(indices.begin(), iterator);

          // include this vector index
          sampleGroup.push_back(index);

          // include the new pedestal index for this sample
          m_pedestals.at(p_crystals.at(i_crystal)).at(index) = p_pedestal;
        }
      } // end loop over new time indices

      // include this sample group
      crystal.sampleGroups.push_back(sampleGroup);

      // update cluster-combined sample list
      crystal.samples = MergeGroups(crystal.sampleGroups);
    } // end if this is a new cluster index

  } // end loop over the crystals
} // end function ModelFunction::ExtendCluster

void ModelFunction::RemoveCluster(const int &p_cluster) {
  // loop over the crystals
  for (std::size_t i_crystal = 0; i_crystal < m_crystals.size(); ++i_crystal) {
    // this crystal
    auto &crystal = m_crystals.at(i_crystal);

    // find index of cluster to be removed
    const std::size_t index = std::distance(crystal.clusters.begin(),
      std::find(crystal.clusters.begin(), crystal.clusters.end(), p_cluster));

    // whether the cluster index exists
    const bool exists = (index != crystal.clusters.size());

    // if the cluster index exists
    if (exists) {
      // remove provided cluster index
      crystal.clusters    .erase(crystal.clusters    .begin() + index);
      crystal.pedestals   .erase(crystal.pedestals   .begin() + index);
      crystal.sampleGroups.erase(crystal.sampleGroups.begin() + index);
    }

    // lower any higher-indexed clusters by one
    std::for_each(crystal.clusters.begin(), crystal.clusters.end(),
      [&](auto &l_cluster) { l_cluster -= (l_cluster > p_cluster ? 1 : 0); });

    // if the cluster index doesn't exist
    if (!exists) {
      // the cluster to be removed isn't in this crystal
      continue;
    }

    // update cluster-combined sample list
    crystal.samples = MergeGroups(crystal.sampleGroups);
  } // end loop over the crystals

  // decrement the global cluster count
  --m_clusterCount;
} // end function ModelFunction::RemoveCluster

std::vector<std::vector<double>> ModelFunction::Residuals(
    const std::vector<double> &p_parameters) const {
  // make a copy of the traces
  std::vector<std::vector<double>> residuals(m_traces.size());

  // loop over the crystals
  for (std::size_t i_crystal = 0; i_crystal < m_traces.size(); ++i_crystal) {
    // solve for the template scales and pedestal
    const auto X = AnalyticalMinimization(p_parameters, i_crystal);

    // loop over the samples
    for (std::size_t i_sample = 0; i_sample < m_traces.at(i_crystal).size(); ++i_sample) {

      // if this residual should be ignored
      if (MaskedIndex(p_parameters, i_crystal, i_sample)) {
        // zero out the residual
        residuals.at(i_crystal).push_back(0);
      } else {
        // initialize model value to the pedestal
        double model = X(m_pedestals.at(i_crystal).at(i_sample));

        // loop over the clusters
        for (std::size_t i = 0; i < m_crystals.at(i_crystal).clusters.size(); ++i) {
          // include contribution to the model value
          model += X(i + m_pedestalCount) * Evaluate(i_crystal,
            m_indices.at(i_crystal).at(i_sample),
            p_parameters.at(m_crystals.at(i_crystal).clusters.at(i)));
        } // end loop over the clusters

        // subtract the model value
        residuals.at(i_crystal).push_back(m_traces.at(i_crystal).at(i_sample) - model);
      } // end if this residual should be ignored

    } // end loop over the samples
  } // end loop over the crystals

  // return the residuals list
  return residuals;
} // end function ModelFunction::Residuals

int ModelFunction::TraceSampleCount(
    const std::vector<double> &p_parameters,
    const int                 &p_crystal) const {
  // initialize this crystal's sample count
  int sampleCount = 0;

  // loop over the samples
  for (std::size_t i_sample = 0; i_sample < m_traces.at(p_crystal).size(); ++i_sample) {
    // if this sample is not near a cluster
    if (!MaskedIndex(p_parameters, p_crystal, i_sample)) {
      // increment the sample count
      ++sampleCount;
    }
  } // end loop over the samples

  // return this crystal's sample count
  return sampleCount;
} // end function ModelFunction::TraceSampleCount

std::vector<std::pair<std::vector<double>, std::vector<double>>> ModelFunction::LinearParameters(
    const std::vector<double> &p_parameters) const {
  // pedestal and template scale for each crystal for each cluster
  std::vector<std::pair<std::vector<double>, std::vector<double>>> parameters(m_clusterCount);

  // loop over the crystals
  for (std::size_t i_crystal = 0; i_crystal < m_crystals.size(); ++i_crystal) {
    // this crystal
    const auto &crystal = m_crystals.at(i_crystal);

    // if no clusters for this crystal
    if (crystal.clusters.size() == 0) {
      // skip this crystal
      continue;
    }

    // solve for the template scales and pedestal
    const auto X = AnalyticalMinimization(p_parameters, i_crystal);

    // loop over the clusters
    for (std::size_t i = 0; i < crystal.clusters.size(); ++i) {
      // append the parameters
      parameters.at(crystal.clusters.at(i)).first .push_back(X(crystal.pedestals.at(i)));
      parameters.at(crystal.clusters.at(i)).second.push_back(X(i + m_pedestalCount));
    }
  } // end loop over the crystals

  // return the parameter list
  return parameters;
} // end function ModelFunction::LinearParameters

Eigen::MatrixXd ModelFunction::AnalyticalMinimization(
    const std::vector<double> &p_parameters,
    const int                 &p_crystal) const {
  // analytical minimization matrices
  // note: L L^T X = L Y
  Eigen::MatrixXd L; // models matrix
  Eigen::MatrixXd Y; // sample matrix

  // the sample length
  const int length = (m_crystals.at(p_crystal).clusters.size() == 0)
    ? m_traces  .at(p_crystal)        .size()
    : m_crystals.at(p_crystal).samples.size();

  // resize matrices for the cluster count and sample length
  L.resize(m_crystals.at(p_crystal).clusters.size() + m_pedestalCount, length);
  Y.resize(length, 1);

  // loop over the columns
  for (int i = 0; i < L.cols(); ++i) {
    // loop over the rows
    for (int j = 0; j < L.rows(); ++j) {

      // update matrix
      L(j, i) = (j < m_pedestalCount)
        ? (m_crystals .at(p_crystal).clusters.size() > 0 &&
           m_pedestals.at(p_crystal).at(m_crystals.at(p_crystal).samples.at(i)) != j)
          ? 0.0
          : 1.0 / Weight(p_crystal)
        : Evaluate(p_crystal,
            m_indices.at(p_crystal).at(m_crystals.at(p_crystal).samples.at(i)),
            p_parameters.at(m_crystals.at(p_crystal).clusters.at(j - m_pedestalCount)))
          / Weight(p_crystal);

    } // end loop over the rows

    // update matrix
    Y(i, 0) = (m_crystals.at(p_crystal).clusters.size() == 0
      ? m_traces.at(p_crystal).at(i)
      : m_traces.at(p_crystal).at(m_crystals.at(p_crystal).samples.at(i)))
      / Weight(p_crystal);
  } // end loop over the columns

  // solve for the linear parameters
  return (L * L.transpose())
    .ldlt()
    .solve(L * Y);
} // end function ModelFunction::AnalyticalMinimization

double ModelFunction::Evaluate(
    const int    &p_crystal,
    const double &p_sample,
    const double &p_peakTime) const {
  // return the model function value
  return ((p_sample - p_peakTime) >= m_templates.at(p_crystal)->GetXmin() &&
          (p_sample - p_peakTime) <= m_templates.at(p_crystal)->GetXmax())
    ? m_templates.at(p_crystal)->Eval(p_sample - p_peakTime) /
      m_templates.at(p_crystal)->Eval(0)
    : 0.0;
} // end function ModelFunction::Evaluate

double ModelFunction::Weight(const int &p_crystal) const {
  // return this crystal's noise level
  return m_noiseLevels.at(p_crystal);
} // end function ModelFunction::Weight
