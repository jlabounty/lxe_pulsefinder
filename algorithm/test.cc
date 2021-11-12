#include <time.h>
#include "stdio.h"
#include <iostream>
#include "standalone/GlobalFitAlgorithm.hh"
// #include "standalone/Print.hh"
// #include "standalone/PositionAlgorithm.hh"
// #include "standalone/ModelFunction.hh"
// #include "standalone/NormalOrder.hh"
#include "standalone/Additions.hh"
// #include "Minuit2/MnPrint.h"

int main()
{
    std::cout << "hello world!" << std::endl;

    // copied from GlobalFit_module.cc
    const Parameters::AlgorithmParameters m_algorithmParameters{
        4,
        4,
        {1,1},
        {1,1},
        {0.2,0.2},
        50.0,
        1.0,
        {200,100,5},
        0.1,
        -1.0,
        0,
        0,
        0,
        true,
        true,
        {8,16},
        {2,1},
        {-1.0,1.0,-2.0,2.0},
        0.1,
        0.5
    };

    const std::vector<std::shared_ptr<const TSpline3>> templates;   // = GetTemplates(4);
    const std::vector<std::shared_ptr<const TSpline3>> pseudotimes; // = GetTemplates(4);

    int caloNum = 1;
    // std::vector<std::shared_ptr<const TSpline3>> m_templates, m_pseudotimes;
    std::vector<double> noiseLevels;


    GlobalFitAlgorithm caloAlgorithm(
        templates,
        pseudotimes,
        noiseLevels,
        m_algorithmParameters.crystalRows,
        m_algorithmParameters.crystalColumns,
        m_algorithmParameters.windowRows,
        m_algorithmParameters.windowColumns,
        m_algorithmParameters.showerConstants,
        m_algorithmParameters.clusterThreshold,
        m_algorithmParameters.artificialDeadTime,
        m_algorithmParameters.maxFunctionCalls,
        m_algorithmParameters.minimizationTolerance,
        m_algorithmParameters.minimizationPrecision,
        m_algorithmParameters.strategyLevel,
        m_algorithmParameters.minimizerPrintLevel,
        m_algorithmParameters.algorithmPrintLevel,
        m_algorithmParameters.refineTimeEstimate,
        m_algorithmParameters.allowNegativeMetrics,
        m_algorithmParameters.fitSamples,
        m_algorithmParameters.positionAlgorithm,
        m_algorithmParameters.seedTweaks,
        m_algorithmParameters.chiSquaredTolerance,
        m_algorithmParameters.limitTolerance);

    caloAlgorithm.ding = 0;
    std::cout << caloAlgorithm.ding << std::endl;

}