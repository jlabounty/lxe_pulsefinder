#include "TFile.h"
#include "TSpline.h"
#include <iostream>


class Parameters{
    public:
        typedef struct {
            const int                 crystalRows;
            const int                 crystalColumns;
            const std::pair<int, int> windowRows;
            const std::pair<int, int> windowColumns;
            const std::vector<double> showerConstants;
            const double              clusterThreshold;
            const double              artificialDeadTime;
            const std::vector<int>    maxFunctionCalls;
            const double              minimizationTolerance;
            const double              minimizationPrecision;
            const int                 strategyLevel;
            const int                 minimizerPrintLevel;
            const int                 algorithmPrintLevel;
            const bool                refineTimeEstimate;
            const bool                allowNegativeMetrics;
            const std::pair<int, int> fitSamples;
            const std::pair<int, int> positionAlgorithm;
            const std::vector<double> seedTweaks;
            const double              chiSquaredTolerance;
            const double              limitTolerance;
        } AlgorithmParameters;
};

std::vector<TSpline3> GetTemplates(int size)
{
    std::vector<TSpline3> ding;
    TFile* f = TFile::Open("/home/jlab/g-2/docker/nearline_sl7/srcs/gm2calo/g2RunTimeFiles/templatesRun3/beamTemplates/calo24BeamTemplates/template5.root");
    TSpline3* spl = (TSpline3*)f->Get("masterSpline");
    for(int i = 0; i < size; i++)
    {
        ding.push_back(*spl);
        std::cout << ding[i].GetName() << std::endl;
    }
    f->Close();
    return ding;
}
