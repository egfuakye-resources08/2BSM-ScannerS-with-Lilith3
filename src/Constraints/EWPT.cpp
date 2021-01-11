#include "ScannerS/Constraints/EWPT.hpp"
#include "BSMPT/minimizer/Minimizer.h"
#include "BSMPT/models/ClassPotentialOrigin.h"
#include "ScannerS/Utilities.hpp"
#include <cstddef>

using namespace std::string_literals;

namespace ScannerS::Constraints::EWPTDetail {

DataMap::Map InitModel(std::shared_ptr<BSMPT::Class_Potential_Origin> model,
                       const std::vector<double> &input) {
  auto counterTerms = model->initModel(input);
  auto labels = model->addLegendCT();
  return Utilities::ZipToMap(std::move(labels), counterTerms);
}

DataMap::Map
TripleHiggsCoups(std::shared_ptr<BSMPT::Class_Potential_Origin> model) {
  model->Prepare_Triple();
  model->TripleHiggsCouplings();
  std::vector<double> tripleHCoups{};
  for (size_t i = 0; i < model->get_NHiggs(); i++) {
    for (size_t j = i; j < model->get_NHiggs(); j++) {
      for (size_t k = j; k < model->get_NHiggs(); k++) {
        tripleHCoups.push_back(
            -model->get_TripleHiggsCorrectionsTreePhysical(i, j, k));
        tripleHCoups.push_back(
            -model->get_TripleHiggsCorrectionsCTPhysical(i, j, k));
        tripleHCoups.push_back(
            -model->get_TripleHiggsCorrectionsCWPhysical(i, j, k));
      }
    }
  }
  auto labels = model->addLegendTripleCouplings();
  return Utilities::ZipToMap(std::move(labels), tripleHCoups);
}

std::pair<DataMap::Map, bool>
NloVev(std::shared_ptr<BSMPT::Class_Potential_Origin> model) {
  std::vector<double> check, nloVev;
  nloVev = BSMPT::Minimizer::Minimize_gen_all(model, 0, check,
                                              model->get_vevTreeMin());
  auto labels = model->addLegendVEV();
  return {Utilities::ZipToMap(std::move(labels), nloVev),
          model->CheckNLOVEV(nloVev)};
}

DataMap::Map
PhaseTransition(std::shared_ptr<BSMPT::Class_Potential_Origin> model) {
  std::vector<double> solution;
  auto EWPTReturn = BSMPT::Minimizer::PTFinder_gen_all(model, 0, 300);
  solution.push_back(EWPTReturn.Tc);
  solution.push_back(EWPTReturn.vc);
  solution.push_back(static_cast<double>(EWPTReturn.StatusFlag));
  solution.insert(solution.end(), EWPTReturn.EWMinimum.begin(),
                  EWPTReturn.EWMinimum.end());
  auto labels = model->addLegendTemp();
  labels.at(0) = "EWPT_T_c"s;
  labels.at(1) = "EWPT_omega_c"s;
  labels.at(2) = "BSMPT_ok"s;
  auto result = Utilities::ZipToMap(std::move(labels), solution);
  return result;
}
} // namespace ScannerS::Constraints::EWPTDetail
