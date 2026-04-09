#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct DataPoint {
  double t9 = 0.0;
  double rate = 0.0;
};

struct FitTerm {
  double a[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

struct FitSummary {
  std::vector<FitTerm> terms;
};

std::vector<FitSummary> gPlotSummaryStorage;
const std::vector<FitSummary> *gPlotSummaries = nullptr;

std::string trim_copy(const std::string &s) {
  const std::string whitespace = " \t\r\n";
  const std::string::size_type first = s.find_first_not_of(whitespace);
  if (first == std::string::npos) {
    return "";
  }
  const std::string::size_type last = s.find_last_not_of(whitespace);
  return s.substr(first, last - first + 1);
}

std::vector<std::string> split_csv(const std::string &text) {
  std::vector<std::string> out;
  std::string item;
  std::istringstream iss(text);
  while (std::getline(iss, item, ',')) {
    item = trim_copy(item);
    if (!item.empty()) {
      out.push_back(item);
    }
  }
  return out;
}

std::string base_name(const std::string &path) {
  const std::string::size_type slash = path.find_last_of("/\\");
  return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

std::vector<DataPoint> load_input_table(const std::string &path) {
  std::ifstream in(path.c_str());
  if (!in) {
    throw std::runtime_error("Failed to open input file: " + path);
  }

  std::vector<DataPoint> data;
  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string stripped = trim_copy(line);
    if (stripped.empty() || stripped[0] == '#') {
      continue;
    }
    std::istringstream iss(stripped);
    DataPoint p;
    if (!(iss >> p.t9 >> p.rate)) {
      throw std::runtime_error("Failed to parse input line " + std::to_string(lineNo) + " in " +
                               path);
    }
    if (p.t9 > 0.0 && p.rate > 0.0) {
      data.push_back(p);
    }
  }
  return data;
}

FitSummary load_fit_summary(const std::string &path) {
  std::ifstream in(path.c_str());
  if (!in) {
    throw std::runtime_error("Failed to open fit summary: " + path);
  }

  FitSummary summary;
  bool inFittedBlock = false;
  std::string line;
  while (std::getline(in, line)) {
    const std::string stripped = trim_copy(line);
    if (stripped.empty()) {
      continue;
    }
    if (stripped[0] == '#') {
      if (stripped == "# Fitted parameters") {
        inFittedBlock = true;
      } else if (inFittedBlock) {
        break;
      }
      continue;
    }
    if (!inFittedBlock) {
      continue;
    }

    std::istringstream iss(stripped);
    FitTerm term;
    std::string kind;
    if (!(iss >> term.a[0] >> term.a[1] >> term.a[2] >> term.a[3] >> term.a[4] >> term.a[5] >>
          term.a[6] >> kind)) {
      throw std::runtime_error("Failed to parse fitted parameter line in " + path);
    }
    summary.terms.push_back(term);
  }

  if (summary.terms.empty()) {
    throw std::runtime_error("No fitted parameter lines found in " + path);
  }
  return summary;
}

double reaclib_rate(double t9, const double *a) {
  const double exponent = a[0] + a[1] * std::pow(t9, -1.0) + a[2] * std::pow(t9, -1.0 / 3.0) +
                          a[3] * std::pow(t9, 1.0 / 3.0) + a[4] * t9 +
                          a[5] * std::pow(t9, 5.0 / 3.0) + a[6] * std::log(t9);
  if (exponent > 700.0) {
    return std::numeric_limits<double>::infinity();
  }
  return std::exp(exponent);
}

double summed_rate(double t9, const FitSummary &summary) {
  double sum = 0.0;
  for (size_t i = 0; i < summary.terms.size(); ++i) {
    const double value = reaclib_rate(t9, summary.terms[i].a);
    if (!std::isfinite(value)) {
      return std::numeric_limits<double>::infinity();
    }
    sum += value;
  }
  return sum;
}

double evaluate_fit_tf1(double *x, double *par) {
  if (gPlotSummaries == nullptr) {
    return 0.0;
  }
  const int index = static_cast<int>(par[0]);
  if (index < 0 || index >= static_cast<int>(gPlotSummaries->size())) {
    return 0.0;
  }
  return summed_rate(x[0], (*gPlotSummaries)[index]);
}

std::vector<double> make_log_grid(double xmin, double xmax, int n = 500) {
  std::vector<double> x;
  x.reserve(n);
  const double logMin = std::log10(xmin);
  const double logMax = std::log10(xmax);
  for (int i = 0; i < n; ++i) {
    const double f = (n == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n - 1);
    x.push_back(std::pow(10.0, logMin + f * (logMax - logMin)));
  }
  return x;
}

}  // namespace

void plot_fit_collection(const char *inputCsv, const char *outputCsv, double xmin = 0.1,
                         double xmax = 10.0) {
  const std::vector<std::string> inputPaths = split_csv(inputCsv);
  const std::vector<std::string> outputPaths = split_csv(outputCsv);
  if (inputPaths.empty()) {
    throw std::runtime_error("No input files were provided");
  }
  if (inputPaths.size() != outputPaths.size()) {
    throw std::runtime_error("input/output list sizes must match");
  }
  if (xmin <= 0.0 || xmax <= xmin) {
    throw std::runtime_error("Require 0 < xmin < xmax");
  }

  gStyle->SetOptStat(0);
  auto *canvas = new TCanvas("c_fit_collection", "Fit collection", 1000, 750);
  canvas->SetLogx();
  canvas->SetLogy();

  const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7, kMagenta + 1};
  const int markers[] = {20, 21, 22, 23, 33, 34};

  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::lowest();
  std::vector<TGraph *> inputGraphs;
  std::vector<TF1 *> fitFuncs;
  std::vector<std::string> inputLabels;
  std::vector<std::string> outputLabels;
  gPlotSummaryStorage.clear();

  for (size_t i = 0; i < inputPaths.size(); ++i) {
    const std::vector<DataPoint> data = load_input_table(inputPaths[i]);
    const FitSummary summary = load_fit_summary(outputPaths[i]);
    gPlotSummaryStorage.push_back(summary);

    std::vector<double> xData;
    std::vector<double> yData;
    xData.reserve(data.size());
    yData.reserve(data.size());
    for (size_t j = 0; j < data.size(); ++j) {
      xData.push_back(data[j].t9);
      yData.push_back(data[j].rate);
      ymin = std::min(ymin, data[j].rate);
      ymax = std::max(ymax, data[j].rate);
    }

    const std::vector<double> xFit = make_log_grid(xmin, xmax, 600);
    for (size_t j = 0; j < xFit.size(); ++j) {
      const double y = summed_rate(xFit[j], gPlotSummaryStorage.back());
      if (std::isfinite(y) && y > 0.0) {
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
      }
    }

    auto *gIn = new TGraph(static_cast<int>(xData.size()), xData.data(), yData.data());

    const int color = colors[i % (sizeof(colors) / sizeof(colors[0]))];
    const int marker = markers[i % (sizeof(markers) / sizeof(markers[0]))];
    gIn->SetLineColor(color);
    gIn->SetMarkerColor(color);
    gIn->SetMarkerStyle(marker);
    gIn->SetMarkerSize(0.7);
    gIn->SetLineWidth(2);

    inputGraphs.push_back(gIn);
    inputLabels.push_back(base_name(inputPaths[i]));
    outputLabels.push_back(base_name(outputPaths[i]));
  }

  gPlotSummaries = &gPlotSummaryStorage;

  if (!std::isfinite(ymin) || !std::isfinite(ymax) || ymin <= 0.0 || ymax <= ymin) {
    throw std::runtime_error("Failed to determine y-range for plotting");
  }

  const double logYMin = std::log10(ymin);
  const double logYMax = std::log10(ymax);
  const double logYSpan = std::max(1e-6, logYMax - logYMin);
  const double yPlotMin = std::pow(10.0, logYMin - 0.08 * logYSpan);
  const double yPlotMax = std::pow(10.0, logYMax + 0.08 * logYSpan);

  inputGraphs[0]->SetTitle("");
  inputGraphs[0]->GetXaxis()->SetLimits(xmin, xmax);
  inputGraphs[0]->GetYaxis()->SetRangeUser(yPlotMin, yPlotMax);
  inputGraphs[0]->GetXaxis()->SetTitle("T_{9} [GK]");
  inputGraphs[0]->GetYaxis()->SetTitle("Reaction rate");
  //inputGraphs[0]->Draw("APL");
  inputGraphs[0]->Draw("AL");

  for (size_t i = 0; i < gPlotSummaryStorage.size(); ++i) {
    const std::string funcName = "fit_func_" + std::to_string(i);
    auto *fitFunc = new TF1(funcName.c_str(), evaluate_fit_tf1, xmin, xmax, 1);
    fitFunc->SetParameter(0, static_cast<double>(i));
    fitFunc->SetNpx(1000);
    fitFunc->SetLineColor(colors[i % (sizeof(colors) / sizeof(colors[0]))]);
    fitFunc->SetLineWidth(3);
    fitFunc->SetLineStyle(2);
    fitFuncs.push_back(fitFunc);
  }

  fitFuncs[0]->Draw("SAME");

  for (size_t i = 1; i < inputGraphs.size(); ++i) {
    //inputGraphs[i]->Draw("PL SAME");
    inputGraphs[i]->Draw("L SAME");
    fitFuncs[i]->Draw("SAME");
  }

  const double y2 = std::min(0.90, 0.12 + 0.10 * static_cast<double>(inputLabels.size()));
  auto *legend = new TLegend(0.52, 0.12, 0.92, y2);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  for (size_t i = 0; i < inputLabels.size(); ++i) {
    //legend->AddEntry(inputGraphs[i], inputLabels[i].c_str(), "lp");
    legend->AddEntry(inputGraphs[i], inputLabels[i].c_str(), "l");
    legend->AddEntry(fitFuncs[i], outputLabels[i].c_str(), "l");
  }
  legend->Draw();

  canvas->Update();
}
