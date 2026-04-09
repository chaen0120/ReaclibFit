#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct RateVariantData {
  std::vector<double> xCentral;
  std::vector<double> yCentral;
  std::vector<double> xLower;
  std::vector<double> yLower;
  std::vector<double> xUpper;
  std::vector<double> yUpper;
};

std::string trim_copy(const std::string &s) {
  const std::string whitespace = " \t\r\n";
  const std::string::size_type first = s.find_first_not_of(whitespace);
  if (first == std::string::npos) {
    return "";
  }
  const std::string::size_type last = s.find_last_not_of(whitespace);
  return s.substr(first, last - first + 1);
}

std::string default_variant_path(const std::string &inputPath, const std::string &suffix) {
  const std::string::size_type slash = inputPath.find_last_of("/\\");
  const std::string dir =
      (slash == std::string::npos) ? "" : inputPath.substr(0, slash + 1);
  const std::string file =
      (slash == std::string::npos) ? inputPath : inputPath.substr(slash + 1);
  const std::string::size_type dot = file.find_last_of('.');
  const std::string stem = (dot == std::string::npos) ? file : file.substr(0, dot);
  return dir + stem + "_" + suffix + ".dat";
}

RateVariantData write_rate_variants_impl(const char *inputPath, const char *centralPath,
                                         const char *lowerPath, const char *upperPath) {
  std::ifstream in(inputPath);
  if (!in) {
    throw std::runtime_error(std::string("Failed to open input file: ") + inputPath);
  }

  std::ofstream centralOut(centralPath);
  std::ofstream lowerOut(lowerPath);
  std::ofstream upperOut(upperPath);
  if (!centralOut || !lowerOut || !upperOut) {
    throw std::runtime_error("Failed to open one of the output files");
  }

  centralOut << std::scientific << std::setprecision(6);
  lowerOut << std::scientific << std::setprecision(6);
  upperOut << std::scientific << std::setprecision(6);

  std::string line;
  int lineNo = 0;
  int written = 0;
  int skipped = 0;
  RateVariantData data;

  while (std::getline(in, line)) {
    ++lineNo;
    const std::string stripped = trim_copy(line);
    if (stripped.empty() || stripped[0] == '#') {
      continue;
    }

    std::istringstream iss(stripped);
    double t9 = 0.0;
    double rate = 0.0;
    double lowerErr = 0.0;
    double upperErr = 0.0;
    if (!(iss >> t9 >> rate >> lowerErr >> upperErr)) {
      std::cerr << "Skipping line " << lineNo
                << ": expected at least 4 numeric columns: T9 rate lowerErr upperErr\n";
      ++skipped;
      continue;
    }

    const double lowerRate = rate - lowerErr;
    const double upperRate = rate + upperErr;

    centralOut << t9 << "\t" << rate << "\n";
    data.xCentral.push_back(t9);
    data.yCentral.push_back(rate);
    if (lowerRate > 0.0) {
      lowerOut << t9 << "\t" << lowerRate << "\n";
      data.xLower.push_back(t9);
      data.yLower.push_back(lowerRate);
    } else {
      std::cerr << "Skipping non-positive lower rate on line " << lineNo << "\n";
      ++skipped;
    }
    if (upperRate > 0.0) {
      upperOut << t9 << "\t" << upperRate << "\n";
      data.xUpper.push_back(t9);
      data.yUpper.push_back(upperRate);
    } else {
      std::cerr << "Skipping non-positive upper rate on line " << lineNo << "\n";
      ++skipped;
    }
    ++written;
  }

  std::cout << "Wrote " << written << " central rows to " << centralPath << "\n";
  std::cout << "Wrote lower/upper rows with non-positive results skipped as needed\n";
  if (skipped > 0) {
    std::cout << "Skipped entries: " << skipped << "\n";
  }
  return data;
}

void draw_rate_variants(const RateVariantData &data) {
  gStyle->SetOptStat(0);

  auto *canvas = new TCanvas("c_rate_variants", "Rate variants", 900, 700);
  canvas->SetLogx();
  canvas->SetLogy();

  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::lowest();
  auto updateRange = [&](const std::vector<double> &values) {
    for (double y : values) {
      if (y > 0.0) {
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
      }
    }
  };
  updateRange(data.yCentral);
  updateRange(data.yLower);
  updateRange(data.yUpper);

  if (!std::isfinite(ymin) || !std::isfinite(ymax) || ymin <= 0.0 || ymax <= ymin) {
    ymin = 1e-30;
    ymax = 1.0;
  }

  const double logYMin = std::log10(ymin);
  const double logYMax = std::log10(ymax);
  const double logYSpan = std::max(1e-6, logYMax - logYMin);
  const double yPlotMin = std::pow(10.0, logYMin - 0.08 * logYSpan);
  const double yPlotMax = std::pow(10.0, logYMax + 0.08 * logYSpan);

  auto *central = new TGraph(static_cast<int>(data.xCentral.size()), data.xCentral.data(),
                             data.yCentral.data());
  central->SetTitle("");
  central->SetLineColor(kBlack);
  central->SetMarkerColor(kBlack);
  central->SetMarkerStyle(20);
  central->SetMarkerSize(0.5);
  central->SetLineWidth(2);
  central->GetXaxis()->SetTitle("T_{9} [GK]");
  central->GetYaxis()->SetTitle("Reaction rate");
  central->GetYaxis()->SetRangeUser(yPlotMin, yPlotMax);
  central->Draw("APL");

  auto *lower = new TGraph(static_cast<int>(data.xLower.size()), data.xLower.data(),
                           data.yLower.data());
  lower->SetLineColor(kBlue + 1);
  lower->SetMarkerColor(kBlue + 1);
  lower->SetMarkerStyle(24);
  lower->SetMarkerSize(0.5);
  lower->SetLineWidth(2);
  lower->Draw("PL SAME");

  auto *upper = new TGraph(static_cast<int>(data.xUpper.size()), data.xUpper.data(),
                           data.yUpper.data());
  upper->SetLineColor(kRed + 1);
  upper->SetMarkerColor(kRed + 1);
  upper->SetMarkerStyle(25);
  upper->SetMarkerSize(0.5);
  upper->SetLineWidth(2);
  upper->Draw("PL SAME");

  auto *legend = new TLegend(0.62, 0.15, 0.90, 0.30);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(central, "central", "lp");
  legend->AddEntry(lower, "lower", "lp");
  legend->AddEntry(upper, "upper", "lp");
  legend->Draw();

  canvas->Update();
}

}  // namespace

void make_rate_variants(const char *inputPath, const char *centralPath = "",
                        const char *lowerPath = "", const char *upperPath = "") {
  const std::string input(inputPath);
  const std::string central =
      std::string(centralPath).empty() ? default_variant_path(input, "central") : centralPath;
  const std::string lower =
      std::string(lowerPath).empty() ? default_variant_path(input, "lower") : lowerPath;
  const std::string upper =
      std::string(upperPath).empty() ? default_variant_path(input, "upper") : upperPath;
  const RateVariantData data =
      write_rate_variants_impl(inputPath, central.c_str(), lower.c_str(), upper.c_str());
  draw_rate_variants(data);
}
