#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <TAxis.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TDecompLU.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TStyle.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct DataPoint {
  double t9 = 0.0;
  double rate = 0.0;
};

struct ParsedTable {
  std::vector<DataPoint> data;
  int skippedNonPositiveRates = 0;
};

std::vector<DataPoint> filterByTmin(const std::vector<DataPoint> &data, double tmin) {
  std::vector<DataPoint> filtered;
  filtered.reserve(data.size());
  for (const auto &p : data) {
    if (p.t9 >= tmin) {
      filtered.push_back(p);
    }
  }
  return filtered;
}

std::string baseName(const std::string &path) {
  const std::string::size_type pos = path.find_last_of("/\\");
  if (pos == std::string::npos) {
    return path;
  }
  return path.substr(pos + 1);
}

std::string trim(const std::string &s) {
  std::string::size_type begin = 0;
  while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin]))) {
    ++begin;
  }
  std::string::size_type end = s.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return s.substr(begin, end - begin);
}

std::string toLower(std::string s) {
  for (char &c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

bool parseBool(const std::string &text) {
  const std::string value = toLower(trim(text));
  if (value == "1" || value == "true" || value == "yes" || value == "on") {
    return true;
  }
  if (value == "0" || value == "false" || value == "no" || value == "off") {
    return false;
  }
  throw std::runtime_error("Invalid boolean value: " + text);
}

std::vector<std::string> splitCommaSeparated(const std::string &text) {
  std::vector<std::string> out;
  std::string current;
  std::istringstream iss(text);
  while (std::getline(iss, current, ',')) {
    out.push_back(trim(current));
  }
  return out;
}

bool hasDirectoryComponent(const std::string &path) {
  return path.find('/') != std::string::npos || path.find('\\') != std::string::npos;
}

bool isAbsolutePath(const std::string &path) {
  return !path.empty() && path[0] == '/';
}

std::string resolveWithDefaultDir(const std::string &path, const std::string &defaultDir) {
  if (path.empty() || isAbsolutePath(path) || hasDirectoryComponent(path)) {
    return path;
  }
  return defaultDir + "/" + path;
}

int parameterIndexFromName(const std::string &name) {
  const std::string low = toLower(trim(name));
  if (low == "a0") return 0;
  if (low == "a1") return 1;
  if (low == "a2") return 2;
  if (low == "a3") return 3;
  if (low == "a4") return 4;
  if (low == "a5") return 5;
  if (low == "a6") return 6;
  throw std::runtime_error("Unknown parameter name: " + name);
}

enum class FitMode {
  ChargedNR,
  NeutronNR,
  NarrowResonance
};

struct FitConfig {
  FitMode mode = FitMode::ChargedNR;
  int orbitalL = 0;
  int nterms = 1;
  bool useEr = false;
  bool fixA1 = false;
  bool fixA2 = false;
  bool fixA3 = false;
  bool fixA4 = false;
  bool fixA5 = false;
  bool fixA6 = false;
  double a1Value = 0.0;
  double a2Value = 0.0;
  double a3Value = 0.0;
  double a4Value = 0.0;
  double a5Value = 0.0;
  double a6Value = 0.0;
};

struct FitResult {
  FitMode mode = FitMode::ChargedNR;
  int nterms = 1;
  double parameters[224] = {0.0};
  double objective = 0.0;
  double rmsLog = 0.0;
  double maxFracErr = 0.0;
  bool converged = false;
  std::vector<std::string> termKinds;
};

struct CompareCurve {
  bool enabled = false;
  int nterms = 1;
  double parameters[224] = {0.0};
  std::vector<std::string> termKinds;
};

enum class TermMode {
  ChargedNR,
  NeutronNR,
  NarrowResonance
};

struct TermSpec {
  TermMode mode = TermMode::ChargedNR;
  double a[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  bool floating[7] = {false, false, false, false, false, false, false};
  std::string label;
};

struct ProgramOptions {
  std::string inputPath;
  std::string termsPath;
  std::string outputPath;
  bool useOutput = false;
  bool noGui = false;
  bool useGlobalReaction = false;
  double A1 = 0.0;
  double A2 = 0.0;
  double Z1 = 0.0;
  double Z2 = 0.0;
  double Threshold = 0.0;
  bool useTmin = false;
  double tmin = 0.0;
  double plotTmin = 0.1;
  double plotTmax = 10.0;
  bool useCompare = false;
  std::string comparePath;
};

struct ReactionGlobals {
  double A1 = 0.0;
  double A2 = 0.0;
  double Z1 = 0.0;
  double Z2 = 0.0;
  double Ared = 0.0;
  double Threshold = 0.0;
};

double reaclibExponent(double t9, const double *par) {
  const double t9m1 = std::pow(t9, -1.0);
  const double t9m13 = std::pow(t9, -1.0 / 3.0);
  const double t913 = std::pow(t9, 1.0 / 3.0);
  const double t9p1 = t9;
  const double t953 = std::pow(t9, 5.0 / 3.0);
  const double lnt9 = std::log(t9);

  return par[0] +
         par[1] * t9m1 +
         par[2] * t9m13 +
         par[3] * t913 +
         par[4] * t9p1 +
         par[5] * t953 +
         par[6] * lnt9;
}

double reaclibRate(double t9, const double *par) {
  const double exponent = reaclibExponent(t9, par);
  if (exponent > 700.0) {
    return std::numeric_limits<double>::infinity();
  }
  return std::exp(exponent);
}

double evaluateReaclibTf1(double *x, double *par) {
  return reaclibRate(x[0], par);
}

double fitRate(double t9, const FitResult &result);

const FitResult *gCurrentFitResult = nullptr;
const CompareCurve *gCurrentCompareCurve = nullptr;

double evaluateFitTf1(double *x, double *) {
  if (gCurrentFitResult == nullptr) {
    return 0.0;
  }
  return fitRate(x[0], *gCurrentFitResult);
}

double evaluateCompareTf1(double *x, double *) {
  if (gCurrentCompareCurve == nullptr) {
    return 0.0;
  }
  double sum = 0.0;
  for (int term = 0; term < gCurrentCompareCurve->nterms; ++term) {
    const double rate = reaclibRate(x[0], gCurrentCompareCurve->parameters + 7 * term);
    if (!std::isfinite(rate)) {
      return std::numeric_limits<double>::infinity();
    }
    sum += rate;
  }
  return sum;
}

double fitRate(double t9, const FitResult &result) {
  double sum = 0.0;
  for (int term = 0; term < result.nterms; ++term) {
    const double *par = result.parameters + 7 * term;
    const double rate = reaclibRate(t9, par);
    if (!std::isfinite(rate)) {
      return std::numeric_limits<double>::infinity();
    }
    sum += rate;
  }
  return sum;
}

std::string modeName(FitMode mode) {
  switch (mode) {
    case FitMode::ChargedNR:
      return "charged_nr";
    case FitMode::NeutronNR:
      return "neutron_nr";
    case FitMode::NarrowResonance:
      return "narrow_resonance";
  }
  return "unknown";
}

FitMode parseMode(const std::string &text) {
  if (text == "charged_nr") {
    return FitMode::ChargedNR;
  }
  if (text == "neutron_nr") {
    return FitMode::NeutronNR;
  }
  if (text == "narrow_resonance") {
    return FitMode::NarrowResonance;
  }
  throw std::runtime_error("Unknown mode: " + text);
}

ParsedTable loadTable(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open input file: " + path);
  }

  ParsedTable table;
  std::string line;
  int lineNo = 0;

  while (std::getline(in, line)) {
    ++lineNo;
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    DataPoint p;
    if (!(iss >> p.t9 >> p.rate)) {
      throw std::runtime_error("Failed to parse line " + std::to_string(lineNo));
    }
    if (p.t9 <= 0.0) {
      throw std::runtime_error("Temperature must be > 0 on line " + std::to_string(lineNo));
    }
    if (p.rate <= 0.0) {
      ++table.skippedNonPositiveRates;
      continue;
    }

    table.data.push_back(p);
  }

  if (table.data.empty()) {
    throw std::runtime_error("No positive-rate data rows were found in: " + path);
  }
  return table;
}

CompareCurve loadCompareCurve(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open compare parameter file: " + path);
  }

  CompareCurve curve;
  curve.enabled = true;
  std::string line;
  int lineNo = 0;
  int termCount = 0;

  auto isNumericToken = [](const std::string &text) {
    char *end = nullptr;
    std::strtod(text.c_str(), &end);
    return end != text.c_str() && *end == '\0';
  };

  while (std::getline(in, line)) {
    ++lineNo;
    const std::string stripped = trim(line);
    if (stripped.empty() || stripped[0] == '#') {
      continue;
    }

    std::istringstream iss(stripped);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }
    if (tokens.empty()) {
      continue;
    }

    if (tokens.size() < 8) {
      throw std::runtime_error("Compare parameter file parse error on line " +
                               std::to_string(lineNo));
    }

    if ((termCount + 1) * 7 > 224) {
      throw std::runtime_error("Compare parameter file contains too many terms");
    }

    for (int i = 0; i < 7; ++i) {
      if (!isNumericToken(tokens[i])) {
        throw std::runtime_error("Expected numeric REACLIB parameter on line " +
                                 std::to_string(lineNo));
      }
      curve.parameters[7 * termCount + i] = std::stod(tokens[i]);
    }

    std::string kind;
    kind = toLower(tokens[7]);
    if (kind != "n" && kind != "r") {
      throw std::runtime_error("Unknown compare term type on line " +
                               std::to_string(lineNo) + ": " + tokens[7]);
    }
    curve.termKinds.push_back(kind);
    ++termCount;
  }

  if (termCount == 0) {
    throw std::runtime_error("No valid compare terms found in: " + path);
  }

  curve.nterms = termCount;
  return curve;
}

constexpr double kReaclibB = 7.8318e9;
constexpr double kReaclibC = 0.08617;
constexpr double kReaclibD = 1.5394e11;
constexpr double kAvogadro = 6.02214076e23;

TermSpec parseTermLine(const std::string &line, int lineNo, const ReactionGlobals &globals) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }
  if (tokens.empty()) {
    throw std::runtime_error("Empty term specification on line " + std::to_string(lineNo));
  }

  TermSpec term;
  const std::string mode = toLower(tokens[0]);
  std::vector<int> extraFloat;

  auto parseFloatToken = [&](const std::string &tok) {
    const std::string::size_type eq = tok.find('=');
    if (eq == std::string::npos) {
      throw std::runtime_error("Invalid float spec on line " + std::to_string(lineNo));
    }
    for (const auto &name : splitCommaSeparated(tok.substr(eq + 1))) {
      extraFloat.push_back(parameterIndexFromName(name));
    }
  };

  if (mode == "nr") {
    if (tokens.size() < 2) {
      throw std::runtime_error("NR line must be: NR S0 [float=...]");
    }
    if (globals.Ared <= 0.0) {
      throw std::runtime_error("NR term requires global A1/A2/Z1/Z2 to be set");
    }
    const double s0 = std::stod(tokens[1]);
    term.mode = TermMode::ChargedNR;
    term.label = "NR";
    term.a[0] =
        std::log(kReaclibB * std::pow((globals.Z1 * globals.Z2) / globals.Ared, 1.0 / 3.0) * s0);
    term.a[1] = 0.0;
    term.a[2] = -4.2486 *
                std::pow(globals.Z1 * globals.Z1 * globals.Z2 * globals.Z2 * globals.Ared,
                         1.0 / 3.0);
    term.a[6] = -2.0 / 3.0;
    term.floating[3] = true;
    term.floating[4] = true;
    term.floating[5] = true;
    for (size_t i = 2; i < tokens.size(); ++i) {
      parseFloatToken(tokens[i]);
    }
  } else if (mode == "n") {
    if (tokens.size() < 3) {
      throw std::runtime_error(
          "n line must be: n l sig0 [float=...]  where sig0 = (sigma*v/E^l) at E=0");
    }
    const int ell = std::stoi(tokens[1]);
    const double sig0 = std::stod(tokens[2]);
    term.mode = TermMode::NeutronNR;
    term.label = "n";
    term.a[0] = std::log(kAvogadro * std::pow(kReaclibC, ell) * std::tgamma(ell + 1.5) /
                         std::tgamma(1.5) * sig0);
    term.a[1] = 0.0;
    term.a[2] = 0.0;
    term.a[6] = static_cast<double>(ell);
    term.floating[3] = true;
    term.floating[4] = true;
    term.floating[5] = true;
    for (size_t i = 3; i < tokens.size(); ++i) {
      parseFloatToken(tokens[i]);
    }
  } else if (mode == "r") {
    if (tokens.size() < 3) {
      throw std::runtime_error("R line must be: R omega_gamma Ex [float=...]");
    }
    if (globals.Ared <= 0.0) {
      throw std::runtime_error("R term requires global A1/A2/Z1/Z2 to be set");
    }
    if (!std::isfinite(globals.Threshold)) {
      throw std::runtime_error("R term requires global Threshold to be set");
    }
    const double omegaGamma = std::stod(tokens[1]);
    const double ex = std::stod(tokens[2]);
    const double er = ex - globals.Threshold;
    term.mode = TermMode::NarrowResonance;
    term.label = "R";
    term.a[0] = std::log(kReaclibD * std::pow(globals.Ared, -1.5) * omegaGamma);
    term.a[1] = -11.6045 * er;
    term.a[2] = 0.0;
    term.a[3] = 0.0;
    term.a[4] = 0.0;
    term.a[5] = 0.0;
    term.a[6] = -1.5;
    for (size_t i = 3; i < tokens.size(); ++i) {
      parseFloatToken(tokens[i]);
    }
  } else {
    throw std::runtime_error("Unknown term mode on line " + std::to_string(lineNo) + ": " +
                             tokens[0]);
  }

  for (int idx : extraFloat) {
    term.floating[idx] = true;
  }
  return term;
}

std::vector<TermSpec> loadTermFile(const std::string &path, const ReactionGlobals &globals) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open term file: " + path);
  }

  std::vector<TermSpec> terms;
  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string stripped = trim(line);
    if (stripped.empty() || stripped[0] == '#') {
      continue;
    }
    terms.push_back(parseTermLine(stripped, lineNo, globals));
  }
  if (terms.empty()) {
    throw std::runtime_error("No valid terms found in: " + path);
  }
  return terms;
}

void applyConfigEntry(ProgramOptions &options, const std::string &key, const std::string &value) {
  const std::string rawKey = trim(key);
  const std::string k = toLower(rawKey);
  const std::string v = trim(value);

  if (k == "input") {
    options.inputPath = v;
  } else if (k == "terms") {
    options.termsPath = v;
  } else if (k == "output") {
    options.outputPath = v;
    options.useOutput = !v.empty();
  } else if (rawKey == "A1" || k == "a1_global" || k == "mass1") {
    options.A1 = std::stod(v);
    options.useGlobalReaction = true;
  } else if (rawKey == "A2" || k == "a2_global" || k == "mass2") {
    options.A2 = std::stod(v);
    options.useGlobalReaction = true;
  } else if (rawKey == "Z1" || k == "z1") {
    options.Z1 = std::stod(v);
    options.useGlobalReaction = true;
  } else if (rawKey == "Z2" || k == "z2") {
    options.Z2 = std::stod(v);
    options.useGlobalReaction = true;
  } else if (rawKey == "Threshold" || k == "threshold" || k == "channel_threshold") {
    options.Threshold = std::stod(v);
    options.useGlobalReaction = true;
  } else if (k == "tmin") {
    options.tmin = std::stod(v);
    options.useTmin = true;
  } else if (k == "plot_tmin" || k == "xmin") {
    options.plotTmin = std::stod(v);
  } else if (k == "plot_tmax" || k == "xmax") {
    options.plotTmax = std::stod(v);
  } else if (k == "compare") {
    options.comparePath = v;
    options.useCompare = !v.empty();
  } else if (k == "use_compare") {
    options.useCompare = parseBool(v);
  } else {
    throw std::runtime_error("Unknown config key: " + key);
  }
}

ProgramOptions loadConfigFile(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open config file: " + path);
  }

  ProgramOptions options;
  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::string stripped = trim(line);
    if (stripped.empty() || stripped[0] == '#') {
      continue;
    }
    const std::string::size_type eq = stripped.find('=');
    if (eq == std::string::npos) {
      throw std::runtime_error("Config parse error on line " + std::to_string(lineNo));
    }
    applyConfigEntry(options, stripped.substr(0, eq), stripped.substr(eq + 1));
  }
  return options;
}

FitConfig buildConfig(FitMode mode, int orbitalL, bool useA2, double a2Value, bool useEr,
                      double erValue) {
  FitConfig cfg;
  cfg.mode = mode;
  cfg.orbitalL = orbitalL;
  cfg.useEr = useEr;

  switch (mode) {
    case FitMode::ChargedNR:
      cfg.fixA1 = true;
      cfg.a1Value = 0.0;
      cfg.fixA6 = true;
      cfg.a6Value = -2.0 / 3.0;
      if (useA2) {
        cfg.fixA2 = true;
        cfg.a2Value = a2Value;
      }
      break;
    case FitMode::NeutronNR:
      cfg.fixA1 = true;
      cfg.a1Value = 0.0;
      cfg.fixA2 = true;
      cfg.a2Value = 0.0;
      cfg.fixA6 = true;
      cfg.a6Value = static_cast<double>(orbitalL);
      break;
    case FitMode::NarrowResonance:
      if (useEr) {
        cfg.fixA1 = true;
        cfg.a1Value = -11.6045 * erValue;
      }
      cfg.fixA2 = true;
      cfg.a2Value = 0.0;
      cfg.fixA3 = true;
      cfg.a3Value = 0.0;
      cfg.fixA4 = true;
      cfg.a4Value = 0.0;
      cfg.fixA5 = true;
      cfg.a5Value = 0.0;
      cfg.fixA6 = true;
      cfg.a6Value = -1.5;
      break;
  }

  return cfg;
}

void seedInitialParameters(const std::vector<DataPoint> &data, const FitConfig &cfg,
                           double *par) {
  for (int i = 0; i < 7; ++i) {
    par[i] = 0.0;
  }

  const DataPoint &mid = data[data.size() / 2];
  par[0] = std::log(mid.rate);

  if (cfg.fixA1) {
    par[1] = cfg.a1Value;
  }
  if (cfg.fixA2) {
    par[2] = cfg.a2Value;
  }
  if (cfg.fixA3) {
    par[3] = cfg.a3Value;
  }
  if (cfg.fixA4) {
    par[4] = cfg.a4Value;
  }
  if (cfg.fixA5) {
    par[5] = cfg.a5Value;
  }
  if (cfg.fixA6) {
    par[6] = cfg.a6Value;
  }
}

class Chi2 {
 public:
  Chi2(std::vector<DataPoint> data, FitConfig config)
      : data_(std::move(data)), config_(config) {}

  double operator()(const double *x) const {
    double par[7];
    for (int i = 0; i < 7; ++i) {
      par[i] = x[i];
    }

    if (config_.fixA1) {
      par[1] = config_.a1Value;
    }
    if (config_.fixA2) {
      par[2] = config_.a2Value;
    }
    if (config_.fixA3) {
      par[3] = config_.a3Value;
    }
    if (config_.fixA4) {
      par[4] = config_.a4Value;
    }
    if (config_.fixA5) {
      par[5] = config_.a5Value;
    }
    if (config_.fixA6) {
      par[6] = config_.a6Value;
    }

    double chi2 = 0.0;
    for (const auto &p : data_) {
      const double model = reaclibRate(p.t9, par);
      const double dlog = std::log(model) - std::log(p.rate);
      chi2 += dlog * dlog;
    }
    return chi2;
  }

 private:
  std::vector<DataPoint> data_;
  FitConfig config_;
};

void printUsage(const char *prog) {
  std::cerr
      << "Usage: " << prog << " input.txt --terms terms.txt [options]\n"
      << "   or: " << prog << " --input input.txt --terms terms.txt [options]\n"
      << "   or: " << prog << " -c exp.conf [options]\n"
      << "Options:\n"
      << "  --config, -c <file>   Read config file (default dir: conf/)\n"
      << "  --input <file>        Read input data file (default dir: input/)\n"
      << "  --terms <file>        Read term file (default dir: conf/)\n"
      << "  --output <file>       Write fit summary (default dir: output/)\n"
      << "  --no-gui              Run fit without opening the ROOT canvas\n"
      << "  --tmin <double>       Fit only data with T9 >= tmin\n"
      << "  --compare <file>      Overlay reference curve (default dir: input/)\n"
      << "  config only: plot_tmin, plot_tmax (default 0.1, 10)\n"
      << "\n"
      << "Input file format:\n"
      << "  T9[GK]   rate\n"
      << "  # comments are allowed\n"
      << "\nTerm file format examples:\n"
      << "  NR S0 [float=a0,a1,a2,a6]\n"
      << "  n l sig0 [float=a0,a1,a2,a6]\n"
      << "  R omega_gamma Ex [float=a0,a1,a2,a3,a4,a5,a6]\n"
      << "\nCompare file format examples:\n"
      << "  a0 a1 a2 a3 a4 a5 a6 n\n"
      << "  a0 a1 a2 a3 a4 a5 a6 r\n";
}

double maxFractionalError(const std::vector<DataPoint> &data, const double *par) {
  double maxErr = 0.0;
  for (const auto &p : data) {
    const double model = reaclibRate(p.t9, par);
    if (!std::isfinite(model)) {
      return std::numeric_limits<double>::infinity();
    }
    const double frac = std::abs(model / p.rate - 1.0);
    if (frac > maxErr) {
      maxErr = frac;
    }
  }
  return maxErr;
}

double rmsLogError(const std::vector<DataPoint> &data, const double *par) {
  double sum = 0.0;
  for (const auto &p : data) {
    const double model = reaclibRate(p.t9, par);
    if (!std::isfinite(model) || model <= 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    const double dlog = std::log(model) - std::log(p.rate);
    sum += dlog * dlog;
  }
  return std::sqrt(sum / static_cast<double>(data.size()));
}

double rateFromPackedTerms(double t9, const std::vector<double> &packed) {
  double sum = 0.0;
  for (size_t i = 0; i < packed.size(); i += 7) {
    const double rate = reaclibRate(t9, packed.data() + static_cast<long>(i));
    if (!std::isfinite(rate)) {
      return std::numeric_limits<double>::infinity();
    }
    sum += rate;
  }
  return sum;
}

class MultiTermChi2 {
 public:
  MultiTermChi2(std::vector<DataPoint> data, std::vector<TermSpec> terms, FitConfig config)
      : data_(std::move(data)), terms_(std::move(terms)), config_(config) {
    for (size_t it = 0; it < terms_.size(); ++it) {
      for (int ip = 0; ip < 7; ++ip) {
        if (terms_[it].floating[ip]) {
          mapping_.push_back({static_cast<int>(it), ip});
        }
      }
    }
  }

  std::vector<double> apply(const double *x) const {
    std::vector<double> packed(7 * terms_.size(), 0.0);
    for (size_t it = 0; it < terms_.size(); ++it) {
      for (int ip = 0; ip < 7; ++ip) {
        packed[7 * it + ip] = terms_[it].a[ip];
      }
    }
    for (size_t i = 0; i < mapping_.size(); ++i) {
      packed[7 * mapping_[i].first + mapping_[i].second] = x[i];
    }
    return packed;
  }

  double operator()(const double *x) const {
    const std::vector<double> packed = apply(x);
    double chi2 = 0.0;
    for (const auto &p : data_) {
      const double model = rateFromPackedTerms(p.t9, packed);
      if (!std::isfinite(model) || model <= 0.0) {
        return std::numeric_limits<double>::infinity();
      }
      const double dlog = std::log(model) - std::log(p.rate);
      chi2 += dlog * dlog;
    }
    return chi2;
  }

  size_t nFloating() const { return mapping_.size(); }

  std::string name(size_t i) const {
    static const char *anames[7] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
    return std::string(anames[mapping_[i].second]) + "_" + std::to_string(mapping_[i].first + 1);
  }

  double initial(size_t i) const { return terms_[mapping_[i].first].a[mapping_[i].second]; }
  int parameterIndex(size_t i) const { return mapping_[i].second; }

 private:
  std::vector<DataPoint> data_;
  std::vector<TermSpec> terms_;
  FitConfig config_;
  std::vector<std::pair<int, int>> mapping_;
};

double fallbackInitialValue(const std::vector<DataPoint> &data, int parameterIndex) {
  if (parameterIndex == 0) {
    double maxLogRate = std::numeric_limits<double>::lowest();
    for (const auto &p : data) {
      if (p.rate > 0.0) {
        maxLogRate = std::max(maxLogRate, std::log(p.rate));
      }
    }
    if (std::isfinite(maxLogRate)) {
      return maxLogRate;
    }
  }
  return 0.0;
}

FitResult runTermSpecsFit(const std::vector<DataPoint> &data, const std::vector<TermSpec> &terms,
                          const FitConfig &cfg) {
  MultiTermChi2 chi2(data, terms, cfg);
  FitResult result;
  result.mode = FitMode::ChargedNR;
  result.nterms = static_cast<int>(terms.size());
  result.termKinds.reserve(terms.size());
  for (const auto &term : terms) {
    if (term.mode == TermMode::NarrowResonance) {
      result.termKinds.push_back("r");
    } else {
      result.termKinds.push_back("n");
    }
  }

  if (chi2.nFloating() == 0) {
    for (size_t it = 0; it < terms.size(); ++it) {
      for (int ip = 0; ip < 7; ++ip) {
        result.parameters[7 * it + ip] = terms[it].a[ip];
      }
    }
    result.converged = true;
  } else {
    auto minimizer =
        std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer(
            "Minuit2", "Migrad"));
    if (!minimizer) {
      throw std::runtime_error("Failed to create ROOT minimizer");
    }
    minimizer->SetMaxFunctionCalls(200000);
    minimizer->SetMaxIterations(200000);
    minimizer->SetTolerance(1e-8);
    minimizer->SetPrintLevel(0);

    ROOT::Math::Functor f(chi2, static_cast<unsigned int>(chi2.nFloating()));
    minimizer->SetFunction(f);
    for (size_t i = 0; i < chi2.nFloating(); ++i) {
      double init = chi2.initial(i);
      if (!std::isfinite(init)) {
        init = fallbackInitialValue(data, chi2.parameterIndex(i));
      }
      minimizer->SetVariable(static_cast<unsigned int>(i), chi2.name(i).c_str(), init, 0.1);
    }
    result.converged = minimizer->Minimize();
    const std::vector<double> packed = chi2.apply(minimizer->X());
    for (size_t i = 0; i < packed.size() && i < 224; ++i) {
      result.parameters[i] = packed[i];
    }
    result.objective = minimizer->MinValue();
  }

  if (chi2.nFloating() == 0) {
    std::vector<double> packed(7 * terms.size(), 0.0);
    for (size_t it = 0; it < terms.size(); ++it) {
      for (int ip = 0; ip < 7; ++ip) {
        packed[7 * it + ip] = terms[it].a[ip];
      }
    }
    double objective = 0.0;
    for (const auto &p : data) {
      const double model = rateFromPackedTerms(p.t9, packed);
      const double dlog = std::log(model) - std::log(p.rate);
      objective += dlog * dlog;
    }
    result.objective = objective;
  }

  double sum = 0.0;
  double maxErr = 0.0;
  std::vector<double> packed(7 * terms.size(), 0.0);
  for (size_t it = 0; it < terms.size(); ++it) {
    for (int ip = 0; ip < 7; ++ip) {
      packed[7 * it + ip] = result.parameters[7 * it + ip];
    }
  }
  for (const auto &p : data) {
    const double model = rateFromPackedTerms(p.t9, packed);
    if (!std::isfinite(model) || model <= 0.0) {
      result.rmsLog = std::numeric_limits<double>::infinity();
      result.maxFracErr = std::numeric_limits<double>::infinity();
      return result;
    }
    const double dlog = std::log(model) - std::log(p.rate);
    sum += dlog * dlog;
    maxErr = std::max(maxErr, std::abs(model / p.rate - 1.0));
  }
  result.rmsLog = std::sqrt(sum / static_cast<double>(data.size()));
  result.maxFracErr = maxErr;
  return result;
}

void writeFitSummary(const std::string &path, const ProgramOptions &options, const FitResult &result,
                     const CompareCurve &compare, int inputPoints, int skippedPoints,
                     int fitPoints) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Failed to open output file: " + path);
  }

  out << std::setprecision(12);
  out << "# Input points: " << inputPoints << "\n";
  if (skippedPoints > 0) {
    out << "# Skipped non-positive rates: " << skippedPoints << "\n";
  }
  if (options.useTmin) {
    out << "# Fit points with T9 >= " << options.tmin << ": " << fitPoints << "\n";
  }
  out << "# Summed terms: " << result.nterms << "\n";
  out << "# Objective (sum of log-residual^2): " << result.objective << "\n";
  out << "# RMS log error: " << result.rmsLog << "\n";
  out << "# Max fractional error: " << std::fixed << std::setprecision(2)
      << 100.0 * result.maxFracErr << std::defaultfloat << " %\n";

  out << "# Fitted parameters\n";
  for (int term = 0; term < result.nterms; ++term) {
    for (int i = 0; i < 7; ++i) {
      if (i > 0) {
        out << "\t";
      }
      out << std::scientific << std::setprecision(6) << result.parameters[7 * term + i];
    }
    if (term < static_cast<int>(result.termKinds.size()) && !result.termKinds[term].empty()) {
      out << "\t" << result.termKinds[term];
    }
    out << std::defaultfloat << "\n";
  }

  if (compare.enabled) {
    out << "# Reference parameters\n";
    for (int term = 0; term < compare.nterms; ++term) {
      for (int i = 0; i < 7; ++i) {
        if (i > 0) {
          out << "\t";
        }
        out << std::scientific << std::setprecision(6) << compare.parameters[7 * term + i];
      }
      if (term < static_cast<int>(compare.termKinds.size()) && !compare.termKinds[term].empty()) {
        out << "\t" << compare.termKinds[term];
      }
      out << std::defaultfloat << "\n";
    }
  }
}

FitResult runChargedNrLinearFit(const std::vector<DataPoint> &data, const FitConfig &cfg) {
  TMatrixD normal(4, 4);
  TVectorD rhs(4);

  for (const auto &p : data) {
    const double basis[4] = {
        1.0,
        std::pow(p.t9, 1.0 / 3.0),
        p.t9,
        std::pow(p.t9, 5.0 / 3.0),
    };
    const double target = std::log(p.rate) - cfg.a1Value * std::pow(p.t9, -1.0) -
                          cfg.a2Value * std::pow(p.t9, -1.0 / 3.0) -
                          cfg.a6Value * std::log(p.t9);

    for (int i = 0; i < 4; ++i) {
      rhs[i] += basis[i] * target;
      for (int j = 0; j < 4; ++j) {
        normal(i, j) += basis[i] * basis[j];
      }
    }
  }

  TDecompLU decomp(normal);
  Bool_t ok = false;
  TVectorD solution = decomp.Solve(rhs, ok);
  if (!ok) {
    throw std::runtime_error("Linear least-squares solve failed for charged_nr");
  }

  FitResult result;
  result.mode = cfg.mode;
  result.nterms = 1;
  result.converged = true;
  result.parameters[0] = solution[0];
  result.parameters[1] = cfg.a1Value;
  result.parameters[2] = cfg.a2Value;
  result.parameters[3] = solution[1];
  result.parameters[4] = solution[2];
  result.parameters[5] = solution[3];
  result.parameters[6] = cfg.a6Value;

  double objective = 0.0;
  for (const auto &p : data) {
    const double dlog = std::log(reaclibRate(p.t9, result.parameters)) - std::log(p.rate);
    objective += dlog * dlog;
  }
  result.objective = objective;
  result.rmsLog = rmsLogError(data, result.parameters);
  result.maxFracErr = maxFractionalError(data, result.parameters);
  return result;
}

class ChargedNrTwoTermChi2 {
 public:
  ChargedNrTwoTermChi2(std::vector<DataPoint> data, FitConfig config)
      : data_(std::move(data)), config_(config) {}

  double operator()(const double *x) const {
    double par1[7] = {x[0], config_.a1Value, config_.a2Value, x[1], x[2], x[3], config_.a6Value};
    double par2[7] = {x[4], config_.a1Value, config_.a2Value, x[5], x[6], x[7], config_.a6Value};

    double chi2 = 0.0;
    for (const auto &p : data_) {
      const double model = reaclibRate(p.t9, par1) + reaclibRate(p.t9, par2);
      if (!std::isfinite(model) || model <= 0.0) {
        return std::numeric_limits<double>::infinity();
      }
      const double dlog = std::log(model) - std::log(p.rate);
      chi2 += dlog * dlog;
    }
    return chi2;
  }

 private:
  std::vector<DataPoint> data_;
  FitConfig config_;
};

FitResult runChargedNrTwoTermFit(const std::vector<DataPoint> &data, const FitConfig &cfg) {
  auto minimizer =
      std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer(
          "Minuit2", "Migrad"));
  if (!minimizer) {
    throw std::runtime_error("Failed to create ROOT minimizer");
  }

  minimizer->SetMaxFunctionCalls(200000);
  minimizer->SetMaxIterations(200000);
  minimizer->SetTolerance(1e-8);
  minimizer->SetPrintLevel(0);

  ChargedNrTwoTermChi2 chi2(data, cfg);
  ROOT::Math::Functor f(chi2, 8);
  minimizer->SetFunction(f);

  const DataPoint &mid = data[data.size() / 2];
  const DataPoint &last = data.back();
  const double init[8] = {
      std::log(std::max(mid.rate * 0.6, 1e-300)), 0.0, 0.0, 0.0,
      std::log(std::max(last.rate * 0.4, 1e-300)), 0.0, 0.0, 0.0};
  const char *names[8] = {"a0_1", "a3_1", "a4_1", "a5_1", "a0_2", "a3_2", "a4_2", "a5_2"};
  for (int i = 0; i < 8; ++i) {
    minimizer->SetVariable(i, names[i], init[i], 0.1);
  }

  FitResult result;
  result.mode = cfg.mode;
  result.nterms = 2;
  result.converged = minimizer->Minimize();
  result.objective = minimizer->MinValue();

  const double *xs = minimizer->X();
  double *p1 = result.parameters;
  double *p2 = result.parameters + 7;
  p1[0] = xs[0];
  p1[1] = cfg.a1Value;
  p1[2] = cfg.a2Value;
  p1[3] = xs[1];
  p1[4] = xs[2];
  p1[5] = xs[3];
  p1[6] = cfg.a6Value;
  p2[0] = xs[4];
  p2[1] = cfg.a1Value;
  p2[2] = cfg.a2Value;
  p2[3] = xs[5];
  p2[4] = xs[6];
  p2[5] = xs[7];
  p2[6] = cfg.a6Value;

  double sum = 0.0;
  double maxErr = 0.0;
  for (const auto &p : data) {
    const double model = fitRate(p.t9, result);
    if (!std::isfinite(model) || model <= 0.0) {
      result.rmsLog = std::numeric_limits<double>::infinity();
      result.maxFracErr = std::numeric_limits<double>::infinity();
      return result;
    }
    const double dlog = std::log(model) - std::log(p.rate);
    sum += dlog * dlog;
    maxErr = std::max(maxErr, std::abs(model / p.rate - 1.0));
  }
  result.rmsLog = std::sqrt(sum / static_cast<double>(data.size()));
  result.maxFracErr = maxErr;
  return result;
}

FitResult runFit(const std::vector<DataPoint> &data, const FitConfig &cfg) {
  if (cfg.mode == FitMode::ChargedNR && cfg.nterms == 2 && cfg.fixA1 && cfg.fixA2 && cfg.fixA6) {
    return runChargedNrTwoTermFit(data, cfg);
  }
  if (cfg.mode == FitMode::ChargedNR && cfg.fixA1 && cfg.fixA2 && cfg.fixA6 &&
      !cfg.fixA3 && !cfg.fixA4 && !cfg.fixA5) {
    return runChargedNrLinearFit(data, cfg);
  }

  double initial[7];
  seedInitialParameters(data, cfg, initial);

  auto minimizer =
      std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer(
          "Minuit2", "Migrad"));
  if (!minimizer) {
    throw std::runtime_error("Failed to create ROOT minimizer");
  }

  minimizer->SetMaxFunctionCalls(100000);
  minimizer->SetMaxIterations(100000);
  minimizer->SetTolerance(1e-8);
  minimizer->SetPrintLevel(0);

  Chi2 chi2(data, cfg);
  ROOT::Math::Functor f(chi2, 7);
  minimizer->SetFunction(f);

  const char *names[7] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
  const double steps[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.01};

  for (int i = 0; i < 7; ++i) {
    minimizer->SetVariable(i, names[i], initial[i], steps[i]);
  }

  if (cfg.fixA1) {
    minimizer->FixVariable(1);
    minimizer->SetVariableValue(1, cfg.a1Value);
  }
  if (cfg.fixA2) {
    minimizer->FixVariable(2);
    minimizer->SetVariableValue(2, cfg.a2Value);
  }
  if (cfg.fixA3) {
    minimizer->FixVariable(3);
    minimizer->SetVariableValue(3, cfg.a3Value);
  }
  if (cfg.fixA4) {
    minimizer->FixVariable(4);
    minimizer->SetVariableValue(4, cfg.a4Value);
  }
  if (cfg.fixA5) {
    minimizer->FixVariable(5);
    minimizer->SetVariableValue(5, cfg.a5Value);
  }
  if (cfg.fixA6) {
    minimizer->FixVariable(6);
    minimizer->SetVariableValue(6, cfg.a6Value);
  }

  FitResult result;
  result.mode = cfg.mode;
  result.nterms = 1;
  result.converged = minimizer->Minimize();
  result.objective = minimizer->MinValue();

  const double *xs = minimizer->X();
  for (int i = 0; i < 7; ++i) {
    result.parameters[i] = xs[i];
  }
  if (cfg.fixA1) {
    result.parameters[1] = cfg.a1Value;
  }
  if (cfg.fixA2) {
    result.parameters[2] = cfg.a2Value;
  }
  if (cfg.fixA3) {
    result.parameters[3] = cfg.a3Value;
  }
  if (cfg.fixA4) {
    result.parameters[4] = cfg.a4Value;
  }
  if (cfg.fixA5) {
    result.parameters[5] = cfg.a5Value;
  }
  if (cfg.fixA6) {
    result.parameters[6] = cfg.a6Value;
  }

  result.rmsLog = rmsLogError(data, result.parameters);
  result.maxFracErr = maxFractionalError(data, result.parameters);
  return result;
}

TCanvas *makePlot(const std::vector<DataPoint> &data, const FitResult &result,
                  const CompareCurve &compare, const std::string &inputLabel,
                  const std::string &compareLabel, double plotXmin, double plotXmax) {
  gStyle->SetOptStat(0);

  std::vector<double> x;
  std::vector<double> yInput;
  x.reserve(data.size());
  yInput.reserve(data.size());

  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::lowest();
  double ratioMaxDeviation = 0.0;

  for (const auto &p : data) {
    x.push_back(p.t9);
    yInput.push_back(p.rate);

    ymin = std::min(ymin, p.rate);
    ymax = std::max(ymax, p.rate);
  }

  const double xmin = plotXmin;
  const double xmax = plotXmax;

  for (const auto &p : data) {
    const double fit = fitRate(p.t9, result);
    if (!std::isfinite(fit) || fit <= 0.0) {
      continue;
    }
    const double ratio = fit / p.rate;
    ratioMaxDeviation = std::max(ratioMaxDeviation, std::abs(ratio - 1.0));
  }
  if (compare.enabled) {
    for (const auto &p : data) {
      double ref = 0.0;
      for (int term = 0; term < compare.nterms; ++term) {
        ref += reaclibRate(p.t9, compare.parameters + 7 * term);
      }
      if (!std::isfinite(ref) || ref <= 0.0) {
        continue;
      }
      const double ratio = ref / p.rate;
      ratioMaxDeviation = std::max(ratioMaxDeviation, std::abs(ratio - 1.0));
    }
  }

  auto *canvas = new TCanvas("c_reaclib", "REACLIB fit", 900, 800);
  auto *topPad = new TPad("pad_top", "pad_top", 0.0, 0.30, 1.0, 1.0);
  auto *bottomPad = new TPad("pad_bottom", "pad_bottom", 0.0, 0.0, 1.0, 0.30);

  topPad->SetLeftMargin(0.13);
  topPad->SetRightMargin(0.05);
  topPad->SetTopMargin(0.05);
  topPad->SetBottomMargin(0.02);
  topPad->SetLogx();
  topPad->SetLogy();

  bottomPad->SetLeftMargin(0.13);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetTopMargin(0.03);
  bottomPad->SetBottomMargin(0.30);
  bottomPad->SetLogx();

  topPad->Draw();
  bottomPad->Draw();

  topPad->cd();
  const double logYMin = std::log10(ymin);
  const double logYMax = std::log10(ymax);
  const double logYSpan = std::max(1e-6, logYMax - logYMin);
  const double yPlotMin = std::pow(10.0, logYMin - 0.08 * logYSpan);
  const double yPlotMax = std::pow(10.0, logYMax + 0.08 * logYSpan);
  auto *inputGraph = new TGraph(static_cast<int>(x.size()), x.data(), yInput.data());
  inputGraph->SetTitle("");
  inputGraph->SetMarkerStyle(20);
  inputGraph->SetMarkerSize(0.8);
  inputGraph->SetLineWidth(2);
  inputGraph->GetXaxis()->SetLimits(xmin, xmax);
  inputGraph->GetYaxis()->SetRangeUser(yPlotMin, yPlotMax);
  inputGraph->GetYaxis()->SetTitle("Reaction rate");
  inputGraph->GetYaxis()->SetTitleSize(0.05);
  inputGraph->GetYaxis()->SetTitleOffset(1.1);
  inputGraph->GetYaxis()->SetLabelSize(0.04);
  inputGraph->GetXaxis()->SetLabelSize(0.0);
  inputGraph->Draw("APL");

  auto *legend = new TLegend(0.58, 0.12, 0.92, 0.32);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(inputGraph, inputLabel.c_str(), "lp");
  gCurrentFitResult = &result;
  auto *fitFunc = new TF1("fit_func", evaluateFitTf1, xmin, xmax, 0);
  fitFunc->SetNpx(1000);
  fitFunc->SetLineColor(kRed + 1);
  fitFunc->SetLineWidth(3);
  fitFunc->Draw("SAME");
  legend->AddEntry(fitFunc, "this fit", "l");
  if (compare.enabled) {
    gCurrentCompareCurve = &compare;
    auto *refFunc = new TF1("ref_func", evaluateCompareTf1, xmin, xmax, 0);
    refFunc->SetNpx(1000);
    refFunc->SetLineColor(kBlue + 1);
    refFunc->SetLineWidth(3);
    refFunc->SetLineStyle(2);
    refFunc->Draw("SAME");
    legend->AddEntry(refFunc, compareLabel.c_str(), "l");
  }
  legend->Draw();

  bottomPad->cd();
  std::vector<double> xRatio;
  std::vector<double> yRatio;
  xRatio.reserve(data.size());
  yRatio.reserve(data.size());
  for (const auto &p : data) {
    const double fit = fitRate(p.t9, result);
    if (std::isfinite(fit) && fit > 0.0) {
      xRatio.push_back(p.t9);
      yRatio.push_back(fit / p.rate);
    }
  }

  auto *ratioGraph = new TGraph(static_cast<int>(xRatio.size()), xRatio.data(), yRatio.data());
  ratioGraph->SetTitle("");
  ratioGraph->SetMarkerStyle(20);
  ratioGraph->SetMarkerSize(0.8);
  ratioGraph->SetLineWidth(2);
  ratioGraph->SetLineColor(kRed + 1);
  ratioGraph->SetMarkerColor(kRed + 1);
  ratioGraph->GetXaxis()->SetLimits(xmin, xmax);
  ratioGraph->GetYaxis()->SetRangeUser(1.0 - std::max(0.05, 1.2 * ratioMaxDeviation),
                                       1.0 + std::max(0.05, 1.2 * ratioMaxDeviation));
  ratioGraph->GetXaxis()->SetTitle("T_{9} [GK]");
  ratioGraph->GetYaxis()->SetTitle("fit/input");
  ratioGraph->GetXaxis()->SetTitleSize(0.11);
  ratioGraph->GetYaxis()->SetTitleSize(0.10);
  ratioGraph->GetXaxis()->SetTitleOffset(1.0);
  ratioGraph->GetYaxis()->SetTitleOffset(0.55);
  ratioGraph->GetXaxis()->SetLabelSize(0.09);
  ratioGraph->GetYaxis()->SetLabelSize(0.08);
  ratioGraph->Draw("APL");
  if (compare.enabled) {
    std::vector<double> xRefRatio;
    std::vector<double> yRefRatio;
    xRefRatio.reserve(data.size());
    yRefRatio.reserve(data.size());
    for (const auto &p : data) {
      double ref = 0.0;
      for (int term = 0; term < compare.nterms; ++term) {
        ref += reaclibRate(p.t9, compare.parameters + 7 * term);
      }
      if (std::isfinite(ref) && ref > 0.0) {
        xRefRatio.push_back(p.t9);
        yRefRatio.push_back(ref / p.rate);
      }
    }
    auto *refRatioGraph =
        new TGraph(static_cast<int>(xRefRatio.size()), xRefRatio.data(), yRefRatio.data());
    refRatioGraph->SetLineColor(kBlue + 1);
    refRatioGraph->SetMarkerColor(kBlue + 1);
    refRatioGraph->SetLineStyle(2);
    refRatioGraph->SetLineWidth(3);
    refRatioGraph->SetMarkerStyle(24);
    refRatioGraph->SetMarkerSize(0.7);
    refRatioGraph->Draw("LP SAME");
  }

  auto *unityLine = new TLine(xmin, 1.0, xmax, 1.0);
  unityLine->SetLineColor(kRed + 1);
  unityLine->SetLineStyle(2);
  unityLine->SetLineWidth(2);
  unityLine->Draw();

  canvas->Update();
  return canvas;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    printUsage(argv[0]);
    return 1;
  }

  try {
    ProgramOptions options;
    int argi = 1;
    if (argi < argc &&
        (std::string(argv[argi]) == "--config" || std::string(argv[argi]) == "-c")) {
      if (argi + 1 >= argc) {
        throw std::runtime_error("--config/-c requires a file path");
      }
      options = loadConfigFile(resolveWithDefaultDir(argv[argi + 1], "conf"));
      argi += 2;
    }

    if (argi < argc && argv[argi][0] != '-') {
      options.inputPath = argv[argi++];
    }

    for (int i = argi; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--config" || arg == "-c") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--config/-c requires a file path");
        }
        ProgramOptions fileOptions = loadConfigFile(resolveWithDefaultDir(argv[++i], "conf"));
        if (options.inputPath.empty()) {
          options.inputPath = fileOptions.inputPath;
        }
        if (!options.useTmin && fileOptions.useTmin) {
          options.tmin = fileOptions.tmin;
          options.useTmin = true;
        }
        if (!options.useCompare && fileOptions.useCompare) {
          options.comparePath = fileOptions.comparePath;
          options.useCompare = true;
        }
        if (options.termsPath.empty()) {
          options.termsPath = fileOptions.termsPath;
        }
        continue;
      } else if (arg[0] != '-') {
        throw std::runtime_error("Unexpected positional argument: " + arg);
      }

      if (arg == "--tmin") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--tmin requires a numeric value");
        }
        options.tmin = std::stod(argv[++i]);
        options.useTmin = true;
      } else if (arg == "--input") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--input requires a file path");
        }
        options.inputPath = resolveWithDefaultDir(argv[++i], "input");
      } else if (arg == "--compare") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--compare requires a file path");
        }
        options.comparePath = resolveWithDefaultDir(argv[++i], "input");
        options.useCompare = true;
      } else if (arg == "--output") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--output requires a file path");
        }
        options.outputPath = resolveWithDefaultDir(argv[++i], "output");
        options.useOutput = true;
      } else if (arg == "--no-gui") {
        options.noGui = true;
      } else if (arg == "--terms") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--terms requires a file path");
        }
        options.termsPath = resolveWithDefaultDir(argv[++i], "conf");
      } else {
        throw std::runtime_error("Unknown option: " + arg);
      }
    }

    if (options.inputPath.empty()) {
      throw std::runtime_error("Input file is required");
    }
    if (options.termsPath.empty()) {
      throw std::runtime_error("Terms file is required");
    }
    options.inputPath = resolveWithDefaultDir(options.inputPath, "input");
    options.termsPath = resolveWithDefaultDir(options.termsPath, "conf");
    if (options.useCompare) {
      options.comparePath = resolveWithDefaultDir(options.comparePath, "input");
    }
    if (options.useOutput) {
      options.outputPath = resolveWithDefaultDir(options.outputPath, "output");
    }

    const std::string inputPath = options.inputPath;
    ReactionGlobals globals;
    if (options.useGlobalReaction) {
      globals.A1 = options.A1;
      globals.A2 = options.A2;
      globals.Z1 = options.Z1;
      globals.Z2 = options.Z2;
      globals.Threshold = options.Threshold;
      if (globals.A1 <= 0.0 || globals.A2 <= 0.0) {
        throw std::runtime_error("Global inputs require positive A1 and A2");
      }
      globals.Ared = globals.A1 * globals.A2 / (globals.A1 + globals.A2);
    }

    const auto table = loadTable(inputPath);
    const auto &data = table.data;
    std::vector<DataPoint> fitData = options.useTmin ? filterByTmin(data, options.tmin) : data;
    if (fitData.empty()) {
      throw std::runtime_error("No positive-rate data rows satisfy the requested tmin");
    }
    const CompareCurve compare =
        options.useCompare ? loadCompareCurve(options.comparePath) : CompareCurve{};
    const std::vector<TermSpec> terms = loadTermFile(options.termsPath, globals);
    FitConfig cfg;
    const FitResult result = runTermSpecsFit(fitData, terms, cfg);

    std::cout << std::setprecision(12);
    std::cout << "Input points: " << data.size() << "\n";
    if (table.skippedNonPositiveRates > 0) {
      std::cout << "Skipped non-positive rates: " << table.skippedNonPositiveRates << "\n";
    }
    if (options.useTmin) {
      std::cout << "Fit points with T9 >= " << options.tmin << ": " << fitData.size() << "\n";
    }
    std::cout << "\nREACLIB form\n";
    std::cout
        << "  rate(T9) = exp(a0 + a1*T9^-1 + a2*T9^(-1/3) + a3*T9^(1/3) + a4*T9"
        << " + a5*T9^(5/3) + a6*ln(T9))\n";
    const char *names[7] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6"};
    if (!result.converged) {
      std::cout << "\nWARNING: " << modeName(result.mode) << " did not fully converge\n";
    }
    std::cout << "\nSummed terms: " << result.nterms << "\n";
    std::cout << "Objective (sum of log-residual^2): " << result.objective << "\n";
    std::cout << "RMS log error: " << result.rmsLog << "\n";
    std::cout << "Max fractional error: " << std::fixed << std::setprecision(2)
              << 100.0 * result.maxFracErr << std::defaultfloat << " %\n";
    for (int term = 0; term < result.nterms; ++term) {
      if (result.nterms > 1) {
        std::cout << "  term " << (term + 1) << "\n";
      }
      for (int i = 0; i < 7; ++i) {
        std::cout << "    " << names[i] << " = " << std::scientific << std::setprecision(6)
                  << result.parameters[7 * term + i] << std::defaultfloat << "\n";
      }
    }
    if (compare.enabled) {
      std::cout << "\nReference parameters\n";
      for (int term = 0; term < compare.nterms; ++term) {
        if (compare.nterms > 1) {
          std::cout << "  term " << (term + 1);
          if (term < static_cast<int>(compare.termKinds.size()) &&
              !compare.termKinds[term].empty()) {
            std::cout << " (" << compare.termKinds[term] << ")";
          }
          std::cout << "\n";
        }
        for (int i = 0; i < 7; ++i) {
          std::cout << "    " << names[i] << " = " << std::scientific << std::setprecision(6)
                    << compare.parameters[7 * term + i] << std::defaultfloat << "\n";
        }
      }
    }
    if (options.useOutput) {
      writeFitSummary(options.outputPath, options, result, compare, static_cast<int>(data.size()),
                      table.skippedNonPositiveRates, static_cast<int>(fitData.size()));
      std::cout << "\nSummary written to: " << options.outputPath << "\n";
    }

    if (!options.noGui) {
      int rootArgc = 1;
      char appName[] = "reaclib_fit";
      char *rootArgv[] = {appName, nullptr};
      TApplication app("reaclib_fit_app", &rootArgc, rootArgv);
      if (options.plotTmin <= 0.0 || options.plotTmax <= options.plotTmin) {
        throw std::runtime_error("Plot range requires 0 < plot_tmin < plot_tmax");
      }
      TCanvas *canvas = makePlot(data, result, compare, baseName(inputPath),
                                 baseName(options.comparePath), options.plotTmin,
                                 options.plotTmax);
      canvas->Draw();
      app.Run();
    }

    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "ERROR: " << ex.what() << "\n";
    return 1;
  }
}
