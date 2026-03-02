#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <numeric>
#include <cmath>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <chrono>
#include <utility>

namespace fs = std::filesystem;

// Functions
void theme() {
    std::cout << "\033[38;2;255;216;102m";
    std::cout << R"(
 ____  _____ __  __    ____            _             _     _____            _             _
|  _ \| ____|  \/  |  / ___|___  _ __ | |_ __ _  ___| |_  | ____|_   ____ _| |_   _  __ _| |_ ___  _ __
| | | |  _| | |\/| | | |   / _ \| '_ \| __/ _` |/ __| __| |  _| \ \ / / _` | | | | |/ _` | __/ _ \| '__|
| |_| | |___| |  | | | |__| (_) | | | | || (_| | (__| |_  | |___ \ V / (_| | | |_| | (_| | || (_) | |
|____/|_____|_|  |_|  \____\___/|_| |_|\__\__,_|\___|\__| |_____| \_/ \__,_|_|\__,_|\__,_|\__\___/|_|
)" << std::endl;
    std::cout << "\033[38;2;255;255;255m";
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Author: Marek Semancik" << std::endl;
    std::cout << "Version: 2026-02-25" << std::endl;
    std::cout << "OS: Windows 10 Pro" << std::endl;
    std::cout << "Editor: CLion" << std::endl;
    std::cout << "Encoding: UTF-8" << std::endl;
    std::cout << "Compiler: g++" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
}

void printUsage() {
    std::cout << "Short usage: " << std::endl;
    std::cout << "To start the application write: " << std::endl;
    std::cout << "<\\.contactDem.exe> < -option> <folder> <initialStep> <dumpStep> <endStep> [trimRatio]" <<std::endl;
    std::cout << "\033[38;2;255;255;255m" << std::endl;
    std::cout << "Parameter settings: " << std::endl;
    std::cout << "\033[38;2;0;255;255m";
    std::cout << "< -option>: < -default> for forceChain evaluation, < -liquid> for cohesion" << std::endl;
    std::cout << "[trimRatio]: optionally parameter for trimmed mean - default is 0.1 (10 %) " << std::endl;
    std::cout << "\033[38;2;255;255;255m" << std::endl;
}
/*
std::string joinPath(const std::string& folder, const std::string& name) {
    if (folder.empty()) return name;
    char last = folder.back();
#ifdef _WIN32
    if (last == '\\' || last == '/') return folder + name;
    return folder + "\\" + name;
#else
    if (last == '/') return folder + name;
    return folder + "/" + name;
#endif
}
*/

std::string joinPath(const std::string& folder, const std::string& name) {
    if (folder.empty()) return name;
    // std::filesystem::path handles (\ vs /) based on OS
    return (fs::path(folder) / name).string();
}

bool fileExists(const std::string& filename) {
    return fs::exists(filename);
}

bool findInitDump(const std::string& name) {
    bool starts = (name.size() >= 5 && name.compare(0, 5, "dump_") == 0);

    bool ends = (name.size() >= 5 && name.compare(name.size() - 5, 5, ".atom") == 0);

    return starts && (name.find("rev_") != std::string::npos) && ends;
}

std::string findRadius(const std::string& name) {
    std::ifstream file(name);
    if (!file.is_open()) return "";

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("ITEM: ATOMS") != std::string::npos) {
            std::istringstream ss(line);
            std::string column;
            std::vector<std::string> header;

            while (ss >> column) {
                header.push_back(column);
            }

            int radiusIndex = -1;
            for (int i = 0; i < header.size(); ++i) {
                if (header[i] == "radius") {
                    radiusIndex = i - 2;
                    break;
                }
            }

            if (radiusIndex != -1 && std::getline(file, line)) {
                std::istringstream data(line);
                std::string value;
                int currentIndex = 0;
                while (data >> value) {
                    if (currentIndex == radiusIndex) return value;
                    currentIndex++;
                }
            }
        }
    }
    return "";
}

// Structures
struct Force {
    double fx;
    double fy;
    double fz;
    double magnitude;
};

struct Positon {
    double x, y, z;
};

// Classes
class Statistics {
public:
    static double mean(const std::vector<double>& data) {
        if (data.empty()) return 0.0;
        double sum = std::accumulate(data.begin(), data.end(), 0.0);
        return sum / data.size();
    }

    static double standardDeviation(const std::vector<double>& data) {
        if (data.empty()) return 0.0;
        double m = mean(data);
        double sumSq = 0.0;
        for (double x : data) sumSq += (x - m) * (x - m);
        return std::sqrt(sumSq / data.size());
    }

    static double median(std::vector<double> data) {
        if (data.empty()) return 0.0;
        size_t n = data.size();
        std::nth_element(data.begin(), data.begin() + n / 2, data.end());
        double med = data[n / 2];
        if (n % 2 == 0) {
            std::nth_element(data.begin(), data.begin() + n / 2 - 1, data.end());
            med = (med + data[n / 2 - 1]) / 2.0;
        }
        return med;
    }

    static double modeBinned(const std::vector<double>& data, double binWidth = 1e-6) {
        if (data.empty()) return 0.0;
        std::unordered_map<long long, int> bins;
        for (double v : data) bins[std::llround(v / binWidth)]++;
        auto maxIt = std::max_element(
            bins.begin(), bins.end(),
                                      [](auto& a, auto& b) { return a.second < b.second; }
        );
        return maxIt->first * binWidth;
    }

    static double truncatedMean(std::vector<double> data, double trimRatio) {
        if (trimRatio < 0.0 || trimRatio >= 0.5) throw std::invalid_argument("trimRatio must be in <0,0.5)");
        std::sort(data.begin(), data.end());
        size_t k = static_cast<size_t>(data.size() * trimRatio);
        double sum = std::accumulate(data.begin() + k, data.end() - k, 0.0);
        return sum / (data.size() - 2 * k);
    }

    static double meanCoordinationNumber(const std::vector<double>& atomIds, int totalParticles) {
        if (atomIds.empty() || totalParticles == 0) return 0.0;
        return static_cast<double>(atomIds.size()) / totalParticles;
    }
};

class ForceChainTime {
private:
    std::string folder;
    long long initialStep, dumpStep, endStep;
    std::string pairsName = "pair_forcechains_";
    std::string wallsName = "wall_forcechains_";
    double deltaT;

public:
    ForceChainTime() = default;
    ForceChainTime(std::string folder_, long long initialStep_, long long dumpStep_, long long endStep_)
    : folder(folder_), initialStep(initialStep_), dumpStep(dumpStep_), endStep(endStep_)
    {}

    std::vector<std::pair<double, double>> idPairs;
    std::vector<Force> forces;
    std::vector<Positon> positionsX;
    std::vector<Positon> positionsY;

    void setDeltaT(double radius) {
        if (radius < 5e-4) {
            deltaT = 1e-5;
        }
        else {
            deltaT = 2e-5;
        }
    }

    void loadFromFile(const std::string& filename) {
        idPairs.clear();

        std::ifstream in(filename);
        if (!in) {
            std::cerr << "ERROR: Cannot open file " << filename << std::endl;
            return;
        }

        std::string dummy;
        for (int i = 0; i < 9; ++i) std::getline(in, dummy);

        double c_forceChain[19];

        while (in >> c_forceChain[1] >> c_forceChain[2] >> c_forceChain[3] >> c_forceChain[4]
                  >> c_forceChain[5] >> c_forceChain[6] >> c_forceChain[7] >> c_forceChain[8]
                  >> c_forceChain[9] >> c_forceChain[10] >> c_forceChain[11] >> c_forceChain[12]
                  >> c_forceChain[13] >> c_forceChain[14] >> c_forceChain[15] >> c_forceChain[16]
                  >> c_forceChain[17] >> c_forceChain[18])
        {
            idPairs.push_back(std::make_pair(c_forceChain[13], c_forceChain[14]));
            positionsX.push_back({c_forceChain[1], c_forceChain[2], c_forceChain[3]});
            positionsY.push_back({c_forceChain[4], c_forceChain[5], c_forceChain[6]});
            forces.push_back({c_forceChain[16],c_forceChain[17],c_forceChain[18], sqrt(pow(c_forceChain[16],2) + pow(c_forceChain[17],2) + pow(c_forceChain[18],2))});
        }
    }

    void evalContact () {

        std::map<std::pair<double,double>, double> pairDeltaT;
        std::map<std::pair<double,double>, double> wallDeltaT;

        long long dumpNumber = initialStep;
        while (dumpNumber <= endStep) {
            std::string outputPairs = joinPath(folder,"time_forcechain_pairs_" + std::to_string(dumpNumber) + ".dump");
            std::ofstream outPairs(outputPairs, std::ios::trunc);
            outPairs << "id[1] id[2] time[3] posX[4] posY[5] posZ[6] posX[7] posY[8] posZ[9] fx[10] fy[11] fz[12] fMag[13]" << std::endl;
            std::string outputWalls = joinPath(folder,"time_forcechain_walls_" + std::to_string(dumpNumber) + ".dump");
            std::ofstream outWalls(outputWalls, std::ios::trunc);
            outWalls << "id[1] id[2] time[3] posX[4] posY[5] posZ[6] posX[7] posY[8] posZ[9] fx[10] fy[11] fz[12] fMag[13]" << std::endl;
            std::string iPairsName = pairsName + std::to_string(dumpNumber) + ".dump";
            std::string iFilePairs = joinPath(folder, iPairsName);

            std::string iWallsName = wallsName + std::to_string(dumpNumber) + ".dump";
            std::string iFileWalls = joinPath(folder, iWallsName);

            std::cout << "File: " << iPairsName << "; " << iWallsName << std::endl;

            ForceChainTime pairs, walls; // Constructor

            pairs.loadFromFile(iFilePairs);
            walls.loadFromFile(iFileWalls);

            std::map<std::pair<double,double>, double> newPairDeltaT;
            std::map<std::pair<double,double>, double> newWallDeltaT;

            for (const auto& p : pairs.idPairs) {
                auto it = pairDeltaT.find(p); // Searches the previous map to see if this pair already exists
                if (it != pairDeltaT.end()) { // A pair was found on the map (iterator is not equal to end)
                    newPairDeltaT[p] = it->second + deltaT; // incrementation of deltaT
                } else {
                    newPairDeltaT[p] = deltaT; // new pair in dump
                }
            }

            for (const auto& p : walls.idPairs) {
                auto it = wallDeltaT.find(p);
                if (it != wallDeltaT.end()) {
                    newWallDeltaT[p] = it->second + deltaT;
                } else {
                    newWallDeltaT[p] = deltaT;
                }
            }

            pairDeltaT = newPairDeltaT;
            wallDeltaT = newWallDeltaT;

            /*
            for (const auto& [pair, dt] : pairDeltaT) {
                outPairs << pair.first << " " << pair.second << " " << dt << std::endl;
            }
            */
            for (size_t i = 0; i < pairs.idPairs.size(); ++i) {
                const auto& p = pairs.idPairs[i];
                double dt = pairDeltaT[p];

                const auto& posX = pairs.positionsX[i];
                const auto& posY = pairs.positionsY[i];
                const auto& f    = pairs.forces[i];

                outPairs << p.first << " " << p.second << " " << dt << " "
                         << posX.x << " " << posX.y << " " << posX.z << " "
                         << posY.x << " " << posY.y << " " << posY.z << " "
                         << f.fx << " " << f.fy << " " << f.fz << " " << f.magnitude
                         << std::endl;
            }
            /*
            for (const auto& [pair, dt] : wallDeltaT) {
                outWalls << pair.first << " " << pair.second << " " << dt << std::endl;
            }
            */
            for (size_t i = 0; i < walls.idPairs.size(); ++i) {
                const auto& p = walls.idPairs[i];
                double dt = wallDeltaT[p]; // deltaT pro daný pár

                const auto& posX = walls.positionsX[i];
                const auto& posY = walls.positionsY[i];
                const auto& f    = walls.forces[i];

                outWalls << p.first << " " << p.second << " " << dt << " "
                         << posX.x << " " << posX.y << " " << posX.z << " "
                         << posY.x << " " << posY.y << " " << posY.z << " "
                         << f.fx << " " << f.fy << " " << f.fz << " " << f.magnitude
                         << std::endl;
            }

            dumpNumber += dumpStep;
        }

    }

};

class DefaultData {
private:
    std::string folder;
    long long initialStep, dumpStep, endStep;
    double trimRatio;
    int totalParticles = 0;
    std::string pairsName = "time_forcechain_pairs_";
    std::string wallsName = "time_forcechain_walls_";
public:
    std::vector<double> contactTimes;
    std::vector<double> atomIds;
    DefaultData() = default;
    DefaultData(const std::string& folder_,
         long long initStep_,
         long long dumpStep_,
         long long endStep_,
         double trimRatio_ = 0.1)
    : folder(folder_), initialStep(initStep_),
      dumpStep(dumpStep_), endStep(endStep_), trimRatio(trimRatio_) {}

    void setParticleInfo(int count) {
        totalParticles = count;
    }

    void loadFromFile(const std::string& filename) {
        contactTimes.clear();
        atomIds.clear();

        std::ifstream in(filename);
        if (!in) return;

        // Skip first row
        std::string dummy;
        std::getline(in, dummy);

        double c_default[14];
        while (in >> c_default[1] >> c_default[2] >> c_default[3]
          >> c_default[4] >> c_default[5] >> c_default[6] >> c_default[7]
          >> c_default[8] >> c_default[9] >> c_default[10] >> c_default[11]
          >> c_default[12] >> c_default[13]) {
            contactTimes.push_back(c_default[3]);
            atomIds.push_back(c_default[1]);
        }
    }
    double meanContact() const { return Statistics::mean(contactTimes); }
    double stdContact() const { return Statistics::standardDeviation(contactTimes); }
    double medianContact() const { return Statistics::median(contactTimes); }
    double modeContact() const { return Statistics::modeBinned(contactTimes); }
    double trimmedMeanContact(double trimRatio) const { return Statistics::truncatedMean(contactTimes, trimRatio); }
    double meanCoordination(int totalParticles) const { return Statistics::meanCoordinationNumber(atomIds, totalParticles); }

        void runTimeEvaluation()
    {
        std::string outputFile = joinPath(folder,"contactStats.txt");
        std::ofstream out(outputFile, std::ios::trunc);
        out << "dumpNumber[1] MeanPairs[2] StdPairs[3] ModePairs[4] MedianPairs[5] TrimmedMeanPairs[6] CvPairs[7] MeanCNPairs[8] "
               "MeanWalls[9] StdWalls[10] ModeWalls[11] MedianWalls[12] TrimmedMeanWalls[13] CvWalls[14] MeanCNWalls[15] "
               "ContactTimeRatio[16]" << std::endl;

        long long dumpNumber = initialStep;
        while(dumpNumber <= endStep) {
            std::string iPairsName = pairsName + std::to_string(dumpNumber) + ".dump";
            std::string iFilePairs = joinPath(folder, iPairsName);

            std::string iWallsName = wallsName + std::to_string(dumpNumber) + ".dump";
            std::string iFileWalls = joinPath(folder, iWallsName);

            std::cout << "File: " << iPairsName << "; " << iWallsName << std::endl;

            DefaultData pairs, walls; // Constructor

            pairs.loadFromFile(iFilePairs);
            walls.loadFromFile(iFileWalls);

            double meanPairs = pairs.meanContact();
            double stdPairs = pairs.stdContact();
            double modePairs = pairs.modeContact();
            double medianPairs = pairs.medianContact();
            double trimmedMeanPairs = pairs.trimmedMeanContact(trimRatio);
            double cvPairs = stdPairs / meanPairs;
            double meanCNPairs = pairs.meanCoordination(totalParticles);

            double meanWalls = walls.meanContact();
            double stdWalls = walls.stdContact();
            double modeWalls = walls.modeContact();
            double medianWalls = walls.medianContact();
            double trimmedMeanWalls = walls.trimmedMeanContact(trimRatio);
            double cvWalls = stdWalls / meanWalls;
            double meanCNWalls = walls.meanCoordination(totalParticles);

            double contactTimeRatio = meanWalls / meanPairs;

            out << dumpNumber << " "
            << meanPairs << " " << stdPairs << " " << modePairs << " " << medianPairs << " " << trimmedMeanPairs << " "
            << cvPairs << " " << meanCNPairs << " "
            << meanWalls << " " << stdWalls << " " << modeWalls << " " << medianWalls << " " << trimmedMeanWalls << " "
            << cvWalls << " " << meanCNWalls << " "
            << contactTimeRatio << "\n";

            dumpNumber += dumpStep;
        }
        out.close();
    }
};

class CohesionData {
private:
    std::string folder;
    long long initialStep, dumpStep, endStep;
    double trimRatio;
    int totalParticles = 0;
    double particleRadius = 0.0;
    std::string pairsName = "liquid_bridges_pairs_";
    std::string wallsName = "liquid_bridges_walls_";

public:
    CohesionData() = default;
    CohesionData(const std::string& folder_,
             long long initStep_,
             long long dumpStep_,
             long long endStep_,
             double trimRatio_ = 0.1)
    : folder(folder_), initialStep(initStep_),
      dumpStep(dumpStep_), endStep(endStep_), trimRatio(trimRatio_) {}

    std::vector<double> contactTimes;
    std::vector<double> atomIds;
    std::vector<double> bridgeVolumes;

    void setParticleInfo(int count, double radius) {
        totalParticles = count;
        particleRadius = radius;
    }

    void loadFromFile(const std::string& filename) {
        contactTimes.clear();
        atomIds.clear();
        bridgeVolumes.clear();

        std::ifstream in(filename);
        if (!in) return;

        std::string dummy;
        for (int i = 0; i < 9; ++i) std::getline(in, dummy);

        double c_liquidTime[15];
        while (in >> c_liquidTime[1] >> c_liquidTime[2] >> c_liquidTime[3] >> c_liquidTime[4] >> c_liquidTime[5] >> c_liquidTime[6]
            >> c_liquidTime[7] >> c_liquidTime[8] >> c_liquidTime[9] >> c_liquidTime[10] >> c_liquidTime[11] >> c_liquidTime[12]
            >> c_liquidTime[13] >> c_liquidTime[14]) {
            contactTimes.push_back(c_liquidTime[14]);
            atomIds.push_back(c_liquidTime[7]);
            bridgeVolumes.push_back(c_liquidTime[10]);
            }
    }

    double meanContact() const { return Statistics::mean(contactTimes); }
    double stdContact() const { return Statistics::standardDeviation(contactTimes); }
    double medianContact() const { return Statistics::median(contactTimes); }
    double modeContact() const { return Statistics::modeBinned(contactTimes); }
    double trimmedMeanContact(double trimRatio) const { return Statistics::truncatedMean(contactTimes, trimRatio); }
    double meanCoordination(int totalParticles) const { return Statistics::meanCoordinationNumber(atomIds, totalParticles); }
    double meanBridgeVolume(double radius, bool isWall=false) const {
        if (bridgeVolumes.empty()) return 0.0;
        double vol = 0.0;
        for (double v : bridgeVolumes) vol += v;
        double base = isWall ? (4.0/3.0*M_PI*std::pow(radius,3) + M_PI*std::pow(radius,2))
        : 2*(4.0/3.0*M_PI*std::pow(radius,3));
        return (vol / bridgeVolumes.size()) / base;
    }
    double medianBridgeVolume(double radius, bool isWall=false) const {
        return Statistics::median(bridgeVolumes) / (isWall ? (4.0/3.0*M_PI*std::pow(radius,3) + M_PI*std::pow(radius,2))
        : 2*(4.0/3.0*M_PI*std::pow(radius,3)));
    }

    void runTimeEvaluation()
    {
        std::string outputFile = joinPath(folder,"contactStatsLiquid.txt");
        std::ofstream out(outputFile, std::ios::trunc);
        out << "dumpNumber[1] MeanPairs[2] StdPairs[3] ModePairs[4] MedianPairs[5] TrimmedMeanPairs[6] CvPairs[7] MeanCNPairs[8] MeanBridgePairs[9] MedianBridgePairs[10] "
               "MeanWalls[11] StdWalls[12] ModeWalls[13] MedianWalls[14] TrimmedMeanWalls[15] CvWalls[16] MeanCNWalls[17] MeanBridgeWalls[18] MedianBridgeWalls[19] "
               "ContactTimeRatio[20]" << std::endl;

        long long dumpNumber = initialStep;
        while(dumpNumber <= endStep) {
            std::string iPairsName = pairsName + std::to_string(dumpNumber) + ".dump";
            std::string iFilePairs = joinPath(folder, iPairsName);

            std::string iWallsName = wallsName + std::to_string(dumpNumber) + ".dump";
            std::string iFileWalls = joinPath(folder, iWallsName);

            std::cout << "File: " << iPairsName << "; " << iWallsName << std::endl;

            CohesionData pairs, walls; // Constructor

            pairs.loadFromFile(iFilePairs);
            walls.loadFromFile(iFileWalls);

            double meanPairs = pairs.meanContact();
            double stdPairs = pairs.stdContact();
            double modePairs = pairs.modeContact();
            double medianPairs = pairs.medianContact();
            double trimmedMeanPairs = pairs.trimmedMeanContact(trimRatio);
            double cvPairs = stdPairs / meanPairs;
            double meanCNPairs = pairs.meanCoordination(totalParticles);
            double meanBridgePairs = pairs.meanBridgeVolume(particleRadius,false);
            double medianBridgePairs = pairs.medianBridgeVolume(particleRadius,false);

            double meanWalls = walls.meanContact();
            double stdWalls = walls.stdContact();
            double modeWalls = walls.modeContact();
            double medianWalls = walls.medianContact();
            double trimmedMeanWalls = walls.trimmedMeanContact(trimRatio);
            double cvWalls = stdWalls / meanWalls;
            double meanCNWalls = walls.meanCoordination(totalParticles);
            double meanBridgeWalls = walls.meanBridgeVolume(particleRadius,true);
            double medianBridgeWalls = walls.medianBridgeVolume(particleRadius,true);

            double contactTimeRatio = meanWalls / meanPairs;

            out << dumpNumber << " "
            << meanPairs << " " << stdPairs << " " << modePairs << " " << medianPairs << " " << trimmedMeanPairs << " "
            << cvPairs << " " << meanCNPairs << " " << meanBridgePairs << " " << medianBridgePairs << " "
            << meanWalls << " " << stdWalls << " " << modeWalls << " " << medianWalls << " " << trimmedMeanWalls << " "
            << cvWalls << " " << meanCNWalls << " " << meanBridgeWalls << " " << medianBridgeWalls << " "
            << contactTimeRatio << "\n";

            dumpNumber += dumpStep;
        }
        out.close();
    }
};

int main(int argc, char* argv[]) {
    theme();
    if(argc < 6) {
        printUsage();
        return 1;
    }

    std::string option = argv[1];
    std::string folder = argv[2];
    long long initialStep = std::stoll(argv[3]);
    long long dumpStep = std::stoll(argv[4]);
    long long endStep = std::stoll(argv[5]);
    double trimRatio = 0.1;
    if(argc > 6) trimRatio = std::stod(argv[6]);

    bool foundInitDump = false;
    std::string particleCount;
    std::string particleRadius;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (!entry.is_regular_file())
            continue;

        if (findInitDump(entry.path().filename().string())) {
            foundInitDump = true;

            std::ifstream file(entry.path());
            if (!file.is_open()) {
                std::cerr << "Invalid dump format for read...\n";
                break;
            }

            std::string line;
            for (int i = 1; i <= 4; ++i) {
                if (!std::getline(file, line)) {
                    line.clear();
                    break;
                }
            }

            particleCount = line;

            particleRadius = findRadius(entry.path().string());
            break;
        }
    }

    if (!foundInitDump) {
        std::cerr << "No dump in format: dump_*rev_*.atom in folder: "
                  << folder << std::endl;
    }
    double radiusDouble = std::stod(particleRadius);
    const std::string outputFile = joinPath(folder, "contactStats.txt");
    std::ofstream out(outputFile, std::ios::trunc);
    std::cout << "Output file path: " << outputFile << std::endl;
    std::cout << "Particle count: " << particleCount << std::endl;
    std::cout << "Particle radius: " << particleRadius << "m" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    try {
        auto start = std::chrono::steady_clock::now();
        if (option == "-liquid")
        {
            CohesionData lq(folder, initialStep, dumpStep, endStep, trimRatio);
            lq.setParticleInfo(std::stoi(particleCount), radiusDouble);
            lq.runTimeEvaluation();

        }
        else if (option == "-default")
        {
            ForceChainTime sim(folder, initialStep, dumpStep, endStep);
            sim.setDeltaT(radiusDouble);
            std::string testFilePairs = joinPath(folder, "time_forcechain_pairs_" + std::to_string(initialStep) + ".dump");
            std::string testFileWalls = joinPath(folder, "time_forcechain_walls_" + std::to_string(initialStep) + ".dump");
            if (fileExists(testFilePairs) || fileExists(testFileWalls)) {
                DefaultData df(folder, initialStep, dumpStep, endStep, trimRatio);
                df.setParticleInfo(std::stoi(particleCount));
                df.runTimeEvaluation();
            }
            else {
                std::cout << "Time_forcechain_ does not exist.\n";
                std::cout << "Starting contact time computing...\n";
                sim.evalContact();
                std::cout << "\033[38;2;0;250;154m";
                std::cout << "Done; Starting evaluation...\n";
                std::cout << "\033[38;2;255;255;255m";
                DefaultData df(folder, initialStep, dumpStep, endStep, trimRatio);
                df.setParticleInfo(std::stoi(particleCount));
                df.runTimeEvaluation();
            }
        }
        auto end = std::chrono::steady_clock::now();
        auto diff = end - start;
        double minutes = std::chrono::duration<double, std::ratio<60>>(diff).count();
        std::cout << "\033[38;2;0;250;154m";
        std::cout << "DONE; execution time: " << minutes << " min" << std::endl;
        std::cout << "\033[38;2;255;255;255m";
        std::cout << "\nPress any key to continue..." << std::endl;
        std::cin.get();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
