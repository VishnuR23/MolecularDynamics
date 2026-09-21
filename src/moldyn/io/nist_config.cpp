#include "moldyn/io/nist_config.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace moldyn {

namespace {

Vec3 parseBoxLine(const std::string& line, const std::string& path) {
    std::istringstream iss(line);
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    if (!(iss >> lx >> ly >> lz)) {
        throw std::runtime_error("readNistConfig: malformed box line in " + path);
    }
    return Vec3{lx, ly, lz};
}

bool isBlank(const std::string& line) {
    return line.find_first_not_of(" \t\r\n") == std::string::npos;
}

}  // namespace

NistConfig readNistConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("readNistConfig: cannot open file " + path);
    }

    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("readNistConfig: empty file " + path);
    }
    const Vec3 lengths = parseBoxLine(line, path);

    if (!std::getline(in, line)) {
        throw std::runtime_error("readNistConfig: missing count line in " + path);
    }
    {
        std::istringstream iss(line);
        long declaredCount = 0;
        if (!(iss >> declaredCount)) {
            throw std::runtime_error("readNistConfig: malformed count line in " + path);
        }
    }

    NistConfig cfg{System(Box(lengths)), {}};
    std::vector<std::string> types;
    bool anyTyped = false;

    // Line 2 is not a trustworthy atom count (it is a molecule count for the
    // SPC/E files), so read atom lines until end of file instead.
    while (std::getline(in, line)) {
        if (isBlank(line)) {
            continue;
        }

        std::istringstream iss(line);
        long index = 0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(iss >> index >> x >> y >> z)) {
            throw std::runtime_error("readNistConfig: malformed atom line in " + path);
        }

        std::string type;
        const bool hasType = static_cast<bool>(iss >> type);

        cfg.system.addAtom(Vec3{x, y, z}, Vec3{0.0, 0.0, 0.0}, 1.0);
        types.push_back(hasType ? type : std::string());
        anyTyped = anyTyped || hasType;
    }

    if (anyTyped) {
        cfg.types = std::move(types);
    }

    return cfg;
}

}  // namespace moldyn
