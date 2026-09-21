#include "doctest/doctest.h"
#include "moldyn/io/nist_config.hpp"
#include <stdexcept>
#include <string>

static std::string dataPath(const std::string& rel) {
    return std::string(MOLDYN_DATA_DIR) + "/" + rel;
}

TEST_CASE("reads NIST LJ configuration 4") {
    auto cfg = moldyn::readNistConfig(dataPath("nist/lj/lj_sample_config_periodic4.txt"));
    CHECK(cfg.system.size() == 30);
    CHECK(cfg.system.box().lengths().x == doctest::Approx(8.0));
    // first atom from the file
    CHECK(cfg.system.position(0).x == doctest::Approx(1.077169909511));
    CHECK(cfg.system.position(0).y == doctest::Approx(-1.020988125886));
    // last atom from the file
    CHECK(cfg.system.position(29).z == doctest::Approx(-1.252452130644));
    CHECK(cfg.types.empty());
}

TEST_CASE("reads NIST LJ configuration 1") {
    auto cfg = moldyn::readNistConfig(dataPath("nist/lj/lj_sample_config_periodic1.txt"));
    CHECK(cfg.system.size() == 800);
    CHECK(cfg.system.box().lengths().x == doctest::Approx(10.0));
}

TEST_CASE("reads SPC/E configuration and keeps atom types") {
    auto cfg = moldyn::readNistConfig(dataPath("nist/spce/spce_sample_config_periodic1.txt"));
    CHECK(cfg.system.size() == 300);            // 100 molecules x 3 sites
    CHECK(cfg.types.size() == 300);
    CHECK(cfg.types[0] == "O");
    CHECK(cfg.types[1] == "H");
    CHECK(cfg.types[2] == "H");
}

TEST_CASE("missing file throws") {
    CHECK_THROWS_AS(moldyn::readNistConfig("/nonexistent/path.txt"), std::runtime_error);
}
