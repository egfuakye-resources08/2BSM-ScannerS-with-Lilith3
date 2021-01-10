#include "ScannerS/DataMap.hpp"
#include "catch.hpp"

using Catch::Matchers::Contains;

TEST_CASE("DataMap", "[datamap][unit]") {

  auto m = ScannerS::DataMap();

  SECTION("Store, retrieve, iterate") {
    m.Store("test", 2);
    m.Store("test1", -0.1);

    REQUIRE(m["test"] == Approx(2));
    REQUIRE(m["test1"] == Approx(-0.1));

    size_t i = 0;
    for (auto [key, value] : m) {
      REQUIRE(m[key] == value);
      ++i;
    }
    REQUIRE(i == 2);
  }

  SECTION("Merge") {
    auto map = ScannerS::DataMap::Map{{"merge1", 10}, {"merge2", 20}};
    m.Merge(std::move(map));
  }

  SECTION("Errors") {
    REQUIRE_THROWS_WITH(m["unknown_key"], "Unknown key unknown_key");

    m.Store("duplicate", 2);
    REQUIRE_THROWS_WITH(m.Store("duplicate", 3),
                        Contains("duplicate already exists"));

    auto map = ScannerS::DataMap::Map{{"entry", 10}, {"duplicate", -1}};
    REQUIRE_THROWS_WITH(m.Merge(std::move(map)),
                        Contains("Entry {duplicate, -1") &&
                            Contains("can't be merged. Exists with value 2"));
  }
}
