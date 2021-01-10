#include "ScannerS/Tools/SushiTables.hpp"

#include "catch.hpp"

TEST_CASE("SusHiTables", "[unit][cxn]") {
  using ScannerS::Tools::SushiTables;
  SushiTables cxns;
  SECTION("LHC7") {
    REQUIRE(cxns.GGe(264, 1.7, 1.1, SushiTables::Collider::LHC7) ==
            Approx(9.141960019748051));
    REQUIRE(cxns.GGe(1342, 0.9, 14.2, SushiTables::Collider::LHC7) ==
            Approx(0.0011328425852545317));

    REQUIRE(cxns.GGo(403, 0.7, 0.8, SushiTables::Collider::LHC7) ==
            Approx(2.9387130405460584));
    REQUIRE(cxns.GGo(77, 0.01, 1.3, SushiTables::Collider::LHC7) ==
            Approx(2.2573097484453486));

    REQUIRE(cxns.BBe(593, 2.2, SushiTables::Collider::LHC7) ==
            Approx(0.000886466074127348));
    REQUIRE(cxns.BBo(111, 1.1, SushiTables::Collider::LHC7) ==
            Approx(0.32261803241240417));
  }
  SECTION("LHC8") {
    REQUIRE(cxns.GGe(672, 2.5, 0.84, SushiTables::Collider::LHC8) ==
            Approx(1.7019383339706202));

    REQUIRE(cxns.GGo(62, 1.8, 2.4, SushiTables::Collider::LHC8) ==
            Approx(553.2567054487868));

    REQUIRE(cxns.BBe(1242, 12.42, SushiTables::Collider::LHC8) ==
            Approx(0.0004906506046017591));

    REQUIRE(cxns.BBe(528, 3.1, SushiTables::Collider::LHC8) ==
            Approx(cxns.BBo(528, 3.1, SushiTables::Collider::LHC8)));
  }
  SECTION("LHC13") {
    REQUIRE(cxns.GGe(125.09, 1, 1, SushiTables::Collider::LHC13) ==
            Approx(43.765292327836285));
    REQUIRE(cxns.GGe(316, 0.6, 1.9, SushiTables::Collider::LHC13) ==
            Approx(3.2990785499718243));

    REQUIRE(cxns.GGo(33, 2.2, 0.3, SushiTables::Collider::LHC13) ==
            Approx(4548.982963775852));
    REQUIRE(cxns.GGo(2013, 12, 34, SushiTables::Collider::LHC13) ==
            Approx(0.1623884157701723));

    REQUIRE(cxns.BBo(99, 1.3, SushiTables::Collider::LHC13) ==
            Approx(1.8792582065151802));
    REQUIRE(cxns.BBe(862, 0.9, SushiTables::Collider::LHC13) ==
            Approx(0.00016610854715806757));
  }
  SECTION("LHC14") {
    REQUIRE(cxns.GGe(2998, 3.8, 6.5, SushiTables::Collider::LHC14) ==
            Approx(0.0007218513));

    REQUIRE(cxns.GGo(758, 0.8, 0.4, SushiTables::Collider::LHC14) ==
            Approx(0.643623542615857));

    REQUIRE(cxns.BBe(2612, 100, SushiTables::Collider::LHC14) ==
            Approx(0.0025486026635684585));

    REQUIRE(cxns.BBe(164, 0.3, SushiTables::Collider::LHC14) ==
            Approx(cxns.BBo(164, 0.3, SushiTables::Collider::LHC14)));
  }
  SECTION("TEV") {
    REQUIRE(cxns.GGe(11, 0.1, 0.5, SushiTables::Collider::TEV) ==
            Approx(113.524065883));

    REQUIRE(cxns.GGo(582, 1.1, 3.6, SushiTables::Collider::TEV) ==
            Approx(0.0026037872));

    REQUIRE(cxns.BBe(141, 2.3, SushiTables::Collider::TEV) ==
            Approx(0.02443212274173205));

    REQUIRE(cxns.BBe(33, 0.7, SushiTables::Collider::TEV) ==
            Approx(cxns.BBo(33, 0.7, SushiTables::Collider::TEV)));
  }
}
