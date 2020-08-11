
#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/place_pwlin.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/util/optional.hpp>

#include <arbornml/arbornml.hpp>
#include <arbornml/nmlexcept.hpp>
#include <arbornml/with_xml.hpp>

#include "../test/gtest.h"
#include "morph_pred.hpp"

using testing::region_eq;

TEST(neuroml, with_xml) {
    // This (hopefully) will not blow up.
    {
        arbnml::with_xml scope;
    }
    {
        arbnml::with_xml scope;
    }
}

// Tests for basic morphology scanning and collection from XML.

TEST(neuroml, morph_badxml) {
    std::string illformed = "<wha?";

    EXPECT_THROW(arbnml::neuroml{illformed}, arbnml::xml_error);
}

TEST(neuroml, morph_none) {
    // No NeuroML doc, with and without declaration:
    {
        std::string empty1 = R"~(<?xml version="1.0" encoding="UTF-8"?><foo/>)~";

        arbnml::neuroml N1(empty1);
        EXPECT_TRUE(N1.cell_ids().empty());
        EXPECT_TRUE(N1.morphology_ids().empty());

        std::string empty2 = "<foo/>";

        arbnml::neuroml N2(empty2);
        EXPECT_TRUE(N2.cell_ids().empty());
        EXPECT_TRUE(N2.morphology_ids().empty());
    }

    // Empty NeuroML doc:
    {
        std::string empty3 =
R"~(<?xml version="1.0" encoding="UTF-8"?>
<neuroml xmlns="http://www.neuroml.org/schema/neuroml2">
</neuroml>)~";

        arbnml::neuroml N3(empty3);
        EXPECT_TRUE(N3.cell_ids().empty());
        EXPECT_TRUE(N3.morphology_ids().empty());
    }
}

TEST(neuroml, morph_ids) {
    // Two top-level morphologies (m1 and m2);
    // cell c3 uses top-level morphology m1;
    // cell c4 uses internally defined morphology m4.
    std::string doc =
R"~(
<neuroml xmlns="http://www.neuroml.org/schema/neuroml2">
<morphology id="m1"/>
<morphology id="m2"/>
<cell id="c3" morphology="m1"/>
<cell id="c4">
    <morphology id="m4"/>
</cell>
</neuroml>
)~";

    using svector = std::vector<std::string>;

    arbnml::neuroml N(doc);

    svector m_ids = N.morphology_ids(); // only top-level!
    std::sort(m_ids.begin(), m_ids.end());
    EXPECT_EQ((svector{"m1", "m2"}), m_ids);

    svector c_ids = N.cell_ids();
    std::sort(c_ids.begin(), c_ids.end());
    EXPECT_EQ((svector{"c3", "c4"}), c_ids);

    arbnml::morphology_data mdata;

    mdata = N.cell_morphology("c4").value();
    EXPECT_EQ("c4", mdata.cell_id);
    EXPECT_EQ("m4", mdata.id);

    mdata = N.cell_morphology("c3").value();
    EXPECT_EQ("c3", mdata.cell_id);
    EXPECT_EQ("m1", mdata.id);

    EXPECT_THROW(N.cell_morphology("mr. bobbins").value(), arb::util::optional_unset_error);
}

TEST(neuroml, simple_morphologies) {
    using namespace arb;

    // Points used in morphology definitions below.

    mpoint p0{1, -2, 3.5, 4};
    mpoint p1{3, -3.5, 4, 4.25};
    mpoint p2{3, -4, 4, 2.25};
    mpoint p3{4.5, -5, 5, 0.25};

    std::string doc =
R"~(
<neuroml xmlns="http://www.neuroml.org/schema/neuroml2">
<morphology id="m1">
    <!-- Just one segment between p0 and p1. -->
    <segment id="0">
        <proximal x="1" y="-2" z="3.5" diameter="8"/>
        <distal x="3" y="-3.5" z="4" diameter="8.5"/>
    </segment>
</morphology>
<morphology id="m2">
    <!-- Two segments, implicit proximal, [p0 p1] [p1 p3]. -->
    <segment id="0">
        <proximal x="1" y="-2" z="3.5" diameter="8"/>
        <distal x="3" y="-3.5" z="4" diameter="8.5"/>
    </segment>
    <segment id="1">
        <parent segment="0"/>
        <distal x="4.5" y="-5" z="5" diameter="0.5"/>
    </segment>
</morphology>
<morphology id="m3">
    <!-- Two segments, explicit proximal (with gap)
         [p0 p1] [p2 p3]. -->
    <segment id="0" name="soma">
        <proximal x="1" y="-2" z="3.5" diameter="8"/>
        <distal x="3" y="-3.5" z="4" diameter="8.5"/>
    </segment>
    <segment id="1">
        <parent segment="0"/>
        <proximal x="3" y="-4" z="4" diameter="4.5"/>
        <distal x="4.5" y="-5" z="5" diameter="0.5"/>
    </segment>
</morphology>
<morphology id="m4">
    <!-- Two segments, meeting at root point p0,
         [p0 p1] and [p0 p3]. -->
    <segment id="0">
        <proximal x="1" y="-2" z="3.5" diameter="8"/>
        <distal x="3" y="-3.5" z="4" diameter="8.5"/>
    </segment>
    <segment id="1">
        <parent segment="0" fractionAlong="0.0"/>
        <distal x="4.5" y="-5" z="5" diameter="0.5"/>
    </segment>
</morphology>
</neuroml>
)~";

    arbnml::neuroml N(doc);

    {
        arbnml::morphology_data m1 = N.morphology("m1").value();
        label_dict labels;
        labels.import(m1.segments, "seg:");
        mprovider P(m1.morphology, labels);

        EXPECT_TRUE(region_eq(P, "seg:0", reg::all()));

        place_pwlin G(P.morphology());
        EXPECT_EQ(p0, G.at(mlocation{0, 0}));
        EXPECT_EQ(p1, G.at(mlocation{0, 1}));
    }

    {
        arbnml::morphology_data m2 = N.morphology("m2").value();
        label_dict labels;
        labels.import(m2.segments, "seg:");
        mprovider P(m2.morphology, labels);

        mextent seg0_extent = thingify(reg::named("seg:0"), P);
        ASSERT_EQ(1u, seg0_extent.size());
        mcable seg0 = seg0_extent.cables()[0];

        mextent seg1_extent = thingify(reg::named("seg:1"), P);
        ASSERT_EQ(1u, seg1_extent.size());
        mcable seg1 = seg1_extent.cables()[0];

        EXPECT_EQ(0u, seg0.branch);
        EXPECT_EQ(0.0, seg0.prox_pos);

        EXPECT_EQ(0u, seg1.branch);
        EXPECT_EQ(seg0.dist_pos, seg1.prox_pos);
        EXPECT_EQ(1.0, seg1.dist_pos);

        place_pwlin G(P.morphology());
        EXPECT_EQ(p0, G.at(prox_loc(seg0)));
        EXPECT_EQ(p1, G.at(dist_loc(seg0)));
        EXPECT_EQ(p1, G.at(prox_loc(seg1)));
        EXPECT_EQ(p3, G.at(dist_loc(seg1)));
    }

    {
        arbnml::morphology_data m3 = N.morphology("m3").value();
        label_dict labels;
        labels.import(m3.segments, "seg:");
        mprovider P(m3.morphology, labels);

        mextent seg0_extent = thingify(reg::named("seg:0"), P);
        ASSERT_EQ(1u, seg0_extent.size());
        mcable seg0 = seg0_extent.cables()[0];

        mextent seg1_extent = thingify(reg::named("seg:1"), P);
        ASSERT_EQ(1u, seg1_extent.size());
        mcable seg1 = seg1_extent.cables()[0];

        EXPECT_EQ(0u, seg0.branch);
        EXPECT_EQ(0.0, seg0.prox_pos);

        EXPECT_EQ(0u, seg1.branch);
        EXPECT_EQ(seg0.dist_pos, seg1.prox_pos);
        EXPECT_EQ(1.0, seg1.dist_pos);

        place_pwlin G(P.morphology());
        auto seg0_segments = G.segments(seg0_extent);
        auto seg1_segments = G.segments(seg1_extent);

        ASSERT_EQ(1u, seg0_segments.size());
        EXPECT_EQ(p0, seg0_segments[0].prox);
        EXPECT_EQ(p1, seg0_segments[0].dist);

        ASSERT_EQ(1u, seg1_segments.size());
        EXPECT_EQ(p2, seg1_segments[0].prox);
        EXPECT_EQ(p3, seg1_segments[0].dist);
    }
}


// also todo: parse neuroml embedded in another XML document.
