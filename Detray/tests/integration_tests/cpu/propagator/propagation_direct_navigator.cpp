// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/direct_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/tracks/trajectories.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <vector>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using vector3 = test::vector3;

constexpr std::size_t cache_size{navigation::default_cache_size};

/// Fixture for Runge-Kutta Propagation with direct navigator and toy detector
class PropagatorWithRkStepperDirectNavigatorToyDetector
    : public ::testing::TestWithParam<std::tuple<scalar, vector3>> {};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepperDirectNavigatorToyDetector, direct_navigator) {
  constexpr scalar tol{2e-3f};

  // Track generator configuration
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  // Toy detector
  using detector_t = host::detector<test::toy_metadata>;
  using surface_t = typename detector_t::surface_type;

  // Runge-Kutta propagation
  using navigator_t =
      caching_navigator<detector_t, cache_size, navigation::print_inspector>;
  using direct_navigator_t = direct_navigator<detector_t>;

  // Constant magnetic field and stepper type
  using bfield_t = bfield::const_field_t<scalar>;
  using stepper_t = rk_stepper<bfield_t::view_t, test_algebra>;

  // Include helix actor to check track position/covariance
  using actor_chain_t = actor_chain<
      actor::parameter_updater<
          test_algebra, actor::pointwise_material_interactor<test_algebra>>,
      actor::surface_sequencer<surface_t>>;

  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;
  using direct_propagator_t =
      propagator<stepper_t, direct_navigator_t, actor_chain_t>;

  // Memory resource
  vecmem::host_memory_resource host_mr;

  // Build toy detector
  toy_det_config<scalar> toy_cfg =
      toy_det_config<scalar>{}.n_brl_layers(4u).n_edc_layers(7u);

  // TODO: Test the collection of portal material
  toy_cfg.use_material_maps(false);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Build mangetic field
  const bfield_t bfield = create_const_field<scalar>(std::get<1>(GetParam()));
  const bfield_t::view_t bfield_view(bfield);

  // Track generator config
  generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.phi_steps(50u).theta_steps(50u);
  trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);

  // Propagation configuration
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(std::get<0>(GetParam()));
  cfg.navigation.search_window = {3u, 3u};
  cfg.navigation.estimate_scattering_noise = false;
  propagator_t p{cfg};
  direct_propagator_t direct_p{cfg};

  // Iterate through uniformly distributed momentum directions
  for (auto track : generator_t{trk_gen_cfg}) {
    // Build actor states: the helix inspector can be shared
    actor::parameter_updater_state<test_algebra> updater_state{cfg, track};
    actor::parameter_updater_state<test_algebra> fw_updater_state{cfg, track};
    actor::parameter_updater_state<test_algebra> bw_updater_state{cfg, track};

    actor::pointwise_material_interactor<test_algebra>::state
        interactor_state{};

    vecmem::data::vector_buffer<surface_t> seqs_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    vecmem::data::vector_buffer<surface_t> seqs_forward_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    vecmem::data::vector_buffer<surface_t> seqs_backward_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    vecmem::copy m_copy;
    m_copy.setup(seqs_buffer)->wait();
    m_copy.setup(seqs_forward_buffer)->wait();
    m_copy.setup(seqs_backward_buffer)->wait();

    vecmem::device_vector<surface_t> seqs_device(seqs_buffer);
    vecmem::device_vector<surface_t> seqs_forward_device(seqs_forward_buffer);
    vecmem::device_vector<surface_t> seqs_backward_device(seqs_backward_buffer);

    actor::surface_sequencer<surface_t>::state sequencer_state(seqs_device);
    actor::surface_sequencer<surface_t>::state sequencer_forward_state(
        seqs_forward_device);
    actor::surface_sequencer<surface_t>::state sequencer_backward_state(
        seqs_backward_device);

    auto actor_states =
        detray::tie(updater_state, interactor_state, sequencer_state);

    propagator_t::state state(track, bfield_view, det);
    navigator_t::state& navigation = state.navigation();

    // Propagate the entire detector
    ASSERT_TRUE(p.propagate(state, actor_states))
        << navigation.inspector().to_string();

    if (!seqs_device.empty()) {
      auto direct_forward_actor_states = detray::tie(
          fw_updater_state, interactor_state, sequencer_forward_state);
      auto direct_backward_actor_states = detray::tie(
          bw_updater_state, interactor_state, sequencer_backward_state);

      direct_propagator_t::state direct_forward_state(track, bfield_view, det,
                                                      seqs_buffer);

      ASSERT_TRUE(direct_p.propagate(direct_forward_state,
                                     direct_forward_actor_states));

      // Check if all surfaces in the sequence are encountered
      ASSERT_TRUE(direct_forward_state.navigation().finished());
      ASSERT_EQ(sequencer_state.sequence().size(),
                sequencer_forward_state.sequence().size());
      for (unsigned int i = 0; i < sequencer_state.sequence().size(); i++) {
        ASSERT_EQ(sequencer_state.sequence().at(i),
                  sequencer_forward_state.sequence().at(i));
      }

      const auto ptc = state.stepping().particle_hypothesis();
      ASSERT_EQ(
          ptc.pdg_num(),
          direct_forward_state.stepping().particle_hypothesis().pdg_num());
      const auto q = ptc.charge();

      // The initial momentum should be higher than the momentum at the
      // last surface
      ASSERT_TRUE(track.p(q) > updater_state.bound_params().p(q));
      ASSERT_NEAR(static_cast<float>(updater_state.bound_params().p(q)),
                  static_cast<float>(fw_updater_state.bound_params().p(q)),
                  static_cast<float>(updater_state.bound_params().p(q)) * tol);

      direct_propagator_t::state direct_backward_state(
          fw_updater_state.bound_params(), bfield_view, det, seqs_buffer);
      direct_backward_state.navigation().set_direction(
          detray::navigation::direction::e_backward);
      direct_backward_state.navigation().reset();

      ASSERT_TRUE(direct_p.propagate(direct_backward_state,
                                     direct_backward_actor_states));
      // Check if all surfaces in the sequence are encountered
      ASSERT_TRUE(direct_backward_state.navigation().finished());
      ASSERT_EQ(sequencer_state.sequence().size(),
                sequencer_backward_state.sequence().size());
      for (unsigned int i = 0u; i < sequencer_state.sequence().size(); i++) {
        unsigned int j = sequencer_state.sequence().size() - 1u - i;
        ASSERT_EQ(sequencer_state.sequence().at(i),
                  sequencer_backward_state.sequence().at(j));
      }

      ASSERT_NEAR(static_cast<float>(track.p(q)),
                  static_cast<float>(bw_updater_state.bound_params().p(q)),
                  static_cast<float>(track.p(q)) * tol);
      ASSERT_TRUE(bw_updater_state.bound_params().p(q) >
                  fw_updater_state.bound_params().p(q));
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_direct_navigator_toy_detector,
    PropagatorWithRkStepperDirectNavigatorToyDetector,
    ::testing::Values(
        std::make_tuple(-100.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                2.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-600.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T})));

/// Test that the parameter updater records the same surface sequence as the
/// surface sequencer actor, which can then be used by the direct navigator
TEST(detray_propagator_direct_navigator, parameter_updater_surface_sequence) {
  // Track generator configuration
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  // Toy detector
  using detector_t = host::detector<test::toy_metadata>;
  using surface_t = typename detector_t::surface_type;

  using navigator_t = caching_navigator<detector_t, cache_size>;
  using direct_navigator_t = direct_navigator<detector_t>;

  // Constant magnetic field and stepper type
  using bfield_t = bfield::const_field_t<scalar>;
  using stepper_t = rk_stepper<bfield_t::view_t, test_algebra>;

  using actor_chain_t = actor_chain<
      actor::parameter_updater<
          test_algebra, actor::pointwise_material_interactor<test_algebra>>,
      actor::surface_sequencer<surface_t>>;

  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;
  using direct_propagator_t =
      propagator<stepper_t, direct_navigator_t, actor_chain_t>;

  // Memory resource
  vecmem::host_memory_resource host_mr;
  vecmem::copy copy;

  // Build toy detector
  toy_det_config<scalar> toy_cfg =
      toy_det_config<scalar>{}.n_brl_layers(4u).n_edc_layers(7u);
  toy_cfg.use_material_maps(false);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Build mangetic field
  const bfield_t bfield = create_const_field<scalar>(vector3{
      0.f * unit<scalar>::T, 0.f * unit<scalar>::T, 2.f * unit<scalar>::T});
  const bfield_t::view_t bfield_view(bfield);

  // Track generator config
  generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.phi_steps(20u).theta_steps(20u);
  trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);

  // Propagation configuration
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(-100.f * unit<scalar>::um);
  cfg.navigation.search_window = {3u, 3u};
  cfg.navigation.estimate_scattering_noise = false;
  propagator_t p{cfg};
  direct_propagator_t direct_p{cfg};

  // Small capacity to test the buffer overflow
  constexpr unsigned int small_capacity{2u};

  for (auto track : generator_t{trk_gen_cfg}) {
    actor::pointwise_material_interactor<test_algebra>::state
        interactor_state{};

    // Reference sequence from the surface sequencer
    vecmem::data::vector_buffer<surface_t> seqs_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    copy.setup(seqs_buffer)->wait();
    vecmem::device_vector<surface_t> seqs_device(seqs_buffer);
    actor::surface_sequencer<surface_t>::state sequencer_state(seqs_device);

    // Let the parameter updater record the surfaces
    std::vector<surface_t> updater_seq(100u);
    actor::parameter_updater_state<test_algebra> updater_state{cfg, track};
    ASSERT_FALSE(updater_state.has_surface_sequence());
    updater_state.set_surface_sequence(
        updater_seq.data(), static_cast<unsigned int>(updater_seq.size()));
    ASSERT_TRUE(updater_state.has_surface_sequence());

    propagator_t::state state(track, bfield_view, det);
    auto actor_states =
        detray::tie(updater_state, interactor_state, sequencer_state);
    ASSERT_TRUE(p.propagate(state, actor_states));

    // Same surfaces in the same order
    const auto& ref_seq = sequencer_state.sequence();
    ASSERT_FALSE(updater_state.surface_sequence_overflow());
    ASSERT_EQ(updater_state.n_recorded_surfaces(), ref_seq.size());
    for (unsigned int i = 0u; i < ref_seq.size(); i++) {
      ASSERT_EQ(updater_seq.at(i), ref_seq.at(i));
    }

    // Some tracks do not encounter any surfaces
    if (ref_seq.empty()) {
      continue;
    }

    // The direct navigator has to follow the recorded sequence
    vecmem::data::vector_buffer<surface_t> seqs_direct_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    copy.setup(seqs_direct_buffer)->wait();
    vecmem::device_vector<surface_t> seqs_direct_device(seqs_direct_buffer);
    actor::surface_sequencer<surface_t>::state sequencer_direct_state(
        seqs_direct_device);
    actor::parameter_updater_state<test_algebra> direct_updater_state{cfg,
                                                                      track};

    const vecmem::data::vector_view<surface_t> updater_seq_view(
        updater_state.n_recorded_surfaces(), updater_seq.data());
    direct_propagator_t::state direct_state(track, bfield_view, det,
                                            updater_seq_view);
    auto direct_actor_states = detray::tie(
        direct_updater_state, interactor_state, sequencer_direct_state);
    ASSERT_TRUE(direct_p.propagate(direct_state, direct_actor_states));
    ASSERT_TRUE(direct_state.navigation().finished());
    ASSERT_EQ(sequencer_direct_state.sequence().size(), ref_seq.size());

    // Overflow: count all surfaces, but do not write past the capacity
    if (ref_seq.size() > small_capacity) {
      std::vector<surface_t> small_seq(small_capacity);
      actor::parameter_updater_state<test_algebra> small_updater_state{cfg,
                                                                       track};
      small_updater_state.set_surface_sequence(small_seq.data(),
                                               small_capacity);
      sequencer_state.sequence().clear();

      propagator_t::state small_state(track, bfield_view, det);
      auto small_actor_states =
          detray::tie(small_updater_state, interactor_state, sequencer_state);
      ASSERT_TRUE(p.propagate(small_state, small_actor_states));

      ASSERT_TRUE(small_updater_state.surface_sequence_overflow());
      ASSERT_EQ(small_updater_state.n_recorded_surfaces(), ref_seq.size());
      // First entries are kept, the last slot holds the last surface
      ASSERT_EQ(small_seq.front(), ref_seq.front());
      ASSERT_EQ(small_seq.back(), ref_seq.back());
    }
  }
}

/// Fixture for Runge-Kutta Propagation with direct navigator and wire chamber
class PropagatorWithRkStepperDirectNavigatorWireChamber
    : public ::testing::TestWithParam<std::tuple<scalar, vector3>> {};

/// Test propagation in a constant magnetic field using a Runge-Kutta stepper
TEST_P(PropagatorWithRkStepperDirectNavigatorWireChamber, direct_navigator) {
  constexpr scalar tol{2e-5f};

  // Track generator configuration
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  // Toy detector
  using detector_t = host::detector<test::wire_chamber_metadata>;
  using surface_t = typename detector_t::surface_type;

  // Default navigator for comparison
  using navigator_t =
      caching_navigator<detector_t, cache_size, navigation::print_inspector>;
  using direct_navigator_t = direct_navigator<detector_t>;

  // Constant magnetic field and stepper type
  using bfield_t = bfield::const_field_t<scalar>;
  using stepper_t = rk_stepper<bfield_t::view_t, test_algebra>;

  // Include helix actor to check track position/covariance
  using actor_chain_t = actor_chain<
      actor::parameter_updater<
          test_algebra, actor::pointwise_material_interactor<test_algebra>>,
      actor::surface_sequencer<surface_t>>;

  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;
  using direct_propagator_t =
      propagator<stepper_t, direct_navigator_t, actor_chain_t>;

  // Memory resource
  vecmem::host_memory_resource host_mr;

  // Build wire chamber
  wire_chamber_config<scalar> wire_cfg{};
  auto [det, names] = build_wire_chamber<test::algebra>(host_mr, wire_cfg);

  // Build mangetic field
  const bfield_t bfield = create_const_field<scalar>(std::get<1>(GetParam()));
  const bfield_t::view_t bfield_view(bfield);

  // Truth track generation
  generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.phi_steps(50u).theta_steps(50u);
  trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);

  // Propagation configuration
  propagation::config cfg{};
  cfg.navigation.intersection.overstep_tolerance =
      static_cast<float>(std::get<0>(GetParam()));
  cfg.navigation.search_window = {5u, 5u};
  cfg.navigation.estimate_scattering_noise = false;
  propagator_t p{cfg};
  direct_propagator_t direct_p{cfg};

  // Iterate through uniformly distributed momentum directions
  for (auto track : generator_t{trk_gen_cfg}) {
    // Build actor states: the helix inspector can be shared
    actor::parameter_updater_state<test_algebra> updater_state{cfg, track};
    actor::parameter_updater_state<test_algebra> fw_updater_state{cfg, track};
    actor::parameter_updater_state<test_algebra> bw_updater_state{cfg, track};

    actor::pointwise_material_interactor<test_algebra>::state
        interactor_state{};

    vecmem::data::vector_buffer<surface_t> seqs_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    vecmem::data::vector_buffer<surface_t> seqs_forward_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    vecmem::data::vector_buffer<surface_t> seqs_backward_buffer{
        100u, host_mr, vecmem::data::buffer_type::resizable};
    vecmem::copy m_copy;
    m_copy.setup(seqs_buffer)->wait();
    m_copy.setup(seqs_forward_buffer)->wait();
    m_copy.setup(seqs_backward_buffer)->wait();

    vecmem::device_vector<surface_t> seqs_device(seqs_buffer);
    vecmem::device_vector<surface_t> seqs_forward_device(seqs_forward_buffer);
    vecmem::device_vector<surface_t> seqs_backward_device(seqs_backward_buffer);

    actor::surface_sequencer<surface_t>::state sequencer_state(seqs_device);
    actor::surface_sequencer<surface_t>::state sequencer_forward_state(
        seqs_forward_device);
    actor::surface_sequencer<surface_t>::state sequencer_backward_state(
        seqs_backward_device);

    auto actor_states =
        detray::tie(updater_state, interactor_state, sequencer_state);

    propagator_t::state state(track, bfield_view, det);
    navigator_t::state& navigation = state.navigation();

    // Propagate the entire detector
    ASSERT_TRUE(p.propagate(state, actor_states))
        << navigation.inspector().to_string();

    if (!seqs_device.empty()) {
      auto direct_forward_actor_states = detray::tie(
          fw_updater_state, interactor_state, sequencer_forward_state);
      auto direct_backward_actor_states = detray::tie(
          bw_updater_state, interactor_state, sequencer_backward_state);

      direct_propagator_t::state direct_forward_state(track, bfield_view, det,
                                                      seqs_buffer);

      ASSERT_TRUE(direct_p.propagate(direct_forward_state,
                                     direct_forward_actor_states));

      // Check if all surfaces in the sequence are encountered
      ASSERT_TRUE(direct_forward_state.navigation().finished());
      ASSERT_EQ(sequencer_state.sequence().size(),
                sequencer_forward_state.sequence().size());
      for (unsigned int i = 0; i < sequencer_state.sequence().size(); i++) {
        ASSERT_EQ(sequencer_state.sequence().at(i),
                  sequencer_forward_state.sequence().at(i));
      }

      const auto ptc = state.stepping().particle_hypothesis();
      ASSERT_EQ(
          ptc.pdg_num(),
          direct_forward_state.stepping().particle_hypothesis().pdg_num());
      const auto q = ptc.charge();

      // The initial momentum should be higher than or equal to the
      // momentum at the last surface
      ASSERT_GE(track.p(q), updater_state.bound_params().p(q));
      ASSERT_NEAR(static_cast<float>(updater_state.bound_params().p(q)),
                  static_cast<float>(fw_updater_state.bound_params().p(q)),
                  static_cast<float>(updater_state.bound_params().p(q)) * tol);

      direct_propagator_t::state direct_backward_state(
          fw_updater_state.bound_params(), bfield_view, det, seqs_buffer);
      direct_backward_state.navigation().set_direction(
          detray::navigation::direction::e_backward);
      direct_backward_state.navigation().reset();

      ASSERT_TRUE(direct_p.propagate(direct_backward_state,
                                     direct_backward_actor_states));

      // Check if all surfaces in the sequence are encountered
      ASSERT_TRUE(direct_backward_state.navigation().finished());

      ASSERT_EQ(sequencer_state.sequence().size(),
                sequencer_backward_state.sequence().size());
      for (unsigned int i = 0u; i < sequencer_state.sequence().size(); i++) {
        unsigned int j = sequencer_state.sequence().size() - 1u - i;
        EXPECT_EQ(sequencer_state.sequence().at(i),
                  sequencer_backward_state.sequence().at(j));
      }

      ASSERT_NEAR(static_cast<float>(track.p(q)),
                  static_cast<float>(bw_updater_state.bound_params().p(q)),
                  static_cast<float>(track.p(q)) * tol);
      ASSERT_GE(bw_updater_state.bound_params().p(q),
                fw_updater_state.bound_params().p(q));
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    detray_propagator_direct_navigator_wire_chamber,
    PropagatorWithRkStepperDirectNavigatorWireChamber,
    ::testing::Values(
        std::make_tuple(-100.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                2.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{0.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-400.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                                1.f * unit<scalar>::T}),
        std::make_tuple(-600.f * unit<scalar>::um,
                        vector3{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
                                1.f * unit<scalar>::T})));
