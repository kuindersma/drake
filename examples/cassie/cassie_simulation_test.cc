/// @file
///
/// Implements a simulation of the Cassie bipedal robot.

#include <memory>

#include <gflags/gflags.h>

#include "drake/common/drake_assert.h"
#include "drake/common/find_resource.h"
#include "drake/common/text_logging.h"
#include "drake/common/text_logging_gflags.h"
#include "drake/examples/cassie/cassie_common.h"
#include "drake/lcm/drake_lcm.h"
#include "drake/lcmt_cassie_state.hpp"
#include "drake/lcmt_viewer_draw.hpp"
#include "drake/manipulation/util/sim_diagram_builder.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_plant/frame_visualizer.h"
#include "drake/multibody/rigid_body_plant/rigid_body_plant.h"
#include "drake/multibody/rigid_body_tree_construction.h"
#include "drake/systems/analysis/simulator.h"
#include "drake/systems/controllers/inverse_dynamics_controller.h"
#include "drake/systems/framework/diagram.h"
#include "drake/systems/framework/diagram_builder.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/lcm/lcm_publisher_system.h"
#include "drake/systems/lcm/lcm_subscriber_system.h"
#include "drake/systems/primitives/constant_vector_source.h"

DEFINE_double(simulation_sec, std::numeric_limits<double>::infinity(),
              "Number of seconds to simulate.");
DEFINE_string(urdf, "", "Name of urdf to load");
DEFINE_bool(visualize_frames, true, "Visualize end effector frames");
DEFINE_double(target_realtime_rate, 1.0,
              "Playback speed.  See documentation for "
              "Simulator::set_target_realtime_rate() for details.");

namespace drake {
namespace examples {
namespace cassie {
namespace {
using manipulation::util::SimDiagramBuilder;
using systems::ConstantVectorSource;
using systems::Context;
using systems::Diagram;
using systems::DiagramBuilder;
using systems::FrameVisualizer;
using systems::RigidBodyPlant;
using systems::Simulator;

int DoMain() {
  drake::lcm::DrakeLcm lcm;
  SimDiagramBuilder<double> builder;

  // Adds a plant.
  RigidBodyPlant<double>* plant = nullptr;
  const char* kModelPath =
      "drake/examples/cassie/models/urdf/"
      "cassie.urdf";
  const std::string urdf =
      (!FLAGS_urdf.empty() ? FLAGS_urdf : FindResourceOrThrow(kModelPath));
  {
    auto tree = std::make_unique<RigidBodyTree<double>>();
    drake::parsers::urdf::AddModelInstanceFromUrdfFile(
        FindResourceOrThrow(
            "drake/examples/cassie/models/urdf/"
            "cassie.urdf"),
        multibody::joints::kRollPitchYaw, nullptr /* weld to frame */,
        tree.get());
    multibody::AddFlatTerrainToWorld(tree.get(), 100., 10.);
    plant = builder.AddPlant(std::move(tree));
  }
  // Creates and adds LCM publisher for visualization.
  builder.AddVisualizer(&lcm);
  builder.get_visualizer()->set_publish_period(1e-2);

  // Contact parameters
  const double kYoungsModulus = 1e8;  // Pa
  const double kDissipation = 5.0;  // s/m
  const double kStaticFriction = 0.9;
  const double kDynamicFriction = 0.5;
  systems::CompliantMaterial default_material;
  default_material.set_youngs_modulus(kYoungsModulus)
      .set_dissipation(kDissipation)
      .set_friction(kStaticFriction, kDynamicFriction);
  plant->set_default_compliant_material(default_material);

  const double kStictionSlipTolerance = 0.01;  // m/s
  const double kContactRadius = 2e-3;  // m
  systems::CompliantContactModelParameters model_parameters;
  model_parameters.characteristic_radius = kContactRadius;
  model_parameters.v_stiction_tolerance = kStictionSlipTolerance;
  plant->set_contact_model_parameters(model_parameters);


  // const RigidBodyTree<double>& tree = plant->get_rigid_body_tree();
  // const int num_joints = tree.get_num_positions();

  systems::DiagramBuilder<double>* base_builder = builder.get_mutable_builder();

  VectorX<double> constant_vector(plant->get_input_port(0).size());
  constant_vector.setZero();
  auto constant_zero_source = base_builder->AddSystem<ConstantVectorSource<double>>(constant_vector);
  constant_zero_source->set_name("zero input");

  // Connects the blank input command
  base_builder->Connect(constant_zero_source->get_output_port(),
                  plant->get_input_port(0));

  auto sys = builder.Build();

  Simulator<double> simulator(*sys);

  lcm.StartReceiveThread();
  simulator.set_publish_every_time_step(true);
  simulator.set_target_realtime_rate(1.0);

  
  simulator.get_mutable_context().get_mutable_continuous_state().
        SetFromVector(RPYCassieFixedPointState());

  simulator.Initialize();

  simulator.StepTo(1.0);

  return 0;
}

}  // namespace
}  // namespace cassie
}  // namespace examples
}  // namespace drake

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  drake::logging::HandleSpdlogGflags();
  return drake::examples::cassie::DoMain();
}
