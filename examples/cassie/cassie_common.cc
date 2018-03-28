#include "drake/examples/cassie/cassie_common.h"

namespace drake {
namespace examples {
namespace cassie {

VectorX<double> RPYCassieFixedPointState() {
  VectorX<double> ret(kRPYCassieDof * 2);
  ret << 0, 0, 0.9342, 0, 0, 0, 0, 0, 0.0057, 0.0057,
          0.6726, 0.6726, -1.4100, -1.4100, -0.0374, -0.0374, 1.6493, 1.6493,-0.0289,-0.0289,
          -1.7479,-1.7479, 0, 0, 0, 0, 0, 0, 0, 0, 
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  return ret;
}

VectorX<double> RPYCassieFixedPointTorque() {
  VectorX<double> ff_torque(kCassieActuators);
  ff_torque << 0,0,0,0,0,0,0,0,0,0;
  return ff_torque;
}

Vector3<double> GetFourBarHipMountPoint() {
  Vector3<double> hip_mount_point(0, 0, 0.0045);
  return hip_mount_point;
}

Vector3<double> GetFourBarHeelMountPoint() {
  Vector3<double> heel_mount_point(0.11877, -0.0001, 0);
  return heel_mount_point;
}


}  // namespace cassie
}  // namespace examples
}  // namespace drake
