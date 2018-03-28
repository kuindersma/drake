#pragma once

#include "drake/common/eigen_types.h"
#include "drake/common/find_resource.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/systems/framework/diagram_builder.h"


#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_plant/rigid_body_plant.h"
#include "drake/multibody/rigid_body_tree_construction.h"



namespace drake {
namespace examples {
namespace cassie {

constexpr int kRPYCassieDof = 22;
constexpr int kCassieActuators = 10;
constexpr float kCassieFourBarDistance = 0.5012;

VectorX<double> RPYCassieFixedPointState();

VectorX<double> RPYCassieFixedPointTorque();

Vector3<double> GetFourBarHipMountPoint();

Vector3<double> GetFourBarHeelMountPoint();


template <typename T>
class CassieTree : public RigidBodyTree<T> {
 public:
  CassieTree();

  size_t getNumPositionConstraints();

  // Overide position constraints to handle 4-bar constraint
  template <typename Scalar>
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> positionConstraints(
    const KinematicsCache<Scalar>& cache) const;

 
	// Overide position constraints to handle 4-bar constraint
	template <typename Scalar>
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> positionConstraintsJacobian(
    const KinematicsCache<Scalar>& cache, bool in_terms_of_qdot) const;


	template <typename Scalar>
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> positionConstraintsJacDotTimesV(
    const KinematicsCache<Scalar>& cache) const;
};




template <typename T>
CassieTree<T>::CassieTree() {
  drake::parsers::urdf::AddModelInstanceFromUrdfFile(
    FindResourceOrThrow(
      "drake/examples/cassie/models/urdf/cassie.urdf"),
        multibody::joints::kRollPitchYaw, nullptr /* weld to frame */,
        this);

   // add terrain
  multibody::AddFlatTerrainToWorld(this, 100., 10.);
}

template <typename T>
size_t CassieTree<T>::getNumPositionConstraints() {
  return 2; // one for each relative distance constraint
}

// Overide position constraints to handle 4-bar constraint
template <typename T>
template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> CassieTree<T>::positionConstraints(
  const KinematicsCache<Scalar>& cache) const {
  
  CheckCacheValidity(cache);
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ret(getNumPositionConstraints(), 1);

  // TODO: Make these class properties. 
  const RigidBody<double>* thigh_left = this->FindBody("thigh_left");
  const RigidBody<double>* thigh_right = this->FindBody("thigh_right");
  const RigidBody<double>* heel_spring_left = this->FindBody("heel_spring_left");
  const RigidBody<double>* heel_spring_right = this->FindBody("heel_spring_right");

  Vector3<double> pos_left = transformPoints(cache, GetFourBarHipMountPoint(), *thigh_left, 0) -
                              transformPoints(cache, GetFourBarHeelMountPoint(), *heel_spring_left, 0);
  Vector3<double> pos_right = transformPoints(cache, GetFourBarHipMountPoint(), *thigh_right, 0) - 
                              transformPoints(cache, GetFourBarHeelMountPoint(), *heel_spring_right, 0);

  ret(0,0) = pos_left.norm() - kCassieFourBarDistance;
  ret(0,0) = pos_right.norm() - kCassieFourBarDistance;

  return ret;
}

 
// Overide position constraints to handle 4-bar constraint
template <typename T>
template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> CassieTree<T>::positionConstraintsJacobian(
    const KinematicsCache<Scalar>& cache, bool in_terms_of_qdot) const {

  CheckCacheValidity(cache);
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> ret(
      getNumPositionConstraints(), in_terms_of_qdot ? this->get_num_positions() : this->get_num_velocities());

  // TODO: Make these class properties. 
  const RigidBody<double>* thigh_left = this->FindBody("thigh_left");
  const RigidBody<double>* thigh_right = this->FindBody("thigh_right");
  const RigidBody<double>* heel_spring_left = this->FindBody("heel_spring_left");
  const RigidBody<double>* heel_spring_right = this->FindBody("heel_spring_right");

  Vector3<double> pos_left = transformPoints(cache, GetFourBarHipMountPoint(), *thigh_left, 0) -
                              transformPoints(cache, GetFourBarHeelMountPoint(), *heel_spring_left, 0);
  Vector3<double> pos_right = transformPoints(cache, GetFourBarHipMountPoint(), *thigh_right, 0) - 
                              transformPoints(cache, GetFourBarHeelMountPoint(), *heel_spring_right, 0);

  auto J_left = transformPointsJacobian(cache, GetFourBarHipMountPoint(), *thigh_left, 0, in_terms_of_qdot) -
                transformPointsJacobian(cache, GetFourBarHeelMountPoint(), *heel_spring_left, 0, in_terms_of_qdot);
  auto J_right = transformPointsJacobian(cache, GetFourBarHipMountPoint(), *thigh_right, 0, in_terms_of_qdot) - 
                transformPointsJacobian(cache, GetFourBarHeelMountPoint(), *heel_spring_right, 0, in_terms_of_qdot);


  ret(0) = (pos_left.transpose() / pos_left.norm()) * J_left;
  ret(1) = (pos_right.transpose() / pos_right.norm()) * J_right;

  return ret;
}


template <typename T>
template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> CassieTree<T>::positionConstraintsJacDotTimesV(
    const KinematicsCache<Scalar>& cache) const {

  CheckCacheValidity(cache);
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ret(getNumPositionConstraints(), 1);

  // TODO: Make these class properties. 
  const RigidBody<double>* thigh_left = this->FindBody("thigh_left");
  const RigidBody<double>* thigh_right = this->FindBody("thigh_right");
  const RigidBody<double>* heel_spring_left = this->FindBody("heel_spring_left");
  const RigidBody<double>* heel_spring_right = this->FindBody("heel_spring_right");

  Vector3<double> pos_left = transformPoints(cache, GetFourBarHipMountPoint(), *thigh_left, 0) -
                              transformPoints(cache, GetFourBarHeelMountPoint(), *heel_spring_left, 0);
  Vector3<double> pos_right = transformPoints(cache, GetFourBarHipMountPoint(), *thigh_right, 0) - 
                              transformPoints(cache, GetFourBarHeelMountPoint(), *heel_spring_right, 0);

  auto J_left = transformPointsJacobianDotTimesV(cache, GetFourBarHipMountPoint(), *thigh_left, 0) -
                transformPointsJacobianDotTimesV(cache, GetFourBarHeelMountPoint(), *heel_spring_left, 0);
  auto J_right = transformPointsJacobianDotTimesV(cache, GetFourBarHipMountPoint(), *thigh_right, 0) - 
                transformPointsJacobianDotTimesV(cache, GetFourBarHeelMountPoint(), *heel_spring_right, 0);


  ret(0) = (pos_left.transpose() / pos_left.norm()) * J_left;
  ret(1) = (pos_right.transpose() / pos_right.norm()) * J_right;



  return ret;
}



}  // namespace cassie
}  // namespace examples
}  // namespace drake

