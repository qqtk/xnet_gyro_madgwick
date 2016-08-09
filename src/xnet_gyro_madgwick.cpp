// get ticks 'msg', publish /imu topic '
// and publish 'the imu/ base_link TransformStamped'msg to tf.
//
#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/Imu.h>
#include <geometry_msgs/Quaternion.h>

#include <string>
// #include "robbase_msg/encoders.h"
#include "robbase_msg/RazorImu.h"
#include "xnet_gyro_madgwick/xdriver.h"
// #include "xnet_gyro_madgwick/xdriver.h"
#define M_PI 3.1415926

typedef struct{
	int gx, gy, gz;
	int ax, ay, az;
} imuRaw_vec_struct;
imuRaw_vec_struct imuRaw_struct;

float vf_gx, vf_gy, vf_gz;
float vf_ax, vf_ay, vf_az;
float vf_roll, vf_pitch, vf_yaw;

int sampleFreq;
// #define sampleFreq	100.0f // sample frequency in Hz
#define G_2_MPSS 1.0
#define betaDef		0.1f // 2 * proportional gain
float beta = betaDef; // 2 * proportional gain (Kp)
float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;	// quaternion of sensor frame relative to auxiliary frame

bool valid_first_tick_flag; ///
double left_ticks, right_ticks;
double left_ticks_prev, right_ticks_prev;
double delta_left_ticks, delta_right_ticks;
ros::Time current_time, last_time;

double self_x=0;
double self_y=0;
double self_th=0;
double base_width, ticks_per_meter;
ros::NodeHandle *private_n;
// IMU algorithm update

// Fast inverse square-root
float invSqrt(float x)
{
    // See: http://en.wikipedia.org/wiki/Fast_inverse_square_root
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
	return y;
}

	/*
	* get Euler angles
	* aerospace sequence, to obtain sensor attitude:
	* 1. rotate around sensor Z plane by yaw
	* 2. rotate around sensor Y plane by pitch
	* 3. rotate around sensor X plane by roll
	*/
void quat2eulerOK(void) {	//float *roll, float *pitch, float *yaw

		vf_yaw =  atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2))* 180.0/M_PI;  // psi 'yaw
		vf_pitch = asin(2 * (q0 * q2 - q3 * q1))* 180.0/M_PI; //theta = pitch
		vf_roll =  (atan2(2 * (q0 * q3 + q1 * q2) , 1 - 2* (q2 * q2 + q3 * q3) ) * 180.0/M_PI *(-1)); // phi 'roll
}

//  Euler angles in radians defined with the Aerospace sequence.
// See Sebastian O.H. Madwick report
// "An efficient orientation filter for inertial and intertial/magnetic sensor arrays" Chapter 2 Quaternion representation
void getEulerErr(float * angles) {
  float q[4]; // quaternion
  // getQ(q);
  angles[0] = atan2(2 * q[1] * q[2] - 2 * q[0] * q[3], 2 * q[0]*q[0] + 2 * q[1] * q[1] - 1) * 180/M_PI; // psi
  angles[1] = asin(2 * q[0] * q[2] - 2 * q[1] * q[3] ) * 180/M_PI; // theta'OK!!
  angles[2] = atan2(2 * q[2] * q[3] - 2 * q[0] * q[1], 2 * q[0] * q[0] + 2 * q[3] * q[3] - 1) * 180/M_PI; // phi
}
// https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
// phi, theta, psi'

// https://en.wikipedia.org/wiki/Quaternion

void getYawPitchRoll(float * ypr) {
  float q[4]; // quaternion
  float gx, gy, gz; // estimated gravity direction
  // getQ(q);
  
  gx = 2 * (q[1]*q[3] - q[0]*q[2]);
  gy = 2 * (q[0]*q[1] + q[2]*q[3]);
  gz = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  
  ypr[0] = atan2(2 * q[1] * q[2] - 2 * q[0] * q[3], 2 * q[0]*q[0] + 2 * q[1] * q[1] - 1) * 180/M_PI;
  ypr[1] = atan(gx / sqrt(gy*gy + gz*gz))  * 180/M_PI;
  ypr[2] = atan(gy / sqrt(gx*gx + gz*gz))  * 180/M_PI;
}


void MadgwickAHRSupdateIMU(float gx, float gy, float gz, float ax, float ay, float az) {
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
	qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
	qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
	qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0f * q0;
		_2q1 = 2.0f * q1;
		_2q2 = 2.0f * q2;
		_2q3 = 2.0f * q3;
		_4q0 = 4.0f * q0;
		_4q1 = 4.0f * q1;
		_4q2 = 4.0f * q2;
		_8q1 = 8.0f * q1;
		_8q2 = 8.0f * q2;
		q0q0 = q0 * q0;
		q1q1 = q1 * q1;
		q2q2 = q2 * q2;
		q3q3 = q3 * q3;

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= beta * s0;
		qDot2 -= beta * s1;
		qDot3 -= beta * s2;
		qDot4 -= beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	q0 += qDot1 * (1.0f / sampleFreq);
	q1 += qDot2 * (1.0f / sampleFreq);
	q2 += qDot3 * (1.0f / sampleFreq);
	q3 += qDot4 * (1.0f / sampleFreq);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
}


void ticksLR_callback(const robbase_msg::encoders& ticks_msg){

    if (valid_first_tick_flag == false){
        left_ticks_prev = ticks_msg.ticks_l;
        right_ticks_prev = ticks_msg.ticks_r;
        valid_first_tick_flag = true;
    }
    //    left_ticks = ticks_msg.lwheelticks;
    left_ticks = ticks_msg.ticks_l;
    right_ticks = ticks_msg.ticks_r;
}
//
int main( int argc, char* argv[] )
{
    ros::init(argc, argv, "xnet_imuRaw_node" );
    ros::NodeHandle nh;
    private_n= new ros::NodeHandle("~");

    // ros::Subscriber ticks_sub = nh.subscribe("/encoder", 20, ticksLR_callback);
    ros::Publisher imu_pub = nh.advertise<sensor_msgs::Imu>("/quatIMU", 20);
    ros::Publisher pubRazorImu = nh.advertise<robbase_msg::RazorImu>("/rpyTFimu", 1);
    ros::Publisher pubRazorEulerImu = nh.advertise<robbase_msg::RazorImu>("/rpyEulerImu", 1);

    valid_first_tick_flag = false; ///
    // nh.param("base_width", base_width, 0.5);
    double gyro_bias;
    if(!private_n->getParam("gyro_bias", gyro_bias)) {
        ROS_WARN("No gyro_bias provided - default: 0");
        gyro_bias = 0.0;
    }

	std::string imu_frame_id_;
    if(!private_n->getParam("imu_frame", imu_frame_id_)) {
        ROS_WARN("No imu_frame provided - default: imu");
        imu_frame_id_ = "imu";
    }

    int imuFreq;
    // nh.param("ticks_per_meter", ticks_per_meter,88000);
    if(!private_n->getParam("imuFreq", imuFreq)) {
        ROS_WARN("No imuFreq provided - default: 24298");
        imuFreq = 50;
    }
    sampleFreq = imuFreq;

    geometry_msgs::TransformStamped imu_transform_msg;
    tf::TransformBroadcaster imu_tf_broadcaster;
    ros::Rate loop_rate(20);
    ROS_INFO("Node base_imuetry started");

    last_time = ros::Time::now();

    while (ros::ok()) {
        double dx, dr, dist, dtheta, d_left, d_right;
        double x, y;
        current_time = ros::Time::now();

        double elapsed_dt = (current_time - last_time).toSec();
	    xdriver_get("imuRaw", &imuRaw_struct, 24);
        // 2ex'scale;
        vf_gx = imuRaw_struct.gx * 0.001;
        vf_gy = imuRaw_struct.gy * 0.001;
        vf_gz = imuRaw_struct.gz * 0.001;
        vf_ax = imuRaw_struct.ax * 0.001;
        vf_ay = imuRaw_struct.ay * 0.001;
        vf_az = imuRaw_struct.az * 0.001;
        // call :: madgwick 'AHRS 'see' .c'.h 'github'a= sivertism/madgwicks_AHRS_algorithm_c--sivertism _y16m2git
        MadgwickAHRSupdateIMU(vf_gx, vf_gy, vf_gz, vf_ax, vf_ay, vf_az);

		// ros::Time current_time = ros::Time::now();

        // publish the /imu topic
        sensor_msgs::Imu imuMsg;
        robbase_msg::RazorImu rz;
        robbase_msg::RazorImu rz_vf;
        imuMsg.header.stamp = current_time;
        imuMsg.header.frame_id = imu_frame_id_; // "imu";

        //set the velocity
        imuMsg.child_frame_id = "base_link";
        imuMsg.angular_velocity.x=vf_gx;
        imuMsg.angular_velocity.y=vf_gy;
        imuMsg.angular_velocity.z=vf_gz;

        imuMsg.linear_acceleration.x=vf_ax * G_2_MPSS;
        imuMsg.linear_acceleration.y=vf_ay * G_2_MPSS;
        imuMsg.linear_acceleration.z=vf_az * G_2_MPSS;

 // 2e' madgwick' Quaternion q0'q1'q2'q3
    double roll, pitch, yaw;
	// We use a quaternion created from yaw
    tf::Quaternion imu_quat = tf::Quaternion(q0, q1, q2, q3 );
    tf::Matrix3x3(quat).getRPY(roll, pitch, yaw);
 // tf Message
      tf::Transform transform_tf;
      transform_tf.setOrigin( tf::Vector3(0.0, 0.0, 0.0) );
      transform_tf.setRotation(imu_quat);

	// publishing 'tf: imu/ base_link
	imu_transform_msg.header.frame_id = "imu"; // imu_frame_id_;
	imu_transform_msg.header.child_frame_id = "base_link";
	imu_transform_msg.header.stamp = current_time;
	// imu_tf_broadcaster.sendTransform(imu_transform_msg);
	imu_tf_broadcaster.sendTransform( (0.0, 0.0, 0.0), (q0, q1, q2, q3), current_Time, "base_link", imu_frame_id_);
   // imu_tf_broadcaster.sendTransform( tf::StampedTransform(transform_tf, ros::Time::now(), "world","laser") );

        // imu message
        // imuMsg.orientation = imu_quat;
        imuMsg.orientation.x=q0; //TODO orientation
        imuMsg.orientation.y=q1;
        imuMsg.orientation.z=q2;
        imuMsg.orientation.w=q3;

        quat2eulerOK();
        rz_vf.header.stamp = current_time;
        rz_vf.header.frame_id = "rpyeuler";
        rz_vf.roll = vf_roll;
        rz_vf.pitch = vf_pitch;
        rz_vf.yaw = vf_yaw;

        rz.header.stamp = current_time;
        rz.header.frame_id = "rpytfval";
        rz.yaw = (float) yaw;
        rz.pitch = (float) pitch;
        rz.roll = (float) roll;
		// publish the /imu topic
        imu_pub.publish(imuMsg);
        pubRazorImu.publish(rz);
        pubRazorEulerImu.publish(rz_vf);

        last_time = current_time;

        loop_rate.sleep();
    }// end.while _ ros::ok '
}// end.main

