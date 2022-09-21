#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <ros/ros.h>
#include <ros/console.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Joy.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/Point.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>
#include <quadrotor_msgs/PositionCommand.h>
#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <quadrotor_msgs/OptimalTimeAllocator.h>
#include <quadrotor_msgs/SpatialTemporalTrajectory.h>
#include <quadrotor_msgs/ReplanCheck.h>

#include <global_planner/utils/a_star.h>
#include <global_planner/utils/backward.hpp>
#include <global_planner/utils/poly_utils.h>
#include <global_planner/utils/bezier_base.h>
#include <global_planner/temporal_optimizer.h>
#include <global_planner/spatial_optimizer.h>
#include <global_planner/ddp_optimizer.h>

#include <msgs/facet3.h>
#include <msgs/polyhedron.h>
#include <msgs/corridor.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>

using namespace std;
using namespace Eigen;

namespace backward {
backward::SignalHandling sh;
}

// simulation param from launch file
double _vis_traj_width;
double _x_size, _y_size, _z_size, _resolution, _inv_resolution;
double _cloud_margin, _minimize_order;
double _MAX_Vel, _MAX_Acc, _MAX_d_Acc, _MAX_Jer;
double _w_Q, _w_snap, _w_terminal, _w_time;
double _w_terminal_line, _w_time_line;
double _w_terminal_zero, _w_time_zero, _w_snap_zero;
double _w_terminal_bez, _w_time_bez, _w_snap_bez;
int _time_power;
int    _traj_order, _max_inf_iter, _max_clu_iter, _iter_max, _iter_max_zero;
string _your_file_path;
int _fast_test;
double _relaxtol;
int _benchmark_poly_num; // read this number of polytope from a recorded corridor

// useful global variables
bool _has_odom  = false;
bool _has_map   = false;
bool _has_traj  = false;

Vector3d _start_pt,  _end_pt;
Vector3d _map_lower, _map_upper;

VectorXd init_allocTime, trr_allocTime, ddp_allocTime;
double trr_compTime, ddp_compTime;

VectorXd Results = VectorXd::Zero(11); // compTime, sum of allocTime, sum of initallocTime

vector<resultsStruct> saveResults; 
int _alg_id;
int _segNum;
double _pt_max_x, _pt_min_x, _pt_max_y, _pt_min_y, _pt_max_z, _pt_min_z;
double _rho;
int _max_x_id, _max_y_id, _max_z_id;
int _traj_id = 1;
int _vis_iter_num = 0;

vector<Vector3d> _manual_path;
decomp_ros_msgs::PolyhedronArray _poly_array_msg;

// ros related
ros::Subscriber _map_sub, _odom_sub, _joy_sub, _target_sub, _start_sub;
ros::Publisher _vis_polytope_pub, _vis_traj_pub, _vis_grid_path_pub, _vis_inf_map_pub;
ros::Publisher _space_time_pub;
ros::Publisher _corridor_rec_pub;
ros::Subscriber _corridor_rec_sub;

ros::Time _odom_time, _traj_time_start, _traj_time_final, _time_update_odom;
nav_msgs::Odometry _odom;

// useful object
Bernstein      * _bezier_basis             = new Bernstein(); 
gridPathFinder * _path_finder              = new gridPathFinder();
polyhedronGenerator * _polyhedronGenerator = new polyhedronGenerator();

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map);
void rcvOdometryCallBack(const nav_msgs::Odometry odom);
void rcvJoyCallBack(const sensor_msgs::Joy joy);
void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg);
void corridorRecCallBack(const msgs::corridor &msg);
void writeCorridorMsg(int path_id, decomp_cvx_space::FlightCorridor & corridor);
void readCorridorMsg(msgs::corridor corridor_rec_msg, decomp_cvx_space::FlightCorridor &corridor, int &path_id );

void startTrajCallBack(const geometry_msgs::PointStamped::ConstPtr &msg);
void fastTrajPlanning(int _benchmark_poly_num);
void trajPlanning(int _benchmark_poly_num);


void initTimeAllocation(decomp_cvx_space::FlightCorridor & corridor);
bool isNext(Vector3d coord_1, Vector3d coord_2);
double getEnergy(timeAllocator * time_allocator, const MatrixXd bezier_coeff, const VectorXd time);

void visCorridor( decomp_ros_msgs::PolyhedronArray & poly_array_msg );
void visGridPath();
void visBezierTrajectory(MatrixXd polyCoeff, VectorXd time, int iter_id);
void visColoredBezierTrajectory(MatrixXd polyCoeff, VectorXd time, int iter_id, Vector3d color);

void visFinalBezierTrajectory(MatrixXd polyCoeff, VectorXd time);
void visFinalColoredBezierTrajectory(MatrixXd polyCoeff, VectorXd time, Vector3d color);
void clearVisualization(int grid_path_mk_num, int traj_mk_iter_num);

quadrotor_msgs::SpatialTemporalTrajectory getSpaceTimeTraj(const timeAllocator * time_allocator, const MatrixXd & bezier_coeff, const VectorXd & range );

fstream path_record;
void rcvWaypointsCallBack(const nav_msgs::Path & wp)
{    
    _end_pt << wp.poses[0].pose.position.x,
               wp.poses[0].pose.position.y,
               wp.poses[0].pose.position.z;
}

bool _WRITE_PATH, _READ_PATH;
bool _WRITE_CORRIDOR, _READ_CORRIDOR;
bool _no_cons;
void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
{
    vector<Vector3d> connect_path;
    if (!_no_cons){

        if (!_READ_CORRIDOR && !_WRITE_CORRIDOR){
            if (abs(msg->pose.position.x)>_x_size/2 || abs(msg->pose.position.y)>_y_size/2){
                _manual_path.clear();
                _polyhedronGenerator->reset();   
                _poly_array_msg.polyhedrons.clear();
                _poly_array_msg.ids.clear();
                _vis_polytope_pub.publish(_poly_array_msg);        
                return;     
            }        
            double zGoal = 1.575;
            Vector3d g(msg->pose.position.x, msg->pose.position.y, zGoal);
            // Vector3d g(7.925, 2.525, zGoal);

            Vector3d current_pt(_odom.pose.pose.position.x, _odom.pose.pose.position.y, _odom.pose.pose.position.z);

            if (_manual_path.size() == 0) {
                _path_finder->AstarSearch(current_pt, g);
            }else{
                _path_finder->AstarSearch(_manual_path.back(), g);
            }
            connect_path = _path_finder->getPath();
            _path_finder->resetMap();

            // add polytopoes at newly found path to connect to starting position
            if( _polyhedronGenerator->corridorInsertGeneration(connect_path, _poly_array_msg) == 1 )
            {
                ROS_WARN("[trr_global_planner] Add extra polytopes to connecting path");
                visCorridor(_poly_array_msg);
            }
        } 

        if (_WRITE_CORRIDOR){
            Vector3d current_pt(-8.0, -8.0, 2.0);
            decomp_cvx_space::FlightCorridor corridor;
            for (int j = 201; j < 202; j++){
                // std::string path = "/home/dji/corridor_bag/path"+ std::to_string(j) +".bag";
                // std::string topics = " /tr_node/corridor";
                // std::string node_name = " __name:=my_record_node";
                // std::string cmd_str = "gnome-terminal -x bash -c 'rosbag record -O " + path + topics + node_name + "'";
                // int start_ret = system(cmd_str.c_str()); // #include <stdlib.h>
                Vector3d g_lst = current_pt;
                double theta_lst = 1.25*M_PI;
                default_random_engine eng;
                eng.seed(j);
                int N = 0;
                // uniform_real_distribution<double> rand_x = uniform_real_distribution<double>( 0.0, 6.0);
                // uniform_real_distribution<double> rand_y = uniform_real_distribution<double>( 0.0, 6.0 );
                // uniform_real_distribution<double> rand_z = uniform_real_distribution<double>( -1.0, 1.0 );
                uniform_real_distribution<double> rand_x = uniform_real_distribution<double>(-_x_size/2.0+1.0, _x_size/2.0-1.0);
                uniform_real_distribution<double> rand_y = uniform_real_distribution<double>(-_y_size/2.0+1.0, _y_size/2.0-1.0 );
                uniform_real_distribution<double> rand_z = uniform_real_distribution<double>(0.5, 1.8);
                uniform_real_distribution<double> rand_rad = uniform_real_distribution<double>(15.0, 25.0);
                uniform_real_distribution<double> rand_theta = uniform_real_distribution<double>(0.5, 2.5);
                while (N <= 105){
                    // double sign_advancedx = 1.0;
                    // if (N%2 == 0){
                    //     sign_advancedx = -1.0;
                    // }
                    // double sign_advancedy = 1.0;
                    // if (g_lst(2) > _z_size/2.0){
                    //     sign_advancedy = -1.0;
                    // }
                    // Vector3d g(sign_advancedx * rand_x(eng) + g_lst(0), sign_advancedy * rand_y(eng) +  g_lst(1), rand_z(eng) + g_lst(2));
                    // g(0) = std::max( -_x_size/2.0+1.0, std::min(g(0), _x_size/2.0-1.0) );
                    // g(1) = std::max( -_y_size/2.0+1.0, std::min(g(1), _y_size/2.0-1.0) );
                    // g(2) = std::max( -_z_size/2.0+0.2, std::min(g(2), _z_size/2.0-0.2) );
                    // Vector3d g(rand_x(eng), rand_y(eng), rand_z(eng));
                    double cur_rad = rand_rad(eng);
                    double cur_theta = theta_lst + rand_theta(eng);
                    Vector3d g(cur_rad * cos(cur_theta), cur_rad * sin(cur_theta), rand_z(eng));


                    cout<< "g" << g.transpose() << endl;
                    if (_manual_path.size() == 0) {
                        _path_finder->AstarSearch(current_pt, g);
                    }else{
                        _path_finder->AstarSearch(_manual_path.back(), g);
                    }
                    connect_path = _path_finder->getPath();
                    _path_finder->resetMap();        
                    if( _polyhedronGenerator->corridorInsertGeneration(connect_path, _poly_array_msg) == 1 )
                    {
                        ROS_WARN("[trr_global_planner] Add extra polytopes to connecting path");
                    }
                    corridor = _polyhedronGenerator->getCorridor();
                    N = corridor.polyhedrons.size();
                    cout << "current number of polyhedrons:" << N << endl;
                    for(auto coord: connect_path) _manual_path.push_back(coord);
                    visGridPath();
                    g_lst = g;
                    theta_lst = cur_theta;
                }
                // system("cd /home/dji/corridor_bag/ && rosbag record /tr_node/corridor");



                writeCorridorMsg(j, corridor);
                ROS_WARN("generated %d -th corridor with %d polytopes", j, corridor.polyhedrons.size());

                _manual_path.clear();
                connect_path.clear();
                corridor.clear();
                _polyhedronGenerator->reset();
            }

        } else {
            for(auto coord: connect_path) _manual_path.push_back(coord);
            vector<int> compList;
            decomp_cvx_space::FlightCorridor corridor;
            corridor = _polyhedronGenerator->getCorridor();
            int nn = corridor.polyhedrons.size();     
            visGridPath();
            fastTrajPlanning(nn);
            trajPlanning(nn);

            // for (int i=1; i<15; i++){
            //     compList.push_back(i*5);
            // }
            // for (int i=0; i<compList.size();i++){
            //     _benchmark_poly_num = compList[i];
            //     fastTrajPlanning(_benchmark_poly_num);
            //     trajPlanning(_benchmark_poly_num);
            // }


            // ROS_WARN("--- comparison ---");
            // printf("\n");
            // printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s","NumofPoly", "initTime", "trrCTime", "trrTime",
            //                                                         "ddpCTime","ddpTime","ddpPos","ddpRtn","ddpTerm");
            // printf("\n");
            // for (int i=0; i<compList.size();i++){
            //     printf("%-12.2d%-12.4g%-12.4g%-12.4g%-12.4g%-12.4g%-12.4g%-12.4g%-12.4g\n", compList[i], initResults[i], trrResults[i](0),trrResults[i](1),
            //                                                         ddpResults[i](0),ddpResults[i](1),ddpResults[i](2),ddpResults[i](3),ddpResults[i](4) );

            // }
            // cout << "init_allocTime = " << init_allocTime.transpose() << endl;
            // cout << "trr_allocTime = " << trr_allocTime.transpose() << endl;
            // cout << "ddp_allocTime = " << ddp_allocTime.transpose() << endl;
            // cout << "init-trr = " << (init_allocTime - trr_allocTime).sum() << endl;
            // cout << "init-ddp = " << (init_allocTime - ddp_allocTime).sum() << endl;
            // cout << "trr_compTime = " << trr_compTime << endl;
            // cout << "ddp_compTime = " << ddp_compTime << endl;
        }
    } else {
        vector<int> compList;
        for (int i=10; i<11; i++){
            compList.push_back(_segNum);
        }
        for (int i=0; i<compList.size();i++){
            _benchmark_poly_num = compList[i];
            if (_traj_order == 5){
            }
            if (_traj_order == 7){
            }

        }
    }


}

decomp_cvx_space::FlightCorridor corridor_recd;
void corridorRecCallBack(const msgs::corridor &msg)
{   
    corridor_recd.clear();
    int path_id;
    readCorridorMsg(msg, corridor_recd, path_id);
    ROS_WARN("received recorded corridor msg path_id: %d with %d polytopes", path_id, corridor_recd.polyhedrons.size());
    if (_READ_CORRIDOR){
        vector<int> compList;
        for (int i=2; i<=64; i++){
        // for (int i=_fast_test; i<=_fast_test; i++){
            compList.push_back(i);
        }
        std::string save_file_name = "/home/dji/comparisons/alg"+ std::to_string(_alg_id) + "path" + std::to_string(path_id);
        // file.open(save_file_name);
        FILE *fout = fopen(save_file_name.c_str(), "w");

        resultsStruct temp;

        for (int i=0; i<compList.size();i++){
            _benchmark_poly_num = compList[i];
            switch(_alg_id){
                case 0:
                fastTrajPlanning(_benchmark_poly_num);
                break; 
                case 1:
                trajPlanning(_benchmark_poly_num);
                break; 

            }


            temp.num_of_poly = _benchmark_poly_num;
            temp.compTime = Results(0);
            temp.allocTime = Results(1);
            temp.initallocTime = Results(2);
            // fprintf(fout, "%d %f %f %f\n",temp.num_of_poly,temp.compTime,temp.allocTime,temp.initallocTime);
            // file<<temp.num_of_poly<<temp.compTime<<temp.allocTime<<temp.initallocTime<<endl;
            fprintf(fout, "%d %f %f %f %f %f %f %f %f %f %f %f\n",temp.num_of_poly, Results(0), Results(1), Results(2), Results(3), Results(4), Results(5), Results(6), Results(7), Results(8), Results(9), Results(10));
        }
        fclose(fout);
        
    }
}

void writeCorridorMsg(int path_id, decomp_cvx_space::FlightCorridor & corridor)
{
    // cout << "writed msg with num of polys" << corridor.polyhedrons.size() << endl;
    msgs::corridor _corridor_rec_msg;
    vector<decomp_cvx_space::Polytope> polyhedrons = corridor.polyhedrons;
    int N = polyhedrons.size();
    for (int i = 0; i < N; i++){
        msgs::polyhedron polyhedron_msg;
        decomp_cvx_space::Polytope pltp = polyhedrons[i];
        polyhedron_msg.center.x = pltp.center(0);
        polyhedron_msg.center.y = pltp.center(1);
        polyhedron_msg.center.z = pltp.center(2);
        polyhedron_msg.seed_coord.x = pltp.seed_coord(0);
        polyhedron_msg.seed_coord.y = pltp.seed_coord(1);
        polyhedron_msg.seed_coord.z = pltp.seed_coord(2);
        
        int num_plane = pltp.planes.size(); // hyperplane num of this polyhedra
        for (int j = 0; j < num_plane; j++){
            msgs::facet3 facet3_msg;
            facet3_msg.a = pltp.planes[j](0);
            facet3_msg.b = pltp.planes[j](1);
            facet3_msg.c = pltp.planes[j](2);
            facet3_msg.d = pltp.planes[j](3);
            polyhedron_msg.coeffPolyhedron.push_back(facet3_msg);
        }
        _corridor_rec_msg.polyhedrons.push_back(polyhedron_msg);
    }
    _corridor_rec_msg.path_id = path_id;
    _corridor_rec_pub.publish(_corridor_rec_msg);
}

void readCorridorMsg(msgs::corridor corridor_rec_msg, decomp_cvx_space::FlightCorridor &corridor, int &path_id )
{
    vector<decomp_cvx_space::Polytope> polyhedrons;
    int N = corridor_rec_msg.polyhedrons.size();
    for (int i = 0; i < N; i++){
        msgs::polyhedron polyhedron_msg = corridor_rec_msg.polyhedrons[i];
        decomp_cvx_space::Polytope pltp;
        pltp.center(0) = polyhedron_msg.center.x;
        pltp.center(1) = polyhedron_msg.center.y;
        pltp.center(2) = polyhedron_msg.center.z;
        pltp.seed_coord(0) = polyhedron_msg.seed_coord.x;
        pltp.seed_coord(1) = polyhedron_msg.seed_coord.y;
        pltp.seed_coord(2) = polyhedron_msg.seed_coord.z;

        int num_plane = polyhedron_msg.coeffPolyhedron.size(); // hyperplane num of this polyhedra

        for (int j = 0; j < num_plane; j++){
            msgs::facet3 facet3_msg = polyhedron_msg.coeffPolyhedron[j];
            Vector4d coeff(facet3_msg.a, facet3_msg.b, facet3_msg.c, facet3_msg.d);
            pltp.planes.push_back(coeff);
        }        
        polyhedrons.push_back(pltp);
    }
    corridor.polyhedrons = polyhedrons;
    path_id = corridor_rec_msg.path_id;
}


void startTrajCallBack(const geometry_msgs::PointStamped::ConstPtr &msg)
{

}

void rcvJoyCallBack(const sensor_msgs::Joy joy)
{   
    if(joy.buttons[7] == 1.0)
    {
        ROS_WARN("[trr_global_planner] Enter in autonomous mode");

        if(_manual_path.size() == 0) return;

        Vector3d current_pt(_odom.pose.pose.position.x, _odom.pose.pose.position.y, _odom.pose.pose.position.z);
        _path_finder->AstarSearch(current_pt, _manual_path.front());
        vector<Vector3d> connect_path = _path_finder->getPath();
        _path_finder->resetMap();

        // add polytopoes at newly found path to connect to starting position
        // if( _polyhedronGenerator->corridorInsertGeneration(connect_path, _poly_array_msg) == 1 )
        // {
        //     ROS_WARN("[trr_global_planner] Add extra polytopes to connecting path");
        //     visCorridor(_poly_array_msg);
        // }
        
        for(auto coord: _manual_path) connect_path.push_back(coord);
        // _manual_path = connect_path;
        visGridPath( );

        // trajPlanning(); 
    }
    else if(joy.buttons[6] == 1.0)
    {
        ROS_WARN("[trr_global_planner] Enter in maunal flight mode");
        clearVisualization(_manual_path.size(), _vis_iter_num);

        // _manual_path.clear();
        _polyhedronGenerator->reset();
        _has_traj = false;
        _vis_iter_num = 0;
    }
}

Vector3i _pose_idx, _pose_idx_lst;
void rcvOdometryCallBack(const nav_msgs::Odometry odom)
{   
    if(!_has_map)
        return;

//    cout<<"odom"<<endl;
    _odom = odom;
    _odom_time = odom.header.stamp;
    _time_update_odom = ros::Time::now();

    _has_odom = true;

    Vector3d pos, pos_round;
    pos(0)  = odom.pose.pose.position.x;
    pos(1)  = odom.pose.pose.position.y;
    pos(2)  = odom.pose.pose.position.z;    

    _pose_idx = _path_finder->coord2gridIndex(pos);
    vector<Vector3d> coord_set;

    if(_manual_path.size() == 0)
    {   
        pos_round = _path_finder->gridIndex2coord(_pose_idx);
        _manual_path.push_back(pos_round);
        _pose_idx_lst = _pose_idx;

        coord_set.push_back(pos_round);
    }
    else if( _pose_idx != _pose_idx_lst ) 
    {   
        if( _path_finder->IndexQuery(_pose_idx) > 0 ) 
            return;
        
        pos_round = _path_finder->gridIndex2coord(_pose_idx);
        
        if(isNext(_manual_path.back(), pos_round) == false)
        {
            _path_finder->AstarSearch(_manual_path.back(), pos_round);
            
            vector<Vector3d> localPath = _path_finder->getPath();

            for(auto coord: localPath)
            {   
                coord_set.   push_back(coord);
                // _manual_path.push_back(coord);
            }

            _path_finder->resetMap();
        }
        else
        {
            coord_set.   push_back(pos_round);
            // _manual_path.push_back(pos_round);
        }    

        _pose_idx_lst = _pose_idx;
    }
    else
        return;

    if( _has_traj || _READ_PATH) return;

    // if( _polyhedronGenerator->corridorIncreGeneration(coord_set, _poly_array_msg) == 1 )
    //     visCorridor(_poly_array_msg);
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map)
{   
    if( _has_map ) return;

    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::PointCloud<pcl::PointXYZ> cloud_inf;
    pcl::fromROSMsg(pointcloud_map, cloud);
    sensor_msgs::PointCloud2 map_inflation;

    if( (int)cloud.points.size() == 0)
        return;

    pcl::PointXYZ pt, pt_inf;

    int inf_step   = round(_cloud_margin * _inv_resolution);
    int inf_step_z = max(1, inf_step / 2);
    for (int idx = 0; idx < (int)cloud.points.size(); idx++)
    {    
        pt = cloud.points[idx];        
        for(int x = -inf_step ; x <= inf_step; x ++ )
        {
            for(int y = -inf_step ; y <= inf_step; y ++ )
            {
                for(int z = -inf_step_z; z <= inf_step_z; z ++ )
                {
                    double inf_x = pt.x + x * _resolution;
                    double inf_y = pt.y + y * _resolution;
                    double inf_z = pt.z + z * _resolution;

                    Vector3d vec_inf(inf_x, inf_y, inf_z);
                    Vector3i idx_inf = _path_finder->coord2gridIndex(vec_inf);

                    // set in obstacle points
                    _path_finder->setObs(inf_x, inf_y, inf_z);
                    _polyhedronGenerator->setObs(idx_inf);

                    // rounding for visualizing the grid map
                    Vector3d coor_round = _path_finder->gridIndex2coord( idx_inf );
                    pt_inf.x = coor_round(0);
                    pt_inf.y = coor_round(1);
                    pt_inf.z = coor_round(2);
                    cloud_inf.points.push_back(pt_inf);
                }
            }
        }
    }

    _polyhedronGenerator->finishMap();

    cloud_inf.width    = cloud_inf.points.size();
    cloud_inf.height   = 1;
    cloud_inf.is_dense = true;

    pcl::toROSMsg(cloud_inf, map_inflation);
    map_inflation.header.frame_id = "/map";
    _vis_inf_map_pub.publish(map_inflation);

    _has_map = true;
}

void initTimeAllocation(decomp_cvx_space::FlightCorridor & corridor)
{   
    VectorXd time;
    time.resize((int)corridor.polyhedrons.size());

    vector<Vector3d> points;
    points.push_back (_start_pt);

    for(int i = 1; i < (int)corridor.polyhedrons.size(); i++)
        points.push_back(corridor.polyhedrons[i].seed_coord);

    points.push_back (_end_pt);

    double _Vel = _MAX_Vel;
    double _Acc = _MAX_Acc;

    for (int k = 0; k < (int)points.size() - 1; k++)
    {
        double dtxyz;
        Vector3d p0   = points[k];        
        Vector3d p1   = points[k + 1];    
        Vector3d d    = p1 - p0;          
        Vector3d v0(0.0, 0.0, 0.0);       
        double D    = d.norm();                  
        double V0   = v0.dot(d / D);             
        double aV0  = fabs(V0);           

        double acct = (_Vel - V0) / _Acc * ((_Vel > V0)?1:-1);
        double accd = V0 * acct + (_Acc * acct * acct / 2) * ((_Vel > V0)?1:-1);
        double dcct = _Vel / _Acc;                                              
        double dccd = _Acc * dcct * dcct / 2;                                   

        if (D < aV0 * aV0 / (2 * _Acc))
        {               
            double t1 = (V0 < 0)?2.0 * aV0 / _Acc:0.0;
            double t2 = aV0 / _Acc;
            dtxyz     = t1 + t2;                 
        }
        else if (D < accd + dccd)
        {
            double t1 = (V0 < 0)?2.0 * aV0 / _Acc:0.0;
            double t2 = (-aV0 + sqrt(aV0 * aV0 + _Acc * D - aV0 * aV0 / 2)) / _Acc;
            double t3 = (aV0 + _Acc * t2) / _Acc;
            dtxyz     = t1 + t2 + t3;    
        }
        else
        {
            double t1 = acct;                              
            double t2 = (D - accd - dccd) / _Vel;
            double t3 = dcct;
            dtxyz     = t1 + t2 + t3;                                                                  
        }

        time(k) = dtxyz; 
        corridor.appendTime( time(k) );   
    }
}

bool isNext(Vector3d coord_1, Vector3d coord_2)
{
    Vector3i index_1 = _path_finder->coord2gridIndex(coord_1);
    Vector3i index_2 = _path_finder->coord2gridIndex(coord_2);

    if( abs(index_1(0) - index_2(0)) <= 1 
     && abs(index_1(1) - index_2(1)) <= 1 
     && abs(index_1(2) - index_2(2)) <= 1 )
        return true;

    return false;
}

double getEnergy(timeAllocator * time_allocator, const MatrixXd bezier_coeff, const VectorXd time)
{
    double jerk_square_inte = 0;
    int traj_seg_num = bezier_coeff.rows();

    for(int i = 0; i < traj_seg_num; i++)
    {
        double jerk_square_inte_i = 0;
        int K = time_allocator->K(i);

        for(int k = 0; k < K; k++)
        {   
            double delta_t; 
            double s_k, s_k_1, a_k, b_k, b_k_1;
            Vector3d velocity_s, acceleration_s, acceleration, acceleration_1, acceleration_2, jerk;

            if( k == 0 )
            {   
                if(i == 0)
                {
                    delta_t = time_allocator->time_acc(i, k);
                    s_k     = time_allocator->s(i, k);
                    s_k_1   = time_allocator->s(i, k + 1);
                    a_k     = time_allocator->a(i, k);
                    b_k     = time_allocator->b(i, k);
                    b_k_1   = time_allocator->b(i, k + 1);

                    velocity_s     = _bezier_basis->getVel(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i));
                    acceleration_s = _bezier_basis->getAcc(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i)) /time(i);
                    acceleration   = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                    jerk = acceleration / delta_t;

                    for(int p = 0; p < 3; p++)
                    {   
                        jerk_square_inte_i += pow(jerk(p),2) * delta_t;
                    }
                }
                else // do nothing
                {
                    jerk_square_inte_i += 0.0;
                }

            }
            else if (k == K-1)
            {   
                if(i == traj_seg_num - 1)
                {
                    delta_t = time_allocator->time(i, k) - time_allocator->time_acc(i, k);
                    s_k     = time_allocator->s(i, k);
                    s_k_1   = time_allocator->s(i, k + 1);
                    a_k     = time_allocator->a(i, k);
                    b_k     = time_allocator->b(i, k);
                    b_k_1   = time_allocator->b(i, k + 1);

                    velocity_s     = _bezier_basis->getVel(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i));
                    acceleration_s = _bezier_basis->getAcc(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i)) /time(i);
                    acceleration   = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                    jerk = acceleration / delta_t;
                    //cout<<"case 3"<<endl;
                }
                else
                {
                    delta_t = time_allocator->time(i, k) - time_allocator->time_acc(i, k) + time_allocator->time_acc(i + 1, 0);
                    s_k     = time_allocator->s(i, k);
                    s_k_1   = time_allocator->s(i, k + 1);
                    a_k     = time_allocator->a(i, k);
                    b_k     = time_allocator->b(i, k);
                    b_k_1   = time_allocator->b(i, k + 1);

                    velocity_s     = _bezier_basis->getVel(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i));
                    acceleration_s = _bezier_basis->getAcc(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i)) /time(i);
                    acceleration_1 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                    s_k     = time_allocator->s(i + 1, 0);
                    s_k_1   = time_allocator->s(i + 1, 1);
                    a_k     = time_allocator->a(i + 1, 0);
                    b_k     = time_allocator->b(i + 1, 0);
                    b_k_1   = time_allocator->b(i + 1, 1);

                    velocity_s     = _bezier_basis->getVel(bezier_coeff, i + 1, (s_k + s_k_1 ) / 2.0 /time(i + 1));
                    acceleration_s = _bezier_basis->getAcc(bezier_coeff, i + 1, (s_k + s_k_1 ) / 2.0 /time(i + 1)) /time(i + 1);
                    acceleration_2 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                    jerk = (acceleration_2 - acceleration_1) / delta_t;
                }

                for(int p = 0; p < 3; p++)
                {   
                    jerk_square_inte_i += pow(jerk(p),2) * delta_t;
                }
            }
            else
            {
                delta_t = time_allocator->time_acc(i, k + 1) - time_allocator->time_acc(i, k);

                s_k     = time_allocator->s(i, k);
                s_k_1   = time_allocator->s(i, k + 1);
                a_k     = time_allocator->a(i, k);
                b_k     = time_allocator->b(i, k);
                b_k_1   = time_allocator->b(i, k + 1);

                velocity_s     = _bezier_basis->getVel(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i));
                acceleration_s = _bezier_basis->getAcc(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i)) /time(i);
                acceleration_1 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                s_k     = time_allocator->s(i, k + 1);
                s_k_1   = time_allocator->s(i, k + 2);
                a_k     = time_allocator->a(i, k + 1);
                b_k     = time_allocator->b(i, k + 1);
                b_k_1   = time_allocator->b(i, k + 2);

                velocity_s     = _bezier_basis->getVel(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i));
                acceleration_s = _bezier_basis->getAcc(bezier_coeff, i, (s_k + s_k_1 ) / 2.0 /time(i)) /time(i);
                acceleration_2 = velocity_s * a_k + acceleration_s * (b_k + b_k_1) / 2.0;

                jerk = (acceleration_2 - acceleration_1) / delta_t;
                for(int p = 0; p < 3; p++)
                {   
                    jerk_square_inte_i += pow(jerk(p),2) * delta_t;
                }
            }
        }

        jerk_square_inte += jerk_square_inte_i;
    }

    return jerk_square_inte;
}

void UpdateTime(decomp_cvx_space::FlightCorridor & corridor, VectorXd bezier_time)
{
    corridor.durations.clear();
    for(int i = 0; i < (int)bezier_time.size(); i++)
        corridor.appendTime(bezier_time(i));
}   

void fastTrajPlanning(int _benchmark_poly_num){
    if( _has_map == false || _has_odom == false) 
        return;

    decomp_cvx_space::FlightCorridor corridor;
    if (_READ_CORRIDOR){
        int N_recd = corridor_recd.polyhedrons.size();
        ROS_WARN("read recorded %d corridors from msg, will take %d of them", N_recd, _benchmark_poly_num);

        if (_benchmark_poly_num > N_recd){
            ROS_WARN("no enough recorded polyhedrons");
            return;
        } else {
            for (int i = 0; i < _benchmark_poly_num; i++)
            {
                corridor.polyhedrons.push_back( corridor_recd.polyhedrons[i]);
            }
            _start_pt = corridor.polyhedrons[0].center;
            _end_pt = corridor.polyhedrons[_benchmark_poly_num-1].center;
        }
    } else {
        _start_pt = _manual_path.front();
        _end_pt   = _manual_path.back();
        corridor = _polyhedronGenerator->getCorridor();
    }

    ROS_WARN("[ddp_global_planner] Start ddp optimization");
    MatrixXd bezier_coeff;
    VectorXd bezier_time;
    
    /*  Initial time allocation  */
    initTimeAllocation(corridor);

    vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
    vector<double> durations = corridor.durations;
    int N = polyhedrons.size();
    init_allocTime = VectorXd::Zero(N);
    for (int i = 0; i < N; i++){
        init_allocTime(i) = durations[i];
    }
    MatrixXd pos = MatrixXd::Zero(2,3);
    MatrixXd vel = MatrixXd::Zero(2,3);
    MatrixXd acc = MatrixXd::Zero(2,3);
    MatrixXd jer = MatrixXd::Zero(2,3);

    pos.row(0) = _start_pt;
    pos.row(1) = _end_pt;    

    vel.row(0) << 0.0, 0.0, 0.0;
    acc.row(0) << 0.0, 0.0, 0.0;
    jer.row(0) << 0.0, 0.0, 0.0;
    
    double delta_s = 0.025;
    int iter_num = 0;
    int iter_max = _iter_max;
    double J, J_lst, J_best;
    J_lst = inf;

    temporalTrajOptimizer * time_optimizer  = new temporalTrajOptimizer();
    spatialTrajOptimizer  * curve_optimizer = new spatialTrajOptimizer();
    timeAllocator         * time_profile_   = new timeAllocator();
    ddpTrajOptimizer * ddp_optimizer = new ddpTrajOptimizer();
    ddpTrajOptimizer * ddp_optimizer0 = new ddpTrajOptimizer();

    time_optimizer->setParam(_traj_order, _bezier_basis->getC(), _bezier_basis->getC_v(), _bezier_basis->getC_a());
    time_optimizer->setType(1);

    VectorXd bezier_range_;
    MatrixXd bezier_coeff_;
    MatrixXd poly_coeff_;

    MatrixXd Qo_u = _bezier_basis->getMQM_u();
    MatrixXd Qo_l = _bezier_basis->getMQM_l();

    // create fake corridor for testing
    decomp_cvx_space::FlightCorridor corridor_fake = corridor;
    int N_fake = corridor_fake.polyhedrons.size();

    int num_cons = 0;
    for (int i = 0; i < N_fake; i++) {
        int num_plane = corridor_fake.polyhedrons[i].planes.size(); // hyperplane num of this polyhedra
        // constraint function
        for (int j = 0; j < num_plane; j++){
            corridor_fake.polyhedrons[i].planes[j].head(3) = Vector3d::Zero(3);
            corridor_fake.polyhedrons[i].planes[j](3) = -1.0; 
        }   
        num_cons += num_plane*6 + 5*3*2 + 4*3*2;     
    }
    cout << "num_cons = " << num_cons << endl;
    double terminalNorm = 0.0;
    MatrixXd initbezCoeff = MatrixXd::Zero(N, 3 * (_traj_order+1));
    
    Vector3d color(0.0,1.0,0.0);
    ddp_compTime = 0.0;
    bool infeas = true;
    bool line_failed = true;


    ROS_WARN("--------------line initialization ---------------------");

    ros::Time time1 = ros::Time::now();
    int rtn0;

    rtn0 = ddp_optimizer0->polyCurveGeneration(
                    corridor, Qo_u, Qo_l, pos, vel, acc, jer, _minimize_order, 
                    _MAX_Vel, _MAX_Acc, _MAX_Jer, initbezCoeff, _w_snap_zero, _w_terminal_zero, _w_time_zero, 1 * _iter_max_zero, infeas, true, false, line_failed, _time_power, false);


    ros::Time time2 = ros::Time::now();
    ROS_WARN("[ddp_global_planner] Time consumation of enforced ddp optimization is: %f",(time2 - time1).toSec() );
    int iter_used0;
    double jerkCost0;
    ddp_compTime += ddp_optimizer0->getCompTime();
    iter_used0 = ddp_optimizer0->getIterUsed();
    jerkCost0 = ddp_optimizer0->getJerkCost(); 
    bezier_coeff_ = ddp_optimizer0->getBezCoeff();
    bezier_range_ = ddp_optimizer0->getPolyTime(); 


    if (rtn0 == 2){
        UpdateTime(corridor, bezier_range_);
    } else {
        ROS_WARN("failed");
    }


    initbezCoeff = bezier_coeff_;
    int rtn = ddp_optimizer->polyCurveGeneration(
                corridor, Qo_u, Qo_l, pos, vel, acc, jer, _minimize_order, 
                _MAX_Vel, _MAX_Acc, _MAX_Jer, initbezCoeff, _w_snap, _w_terminal, _w_time, _iter_max, infeas, false, false, line_failed, _time_power, false); //todo
    ddp_compTime += ddp_optimizer->getCompTime();
    int iter_used = ddp_optimizer->getIterUsed();
    double jerkCost = ddp_optimizer->getJerkCost(); 
    terminalNorm = ddp_optimizer->getTerminalNorm();


    color << 0.0,0.0,1.0;
    bezier_coeff_ = ddp_optimizer->getBezCoeff(); 
    bezier_range_ = ddp_optimizer->getPolyTime(); 
    visFinalColoredBezierTrajectory(bezier_coeff_, bezier_range_, color); 
    ddp_allocTime = bezier_range_;
    bool allocTimePositive = true;
    for (int i=0; i<ddp_allocTime.size(); i++){
        if (ddp_allocTime(i) < 0.0){
            allocTimePositive = false;
        }
    }


    Results << ddp_compTime, ddp_allocTime.sum(), init_allocTime.sum(), rtn0, iter_used0, jerkCost0, rtn, iter_used, jerkCost, terminalNorm, 0.0;
    if (!allocTimePositive){
        Results(10) = -1.0;
    }

    
    delete curve_optimizer;
    delete time_optimizer;
    delete time_profile_;
    delete ddp_optimizer;
    delete ddp_optimizer0;
}


void trajPlanning(int _benchmark_poly_num)
{   
    if( _has_map == false || _has_odom == false) 
        return;

    decomp_cvx_space::FlightCorridor corridor;
    if (_READ_CORRIDOR){
        int N_recd = corridor_recd.polyhedrons.size();
        ROS_WARN("read recorded %d corridors from msg, will take %d of them", N_recd, _benchmark_poly_num);

        if (_benchmark_poly_num > N_recd){
            ROS_WARN("no enough recorded polyhedrons");
            return;
        } else {
            for (int i = 0; i < _benchmark_poly_num; i++)
            {
                corridor.polyhedrons.push_back( corridor_recd.polyhedrons[i]);
            }
            _start_pt = corridor.polyhedrons[0].center;
            _end_pt = corridor.polyhedrons[_benchmark_poly_num-1].center;
        }
    } else {
        _start_pt = _manual_path.front();
        _end_pt   = _manual_path.back();
        corridor = _polyhedronGenerator->getCorridor();
    }

    std::cout<<"start pt is "<<_start_pt<<"\n";
    std::cout<<"end pt is "<<_end_pt<<"\n";

    ROS_WARN("[trr_global_planner] Start coordinate descent-based spatial-temporal joint optimization");
    MatrixXd bezier_coeff;
    VectorXd bezier_time;
    
    /*  Initial time allocation  */
    initTimeAllocation(corridor);
    vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
    vector<double> durations = corridor.durations;
    int N = polyhedrons.size();
    init_allocTime = VectorXd::Zero(N);
    for (int i = 0; i < N; i++){
        init_allocTime(i) = durations[i];
    }

    MatrixXd pos = MatrixXd::Zero(2,3);
    MatrixXd vel = MatrixXd::Zero(2,3);
    MatrixXd acc = MatrixXd::Zero(2,3);

    pos.row(0) = _start_pt;
    pos.row(1) = _end_pt;    

    vel.row(0) << 0.0, 0.0, 0.0;
    acc.row(0) << 0.0, 0.0, 0.0;

    double delta_s = 0.025;
    int iter_num = 0;
    int iter_max = _iter_max;
    double J, J_lst, J_best;
    J_lst = inf;

    temporalTrajOptimizer * time_optimizer  = new temporalTrajOptimizer();
    spatialTrajOptimizer  * curve_optimizer = new spatialTrajOptimizer();
    timeAllocator         * time_profile_   = new timeAllocator();

    time_optimizer->setParam(_traj_order, _bezier_basis->getC(), _bezier_basis->getC_v(), _bezier_basis->getC_a());
    time_optimizer->setType(1);

    VectorXd bezier_range_;
    MatrixXd bezier_coeff_;

    MatrixXd Qo_u = _bezier_basis->getMQM_u();
    MatrixXd Qo_l = _bezier_basis->getMQM_l();
    trr_compTime = 0.0;
    double energy = 0.0;
    while(iter_num < iter_max)
    {
        ros::Time time_before_spatial_optimization = ros::Time::now();
        
        int error_code = 
        curve_optimizer->bezierCurveGeneration(
                corridor, Qo_u, Qo_l, pos, vel, acc, _traj_order, _minimize_order, _MAX_Vel, _MAX_Acc);

        ros::Time time_after_spatial_optimization = ros::Time::now();
        // ROS_WARN("[trr_global_planner] Time consumation of spatial optimization is: %f",(time_after_optimization - time_before_optimization).toSec() );
        trr_compTime += (time_after_spatial_optimization - time_before_spatial_optimization).toSec();

        bezier_coeff = curve_optimizer->getPolyCoeff();
        bezier_time  = curve_optimizer->getPolyTime();

        VectorXd bezier_range = bezier_time;

        if( error_code != 0 ){
            trr_allocTime = bezier_range_;
            Results << trr_compTime, trr_allocTime.sum(), init_allocTime.sum(), error_code, 0.0, 
                        0.0, 0.0, 0.0, 0.0, 0.0, 1.0 * iter_num;
            // ROS_WARN("[trr_global_planner] Cannot find a safe solution, somthing wrong with the solver");  
            return;
        }
        else{   
            ros::Time time_before_time_optimization = ros::Time::now();
            time_optimizer->timeProfileGeneration( bezier_coeff, bezier_time, _MAX_Vel, _MAX_Acc, _MAX_d_Acc, delta_s, 0.0 ); // 0.0 for minimizing the time, no regulizer in control cost
            // bezier_time already returns the optimized time
            ros::Time time_after_time_optimization = ros::Time::now();
            trr_compTime += (time_after_time_optimization - time_before_time_optimization).toSec();

            timeAllocator * time_profile = time_optimizer->getTimeProfile();
        
            energy = getEnergy(time_profile, bezier_coeff, bezier_range);
            J = energy + _rho * bezier_time.sum();
            
            // cout<<"Energy cost is: "<<energy<<", Time duration is: "<<bezier_time.sum()<<endl;
            if((pow((J - J_lst), 2) < J_lst * 0.01) || J > J_lst)
                break;
            else{
                J_lst  = J;
                J_best = J;
                bezier_coeff_ = bezier_coeff;
                bezier_range_ = bezier_range;
                time_profile_ = time_profile;
            }
            
            // visBezierTrajectory(bezier_coeff, bezier_range, iter_num);
            UpdateTime(corridor, bezier_time);
        }

        iter_num ++;
    }

    // ROS_WARN("[trr_global_planner] find a solution after iteration: %d, the optimal cost is: %f", iter_num, J_best);

    // visFinalBezierTrajectory(bezier_coeff_, bezier_range_); 

    _traj_time_start = _odom_time + ( ros::Time::now() - _time_update_odom);
    _traj_time_final = _traj_time_start;

    for(int i = 0; i < time_profile_->time.rows(); i++){   
        int K = time_profile_->K(i); // K is number of discretized intervel in each segment
        _traj_time_final += ros::Duration(time_profile_->time(i, K - 1)); 
    }

    quadrotor_msgs::SpatialTemporalTrajectory space_time_msgs 
        = getSpaceTimeTraj(time_profile_, bezier_coeff_, bezier_range_);

    // _space_time_pub.publish(space_time_msgs);

    _traj_id ++;
    _has_traj = true;
    
    cout<<"publish global planning"<<endl;
    trr_allocTime = bezier_range_;

    Results << trr_compTime, trr_allocTime.sum(), init_allocTime.sum(), 0.0, energy, 
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0*iter_num;
    
    delete curve_optimizer;
    delete time_optimizer;
    delete time_profile_;
}




int main(int argc, char** argv)
{
    ros::init(argc, argv, "trr_global_planner_node");
    ros::NodeHandle nh("~");

    nh.param("write_path", _WRITE_PATH, false);
    nh.param("read_path",  _READ_PATH,  false);

    nh.param("write_corridor", _WRITE_CORRIDOR, false);
    nh.param("read_corridor",  _READ_CORRIDOR,  false);
    nh.param("no_cons",  _no_cons,  false);

    nh.param("map/map_margin", _cloud_margin, 0.25);
    nh.param("map/resolution", _resolution,   0.2 );
    
    nh.param("map/x_size",         _x_size, 50.0);
    nh.param("map/y_size",         _y_size, 50.0);
    nh.param("map/z_size",         _z_size, 5.0 );
    
    nh.param("planning/rho_time",  _rho,       1.0);
    nh.param("planning/max_vel",   _MAX_Vel,   1.0);
    nh.param("planning/max_acc",   _MAX_Acc,   1.0);
    nh.param("planning/max_jer",   _MAX_Jer,   5.0);
    nh.param("planning/max_d_acc", _MAX_d_Acc, 1.0);
    nh.param("planning/w_Q",    _w_Q,    1.0e0);
    nh.param("planning/w_snap",    _w_snap,    1.0e-5);
    nh.param("planning/w_terminal",_w_terminal,1.0e-3);
    nh.param("planning/w_time",    _w_time,    1.0e-2);
    nh.param("planning/w_terminal_line",_w_terminal_line,1.0e-3);
    nh.param("planning/w_time_line",    _w_time_line,    1.0e-2);
    nh.param("planning/w_snap_zero",    _w_snap_zero,    1.0e-5);
    nh.param("planning/w_terminal_zero",_w_terminal_zero,1.0e0);
    nh.param("planning/w_time_zero",    _w_time_zero,    1.0e0);
    nh.param("planning/time_power",    _time_power,    2);

    nh.param("planning/w_terminal_bez",_w_terminal_bez,1.0e-3);
    nh.param("planning/w_time_bez",    _w_time_bez,    1.0e-2);
    nh.param("planning/w_snap_bez",    _w_snap_bez,    1.0e-2);

    nh.param("planning/alg_id",    _alg_id,    0);
    nh.param("planning/segNum",    _segNum,    10);
    nh.param("planning/fast_test",    _fast_test,    20);
    nh.param("planning/relaxtol",    _relaxtol,    0.0);

    nh.param("planning/max_inf_iter", _max_inf_iter,        10  );
    nh.param("planning/max_clu_iter", _max_clu_iter,        50  );
    nh.param("planning/iter_max",         _iter_max,        50  );
    nh.param("planning/iter_max_zero",    _iter_max_zero,        50  );

    // nh.param("planning/benchmark_poly_num",  _benchmark_poly_num,  3);

    nh.param("optimization/min_order",  _minimize_order,    3.0 );
    nh.param("optimization/poly_order", _traj_order,        10  );

    nh.param("vis/vis_traj_width",     _vis_traj_width,     0.15);

    nh.param("your_file_path",         _your_file_path,     string("")  );

    _map_sub  = nh.subscribe( "map",       1, rcvPointCloudCallBack );
    _odom_sub = nh.subscribe( "odometry",  1, rcvOdometryCallBack);
    _joy_sub  = nh.subscribe( "joystick",  1, rcvJoyCallBack );
    _target_sub  = nh.subscribe( "/move_base_simple/goal",  1, targetCallBack );
    _start_sub  = nh.subscribe( "/clicked_point",  1, startTrajCallBack );
    _corridor_rec_sub  = nh.subscribe("corridor", 20, corridorRecCallBack); 

    // for visualization of the planning results
    _vis_traj_pub      = nh.advertise<visualization_msgs::Marker>("trajectory_vis", 1);    
    _vis_grid_path_pub = nh.advertise<visualization_msgs::MarkerArray>("grid_path_vis", 1);
    _vis_inf_map_pub   = nh.advertise<sensor_msgs::PointCloud2>("inflation_map", 10);
    _vis_polytope_pub  = nh.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedron_corridor_mesh", 1, true);
    _space_time_pub    = nh.advertise<quadrotor_msgs::SpatialTemporalTrajectory>("space_time_traj", 10); 
    _corridor_rec_pub  = nh.advertise<msgs::corridor>("corridor", 1); 

    _map_lower << -_x_size/2.0, -_y_size/2.0, 0.0;
    _map_upper << +_x_size/2.0, +_y_size/2.0, _z_size;

    _poly_array_msg.header.frame_id = "/map";
    _vis_polytope_pub.publish(_poly_array_msg); 

    _pt_max_x = + _x_size / 2.0; _pt_min_x = - _x_size / 2.0;
    _pt_max_y = + _y_size / 2.0; _pt_min_y = - _y_size / 2.0; 
    _pt_max_z = + _z_size;       _pt_min_z = 0.0;
    
    _resolution = max(0.1, _resolution); // In case a too fine resolution crashes the CUDA code.
    _inv_resolution = 1.0 / _resolution;
    _max_x_id = (int)(_x_size * _inv_resolution);
    _max_y_id = (int)(_y_size * _inv_resolution);
    _max_z_id = (int)(_z_size * _inv_resolution);

    _bezier_basis = new Bernstein(_minimize_order);
    _bezier_basis ->setFixedOrder(_traj_order);

    _path_finder  = new gridPathFinder(_max_x_id, _max_y_id, _max_z_id);
    _path_finder  -> initGridMap(_resolution, _map_lower, _map_upper);
    
    _polyhedronGenerator->initialize(false, true, true, // method in generating corridor
        _max_x_id, _max_y_id, _max_z_id, _map_lower, _map_upper, _resolution, _inv_resolution, // map information
        _max_inf_iter, _max_clu_iter); // max infaltion/clustering num

    if(_WRITE_PATH){
        path_record.open(_your_file_path);
    }

    ros::Rate rate(100);
    bool status = ros::ok();
    while(status){
        ros::spinOnce();  
        status = ros::ok();
        rate.sleep();
    }

    delete _bezier_basis;
    delete _path_finder;

    return 0;
}

quadrotor_msgs::SpatialTemporalTrajectory getSpaceTimeTraj(const timeAllocator * time_allocator, const MatrixXd & bezier_coeff, const VectorXd & range )
{
    quadrotor_msgs::SpatialTemporalTrajectory space_time_traj;
    space_time_traj.action = quadrotor_msgs::OptimalTimeAllocator::ACTION_ADD;
    space_time_traj.header.frame_id = "/trr_global_planner";
    space_time_traj.trajectory_id = _traj_id;

    space_time_traj.header.stamp = _traj_time_start; 
    space_time_traj.start_time   = _traj_time_start; 
    space_time_traj.final_time   = _traj_time_final; 

    int seg_num = time_allocator->seg_num;
    space_time_traj.s_step = time_allocator->s_step; 
    space_time_traj.K_max  = time_allocator->K_max;

    for (int i = 0; i < seg_num; i++ )
    {     
        space_time_traj.K.push_back(time_allocator->K(i));

        for(int j = 0; j < time_allocator->K(i) + 1; j++ )
        {   
            if( j < time_allocator->K(i) )
            {
                space_time_traj.a.push_back( time_allocator->a(i, j) );
                space_time_traj.time.push_back    ( time_allocator->time(i, j) );
                space_time_traj.time_acc.push_back( time_allocator->time_acc(i, j) );
            }

            space_time_traj.b.push_back( time_allocator->b(i, j) );
            space_time_traj.s.push_back( time_allocator->s(i, j) );
        }
    }

    /*  stack the spatial curve  */
    seg_num = range.size();
    space_time_traj.num_segment = seg_num;

    for(int i = 0; i < seg_num; i++ )
    {    
        int poly_num1d = _traj_order + 1;
        for(int j =0; j < poly_num1d; j++)
        { 
            space_time_traj.coef_x.push_back(bezier_coeff(i,                  j));
            space_time_traj.coef_y.push_back(bezier_coeff(i,     poly_num1d + j));
            space_time_traj.coef_z.push_back(bezier_coeff(i, 2 * poly_num1d + j));
        }
        space_time_traj.range.push_back(range(i));
    }

    space_time_traj.start_yaw = 0.0;
    space_time_traj.final_yaw = 0.0;
    space_time_traj.trajectory_id = _traj_id;

    return space_time_traj;
}

geometry_msgs::Point Vector2Point(Vector3d vec)
{
    geometry_msgs::Point pt;
    pt.x = vec(0);   
    pt.y = vec(1);   
    pt.z = vec(2);   

    return pt;
}

void clearVisualization(int grid_path_mk_num, int traj_mk_iter_num)
{
    // 1. Clear the Grid Path :
    visualization_msgs::MarkerArray grid_vis; 
    visualization_msgs::Marker mk;
    mk.header.frame_id = "map";
    mk.header.stamp = ros::Time::now();
    mk.ns = "trr_global_planner/grid_path";
    mk.type = visualization_msgs::Marker::CUBE;
    mk.action = visualization_msgs::Marker::DELETE;
    
    for(int i = 0; i < grid_path_mk_num; i++)
    {
        mk.id = i;
        grid_vis.markers.push_back(mk);
    }

    _vis_grid_path_pub.publish(grid_vis);

    // 2. Clear the Bezier Curves
    cout<<"bezier curves visualization num: "<<_vis_iter_num<<endl;
    for(int i = 0; i < _vis_iter_num; i++)
    {
        visualization_msgs::Marker traj_vis;

        traj_vis.header.stamp       = ros::Time::now();
        traj_vis.header.frame_id    = "map";
        traj_vis.ns = "trr_global_planner/trajectory" + to_string(i);
        traj_vis.id = 0;
        traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
        traj_vis.action = visualization_msgs::Marker::DELETE;

        for(int k = 0; k < 100; k++)
            _vis_traj_pub.publish(traj_vis);
    }

    // 3. Clear the polyhedrons
    _poly_array_msg.polyhedrons.clear();
    _poly_array_msg.ids.clear();
    _vis_polytope_pub.publish(_poly_array_msg);
}

void visCorridor( decomp_ros_msgs::PolyhedronArray & poly_array_msg )
{
    _vis_polytope_pub.publish(poly_array_msg);
}

void visBezierTrajectory(MatrixXd polyCoeff, VectorXd time, int iter_id)
{   
    visualization_msgs::Marker traj_vis;

    _vis_iter_num ++;
    traj_vis.header.stamp       = ros::Time::now();
    traj_vis.header.frame_id    = "map";

    traj_vis.ns = "trr_global_planner/trajectory" + to_string(iter_id);
    traj_vis.id = 0;
    traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;

    traj_vis.action = visualization_msgs::Marker::ADD;
    traj_vis.scale.x = _vis_traj_width;
    traj_vis.scale.y = _vis_traj_width;
    traj_vis.scale.z = _vis_traj_width;
    traj_vis.pose.orientation.x = 0.0;
    traj_vis.pose.orientation.y = 0.0;
    traj_vis.pose.orientation.z = 0.0;
    traj_vis.pose.orientation.w = 1.0;

    traj_vis.color.a = 1.0;
    traj_vis.color.r = min(1.0, (iter_id + 1) / 5.0);
    traj_vis.color.g = 0.0;//max(0.0, 1 - rgb / 5.0);
    traj_vis.color.b = 0.0;

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    traj_vis.points.clear();

    Vector3d state;
    geometry_msgs::Point pt;

    int segment_num  = polyCoeff.rows();
    for(int i = 0; i < segment_num; i++ ){
        for (double t = 0.0; t < 1.0; t += 0.1 / time(i), count += 1){
            state = _bezier_basis->getPosFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = time(i) * state(0);
            cur(1) = pt.y = time(i) * state(1);
            cur(2) = pt.z = time(i) * state(2);
            traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    _vis_traj_pub.publish(traj_vis);
}

void visColoredBezierTrajectory(MatrixXd polyCoeff, VectorXd time, int iter_id, Vector3d color)
{   
    visualization_msgs::Marker traj_vis;

    _vis_iter_num ++;
    traj_vis.header.stamp       = ros::Time::now();
    traj_vis.header.frame_id    = "map";

    traj_vis.ns = "trr_global_planner/ddptrajectory" + to_string(iter_id);
    traj_vis.id = 0;
    traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;

    traj_vis.action = visualization_msgs::Marker::ADD;
    traj_vis.scale.x = _vis_traj_width;
    traj_vis.scale.y = _vis_traj_width;
    traj_vis.scale.z = _vis_traj_width;
    traj_vis.pose.orientation.x = 0.0;
    traj_vis.pose.orientation.y = 0.0;
    traj_vis.pose.orientation.z = 0.0;
    traj_vis.pose.orientation.w = 1.0;

    traj_vis.color.a = 1.0;
    // traj_vis.color.r = min(1.0, (iter_id + 1) / 5.0);
    // traj_vis.color.g = 0.0;//max(0.0, 1 - rgb / 5.0);
    // traj_vis.color.b = 0.0;
    traj_vis.color.r = color(0);
    traj_vis.color.g = color(1);
    traj_vis.color.b = color(2);
    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    traj_vis.points.clear();

    Vector3d state;
    geometry_msgs::Point pt;

    int segment_num  = polyCoeff.rows();
    for(int i = 0; i < segment_num; i++ ){
        if(time(i) < 0){
            ROS_WARN("time less than 0, returning");
            return;
        }
        for (double t = 0.0; t < 1.0; t += 0.1 / time(i), count += 1){
            state = _bezier_basis->getPosFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = time(i) * state(0);
            cur(1) = pt.y = time(i) * state(1);
            cur(2) = pt.z = time(i) * state(2);
            traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    _vis_traj_pub.publish(traj_vis);
}

void visFinalBezierTrajectory(MatrixXd polyCoeff, VectorXd time)
{   
    visualization_msgs::Marker traj_vis;

    traj_vis.header.stamp       = ros::Time::now();
    traj_vis.header.frame_id    = "map";

    traj_vis.ns = "trr_global_planner/trajectory_final";
    traj_vis.id = 0;
    traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    
    traj_vis.action = visualization_msgs::Marker::ADD;
    traj_vis.scale.x = _vis_traj_width;
    traj_vis.scale.y = _vis_traj_width;
    traj_vis.scale.z = _vis_traj_width;
    traj_vis.pose.orientation.x = 0.0;
    traj_vis.pose.orientation.y = 0.0;
    traj_vis.pose.orientation.z = 0.0;
    traj_vis.pose.orientation.w = 1.0;

    traj_vis.color.a = 1.0;
    traj_vis.color.r = 1.0;
    traj_vis.color.g = 0.0;
    traj_vis.color.b = 0.0;

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    traj_vis.points.clear();

    Vector3d state;
    geometry_msgs::Point pt;

    int segment_num  = polyCoeff.rows();
    for(int i = 0; i < segment_num; i++ ){
        for (double t = 0.0; t < 1.0; t += 0.01 / time(i), count += 1){
            state = _bezier_basis->getPosFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = time(i) * state(0);
            cur(1) = pt.y = time(i) * state(1);
            cur(2) = pt.z = time(i) * state(2);
            traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    _vis_traj_pub.publish(traj_vis);
}

void visFinalColoredBezierTrajectory(MatrixXd polyCoeff, VectorXd time, Vector3d color)
{   
    visualization_msgs::Marker traj_vis;

    traj_vis.header.stamp       = ros::Time::now();
    traj_vis.header.frame_id    = "map";

    traj_vis.ns = "trr_global_planner/ddptrajectory_final";
    traj_vis.id = 0;
    traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    
    traj_vis.action = visualization_msgs::Marker::ADD;
    traj_vis.scale.x = _vis_traj_width;
    traj_vis.scale.y = _vis_traj_width;
    traj_vis.scale.z = _vis_traj_width;
    traj_vis.pose.orientation.x = 0.0;
    traj_vis.pose.orientation.y = 0.0;
    traj_vis.pose.orientation.z = 0.0;
    traj_vis.pose.orientation.w = 1.0;

    traj_vis.color.a = 1.0;
    traj_vis.color.r = color(0);
    traj_vis.color.g = color(1);
    traj_vis.color.b = color(2);

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    traj_vis.points.clear();

    Vector3d state;
    geometry_msgs::Point pt;

    int segment_num  = polyCoeff.rows();
    for(int i = 0; i < segment_num; i++ ){
        if(time(i) < 0){
            ROS_WARN("time less than 0, returning");
            return;
        }
        for (double t = 0.0; t < 1.0; t += 0.2 / time(i), count += 1){
            state = _bezier_basis->getPosFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = time(i) * state(0);
            cur(1) = pt.y = time(i) * state(1);
            cur(2) = pt.z = time(i) * state(2);
            traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    _vis_traj_pub.publish(traj_vis);
}

void visGridPath( )
{   
    visualization_msgs::MarkerArray grid_vis; 
    visualization_msgs::Marker mk;
    mk.header.frame_id = "map";
    mk.header.stamp = ros::Time::now();
    mk.ns = "trr_global_planner/grid_path";
    mk.type = visualization_msgs::Marker::CUBE;
    mk.action = visualization_msgs::Marker::ADD;

    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;

    mk.color.a = 1.0;
    mk.color.r = 1.0;
    mk.color.g = 1.0;
    mk.color.b = 0.0;
    
    int idx = 0;
    for(int i = 0; i < int(_manual_path.size()); i++)
    {
        mk.id = idx;
        mk.pose.position.x = _manual_path[i](0); 
        mk.pose.position.y = _manual_path[i](1); 
        mk.pose.position.z = _manual_path[i](2);  

        mk.scale.x = _resolution;
        mk.scale.y = _resolution;
        mk.scale.z = _resolution;

        idx ++;
        grid_vis.markers.push_back(mk);
    }

    _vis_grid_path_pub.publish(grid_vis);
}