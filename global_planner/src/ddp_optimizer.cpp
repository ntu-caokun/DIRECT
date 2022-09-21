#include <global_planner/ddp_optimizer.h>
using namespace std;    
using namespace Eigen;

int ddpTrajOptimizer::polyCurveGeneration( 
        const decomp_cvx_space::FlightCorridor &corridor,
        const MatrixXd &MQM_u,
        const MatrixXd &MQM_l,
        const MatrixXd &pos,
        const MatrixXd &vel,
        const MatrixXd &acc,
        const MatrixXd &jer,
        const double minimize_order,
        const double maxVel, 
        const double maxAcc,
        const double maxJer,
        Eigen::MatrixXd initbezCoeff,
        const double w_snap,
        const double w_terminal,
        const double w_time,
        const int iter_max,
        bool &infeas,
        bool zero_init_flag, 
        bool line_init_flag,
        bool &line_failed,
        int time_power,
        bool minvo_flag) 
{   

    ros::Time time_before_optimization = ros::Time::now();

    // ROS_WARN("~~~entered polyCurvegeneration~~~");
    int rtn = 0;
    vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
    vector<double> durations = corridor.durations;
    int N = polyhedrons.size();
    int traj_order = 5;
    int dim = 3;
    int sys_order = (traj_order+1)/2;
    int num_ctrlP = traj_order + 1;

    alg.maxiter = iter_max;
    alg.tol = 1.0e-7;
    alg.infeas = infeas;

    vector<double> steps;

    // fwdPass fp;
    fp.N = N;
    fp.traj_order = traj_order;
    fp.dim = dim;
    fp.sys_order = sys_order;
    fp.num_ctrlP = num_ctrlP;
    fp.w_snap = w_snap;
    fp.Rtime = w_time;
    fp.time_power = time_power;
    fp.maxVel = maxVel;
    fp.maxAcc = maxAcc;
    fp.minvo_enabled = minvo_flag;
    fp.reg_exp_base = 1.6;
    if (!zero_init_flag) fp.reg_exp_base = 4.0;
    if (fp.minvo_enabled){
        fp.polyt2minvotau6_ << 1.0, -0.06471861202, -0.03728008486, -0.02577637794, -0.02027573243,  -0.01678273037,
                            1.0,  0.03314986096, -0.06548114211, -0.05530463802, -0.04362718953,  -0.03671639115,
                            1.0,   0.3375528997,  0.05836232552, -0.02920033165, -0.04690387913,  -0.04376447947,
                            1.0,   0.6624471003,   0.3832565261,   0.1916286091,  0.06985980172, -0.002892843108,
                            1.0,    0.966850139,    0.868219136,   0.7594116288,   0.6521050661,    0.5510660979,
                            1.0,    1.064718612,    1.092157139,    1.108091959,    1.118023718,     1.123960059;     
        fp.polyt2minvotau_v6_ <<0, 1.0, -0.1423379297, -0.1332742327, -0.1242105357,  -0.126304257,
                            0, 1.0,  0.1887439858, -0.1831318297, -0.2466606848, -0.2321393311,
                            0, 1.0,           1.0,  0.5411016575, 0.08220331498, -0.2433474658,
                            0, 1.0,   1.811256014,   2.250636213,   2.381669451,   2.282405938,
                            0, 1.0,    2.14233793,   3.293739556,   4.445141183,   5.585385392;      
        fp.polyt2minvotau_a6_ <<0, 0, 2.0, -0.4472869252, -0.6133793313, -0.6406553622,
                            0, 0, 2.0,   1.223711659, -0.5552714346,  -1.854618819,
                            0, 0, 2.0,   4.776288341,    6.54988193,   6.841145057,
                            0, 0, 2.0,   6.447286925,   13.17576837,   22.04662796;    
    } else {
        fp.polyt2minvotau6_ <<      1.0,   0,   0,   0,   0,   0,
                                1.0, 0.2,   0,   0,   0,   0,
                                1.0, 0.4, 0.1,   0,   0,   0,
                                1.0, 0.6, 0.3, 0.1,   0,   0,
                                1.0, 0.8, 0.6, 0.4, 0.2,   0,
                                1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

        fp.polyt2minvotau_v6_  <<      0, 1.0,   0,   0,   0,   0,
                                0, 1.0, 0.5,   0,   0,   0,
                                0, 1.0, 1.0, 0.5,   0,   0,
                                0, 1.0, 1.5, 1.5, 1.0,   0,
                                0, 1.0, 2.0, 3.0, 4.0, 5.0;

        fp.polyt2minvotau_a6_ <<      0, 0, 2.0,   0,    0,    0,
                                0, 0, 2.0, 2.0,    0,    0,
                                0, 0, 2.0, 4.0,  4.0,    0,
                                0, 0, 2.0, 6.0, 12.0, 20.0;      
    }   

    if (traj_order == 7){
        fp.Ek_inv = {1,1,1/2.0,1/6.0};
    }
    if (traj_order == 5){
        fp.Ek_inv = {1,1,1/2.0};
    }
    VectorXd xf = VectorXd::Zero(dim*sys_order);
    if (sys_order == 4){
        xf << pos.row(1).transpose(), vel.row(1).transpose(), acc.row(1).transpose(), jer.row(1).transpose();
    }
    if (sys_order == 3){
        xf << pos.row(1).transpose(), vel.row(1).transpose(), acc.row(1).transpose();
    }
    fp.x_d = xf;
    fp.setPmat(w_terminal*MatrixXd::Identity(sys_order*dim, sys_order*dim)); // default is identity

    // initialization
    VectorXd x0 = VectorXd::Zero(dim*sys_order);
    if (sys_order == 4){
        x0 << pos.row(0).transpose(), vel.row(0).transpose(), acc.row(0).transpose(), jer.row(0).transpose();
    }
    if (sys_order == 3){
        x0 << pos.row(0).transpose(), vel.row(0).transpose(), acc.row(0).transpose();
    }
    fp.x.emplace_back(x0);

    for (int i = 0; i < N; i++){
        fp.x.emplace_back( VectorXd::Zero(dim*sys_order) );
        VectorXd ui = VectorXd::Zero(dim*sys_order+1); 
        ui.tail(1) << durations[i];
        fp.u.emplace_back(ui);
        decomp_cvx_space::Polytope pltp = polyhedrons[i];
        int num_plane = pltp.planes.size(); // hyperplane num of this polyhedra

        fp.fx.emplace_back(MatrixXd::Zero(sys_order*dim, sys_order*dim));
        fp.fu.emplace_back(MatrixXd::Zero(sys_order*dim, sys_order*dim+1));
        fp.fxx.emplace_back(MatrixXd::Zero(sys_order*dim * sys_order*dim, sys_order*dim)); // write tensor as matrix
        fp.fxu.emplace_back(MatrixXd::Zero(sys_order*dim * (sys_order*dim+1), sys_order*dim));
        fp.fuu.emplace_back(MatrixXd::Zero( sys_order*dim * (sys_order*dim+1), sys_order*dim+1));

        // fp.q is initialized in initialroll
        fp.qx.emplace_back(MatrixXd::Zero(sys_order*dim, 1));
        fp.qu.emplace_back(MatrixXd::Zero(sys_order*dim+1, 1));
        fp.qxx.emplace_back( MatrixXd::Zero(sys_order*dim, sys_order*dim) );
        fp.qxu.emplace_back( MatrixXd::Zero(sys_order*dim, sys_order*dim+1) );
        fp.quu.emplace_back( MatrixXd::Zero(sys_order*dim+1, sys_order*dim+1) );

        int num_cons = num_plane*num_ctrlP + (num_ctrlP-1)*dim*2 + (num_ctrlP-2)*dim*2 + 1;
        fp.c.emplace_back(MatrixXd::Zero(num_cons, 1)); 
        fp.cx.emplace_back(MatrixXd::Zero(num_cons, sys_order*dim)); 
        fp.cu.emplace_back(MatrixXd::Zero(num_cons, sys_order*dim+1)); 

        fp.s.emplace_back( 1.0e-1 * VectorXd::Ones(num_cons) );
        fp.y.emplace_back( 0.01 * VectorXd::Ones(num_cons));
        fp.mu.emplace_back( fp.s[i].array() * fp.y[i].array() );

        bp.ky.emplace_back( VectorXd::Zero(num_cons) );
        bp.Ky.emplace_back( MatrixXd::Zero(num_cons, sys_order*dim) );
        bp.ks.emplace_back( VectorXd::Zero(num_cons) );
        bp.Ks.emplace_back( MatrixXd::Zero(num_cons, sys_order*dim) );
        bp.ku.emplace_back( VectorXd::Zero(sys_order*dim+1) );
        bp.Ku.emplace_back( MatrixXd::Zero(sys_order*dim+1, sys_order*dim) );
    }


    PolyTime = VectorXd::Zero(N);
    for(int i = 0; i < N; i++ )
    {   PolyTime(i) = durations[i]; }

    BezCoeff = initbezCoeff; // this is N * (num_ctrlP * dim), each row is [x num_ctrlP, y num_ctrlP, z num_ctrlP]

    // convert to dim * num_ctrlP, i.e., [x y z]*8
    for (int i = 0; i < N; i++){
        VectorXd tempv = BezCoeff.row(i); 
        Map<MatrixXd> BezCoeffi_mat(tempv.data(), fp.num_ctrlP, fp.dim);
        MatrixXd tempm = BezCoeffi_mat.transpose();
        Map<RowVectorXd> BezCoeffi(tempm.data(), BezCoeffi_mat.size()); 
        BezCoeff.row(i) =  BezCoeffi;
    }

    VectorXd initPolyTime = PolyTime;
    MatrixXd initpolyCoeff = bez2polyFunc();

    fp.barEk_inv = VectorXd::Zero(sys_order*dim);
    for (int i = 0; i < sys_order; i++){
        for (int j = 0; j < dim; j++){
            fp.barEk_inv(i*dim + j) = fp.Ek_inv[i];
        }
    }    

    // initialize fp.u
    if (!zero_init_flag){
        if (!line_init_flag){  // init from passed values
            for (int i = 0; i < N; i++){
                fp.u[i].head(dim*sys_order) = initpolyCoeff.row(i).tail(dim*sys_order);
            }
        } else {
            vector<Vector3d> points;
            points.push_back (pos.row(0).transpose());
            for(int i = 1; i < (int)corridor.polyhedrons.size(); i++) // from 1
                points.push_back(corridor.polyhedrons[i].seed_coord);
            Vector3d endpoint_pert = pos.row(1).transpose();
            points.push_back (endpoint_pert);
            cout << "enter line init"<< endl;

            for (int l = 0; l < N; l++){
                bool constraint_vio_flag = true;
                int increase_time_count = 0;
                while (constraint_vio_flag && increase_time_count <= 4){
                    double Tk = (fp.u[l].tail(1))(0);
                    double Tk2 = Tk * Tk;
                    double Tk3 = Tk2 * Tk;
                    double Tk4 = Tk3 * Tk;
                    double Tk5 = Tk4 * Tk;
                    double Tk6 = Tk5 * Tk;
                    double Tk7 = Tk6 * Tk;
                    MatrixXd tempFk(sys_order,sys_order);
                    MatrixXd tempGkinv(sys_order,sys_order); 
                    tempFk <<   1.0,   Tk, Tk2/2.0, 
                                0.0,  1.0,      Tk, 
                                0.0,  0.0,     1.0;   
                    tempGkinv << 10.0/Tk3, -4.0/Tk2,   0.5/Tk,
                                -15.0/Tk4,  7.0/Tk3, -1.0/Tk2,
                                6.0/Tk5, -3.0/Tk4,  0.5/Tk3;  
                    MatrixXd tempbarFk=MatrixXd::Zero(dim * sys_order, dim * sys_order);
                    MatrixXd tempbarGkinv=MatrixXd::Zero(dim * sys_order, dim * sys_order);
                    for (int i = 0; i < sys_order; i++){
                        for (int j = 0; j < sys_order; j++){
                            for (int k = 0; k < dim; k++){
                                tempbarFk(i*dim+k, j*dim+k) = tempFk(i,j);
                                tempbarGkinv(i*dim+k, j*dim+k) = tempGkinv(i,j);
                            }
                        }
                    }          
                    VectorXd xnext = VectorXd::Zero(dim*sys_order);
                    VectorXd xcur = VectorXd::Zero(dim*sys_order);
                    xnext.head(dim) = points[l+1];
                    xcur.head(dim) = points[l];
                    fp.u[l].head(dim*sys_order) = tempbarGkinv * (xnext - tempbarFk * xcur);
                    VectorXd cons = fp.computecminvo(xcur, fp.u[l], l, corridor);

                    if ((cons.array() < 0).all()) {
                        constraint_vio_flag = false;
                    } else {
                        fp.u[l].tail(1) << 2 * Tk;
                        cout << "time:"<<l<<"increase time to" << 2*Tk << endl;
                        increase_time_count++;
                    }
                }
            }
        }

    }
    // ROS_WARN("initialization finished");
    fp.initialroll(corridor);
    double initCost = fp.cost;

    if (line_init_flag){ // for line initilization, will reset infeas
        int count = 0;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < fp.c[i].size(); j++){
                if (fp.c[i](j) > 0){
                    count++;
                }
            }
        }
        if (count == 0){
            alg.infeas = false;
            ROS_WARN("line_init satisfies constraints, thus reset alg.infeas to false");
        }
        cout << "current alginfeas="<<alg.infeas<<endl; 
    }

    PolyCoeff = MatrixXd::Zero(N, num_ctrlP * dim );
    PolyTime  = VectorXd::Zero(N);

    sysparam2polyFunc();
    PolyCoeffTraj.emplace_back(PolyCoeff);
    PolyTimeTraj.emplace_back(PolyTime);
    costTraj.emplace_back(fp.cost);
    costqTraj.emplace_back(fp.costq);

 
    alg.mu = fp.cost/fp.N/fp.s[0].size();
    fp.resetfilter(alg);
    bp.resetreg();
    if (line_init_flag){
        bp.initreg(10.0);
    }

    compTime = 0.0;
    int iter = 0;
    int bp_no_upd_count = 0; 
    int no_upd_count = 0;
    int bp_no_upd_count_max = 20;
    int opt_no_upd_count = 0;
    int opt_no_upd_count_max = 5;
    for (iter = 0; iter < alg.maxiter; iter++) {
        // cout << "bpcount";
        while (true) {
            backwardpass(corridor);
            if (!bp.failed) {break;}

            // in case dead loop in bp
            if (bp.reg == 24 && bp.failed){
                bp_no_upd_count++;
            } else {
                bp_no_upd_count = 0;
            }
            if (bp_no_upd_count > bp_no_upd_count_max){
                break;
            }           
        }
        forwardpass(corridor);

        // for record
        PolyCoeff = MatrixXd::Zero(N, num_ctrlP * dim );
        PolyTime  = VectorXd::Zero(N);
        sysparam2polyFunc();
        bool timePosiInd = true;
        for (int i = 0; i < N; i++){
            if (PolyTime(i) < 0) timePosiInd = false;
        }
        if (!timePosiInd){
            rtn = -3.0;
            ROS_WARN("~~~negative time~~~");
            cout  << "iter:" << iter << "timealloc:" << PolyTime.transpose() << endl;
            break;
        }
        PolyCoeffTraj.emplace_back(PolyCoeff);
        PolyTimeTraj.emplace_back(PolyTime);
        costTraj.emplace_back(fp.cost);
        costqTraj.emplace_back(fp.costq);
        steps.emplace_back(fp.stepsize);


        //  -----------termination conditions---------------
        if (std::max(bp.opterr, alg.mu)<=alg.tol){
            ROS_WARN("~~~Optimality reached~~~");
            break;
        }
        
        if (bp.opterr <= 0.2*alg.mu) {
            alg.mu = std::max(alg.tol/10.0, std::min(0.2*alg.mu, pow(alg.mu, 1.2) ) );
            fp.resetfilter(alg);
            bp.resetreg();
        }

        int count = 0;
        VectorXd countv = VectorXd::Zero(N);
        for (int i = 0; i < N; i++){
            for (int j = 0; j < fp.c[i].size(); j++){
                if (fp.c[i](j) >= 2.0e-4){
                    count++;
                    countv(i)++;
                }
            }
        }
        if (count == 0){

            if (zero_init_flag){
                infeas = false;
                rtn = 2;
                ROS_WARN("~~~zero init step 1 find feasible~~~");
                break;
            }

            if (!zero_init_flag && !line_init_flag){
                if ((pow((fp.costq - costqTraj.end()[-2] ), 2) < costqTraj.end()[-2] * 1.0e-2)) // keep opterr or not???
                {
                    opt_no_upd_count ++;
                } else {
                    opt_no_upd_count = 0;
                }

                
                if ((pow((fp.cost - costTraj.end()[-2] ), 2) < costTraj.end()[-2] * 1.0e-2) && bp.opterr < 5.0e1){
                        rtn = 1;
                        ROS_WARN("~~~zero init step 2 find feasible and terminated with some opt conditions~~~");
                        break;
                    }            
            }

            if (line_init_flag) {
                if ((pow((fp.cost - costTraj.end()[-2] ), 2) < costTraj.end()[-2] * 0.01))
                {
                        line_failed = false;
                        ROS_WARN("~~~line init find feasible and terminated with some opt conditions~~~");
                        break;
                } 
            }

        }

        if (bp_no_upd_count > bp_no_upd_count_max){
            rtn = -4.0;
            ROS_WARN("~~~ bp no update, terminate prematurely ~~~");
            break;
        }    

        if (line_init_flag) {
            // premature termination for non-update case
            if (fp.stepsize < 1.0e-6){
                no_upd_count++;
            } else {
                no_upd_count = 0;
            }
            if (no_upd_count > 100){
                ROS_WARN("~~~ line init no update, terminate prematurely ~~~");
                break;
            }
        }
   
        
    }    

    ros::Time time_after_optimization = ros::Time::now();
    compTime = (time_after_optimization - time_before_optimization).toSec();
    ROS_WARN("[ddp_global_planner] Time consumation of %d ddp optimization is: %f", iter, compTime);
    
    iter_used = iter;
    fp.finalroll();
    jerkCost = fp.jerkCost;
    // for output and plot
    PolyCoeff = MatrixXd::Zero(N, num_ctrlP * dim );
    PolyTime  = VectorXd::Zero(N);
    ddpobj = 0.0;
    sysparam2polyFunc();
    ddpobj = fp.cost;
    BezCoeff = poly2bezFunc();

    // convert  dim * num_ctrlP back to num_ctrlP * dim to plot
    for (int i = 0; i < N; i++){
        VectorXd tempv = BezCoeff.row(i); 
        Map<MatrixXd> BezCoeffi_mat(tempv.data(), fp.dim, fp.num_ctrlP);
        MatrixXd tempm = BezCoeffi_mat.transpose();
        Map<RowVectorXd> BezCoeffi(tempm.data(), BezCoeffi_mat.size());
        BezCoeff.row(i) =  BezCoeffi;
    }   
    return rtn;
}

void ddpTrajOptimizer::backwardpass(const decomp_cvx_space::FlightCorridor &corridor)
{
    // ROS_WARN("-------backwardpass Once");
    int N = fp.N;
    int sys_order = fp.sys_order;
    int dim = fp.dim;
    Vector2d dV(0.0,0.0);
    double c_err = 0.0;
    double mu_err = 0.0;
    double Qu_err = 0.0;

    // cout << "fpstep" << fp.step << endl;
    if (fp.failed || bp.failed){
        bp.reg = bp.reg + 1.0;
    } else {
        if (fp.step == 0){ // should be first element
            bp.reg = bp.reg - 1.0;
            // ROS_WARN("entered reduce");
        } else {
            if (fp.step <= 3) {
                bp.reg = bp.reg;
            } else {
                bp.reg = bp.reg + 1.0;
            }
        }
    }


    if (bp.reg < 0.0){
        bp.reg = 0.0;
    } else {
        if (bp.reg > 24.0) {
            bp.reg = 24.0;
        }
    }

    if(!fp.failed){  // only recompute after fp success
        fp.computeall(corridor);
    }
    VecOfVecXd x = fp.getx(); 
    VecOfVecXd u = fp.getu();
    VecOfVecXd c = fp.getc();
    VecOfVecXd y = fp.gety();
    VecOfVecXd s = fp.gets();
    VecOfVecXd mu = fp.getmu();  

    // double V = fp.getp();
    VectorXd Vx = fp.getpx();
    MatrixXd Vxx = fp.getpxx();

    VecOfMatXd fx = fp.getfx();
    VecOfMatXd fu = fp.getfu();
    VecOfMatXd fxx = fp.getfxx();
    VecOfMatXd fxu = fp.getfxu();
    VecOfMatXd fuu = fp.getfuu();

    VecOfMatXd qx = fp.getqx();
    VecOfMatXd qu = fp.getqu();

    VecOfMatXd qxx = fp.getqxx();
    VecOfMatXd qxu = fp.getqxu();
    VecOfMatXd quu = fp.getquu();

    VecOfMatXd cx = fp.getcx();
    VecOfMatXd cu = fp.getcu();


    for (int i = N-1; i>=0; i--){
        MatrixXd Qx = qx[i] + cx[i].transpose() * s[i] + fx[i].transpose() * Vx;
        MatrixXd Qu = qu[i] + cu[i].transpose() * s[i] + fu[i].transpose() * Vx; // (5b)


        // tensor contraction
        MatrixXd tempm1, tempm2;
        tempm1 = MatrixXd::Zero(sys_order*dim, sys_order*dim+1);
        tempm2 = MatrixXd::Zero(sys_order*dim+1, sys_order*dim+1);

        MatrixXd fxiVxx = fx[i].transpose() * Vxx;
        MatrixXd Qxx = qxx[i] + fxiVxx * fx[i];
        MatrixXd Qxu = qxu[i] + fxiVxx * fu[i] + tempm1;
        MatrixXd Quu = quu[i] + fu[i].transpose() * Vxx * fu[i] + tempm2;  // (5c-5e)
        Quu = 0.5 * (Quu + Quu.transpose());

        MatrixXd kK;
        MatrixXd ku;
        MatrixXd Ku;
        VectorXd r;

        // DiagonalMatrix<double, Dynamic> S = s[i].asDiagonal();
        MatrixXd Quu_reg = Quu + (pow(fp.reg_exp_base, bp.reg) - 1) * MatrixXd::Identity(sys_order * dim +1, sys_order * dim + 1);

        // MatrixXd Quu_reg = Quu + quu[i] * (pow(1.6, bp.reg) - 1);
        if (alg.infeas)
        {
            // ROS_WARN("--- entered here ----");
            r = (s[i].array() * y[i].array()).array() - alg.mu;
            VectorXd rhat = s[i].array() * (c[i] + y[i]).array() - r.array();
            VectorXd yinv = y[i].array().inverse();
            VectorXd tempv1 = s[i].array() * yinv.array();
            DiagonalMatrix<double, Dynamic>  SYinv = tempv1.asDiagonal(); // y is vector
            MatrixXd cuitSYinvcui = cu[i].transpose() * SYinv * cu[i];
            MatrixXd SYinvcxi = SYinv * cx[i];

            LLT<MatrixXd> lltofQuuReg(Quu_reg + cuitSYinvcui); // compute the Cholesky decomposition 


            if (lltofQuuReg.info() == Eigen::NumericalIssue)
            {
                bp.failed = true;
                bp.opterr = std::numeric_limits<double>::infinity();
                return;
            }


            VectorXd tempv2 = (yinv.array() * rhat.array());
            Qu.noalias() += cu[i].transpose() * tempv2;
            MatrixXd tempQux = Qxu.transpose() + cu[i].transpose() * SYinvcxi;

            MatrixXd tempm(Qu.rows(), Qu.cols()+tempQux.cols());
            tempm << Qu, tempQux;

            kK = - lltofQuuReg.solve(tempm);

            ku = kK.col(0);
            Ku = kK.rightCols(kK.cols() - 1);
            VectorXd cuiku = cu[i]*ku;
            MatrixXd cxiPluscuiKu = cx[i] + cu[i]*Ku;

            bp.ks[i] = yinv.array() * (rhat.array()+s[i].array()*cuiku.array());

            bp.Ks[i] = SYinv * cxiPluscuiKu;
            bp.ky[i] = - (c[i] + y[i]) - cuiku;
            bp.Ky[i] = -cxiPluscuiKu;

            Quu = Quu + cuitSYinvcui;
            Qxu = tempQux.transpose(); //Qxu + cx[i].transpose() * SYinvcui;
            Qxx.noalias() += cx[i].transpose() * SYinvcxi;

            Qx.noalias() += cx[i].transpose() * tempv2;


        } else {
            // ROS_WARN("~~~entered bwd iteration~~~");
            r = (s[i].array() * c[i].array()).array() + alg.mu;
            // r = (S * c[i]).array() + alg.mu;
            VectorXd cinv = c[i].array().inverse();
            VectorXd tempv1 = (s[i].array() * cinv.array());
            DiagonalMatrix<double, Dynamic>  SCinv = tempv1.asDiagonal(); // y is vector

            MatrixXd SCinvcui = SCinv * cu[i];
            MatrixXd SCinvcxi = SCinv * cx[i];
            MatrixXd cuitSCinvcui = cu[i].transpose() * SCinvcui;
            LLT<MatrixXd> lltofQuuReg(Quu_reg - cuitSCinvcui); // compute the Cholesky decomposition 

            // ROS_WARN("~~~LLT~~~");
            if (lltofQuuReg.info() == Eigen::NumericalIssue) 
            {   
                bp.failed = true;
                bp.opterr = std::numeric_limits<double>::infinity();
                return;
            }
            VectorXd tempv2 = (cinv.array() * r.array());
            Qu.noalias() -= cu[i].transpose() * tempv2; // (12b)            
            MatrixXd tempQux = Qxu.transpose() - cu[i].transpose() * SCinvcxi;
            MatrixXd temp(Qu.rows(), Qu.cols()+tempQux.cols());
            temp << Qu, tempQux;

            kK = - lltofQuuReg.solve(temp);
            ku = kK.col(0);
            Ku = kK.rightCols(kK.cols() - 1);
            VectorXd cuiku = cu[i]*ku;
            bp.ks[i] = -(cinv.array() * (r.array()+s[i].array()*cuiku.array()));
            bp.Ks[i] = - (SCinv * (cx[i] + cu[i] * Ku)); // (11) checked
            bp.ky[i] = MatrixXd::Zero(c[i].size(), 1);
            bp.Ky[i] = MatrixXd::Zero(c[i].size(), fp.sys_order*fp.dim);;        
            Quu = Quu - cuitSCinvcui; // (12e)
            Qxu = tempQux.transpose(); //Qxu - cx[i].transpose() * SCinvcui; // (12d)
            Qxx.noalias() -= cx[i].transpose() * SCinvcxi; // (12c)
            Qx.noalias() -= cx[i].transpose() * tempv2; // (12a)
        }
        dV(0) = dV(0) + (ku.transpose() * Qu)(0);

        MatrixXd QxuKu = Qxu * Ku;
        MatrixXd KutQuu = Ku.transpose() * Quu;

        dV(1) = dV(1) + (0.5 * ku.transpose() * Quu * ku)(0);
        Vx = Qx + Ku.transpose() * Qu + KutQuu * ku + Qxu * ku; // (btw 11-12)
        Vxx = Qxx + QxuKu.transpose() + QxuKu + KutQuu * Ku; // (btw 11-12)
        Vxx = 0.5 * ( Vxx + Vxx.transpose() ); // for symmetry

        bp.ku[i] = ku;
        bp.Ku[i] = Ku;

        Qu_err = std::max(Qu_err, Qu.lpNorm<Infinity>()  );
        mu_err = std::max(mu_err, r.lpNorm<Infinity>()  );
        if (alg.infeas){
            c_err=max(c_err, (c[i]+y[i]).lpNorm<Infinity>() );
        }
    }

    bp.failed = false;
    bp.opterr = std::max(std::max( Qu_err, c_err ), mu_err);
    bp.dV = dV;

}


void ddpTrajOptimizer::forwardpass(const decomp_cvx_space::FlightCorridor &corridor)  
{
    // ROS_WARN("------entered forwardpass");
    int N=fp.N;

    VecOfVecXd xold=fp.getx();
    VecOfVecXd uold=fp.getu();
    VecOfVecXd yold=fp.gety();
    VecOfVecXd sold=fp.gets();
    VecOfVecXd cold=fp.getc();

    VecOfVecXd xnew = xold;
    VecOfVecXd unew = uold;
    VecOfVecXd cnew = cold;
    VecOfVecXd ynew = yold;
    VecOfVecXd snew = sold; //initilization
    double cost;
    double costq;
    double logcost;   
    VectorXd qnew = VectorXd::Zero(N);
    double stepsize;
    double err;
    double tau = std::max(0.99, 1-alg.mu);
    VectorXd steplist = pow(2.0, ArrayXd::LinSpaced(11, -10, 0).reverse() );
    int step;
    bool failed;
    // ros::Time time_fp_start = ros::Time::now();
    for (step = 0; step < steplist.size(); step++) {
        failed = false;
        stepsize = steplist(step);
        xnew[0] = xold[0];
        if (alg.infeas){
            for (int i = 0; i < N; i++) {
                ynew[i] = yold[i] + stepsize*bp.ky[i]+bp.Ky[i]*(xnew[i]-xold[i]);
                snew[i] = sold[i] + stepsize*bp.ks[i]+bp.Ks[i]*(xnew[i]-xold[i]);

                if (    (ynew[i].array()<(1-tau)*yold[i].array()).any() || 
                        (snew[i].array()<(1-tau)*sold[i].array()).any()   ){
                    failed = true;
                    break;
                }
                
                unew[i] = uold[i] + stepsize*bp.ku[i]+bp.Ku[i]*(xnew[i]-xold[i]);
                xnew[i+1] = fp.computenextx(xnew[i], unew[i]);
            }
        } else {
            for (int i = 0; i < N; i++) {
                snew[i] = sold[i] + stepsize*bp.ks[i]+bp.Ks[i]*(xnew[i]-xold[i]);
                unew[i] = uold[i] + stepsize*bp.ku[i]+bp.Ku[i]*(xnew[i]-xold[i]);
                cnew[i] = fp.computecminvo(xnew[i], unew[i], i, corridor);


                if (    (cnew[i].array()>(1-tau)*cold[i].array()).any() ||  
                        (snew[i].array()<(1-tau)*sold[i].array()).any()   ){
                    failed = true;
                    break;
                }

                xnew[i+1] = fp.computenextx(xnew[i], unew[i]);
            }
        }
        
        if (failed) {
            continue;
        } else {
            for (int i = 0; i < N; i++) {
                qnew[i] = fp.computeq(xnew[i], unew[i]);
            }

            cost = qnew.sum() + fp.computep(xnew[N], fp.x_d);
            costq = qnew.sum();
            logcost = cost;   
            err = 0.0;            
            if (alg.infeas) {
                for (int i = 0; i < N; i++) {
                    logcost -= alg.mu * ynew[i].array().log().sum();
                    cnew[i] = fp.computecminvo(xnew[i], unew[i], i, corridor);

                    err += (cnew[i]+ynew[i]).lpNorm<1>(); 
                }
                err = std::max(alg.tol, err);
            } else {
                for (int i = 0; i < N; i++) {
                    cnew[i] = fp.computecminvo(xnew[i], unew[i], i, corridor);
                    logcost -= alg.mu * (-cnew[i]).array().log().sum();
                }
                err=0.0;
            }


            Vector2d candidate(logcost, err);
            ros::Time _time = ros::Time::now();


            std::vector<int> columnidtokeep;
            for (int i=0; i<fp.filter.cols(); i++){
                if (candidate(0)>=fp.filter(0, i) && candidate(1)>=fp.filter(1, i)){
                    failed=true;
                    break;                    
                } else if (candidate(0)>fp.filter(0, i) || candidate(1)>fp.filter(1, i)){
                    columnidtokeep.push_back(i);
                }
            }
            if (failed) continue;  
            MatrixXd tempm(2,columnidtokeep.size());
            for (int i=0; i<columnidtokeep.size(); i++){
                tempm.col(i) = fp.filter.col(columnidtokeep[i]);
            }
            fp.filter.resize(2, tempm.cols() + 1);
            fp.filter << tempm, candidate;            
            break;
        }
    }
    if (failed){
        fp.failed=true;
        fp.stepsize=0.0;
    } else {
        fp.cost=cost;
        fp.costq = costq;
        fp.logcost=logcost;
        fp.x=xnew;
        fp.u=unew;
        fp.y=ynew;
        fp.s=snew;
        fp.c=cnew;
        fp.q=qnew;
        fp.err=err;
        fp.stepsize=stepsize;
        fp.step=step;
        fp.failed=false;
    }
}


// convert bezier coefficient to poly coefficient
MatrixXd ddpTrajOptimizer::bez2polyFunc(){
    // ROS_WARN("~~~entered bez2polyFunc~~~");
    int N = PolyTime.size();
    PolyCoeff = MatrixXd::Zero(BezCoeff.rows(), BezCoeff.cols());
    for (int i = 0; i < N; i++){
        fp.t2tau( PolyTime(i) );
        fp.beztau2polyt = fp.poly2bez * fp.t2tauMat;
        VectorXd tempv = PolyTime(i) * BezCoeff.row(i);
        Map<MatrixXd> BezCoeffi_mat(tempv.data(), fp.dim, fp.num_ctrlP);
        MatrixXd PolyCoeffi_mat = BezCoeffi_mat * fp.beztau2polyt; 
        Map<RowVectorXd> PolyCoeffi(PolyCoeffi_mat.data(), PolyCoeffi_mat.size());               
        PolyCoeff.row(i) =  PolyCoeffi;
    }
    return PolyCoeff;
}  

// convert poly coefficient to bezier coefficient
MatrixXd ddpTrajOptimizer::poly2bezFunc(){
    int N = PolyTime.size();
    BezCoeff = MatrixXd::Zero(PolyCoeff.rows(), PolyCoeff.cols());
    for (int i = 0; i < N; i++){
        fp.t2tau( PolyTime(i) );
        fp.polyt2beztau = (fp.poly2bez * fp.t2tauMat).inverse();
        VectorXd tempv = 1.0 / PolyTime(i) * PolyCoeff.row(i);
        Map<MatrixXd> PolyCoeffi_mat(tempv.data(), fp.dim, fp.num_ctrlP);
        MatrixXd BezCoeffi_mat = PolyCoeffi_mat * fp.polyt2beztau;
        Map<RowVectorXd> BezCoeffi(BezCoeffi_mat.data(), BezCoeffi_mat.size());
        BezCoeff.row(i) =  BezCoeffi;
    }
    return BezCoeff;
}        

void ddpTrajOptimizer::sysparam2polyFunc(){
    int N = PolyTime.size();
    for (int i = 0; i < N; i++){
        VectorXd u_temp = fp.u[i];
        VectorXd x_temp = fp.x[i];   
        PolyTime(i) = (u_temp.tail(1))(0);  
        PolyCoeff.row(i).head(fp.sys_order*fp.dim) = fp.barEk_inv.array() * x_temp.array();
        PolyCoeff.row(i).tail(fp.sys_order*fp.dim) = u_temp.head(fp.sys_order*fp.dim);
    } 
}
        
// for call in forward pass or in the function of class fwdPass
void fwdPass::setPmat()
{
    Pmat = MatrixXd::Identity(dim*sys_order,dim*sys_order);
}

void fwdPass::setPmat(MatrixXd Pmat_input)
{
    Pmat = Pmat_input;
}

void fwdPass::time2barFkbarGk(double time)
{
    double Tk = time;
    double Tk2 = Tk * Tk;
    double Tk3 = Tk2 * Tk;
    double Tk4 = Tk3 * Tk;
    double Tk5 = Tk4 * Tk;
    double Tk6 = Tk5 * Tk;
    double Tk7 = Tk6 * Tk;

    MatrixXd Fk(sys_order,sys_order);
    MatrixXd Gk(sys_order,sys_order);

    if (sys_order == 4){
        // fourth order system
        Fk <<   1.0,   Tk, Tk2/2.0, Tk3/6.0, 
                0.0,  1.0,      Tk, Tk2/2.0, 
                0.0,  0.0,     1.0,      Tk,
                0.0,  0.0,     0.0,     1.0; 

        Gk <<      Tk4,    Tk5,     Tk6,     Tk7,
                4*Tk3,  5*Tk4,   6*Tk5,   7*Tk6,
                12*Tk2, 20*Tk3,  30*Tk4,  42*Tk5, 
                24*Tk, 60*Tk2, 120*Tk3, 210*Tk4;  
    }

    if (sys_order == 3){
        // third order system
        Fk <<   1.0,   Tk, Tk2/2.0, 
                0.0,  1.0,      Tk, 
                0.0,  0.0,     1.0; 

        Gk <<   Tk3,   Tk4,    Tk5,
                3*Tk2, 4*Tk3,  5*Tk4,
                6*Tk, 12*Tk2, 20*Tk3;  
    }

    // three dimensional
    barFk=MatrixXd::Zero(dim * sys_order, dim * sys_order);
    barGk=MatrixXd::Zero(dim * sys_order, dim * sys_order);
    for (int i = 0; i < sys_order; i++){
        for (int j = 0; j < sys_order; j++){
            for (int k = 0; k < dim; k++){
                barFk(i*dim+k, j*dim+k) = Fk(i,j);
            }
        }
    }
    for (int i = 0; i < sys_order; i++){
        for (int j = 0; j < sys_order; j++){
            for (int k = 0; k < dim; k++){
                barGk(i*dim+k, j*dim+k) = Gk(i,j);
            }
        }
    }
}

void fwdPass::time2barFkprimebarGkprime(double time)
{
    double Tk = time;
    double Tk2 = Tk * Tk;
    double Tk3 = Tk2 * Tk;
    double Tk4 = Tk3 * Tk;
    double Tk5 = Tk4 * Tk;
    double Tk6 = Tk5 * Tk;
    // double Tk7 = Tk6 * Tk;


    MatrixXd Fkprime(sys_order,sys_order);
    MatrixXd Gkprime(sys_order,sys_order);
    MatrixXd Fkpprime(sys_order,sys_order);
    MatrixXd Gkpprime(sys_order,sys_order);

    if (sys_order == 4){
        // fourth order system
        Fkprime <<  0,  1,Tk, Tk2/2.0,
                    0,  0, 1, Tk,
                    0,  0, 0, 1,
                    0,  0, 0, 0;
        Gkprime <<  4*Tk3, 5*Tk4, 6*Tk5, 7*Tk6, 
                    12*Tk2, 20*Tk3, 30*Tk4, 42*Tk5, 
                    24*Tk, 60*Tk2, 120*Tk3, 210*Tk4,
                    24,    120*Tk, 360*Tk2, 840*Tk3;
        // Fkpprime << 0, 0, 1, Tk,
        //             0, 0, 0, 1,
        //             0, 0, 0, 0,
        //             0, 0, 0, 0;
        // Gkpprime << 12*Tk2, 20*Tk3, 30*Tk4,  42*Tk5, 
        //             24*Tk, 60*Tk2, 120*Tk3, 210*Tk4,
        //             24,    120*Tk, 360*Tk2, 840*Tk3,
        //             0,        120,  720*Tk,2520*Tk2;           
    }

    if (sys_order == 3){
        // third order system
        Fkprime <<  0,  1,Tk,
                    0,  0, 1,
                    0,  0, 0;
        Gkprime <<  3*Tk2,  4*Tk3,  5*Tk4,
                    6*Tk,  12*Tk2, 20*Tk3,
                    6,    24*Tk, 60*Tk2;
        // Fkpprime <<  0, 0, 1,
        //              0, 0, 0,
        //              0, 0, 0;
        // Gkpprime <<  6*Tk,  12*Tk2,  20*Tk3,
        //                 6,   24*Tk,  60*Tk2,
        //                 0,      24,  120*Tk;          
    }

    // three dimensional
    barFkprime=MatrixXd::Zero(dim * sys_order, dim * sys_order);
    barGkprime=MatrixXd::Zero(dim * sys_order, dim * sys_order);
    for (int i = 0; i < sys_order; i++){
        for (int j = i+1; j < sys_order; j++){
            for (int k = 0; k < dim; k++){
                barFkprime(i*dim+k, j*dim+k) = Fkprime(i,j);

            }
        }
    }
    for (int i = 0; i < sys_order; i++){
        for (int j = 0; j < sys_order; j++){
            for (int k = 0; k < dim; k++){
                barGkprime(i*dim+k, j*dim+k) = Gkprime(i,j);
            }
        }
    }
}

void fwdPass::time2barR(double time)
{
    double Tk = time;  // take the last one as time
    double Tk2 = Tk * Tk;
    double Tk3 = Tk2 * Tk;
    double Tk4 = Tk3 * Tk;
    double Tk5 = Tk4 * Tk;
    double Tk6 = Tk5 * Tk;
    double Tk7 = Tk6 * Tk;                 
    MatrixXd R(sys_order,sys_order);
    MatrixXd Rprime(sys_order,sys_order);
    MatrixXd Rpprime(sys_order,sys_order);
    if (sys_order == 4){
        R <<   576*Tk,  1440*Tk2,  2880*Tk3,   5040*Tk4,
            1440*Tk2,  4800*Tk3, 10800*Tk4,  20160*Tk5,
            2880*Tk3, 10800*Tk4, 25920*Tk5,  50400*Tk6,
            5040*Tk4, 20160*Tk5, 50400*Tk6, 100800*Tk7;
        Rprime <<       576,  2880*Tk,  8640*Tk2,   20160*Tk3,
                2880*Tk,  14400*Tk2, 43200*Tk3,  100800*Tk4,
                8640*Tk2, 43200*Tk3, 129600*Tk4,  302400*Tk5,
            20160*Tk3, 100800*Tk4, 302400*Tk5,  705600*Tk6;
        Rpprime <<           0,  2880,  17280*Tk,   60480*Tk2,
                    2880,  28800*Tk, 129600*Tk2,  403200*Tk3,
            17280*Tk, 129600*Tk2, 518400*Tk3,  1512000*Tk4,
            60480*Tk2, 403200*Tk3, 1512000*Tk4, 4233600*Tk5;   
    }
    if (sys_order == 3){
        R <<        36*Tk,  72*Tk2, 120*Tk3,
                   72*Tk2, 192*Tk3, 360*Tk4,
                  120*Tk3, 360*Tk4, 720*Tk5;
        Rprime <<              36,    144*Tk,  360*Tk2,
                            144*Tk,  576*Tk2, 1440*Tk3,
                           360*Tk2, 1440*Tk3, 3600*Tk4;
        Rpprime <<            0,       144,     720*Tk,
                            144,   1152*Tk,   4320*Tk2,
                         720*Tk,  4320*Tk2,  14400*Tk3;  
    }
    barR = MatrixXd::Zero(dim * sys_order,dim * sys_order); // for first 12 elements
    barRprime = MatrixXd::Zero(dim * sys_order,dim * sys_order);   
    barRpprime = MatrixXd::Zero(dim * sys_order,dim * sys_order);   

    for (int i = 0; i < sys_order; i++){
        for (int j = 0; j < sys_order; j++){
            for (int k = 0; k < dim; k++){
                barR(i*dim+k, j*dim+k) = R(i,j);
                barRprime(i*dim+k, j*dim+k) = Rprime(i,j);
                barRpprime(i*dim+k, j*dim+k) = Rpprime(i,j);

            }
        }
    }
}

// monomial basis P(tau) = t2tauMat * P(t)
void fwdPass::t2tau(double time)
{
    double Tk = 1.0 / time;
    double Tk2 = Tk * Tk;
    double Tk3 = Tk2 * Tk;
    double Tk4 = Tk3 * Tk;
    double Tk5 = Tk4 * Tk;
    double Tk6 = Tk5 * Tk;
    double Tk7 = Tk6 * Tk;                 
    VectorXd temp(num_ctrlP);

    if (traj_order == 7){
        temp << 1.0, Tk, Tk2, Tk3, Tk4, Tk5, Tk6, Tk7;
    }
    if (traj_order == 5){
        temp << 1.0, Tk, Tk2, Tk3, Tk4, Tk5;
    }
    t2tauMat = temp.asDiagonal();

    MatrixXd tempm( num_ctrlP , num_ctrlP);

    if (traj_order == 7){
        tempm  <<    1,    0,    0,    0,    0,   0,   0,   0,
				    -7,    7,    0,    0,    0,   0,   0,   0,
				    21,  -42,   21,    0,    0,   0,   0,   0,
				   -35,  105, -105,   35,    0,   0,   0,   0, 
				    35, -140,  210, -140,   35,   0,   0,   0,
				   -21,  105, -210,  210, -105,  21,   0,   0,
				     7,  -42,  105, -140,  105, -42,   7,   0,
				    -1,    7,  -21,   35,  -35,  21,  -7,   1;
    }
    if (traj_order == 5){
        tempm  <<    1,   0,   0,   0,  0,  0,
					-5,   5,   0,   0,  0,  0,
					10, -20,  10,   0,  0,  0,
				   -10,  30, -30,  10,  0,  0,
				     5, -20,  30, -20,  5,  0,
				    -1,   5, -10,  10, -5,  1;
    }


    poly2bez = tempm.transpose(); // Bezier basis B(tau) = poly2bez * P(tau)
}

VectorXd fwdPass::computenextx(VectorXd xnew, VectorXd unew)
{
    time2barFkbarGk((unew.tail(1))(0));
    VectorXd xnext = barFk.triangularView<Upper>() * xnew + barGk * unew.head(sys_order*dim);
    return xnext;
}

VectorXd fwdPass::computec( VectorXd xnew, 
                            VectorXd unew, 
                            int time_id, 
                            const decomp_cvx_space::FlightCorridor &corridor)
{
        VectorXd x_temp = xnew;
        VectorXd u_temp = unew;        
        MatrixXd polyCoeff = MatrixXd::Zero(num_ctrlP, dim);
        for (int j = 0; j<dim; j++){
            for (int k = 0; k<sys_order; k++)
                polyCoeff(k,j) = x_temp(k*dim+j) * Ek_inv[k];
            for (int k = sys_order; k<num_ctrlP; k++)
                polyCoeff(k,j) = u_temp((k-sys_order)*dim+j);  // second part of polyCoeff is exactly the control input                             
        }

        t2tau((u_temp.tail(1))(0));
        double Tk = (u_temp.tail(1))(0);
        double Tk2 = Tk * Tk;
        double Tk3 = Tk2 * Tk;
        double Tk4 = Tk3 * Tk;
        double Tk5 = Tk4 * Tk;
        double Tk6 = Tk5 * Tk;

        MatrixXd posCoeff;        
        polyt2beztau = MatrixXd::Zero(num_ctrlP, num_ctrlP);
        
        if (num_ctrlP == 6){
            polyt2beztau <<  1,        0,           0,          0,      0,    0,
                            1,     Tk/5.0,           0,          0,      0,    0,
                            1, (2*Tk)/5.0,     Tk2/10.0,          0,      0,    0,
                            1, (3*Tk)/5.0, (3*Tk2)/10.0,    Tk3/10.0,      0,    0,
                            1, (4*Tk)/5.0,  (3*Tk2)/5.0, (2*Tk3)/5.0, Tk4/5.0,    0,
                            1,       Tk,            Tk2,       Tk3,      Tk4,    Tk5;
        }
        if (num_ctrlP == 8){
            polyt2beztau <<      1,        0,            0,           0,          0,          0,      0,    0,
                                1,     Tk/7.0,            0,           0,          0,          0,      0,    0,
                                1, (2*Tk)/7.0,      Tk2/21.0,           0,          0,          0,      0,    0,
                                1, (3*Tk)/7.0,       Tk2/7.0,     Tk3/35.0,          0,          0,      0,    0,
                                1, (4*Tk)/7.0,   (2*Tk2)/7.0, (4*Tk3)/35.0,    Tk4/35.0,          0,      0,    0,
                                1, (5*Tk)/7.0, (10*Tk2)/21.0,  (2*Tk3)/7.0,     Tk4/7.0,    Tk5/21.0,      0,    0,
                                1, (6*Tk)/7.0,   (5*Tk2)/7.0,  (4*Tk3)/7.0, (3*Tk4)/7.0, (2*Tk5)/7.0, Tk6/7.0,    0;
        }

        posCoeff =  polyt2beztau.triangularView<Lower>() * polyCoeff;
        vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
        decomp_cvx_space::Polytope pltp = polyhedrons[time_id];
        int num_plane = pltp.planes.size(); // hyperplane num of this polyhedra

        // constraint function
        VectorXd c_temp = VectorXd::Zero(num_plane * num_ctrlP);
        for (int j = 0; j < num_plane; j++){
            for (int k = 0; k < num_ctrlP; k++){
                c_temp(j * num_ctrlP + k) =    pltp.planes[j](0) * posCoeff(k, 0) 
                                            +  pltp.planes[j](1) * posCoeff(k, 1) 
                                            +  pltp.planes[j](2) * posCoeff(k, 2) 
                                            +  pltp.planes[j](3) - 1.0e-10; // ax+by+cz+d <= 0
            }

        }
        return c_temp;
}

VectorXd fwdPass::computecminvo( const VectorXd& x_temp, 
                            const VectorXd& u_temp, 
                            int time_id, 
                            const decomp_cvx_space::FlightCorridor &corridor)
{         
        MatrixXd polyCoeff = MatrixXd::Zero(num_ctrlP, dim);
        for (int j = 0; j<dim; j++){
            for (int k = 0; k<sys_order; k++)
                polyCoeff(k,j) = x_temp(k*dim+j) * Ek_inv[k];
            for (int k = sys_order; k<num_ctrlP; k++)
                polyCoeff(k,j) = u_temp((k-sys_order)*dim+j);  // second part of polyCoeff is exactly the control input                             
        }

        t2tau((u_temp.tail(1))(0));

        MatrixXd posCoeff;
        VectorXd Tkv(7);
        Tkv(0) = (u_temp.tail(1))(0);
        for (int i = 1; i<7; i++){
            Tkv(i) = Tkv(i-1) * Tkv(0);
        }
        // DiagonalMatrix<double, Dynamic> TkD = Tkv.asDiagonal();
        polyt2minvotau = MatrixXd::Zero(num_ctrlP, num_ctrlP);            
        if (num_ctrlP == 6){
            // ros::Time computeminvo_startm = ros::Time::now();
            polyt2minvotau.col(0) = VectorXd::Ones(num_ctrlP);
            for (int i=1; i<num_ctrlP; i++){
                polyt2minvotau.col(i) = polyt2minvotau6_.col(i)*Tkv(i-1);
            }

        }
        if (num_ctrlP == 8){
            polyt2minvotau <<    1.0, -0.05550234532*Tkv(0), -0.03015649917*Tkv(1), -0.02126221895*Tkv(2), -0.01693417728*Tkv(3), -0.01426182871*Tkv(4), -0.01240847743*Tkv(5), -0.01102231757*Tkv(6),
                                    1.0, -0.01444679377*Tkv(0), -0.05377370182*Tkv(1), -0.03962596389*Tkv(2), -0.03098345174*Tkv(3), -0.02595197185*Tkv(4), -0.02265469217*Tkv(5), -0.02032295522*Tkv(6),
                                    1.0,   0.1624020464*Tkv(0), -0.01852864608*Tkv(1), -0.03991647485*Tkv(2), -0.03471412551*Tkv(3), -0.02817909967*Tkv(4), -0.02334086394*Tkv(5), -0.01983518262*Tkv(6),
                                    1.0,   0.3754136457*Tkv(0),  0.09567905677*Tkv(1), -0.00828476061*Tkv(2), -0.03998426794*Tkv(3), -0.04550931887*Tkv(4), -0.04303904825*Tkv(5), -0.03898747361*Tkv(6),
                                    1.0,   0.6245863543*Tkv(0),   0.3448517654*Tkv(1),   0.1690809939*Tkv(2),  0.06557453238*Tkv(3), 0.008157924523*Tkv(4), -0.02134796454*Tkv(5), -0.03470825213*Tkv(6),
                                    1.0,   0.8375979536*Tkv(0),    0.656667261*Tkv(1),   0.4971243973*Tkv(2),   0.3641717116*Tkv(3),   0.2564765275*Tkv(4),   0.1710093785*Tkv(5),   0.1043765621*Tkv(6),
                                    1.0,    1.014446794*Tkv(0),   0.9751198857*Tkv(1),   0.9216452397*Tkv(2),    0.862665368*Tkv(3),   0.8017913027*Tkv(4),    0.740899876*Tkv(5),   0.6810992625*Tkv(6),
                                    1.0,    1.055502345*Tkv(0),    1.080848191*Tkv(1),    1.097299757*Tkv(2),    1.109185085*Tkv(3),    1.118159867*Tkv(4),    1.125060799*Tkv(5),    1.130372771*Tkv(6);
        }
        posCoeff =  polyt2minvotau * polyCoeff;
       
        vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
        decomp_cvx_space::Polytope pltp = polyhedrons[time_id];
        int num_plane = pltp.planes.size(); // hyperplane num of this polyhedra

        // try to be consistent with paper
        VectorXd c_temp = VectorXd::Zero(num_plane * num_ctrlP); // vec( [h, h, ..., h]_{size s_{k} * (n+1)} )
        for (int j = 0; j < num_ctrlP; j++){
            for (int k = 0; k < num_plane; k++){
                c_temp(j * num_plane + k) =    pltp.planes[k](0) * posCoeff(j, 0) 
                                            +  pltp.planes[k](1) * posCoeff(j, 1) 
                                            +  pltp.planes[k](2) * posCoeff(j, 2) 
                                            +  pltp.planes[k](3); // ax+by+cz+d <= 0
            }

        }


        MatrixXd barEkinv = MatrixXd::Zero(sys_order * dim, sys_order * dim);
        for (int j = 0; j < sys_order; j++){
                for (int l = 0; l < dim; l++){
                    barEkinv(j*dim+l, j*dim+l) = Ek_inv[j];
                }
        }   

        VectorXd tempv = MatrixXd::Zero(num_ctrlP * dim, 1);
        tempv.head(sys_order * dim) = barEkinv * x_temp;
        tempv.tail(sys_order * dim) = u_temp.head(sys_order * dim);

        MatrixXd polyt2minvotau_v = MatrixXd::Zero(num_ctrlP-1, num_ctrlP);            
        if (num_ctrlP == 6){
            polyt2minvotau_v.col(1) = polyt2minvotau_v6_.col(1);

            for (int i=2; i<num_ctrlP; i++){
                polyt2minvotau_v.col(i) = polyt2minvotau_v6_.col(i)*Tkv(i-2);
            }            
        }
        if (num_ctrlP == 8){
            double Tk = (u_temp.tail(1))(0);
            double Tk2 = Tk * Tk;
            double Tk3 = Tk2 * Tk;
            double Tk4 = Tk3 * Tk;
            double Tk5 = Tk4 * Tk;
            double Tk6 = Tk5 * Tk;
            double Tk7 = Tk6 * Tk;            
            polyt2minvotau_v <<      0, 1.0,   -0.120457855*Tk,   -0.1012670775*Tk2, -0.09542334482*Tk3, -0.09560892689*Tk4, -0.09741102308*Tk5, -0.09997532547*Tk6,
                                     0, 1.0, 0.009271294098*Tk,    -0.172475818*Tk2,  -0.1741984982*Tk3,  -0.1660275472*Tk4,  -0.1618718157*Tk5,  -0.1592765356*Tk6,
                                     0, 1.0,    0.456193493*Tk, 0.0006682055768*Tk2,  -0.1813536712*Tk3,  -0.2321689979*Tk4,   -0.242203485*Tk5,  -0.2450315643*Tk6,
                                     0, 1.0,                Tk,    0.6062434344*Tk2,   0.2124868687*Tk3, -0.06539993789*Tk4,  -0.2274169855*Tk5,   -0.306569593*Tk6,
                                     0, 1.0,    1.543806507*Tk,     1.632087726*Tk2,    1.446865535*Tk3,    1.119346483*Tk4,   0.7399562815*Tk5,   0.3655462079*Tk6,
                                     0, 1.0,    1.990728706*Tk,       2.7997103*Tk2,    3.428667461*Tk3,    3.887493822*Tk4,    4.190098233*Tk5,    4.352844312*Tk6,
                                     0, 1.0,    2.120457855*Tk,     3.260106488*Tk2,    4.413102165*Tk3,    5.573415573*Tk4,     6.73663391*Tk5,    7.899198684*Tk6;
        }

        barpoly2minvotau_v = MatrixXd::Zero((num_ctrlP-1)*dim, num_ctrlP * dim);
        for (int j = 0; j < num_ctrlP - 1; j++){
            for (int k = 1; k < num_ctrlP; k++){  // since column 0 is zero
                for (int l = 0; l < dim; l++){
                    barpoly2minvotau_v(j*dim+l, k*dim+l) = polyt2minvotau_v(j,k);}
            }
        }

        VectorXd tempcv = barpoly2minvotau_v * tempv;
        VectorXd c_v = VectorXd::Zero( (num_ctrlP - 1) * 2 * dim);
        c_v << tempcv.array() - maxVel, - tempcv.array() - maxVel;

        MatrixXd polyt2minvotau_a = MatrixXd::Zero(num_ctrlP-2, num_ctrlP);            
        if (num_ctrlP == 6){

            polyt2minvotau_a.col(2) = polyt2minvotau_a6_.col(2);

            for (int i=3; i<num_ctrlP; i++){
                polyt2minvotau_a.col(i) = polyt2minvotau_a6_.col(i)*Tkv(i-3);
            }       
               
        }
        if (num_ctrlP == 8){
            double Tk = (u_temp.tail(1))(0);
            double Tk2 = Tk * Tk;
            double Tk3 = Tk2 * Tk;
            double Tk4 = Tk3 * Tk;
            double Tk5 = Tk4 * Tk;
            double Tk6 = Tk5 * Tk;
            double Tk7 = Tk6 * Tk;            
            polyt2minvotau_a <<      0, 0, 2.0, -0.3883116721*Tk, -0.4473610183*Tk2, -0.5155275588*Tk3, -0.6082719729*Tk4, -0.7048746756*Tk5,
                                     0, 0, 2.0,  0.1988991658*Tk, -0.7857737053*Tk2,   -1.10609276*Tk3,  -1.308815686*Tk4,  -1.542088428*Tk5,
                                     0, 0, 2.0,   2.025317398*Tk,  0.7003479062*Tk2, -0.5840066329*Tk3,  -1.407116374*Tk4,  -1.838108138*Tk5,
                                     0, 0, 2.0,   3.974682602*Tk,   4.599078313*Tk2,   3.832572181*Tk3,   2.095794051*Tk4, -0.1214994105*Tk5,
                                     0, 0, 2.0,   5.801100834*Tk,   10.41862963*Tk2,   15.18823258*Tk3,   19.56315198*Tk4,   23.14477611*Tk5,
                                     0, 0, 2.0,   6.388311672*Tk,   13.10588567*Tk2,   22.16183919*Tk3,   33.54071155*Tk4,   47.20632248*Tk5;
        }

        barpoly2minvotau_a = MatrixXd::Zero((num_ctrlP-2)*dim, num_ctrlP * dim);
        for (int j = 0; j < num_ctrlP - 2; j++){
            for (int k = 2; k < num_ctrlP; k++){ // column 0, 1 are zeros
                for (int l = 0; l < dim; l++){
                    barpoly2minvotau_a(j*dim+l, k*dim+l) = polyt2minvotau_a(j,k);}
            }
        }

        VectorXd tempca = barpoly2minvotau_a * tempv;
        VectorXd c_a = VectorXd::Zero( (num_ctrlP - 2) * 2 * dim);
        c_a << tempca.array() - maxAcc, - tempca.array() - maxAcc;

        VectorXd c_rtn = VectorXd::Zero( num_plane * num_ctrlP +  (num_ctrlP - 1) * 2 * dim + (num_ctrlP - 2) * 2 * dim + 1);
        c_rtn << c_temp, c_v, c_a, -(u_temp.tail(1))(0) + 0.3;
        // ROS_WARN("compute minvo causes %f s", (ros::Time::now()-computeminvo_start).toSec());
        if (!minvo_enabled){
            c_rtn = c_rtn.array() - 2.0e-4;
        }
        return c_rtn;
}



double fwdPass::computep(VectorXd& xnew, VectorXd& xnew_d) 
{   
    return 0.5 * ((xnew - xnew_d).transpose() * Pmat * (xnew - xnew_d))(0);
}

double fwdPass::computeq(VectorXd& xnew, VectorXd& unew)
{
    time2barR((unew.tail(1))(0));
     if (time_power == 2){
        return 0.5 * w_snap * (unew.head(sys_order*dim).transpose() * barR * unew.head(sys_order*dim))(0) 
                + 0.5 * (unew.tail(1))(0) * Rtime * (unew.tail(1))(0);
    }
    if (time_power == 1){
        return 0.5 * w_snap * (unew.head(sys_order*dim).transpose() * barR * unew.head(sys_order*dim))(0) 
                + 0.5 * Rtime * (unew.tail(1))(0);
    }   
} 


// values and derivatives along entire traj
void fwdPass::computeall(const decomp_cvx_space::FlightCorridor &corridor)
{
    computeprelated();
    computefrelated();
    computeqrelated();
    computecrelatedminvo(corridor);

}

void fwdPass::computeprelated()
{
    p = computep(x.back(), x_d); //reuse computep with argument
    px = Pmat * (x.back() - x_d);
    pxx = Pmat;
}

void fwdPass::computefrelated() 
{   
    for (int i = 0; i < N; i++){
        VectorXd u_temp = u[i];
        time2barFkbarGk((u_temp.tail(1))(0));
        fx[i] = barFk;
        time2barFkprimebarGkprime((u_temp.tail(1))(0));
        MatrixXd fgradt = barFkprime.triangularView<StrictlyUpper>() * x[i] + barGkprime * u_temp.head(sys_order*dim); // gradient w.r.t. time
        fu[i] = MatrixXd::Zero(sys_order*dim, sys_order*dim+1);
        fu[i] << barGk, fgradt;
    }
}

void fwdPass::computeqrelated()
{
    for (int i = 0; i < N; i++){
        VectorXd u_temp = u[i];
        time2barR((u_temp.tail(1))(0));
        qx[i] = MatrixXd::Zero(sys_order*dim, 1);
        qxx[i] = MatrixXd::Zero(sys_order*dim, sys_order*dim);
        qxu[i] = MatrixXd::Zero(sys_order*dim, sys_order*dim+1);

        qu[i] = MatrixXd::Zero(sys_order*dim+1, 1);
        quu[i] = MatrixXd::Zero(sys_order*dim+1, sys_order*dim+1);
        if (time_power == 2){
            qu[i] << w_snap * barR * u_temp.head(sys_order*dim), 
                    Rtime * u_temp.tail(1) 
                    +  0.5 * w_snap * u_temp.head(sys_order*dim).transpose() * barRprime * u_temp.head(sys_order*dim);
            quu[i] << w_snap * barR, w_snap * barRprime * u_temp.head(sys_order*dim),
                    w_snap * u_temp.head(sys_order*dim).transpose() * barRprime, Rtime + 
                    0.5 * w_snap * (u_temp.head(sys_order*dim).transpose() * barRpprime * u_temp.head(sys_order*dim))(0); 
        }
        if (time_power == 1){
            qu[i] << w_snap * barR * u_temp.head(sys_order*dim), 
                    0.5 * Rtime 
                    +  0.5 * w_snap * u_temp.head(sys_order*dim).transpose() * barRprime * u_temp.head(sys_order*dim);
            quu[i] << w_snap * barR, w_snap * barRprime * u_temp.head(sys_order*dim),
                    w_snap * u_temp.head(sys_order*dim).transpose() * barRprime,
                    0.5 * w_snap * (u_temp.head(sys_order*dim).transpose() * barRpprime * u_temp.head(sys_order*dim))(0); 
        }

    }

}

void fwdPass::computecrelated( const decomp_cvx_space::FlightCorridor &corridor )
{
    vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
    for (int i = 0; i < N; i++){
        decomp_cvx_space::Polytope pltp = polyhedrons[i];
        int num_plane = pltp.planes.size(); // hyperplane num of this polyhedra
        
        VectorXd x_temp = x[i];
        VectorXd u_temp = u[i];
        c[i] = computec(x_temp, u_temp, i, corridor);  // reuse computec    

        MatrixXd barpoly2beztau = MatrixXd::Zero(num_ctrlP * dim, num_ctrlP * dim);
        Matrix3d eye3 = MatrixXd::Identity(dim, dim);
        for (int j = 0; j < dim; j++){
            for (int k = 0; k < dim; k++){
                barpoly2beztau.block(j*num_ctrlP,k*num_ctrlP,num_ctrlP,num_ctrlP) = eye3(j,k) * polyt2beztau;
            }
        }
        MatrixXd hatA = MatrixXd::Zero(num_plane * num_ctrlP, dim * num_ctrlP);
        for (int k = 0; k < num_plane; k++){
            for (int l = 0; l < num_ctrlP; l++){
                    Vector3d tempv = pltp.planes[k].head(dim);
                    hatA.block(k*num_ctrlP+l, l*dim, 1, dim) = tempv.transpose();
            }
        }        
        MatrixXd barEkinv = MatrixXd::Zero(sys_order * dim, sys_order * dim);
        for (int j = 0; j < sys_order; j++){
                for (int l = 0; l < dim; l++){
                    barEkinv(j*dim+l, j*dim+l) = Ek_inv[j];
                }
        }   
        MatrixXd hatAbarpoly2beztau = hatA * barpoly2beztau * comMat.transpose();
        // block diagonal matrix
        MatrixXd temp = MatrixXd::Zero(num_ctrlP * dim, num_ctrlP * dim);
        temp.block(0,0,sys_order * dim, sys_order * dim) = barEkinv;
        temp.block(sys_order * dim,sys_order * dim, sys_order * dim, sys_order * dim) = MatrixXd::Identity(sys_order * dim,sys_order * dim);
        // concatenated vector
        VectorXd tempv = MatrixXd::Zero(num_ctrlP * dim, 1);
        tempv.head(sys_order * dim) = barEkinv * x_temp;
        tempv.tail(sys_order * dim) = u_temp.head(sys_order * dim);
        
        cx[i] = hatAbarpoly2beztau * temp.leftCols(sys_order * dim);

        // computing cu[i]
        MatrixXd poly2beztau_dt(num_ctrlP, num_ctrlP);
        double Tk = (u_temp.tail(1))(0);
        double Tk2 = Tk * Tk;
        double Tk3 = Tk2 * Tk;
        double Tk4 = Tk3 * Tk;
        double Tk5 = Tk4 * Tk;
        double Tk6 = Tk5 * Tk;

        // computed by matlab
        if (num_ctrlP == 8){
            poly2beztau_dt <<       0,   0,          0,            0,           0,           0,          0,      0,
                                     0, 1/7.0,          0,            0,           0,           0,          0,      0,
                                    0, 2/7.0,  (2*Tk)/21.0,            0,           0,           0,          0,      0,
                                     0, 3/7.0,   (2*Tk)/7.0,  (3*Tk2)/35.0,           0,           0,          0,      0,
                                     0, 4/7.0,   (4*Tk)/7.0, (12*Tk2)/35.0, (4*Tk3)/35.0,           0,          0,      0,
                                     0, 5/7.0, (20*Tk)/21.0,   (6*Tk2)/7.0,  (4*Tk3)/7.0, (5*Tk4)/21.0,          0,      0,
                                     0, 6/7.0,  (10*Tk)/7.0,  (12*Tk2)/7.0, (12*Tk3)/7.0, (10*Tk4)/7.0, (6*Tk5)/7.0,      0,
                                     0,   1,       2*Tk,       3*Tk2,      4*Tk3,      5*Tk4,     6*Tk5, 7*Tk6;
        }
        if (num_ctrlP == 6){
            poly2beztau_dt <<      0,   0,        0,           0,          0,      0,
                                 0, 1/5.0,        0,           0,          0,      0,
                                 0, 2/5.0,     Tk/5.0,           0,          0,    0,
                                 0, 3/5.0, (3*Tk)/5.0, (3*Tk2)/10.0,          0,   0,
                                 0, 4/5.0, (6*Tk)/5.0,  (6*Tk2)/5.0, (4*Tk3)/5.0,  0,
                                 0,   1,     2*Tk,      3*Tk2,     4*Tk3, 5*Tk4;
        }
        MatrixXd barpoly2beztau_dt = MatrixXd::Zero(num_ctrlP * dim, num_ctrlP * dim);

        for (int j = 0; j < dim; j++){
            for (int k = 0; k < dim; k++){
                barpoly2beztau_dt.block(j*num_ctrlP,k*num_ctrlP,num_ctrlP,num_ctrlP) = eye3(j,k) * poly2beztau_dt;
            }
        }
        cu[i] = MatrixXd::Zero(num_ctrlP * num_plane, sys_order * dim  + 1);
        cu[i] << hatAbarpoly2beztau * temp.rightCols(sys_order * dim), hatA * 
        barpoly2beztau_dt.triangularView<Lower>() * comMat.transpose() * tempv;

    }
}

void fwdPass::computecrelatedminvo( const decomp_cvx_space::FlightCorridor &corridor )
{
    vector<decomp_cvx_space::Polytope> polyhedrons  = corridor.polyhedrons;
    for (int i = 0; i < N; i++){
        decomp_cvx_space::Polytope pltp = polyhedrons[i];
        int num_plane = pltp.planes.size(); // hyperplane num of this polyhedra
        
        VectorXd x_temp = x[i];
        VectorXd u_temp = u[i];
        c[i] = computecminvo(x_temp, u_temp, i, corridor);  // reuse computecminvo    

        MatrixXd hatAbarpoly2minvotau = MatrixXd::Zero(num_ctrlP * num_plane, num_ctrlP * dim);  

        for (int j = 0; j < num_ctrlP; j++){
            for (int k = 0; k < num_ctrlP; k++){
                for (int ld = 0; ld < num_plane; ld++){
                    Vector3d tempv = pltp.planes[ld].head(dim);
                    hatAbarpoly2minvotau.block(j*num_plane+ld,k*dim,1,dim) =  polyt2minvotau(j,k) * tempv.transpose();
                }
            }
        }

        MatrixXd barEkinv = MatrixXd::Zero(sys_order * dim, sys_order * dim);
        for (int j = 0; j < sys_order; j++){
                for (int l = 0; l < dim; l++){
                    barEkinv(j*dim+l, j*dim+l) = Ek_inv[j];
                }
        }   

        // block diagonal matrix
        MatrixXd temp = MatrixXd::Zero(num_ctrlP * dim, num_ctrlP * dim);
        temp.block(0,0,sys_order * dim, sys_order * dim) = barEkinv;
        temp.block(sys_order * dim,sys_order * dim, sys_order * dim, sys_order * dim) = MatrixXd::Identity(sys_order * dim,sys_order * dim);
        // concatenated vector
        VectorXd tempv = MatrixXd::Zero(num_ctrlP * dim, 1);
        tempv.head(sys_order * dim) = barEkinv * x_temp;
        tempv.tail(sys_order * dim) = u_temp.head(sys_order * dim);
          
        cx[i] = MatrixXd::Zero(c[i].size(), sys_order * dim);
        MatrixXd tempcxv = barpoly2minvotau_v * temp.leftCols(sys_order * dim);
        MatrixXd tempcxa = barpoly2minvotau_a * temp.leftCols(sys_order * dim);
        cx[i] << hatAbarpoly2minvotau * temp.leftCols(sys_order * dim),
                  tempcxv,
                 -tempcxv,
                  tempcxa,
                 -tempcxa,
                 MatrixXd::Zero(1, sys_order * dim);

        // computing cu[i]
        MatrixXd poly2minvotau_dt(num_ctrlP, num_ctrlP);
        MatrixXd poly2minvotau_v_dt(num_ctrlP-1, num_ctrlP);
        MatrixXd poly2minvotau_a_dt(num_ctrlP-2, num_ctrlP);
  
        double Tk = (u_temp.tail(1))(0);
        double Tk2 = Tk * Tk;
        double Tk3 = Tk2 * Tk;
        double Tk4 = Tk3 * Tk;
        double Tk5 = Tk4 * Tk;
        double Tk6 = Tk5 * Tk;

        // computed by matlab
        if (num_ctrlP == 8){
            poly2minvotau_dt <<      0, -0.05550234532, -0.06031299835*Tk, -0.06378665684*Tk2, -0.06773670913*Tk3, -0.07130914354*Tk4, -0.07445086459*Tk5, -0.07715622296*Tk6,
                                     0, -0.01444679377,  -0.1075474036*Tk,  -0.1188778917*Tk2,  -0.1239338069*Tk3,  -0.1297598593*Tk4,   -0.135928153*Tk5,  -0.1422606866*Tk6,
                                     0,   0.1624020464, -0.03705729217*Tk,  -0.1197494245*Tk2,   -0.138856502*Tk3,  -0.1408954984*Tk4,  -0.1400451837*Tk5,  -0.1388462783*Tk6,
                                     0,   0.3754136457,   0.1913581135*Tk, -0.02485428183*Tk2,  -0.1599370718*Tk3,  -0.2275465944*Tk4,  -0.2582342895*Tk5,  -0.2729123153*Tk6,
                                     0,   0.6245863543,   0.6897035308*Tk,   0.5072429816*Tk2,   0.2622981295*Tk3,  0.04078962261*Tk4,  -0.1280877873*Tk5,  -0.2429577649*Tk6,
                                     0,   0.8375979536,    1.313334522*Tk,    1.491373192*Tk2,    1.456686846*Tk3,    1.282382638*Tk4,    1.026056271*Tk5,   0.7306359347*Tk6,
                                     0,    1.014446794,    1.950239771*Tk,    2.764935719*Tk2,    3.450661472*Tk3,    4.008956514*Tk4,    4.445399256*Tk5,    4.767694838*Tk6,
                                     0,    1.055502345,    2.161696383*Tk,    3.291899272*Tk2,    4.436740339*Tk3,    5.590799333*Tk4,    6.750364793*Tk5,    7.912609399*Tk6;

            poly2minvotau_v_dt <<    0, 0,   -0.120457855,  -0.2025341549*Tk, -0.2862700344*Tk2, -0.3824357076*Tk3, -0.4870551154*Tk4, -0.5998519528*Tk5,
                                     0, 0, 0.009271294098,  -0.3449516361*Tk, -0.5225954946*Tk2, -0.6641101889*Tk3, -0.8093590783*Tk4, -0.9556592135*Tk5,
                                     0, 0,    0.456193493, 0.001336411154*Tk, -0.5440610136*Tk2, -0.9286759916*Tk3,  -1.211017425*Tk4,  -1.470189386*Tk5,
                                     0, 0,            1.0,    1.212486869*Tk,  0.6374606062*Tk2, -0.2615997515*Tk3,  -1.137084927*Tk4,  -1.839417558*Tk5,
                                     0, 0,    1.543806507,    3.264175453*Tk,   4.340596606*Tk2,   4.477385934*Tk3,   3.699781407*Tk4,   2.193277248*Tk5,
                                     0, 0,    1.990728706,    5.599420599*Tk,   10.28600238*Tk2,   15.54997529*Tk3,   20.95049117*Tk4,   26.11706587*Tk5,
                                     0, 0,    2.120457855,    6.520212975*Tk,   13.23930649*Tk2,   22.29366229*Tk3,   33.68316955*Tk4,    47.3951921*Tk5;


            poly2minvotau_a_dt <<    0, 0, 0, -0.3883116721, -0.8947220367*Tk, -1.546582676*Tk2, -2.433087892*Tk3,  -3.524373378*Tk4,
                                     0, 0, 0,  0.1988991658,  -1.571547411*Tk, -3.318278281*Tk2, -5.235262744*Tk3,  -7.710442142*Tk4,
                                     0, 0, 0,   2.025317398,   1.400695812*Tk, -1.752019899*Tk2, -5.628465496*Tk3,  -9.190540688*Tk4,
                                     0, 0, 0,   3.974682602,   9.198156626*Tk,  11.49771654*Tk2,  8.383176206*Tk3, -0.6074970527*Tk4,
                                     0, 0, 0,   5.801100834,   20.83725926*Tk,  45.56469773*Tk2,  78.25260793*Tk3,   115.7238806*Tk4,
                                     0, 0, 0,   6.388311672,   26.21177134*Tk,  66.48551756*Tk2,  134.1628462*Tk3,   236.0316124*Tk4;
            
        }
        if (num_ctrlP == 6){
            poly2minvotau_dt <<  0, -0.06471861202, -0.07456016972*Tk, -0.07732913382*Tk2, -0.08110292972*Tk3, -0.08391365186*Tk4,
                                 0,  0.03314986096,  -0.1309622842*Tk,  -0.1659139141*Tk2,  -0.1745087581*Tk3,  -0.1835819558*Tk4,
                                 0,   0.3375528997,    0.116724651*Tk, -0.08760099494*Tk2,  -0.1876155165*Tk3,  -0.2188223973*Tk4,
                                 0,   0.6624471003,   0.7665130522*Tk,   0.5748858272*Tk2,   0.2794392069*Tk3, -0.01446421554*Tk4,
                                 0,    0.966850139,    1.736438272*Tk,    2.278234886*Tk2,    2.608420264*Tk3,    2.755330489*Tk4,
                                 0,    1.064718612,    2.184314278*Tk,    3.324275878*Tk2,    4.472094873*Tk3,    5.619800295*Tk4;

            poly2minvotau_v_dt <<    0, 0, -0.1423379297, -0.2665484655*Tk, -0.3726316072*Tk2, -0.5052170278*Tk3,
                                     0, 0,  0.1887439858, -0.3662636595*Tk, -0.7399820545*Tk2, -0.9285573245*Tk3,
                                     0, 0,           1.0,   1.082203315*Tk,  0.2466099449*Tk2, -0.9733898632*Tk3,
                                     0, 0,   1.811256014,   4.501272426*Tk,   7.145008354*Tk2,   9.129623752*Tk3,
                                     0, 0,    2.14233793,   6.587479113*Tk,   13.33542355*Tk2,   22.34154157*Tk3;

            poly2minvotau_a_dt <<    0, 0, 0, -0.4472869252, -1.226758663*Tk, -1.921966087*Tk2,
                                     0, 0, 0,   1.223711659, -1.110542869*Tk, -5.563856457*Tk2,
                                     0, 0, 0,   4.776288341,  13.09976386*Tk,  20.52343517*Tk2,
                                     0, 0, 0,   6.447286925,  26.35153674*Tk,  66.13988387*Tk2;
        }
        MatrixXd hatAbarpoly2minvotau_dt = MatrixXd::Zero(num_ctrlP * num_plane, num_ctrlP * dim);
        MatrixXd barpoly2minvotau_v_dt = MatrixXd::Zero((num_ctrlP-1) * dim, num_ctrlP * dim);
        MatrixXd barpoly2minvotau_a_dt = MatrixXd::Zero((num_ctrlP-2) * dim, num_ctrlP * dim);

        for (int j = 0; j < num_ctrlP; j++){
            for (int k = 0; k < num_ctrlP; k++){
                for (int ld = 0; ld < num_plane; ld++){
                    Vector3d tempv = pltp.planes[ld].head(dim);
                    hatAbarpoly2minvotau_dt.block(j*num_plane+ld,k*dim,1,dim) =  poly2minvotau_dt(j,k) * tempv.transpose();
                }
            }
        }
      
        for (int j = 0; j < num_ctrlP-1; j++){
            for (int k = 2; k < num_ctrlP; k++){
                for (int l = 0; l < dim; l++){
                    barpoly2minvotau_v_dt(j*dim+l, k*dim+l) = poly2minvotau_v_dt(j,k);
                }
            }
        }
        for (int j = 0; j < num_ctrlP-2; j++){
            for (int k = 3; k < num_ctrlP; k++){
                for (int l = 0; l < dim; l++){
                    barpoly2minvotau_a_dt(j*dim+l, k*dim+l) = poly2minvotau_a_dt(j,k);
                }
            }
        }
        cu[i] = MatrixXd::Zero(c[i].size(), sys_order * dim  + 1);
  
        MatrixXd tempcuv = barpoly2minvotau_v * temp.rightCols(sys_order * dim);
        VectorXd tempcuv_time = barpoly2minvotau_v_dt * tempv;
        MatrixXd tempcua = barpoly2minvotau_a * temp.rightCols(sys_order * dim);
        VectorXd tempcua_time = barpoly2minvotau_a_dt * tempv;

        cu[i] << hatAbarpoly2minvotau * temp.rightCols(sys_order * dim), hatAbarpoly2minvotau_dt * tempv,
                 tempcuv, tempcuv_time,
                 -tempcuv, -tempcuv_time,
                 tempcua, tempcua_time,
                 -tempcua, -tempcua_time,
                 MatrixXd::Zero(1, sys_order * dim), -1.0;

    }
}



void fwdPass::initialroll(const decomp_cvx_space::FlightCorridor  &corridor)
{
    q = VectorXd::Zero(N);
    for (int i = 0; i < N; i++) {
        VectorXd x_temp = x[i];
        VectorXd u_temp = u[i];
        c[i] = computecminvo(x_temp, u_temp, i, corridor);
        q(i) = computeq(x_temp, u_temp);  //  compute cost then used in resetfilter
        x[i+1] = computenextx(x_temp, u_temp);
    }
    cost = q.sum() + computep(x[N], x_d);
    costq = q.sum();
}



void fwdPass::finalroll()
{
    jerkCost = VectorXd::Zero(N);
    for (int i = 0; i < N; i++) {
        VectorXd x_temp = x[i];
        VectorXd u_temp = u[i];
        time2barR((u_temp.tail(1))(0));

        jerkCost(i) = (u_temp.head(sys_order*dim).transpose() * barR * u_temp.head(sys_order*dim))(0); 
    }
}

void fwdPass::resetfilter( algParam alg)
{
    logcost = cost;   
    err = 0.0;
    if (alg.infeas) {
        for (int i = 0; i < N; i++) {
            logcost -= alg.mu * y[i].array().log().sum();
            err += (c[i]+y[i]).lpNorm<1>();
        }
        if (err < alg.tol){
            err = 0.0;
        }

    } else {
        for (int i = 0; i < N; i++) {
            logcost -= alg.mu * (-c[i]).array().log().sum();
            err = 0.0;
        }
    }

    filter = MatrixXd::Zero(2,1);
    Vector2d tempv(logcost, err);

    filter = tempv;
    step = 0;
    failed = false;
}

void fwdPass::removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void bwdPass::resetreg()
{
    reg = 0.0;
    failed = false;
    recovery = 0;
}

void bwdPass::initreg(double regvalue)
{
    reg = regvalue;
    failed = false;
    recovery = 0;
}

