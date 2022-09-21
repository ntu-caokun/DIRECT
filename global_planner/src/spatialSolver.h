class IndoorOptProblem
{
    private:
        int floor;
        Eigen::Vector3d x0_pack;
        Eigen::Vector3d xf_pack;
        double tfweight;
        int obj_order;
        int poly_order;
        int connect_order;
        bool verbose;
        double margin;
        bool is_limit_vel;
        bool is_limit_acc;
        double vel_limit;
        double acc_limit;

        int boxes; // need change

        int num_box;
        Eigen::VectorXd room_time;
        Eigen::MatrixXd bz;
        Eigen::MatrixXd bzM;
        Eigen::MatrixXd MQM;

        double abs_obj_tol = 1.0e-3;
        double rel_obj_tol = 1.0e-3;
        double grad_tol = 1.0e-3;

        vector<double> timeProfile;
        vector<Eigen::Vector2d> log;

    public:
        IndoorOptProblem(double tgpInput, double tfweightInput, int connect_orderInput, bool verboseInput){
            
        };
        ~IndoorOptProblem(){};

        void set_tfweight(double weight){
                tfweight = weight;
            }

        void  construct_prob(Eigen::Vector3d x0_pack, Eigen::Vector3d xf_pack, int poly_order, int obj_order, int connect_order){
            return;
        }

    // def set_x0_pack_value(self, *args):
    //     """Set the contents of x0pack"""
    //     ff = zip(self.x0_pack, args)
    //     print(ff)
    //     for tmp, val in zip(self.x0_pack, args):
    //         tmp[:] = val[:]

    // def set_xf_pack_value(self, *args):
    //     """Set the contents of xfpack"""
    //     for tmp, val in zip(self.xf_pack, args):
    //         tmp[:] = val[:]

    // def solve_with_room_time(self, rm_time):
    //     """Specify room time and solve the problem"""
    //     raise NotImplementedError("Subclass should implement solve_with_room_time function")

    void get_output_coefficients(Eigen::VectorXd &room_time_save, Eigen::MatrixXd &poly_coef){
        // int n_room = num_box;
        
        // poly_coef = np.zeros([n_room, self.poly_order + 1, 3])
        
        // # coefficients in bezier and scaled form
        // mat_x = np.reshape(self.sol, (n_room, 3, self.poly_order + 1))
        // mat_x = np.transpose(mat_x, (0, 2, 1))

        // # change coefficients to monomial and unscaled form
        // for i in range(n_room):
        //     poly_coef[i, :, :] = self.bzM.dot(mat_x[i]) * self.room_time[i]

        // return self.room_time.copy(), poly_coef
        }

    void get_coef_matrix(Eigen::MatrixXd &coef_mat){
        // """Return coefficients"""
        // coef_mat = np.reshape(self.sol, (self.num_box, 3, self.poly_order + 1))
        // coef_mat = np.transpose(coef_mat, (0, 2, 1))
        // #import pdb; pdb.set_trace()
        // for i in range(self.num_box):
        //     coef_mat[i] *= self.room_time[i]
        // return coef_mat
    }
    // def get_coef_matrix2(self):
    //     """
    //     Return Bezier coefficients
    //     """
    //     coef_mat = np.reshape(self.sol, (self.poly_order + 1, 3, self.num_box), order='F')
    //     coef_mat = coef_mat.transpose((0,2,1))
    //     #import pdb;pdb.set_trace()
    //     for i in range(self.num_box):
    //         coef_mat[:,i,:] *= self.room_time[i]
        
    //     return coef_mat


    // def from_coef_matrix(self, mat_in):
    //     """Assign values to sol based on the input coefficient matrix."""
    //     self.sol = np.zeros(mat_in.size)
    //     coef_mat = np.reshape(self.sol, (self.num_box, 3, self.poly_order + 1))
    //     for i in range(self.num_box):
    //         coef_mat[i] = mat_in[i].T / self.room_time[i]
    //     with np.printoptions(precision=4, linewidth=10000):
    //         print(self.sol)

    // def get_output_path(self, n):
    //     """Get a path for output that is linspace in time.

    //     :param n: int, the number of nodes for a path
    //     :return: float, the total time for this problem
    //     :return: ndarray, (n, 2) the optimal path
    //     """
    //     cum_sum_time = np.cumsum(self.room_time)
    //     output = np.zeros((n, 3))  # the output trajectory
    //     sample_time = np.linspace(0, cum_sum_time[-1], n)
    //     n_room = self.num_box

    //     # get all the coef of polynomials
    //     t, all_poly_coeffs = self.get_output_coefficients()
    //     for i in range(n_room):
    //         # poly_coef = self.bzM.dot(mat_x[i]) * self.room_time[i]
    //         poly_coef = all_poly_coeffs[i,:,:]
    //         if i == 0:
    //             t_mask = sample_time <= cum_sum_time[0]
    //             use_s = sample_time[t_mask] / cum_sum_time[0]
    //         else:
    //             t_mask = (sample_time > cum_sum_time[i - 1]) & (sample_time <= cum_sum_time[i])
    //             use_s = (sample_time[t_mask] - cum_sum_time[i - 1]) / self.room_time[i]
    //         output[t_mask, 0] = np.polyval(poly_coef[:, 0][::-1], use_s)
    //         output[t_mask, 1] = np.polyval(poly_coef[:, 1][::-1], use_s)
    //         output[t_mask, 2] = np.polyval(poly_coef[:, 2][::-1], use_s)
    //     return cum_sum_time[-1], output

    void get_gradient(){
        ROS_WARN("NotImplementedError");

    };

    // def get_gradient_fd(self, h=1e-6):
    //     """Use forward finite difference to approximate gradients."""
    //     grad = np.zeros(self.num_box)
    //     obj0 = self.obj
    //     origin_time = self.room_time.copy()
    //     for i in range(self.num_box):
    //         try_time = origin_time.copy()
    //         try_time[i] += h
    //         self.solve_with_room_time(try_time)
    //         grad[i] = (self.obj - obj0) / h
    //     grad += self.tfweight
    //     return grad

    // def get_gradient_mellinger(self, h=1e-6):
    //     """
    //     Use finite difference described in:
    //     http://www-personal.acfr.usyd.edu.au/spns/cdm/papers/Mellinger.pdf
    //     to approximate gradients.
    //     """
    //     grad = np.zeros(self.num_box)
    //     obj0 = self.obj
    //     origin_time = self.room_time.copy()
    //     for i in range(self.num_box):
    //         m = self.num_box
    //         gi = -1/(m-1) * np.ones(m)
    //         gi[i] = 1
    //         try_time = origin_time.copy()
    //         try_time += h * gi
    //         # print("In gg Mellinger: gi", gi," try_time: ", try_time)
    //         self.solve_with_room_time(try_time)
    //         grad[i] = (self.obj - obj0) / h
    //     grad += self.tfweight
    //     return grad

    void refine_time_by_backtrack(double alpha0, double h, double c, double tau, 
                                  int max_iter, int j_iter, bool logInput, bool timeProfileInput, bool adaptiveLineSearch)
    {                          
        if (logInput){
            Vector2d vtemp;
            vtemp << obj, 0;
            log.pushback(vtemp);
        }

        ros::Time t0 = ros::Time::now();
        int major_iteration = 0;
        int num_prob_solve = 0;
        int converge_reason;
        if (num_box == 1 && tfweight == 0){
            major_iteration = 0;
            num_prob_solve = 0;
            ros::Time t1 = ros::Time::now();
            double time_cost = (t1 - t0).toSec();
            converge_reason = -2;
            return;             
        }

        int n_room = num_box;
        VectorXd t_now = room_time;

        bool converged = false;
        converge_reason = -1;
        num_prob_solve = 0;
        for (int i; i < max_iter; i++){
            if (verbose){
                ROS_WARN("not implemented");
            }
                
           
            double obj0 = obj;
            bool is_okay = true;
            
            if (timeProfileInput){
               ros::Time tBeforeGrad = ros::Time::now();}

            VectorXd grad = get_gradient();
            
            if (timeProfileInput){
                ros::Time tAfterGrad = ros::Time::now();
                timeProfile.pushback( (tAfterGrad - tBeforeGrad).toSec() ); }


            if (tfweight == 0){
                VectorXd normal = VectorXd::Ones(n_room) / std::sqrt(n_room);
                grad = grad - grad * normal * normal;
            }

            if (grad.norm() < grad_tol) {
                if (verbose) ROS_WARN("Gradient too small");
                converged = true;
                converge_reason = -3;
                break;
            }

            double m = - grad.norm();
            VectorXd p = grad / m;
            double alpha_max = (-t_now.array() / p.array()).maxCoeff() - 1.0e-6;
            double alpha;
            if (alpha_max > 0){
                alpha = min(alpha_max, alpha0);
            }else{
                alpha = alpha0;
            }
            double t = -c * m;

            bool alpha_found = false;

            if (timeProfileInput){
               ros::Time tBeforeAlpha = ros::Time::now();}

            for (int j=0; j <j_iter; j++){
                if (verbose) ROS_WARN("Search alpha step");
                VectorXd candid_time = t_now + alpha * p;
                if (adaptiveLineSearch){
                    if (alpha < 1.0e-4)
                        // if self.verbose:
                        //     print_yellow('Stop line search because alpha is too small')
                        break;
                }

                if ( (candid_time < 1.0e-6).any() ) { // ???
                    alpha = tau * alpha;
                    continue;
                }

                solve_with_room_time(candid_time);
                num_prob_solve += 1;

                if (obj < 0) is_solved = false;
                if (!is_solved){
                    alpha = tau * alpha;
                    continue; 
                }               
                double objf = obj;
                
                if (objf < 0){
                    ROS_WARN("Negative Objective!");
                    alpha = tau * alpha;
                    continue;
                }
                // if self.verbose:
                //     print('\talpha ', alpha, ' obj0 ', obj0, ' objf ', objf)
                if (obj0 - objf >= alpha * t || obj0 - objf >= 0.1 * obj0){
                    alpha_found = true;
                    if (adaptiveLineSearch){
                        if (j == 0)
                        {   
                            alpha0 = 1.5 * alpha;
                            // if self.verbose:                    
                            //      print('alpha growes from %.3f to %.3f' % (alpha, alpha0))
                        }else{
                            alpha0 = alpha;
                            // if self.verbose:
                            //     print('alpha set to %.3f' % alpha0)
                        }
                    }
                    break;
                }else{
                    alpha = tau * alpha;
                }
            }           

            if (timeProfileInput){
               ros::Time tAfterAlpha = ros::Time::now();
               timeProfile.pushback((tAfterAlpha - tBeforeAlpha).toSec());
            }

            if (verbose){

            }
                // if alpha_found:
                //     print('We found alpha = %f' % alpha)
                // else:
                //     print('Fail to find alpha, use a conservative %f' % alpha)

            if (!alpha_found){
                converge_reason = -4;
                is_okay = true;
                converged = false;
                room_time = t_now;
                obj = obj0;
                if (logInput){
                    ros::Time ttemp = ros::Time::now();
                    double duration = (ttemp - t0).toSec();
                    Vector2d vtemp;
                    vtemp << obj0, duration;
                    log.pushback(vtemp);
                }
                break;
            }

            t_now = candid_time;
            // if self.verbose:
            //     print('obj0 = ', obj0, 'objf = ', objf)
            if (logInput){
                ros::Time ttemp = ros::Time::now();
                double duration = (ttemp - t0).toSec();
                vtemp << objf, duration;
                // if self.verbose == True:
                //     print("Logging: obj: %f, T: %f" % (objf, duration))
                log.pushback(vtemp);
            }       
            if (abs(objf - obj0) < abs_obj_tol){
                // if self.verbose:
                //     print('Absolute obj improvement too small')
                converged = true;
                converge_reason = 1;
                break;
            } else {
                if ( abs(objf - obj0) / abs(obj0) < rel_obj_tol){
                    // if self.verbose:
                    //     print('Relative obj improvement too small')
                    converged = true;
                    converge_reason = 2;
                    break;
                }
            }

        }
        major_iteration = i;
        ros::Time ttemp = ros::Time::now();
        double time_cost = (ttemp - t0).toSec();
    }
}

class IndoorQPProblem(IndoorOptProblem):
    """Formulate the indoor navigation problem explicitly as QP so we can either use mosek or osqp to solve it.

    I will manually maintain those matrices and hope it is more efficient.
    """
    def __init__(self, tgp, tfweight=0, connect_order=2, verbose=False):
        IndoorOptProblem.__init__(self, tgp, tfweight, connect_order, verbose)
        self.h_type = 'F'

    def update_prob(self):
        """Just update the problem since we are changing pretty fast.

        This function assume you might have change in room sequence so it reconstructs things. Take care with this.
        """
        self.construct_prob(self.x0_pack, self.xf_pack, self.poly_order, self.obj_order, self.connect_order)

    def construct_prob(self, x0_pack, xf_pack, poly_order, obj_order, connect_order):
        """Construct a problem."""
        # construct the problem using Fei's code
        self.construct_P()
        self.construct_A()

    def construct_P(self):
        # P is fixed and we do not alter it afterwards, so let's keep going
        pval, prow, pcol = construct_P(self.obj_order, self.num_box, self.poly_order, self.room_time, self.MQM, self.h_type)
        self.sp_P = coo_matrix((pval, (prow, pcol)))  # ugly hack since osqp only support upper triangular part or full
        self.n_var = self.sp_P.shape[0]
        # self.qp_P = spmatrix(sp_P.data, sp_P.row, sp_P.col)
        self.qp_q = np.zeros(self.n_var)

    def construct_A(self):
        lincon = construct_A(
                self.floor.getCorridor(),
                self.MQM,
                self.floor.position.copy(order='F'),
                self.floor.velocity.copy(order='F'),
                self.floor.acceleration.copy(order='F'),
                self.floor.maxVelocity,
                self.floor.maxAcceleration,
                self.floor.trajectoryOrder,
                self.floor.minimizeOrder,
                self.floor.margin,
                self.floor.doLimitVelocity,
                self.floor.doLimitAcceleration)
        self.xlb = lincon.xlb
        self.xub = lincon.xub
        self.clb = lincon.clb
        self.cub = lincon.cub
        # we need more
        self.sp_A = coo_matrix((lincon.aval, (lincon.arow, lincon.acol)))
        self.n_con = self.sp_A.shape[0]
        if self.verbose > 1:
            print('n_con', self.n_con)
            print('n_var', self.sp_A.shape[1])
            print("A has %d nnz" % lincon.aval.shape[0])

    def eval_cost_constr(self, mat_in):
        """Pass in a coefficient matrix, see results."""
        self.from_coef_matrix(mat_in)  # update self.sol
        self.update_prob()
        cost = 0.5 * self.sol.dot(self.sp_P.dot(self.sol))
        Ax = self.sp_A.dot(self.sol)
        # equality part
        error_clb = np.minimum(Ax - self.clb, 0)
        error_cub = np.minimum(-Ax + self.cub, 0)
        error_lb = np.minimum(self.sol - self.xlb, 0)
        error_ub = np.minimum(self.xub - self.sol, 0)
        with np.printoptions(precision=4, linewidth=10000):
            print('cost %f' % cost)
            print('error_eq', np.minimum(error_clb, error_cub))
            print('error_ieq', np.minimum(error_lb, error_ub))

    def get_gradient(self, sol, lmdy, lmdz):
        pgrad = gradient_from_P(self.obj_order, self.num_box, self.poly_order, self.room_time, self.MQM, sol)
        agrad = gradient_from_A(
                    self.floor.getCorridor(),
                    self.MQM,
                    self.floor.position.copy(order='F'),
                    self.floor.velocity.copy(order='F'),
                    self.floor.acceleration.copy(order='F'),
                    self.floor.maxVelocity,
                    self.floor.maxAcceleration,
                    self.floor.trajectoryOrder,
                    self.floor.minimizeOrder,
                    self.floor.margin,
                    self.floor.doLimitVelocity,
                    self.floor.doLimitAcceleration,
                    sol,
                    lmdy,
                    lmdz)
        if self.verbose > 1:
            print('pgrad', pgrad)
            print('agrad', agrad)
        return pgrad + agrad + self.tfweight

    def solve_once(self):
        """Solve the original problem once."""
        raise NotImplementedError

    def solve_with_room_time(self, rm_time):
        self.room_time[:] = rm_time
        self.floor.updateCorridorTime(self.room_time)
        return self.solve_once()


class IndoorQPProblemMOSEK(IndoorQPProblem):
    """Use mosek solver to solve this problem in cvxopt interface or not"""
    def __init__(self, tgp, tfweight=0, connect_order=2, verbose=False):
        IndoorQPProblem.__init__(self, tgp, tfweight, connect_order, verbose)
        self.h_type = "L"

    def solve_once(self):
        self.update_prob()
        # set up A
        A_sp = self.sp_A.tocsc()
        colptr, asub, acof = A_sp.indptr, A_sp.indices, A_sp.data
        aptrb, aptre = colptr[:-1], colptr[1:]
        # set up bounds on x
        bkx = self.n_var * [mosek.boundkey.ra]
        bkc = self.n_con * [mosek.boundkey.ra]
        with mosek.Env() as env:
            with env.Task(0, 1) as task:
                task.inputdata(self.n_con, self.n_var, self.qp_q.tolist(), 0.0,
                                list(aptrb), list(aptre), list(asub), list(acof),
                                bkc, self.clb.tolist(), self.cub.tolist(), 
                                bkx, self.xlb.tolist(), self.xub.tolist()
                                )
                # set up lower triangular part of P
                task.putqobj(self.sp_P.row.tolist(), self.sp_P.col.tolist(), self.sp_P.data.tolist())
                task.putobjsense(mosek.objsense.minimize)
                task.optimize()
                solsta = task.getsolsta(mosek.soltype.itr)
                x = self.n_var * [0.0]
                task.getsolutionslice(mosek.soltype.itr, mosek.solitem.xx, 0, self.n_var, x)
                x = np.array(x)
                # get dual variables on linear constraints
                zu, zl = self.n_con * [0.0], self.n_con * [0.0]
                task.getsolutionslice(mosek.soltype.itr, mosek.solitem.suc, 0, self.n_con, zu)
                task.getsolutionslice(mosek.soltype.itr, mosek.solitem.slc, 0, self.n_con, zl)
                z = np.array(zu) - np.array(zl)
                # get dual variables on variable bounds
                yu, yl = self.n_var * [0.0], self.n_var * [0.0]
                task.getsolutionslice(mosek.soltype.itr, mosek.solitem.sux, 0, self.n_var, yu)
                task.getsolutionslice(mosek.soltype.itr, mosek.solitem.slx, 0, self.n_var, yl)
                y = np.array(yu) - np.array(yl)
                if self.verbose:
                    print("Solving status", solsta)
                if solsta == mosek.solsta.optimal: #solsta == mosek.solsta.near_optimal: near_optimal is longer valid in Mosek 9.0
                    self.is_solved = True
                    self.obj = task.getprimalobj(mosek.soltype.itr) + self.tfweight * np.sum(self.room_time)
                    self.sol = x
                    self.lmdy = z
                    self.lmdz = y
                    return solsta, x, z, y
                else:
                    self.is_solved = False
                    self.obj = np.inf
                    #print("Mosek Failed, solsta: ", solsta)
                    return solsta, None, None, None

    def get_gradient(self):
        return IndoorQPProblem.get_gradient(self, self.sol, self.lmdy, self.lmdz)


def solveProblem():
    """Test the backtrack line search with IP solver."""
    prob = 11
    if len(sys.argv) > 1:
        prob = int(sys.argv[1])

    print_purple("Testing on problem: %d" % prob)

    tgp = loadTGP('dataset/tgp_%d.tgp' % prob)  # every time you have to reload from hard disk since it is modified before
    initial_time_allocation = np.array([box.t for box in tgp.getCorridor()])
    
    # we do not limit velocity in the open source implementation
    tgp.doLimitVelocity = False
    # we minimize jerk
    tgp.minimizeOrder = 3
    # we use 6-th order piecewise Bezier spline
    tgp.trajectoryOrder = 6
    
    #print_green("Use Mosek + Adaptive line search")
    
    solver = IndoorQPProblemMOSEK(tgp, verbose=False)
    
    ts1 = time.time()
    solver.solve_once()
    tf1 = time.time()
    initial_obj = solver.obj
    mosek_once_time = tf1 - ts1
   
    t_before_opt, coeff_before_opt = solver.get_output_coefficients()
    solver.grad_method = 'ours'

    # extract initial trajectory
    poly_coef = solver.get_coef_matrix().transpose((1,0,2))
    #import pdb; pdb.set_trace()
    break_points = np.insert(np.cumsum(solver.room_time), 0, 0.0)
    initial_trajectory = BPoly(poly_coef, break_points)
    tt_initial = np.linspace(0.0, break_points[-1], 100)
    
    ts2 = time.time()
    is_okay, converged = solver.refine_time_by_backtrack(max_iter=100, log=True, adaptiveLineSearch=True)
    tf2 = time.time()

    #print(solver.room_time)

    computation_time = (tf2 - ts2 + mosek_once_time) * 1000


    #print('mosek first computation time', mosek_once_time, 'mosek computation time:', mosek_time)
    mosek_convergence = solver.log
    mosek_convergence[1::2] += mosek_once_time

    t_after_opt, coeff_after_opt = solver.get_output_coefficients()

    result = [["Solver", "Mosek"], 
                ["Solved?", is_okay],
                #["Converged?", converged],
                ["Initial Cost", round(initial_obj, 3)],
                ["Final Cost", round(solver.obj,3)], 
                ["Solve Time [ms]", round(computation_time, 2)],
                ["# Major Iterations", solver.major_iteration], 
                ["# Function Evaluation", solver.num_prob_solve]
                ]
    print("Results")
    print(tabulate(result, tablefmt="psql", stralign="right", numalign="center"))
    
    print("Initial time allocation:\n", np.round(initial_time_allocation, 2))
    print("Final time allocation:\n", np.round(solver.room_time, 2))

    if DO_PLOT_RESULTS:
        poly_coef = solver.get_coef_matrix2()
        break_points = np.insert(np.cumsum(t_after_opt), 0, 0.0)
        final_trajectory = BPoly(poly_coef, break_points)
        tt_final = np.linspace(0.0, break_points[-1], 100)
        
        # plot 3D trajectory
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        boxes = tgp.getCorridor()
        
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        for i in range(len(boxes)):
            vertices = boxes[i].vertex
            ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], s=1)
            Z = vertices            
            verts = [[Z[0],Z[1],Z[2],Z[3]],
                     [Z[4],Z[5],Z[6],Z[7]], 
                     [Z[0],Z[1],Z[5],Z[4]], 
                     [Z[2],Z[3],Z[7],Z[6]], 
                     [Z[1],Z[2],Z[6],Z[5]],
                     [Z[4],Z[7],Z[3],Z[0]], 
                     [Z[2],Z[3],Z[7],Z[6]]]

            # plot safe corridor
            pc = Poly3DCollection(verts, alpha = 0.0, facecolor='gray', linewidths=0.1, edgecolors='red')
            ax.add_collection3d(pc)

        ax.plot(initial_trajectory(tt_initial)[:,0], initial_trajectory(tt_initial)[:,1], initial_trajectory(tt_initial)[:,2], label="Before refinement")
        ax.plot(final_trajectory(tt_final)[:,0], final_trajectory(tt_final)[:,1], final_trajectory(tt_final)[:,2], label="After refinement")
        #import pdb; pdb.set_trace()
        ax.scatter(tgp.position[0,0], tgp.position[0,1], tgp.position[0,2], marker="*", s=20, label="Start")
        ax.scatter(tgp.position[1,0], tgp.position[1,1], tgp.position[1,2], marker="o", s=20, label="Goal")
        set_axes_equal(ax)
        ax.set_axis_off()
        ax.legend()
        ax.set_title("3D trajectory")

        # plot velocity/acceleration
        fig, ax = plt.subplots(3, 1, figsize=(6,7))
        for i in range(3):
            ax[i].plot(tt_initial, initial_trajectory(tt_initial, 1)[:, i], '-.', label="Before refinement")
            ax[i].plot(tt_final, final_trajectory(tt_final, 1)[:, i], label="After refinement")
            if i == 0:
                ax[i].set_title("X Velocity")
            elif i == 1:
                ax[i].set_title("Y Velocity")
            elif i == 2:
                ax[i].set_title("Z Velocity")
            else:
                pass
            ax[i].legend()
            ax[i].grid()
        plt.tight_layout()
        

        # plot acceleration
        fig, ax = plt.subplots(3, 1, figsize=(6,7))
        for i in range(3):
            ax[i].plot(tt_initial, initial_trajectory(tt_initial, 2)[:, i], '-.', label="Before refinement")
            ax[i].plot(tt_final, final_trajectory(tt_final, 2)[:, i], label="After refinement")
            if i == 0:
                ax[i].set_title("X Acceleration")
            elif i == 1:
                ax[i].set_title("Y Acceleration")
            elif i == 2:
                ax[i].set_title("Z Acceleration")
            else:
                pass
            ax[i].legend()
            ax[i].grid()
        plt.tight_layout()



        plt.show()





def main():
    
    solveProblem()

