

class IndoorOptProblem(object):
    """Convex optimization around corridor approach for drone trajectory optimization.

    Currently it support problem where each corridor has been specified a time, later we will let it be free.
    """
    def __init__(self, tgp, tfweight=0, connect_order=2, verbose=False):
        """Constructor for this problem class.

        Parameters
        ----------
        tgp: a PyTGProblem object.
        tfweight: float, the weight on transfer time
        connect_order: to which order of derivative at connection do we guarantee continuity.
        verbose: bool, if the solver is verbose
        """
        self.floor = tgp  # I shall still use the floor name for convenience
        self.x0_pack = [tgp.position[0], tgp.velocity[0], tgp.acceleration[0]]
        self.xf_pack = [tgp.position[1], tgp.velocity[1], tgp.acceleration[1]]
        self.tfweight = tfweight
        self.obj_order = tgp.minimizeOrder
        self.poly_order = tgp.trajectoryOrder
        self.connect_order = connect_order  # guarantee acceleration continuity
        self.verbose = verbose
        self.margin = tgp.margin
        self.is_limit_vel = tgp.doLimitVelocity
        self.is_limit_acc = tgp.doLimitAcceleration
        self.vel_limit = tgp.maxVelocity  # component wise limit on velocity
        self.acc_limit = tgp.maxAcceleration  # component wise limit on acceleration
        # get info on environment
        self.boxes = tgp.getCorridor()
        self.num_box = len(self.boxes)
        self.room_time = np.array([box.t for box in self.boxes])
        bz = Bezier(self.poly_order, self.poly_order, self.obj_order)
        self.bz = bz
        self.bzM = self.bz.M()[self.poly_order]  # use this for some output stuff
        self.MQM = self.bz.MQM()[self.poly_order]
        if verbose:
            print('has %d boxes' % self.num_box)
            print('init time ', self.room_time)
        # some default settings such as how close you can reach the corner
        self.abs_obj_tol = 1e-3  # this controls absolute objective tolerance and should not set small for problem with small obj
        self.rel_obj_tol = 1e-3  # this controls relative objective tolerance
        self.grad_tol = 1e-3  # this controls gradient tolerance

    def set_tfweight(self, weight):
        """Set weight on time"""
        self.tfweight = weight

    def construct_prob(self, x0_pack, xf_pack, poly_order, obj_order, connect_order):
        """Construct a problem."""
        pass

    def set_x0_pack_value(self, *args):
        """Set the contents of x0pack"""
        ff = zip(self.x0_pack, args)
        print(ff)
        for tmp, val in zip(self.x0_pack, args):
            tmp[:] = val[:]

    def set_xf_pack_value(self, *args):
        """Set the contents of xfpack"""
        for tmp, val in zip(self.xf_pack, args):
            tmp[:] = val[:]

    def solve_with_room_time(self, rm_time):
        """Specify room time and solve the problem"""
        raise NotImplementedError("Subclass should implement solve_with_room_time function")

    def get_output_coefficients(self, ):
        """
        Get the monomial coefficients representing a piece-wise polynomial trajectory
        return ndarray, (s, 1), time allocated for each segment 
        return ndarray, (s, o + 1, d), polynomial coefficients, where s is the number of segments, 
        o is the order of trajectory, d is the number of dimensions, which is 3
        """
        
        n_room = self.num_box
        
        poly_coef = np.zeros([n_room, self.poly_order + 1, 3])
        
        # coefficients in bezier and scaled form
        mat_x = np.reshape(self.sol, (n_room, 3, self.poly_order + 1))
        mat_x = np.transpose(mat_x, (0, 2, 1))

        # change coefficients to monomial and unscaled form
        for i in range(n_room):
            poly_coef[i, :, :] = self.bzM.dot(mat_x[i]) * self.room_time[i]

        return self.room_time.copy(), poly_coef

    def get_coef_matrix(self):
        """Return coefficients"""
        coef_mat = np.reshape(self.sol, (self.num_box, 3, self.poly_order + 1))
        coef_mat = np.transpose(coef_mat, (0, 2, 1))
        #import pdb; pdb.set_trace()
        for i in range(self.num_box):
            coef_mat[i] *= self.room_time[i]
        return coef_mat

    def get_coef_matrix2(self):
        """
        Return Bezier coefficients
        """
        coef_mat = np.reshape(self.sol, (self.poly_order + 1, 3, self.num_box), order='F')
        coef_mat = coef_mat.transpose((0,2,1))
        #import pdb;pdb.set_trace()
        for i in range(self.num_box):
            coef_mat[:,i,:] *= self.room_time[i]
        
        return coef_mat


    def from_coef_matrix(self, mat_in):
        """Assign values to sol based on the input coefficient matrix."""
        self.sol = np.zeros(mat_in.size)
        coef_mat = np.reshape(self.sol, (self.num_box, 3, self.poly_order + 1))
        for i in range(self.num_box):
            coef_mat[i] = mat_in[i].T / self.room_time[i]
        with np.printoptions(precision=4, linewidth=10000):
            print(self.sol)

    def get_output_path(self, n):
        """Get a path for output that is linspace in time.

        :param n: int, the number of nodes for a path
        :return: float, the total time for this problem
        :return: ndarray, (n, 2) the optimal path
        """
        cum_sum_time = np.cumsum(self.room_time)
        output = np.zeros((n, 3))  # the output trajectory
        sample_time = np.linspace(0, cum_sum_time[-1], n)
        n_room = self.num_box

        # get all the coef of polynomials
        t, all_poly_coeffs = self.get_output_coefficients()
        for i in range(n_room):
            # poly_coef = self.bzM.dot(mat_x[i]) * self.room_time[i]
            poly_coef = all_poly_coeffs[i,:,:]
            if i == 0:
                t_mask = sample_time <= cum_sum_time[0]
                use_s = sample_time[t_mask] / cum_sum_time[0]
            else:
                t_mask = (sample_time > cum_sum_time[i - 1]) & (sample_time <= cum_sum_time[i])
                use_s = (sample_time[t_mask] - cum_sum_time[i - 1]) / self.room_time[i]
            output[t_mask, 0] = np.polyval(poly_coef[:, 0][::-1], use_s)
            output[t_mask, 1] = np.polyval(poly_coef[:, 1][::-1], use_s)
            output[t_mask, 2] = np.polyval(poly_coef[:, 2][::-1], use_s)
        return cum_sum_time[-1], output

    def get_gradient(self):
        raise NotImplementedError

    def get_gradient_fd(self, h=1e-6):
        """Use forward finite difference to approximate gradients."""
        grad = np.zeros(self.num_box)
        obj0 = self.obj
        origin_time = self.room_time.copy()
        for i in range(self.num_box):
            try_time = origin_time.copy()
            try_time[i] += h
            self.solve_with_room_time(try_time)
            grad[i] = (self.obj - obj0) / h
        grad += self.tfweight
        return grad

    def get_gradient_mellinger(self, h=1e-6):
        """
        Use finite difference described in:
        http://www-personal.acfr.usyd.edu.au/spns/cdm/papers/Mellinger.pdf
        to approximate gradients.
        """
        grad = np.zeros(self.num_box)
        obj0 = self.obj
        origin_time = self.room_time.copy()
        for i in range(self.num_box):
            m = self.num_box
            gi = -1/(m-1) * np.ones(m)
            gi[i] = 1
            try_time = origin_time.copy()
            try_time += h * gi
            # print("In gg Mellinger: gi", gi," try_time: ", try_time)
            self.solve_with_room_time(try_time)
            grad[i] = (self.obj - obj0) / h
        grad += self.tfweight
        return grad

    def refine_time_by_backtrack(self, alpha0=0.175, h=1e-5, c=0.2, tau=0.2, max_iter=50, j_iter=5, log=False, timeProfile=False, adaptiveLineSearch=False):
        """Use backtrack line search to refine time. We fix total time to make things easier.

        Parameters
        ----------
        alpha0: float, initial step length 0.175 and 0.375 are found to be very good.
        h: float, step size for finding gradient using forward differentiation
        c: float, the objective decrease parameter
        tau: float, the step length shrink parameter
        max_iter: int, maximum iteration for gradient descent
        j_iter: int, maximum iteration for finding alpha
        abs_tol: float, absolute objective tolerance
        rel_tol: float, Relative objective tolerance
        Returns
        -------
        is_okay: bool, indicates if some exception occurs
        converged: bool, indicates if the algorithm converges
        """

        if log == True:
            self.log = np.array([])
            self.log = np.append(self.log, [self.obj, 0])

        if timeProfile == True:
            self.timeProfile = np.array([])

        t0 = time.time()
        self.major_iteration = 0
        self.num_prob_solve = 0
        if self.num_box == 1 and self.tfweight == 0:
            self.major_iteration = 0
            self.num_prob_solve = 0
            self.time_cost = time.time() - t0
            self.converge_reason = 'No need to refine'
            return True, True
        n_room = self.num_box
        t_now = self.room_time.copy()

        converged = False
        converge_reason = 'Not converged'
        num_prob_solve = 0  # record number of problems being solved for a BTLS
        for i in range(max_iter):
            if self.verbose:
                print_green('Iteration %d' % i)
           
            obj0 = self.obj
            is_okay = True
            
            if timeProfile == True:
                tBeforeGrad = time.time()

            # choose a method to calculate gradient
            if self.grad_method == 'ours':
                grad = self.get_gradient()
            elif self.grad_method == 'fd':
                grad = self.get_gradient_fd(h=0.25*1e-6)
            elif self.grad_method == 'mel':
                grad = self.get_gradient_mellinger(h=0.25*1e-6)
            else:
                print_red("No this grad method!")
            
            if timeProfile == True:
                tAfterGrad = time.time()
                self.timeProfile = np.append(self.timeProfile, [tAfterGrad - tBeforeGrad])

            #print_red('grad_an ', grad, ' grad_fd ', grad_fd, 'grad_mel', grad_mel)
            if self.tfweight == 0:
                # get projected gradient, the linear manifold is \sum x_i = 0; if tfweight=0, we fix total time
                normal = np.ones(n_room) / np.sqrt(n_room)  # normal direction
                grad = grad - grad.dot(normal) * normal  # projected gradient descent direction
            # print_red('grad_an ', grad / np.linalg.norm(grad), ' grad_fd ', grad_fd / np.linalg.norm(grad_fd), 'grad_mel', grad_mel / np.linalg.norm(grad_mel))
            # print_yellow('sum_an ', np.sum(grad), 'sum_fd ', np.sum(grad_fd), 'sum_mel ', np.sum(grad_mel))
            if np.linalg.norm(grad) < self.grad_tol:  # break out if gradient is too small
                if self.verbose:
                    print('Gradient too small')
                converged = True
                converge_reason = 'Small gradient'
                break
            if self.verbose:
                print_green('At time ', t_now, ' grad is ', grad)
            m = -np.linalg.norm(grad)
            p = grad / m  # p is the descending direction
            # use a maximum alpha that makes sure time are always positive
            alpha_max = np.amax(-t_now / p) - 1e-6  # so I still get non-zero things
            if alpha_max > 0:
                alpha = min(alpha_max, alpha0)
            else:
                alpha = alpha0
            t = -c * m

            # find alpha
            alpha_found = False

            if timeProfile == True:
                tBeforeAlpha = time.time()

            for j in range(j_iter):
                if self.verbose:
                    print('Search alpha step %d, alpha = %f' % (j, alpha))
                candid_time = t_now + alpha * p

                # lower bound on the alpha
                if adaptiveLineSearch == True:
                    if alpha < 1e-4:
                        if self.verbose:
                            print_yellow('Stop line search because alpha is too small')
                        break
                
                # make sure that time will not go too small
                if np.any(candid_time < 1e-6):
                    alpha = tau * alpha
                    continue
                if self.verbose:
                    print('Try time: ', candid_time)
                self.solve_with_room_time(candid_time)
                num_prob_solve += 1

                # in case objective somehow falls below 0
                if self.obj < 0:
                    self.is_solved = False
                if not self.is_solved:
                    alpha = tau * alpha  # decrease step length
                    continue                
                objf = self.obj
                
                # in case objective somehow falls below 0
                if objf < 0:
                    print_yellow("Negative Objective!", self.obj)
                    print(self.solve_with_room_time(candid_time))
                    import pdb; pdb.set_trace()
                    alpha = tau * alpha
                    continue

                if self.verbose:
                    print('\talpha ', alpha, ' obj0 ', obj0, ' objf ', objf)
                if obj0 - objf >= alpha * t or obj0 - objf >= 0.1 * obj0:  # either backtrack or decrease sufficiently
                    alpha_found = True
                    # if use adaptive line search
                    if adaptiveLineSearch == True:
                        # adaptive line search, increase the initial alpha is alpha is good enough for the first ime
                        if j == 0:
                            alpha0 = 1.5 * alpha
                            if self.verbose:                    
                                 print('alpha growes from %.3f to %.3f' % (alpha, alpha0))
                        # set the initial alpha for the next iteration the same as this time's alpha 
                        else:
                            alpha0 = alpha
                            if self.verbose:
                                print('alpha set to %.3f' % alpha0)
                    break
                else:
                    alpha = tau * alpha  # decrease step length
                

            if timeProfile == True:
                tAfterAlpha = time.time()
                self.timeProfile = np.append(self.timeProfile, [tAfterAlpha - tBeforeAlpha])

            if self.verbose:
                if alpha_found:
                    print('We found alpha = %f' % alpha)
                else:
                    print('Fail to find alpha, use a conservative %f' % alpha)
            if not alpha_found:  # for case where alpha is not found
                converge_reason = 'Cannot find step size alpha'
                is_okay = True
                converged = False
                # roll back to t_now
                self.room_time = t_now
                self.obj = obj0
                if log == True:
                    duration = time.time() - t0
                    self.log = np.append(self.log, [obj0, duration])
                break

            # adaptive line search

            # ready to update time now and check convergence
            t_now = candid_time  # this is the alpha we desire
            if self.verbose:
                print('obj0 = ', obj0, 'objf = ', objf)
            if log == True:
                duration = time.time() - t0
                # we log the objective and duration for each major iteration
                if self.verbose == True:
                    print("Logging: obj: %f, T: %f" % (objf, duration))
                self.log = np.append(self.log, [objf, duration])
            if abs(objf - obj0) < self.abs_obj_tol:
                if self.verbose:
                    print('Absolute obj improvement too small')
                converged = True
                converge_reason = 'Absolute cost'
                break
            elif abs(objf - obj0) / abs(obj0) < self.rel_obj_tol:
                if self.verbose:
                    print('Relative obj improvement too small')
                converged = True
                converge_reason = 'Relative cost'
                break
        self.major_iteration = i
        self.num_prob_solve = num_prob_solve
        self.time_cost = time.time() - t0
        self.converge_reason = converge_reason
        return is_okay, converged



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

