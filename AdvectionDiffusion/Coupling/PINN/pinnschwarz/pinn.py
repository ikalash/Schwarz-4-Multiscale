import tensorflow as tf
import numpy as np

# A general class for constructing neural network architecture
class PINN_Architecture(tf.keras.Model):
    def __init__(
        self,
        xl=0,
        xr=1,
        num_hidden_layers=2,
        num_neurons_per_layer=20,
        output_dim=1,
        activation=tf.keras.activations.swish,
        kernel_initializer=tf.keras.initializers.GlorotNormal(seed=0),
        bias_initializer=tf.keras.initializers.RandomNormal(mean=0.4, stddev=0.1, seed=0),
        **kwargs,
    ):
        # Intialize superclass with its default parameter signature
        super().__init__(**kwargs)

        # Store model subdomain bounds
        self.xl = xl
        self.xr = xr

        # Initialize storage for BCs according to adjacent models
        self.u_gamma = np.array([0, 0])
        self.du_gamma = np.array([0, 0])

        # Store hyperparameters
        self.num_hidden_layers = num_hidden_layers
        self.output_dim = output_dim

        # Define NN architecture
        self.hidden = [
            tf.keras.layers.Dense(
                num_neurons_per_layer,
                activation=tf.keras.activations.get(activation),
                kernel_initializer=kernel_initializer,
            )  # ,
            # bias_initializer=bias_initializer)
            for _ in range(self.num_hidden_layers)
        ]
        self.out = tf.keras.layers.Dense(output_dim)

    # Mimic functionality of model(x)
    def call(self, X):
        # Forward-pass through neural network.
        Z = self.hidden[0](X)
        for i in range(1, self.num_hidden_layers):
            Z = self.hidden[i](Z)
        return self.out(Z)


# A general class for FD models on subdomains
class FD_1D_Steady:
    def __init__(self, X, BC, pde):
        self.X = X
        n_FD = len(X)
        xl = X[0]
        xr = X[-1]

        h = X[1] - X[0]

        nu = pde.nu
        beta = pde.beta
        order = pde.order

        a = -nu / (h**2)
        b = (2 * nu) / (h**2)
        c = -(nu / (h**2))

        if order == 1:
            b += beta / h
            c += -beta / h
            d = 0.0
        elif order == 2:
            b += 3 / 2 * beta / h
            c += -2 * beta / h
            d = 1 / 2 * beta / h
        else:
            raise ValueError(f"Invalid order: {order}")

        A = (
            np.diagflat([b] * (n_FD))
            + np.diagflat([c] * (n_FD - 1), -1)
            + np.diagflat([a] * (n_FD - 1), 1)
            + np.diagflat([d] * (n_FD - 2), -2)
        )

        if xr == BC[1]:
            y = np.ones((A.shape[0] - 1, 1))
            y[0] += -d * np.random.rand(1) - c * np.random.rand(1)
            y[1] += -d * np.random.rand(1)

            u_FD = np.linalg.solve(A[:-1, :-1], y)
            self.u = np.hstack((u_FD.flatten(), pde.f(xr)))

        elif xl == BC[0]:
            y = np.ones((A.shape[0] - 1, 1))
            y[-1] += -a * np.random.rand(1)

            u_FD = np.linalg.solve(A[:-1, :-1], y)
            self.u = np.hstack((pde.f(xl), u_FD.flatten()))

        else:
            y = np.ones((A.shape[0], 1))
            y[0] += -d * np.random.rand(1) - c * np.random.rand(1)
            y[1] += -d * np.random.rand(1)
            y[-1] += -a * np.random.rand(1)

            u_FD = np.linalg.solve(A, y)
            self.u = np.squeeze(u_FD)

        # save FD matrices for updating model
        self.A = A
        self.coeff = (a, b, c, d)

        # Initialize storage for BCs according to adjacent models
        self.u_gamma = np.array([0, 0])

    # Mimic functionality of model(x)
    def __call__(self, x):
        return np.interp(x, self.X, self.u)

class PINN_Schwarz_Steady():
    def __init__(self, pde, model_r, model_i, SDBC, X_r, X_b, alpha, snap, lamb_xb, BC_type, Schwarz_Iteration):
        # Store PDE
        self.pde = pde

        # Store models
        self.model_r = model_r
        self.model_i = model_i

        # Check for SDBC enforcement
        self.sdbc = SDBC

        # set BC treatment mixing boolean
        self.BC_type = BC_type

        # Store internal collocation points
        self.x = X_r

        # Store boundary points
        self.xb = X_b

        # store lambda
        self.lamb_xb = lamb_xb

        # Store Schwarz iteration
        self.sch_iter = Schwarz_Iteration

        # Store snapshot points if applicable
        if snap:
            npxs = np.linspace(float(model_r.xl), float(model_r.xr), num=snap, endpoint=False)[1:]
            self.xs = tf.constant(npxs, shape=(snap - 1, 1), dtype="float64")

        self.snap = snap

        # Store loss scaling coefficient
        self.a = alpha

        # scalar modifier for tanh enforcement functions
        self.m = 1

        # NN model loss and FD error storage
        self.loss = 0
        self.err = 0

    def get_u_hat(self, x, model):
        # Check if only schwarz boundaries are enforced
        if not self.sdbc[0]:
            bounds = [model.xl, model.xr]

            # Scaling functions to selectively enforce NN outputs to 0 at boundaries
            v_lr = [tf.math.tanh(self.m * (x - bounds[0])), tf.math.tanh(self.m * (bounds[1] - x))]

            # scaling functions to smooth enforcement of Schwarz BCs
            scale = [10 ** (-10 * (x - bounds[0])), 10 ** (10 * (x - bounds[1]))]

            # determine v_x based on boundary type of input model
            v_x = 1
            for i, b in enumerate(bounds):
                if not (b == 0 or b == 1):
                    v_x *= v_lr[i]

            # scale model outputs
            u_hat = v_x * model(x)

            # add BCs according to adjacent models
            for i, u_g in enumerate(model.u_gamma):
                u_hat = u_hat + scale[i] * u_g

            return u_hat

        # Check if only system boundaries are enforced
        elif not self.sdbc[1]:
            return tf.math.tanh(self.m * (1 - x)) * tf.math.tanh(self.m * x) * model(x)

        # Else both are strong boundaries
        else:
            bounds = [model.xl, model.xr]
            v_x = self.v_x(x, bounds)
            phi_x = self.phi_x(x, bounds)
            psi_x = self.psi_x(x, bounds)

            if self.BC_type == 'DD':
                u_hat = v_x*model(x) + phi_x*model.u_gamma[0] + psi_x*model.u_gamma[1]
            
            elif self.BC_type == 'DN':
                #  THESE SCALING FUNCTIONS NEED TO BE REDERIVED
                if  not self.model_i[1]:
                    # in final subdomain, pure dirichlet, 
                    # NOT THIS SIMPLE HOWEVER SINCE PHI AND PSI DONT SUPPORT RECOVERY OF Gi AND Gi+1  AT BOUNDARIES
                    u_hat = v_x*model(x) + phi_x*model.u_gamma[0] + psi_x*model.u_gamma[1]
                else:
                    u_hat = v_x*model(x) + phi_x*model.u_gamma[0] + psi_x*model.du_gamma[1]
            
            elif self.BC_type == 'RR':
                u_xt = self.u_x_theta(x, bounds, model)
                if bounds[0] == 0:
                    g = model.u_gamma[0]
                else:
                    g = model.u_gamma[1]
                    
                u_hat = v_x*model(x) + phi_x*g + psi_x*u_xt

        return u_hat

    # scaling functions
    def v_x(self, x, bounds):
        if self.BC_type == 'DD':
            v_x = tf.math.tanh(self.m*(x - bounds[0]))*tf.math.tanh(self.m*(bounds[1] - x))
        elif self.BC_type == 'DN':
            v_x = tf.math.exp((x - bounds[0]) * (x - bounds[1])) - 1 - (x - bounds[0]) * (x - bounds[1])
        elif self.BC_type == 'RR':
            phi1, phi2 = Robin_phis(bounds, x)
            # compute phi according to which sudomain is being solved (only 2 subdomains for RBC strong)
            if bounds[0] == 0:
                v_x = phi1 * (phi2**2)
            else:
                v_x = phi2 * (phi1**2)
        return v_x

    def phi_x(self, x, bounds):
        if self.BC_type == 'DD':
            phi_x = 10**(-10*(x-bounds[0]))
        elif self.BC_type == 'DN':
            phi_x = phi_transform(bounds, x)
        elif self.BC_type == 'RR':
            _, w2, _, w4 = Robin_ws(bounds, x)
            if bounds[0] == 0:
                phi_x = w2
            else:
                phi_x = w4

        return phi_x
    
    def psi_x(self, x, bounds):
        if self.BC_type == 'DD':
            psi_x = 10**(10*(x-bounds[1]))
        elif self.BC_type == 'DN':
            psi_x = psi_transform(bounds, x)
        elif self.BC_type == 'RR':
            w1, _, w3, _ = Robin_ws(bounds, x)
            if bounds[0] == 0:
                psi_x = w1
            else:
                psi_x = w3

        return psi_x
    
    def u_x_theta(self, x, bounds, model):
        # function to compute u(x;theta) for RBC depending on domain, only supports 2 subdomains
        
        phi1, phi2 = Robin_phis(bounds, x)
        c = 8 # user defined weight on Dirichlet term
        Du = self.gradient_NN(model, x)

        if bounds[0] == 0:
            Dn = 1 # directional derivative term
            h = Dn * model.du_gamma[1] + c * model.u_gamma[1]
            u = (1 + phi2 * (c + Dn * Du)) * model(x) - phi2 * h
        else:
            Dn = -1
            h = Dn * model.du_gamma[0] + c * model.u_gamma[0]
            u = (1 + phi1 * (c - Dn * Du)) * model(x) - phi1 * h
        return u

    def gradient_NN(self,model,x):
        # function to calculate the gradient of the nueral network

        with tf.GradientTape(persistent=True) as tape:
            # Watch variable x during this GradientTape
            tape.watch(x)

            # Compute current values u(x) with strongly enforced BCs
            u = model(x)

            # Store first derivative
            u_x = tape.gradient(u, x)

        return u_x
    
    @tf.function

    def get_gradients(self,model,x):
        # function to calculate the first and second derivative or either uhat or NN 

                # Compute current values u(x) with strongly enforced BCs
                if any(self.sdbc):
                    u = self.get_u_hat(x, self.model_r)
                else:
                    u = self.model_r(x)

            # Compute current values u(x) with strongly enforced BCs
            if any(self.sdbc):
                u = self.get_u_hat(x, model)
            else:
                u = model(x)

            # Store first derivative
            u_x = tape.gradient(u, x)

        # Store second derivative
        u_xx = tape.gradient(u_x, x)

        del tape
        return u_x, u_xx

    def get_residual(self,model, x):

        u_x, u_xx = self.get_gradients(model,x)

        return self.pde.f_r(u_x, u_xx)
        
    # Enforce system and/or Schwarz boundaries strongly
    def loss_strong(self, x):
        # Compute phi_r
        r = self.get_residual(self.model_r,x)
        phi_r = self.a * tf.reduce_mean(tf.square(r))

        # Initialize loss with residual loss function
        loss = phi_r

        phi_b = 0
        phi_i = 0

        # Check if only system boundaries are enforced
        if not self.sdbc[0]:
            # Calculate system boundary loss for current model if applicable
            for i,model in enumerate(self.model_i):

                # skip if has a neighboring subdomain
                if model:
                    continue

                b = self.xb[i]

                u_hat = self.get_u_hat(b, self.model_r)

                phi_b += (1 - self.a) * tf.reduce_mean(tf.square(self.pde.f_b(b) - u_hat))

        # Check if only Schwarz boundaries are enforced
        elif not self.sdbc[1]:
            # Calculate Schwarz boundary loss for current model if applicable
            for i,model in enumerate(self.model_i):

                # skip if is physical boundary
                if not model:
                    continue

                b = self.xb[i]

                u_pred1 = self.get_u_hat(b, self.model_r)

                if isinstance(model[0], FD_1D_Steady):
                    u_pred2 = model[0](b)
                else:
                    u_pred2 = self.get_u_hat(b, model[0])

                phi_i += (1 - self.a) * tf.reduce_mean(tf.square(u_pred1 - u_pred2))

        phi_s = 0
        if self.snap:
            # calculate snapshot data loss
            phi_s = (1 - self.a) * tf.reduce_mean(
                tf.square(self.get_u_hat(self.xs, self.model_r) - self.pde.f(self.xs))
            )

        # Add phi_b, phi_i, and phi_s to the loss
        loss += phi_b + phi_i + phi_s

        return loss, phi_r, phi_b, phi_i, phi_s
    
    def loss_weak(self, x):
        # this calculates the value of the loss function for Dirichlet-Neuman coupling with relaxation
        # Compute phi_r
        r = self.get_residual(self.model_r,x)
        phi_r = self.a * tf.reduce_mean(tf.square(r))

        # Initialize loss with residual loss function
        loss = phi_r

        phi_b = 0
        phi_i = 0
        
        # set indices of BCs (-10 for non-applicable BC types)
        Dir_pts = []
        Neum_pts = []
        Robin_pts = []
        if self.BC_type == 'DN':
            Dir_pts, Neum_pts = self.Dirch_Neuman_pts()
        elif self.BC_type == 'DD':
            Dir_pts = [0, 1]
        elif self.BC_type == 'RR':
            Robin_pts = [0, 1]

        # step through boundaries, enforcing corresponding BCs
        boundary_pt = -1 # counter
        for i, model in enumerate(self.model_i):
            boundary_pt += 1
            b = self.xb[i]

            # Calculate boundary loss for current model if applicable
            if not model:
                u_pred = self.model_r(b)
                phi_b += (1 - self.a) * tf.reduce_mean(tf.square(self.pde.f_b(b) - u_pred))

            else:
                # compute function and gradient values at boundary
                u_pred1 = self.model_r(b)
                if self.BC_type == 'DN':
                    u_pred2 = self.lamb_xb
                else:
                    u_pred2 = model[0](b)   

                # apply relevant BC enforcement
                if boundary_pt in Neum_pts:
                    # calculate interface loss for current model if applicable at front interface
                    # apply Neuman conditions on left
                    du_pred1, _ = self.get_gradients(self.model_r,b)
                    du_pred2, _ = self.get_gradients(model[0],b)
                    phi_i += (1 - self.a) * tf.reduce_mean(tf.square(du_pred1 - du_pred2))

                elif boundary_pt in Dir_pts:
                    # Calculate interface loss for current model if applicable at back interface of subdomain
                    # apply Dirichlet conditions
                    phi_i += (1 - self.a) * tf.reduce_mean(tf.square(u_pred1 - u_pred2))

                elif boundary_pt in Robin_pts:
                    # apply robin-robin conditions
                    L = 10
                    du_pred1, _ = self.get_gradients(self.model_r,b)
                    du_pred2, _ = self.get_gradients(model[0],b)
                    # compute normal directions
                    n1 = (-1)**(boundary_pt+1)
                    n2 = (-1)**(boundary_pt)
                    phi_i += (1 - self.a) * tf.reduce_mean(tf.square((L*u_pred1 + n1*du_pred1) - (L*u_pred2 + n2*du_pred2)))
        
        phi_s = 0
        if self.snap:
            # calculate snapshot data loss
            phi_s = (1 - self.a) * tf.reduce_mean(tf.square(self.model_r(self.xs) - self.pde.f(self.xs)))

        # Add phi_b, phi_i, and phi_s to the loss
        loss += phi_b + phi_i + phi_s
        
        return loss, phi_r, phi_b, phi_i, phi_s

    def Dirch_Neuman_pts(self):
        # function returns the x point at which the BC is dirchlet and at which it is neumann depending on iteration
        if self.sch_iter % 2 == 0:
        # even iterations
            Dir_index = [0]
            Neum_index = [1]
        else:
            Dir_index = [1]
            Neum_index = [0]

        return Dir_index, Neum_index
        
    @tf.function
    def get_gradient_trainable(self, x):
        with tf.GradientTape(persistent=True) as tape:

            # This tape is for derivatives with respect to trainable variables
            tape.watch(self.model_r.trainable_variables)
            if any(self.sdbc):
                loss, _, _, _, _ = self.loss_strong(x)
            else:
                loss, _, _, _, _ = self.loss_weak(x)

        g = tape.gradient(loss, self.model_r.trainable_variables)

        return g

    def FD_update(self):
        a, b, c, d = self.model_r.coeff
        model = self.model_i
        A = self.model_r.A

        if model[0] and model[1]:
            f_NN = np.ones((A.shape[0], 1))

            if any(self.sdbc) and not isinstance(model[0][0], FD_1D_Steady):
                f_NN[0] = f_NN[0] + (
                    -d * self.model_r.u_gamma[0] - c * self.get_u_hat(tf.reshape(self.x[1], shape=(1, 1)), model[0][0])
                )
                f_NN[1] = f_NN[1] + (-d * self.get_u_hat(tf.reshape(self.x[1], shape=(1, 1)), model[0][0]))

            else:
                f_NN[0] = f_NN[0] + (
                    -d * model[0][0](tf.reshape(self.x[0], shape=(1, 1)))
                    - c * model[0][0](tf.reshape(self.x[1], shape=(1, 1)))
                )
                f_NN[1] = f_NN[1] + (-d * model[0][0](tf.reshape(self.x[1], shape=(1, 1))))

            if any(self.sdbc) and not isinstance(model[1][0], FD_1D_Steady):
                f_NN[-1] = f_NN[-1] + (-a * self.model_r.u_gamma[1])
            else:
                f_NN[-1] = f_NN[-1] + (-a * model[1][0](tf.reshape(self.x[-1], shape=(1, 1))))

            u_FD = np.squeeze(np.linalg.solve(A, f_NN))

        elif model[1]:
            f_NN = np.ones((A.shape[0] - 1, 1))
            A = A[:-1, :-1]

            if any(self.sdbc) and not isinstance(model[1][0], FD_1D_Steady):
                f_NN[-1] = f_NN[-1] + (-a * self.model_r.u_gamma[1])
            else:
                f_NN[-1] = f_NN[-1] + (-a * model[1][0](tf.reshape(self.x[-1], shape=(1, 1))))

            u_FD = np.linalg.solve(A, f_NN)

            u_FD = np.hstack((self.pde.f(self.x[0]), u_FD.flatten()))

        elif model[0]:
            f_NN = np.ones((A.shape[0] - 1, 1))
            A = A[:-1, :-1]

            if any(self.sdbc) and not isinstance(model[0][0], FD_1D_Steady):
                f_NN[0] = f_NN[0] + (
                    -d * self.model_r.u_gamma[0] - c * self.get_u_hat(tf.reshape(self.x[1], shape=(1, 1)), model[0][0])
                )
                f_NN[1] = f_NN[1] + (-d * self.get_u_hat(tf.reshape(self.x[1], shape=(1, 1)), model[0][0]))

            else:
                f_NN[0] = f_NN[0] + (
                    -d * model[0][0](tf.reshape(self.x[0], shape=(1, 1)))
                    - c * model[0][0](tf.reshape(self.x[1], shape=(1, 1)))
                )
                f_NN[1] = f_NN[1] + (-d * model[0][0](tf.reshape(self.x[1], shape=(1, 1))))

            u_FD = np.linalg.solve(A, f_NN)

            u_FD = np.hstack((u_FD.flatten(), self.pde.f(self.x[-1])))

        else:
            u_FD = np.squeeze(self.pde.f(self.x))

        # Update u for current model
        self.model_r.u = u_FD

    def solve(self, optimizer, numEpochs):
        @tf.function
        def train_step(x):
            # Retrieve loss gradient w.r.t. trainable variables
            grad_theta = self.get_gradient_trainable(x)

            # Perform gradient descent step
            optimizer.apply_gradients(zip(grad_theta, self.model_r.trainable_variables))

        # Split data into training batches
        # train_dataset = tf.data.Dataset.from_tensor_slices((self.x,))
        # train_dataset = train_dataset.shuffle(buffer_size=1024).batch(batch_size)

        # If current model is FOM, update interface boundaries with adjacent NN models
        if isinstance(self.model_r, FD_1D_Steady):
            self.FD_update()
            self.err = np.square(self.model_r.u - self.pde.f(self.x)).mean()
        else:
            # Iterate training
            for i in range(numEpochs):
                # Train on each batch
                # for (x_batch_train,) in train_dataset:
                # train_step(x_batch_train)
                train_step(self.x)

            # Compute loss for full dataset to track training progress
            if any(self.sdbc):
                self.loss, self.phi_r, self.phi_b, self.phi_i, self.phi_s = self.loss_strong(self.x)
            else:
                self.loss, self.phi_r, self.phi_b, self.phi_i, self.phi_s = self.loss_weak(self.x)

    # helper functions

def phi_transform(bounds,x):
     # this function computes the coefficients a1-a5 for the phi transformation function to strongly enforce Dirichlet-Neuman BC
     # input is the domain bounds in a single list and the output is coeffs a1-a5 in a single list

    # set gamma_i and gamma_i+1
    gi = bounds[0]
    gip1 = bounds[1]

    a1 = -1/2*(3*gi*gip1**2 + gi**3 + 2 - 3*gip1*gi**2 - gip1**3) / (gi**4 + 2*gi*gip1**3 - gip1**4 - 2*gip1*gi**3)
    a2 = 1
    a3 = -1/2*(gi**5 + gip1*gi**4 - 8*gip1**2*gi**3 + 8*gi**2*gip1**3 - 4*gi**2 - gi*gip1**4 - 4*gi*gip1 - gip1**5 - 4*gip1**2) / (gi**4 + 2*gi*gip1**3 - gip1**4 - 2*gip1*gi**3)
    a4 = (gi**3 - 3*gip1*gi**2 + 3*gi*gip1**2 - 4 - gip1**3) * gi*gip1 / (gi**3 - 3*gip1*gi**2 + 3*gi*gip1**2 - gip1**3)
    a5 = -1/2*gip1**2*(3*gip1**2*gi**3 + 2*gip1**2 - gi**2*gip1**3 - 3*gip1*gi**4 + gi**5 - 4*gi**2 - 4*gi*gip1) / (gi**4 + 2*gi*gip1**3 - gip1**4 - 2*gip1*gi**3);
    
    phi = a1*(x**4) + a2*(x**3) + a3*(x**2) + a4*x + a5

    return phi

def psi_transform(bounds,x):
    # this fxn computes the value of the psi transform function for a specific x value
    
    # set gamma_i and gamma_i+1
    gi = bounds[0]
    gip1 = bounds[1]

    b1 = -1/2*(gi**2 - 2*gip1*gi - 1 + gip1**2) / (gi**3 - gip1*gi**2 - gi*gip1**2 + gip1**3)
    b2 = 1
    b3 = -1/2*(gi**4 + 2*gip1*gi**3 + 3*gi**2 - 6*gi**2*gip1**2 + 2*gi*gip1**3 + 2*gip1*gi + gip1**4 + gip1**2) / (gi**3 - gip1*gi**2 - gi*gip1**2 +gip1**3)
    b4 = (gip1*gi**2 - 2*gi*gip1**2 + gi + gip1 + gip1**3)*gi / (gip1**2 - 2*gip1*gi + gi**2)
    b5 = -1/2*gip1*gi**2*(gip1*gi**2 - 2*gi*gip1**2 + 2*gi + gip1**3 +gip1) / (gi**3 - gip1*gi**2 - gi*gip1**2 + gip1**3) 
    
    psi = b1*(x**4) + b2*(x**3) + b3*(x**2) + b4*x + b5
    
    return psi

def Robin_phis(bounds, x):
    # helper function to compute phi scaling functions for RBC strong enforcement
    
    phi1 = x - bounds[0]
    phi2 = bounds[1] - x

    return phi1, phi2

def Robin_ws(bounds, x):
    # helper function to compute the w's for RBC strong enforcement
    phi1, phi2 = Robin_phis(bounds, x)

    w1 = phi1 / (phi1 + phi2**2)
    w2 = phi2**2 / (phi1 + phi2**2)
    w3 = phi2 / (phi2 +phi1)
    w4 = phi1 / (phi2 + phi1**2)

    return w1, w2, w3, w4