import numpy as np
import pysindy as ps
from pysindy.utils import AxesArray, comprehend_axes
from pysindy.feature_library.base import BaseFeatureLibrary
from pysindy.differentiation import FiniteDifference

def gl_weights(order, n):
    w = np.zeros(n)
    w[0] = 1.0
    for m in range(1, n):
        w[m] = w[m-1] * (1.0 - (order + 1.0) / m)
    return w

def gl_derivative_time(u, dt, order):
    """
    Computes temporal Grünwald-Letnikov fractional derivative along axis 0.
    u shape: (nt, nx, n_vars)
    """
    nt, nx, n_vars = u.shape
    w = gl_weights(order, nt)
    u_diff = np.zeros_like(u)
    for k in range(nt):
        w_slice = w[:k+1][::-1][:, np.newaxis, np.newaxis]
        u_diff[k] = np.sum(u[:k+1] * w_slice, axis=0)
    return u_diff / (dt**order)

def gl_derivative_space(u, dx, order):
    """
    Computes spatial left-sided Grünwald-Letnikov fractional derivative along axis 1.
    u shape: (nt, nx, n_vars) or (nx, n_vars)
    """
    if u.ndim == 3:
        nt, nx, n_vars = u.shape
        w = gl_weights(order, nx)
        u_diff = np.zeros_like(u)
        for i in range(nx):
            w_slice = w[:i+1][::-1][np.newaxis, :, np.newaxis]
            u_diff[:, i, :] = np.sum(u[:, :i+1, :] * w_slice, axis=1)
        return u_diff / (dx**order)
    else:
        nx, n_vars = u.shape
        w = gl_weights(order, nx)
        u_diff = np.zeros_like(u)
        for i in range(nx):
            w_slice = w[:i+1][::-1][:, np.newaxis]
            u_diff[i, :] = np.sum(u[:i+1, :] * w_slice, axis=0)
        return u_diff / (dx**order)

class FractionalPDELibrary(BaseFeatureLibrary):
    """
    A custom PySINDy feature library that generates:
    - Polynomials of states (up to degree 2)
    - Integer-order spatial derivatives (d=1, 2)
    - Fractional-order spatial derivatives (specified in beta_candidates)
    - Interaction terms (polynomials * derivatives)
    """
    def __init__(self, dx, beta_candidates, max_poly_degree=2, include_interaction=True):
        super().__init__()
        self.dx = dx
        self.beta_candidates = beta_candidates
        self.max_poly_degree = max_poly_degree
        self.include_interaction = include_interaction
        self.n_features_in_ = 3
        
    def fit(self, x, y=None):
        self.n_features_in_ = x[0].shape[-1] if isinstance(x, list) else x.shape[-1]
        names = self.get_feature_names()
        self.n_output_features_ = len(names)
        return self
        
    def get_feature_names(self, input_features=None):
        if input_features is None:
            input_features = ['p', 'c', 'n']
            
        names = []
        names.append('1')
        
        # Degree 1 polynomials
        for f in input_features:
            names.append(f)
            
        # Degree 2 polynomials
        if self.max_poly_degree >= 2:
            for i in range(len(input_features)):
                for j in range(i, len(input_features)):
                    names.append(f"{input_features[i]}_{input_features[j]}")
                    
        # Integer spatial derivatives
        deriv_names = []
        for f in input_features:
            deriv_names.append(f"{f}_x")
            deriv_names.append(f"{f}_xx")
            
        # Fractional spatial derivatives
        frac_names = []
        for beta in self.beta_candidates:
            beta_str = str(beta).replace('.', '_')
            for f in input_features:
                frac_names.append(f"D_{beta_str}_{f}")
                
        # Add derivatives to names
        for name in deriv_names + frac_names:
            names.append(name)
            
        # Add interaction terms (state * derivative)
        if self.include_interaction:
            for f in input_features:
                for d in deriv_names + frac_names:
                    names.append(f"{f}*{d}")
                    
        return names
        
    def transform(self, x_list):
        was_list = isinstance(x_list, list)
        if not was_list:
            x_list = [x_list]
            
        xp_list = []
        for x in x_list:
            if x.ndim == 2:
                nx, n_vars = x.shape
                x_3d = x[np.newaxis, :, :]
            elif x.ndim == 3:
                nt, nx, n_vars = x.shape
                x_3d = x
            else:
                raise ValueError("Expected input x to be 2D or 3D")
                
            nt_curr, nx_curr, n_vars_curr = x_3d.shape
            
            # Integer spatial derivatives (using finite difference along spatial axis=1)
            fd1 = FiniteDifference(d=1, axis=1)
            fd2 = FiniteDifference(d=2, axis=1)
            x_x = fd1._differentiate(x_3d, self.dx)
            x_xx = fd2._differentiate(x_3d, self.dx)
            
            # Fractional spatial derivatives
            frac_derivs = []
            for beta in self.beta_candidates:
                frac_derivs.append(gl_derivative_space(x_3d, self.dx, beta))
                
            # Construct feature matrix columns
            cols = []
            cols.append(np.ones((nt_curr, nx_curr, 1)))  # constant
            
            # Polynomials of degree 1
            for v_idx in range(n_vars_curr):
                cols.append(x_3d[..., v_idx:v_idx+1])
                
            # Polynomials of degree 2
            if self.max_poly_degree >= 2:
                for i in range(n_vars_curr):
                    for j in range(i, n_vars_curr):
                        cols.append((x_3d[..., i] * x_3d[..., j])[..., np.newaxis])
                        
            # All derivatives
            deriv_cols = []
            for v_idx in range(n_vars_curr):
                deriv_cols.append(x_x[..., v_idx:v_idx+1])
                deriv_cols.append(x_xx[..., v_idx:v_idx+1])
            for fd_arr in frac_derivs:
                for v_idx in range(n_vars_curr):
                    deriv_cols.append(fd_arr[..., v_idx:v_idx+1])
                    
            for col in deriv_cols:
                cols.append(col)
                
            # Interaction terms
            if self.include_interaction:
                for v_idx in range(n_vars_curr):
                    input_col = x_3d[..., v_idx:v_idx+1]
                    for d_col in deriv_cols:
                        cols.append(input_col * d_col)
                        
            xp = np.concatenate(cols, axis=-1)
            if x.ndim == 2:
                xp = xp[0]
                
            # Wrap in AxesArray for PySINDy compatibility
            xp = AxesArray(xp, comprehend_axes(xp))
            xp_list.append(xp)
            
        return xp_list if was_list else xp_list[0]

def build_fractional_pde_model(u, dt, dx, beta_x_candidates, x_dot, threshold=0.05):
    """
    Instantiates FractionalPDELibrary and fits SINDy using pre-computed temporal derivative x_dot.
    """
    library = FractionalPDELibrary(dx=dx, beta_candidates=beta_x_candidates, include_interaction=False)
    optimizer = ps.STLSQ(threshold=threshold, alpha=1e-4, normalize_columns=True)
    model = ps.SINDy(feature_library=library, optimizer=optimizer)
    model.fit(u, t=dt, x_dot=x_dot)
    return model, library
