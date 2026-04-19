import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import digamma
from types import SimpleNamespace
from tqdm import tqdm

# =====================================================================
# FUTURE IMPORTS:
# Once translated, you will import your missing functions here so 
# cal_h and cal_jpdf_hist can use them.
# from histcn_py import histcn
# from design_B3_py import design_B3
# =====================================================================

def parList(Opts, RandV, isNorm):
    """Assign distribution to different input parameters."""
    dist_v = Opts.distType
    nVar = RandV.nVar
    vNominal = RandV.vNominal.flatten()
    CoV = RandV.CoV.flatten()

    J = np.eye(nVar * 2, nVar * 2) 

    if isinstance(dist_v, str):
        dist_v = [dist_v] * nVar

    ListPar = np.zeros((nVar, 6)) 
    
    for ii in range(nVar):
        dist = dist_v[ii]
        parMean = vNominal[ii]
        parStd = vNominal[ii] * CoV[ii]

        if parMean == 0:
            parStd = CoV[ii]

        if dist.lower() == 'normal':
            distTypeMultiplier = 1
            parA = parMean
            parB = parStd

        elif dist.lower() == 'lognormal':
            distTypeMultiplier = 2
            mu = np.log(parMean**2 / np.sqrt(parMean**2 + parStd**2))
            sigma = np.sqrt(np.log(1 + parStd**2 / parMean**2))
            parA = mu
            parB = sigma

        elif dist.lower() == 'gamma':
            distTypeMultiplier = 3
            k = (parMean / parStd)**2        
            theta = parStd**2 / parMean      
            parA = k
            parB = theta

        elif dist.lower() == 'gumbel':
            distTypeMultiplier = 4
            mu = parMean - np.sqrt(6) * np.euler_gamma * parStd / np.pi 
            sigma = np.sqrt(6) * parStd / np.pi 
            parA = mu
            parB = sigma
        else:
            raise ValueError(f"Unknown distribution: {dist}")

        if isNorm in [2, 3]:
            J = getJ(parMean, parStd, nVar, dist, ii, J)

        ListPar[ii, :] = [ii + 1, vNominal[ii], CoV[ii], distTypeMultiplier, parA, parB]

    parC = np.full((nVar, 1), np.nan)
    parD = np.full((nVar, 1), np.nan)

    ListPar = np.hstack((ListPar, parC, parD))

    return ListPar, J


def getJ(m, s, nPar, dist, index, J):
    """Reparameterize Fisher with respect to mean (m) and standard deviation (s)."""
    dist = dist.lower()
    
    if dist == 'normal':
        mu_m, mu_s = 1, 0
        sigma_m, sigma_s = 0, 1

    elif dist == 'gamma':
        mu_m = 2 * m / s**2
        mu_s = -2 * m**2 / s**3
        sigma_m = -s**2 / m**2
        sigma_s = 2 * s / m

    elif dist == 'lognormal':
        mu_m = 2 / m - m / (m**2 + s**2)
        mu_s = -s / (m**2 + s**2)
        sigma_m = (np.log(1 + s**2 / m**2))**(-0.5) * (-s**2 / m) / (m**2 + s**2)
        sigma_s = (np.log(1 + s**2 / m**2))**(-0.5) * s / (m**2 + s**2)

    elif dist == 'gumbel':
        mu_m = 1
        mu_s = -np.sqrt(6) * np.euler_gamma / np.pi
        sigma_m = 0
        sigma_s = np.sqrt(6) / np.pi

    Jrepa = J.copy()

    Jrepa[index, index] = mu_m
    Jrepa[index + nPar, index + nPar] = sigma_s
    Jrepa[index, index + nPar] = mu_s
    Jrepa[index + nPar, index] = sigma_m

    return Jrepa


def parSampling(ListPar, nPar, Opts):
    """Generate Monte Carlo samples and parameter sensitivities."""
    nS = Opts.nSampMC

    xS = SimpleNamespace()
    xS.samp = np.zeros((nS, nPar))
    xS.senA = np.zeros((nS, nPar))
    xS.senB = np.zeros((nS, nPar))
    xS.senC = np.zeros((nS, nPar))
    xS.senD = np.zeros((nS, nPar))

    ParSen = [[None, None] for _ in range(nPar)]

    for ii in range(nPar):
        Par = ListPar[ii, :]
        dist_type = Par[3]
        
        if dist_type == 1:
            dist = 'Normal'
            samp = stats.norm.rvs(loc=Par[4], scale=Par[5], size=nS)
        elif dist_type == 2:
            dist = 'Lognormal'
            samp = stats.lognorm.rvs(s=Par[5], scale=np.exp(Par[4]), size=nS)
        elif dist_type == 3:
            dist = 'Gamma'
            samp = stats.gamma.rvs(a=Par[4], scale=Par[5], size=nS)
        elif dist_type == 4:
            dist = 'gev'
            samp = stats.gumbel_r.rvs(loc=Par[4], scale=Par[5], size=nS)
        elif dist_type == 5:
            dist = 'GP'
            samp = np.zeros(nS) 
        else:
            dist = 'Normal'
            mu = Par[1]
            st = Par[1] / 1e4
            samp = stats.norm.rvs(loc=mu, scale=st, size=nS)
            ListPar[ii, 4] = mu
            ListPar[ii, 5] = st
            Par[4] = mu
            Par[5] = st

        if dist.lower() == 'normal':
            mu, st = Par[4], Par[5]
            senA = (samp - mu) / st**2
            senB = -1/st + (samp - mu)**2 / st**3
            senC = np.full(nS, np.nan)
            senD = np.full(nS, np.nan)
            
        elif dist.lower() == 'lognormal':
            mu, st = Par[4], Par[5]
            senA = (np.log(samp) - mu) / st**2
            senB = -1/st + (np.log(samp) - mu)**2 / st**3
            senC = np.full(nS, np.nan)
            senD = np.full(nS, np.nan)
            
        elif dist.lower() == 'gamma':
            k, theta = Par[4], Par[5]
            senA = np.log(samp / theta) - digamma(k)
            senB = -k/theta + samp / theta**2
            senC = np.full(nS, np.nan)
            senD = np.full(nS, np.nan)
            
        elif dist.lower() == 'gev':
            mu, sigma = Par[4], Par[5]
            z = (samp - mu) / sigma
            senA = (-1 + np.exp(-z)) * (-1 / sigma)
            senB = (-1 + np.exp(-z)) * (-z / sigma) - 1 / sigma
            senC = np.full(nS, np.nan)
            senD = np.full(nS, np.nan)
            
        elif dist.lower() == 'uniform':
            a, b = Par[4], Par[5]
            senA = (1 / (b - a)) * (samp != 0)
            senB = -(1 / (b - a)) * (samp != 0)
            senC = np.full(nS, np.nan)
            senD = np.full(nS, np.nan)

        xS.samp[:, ii] = samp
        xS.senA[:, ii] = senA
        xS.senB[:, ii] = senB
        xS.senC[:, ii] = senC
        xS.senD[:, ii] = senD

        ParSen[ii][0] = senA
        ParSen[ii][1] = senB

    return xS, ListPar, ParSen


def cal_jpdf_hist(y, xS, Ny):
    """Estimates the joint pdf using the method of histogram."""
    Ns, Ne = y.shape
    nbins = Ny - 1
    edge_v = [nbins] * Ne
    
    # Will call the python version of histcn once imported
    epdf, epdf_dp, edges, h_bin = histcn(y, xS, *edge_v)
    
    y_v = np.array(edges).T

    yjpdf = SimpleNamespace()
    yjpdf.p_y = epdf
    yjpdf.dp_y = epdf_dp
    yjpdf.y_v = y_v.T
    
    return yjpdf



def cal_jFisher(yjpdf_data, nPar):
        print("   [DEBUG] 1. Extracting PDF data...", flush=True)
        p_y = np.nan_to_num(yjpdf_data.p_y, nan=0.0)
        dp_y = yjpdf_data.dp_y
        
        # SAFEGUARD: Ensure y_v is (Bins, Dimensions) so we don't build a 30-D universe
        y_v = np.array(yjpdf_data.y_v)
        if y_v.shape[0] < y_v.shape[1]:
            y_v = y_v.T
        Ny, Ne = y_v.shape
        print(f"   [DEBUG] 1b. Grid check: {Ny} bins, {Ne} dimensions.", flush=True)

        print("   [DEBUG] 2. Computing grid differentials safely...", flush=True)
        dy_v = np.vstack([np.zeros((1, Ne)), np.diff(y_v, axis=0)])
        
        # Bypass meshgrid entirely using memory-safe broadcasting
        dynD = np.ones_like(p_y)
        for i in range(Ne):
            shape = [1] * Ne
            shape[i] = Ny
            dy_1d = dy_v[:, i].reshape(shape)
            dynD = dynD * dy_1d

        print("   [DEBUG] 3. Flattening arrays...", flush=True)
        p_y_flat = p_y.flatten()
        dynD_flat = dynD.flatten()

        print("   [DEBUG] 4. Building sensitivity matrix...", flush=True)
        DP = np.zeros((nPar * 4, p_y_flat.size))
        for ii in range(nPar * 4):
            idx1 = ii % nPar
            idx2 = ii // nPar
            DP[ii, :] = np.nan_to_num(dp_y[idx1][idx2], nan=0.0).flatten()

        print("   [DEBUG] 5. Calculating safe weights...", flush=True)
        safe_p_y = np.where(p_y_flat > 1e-15, p_y_flat, np.inf)
        weight = dynD_flat / safe_p_y

        print("   [DEBUG] 6. Matrix Math (Mac Safe)...", flush=True)
        DP_safe = np.ascontiguousarray(DP)
        weighted_DP = np.ascontiguousarray(DP_safe * weight[np.newaxis, :])
        F = weighted_DP @ DP_safe.T

        print("   [DEBUG] 7. Fisher complete!", flush=True)
        return F



def parTran(Fraw, ListPar, parJ, isNorm=1):
    """Conducts re-parameterization and normalization using matrix transformations."""
    Fn = parJ.T @ Fraw @ parJ 

    if isNorm == 1: 
        parA = ListPar[:, 4]
        parB = ListPar[:, 5]
        b_v = np.diag(np.concatenate([parA, parB]))
        Fn = b_v @ Fn @ b_v 

    elif isNorm == 2: 
        parMean = ListPar[:, 1]
        parStd = parMean * ListPar[:, 2]
        b_v = np.diag(np.concatenate([parMean, parStd]))
        Fn = b_v @ Fn @ b_v 

    elif isNorm == 3: 
        parMean = ListPar[:, 1]
        parCoV = ListPar[:, 2]
        parStd = np.zeros_like(parMean)
        
        for i in range(len(parMean)):
            if parMean[i] == 0:
                parStd[i] = parCoV[i]
            else:
                parStd[i] = np.abs(parMean[i] * parCoV[i])
                
        b_v = np.diag(np.concatenate([parStd, parStd]))
        Fn = b_v @ Fn @ b_v 
        
    return Fn, b_v


def cal_h(xS, Opts):
    """Evaluate the blackbox function over the generated samples."""
    nS = Opts.nSampMC
    
    # Safely get the target function from the global namespace of this module
    hfunction = globals().get(Opts.funName)

    if Opts.funName == 'trial':
        y = np.sum(xS.samp, axis=1)
        ys = None 
    elif Opts.funName == 'design_dummy':
        y = xS.samp
        ys = None
    else:
        if hfunction is None:
            raise NameError(f"Function '{Opts.funName}' is not defined or imported in teds_utils.py")

        yout0 = hfunction(xS.samp[0, :], Opts)
        Ne = len(yout0['y'])
        
        y = np.zeros((nS, Ne))
        ys = [None] * nS

        for ii in tqdm(range(nS), desc="Evaluating Blackbox Function"):
            yout = hfunction(xS.samp[ii, :], Opts)
            y[ii, :] = yout['y']
            ys[ii] = yout['ys']

    h_Results = SimpleNamespace()
    h_Results.y = y
    h_Results.ys = ys
    return h_Results

# =========================================================================
# Blackbox Model & Plotting (Translated from design_B3.m)
# =========================================================================

def solve4frf_B3(xS, ifocus=None):
    """Solve for Frequency Response Function (FRF) for a 3-DOF system."""
    Ndof = 3
    
    # 0-based indexing for Python (xS[0] replaces xS(1))
    m1, m2, m3 = xS[0], xS[1], xS[2]
    k1, k2, k3 = xS[3], xS[4], xS[5]
    c1, c2, c3 = xS[6], xS[7], xS[8]
    F1 = 1.0  # xS[9] was commented out in MATLAB
    
    # Stiffness, Damping and Mass matrices
    K = np.array([
        [k1 + k2, -k2,       0],
        [-k2,      k2 + k3, -k3],
        [0,       -k3,       k3]
    ])
    
    C = np.array([
        [c1 + c2, -c2,       0],
        [-c2,      c2 + c3, -c3],
        [0,       -c3,       c3]
    ])
    
    M = np.diag([m1, m2, m3])
    
    # Eigen analysis (using eigh since K and M are symmetric positive definite)
    D_om, V = la.eigh(K, M)
    
    # Sort eigenvalues and eigenvectors
    idx = np.argsort(D_om)
    D_om = D_om[idx]
    V = V[:, idx]
    
    # Natural frequencies
    # Using np.abs to prevent warning for tiny numerical negative zeros
    omn = np.sqrt(np.abs(D_om)) 
    fn = omn / (2 * np.pi)
    
    # Mass normalization
    for ii in range(Ndof):
        V[:, ii] = V[:, ii] / np.sqrt(np.dot(V[:, ii].T, np.dot(M, V[:, ii])))
        
    # Forced response
    if ifocus is None:
        Nf = 100
        f_up = 1.2 * np.max(fn)
        f_v = np.linspace(0, f_up, Nf)
    else:
        Nf = 1
        # ifocus is 1-based from MATLAB logic, so subtract 1 for Python
        f_v = np.array([fn[ifocus - 1]]) 
        
    Fb = np.array([F1, 0, 0])
    recp_v = np.zeros((Ndof, Nf), dtype=complex)
    
    for ii in range(Nf):
        om = 2 * np.pi * f_v[ii]
        D_b = -om**2 * M + 1j * om * C + K  # 1j is Python's imaginary unit
        # Solve the linear system (equivalent to D_b \ Fb in MATLAB)
        recp_v[:, ii] = np.linalg.solve(D_b, Fb)
        
    return f_v, recp_v


def design_B3(xS, Opts):
    """Calculate the FRF and shear forces."""
    ifocus = 1
    f_v, recp_v = solve4frf_B3(xS, ifocus)
    
    # Get stiffness and damping parameters again
    k1, k2, k3 = xS[3], xS[4], xS[5]
    c1, c2, c3 = xS[6], xS[7], xS[8]
    
    # Responses at each floor
    y1 = recp_v[0, :]
    y2 = recp_v[1, :]
    y3 = recp_v[2, :]
    
    # Calculate shear forces
    s1 = k1 * y1 + c1 * 1j * (2 * np.pi * f_v) * y1
    s2 = k2 * (y2 - y1) + c2 * 1j * (2 * np.pi * f_v) * (y2 - y1)
    s3 = k3 * (y3 - y2) + c3 * 1j * (2 * np.pi * f_v) * (y3 - y2)
    
    # Combine outputs (flatten to ensure 1D array output for cal_h)
    recp_abs = np.abs(recp_v).T
    s_abs = np.column_stack([np.abs(s1), np.abs(s2), np.abs(s3)])
    
    y_out = np.hstack([recp_abs, s_abs]).flatten()
    
    return {'y': y_out, 'ys': f_v}


def display_jpdf_design_B3(y, yjpdf, D_e, V_e, nPar, Opts):
    """Plotting function for Fisher Eigenvalues and Eigenvectors."""
    varName = ['m_1', 'm_2', 'm_3', 'k_1', 'k_2', 'k_3', 'c_1', 'c_2', 'c_3', 'F1']
    
    lambda_val = np.diag(D_e)
    x_indices = np.arange(1, nPar * 2 + 1)
    
    # 1. Plot Fisher eigenvalues
    plt.figure(figsize=(8, 5))
    plt.bar(x_indices, lambda_val)
    plt.xlabel('Index of Fisher Eigenvalue', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.title('Fisher Eigenvalues')
    plt.tight_layout()
    plt.show()

    displayMode = Opts.displayMode
    labels = varName + varName  # Concatenate lists for Mean and Std Dev
    
    if displayMode == 1:
        # Plot eigenvectors (first 4)
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
        
        for ii in range(4):
            ax = axes[ii]
            bars = ax.bar(x_indices, V_e[:, ii])
            ax.set_ylim([-1, 1])
            ax.tick_params(axis='both', which='major', labelsize=14)
            
            # X-ticks for Mean and Std Dev regions
            ax.set_xticks([round(nPar / 2), nPar + round(nPar / 2)])
            ax.set_xticklabels(['Mean', 'Std Dev'])
            
            # Add text labels on top/bottom of bars
            for bar, label in zip(bars, labels):
                yval = bar.get_height()
                offset = 0.05 if yval > 0 else -0.05
                va = 'bottom' if yval > 0 else 'top'
                if yval != 0:  # Only label visible bars
                    ax.text(bar.get_x() + bar.get_width()/2, yval + offset, 
                            label, ha='center', va=va, fontsize=10)
            
            ax.set_title(f'No. {ii+1} Fisher Eigenvector [$\\lambda_{ii+1}$ = {lambda_val[ii]:.1e}]')
            
        plt.tight_layout()
        plt.show()
        
    elif displayMode == 2:
        # Plot only the 1st eigenvector
        plt.figure(figsize=(8, 5))
        bars = plt.bar(x_indices, V_e[:, 0])
        plt.ylim([-1, 1])
        plt.tick_params(axis='both', which='major', labelsize=14)
        
        plt.xticks([round(nPar / 2), nPar + round(nPar / 2)], ['Mean', 'Std Dev'])
        
        for bar, label in zip(bars, labels):
            yval = bar.get_height()
            offset = 0.05 if yval > 0 else -0.05
            va = 'bottom' if yval > 0 else 'top'
            if yval != 0:
                plt.text(bar.get_x() + bar.get_width()/2, yval + offset, 
                         label, ha='center', va=va, fontsize=10)
                
        plt.title(f'No. 1 Fisher Eigenvector [$\\lambda_1$ = {lambda_val[0]:.1e}]')
        plt.tight_layout()
        plt.show()


# =========================================================================
# N-Dimensional Histogram & Sensitivities (Translated from histcn.m)
# =========================================================================

def histcn(X, xS, *args, **kwargs):
    """
    Compute n-dimensional histogram and sensitivities.
    Equivalent to the modified MATLAB histcn function.
    """
    if X.ndim > 2:
        raise ValueError("histcn: X requires to be an (M x N) array of M points in R^N")
        
    DEFAULT_NBINS = 32
    
    # Process kwargs (AccumData, Fun) if passed (though not used in main script)
    AccumData = kwargs.get('AccumData', None)
    
    nd = X.shape[1]
    edges = list(args)
    
    if nd < len(edges):
        nd = len(edges)
    else:
        edges.extend([DEFAULT_NBINS] * (nd - len(edges)))
        
    # Build the edges for numpy's histogramdd
    np_edges = []
    h = np.zeros(nd)
    
    for d in range(nd):
        ed = edges[d]
        Xd = X[:, d]
        
        if ed is None:
            ed = DEFAULT_NBINS
            
        if np.isscalar(ed):
            # automatic linear subdivision
            ed_arr = np.linspace(np.min(Xd), np.max(Xd), int(ed) + 1)
        else:
            ed_arr = np.asarray(ed)
            
        np_edges.append(ed_arr)
        h[d] = ed_arr[1] - ed_arr[0]
        edges[d] = ed_arr

    # np.histogramdd returns the count and the edges
    count, edges_out = np.histogramdd(X, bins=np_edges)
    
    # To replicate the exact loc array and accumarray for sensitivities,
    # we need to find which bin each point falls into.
    # np.digitize returns 1-based indices. 
    # Bin i corresponds to edges_out[d][i-1] <= x < edges_out[d][i]
    loc = np.zeros(X.shape, dtype=int)
    for d in range(nd):
        # right=False means bin is [a, b)
        loc[:, d] = np.digitize(X[:, d], edges_out[d], right=False)
        
    # Points falling exactly on the rightmost edge are placed in the last bin
    # to match MATLAB's histc behavior
    for d in range(nd):
        on_right_edge = X[:, d] == edges_out[d][-1]
        loc[on_right_edge, d] = len(edges_out[d]) - 1

    # Check bounds (1 <= loc <= number of bins)
    sz = np.array([len(e) - 1 for e in edges_out])
    hasdata = np.all((loc > 0) & (loc <= sz), axis=1)
    
    # Convert to 0-based indexing for Python mapping
    valid_locs = loc[hasdata, :] - 1 
    
    # Sensitivities (The custom addition)
    Ns, nPar = xS.samp.shape
    epdf_dp = [[None for _ in range(4)] for _ in range(nPar)]
    
    # Calculate the product of the bin widths
    h_prod = np.prod(h)
    
    # If there is valid data, perform the accumarray logic
    if np.any(hasdata):
        # Convert multi-dimensional indices into a 1D index for bincount
        flat_indices = np.ravel_multi_index(valid_locs.T, sz)
        
        for kk in range(nPar):
            # For each sensitivity column, map the values into the flat bins
            # and then reshape back to the N-dimensional histogram shape
            senA_vals = xS.senA[hasdata, kk]
            sumA = np.bincount(flat_indices, weights=senA_vals, minlength=np.prod(sz))
            epdf_dp[kk][0] = sumA.reshape(sz) / Ns / h_prod
            
            senB_vals = xS.senB[hasdata, kk]
            sumB = np.bincount(flat_indices, weights=senB_vals, minlength=np.prod(sz))
            epdf_dp[kk][1] = sumB.reshape(sz) / Ns / h_prod
            
            senC_vals = xS.senC[hasdata, kk]
            sumC = np.bincount(flat_indices, weights=senC_vals, minlength=np.prod(sz))
            epdf_dp[kk][2] = sumC.reshape(sz) / Ns / h_prod
            
            senD_vals = xS.senD[hasdata, kk]
            sumD = np.bincount(flat_indices, weights=senD_vals, minlength=np.prod(sz))
            epdf_dp[kk][3] = sumD.reshape(sz) / Ns / h_prod
    else:
        # If no data fell into the bins, return empty arrays of the correct shape
        empty_hist = np.zeros(sz)
        for kk in range(nPar):
            epdf_dp[kk][0] = empty_hist
            epdf_dp[kk][1] = empty_hist
            epdf_dp[kk][2] = empty_hist
            epdf_dp[kk][3] = empty_hist

    epdf = count / Ns / h_prod
    
    # MATLAB returned multiple things. The main script expects: epdf, epdf_dp, edges, h
    # We will return the edges excluding the last point to match MATLAB's output format 
    # expected by `cal_jpdf_hist` where it does `y_v= cell2mat(edges.')`
    
    edges_mid = [e[:-1] for e in edges_out]

    return epdf, epdf_dp, edges_mid, h

# =========================================================================
# (4) KPI-based Sensitivity Analysis (Translated from calSen_KPI.m)
# =========================================================================

def calSen_KPI(y, yExLevel, isPrctile, nQoI, xS):
    """
    Calculates sensitivity for failure probability (KPI-based).
    """
    Ns, N_QoI = y.shape
    
    # Initialize arrays
    pF = np.zeros((Ns, N_QoI), dtype=bool)
    pFMean = np.zeros(N_QoI)
    pFSenC = []  # Standard Python list acts perfectly as a MATLAB cell array
    
    # Handle the repmat equivalent for threshold levels
    if nQoI <= N_QoI:
        if np.isscalar(yExLevel) or len(np.atleast_1d(yExLevel)) == 1:
            yExLevel_arr = np.full(N_QoI, float(yExLevel))
        else:
            yExLevel_arr = np.array(yExLevel)
    else:
        yExLevel_arr = np.array(yExLevel)

    _, N_UPar = xS.samp.shape
    
    for ii in range(N_QoI):
        # Add small noise to randomise the percentile every time the function is called
        PeakThresholdFactor = yExLevel_arr[ii] + (np.random.rand() - 0.5)
        
        if isPrctile == 1:
            # np.nanpercentile automatically ignores NaNs and works on a 0-100 scale
            yPeakThreshold = np.nanpercentile(y[:, ii], PeakThresholdFactor)
        else:
            # Set absolute threshold for failure
            yPeakThreshold = PeakThresholdFactor
            
        # Indicator for failure 
        # (Python broadcasting automatically handles what repmat did in MATLAB!)
        pF[:, ii] = y[:, ii] >= yPeakThreshold
        
        # Unconditional probability of failure 
        pFMean[ii] = np.nanmean(pF[:, ii])
        
        # --- SENSITIVITY CALCULATION (Vectorized for Python Speed) ---
        pFSen = np.zeros((2, N_UPar))
        
        # Convert boolean failure indicator to float for math operations
        pF_col = pF[:, ii].astype(float)
        
        # Vectorized multiplication and mean along the columns (axis=0)
        pFSen[0, :] = np.nanmean(pF_col[:, np.newaxis] * xS.senA, axis=0)
        pFSen[1, :] = np.nanmean(pF_col[:, np.newaxis] * xS.senB, axis=0)
        
        # Append to our "cell array"
        pFSenC.append(pFSen)
        
    return pF, pFMean, pFSenC


# =========================================================================
# KPI Sensitivity Visualization
# =========================================================================
import matplotlib.pyplot as plt
import numpy as np

def display_kpi_sensitivity(r, pFMean, pFInex, nPar, V_e):
    """
    Plots the sensitivity of failure probability and its projection 
    onto the Fisher Information Matrix eigenvectors.
    """
    varName = ['m_1', 'm_2', 'm_3', 'k_1', 'k_2', 'k_3', 'c_1', 'c_2', 'c_3', 'F1']
    labels = varName + varName  # Concatenate for Mean and Std Dev
    
    # Extract the specific column to plot
    r_plot = r[:, pFInex]
    x_indices = np.arange(1, nPar * 2 + 1)
    
    # Create the figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # ---------------------------------------------------------
    # Subplot 1: Sensitivity of Failure Probability
    # ---------------------------------------------------------
    bars = ax1.bar(x_indices, r_plot, color='steelblue')
    
    title_str = f"Sensitivity of Failure Probability [Pf = {pFMean[pFInex]:.4g}]"
    ax1.set_title(title_str, loc='left', fontsize=14, fontweight='bold')
    
    # Formatting X-axis
    ax1.set_xticks([round(nPar / 2), nPar + round(nPar / 2)])
    ax1.set_xticklabels(['Mean', 'Std Dev'])
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.set_ylabel('r [-]', fontsize=14)
    
    # Dynamic Y-limits based on data
    y_min, y_max = np.min(r_plot), np.max(r_plot)
    ax1.set_ylim([y_min - 2, y_max + 5])
    
    # Add text labels on top/bottom of bars
    for bar, label in zip(bars, labels):
        yval = bar.get_height()
        # Only display text if the bar has a positive value (matching MATLAB logic)
        if yval > 0:
            ax1.text(bar.get_x() + bar.get_width()/2, yval + 0.1, 
                     label, ha='center', va='bottom', fontsize=11)

    # ---------------------------------------------------------
    # Subplot 2: Projection onto FIM Eigenvectors
    # ---------------------------------------------------------
    # Calculate projection 's'
    rc = r_plot
    s = (rc.T @ V_e) / np.linalg.norm(rc)
    
    ax2.bar(x_indices, np.abs(s), color='darkorange')
    
    ax2.set_title('Projection onto FIM EigVector', loc='left', fontsize=14, fontweight='bold')
    
    # Formatting
    ax2.set_ylim([0, 1])
    ax2.set_xlim([0.5, nPar * 2 + 0.5])
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    plt.tight_layout()
    plt.show()