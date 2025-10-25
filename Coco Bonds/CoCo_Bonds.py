"""
AT1 CoCo Bond Monte Carlo Pricing Demo
--------------------------------------

Educational example showing how to price an Additional Tier 1 (AT1)
Contingent Convertible (CoCo) bond using correlated stochastic processes
for interest rates, credit risk, and equity (CET1 proxy).

These instruments are hybrid: part bond, part equity option.
They can be written down or converted to equity when the issuer's
capital ratio (CET1) falls below a regulatory trigger.

We model three correlated risk factors:
  1. Short rate (Hull–White model) — drives discounting / IR risk
  2. Credit hazard rate (CIR model) — drives survival probability / credit risk
  3. Equity ratio (GBM) — proxy for CET1 capital ratio / trigger dynamics

Conversion logic can follow:
  - WRITE-DOWN of principal (CONVERSION_TYPE = "write_down")
  - CONVERSION TO EQUITY (CONVERSION_TYPE = "convert_to_equity")

The goal: to simulate expected discounted cashflows under these correlated risks.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# DEAL PARAMETERS
# ============================================================

NOTIONAL = 100.0             # Bond face value
COUPON_RATE = 0.05           # 5% annual coupon
COUPON_FREQ = 2              # Semiannual coupons (2 per year)
T_MATURITY = 10.0            # Bond maturity (years)
CALL_DATE = 5.0              # First call date
CALL_PRICE = 100.0           # Redemption price if called

# Conversion settings
CONVERSION_TYPE = "write_down"  # "write_down" or "convert_to_equity"
WRITE_DOWN_RATIO = 1.0          # 1.0 = full write-down, 0 = no write-down

# Conversion to Equity (only relevant if CONVERSION_TYPE = "convert_to_equity")
INITIAL_SHARE_PRICE = 100.0     # Current share price
CONVERSION_PRICE = 50.0         # Price per share when converted
CONVERSION_RATIO = NOTIONAL / CONVERSION_PRICE  # Shares received per bond

# CET1 trigger (regulatory capital threshold)
CET1_TRIGGER_PERCENT = 0.07     # Trigger level (e.g., 7%)
CET1_CURRENT_PERCENT = 0.12     # Current CET1 (e.g., 12%)
EQ_BARRIER = CET1_TRIGGER_PERCENT / CET1_CURRENT_PERCENT  # normalized (~0.583)

# ============================================================
# MODEL PARAMETERS
# ============================================================

# --- Hull–White Interest Rate Model (Short rate) ---
# dr = a*(r_mean - r)*dt + sigma*dW_r
# Captures mean-reverting dynamics of the short-term rate.
HW_R0 = 0.02       # Starting rate (2%)
HW_A = 0.03        # Mean reversion speed
HW_SIGMA = 0.01    # Volatility of the short rate

# --- CIR Model (Hazard rate / credit intensity) ---
# dh = k*(theta - h)*dt + sigma*sqrt(h)*dW_h
# Ensures hazard rate remains non-negative.
CIR_L0 = 0.015     # Initial hazard rate (1.5%)
CIR_K = 0.3        # Mean reversion speed
CIR_THETA = 0.02   # Long-term mean hazard rate
CIR_SIGMA = 0.005  # Volatility of hazard rate

# --- GBM for Equity (CET1 proxy) ---
# dS/S = mu*dt + sigma*dW_e
# Tracks normalized CET1 ratio as an equity-like process.
EQ_SPOT = 1.0
EQ_MU = 0.00       # Drift (expected growth of CET1)
EQ_SIGMA = 0.15    # Volatility of CET1 proxy

# --- Correlations between drivers ---
RHO_R_H = 0.5      # Rate ↔ Hazard (credit/IR correlation)
RHO_R_E = 0.3      # Rate ↔ Equity (macro link)
RHO_H_E = 0.3      # Hazard ↔ Equity (credit-equity link)

# --- Simulation controls ---
DT = 1 / 12        # Monthly time step
T_STEPS = int(T_MATURITY / DT)
N_PATHS = 10000
ANTITHETIC = False # Use antithetic variance reduction
SEED = 42

# Derived parameters
np.random.seed(SEED)
CHOLESKY = np.linalg.cholesky(np.array([
    [1.0, RHO_R_H, RHO_R_E],
    [RHO_R_H, 1.0, RHO_H_E],
    [RHO_R_E, RHO_H_E, 1.0]
]))
EQ_DRIFT = EQ_MU - 0.5 * EQ_SIGMA ** 2

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================


def simulate_paths(num_paths=N_PATHS, antithetic=ANTITHETIC):
    """
    Simulate correlated Hull–White, CIR, and GBM paths.
    If antithetic is True, generate half the paths and mirror them after.
    Returns full array of N_PATHS x (T_STEPS+1).
    """
    base_paths = num_paths // 2 if antithetic else num_paths

    rates = np.zeros((base_paths, T_STEPS + 1))
    hazards = np.zeros((base_paths, T_STEPS + 1))
    equities = np.zeros((base_paths, T_STEPS + 1))

    rates[:, 0] = HW_R0
    hazards[:, 0] = CIR_L0
    equities[:, 0] = EQ_SPOT
    sqrt_dt = np.sqrt(DT)

    for t in range(T_STEPS):
        z = np.random.normal(size=(base_paths, 3))
        dW = z @ CHOLESKY.T

        # --- Hull–White (mean-reverting short rate) ---
        dr = HW_A * (HW_R0 - rates[:, t]) * DT + HW_SIGMA * dW[:, 0] * sqrt_dt
        rates[:, t + 1] = rates[:, t] + dr

        # --- CIR (hazard rate, positive process) ---
        sqrt_h = np.sqrt(np.maximum(hazards[:, t], 0))
        dh = CIR_K * (CIR_THETA - hazards[:, t]) * DT + CIR_SIGMA * sqrt_h * dW[:, 1] * sqrt_dt
        hazards[:, t + 1] = np.maximum(hazards[:, t] + dh, 0)

        # --- GBM (Equity for CET1 proxy) ---
        dE = (EQ_MU - 0.5 * EQ_SIGMA**2) * DT + EQ_SIGMA * sqrt_dt * dW[:, 2]
        equities[:, t + 1] = equities[:, t] * np.exp(dE)

    # Antithetic mirroring
    if antithetic:
        rates = np.vstack([rates, rates[:, ::-1]])
        hazards = np.vstack([hazards, hazards[:, ::-1]])
        equities = np.vstack([equities, 2*EQ_SPOT - equities[:, ::-1]])  # mirror around initial equity

    return rates, hazards, equities

# ============================================================
# PRICING FUNCTION
# ============================================================

def price_path(rates, hazard, equity, dt, coupon, notional, pay_times,
               call_date, call_price, eq_barrier, conv_ratio, conv_type, write_down_ratio):
    """Compute discounted PV for a single path, accounting for call and conversion events."""

    df_path = np.exp(-np.cumsum((rates[:-1] + rates[1:]) / 2 * dt))
    surv_path = np.exp(-np.cumsum((hazard[:-1] + hazard[1:]) / 2 * dt))
    n_pay = len(pay_times)

    pv = 0.0
    called = False
    converted = False

    for i, t in enumerate(pay_times):
        pv += coupon * surv_path[i] * df_path[i]

        if abs(t - call_date) < 1e-9:
            continuation_value = sum(coupon * surv_path[j] * df_path[j] for j in range(i + 1, n_pay))
            continuation_value += notional * surv_path[-1] * df_path[-1]
            if continuation_value > call_price:
                pv += call_price * surv_path[i] * df_path[i]
                called = True
                return pv, called, converted

        if equity[i] < eq_barrier:
            converted = True
            if conv_type == "write_down":
                pv += notional * (1 - write_down_ratio) * surv_path[i] * df_path[i]
            elif conv_type == "convert_to_equity":
                share_price = INITIAL_SHARE_PRICE * equity[i]
                pv += conv_ratio * share_price * surv_path[i] * df_path[i]
            return pv, called, converted

    pv += notional * surv_path[-1] * df_path[-1]
    return pv, called, converted

# ============================================================
# DIAGNOSTICS AND VISUALIZATION
# ============================================================

def run_diagnostics(df, convergence_data):
    """Print event breakdown and convergence diagnostics."""
    summary = df.groupby('Event')['PV'].agg(['count', 'mean', 'std'])
    summary['Probability'] = summary['count'] / len(df)
    print("\nEvent Summary:")
    print(summary[['count', 'Probability', 'mean', 'std']].round(3))

    if convergence_data:
        conv_df = pd.DataFrame(convergence_data, columns=['Paths', 'MeanPV'])
        print("\nConvergence Table:")
        print(conv_df.to_string(index=False))

        plt.figure(figsize=(10, 5))
        plt.plot(conv_df['Paths'], conv_df['MeanPV'], marker='o')
        plt.title("Monte Carlo Convergence (Mean PV vs Paths)")
        plt.xlabel("Paths Simulated")
        plt.ylabel("Mean PV")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.show()

def plot_curves(rates, hazards, equities, num_paths_plot=30):
    """
    Plot Monte Carlo paths for interest rates, hazard rates, and CET1 proxy
    as three separate figures for clarity.
    """
    num_paths_plot = min(num_paths_plot, rates.shape[0])
    num_steps = rates.shape[1]
    time_grid = np.linspace(0, T_MATURITY, num_steps)

    # --- Hull–White Short Rate Paths ---
    plt.figure(figsize=(10, 4))
    for i in range(num_paths_plot):
        plt.plot(time_grid, rates[i, :], alpha=0.3, linewidth=0.8)
    plt.title("Interest Rate Paths (Hull–White Model)")
    plt.xlabel("Time (Years)")
    plt.ylabel("Short Rate")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

    # --- CIR Hazard Rate Paths ---
    plt.figure(figsize=(10, 4))
    for i in range(num_paths_plot):
        plt.plot(time_grid, hazards[i, :], alpha=0.3, linewidth=0.8, color='orange')
    plt.title("Hazard Rate Paths (CIR Model)")
    plt.xlabel("Time (Years)")
    plt.ylabel("Hazard Rate")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

    # --- CET1 Proxy (Equity GBM) Paths ---
    plt.figure(figsize=(10, 4))
    for i in range(num_paths_plot):
        plt.plot(time_grid, equities[i, :], alpha=0.3, linewidth=0.8, color='green')
    plt.axhline(EQ_BARRIER, color='red', linestyle='--', label='Conversion Barrier')
    plt.title("Equity CET1 Proxy Paths (GBM)")
    plt.xlabel("Time (Years)")
    plt.ylabel("Equity Level (Normalized)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

# ============================================================
# SIMULATION DRIVER
# ============================================================

def run_simulation(num_paths=N_PATHS, use_antithetic=ANTITHETIC, show_diagnostics=True, conv_step=1000):
    """
    Full Monte Carlo run with diagnostics.
    """
    dt = DT
    coupon = COUPON_RATE / COUPON_FREQ * NOTIONAL
    pay_times = np.arange(1 / COUPON_FREQ, T_MATURITY + 1e-9, 1 / COUPON_FREQ)

    rates, hazards, equities = simulate_paths(num_paths=num_paths, antithetic=use_antithetic)

    pv_list, called_list, converted_list, conv_list = [], [], [], []

    for i in range(num_paths):
        pv, called, converted = price_path(
            rates[i], hazards[i], equities[i], dt,
            coupon, NOTIONAL, pay_times,
            CALL_DATE, CALL_PRICE, EQ_BARRIER,
            CONVERSION_RATIO, CONVERSION_TYPE, WRITE_DOWN_RATIO
        )
        pv_list.append(pv)
        called_list.append(called)
        converted_list.append(converted)

        if (i + 1) % conv_step == 0:
            conv_list.append((i + 1, np.mean(pv_list)))

    df = pd.DataFrame({
        'PV': pv_list,
        'Called': called_list,
        'Converted': converted_list
    })
    df['Event'] = np.where(df['Converted'], 'Conversion', np.where(df['Called'], 'Called', 'Alive'))

    price = np.mean(pv_list)
    print(f"\n=== AT1 Convertible Bond Price: {price:.2f} ===")

    if show_diagnostics:
        plot_curves(rates, hazards, equities)
        run_diagnostics(df, conv_list)

    return df, rates, hazards, equities

# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == "__main__":
    run_simulation(num_paths=15000, use_antithetic=False, show_diagnostics=True)
