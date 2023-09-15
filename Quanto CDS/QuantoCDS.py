import math

class CreditDefaultSwap:
    def __init__(self, is_buy, notional, cds_spread, frequency, time_to_maturity, recovery_rate, hazard_rate, zero_rate, fx_jump=0.0, correlation=0.0, credit_vol=0.0, fx_vol=0.0):
        self.is_buy = is_buy
        self.notional = notional
        self.cds_spread = cds_spread
        self.frequency = frequency
        self.time_to_maturity = time_to_maturity
        self.recovery_rate = recovery_rate
        self.zero_rate = zero_rate
        self.nCoupons = int(time_to_maturity/frequency)
        self.LGD = notional * (1.0 - recovery_rate)
        self.fx_jump = fx_jump
        self.correlation = correlation
        self.credit_vol = credit_vol
        self.fx_vol = fx_vol
        # Apply Quanto Adjustment,if required
        # See Minqiang et al (2018)
        self.hazard_rate = (1.0 + fx_jump) * (1 + 0.5 * correlation * credit_vol * fx_vol * time_to_maturity) * hazard_rate
    
    # Marginal Default i.e. prob default in a given coupon period
    def probDefault(self, t):
        if (t<=self.frequency):
            return 1.0 - math.exp(-self.hazard_rate * t)
        else:
            return math.exp(-self.hazard_rate * (t-self.frequency)) - math.exp(-self.hazard_rate * t)
    
    def probSurvive(self, t):
        return math.exp(-self.hazard_rate * t)
    
    def discFactor(self, t):
        return math.exp(-self.zero_rate * t)
    
    def premiumCoupon(self):
        n = self.notional
        s = self.cds_spread
        tau = self.frequency
        cpn = n * s * tau    
        return cpn
    
    # Premium Paid on Survival
    def premiumPV(self):
        premium_pv = 0
        cpn = self.premiumCoupon()
        
        for i in range(1, self.nCoupons + 1):
            t = i * self.frequency
            pSurvive = self.probSurvive(t)
            df = self.discFactor(t)
            coupon_pv = cpn * pSurvive * df
            #Diagnostics
            #print("Premium PV({0}) : {1:0.2f}".format(t, coupon_pv))
            premium_pv += coupon_pv

        # buy protection i.e. pay premium
        if self.is_buy:
            return -1.0 * premium_pv
        else:
            return premium_pv
    
    # Accrued Interest Paid on Default
    def accruedInterestPV(self):
        accrued_pv = 0
        # Model assumes default mid-period i.e. pay only half coupon
        cpn = self.premiumCoupon() * 0.5 
        
        for i in range(1, self.nCoupons + 1):
            t = i * self.frequency
            pDefault = self.probDefault(t)
            df = self.discFactor(t)
            coupon_pv = cpn * pDefault * df
            #Diagnostics
            #print("Accrued PV({0}) : {1:0.2f}".format(t, coupon_pv))
            accrued_pv += coupon_pv

        # buy protection i.e. pay premium accrued on default
        if self.is_buy:
            return -1.0 * accrued_pv
        else:
            return accrued_pv
    
    # Protection Paid on Default
    def protectionPV(self):
        protection_pv = 0
       
        for i in range(1, self.nCoupons + 1):
            t = i * self.frequency
            pDefault = self.probDefault(t)
            df = self.discFactor(t)
            lgd_pv = self.LGD * pDefault * df
            #Diagnostics
            #print("Protection PV({0}) : {1:0.2f}".format(t, lgd_pv))
            protection_pv += lgd_pv
        
        # buy protection i.e. receive protection
        if self.is_buy:
            return protection_pv
        else:
            return -1.0 * protection_pv
    
    def cdsPV(self):
        return self.premiumPV() + self.accruedInterestPV() + self.protectionPV()

# Example usage
is_buy = True           # Buy Protection
notional = 10000000     # $10,000,000
cds_spread = 0.01       # 100 basis points
frequency = 0.25        # Quarterly Coupons
time_to_maturity = 5.0  # 5 years
recovery_rate = 0.25    # 25% recovery rate
hazard_rate = 0.022947  # 2.2947%
zero_rate = 0.14        # 14%

# Vanilla CDS
cds = CreditDefaultSwap(is_buy, notional, cds_spread, frequency, time_to_maturity, recovery_rate, hazard_rate, zero_rate, 0.0, 0.0, 0.0, 0.0)

print("VANILLA CDS")
print("Premium PV: {:0.2f}".format(cds.premiumPV()))
print("Accrued PV: {:0.2f}".format(cds.accruedInterestPV()))
print("Protection PV: {:0.2f}".format(cds.protectionPV()))
print("CDS PV: {:0.2f}".format(cds.cdsPV()))
print("")

# Quanto CDS
fx_jump = -0.40     # -40% FX crash risk on default
correlation = -0.25 # -0.25 correlation between underlying credit and FX
credit_vol = 0.25   # 25% credit volatility in Hazard Rate
fx_vol = 0.1        # 10% FX volatility

qcds = CreditDefaultSwap(is_buy, notional, cds_spread, frequency, time_to_maturity, recovery_rate, hazard_rate, zero_rate, fx_jump, correlation, credit_vol, fx_vol)

print("QUANTO CDS")
print("Premium PV: {:0.2f}".format(qcds.premiumPV()))
print("Accrued PV: {:0.2f}".format(qcds.accruedInterestPV()))
print("Protection PV: {:0.2f}".format(qcds.protectionPV()))
print("CDS PV: {:0.2f}".format(qcds.cdsPV()))
