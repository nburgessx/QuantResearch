/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */

#include "funcdef.h"


#define EXCEL_FUNCTION(name, description) \
    TFuncDef name##_FN_DEF = \
    { \
    #name, \
    sizeof(##name##_PARAMS) / sizeof(TParamDef), \
    (TParamDef*) ##name##_PARAMS, \
    description \
    }


/******************************
 * Generic addin functions.
 ******************************/

TFuncDef Version_FN_DEF          = { "Version"         , 0, NULL, "Get library version. " };
TFuncDef ErrorLogStatus_FN_DEF   = { "ErrorLogStatus"  , 0, NULL, "Get error logging status. " };
TFuncDef ErrorLogContents_FN_DEF = { "ErrorLogContents", 0, NULL, "Get last 20 messages within error log. " };
TFuncDef ErrorLogFilename_FN_DEF = { "ErrorLogFilename", 0, NULL, "Get error log filename. " };

/***** SetErrorLogFilename *****/
static TParamDef SetErrorLogFilename_PARAMS[] =
{
    {"Filename", "Name of alternative error log file"},
    {"Append"  , "Should the file be appended to if it already exists "}
};

EXCEL_FUNCTION(SetErrorLogFilename, "Set error log filename. ");


/***** SetErrorLogStatus *****/
static TParamDef SetErrorLogStatus_PARAMS[] =
{
    {"IsOn", "Desired state "}
};

EXCEL_FUNCTION(SetErrorLogStatus, "Set error logging status. ");


/******************************
 * Domain functions.
 ******************************/

/***** LoadHolidays *****/
static TParamDef LoadHolidays_PARAMS[] =
{
    {"Name"    , "Name of holiday calendar"},
    {"Filename", "Name of file containing holiday data "}
};

EXCEL_FUNCTION(LoadHolidays, "Load holiday calendar from file. ");


/***** IRZeroCurveBuild *****/
static TParamDef IRZeroCurveBuild_PARAMS[] =
{
    {"Value Date", "Date for which the PV is calculated"},
    {"Types"     , NULL},
    {"End Dates" , NULL},
    {"Rates"     , NULL},
    {"MMDCC"     , NULL},
    {"Fixed IVL" , "Interval between fixed coupon payments"},
    {"Float IVL" , "Interval between floating coupon payments"},
    {"Fixed DCC" , "Day count convention for fixed coupon payments"},
    {"Float DCC" , "Day count convention for floating coupon payments"},
    {"Swap BDC"  , "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"  , "Calendar used when adjusting coupon dates"},
    {"name"      , "Name of curve object "}
};

EXCEL_FUNCTION(IRZeroCurveBuild, "Bootstrap IR zero curve from cash and swap rates. ");


/***** IRZeroCurveMake *****/
static TParamDef IRZeroCurveMake_PARAMS[] =
{
    {"Base Date" , "Value date for zero curve. This is the date all the rates start at"},
    {"Dates"     , NULL},
    {"Rates"     , NULL},
    {"Basis"     , NULL},
    {"DCC"       , "Defaults to Act/365F"},
    {"name"      , "Name of curve object "}
};

EXCEL_FUNCTION(IRZeroCurveMake, "Recreate IR zero curve from dates and rates. ");


/***** CleanSpreadCurveBuild *****/
static TParamDef CleanSpreadCurveBuild_PARAMS[] =
{
    {"Today"             , "Risk starts at the end of today"},
    {"Start Date"        , "Date when CDS started"},
    {"Stepin Date"       , "Date when new protection begins"},
    {"CashSettle Date"   , "Date when payment made"},

    {"End Dates"         , "End date for each benchmark instrument"},
    {"Coupon Rates"      , "Coupon rate for each benchmark instrument"},
    {"Include Flags"     , "Flags to include/exclude particular benchmarks (NULL = include all)"},
    {"Pay Acc On Default", "Should accrued interest be paid on default?"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates"},
    {"Discount Curve"    , "Interest rate discount curve"},
    {"Recovery Rate"     , "Assumed recovery rate in case of default"},
    {"name"              , "Name of curve object "}
};

EXCEL_FUNCTION(CleanSpreadCurveBuild, NULL);


/***** DiscountFactor *****/
static TParamDef DiscountFactor_PARAMS[] =
{
    {"Curve"             , "Discount/clean spread curve"},
    {"Date"              , "Interpolation date "}
};

EXCEL_FUNCTION(DiscountFactor, "Interpolates a discount factor (or survival probablility) from an interest rate discount (or clean spread) curve. ");


/***** DatesAndRates *****/
static TParamDef DatesAndRates_PARAMS[] =
{
    {"Curve"             , "Discount/clean spread curve"}
};

EXCEL_FUNCTION(DatesAndRates, "Get critical dates and rates (or survival probablility) from an interest rate discount (or clean spread) curve. ");


/***** ParSpreadFlat *****/
static TParamDef ParSpreadFlat_PARAMS[] =
{
    {"Today"             , "Risk starts at the end of today"},
    {"CashSettle Date"   , "Date for which the PV is calculated"},
    {"Benchmark St Date" , "Benchmark CDS start date for internal bootstrapping of credit curve"},
    {"StepIn Date"       , "Date when step-in begins"},
    {"Start Date"        , "Date when protection begins"},
    {"End Date"          , "Date when protection ends (end of day)"},
    {"Coupon Rate"       , "Fixed coupon rate (a.k.a. deal spread) for the fee leg"},
    {"Pay Acc On Default", "Should accrued interest be paid on default?"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates."},
    {"Discount Curve"    , "Interest rate discount curve"},
    {"Upfront Charge"    , "present value of CDS"},
    {"Recovery Rate"     , "Assumed recovery rate in case of default"},
    {"Clean Price"       , "Is the price expressed as a clean price? "}
};

EXCEL_FUNCTION(ParSpreadFlat, "Return a single par spread which prices the CDS");


/***** UpfrontFlat *****/
static TParamDef UpfrontFlat_PARAMS[] =
{
    {"Today"             , "Risk starts at the end of today"},
    {"CashSettle Date"   , "Date for which the PV is calculated"},
    {"Benchmark St Date" , "Benchmark CDS start date for internal bootstrapping of credit curve"},
    {"StepIn Date"       , "Date when step-in begins"},
    {"Start Date"        , "Date when protection begins"},
    {"End Date"          , "Date when protection ends (end of day)"},
    {"Coupon Rate"       , "Fixed coupon rate (a.k.a deal spread) for the fee leg"},
    {"Pay Acc On Default", "Should accrued interest be paid on default?"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates"},
    {"Discount Curve"    , "Interest rate discount curve"},
    {"Market Par Spread" , "Par spread used bootstrapping credit clean curve"},
    {"Recovery Rate"     , "Assumed recovery rate in case of default"},
    {"Clean Price"       , "Is the price expressed as a clean price? "}
};

EXCEL_FUNCTION(UpfrontFlat, "Return the upfront charge of the CDS assuming a flat credit curve");


/***** CdsPrice *****/
static TParamDef CdsPrice_PARAMS[] =
{
    {"Today"             , "Risk starts at the end of today"},
    {"CashSettle Date"   , "Date for which the PV is calculated"},
    {"StepIn Date"       , "Date when step-in begins"},
    {"Start Date"        , "Date when protection begins"},
    {"End Date"          , "Date when protection ends (end of day)"},
    {"Coupon Rate"       , "Fixed coupon rate (a.k.a. spread) for the fee leg"},
    {"Pay Acc On Default", "Should accrued interest be paid on default?"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates"},
    {"Discount Curve"    , "Interest rate discount curve"},
    {"Spread Curve"      , "Credit clean spread curve"},
    {"Recovery Rate"     , "Assumed recovery rate in case of default"},
    {"Clean Price"       , "Is the price expressed as a clean price? "}
};

EXCEL_FUNCTION(CdsPrice, "Computes the price for a vanilla CDS. ");


/***** ParSpreads *****/
static TParamDef ParSpreads_PARAMS[] =
{
    {"Today"             , "Risk starts at the end of today"},
    {"StepIn Date"       , "Date when step in happens"},
    {"Start Date"        , "Date when CDS became/becomes effective"},
    {"End Dates"         , "End date for each benchmark instrument"},
    {"Pay Acc On Default", "Should accrued interest be paid on default?"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates"},
    {"Discount Curve"    , "Interest rate discount curve"},
    {"Spread Curve"      , "Credit clean spread curve"},
    {"Recovery Rate"     , "Assumed recovery rate in case of default"}
};

EXCEL_FUNCTION(ParSpreads, NULL);


/***** FeeLegFlows *****/
static TParamDef FeeLegFlows_PARAMS[] =
{
    {"Start Date"        , "Date when protection begins"},
    {"End Date"          , "End date for fee leg"},
    {"Coupon Rate"       , "Fixed coupon rate (a.k.a. spread) for the fee leg"},
    {"Notional"          , "Notional principal for fee leg"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates"}
};

EXCEL_FUNCTION(FeeLegFlows, "Return par spreads for a list of CDS maturity dates");


/***** DefaultedCDS *****/
static TParamDef DefaultAccrual_PARAMS[] =
{
    {"Trade Date"        , "Trade Date of CDS"},
    {"EDD Date"          , "Event Determination Date"},
    {"Effective Date"    , "The Effective Date of the CDS"},
    {"Maturity Date"     , "The Maturity Date of the CDS"},
    {"Coupon Rate"       , "Fixed coupon rate (a.k.a. spread) for the fee leg"},
    {"Notional"          , "Notional principal for fee leg"},
    {"Coupon Interval"   , "Interval between coupon payments"},
    {"Stub Type"         , "If the startDate and endDate are not on cycle determines location of coupon dates"},
    {"Payment DCC"       , "Day count convention for coupon payment"},
    {"Bad Day Convention", "Bad day convention for adjusting coupon payment dates"},
    {"Holidays"          , "Calendar used when adjusting coupon dates"}
};

EXCEL_FUNCTION(DefaultAccrual, "Computes accrual of defaulted CDS. ");


/* Array of pointers to all function definitions */
TFuncDef *gtoFuncDefList[] =
{
    &Version_FN_DEF,
    &ErrorLogContents_FN_DEF,
    &ErrorLogStatus_FN_DEF,
    &ErrorLogFilename_FN_DEF,
    &SetErrorLogFilename_FN_DEF,
    &SetErrorLogStatus_FN_DEF,
    &LoadHolidays_FN_DEF,
    &IRZeroCurveBuild_FN_DEF,
    &IRZeroCurveMake_FN_DEF,
    &CleanSpreadCurveBuild_FN_DEF,
    &DiscountFactor_FN_DEF,
    &DatesAndRates_FN_DEF,
    &ParSpreadFlat_FN_DEF,
    &UpfrontFlat_FN_DEF,
    &CdsPrice_FN_DEF,
    &ParSpreads_FN_DEF,
    &FeeLegFlows_FN_DEF,
    &DefaultAccrual_FN_DEF
};

size_t gtoFuncDefCount = 18;
