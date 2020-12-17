(window.webpackJsonp=window.webpackJsonp||[]).push([[13],{381:function(e,v,t){"use strict";t.r(v);var _=t(42),o=Object(_.a)({},(function(){var e=this,v=e.$createElement,t=e._self._c||v;return t("ContentSlotsDistributor",{attrs:{"slot-key":e.$parent.slotKey}},[t("h1",{attrs:{id:"analysis"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#analysis"}},[e._v("#")]),e._v(" Analysis")]),e._v(" "),t("h2",{attrs:{id:"analyzeglucoseprofile"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#analyzeglucoseprofile"}},[e._v("#")]),e._v(" analyzeGlucoseProfile")]),e._v(" "),t("div",{staticClass:"language-MATLAB extra-class"},[t("pre",{pre:!0,attrs:{class:"language-matlab"}},[t("code",[t("span",{pre:!0,attrs:{class:"token keyword"}},[e._v("function")]),e._v(" results "),t("span",{pre:!0,attrs:{class:"token operator"}},[e._v("=")]),e._v(" "),t("span",{pre:!0,attrs:{class:"token function"}},[e._v("analyzeGlucoseProfile")]),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("(")]),e._v("data"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(")")]),e._v("\n")])])]),t("p",[e._v("Function that computes the glycemic outcomes of a glucose profile.")]),e._v(" "),t("h3",{attrs:{id:"input"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#input"}},[e._v("#")]),e._v(" Input")]),e._v(" "),t("ul",[t("li",[t("strong",[e._v("data: timetable (required)")]),e._v(" "),t("br"),e._v("\nA timetable with columns "),t("code",[e._v("Time")]),e._v(" and "),t("code",[e._v("glucose")]),e._v(" containing the glucose data to analyze (mg/dl).")])]),e._v(" "),t("h3",{attrs:{id:"output"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#output"}},[e._v("#")]),e._v(" Output")]),e._v(" "),t("ul",[t("li",[t("strong",[e._v("results: structure")]),e._v(" "),t("br"),e._v("\nA structure with fields containing the computed metrics and stats in the glucose profile, i.e.:\n"),t("ul",[t("li",[t("code",[e._v("variabilityMetrics")]),e._v(": a structure with fields containing the values of the computed variability metrics (i.e., {"),t("code",[e._v("aucGlucose")]),e._v(", "),t("code",[e._v("CVGA")]),e._v(", "),t("code",[e._v("cvGlucose")]),e._v(", "),t("code",[e._v("efIndex")]),e._v(", "),t("code",[e._v("gmi")]),e._v(", "),t("code",[e._v("iqrGlucose")]),e._v(", "),t("code",[e._v("jIndex")]),e._v(", "),t("code",[e._v("mageIndex")]),e._v(", "),t("code",[e._v("magePlusIndex")]),e._v(", "),t("code",[e._v("mageMinusIndex")]),e._v(", "),t("code",[e._v("meanGlucose")]),e._v(", "),t("code",[e._v("medianGlucose")]),e._v(", "),t("code",[e._v("rangeGlucose")]),e._v(", "),t("code",[e._v("sddmIndex")]),e._v(", "),t("code",[e._v("sdwIndex")]),e._v(", "),t("code",[e._v("stdGlucose")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("riskMetrics")]),e._v(": a structure with fields containing the values of the computed risk metrics (i.e., {"),t("code",[e._v("adrr")]),e._v(", "),t("code",[e._v("bgri")]),e._v(", "),t("code",[e._v("hbgi")]),e._v(", "),t("code",[e._v("lbgi")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("dataQualityMetrics")]),e._v(": a structure with fields containing the values of the computed data quality metrics (i.e., {"),t("code",[e._v("missingGlucosePercentage")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("timeMetrics")]),e._v(": a structure with fields containing the values of the computed time related metrics (i.e., {"),t("code",[e._v("timeInHyperglycemia")]),e._v(", "),t("code",[e._v("timeInSevereHyperglycemia")]),e._v(", "),t("code",[e._v("timeInHypoglycemia")]),e._v(", "),t("code",[e._v("timeInSevereHypoglycemia")]),e._v(", "),t("code",[e._v("timeInTarget")]),e._v(", "),t("code",[e._v("timeInTightTarget")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("glycemicTransformationMetrics")]),e._v(": a structure with fields containing the values of the computed glycemic transformed metrics (i.e., {"),t("code",[e._v("gradeScore")]),e._v(", "),t("code",[e._v("gradeEuScore")]),e._v(", "),t("code",[e._v("gradeHyperScore")]),e._v(", "),t("code",[e._v("gradeHypoScore")]),e._v(", "),t("code",[e._v("hypoIndex")]),e._v(", "),t("code",[e._v("hyperIndex")]),e._v(", "),t("code",[e._v("igc")]),e._v(", "),t("code",[e._v("mrIndex")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("eventMetrics")]),e._v(": a structure with fields containing the values of the computed event related metrics (i.e., {"),t("code",[e._v("gradeScore")]),e._v(", "),t("code",[e._v("gradeEuScore")]),e._v(", "),t("code",[e._v("gradeHyperScore")]),e._v(", "),t("code",[e._v("gradeHypoScore")]),e._v(", "),t("code",[e._v("hypoIndex")]),e._v(", "),t("code",[e._v("hyperIndex")]),e._v(", "),t("code",[e._v("igc")]),e._v(", "),t("code",[e._v("mrIndex")]),e._v("}) of the metrics for each glucose profile.")])])])]),e._v(" "),t("h3",{attrs:{id:"preconditions"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#preconditions"}},[e._v("#")]),e._v(" Preconditions")]),e._v(" "),t("ul",[t("li",[t("code",[e._v("data")]),e._v(" must be a timetable having an homogeneous time grid;")]),e._v(" "),t("li",[t("code",[e._v("data")]),e._v(" must contain a column named "),t("code",[e._v("Time")]),e._v(" and another named "),t("code",[e._v("glucose")]),e._v(".")])]),e._v(" "),t("h3",{attrs:{id:"reference"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#reference"}},[e._v("#")]),e._v(" Reference")]),e._v(" "),t("ul",[t("li",[e._v("None")])]),e._v(" "),t("div",{staticClass:"custom-block warning"},[t("p",{staticClass:"custom-block-title"},[e._v("WARNING")]),e._v(" "),t("p",[e._v("Currently "),t("code",[e._v("analyzeGlucoseProfile")]),e._v(" is not CI tested.")])]),e._v(" "),t("h2",{attrs:{id:"analyzeonearm"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#analyzeonearm"}},[e._v("#")]),e._v(" analyzeOneArm")]),e._v(" "),t("div",{staticClass:"language-MATLAB extra-class"},[t("pre",{pre:!0,attrs:{class:"language-matlab"}},[t("code",[t("span",{pre:!0,attrs:{class:"token keyword"}},[e._v("function")]),e._v(" results "),t("span",{pre:!0,attrs:{class:"token operator"}},[e._v("=")]),e._v(" "),t("span",{pre:!0,attrs:{class:"token function"}},[e._v("analyzeOneArm")]),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("(")]),e._v("arm"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(")")]),e._v("\n")])])]),t("p",[e._v("Function that computes the glycemic outcomes of one arm.")]),e._v(" "),t("h3",{attrs:{id:"input-2"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#input-2"}},[e._v("#")]),e._v(" Input")]),e._v(" "),t("ul",[t("li",[t("strong",[e._v("arm: cell array of timetable (required)")]),e._v(" "),t("br"),e._v("\nA cell array of timetables containing the glucose data of the first arm. Each timetable corresponds to a patient and contains a column "),t("code",[e._v("Time")]),e._v(" and a column "),t("code",[e._v("glucose")]),e._v(" containg the glucose recordings (in mg/dl).")])]),e._v(" "),t("h3",{attrs:{id:"output-2"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#output-2"}},[e._v("#")]),e._v(" Output")]),e._v(" "),t("ul",[t("li",[t("strong",[e._v("results: structure")]),e._v(" "),t("br"),e._v("\nA structure with fields containing the computed metrics and stats in the arm, i.e.:\n"),t("ul",[t("li",[t("code",[e._v("variabilityMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed variability metrics (i.e., {"),t("code",[e._v("aucGlucose")]),e._v(", "),t("code",[e._v("CVGA")]),e._v(", "),t("code",[e._v("cvGlucose")]),e._v(", "),t("code",[e._v("efIndex")]),e._v(", "),t("code",[e._v("gmi")]),e._v(", "),t("code",[e._v("iqrGlucose")]),e._v(", "),t("code",[e._v("jIndex")]),e._v(", "),t("code",[e._v("mageIndex")]),e._v(", "),t("code",[e._v("magePlusIndex")]),e._v(", "),t("code",[e._v("mageMinusIndex")]),e._v(", "),t("code",[e._v("meanGlucose")]),e._v(", "),t("code",[e._v("medianGlucose")]),e._v(", "),t("code",[e._v("rangeGlucose")]),e._v(", "),t("code",[e._v("sddmIndex")]),e._v(", "),t("code",[e._v("sdwIndex")]),e._v(", "),t("code",[e._v("stdGlucose")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("riskMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed risk metrics (i.e., {"),t("code",[e._v("adrr")]),e._v(", "),t("code",[e._v("bgri")]),e._v(", "),t("code",[e._v("hbgi")]),e._v(", "),t("code",[e._v("lbgi")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("dataQualityMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed data quality metrics (i.e., {"),t("code",[e._v("missingGlucosePercentage")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("timeMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed time related metrics (i.e., {"),t("code",[e._v("timeInHyperglycemia")]),e._v(", "),t("code",[e._v("timeInSevereHyperglycemia")]),e._v(", "),t("code",[e._v("timeInHypoglycemia")]),e._v(", "),t("code",[e._v("timeInSevereHypoglycemia")]),e._v(", "),t("code",[e._v("timeInTarget")]),e._v(", "),t("code",[e._v("timeInTightTarget")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("glycemicTransformationMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed glycemic transformed metrics (i.e., {"),t("code",[e._v("gradeScore")]),e._v(", "),t("code",[e._v("gradeEuScore")]),e._v(", "),t("code",[e._v("gradeHyperScore")]),e._v(", "),t("code",[e._v("gradeHypoScore")]),e._v(", "),t("code",[e._v("hypoIndex")]),e._v(", "),t("code",[e._v("hyperIndex")]),e._v(", "),t("code",[e._v("igc")]),e._v(", "),t("code",[e._v("mrIndex")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("eventMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed event related metrics (i.e., {"),t("code",[e._v("gradeScore")]),e._v(", "),t("code",[e._v("gradeEuScore")]),e._v(", "),t("code",[e._v("gradeHyperScore")]),e._v(", "),t("code",[e._v("gradeHypoScore")]),e._v(", "),t("code",[e._v("hypoIndex")]),e._v(", "),t("code",[e._v("hyperIndex")]),e._v(", "),t("code",[e._v("igc")]),e._v(", "),t("code",[e._v("mrIndex")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])])])])]),e._v(" "),t("h3",{attrs:{id:"preconditions-2"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#preconditions-2"}},[e._v("#")]),e._v(" Preconditions")]),e._v(" "),t("ul",[t("li",[t("code",[e._v("arm")]),e._v(" must be a cell array containing timetables;")]),e._v(" "),t("li",[e._v("Each timetable in "),t("code",[e._v("arm")]),e._v(" must have a column names "),t("code",[e._v("Time")]),e._v(" and a\ncolumn named "),t("code",[e._v("glucose")]),e._v(".")]),e._v(" "),t("li",[e._v("Each timetable in "),t("code",[e._v("arm")]),e._v(" must have an homogeneous time grid.")])]),e._v(" "),t("h3",{attrs:{id:"reference-2"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#reference-2"}},[e._v("#")]),e._v(" Reference")]),e._v(" "),t("ul",[t("li",[e._v("None")])]),e._v(" "),t("div",{staticClass:"custom-block warning"},[t("p",{staticClass:"custom-block-title"},[e._v("WARNING")]),e._v(" "),t("p",[e._v("Currently "),t("code",[e._v("analyzeOneArm")]),e._v(" is not CI tested.")])]),e._v(" "),t("h2",{attrs:{id:"comparetwoarms"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#comparetwoarms"}},[e._v("#")]),e._v(" compareTwoArms")]),e._v(" "),t("div",{staticClass:"language-MATLAB extra-class"},[t("pre",{pre:!0,attrs:{class:"language-matlab"}},[t("code",[t("span",{pre:!0,attrs:{class:"token keyword"}},[e._v("function")]),e._v(" "),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("[")]),e._v("results"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v(" stats"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("]")]),e._v(" "),t("span",{pre:!0,attrs:{class:"token operator"}},[e._v("=")]),e._v(" "),t("span",{pre:!0,attrs:{class:"token function"}},[e._v("compareTwoArms")]),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("(")]),e._v("arm1"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("arm2"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("isPaired"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("alpha"),t("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(")")]),e._v("\n")])])]),t("p",[e._v("Function that compares the glycemic outcomes of two arms.")]),e._v(" "),t("h3",{attrs:{id:"inputs"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#inputs"}},[e._v("#")]),e._v(" Inputs")]),e._v(" "),t("ul",[t("li",[t("strong",[e._v("arm1: cell array of timetable (required)")]),e._v(" "),t("br"),e._v("\nA cell array of timetables containing the glucose data of the first arm. Each timetable corresponds to a patient and contains a column "),t("code",[e._v("Time")]),e._v(" and a column "),t("code",[e._v("glucose")]),e._v(" containg the glucose recordings (in mg/dl);")]),e._v(" "),t("li",[t("strong",[e._v("arm2: cell array of timetable (required)")]),e._v(" "),t("br"),e._v("\nA cell array of timetables containing the glucose data of the second arm. Each timetable corresponds to a patient and contains a column "),t("code",[e._v("Time")]),e._v(" and a column "),t("code",[e._v("glucose")]),e._v(" containg the glucose recordings (in mg/dl);")]),e._v(" "),t("li",[t("strong",[e._v("isPaired: integer (required)")]),e._v(" "),t("br"),e._v("\nA numeric flag defining whether to run paired or unpaired analysis. Commonly paired tests are performed when data of the same patients are present in both arms, unpaired otherwise;")]),e._v(" "),t("li",[t("strong",[e._v("alpha: double (required)")]),e._v(" "),t("br"),e._v("\nA double representing the significance level to use.")])]),e._v(" "),t("h3",{attrs:{id:"outputs"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#outputs"}},[e._v("#")]),e._v(" Outputs")]),e._v(" "),t("ul",[t("li",[t("strong",[e._v("results: structure")]),e._v(" "),t("br"),e._v("\nA structure with field "),t("code",[e._v("arm1")]),e._v(" and "),t("code",[e._v("arm2")]),e._v(", that are two structures with field containing the computed metrics in the two arms, i.e.:\n"),t("ul",[t("li",[t("code",[e._v("variabilityMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed variability metrics (i.e., {"),t("code",[e._v("aucGlucose")]),e._v(", "),t("code",[e._v("CVGA")]),e._v(", "),t("code",[e._v("cvGlucose")]),e._v(", "),t("code",[e._v("efIndex")]),e._v(", "),t("code",[e._v("gmi")]),e._v(", "),t("code",[e._v("iqrGlucose")]),e._v(", "),t("code",[e._v("jIndex")]),e._v(", "),t("code",[e._v("mageIndex")]),e._v(", "),t("code",[e._v("magePlusIndex")]),e._v(", "),t("code",[e._v("mageMinusIndex")]),e._v(", "),t("code",[e._v("meanGlucose")]),e._v(", "),t("code",[e._v("medianGlucose")]),e._v(", "),t("code",[e._v("rangeGlucose")]),e._v(", "),t("code",[e._v("sddmIndex")]),e._v(", "),t("code",[e._v("sdwIndex")]),e._v(", "),t("code",[e._v("stdGlucose")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("riskMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed risk metrics (i.e., {"),t("code",[e._v("adrr")]),e._v(", "),t("code",[e._v("bgri")]),e._v(", "),t("code",[e._v("hbgi")]),e._v(", "),t("code",[e._v("lbgi")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("dataQualityMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed data quality metrics (i.e., {"),t("code",[e._v("missingGlucosePercentage")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("timeMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed time related metrics (i.e., {"),t("code",[e._v("timeInHyperglycemia")]),e._v(", "),t("code",[e._v("timeInSevereHyperglycemia")]),e._v(", "),t("code",[e._v("timeInHypoglycemia")]),e._v(", "),t("code",[e._v("timeInSevereHypoglycemia")]),e._v(", "),t("code",[e._v("timeInTarget")]),e._v(", "),t("code",[e._v("timeInTightTarget")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("glycemicTransformationMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed glycemic transformed metrics (i.e., {"),t("code",[e._v("gradeScore")]),e._v(", "),t("code",[e._v("gradeEuScore")]),e._v(", "),t("code",[e._v("gradeHyperScore")]),e._v(", "),t("code",[e._v("gradeHypoScore")]),e._v(", "),t("code",[e._v("hypoIndex")]),e._v(", "),t("code",[e._v("hyperIndex")]),e._v(", "),t("code",[e._v("igc")]),e._v(", "),t("code",[e._v("mrIndex")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])]),e._v(" "),t("li",[t("code",[e._v("eventMetrics")]),e._v(": a structure with fields:\n"),t("ul",[t("li",[t("code",[e._v("values")]),e._v(": a vector containing the values of the computed event related metrics (i.e., {"),t("code",[e._v("gradeScore")]),e._v(", "),t("code",[e._v("gradeEuScore")]),e._v(", "),t("code",[e._v("gradeHyperScore")]),e._v(", "),t("code",[e._v("gradeHypoScore")]),e._v(", "),t("code",[e._v("hypoIndex")]),e._v(", "),t("code",[e._v("hyperIndex")]),e._v(", "),t("code",[e._v("igc")]),e._v(", "),t("code",[e._v("mrIndex")]),e._v("}) of the metrics for each glucose profile;")]),e._v(" "),t("li",[t("code",[e._v("mean")]),e._v(": the mean of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("median")]),e._v(": the median of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("std")]),e._v(": the standard deviation of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc5")]),e._v(": the 5th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc25")]),e._v(": the 25th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc75")]),e._v(": the 75th percentile of "),t("code",[e._v("values")]),e._v(";")]),e._v(" "),t("li",[t("code",[e._v("prc95")]),e._v(": the 95th percentile of "),t("code",[e._v("values")]),e._v(";")])])])])]),e._v(" "),t("li",[t("strong",[e._v("stats: structure")]),e._v(" "),t("br"),e._v("\nA structure that contains for each of the considered metrics the result of the statistical test with field "),t("code",[e._v("p")]),e._v(" (p-value value) and "),t("code",[e._v("h")]),e._v(" null hypothesis accepted or rejcted. Statistical tests are:\n"),t("ul",[t("li",[t("em",[e._v("t-test")]),e._v(" if the test "),t("code",[e._v("isPaired")]),e._v(" and the samples are both gaussian distributed (checked with the Lilliefors test);")]),e._v(" "),t("li",[t("em",[e._v("unpaired t-test")]),e._v(" if the test not "),t("code",[e._v("isPaired")]),e._v(" and the samples are both gaussian distributed (checked with the Lilliefors test);")]),e._v(" "),t("li",[t("em",[e._v("Wilcoxon rank test")]),e._v(" if the test "),t("code",[e._v("isPaired")]),e._v(" and at least one of the samples is not gaussian distributed (checked with the Lilliefors test);")]),e._v(" "),t("li",[t("em",[e._v("Mann-Whitney U-test")]),e._v(" if the test not "),t("code",[e._v("isPaired")]),e._v(" and at least one of the samples is not gaussian distributed (checked with the Lilliefors test).")])])])]),e._v(" "),t("h3",{attrs:{id:"preconditions-3"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#preconditions-3"}},[e._v("#")]),e._v(" Preconditions")]),e._v(" "),t("ul",[t("li",[t("code",[e._v("arm1")]),e._v(" must be a cell array containing timetables;")]),e._v(" "),t("li",[t("code",[e._v("arm2")]),e._v(" must be a cell array containing timetables;")]),e._v(" "),t("li",[e._v("Each timetable in "),t("code",[e._v("arm1")]),e._v(" and "),t("code",[e._v("arm2")]),e._v(" must have a column names "),t("code",[e._v("Time")]),e._v(" and a\ncolumn named "),t("code",[e._v("glucose")]),e._v(".")]),e._v(" "),t("li",[e._v("Each timetable in "),t("code",[e._v("arm1")]),e._v(" and "),t("code",[e._v("arm2")]),e._v(" must have an homogeneous time grid;")]),e._v(" "),t("li",[t("code",[e._v("isPaired")]),e._v(" can be 0 or 1.")])]),e._v(" "),t("h3",{attrs:{id:"reference-3"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#reference-3"}},[e._v("#")]),e._v(" Reference")]),e._v(" "),t("ul",[t("li",[e._v('Lilliefors et al., "On the Kolmogorov-Smirnov test for normality with mean and variance unknown," Mathematics, vol. 62, 1967, pp. 399–402. DOI: 10.1080/01621459.1967.10482916.')])]),e._v(" "),t("div",{staticClass:"custom-block warning"},[t("p",{staticClass:"custom-block-title"},[e._v("WARNING")]),e._v(" "),t("p",[e._v("Currently "),t("code",[e._v("compareTwoArms")]),e._v(" is not CI tested.")])])])}),[],!1,null,null,null);v.default=o.exports}}]);