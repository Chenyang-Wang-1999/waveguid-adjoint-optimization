/*
 * cross_waveguide_test.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 19 2026, 00:24 by COMSOL 5.5.0.359. */
public class cross_waveguide_test {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/home/wangchenyang/654/research-data/CouplingRings/COMSOL");

    model.label("cross_waveguide_test.mph");

    model.param().set("L_unit", "1[um]");
    model.param().set("T_unit", "L_unit/c_const");
    model.param().set("P_unit", "1[W]");
    model.param().set("E_unit", "P_unit*T_unit");
    model.param().set("S_unit", "P_unit/L_unit^2");
    model.param().group().create("par2");
    model.param("par2").set("L", "10");
    model.param("par2").set("w", "1");
    model.param("par2").set("f_sweep", "0.1");
    model.param().group().create("par3");
    model.param("par3").set("e1x", "cos(pi/6)");
    model.param("par3").set("e1y", "sin(pi/6)");
    model.param("par3").set("e2x", "-sin(pi/6)");
    model.param("par3").set("e2y", "cos(pi/6)");
    model.param().group().create("par4");
    model.param("par4").set("epsilon_wg", "13");
    model.param("par4").set("epsilon_amb", "1");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 2);

    model.result().table().create("tbl1", "Table");
    model.result().table().create("tbl2", "Table");
    model.result().table().create("tbl3", "Table");
    model.result().table().create("tbl4", "Table");
    model.result().table().create("tbl5", "Table");
    model.result().table().create("tbl6", "Table");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").create("pol1", "Polygon");
    model.component("comp1").geom("geom1").feature("pol1").set("x", "L*cos(2*pi/6*range(0,1,6))*L_unit");
    model.component("comp1").geom("geom1").feature("pol1").set("y", "L*sin(2*pi/6*range(0,1,6))*L_unit");
    model.component("comp1").geom("geom1").create("r1", "Rectangle");
    model.component("comp1").geom("geom1").feature("r1").set("rot", 30);
    model.component("comp1").geom("geom1").feature("r1").set("base", "center");
    model.component("comp1").geom("geom1").feature("r1")
         .set("size", new String[]{"2*L*cos(pi/6)*L_unit", "w*L_unit"});
    model.component("comp1").geom("geom1").create("r2", "Rectangle");
    model.component("comp1").geom("geom1").feature("r2").set("rot", -30);
    model.component("comp1").geom("geom1").feature("r2").set("base", "center");
    model.component("comp1").geom("geom1").feature("r2")
         .set("size", new String[]{"2*L*cos(pi/6)*L_unit", "w*L_unit"});
    model.component("comp1").geom("geom1").create("ls1", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls1").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls1")
         .set("coord1", new String[]{"L_unit*(L*cos(pi/6)*e1x + 3*w*e2x)", "L_unit*(L*cos(pi/6)*e1y + 3*w*e2y)"});
    model.component("comp1").geom("geom1").feature("ls1").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls1")
         .set("coord2", new String[]{"L_unit*(L*cos(pi/6)*e1x - 3*w*e2x)", "L_unit*(L*cos(pi/6)*e1y - 3*w*e2y)"});
    model.component("comp1").geom("geom1").create("rot1", "Rotate");
    model.component("comp1").geom("geom1").feature("rot1").set("keep", true);
    model.component("comp1").geom("geom1").feature("rot1").set("rot", "120, 180,-60");
    model.component("comp1").geom("geom1").feature("rot1").selection("input").set("ls1");
    model.component("comp1").geom("geom1").run();

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat1").selection().set(1, 4, 5, 9);
    model.component("comp1").material("mat1").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat2").selection().set(2, 3, 6, 7, 8);
    model.component("comp1").material("mat2").propertyGroup().create("RefractiveIndex", "Refractive index");

    model.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1");
    model.component("comp1").physics("ewfd").create("port1", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port1").selection().set(27, 29, 31);
    model.component("comp1").physics("ewfd").create("port2", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port2").selection().set(4, 7, 11);
    model.component("comp1").physics("ewfd").create("port3", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port3").selection().set(28, 30, 32);
    model.component("comp1").physics("ewfd").create("port4", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port4").selection().set(3, 5, 9);
    model.component("comp1").physics("ewfd").create("sctr1", "Scattering", 1);
    model.component("comp1").physics("ewfd").feature("sctr1").selection().set(1, 2, 13, 14, 15, 16, 25, 26, 33, 34);

    model.component("comp1").mesh("mesh1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").create("size2", "Size");
    model.component("comp1").mesh("mesh1").create("size3", "Size");
    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").feature("size1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("size1").selection().set(2, 3, 6, 7, 8);
    model.component("comp1").mesh("mesh1").feature("size2").selection().geom("geom1", 1);
    model.component("comp1").mesh("mesh1").feature("size2").selection()
         .set(3, 4, 5, 7, 9, 11, 27, 28, 29, 30, 31, 32);
    model.component("comp1").mesh("mesh1").feature("size3").selection().geom("geom1", 1);
    model.component("comp1").mesh("mesh1").feature("size3").selection()
         .set(3, 4, 5, 7, 9, 11, 27, 28, 29, 30, 31, 32);

    model.result().table("tbl2").comments("Line Integration 1 {int1}");
    model.result().table("tbl3").comments("Line Integration 3 {int3}");
    model.result().table("tbl4").comments("Line Integration 2 {int2}");
    model.result().table("tbl5").comments("Line Integration 4 {int4}");
    model.result().table("tbl6").comments("Global Evaluation 3 {gev3}");

    model.component("comp1").view("view1").axis().set("xmin", -2.6687093850341626E-5);
    model.component("comp1").view("view1").axis().set("xmax", 2.6687093850341626E-5);
    model.component("comp1").view("view1").axis().set("ymin", -1.3129599210515153E-5);
    model.component("comp1").view("view1").axis().set("ymax", 1.3129599210515153E-5);

    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex")
         .set("n", new String[]{"sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)"});
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex")
         .set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n_symmetry", "0");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex")
         .set("n", new String[]{"sqrt(epsilon_wg)", "0", "0", "0", "sqrt(epsilon_wg)", "0", "0", "0", "sqrt(epsilon_wg)"});
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex")
         .set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("n_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("ki_symmetry", "0");

    model.component("comp1").physics("ewfd").prop("components").set("components", "outofplane");
    model.component("comp1").physics("ewfd").feature("port1").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port1").set("PortFreq", "f_sweep/T_unit");
    model.component("comp1").physics("ewfd").feature("port1")
         .set("LambdaPort", "-278.8560239932588-1649804.531029822i");
    model.component("comp1").physics("ewfd").feature("port2").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port2").set("PortFreq", "f_sweep/T_unit");
    model.component("comp1").physics("ewfd").feature("port2")
         .set("LambdaPort", "-278.85602399328275-1649804.5310298116i");
    model.component("comp1").physics("ewfd").feature("port3").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port3").set("PortFreq", "f_sweep/T_unit");
    model.component("comp1").physics("ewfd").feature("port3")
         .set("LambdaPort", "-278.85602399342883-1649804.5310298111i");
    model.component("comp1").physics("ewfd").feature("port4").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port4").set("PortFreq", "f_sweep/T_unit");
    model.component("comp1").physics("ewfd").feature("port4")
         .set("LambdaPort", "-278.8560239932169-1649804.5310298244i");

    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "1.2500E-6");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "6.0000E-9");
    model.component("comp1").mesh("mesh1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size1").set("hmax", "3.4670E-7");
    model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hmin", "6.0000E-9");
    model.component("comp1").mesh("mesh1").feature("size1").set("hminactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hcurveactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hgradactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size2").set("hmax", "6.2500E-7");
    model.component("comp1").mesh("mesh1").feature("size2").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hmin", "6.0000E-9");
    model.component("comp1").mesh("mesh1").feature("size2").set("hminactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hcurveactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hgradactive", true);
    model.component("comp1").mesh("mesh1").feature("size3").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size3").set("hmax", "w*L_unit/10");
    model.component("comp1").mesh("mesh1").feature("size3").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").run();

    model.study().create("std1");
    model.study("std1").create("bma", "BoundaryModeAnalysis");
    model.study("std1").create("bma1", "BoundaryModeAnalysis");
    model.study("std1").create("bma2", "BoundaryModeAnalysis");
    model.study("std1").create("bma3", "BoundaryModeAnalysis");
    model.study("std1").create("freq", "Frequency");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("e1", "Eigenvalue");
    model.sol("sol1").create("su1", "StoreSolution");
    model.sol("sol1").create("st2", "StudyStep");
    model.sol("sol1").create("v2", "Variables");
    model.sol("sol1").create("e2", "Eigenvalue");
    model.sol("sol1").create("su2", "StoreSolution");
    model.sol("sol1").create("st3", "StudyStep");
    model.sol("sol1").create("v3", "Variables");
    model.sol("sol1").create("e3", "Eigenvalue");
    model.sol("sol1").create("su3", "StoreSolution");
    model.sol("sol1").create("st4", "StudyStep");
    model.sol("sol1").create("v4", "Variables");
    model.sol("sol1").create("e4", "Eigenvalue");
    model.sol("sol1").create("su4", "StoreSolution");
    model.sol("sol1").create("st5", "StudyStep");
    model.sol("sol1").create("v5", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("e1").create("d1", "Direct");
    model.sol("sol1").feature("e2").create("d1", "Direct");
    model.sol("sol1").feature("e3").create("d1", "Direct");
    model.sol("sol1").feature("e4").create("d1", "Direct");
    model.sol("sol1").feature("s1").create("p1", "Parametric");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");

    model.result().numerical().create("gev1", "EvalGlobal");
    model.result().numerical().create("int1", "IntLine");
    model.result().numerical().create("gev2", "EvalGlobal");
    model.result().numerical().create("int2", "IntLine");
    model.result().numerical().create("int3", "IntLine");
    model.result().numerical().create("int4", "IntLine");
    model.result().numerical().create("gev3", "EvalGlobal");
    model.result().numerical("gev1").set("probetag", "none");
    model.result().numerical("int1").selection().set(27, 29, 31);
    model.result().numerical("int1").set("probetag", "none");
    model.result().numerical("gev2").set("probetag", "none");
    model.result().numerical("int2").selection().set(4, 7, 11);
    model.result().numerical("int2").set("probetag", "none");
    model.result().numerical("int3").selection().set(28, 30, 32);
    model.result().numerical("int3").set("probetag", "none");
    model.result().numerical("int4").selection().set(3, 5, 9);
    model.result().numerical("int4").set("probetag", "none");
    model.result().numerical("gev3").set("probetag", "none");
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").create("surf1", "Surface");

    model.study("std1").feature("bma").set("modeFreq", "f_sweep/T_unit");
    model.study("std1").feature("bma").set("shift", "3");
    model.study("std1").feature("bma1").set("PortName", "2");
    model.study("std1").feature("bma1").set("modeFreq", "f_sweep/T_unit");
    model.study("std1").feature("bma1").set("shift", "3");
    model.study("std1").feature("bma2").set("PortName", "3");
    model.study("std1").feature("bma2").set("modeFreq", "f_sweep/T_unit");
    model.study("std1").feature("bma2").set("shift", "3");
    model.study("std1").feature("bma3").set("PortName", "4");
    model.study("std1").feature("bma3").set("modeFreq", "f_sweep/T_unit");
    model.study("std1").feature("bma3").set("shift", "3");
    model.study("std1").feature("freq").set("plist", "f_sweep/T_unit");
    model.study("std1").feature("freq").set("ngen", 5);

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("e1").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e1").set("neigs", 1);
    model.sol("sol1").feature("e1").set("shift", "3");
    model.sol("sol1").feature("e1").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e1").feature("d1").label("Suggested Direct Solver (ewfd)");
    model.sol("sol1").feature("st2").set("studystep", "bma1");
    model.sol("sol1").feature("v2").set("initmethod", "sol");
    model.sol("sol1").feature("v2").set("initsol", "sol1");
    model.sol("sol1").feature("v2").set("solnum", "auto");
    model.sol("sol1").feature("v2").set("notsolmethod", "sol");
    model.sol("sol1").feature("v2").set("notsol", "sol1");
    model.sol("sol1").feature("v2").set("notsolnum", "auto");
    model.sol("sol1").feature("e2").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e2").set("neigs", 1);
    model.sol("sol1").feature("e2").set("shift", "3");
    model.sol("sol1").feature("e2").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e2").feature("d1").label("Suggested Direct Solver (ewfd)");
    model.sol("sol1").feature("st3").set("studystep", "bma2");
    model.sol("sol1").feature("v3").set("initmethod", "sol");
    model.sol("sol1").feature("v3").set("initsol", "sol1");
    model.sol("sol1").feature("v3").set("solnum", "auto");
    model.sol("sol1").feature("v3").set("notsolmethod", "sol");
    model.sol("sol1").feature("v3").set("notsol", "sol1");
    model.sol("sol1").feature("v3").set("notsolnum", "auto");
    model.sol("sol1").feature("e3").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e3").set("neigs", 1);
    model.sol("sol1").feature("e3").set("shift", "3");
    model.sol("sol1").feature("e3").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e3").feature("d1").label("Suggested Direct Solver (ewfd)");
    model.sol("sol1").feature("st4").set("studystep", "bma3");
    model.sol("sol1").feature("v4").set("initmethod", "sol");
    model.sol("sol1").feature("v4").set("initsol", "sol1");
    model.sol("sol1").feature("v4").set("solnum", "auto");
    model.sol("sol1").feature("v4").set("notsolmethod", "sol");
    model.sol("sol1").feature("v4").set("notsol", "sol1");
    model.sol("sol1").feature("v4").set("notsolnum", "auto");
    model.sol("sol1").feature("e4").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e4").set("neigs", 1);
    model.sol("sol1").feature("e4").set("shift", "3");
    model.sol("sol1").feature("e4").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e4").feature("d1").label("Suggested Direct Solver (ewfd)");
    model.sol("sol1").feature("st5").set("studystep", "freq");
    model.sol("sol1").feature("v5").set("initmethod", "sol");
    model.sol("sol1").feature("v5").set("initsol", "sol1");
    model.sol("sol1").feature("v5").set("solnum", "auto");
    model.sol("sol1").feature("v5").set("notsolmethod", "sol");
    model.sol("sol1").feature("v5").set("notsol", "sol1");
    model.sol("sol1").feature("v5").set("notsolnum", "auto");
    model.sol("sol1").feature("v5").set("clistctrl", new String[]{"p1"});
    model.sol("sol1").feature("v5").set("cname", new String[]{"freq"});
    model.sol("sol1").feature("v5").set("clist", new String[]{"f_sweep/T_unit"});
    model.sol("sol1").feature("s1").set("control", "freq");
    model.sol("sol1").feature("s1").set("stol", 0.01);
    model.sol("sol1").feature("s1").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("s1").feature("p1").set("control", "freq");
    model.sol("sol1").feature("s1").feature("p1").set("pname", new String[]{"freq"});
    model.sol("sol1").feature("s1").feature("p1").set("plistarr", new String[]{"f_sweep/T_unit"});
    model.sol("sol1").feature("s1").feature("p1").set("punit", new String[]{"THz"});
    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");
    model.sol("sol1").feature("s1").feature("p1").set("preusesol", "auto");
    model.sol("sol1").feature("s1").feature("d1").label("Suggested Direct Solver (ewfd)");
    model.sol("sol1").runAll();

    model.result().numerical("gev1").label("Reflectance, Transmittance, and Absorptance (ewfd)");
    model.result().numerical("gev1").set("table", "tbl1");
    model.result().numerical("gev1")
         .set("expr", new String[]{"ewfd.Rport_1", "ewfd.Rtotal", "ewfd.Tport_2", "ewfd.Tport_3", "ewfd.Tport_4", "ewfd.Ttotal", "ewfd.RTtotal", "ewfd.Atotal"});
    model.result().numerical("gev1").set("unit", new String[]{"1", "1", "1", "1", "1", "1", "1", "1"});
    model.result().numerical("gev1")
         .set("descr", new String[]{"Reflectance, port 1", "Total reflectance", "Transmittance, port 2", "Transmittance, port 3", "Transmittance, port 4", "Total transmittance", "Total reflectance and transmittance", "Absorptance"});
    model.result().numerical("int1").set("table", "tbl2");
    model.result().numerical("int1")
         .set("expr", new String[]{"(epsilon0_const*ewfd.epsilonrxx*ewfd.normtEmode_1^2+ewfd.murxx*mu0_const*ewfd.normtHmode_1^2)/(E_unit/L_unit^2)/4", ""});
    model.result().numerical("int1").set("unit", new String[]{"1", ""});
    model.result().numerical("int1").set("descr", new String[]{"", ""});
    model.result().numerical("int2").set("table", "tbl4");
    model.result().numerical("int2").set("expr", new String[]{"ewfd.nPoav/(S_unit*L_unit)/3.126339348168692E-6"});
    model.result().numerical("int2").set("unit", new String[]{"1"});
    model.result().numerical("int2").set("descr", new String[]{""});
    model.result().numerical("int3").set("table", "tbl3");
    model.result().numerical("int3").set("expr", new String[]{"ewfd.nPoav/(S_unit*L_unit)/3.126339348168692E-6"});
    model.result().numerical("int3").set("unit", new String[]{"1"});
    model.result().numerical("int3").set("descr", new String[]{""});
    model.result().numerical("int4").set("table", "tbl5");
    model.result().numerical("int4").set("expr", new String[]{"ewfd.nPoav/(S_unit*L_unit)/3.126339348168692E-6"});
    model.result().numerical("int4").set("unit", new String[]{"1"});
    model.result().numerical("int4").set("descr", new String[]{""});
    model.result().numerical("gev3").set("table", "tbl6");
    model.result().numerical("gev3").set("expr", new String[]{"ewfd.S11", "ewfd.S21", "ewfd.S31", "ewfd.S41"});
    model.result().numerical("gev3").set("unit", new String[]{"1", "1", "1", "1"});
    model.result().numerical("gev3").set("descr", new String[]{"S11", "S21", "S31", "S41"});
    model.result().numerical("gev1").setResult();
    model.result().numerical("int1").setResult();
    model.result().numerical("int2").setResult();
    model.result().numerical("int3").setResult();
    model.result().numerical("int4").setResult();
    model.result().numerical("gev3").setResult();
    model.result("pg1").label("Electric Field (ewfd)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").feature("surf1").set("colortable", "RainbowLight");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("resolution", "normal");

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
