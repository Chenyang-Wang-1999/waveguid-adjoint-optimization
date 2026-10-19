/*
 * CouplingTriangles_double_triangle_wide_bench.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 19 2026, 19:22 by COMSOL 5.5.0.359. */
public class CouplingTriangles_double_triangle_wide_bench {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/home/wangchenyang/654/research-data/FEM/adjoint_optimization/comsol");

    model.label("CouplingTriangles_double_triangle_wide_bench.mph");

    model.param().set("r_ex", "L_unit");
    model.param().set("d_ring", "0.8*L_unit");
    model.param().set("w_ring", "0.1*L_unit");
    model.param().set("a", "(2*r_ex + d_ring)*sqrt(3)");
    model.param().set("h_ring", "0.1*a");
    model.param().set("r_in", "r_ex - w_ring");
    model.param().set("d_triangle", "0.08*L_unit");
    model.param().set("rho_vertex", "0.085*a");
    model.param().set("l_triangle", "a/sqrt(3) - r_ex - rho_vertex");
    model.param().set("rho_vertex_in", "rho_vertex - d_triangle");
    model.param().set("r_air", "8*r_ex");
    model.param().set("approx_f", "c_const/a");
    model.param().set("inverse_r_ex", "5[um]/r_ex/cos(pi/6)");
    model.param().set("rho_wg_cen", "(rho_vertex+rho_vertex_in)/2");
    model.param().group().create("par2");
    model.param("par2").set("n_air", "1");
    model.param("par2").set("n_wg", "3.4");
    model.param().group().create("par3");
    model.param("par3").set("nu0", "1.552169");
    model.param("par3").set("d_coupling", "w_ring", "waveguide ring coupling changed");
    model.param("par3").set("f_sweep", "1");
    model.param("par3").set("w_port", "6*w_ring");
    model.param("par3").set("d_ring_tri", "0.5*w_ring");
    model.param().group().create("par4");
    model.param("par4").set("L_unit", "1[um]");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 2);

    model.result().table().create("tbl1", "Table");

    model.component("comp1").func().create("an1", "Analytic");
    model.component("comp1").func("an1").active(false);
    model.component("comp1").func("an1").set("funcname", "intair");
    model.component("comp1").func("an1").set("expr", "intall(x) - intwg(x)");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").create("c1", "Circle");
    model.component("comp1").geom("geom1").feature("c1").set("pos", new String[]{"- a/2", "- a/2/sqrt(3)"});
    model.component("comp1").geom("geom1").feature("c1").set("r", "r_ex");
    model.component("comp1").geom("geom1").create("c2", "Circle");
    model.component("comp1").geom("geom1").feature("c2").set("pos", new String[]{"- a/2", "- a/2/sqrt(3)"});
    model.component("comp1").geom("geom1").feature("c2").set("r", "r_in");
    model.component("comp1").geom("geom1").create("r1", "Rectangle");
    model.component("comp1").geom("geom1").feature("r1")
         .set("pos", new String[]{"- a/2-r_ex-w_ring/2-d_coupling", "0"});
    model.component("comp1").geom("geom1").feature("r1").set("base", "center");
    model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"w_ring", "r_air"});
    model.component("comp1").geom("geom1").create("dif1", "Difference");
    model.component("comp1").geom("geom1").feature("dif1").selection("input").set("c1");
    model.component("comp1").geom("geom1").feature("dif1").selection("input2").set("c2");
    model.component("comp1").geom("geom1").create("rot1", "Rotate");
    model.component("comp1").geom("geom1").feature("rot1").set("keep", true);
    model.component("comp1").geom("geom1").feature("rot1").setIndex("rot", "120", 0);
    model.component("comp1").geom("geom1").feature("rot1").selection("input").set("dif1");
    model.component("comp1").geom("geom1").create("c3", "Circle");
    model.component("comp1").geom("geom1").feature("c3").set("pos", new String[]{"0", "l_triangle-d_ring_tri"});
    model.component("comp1").geom("geom1").feature("c3").set("rot", 30);
    model.component("comp1").geom("geom1").feature("c3").set("r", "rho_vertex");
    model.component("comp1").geom("geom1").feature("c3").set("angle", 120);
    model.component("comp1").geom("geom1").create("c4", "Circle");
    model.component("comp1").geom("geom1").feature("c4").set("pos", new String[]{"0", "l_triangle-d_ring_tri"});
    model.component("comp1").geom("geom1").feature("c4").set("rot", 30);
    model.component("comp1").geom("geom1").feature("c4").set("r", "rho_vertex_in");
    model.component("comp1").geom("geom1").feature("c4").set("angle", 120);
    model.component("comp1").geom("geom1").create("dif2", "Difference");
    model.component("comp1").geom("geom1").feature("dif2").selection("input").set("c3");
    model.component("comp1").geom("geom1").feature("dif2").selection("input2").set("c4");
    model.component("comp1").geom("geom1").create("rot2", "Rotate");
    model.component("comp1").geom("geom1").feature("rot2").set("keep", true);
    model.component("comp1").geom("geom1").feature("rot2").set("rot", "range(60,60,300)");
    model.component("comp1").geom("geom1").feature("rot2").selection("input").set("dif2");
    model.component("comp1").geom("geom1").create("ls1", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex1").set("dif2(1)", 1);
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex2").set("rot2(2)", 6);
    model.component("comp1").geom("geom1").create("ls2", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls2").selection("vertex1").set("dif2(1)", 2);
    model.component("comp1").geom("geom1").feature("ls2").selection("vertex2").set("rot2(2)", 5);
    model.component("comp1").geom("geom1").create("ls3", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls3").selection("vertex1").set("dif2(1)", 5);
    model.component("comp1").geom("geom1").feature("ls3").selection("vertex2").set("rot2(4)", 2);
    model.component("comp1").geom("geom1").create("ls4", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls4").selection("vertex1").set("dif2(1)", 6);
    model.component("comp1").geom("geom1").feature("ls4").selection("vertex2").set("rot2(4)", 1);
    model.component("comp1").geom("geom1").create("ls11", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls11").selection("vertex1").set("rot2(1)", 6);
    model.component("comp1").geom("geom1").feature("ls11").selection("vertex2").set("rot2(5)", 1);
    model.component("comp1").geom("geom1").create("ls12", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls12").selection("vertex1").set("rot2(1)", 5);
    model.component("comp1").geom("geom1").feature("ls12").selection("vertex2").set("rot2(5)", 2);
    model.component("comp1").geom("geom1").create("ls13", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls13").selection("vertex1").set("rot2(1)", 1);
    model.component("comp1").geom("geom1").feature("ls13").selection("vertex2").set("rot2(3)", 6);
    model.component("comp1").geom("geom1").create("ls14", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls14").selection("vertex1").set("rot2(1)", 2);
    model.component("comp1").geom("geom1").feature("ls14").selection("vertex2").set("rot2(3)", 5);
    model.component("comp1").geom("geom1").create("ls15", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls15").selection("vertex1").set("rot2(3)", 2);
    model.component("comp1").geom("geom1").feature("ls15").selection("vertex2").set("rot2(5)", 5);
    model.component("comp1").geom("geom1").create("ls16", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls16").selection("vertex1").set("rot2(5)", 6);
    model.component("comp1").geom("geom1").feature("ls16").selection("vertex2").set("rot2(3)", 1);
    model.component("comp1").geom("geom1").create("ls17", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls17").selection("vertex1").set("rot2(4)", 5);
    model.component("comp1").geom("geom1").feature("ls17").selection("vertex2").set("rot2(2)", 2);
    model.component("comp1").geom("geom1").create("ls18", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls18").selection("vertex1").set("rot2(2)", 1);
    model.component("comp1").geom("geom1").feature("ls18").selection("vertex2").set("rot2(4)", 6);
    model.component("comp1").geom("geom1").create("mir1", "Mirror");
    model.component("comp1").geom("geom1").feature("mir1").set("keep", true);
    model.component("comp1").geom("geom1").feature("mir1").selection("input").set("r1");
    model.component("comp1").geom("geom1").create("r2", "Rectangle");
    model.component("comp1").geom("geom1").feature("r2").set("base", "center");
    model.component("comp1").geom("geom1").feature("r2").set("size", new String[]{"2*r_air", "r_air"});
    model.component("comp1").geom("geom1").create("ls7", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls7").label("port 4");
    model.component("comp1").geom("geom1").feature("ls7").set("selresult", true);
    model.component("comp1").geom("geom1").feature("ls7").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls7")
         .set("coord1", new String[]{"- a/2-r_ex-w_ring/2-d_coupling-w_port/2", "r_air/2"});
    model.component("comp1").geom("geom1").feature("ls7").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls7")
         .set("coord2", new String[]{"- a/2-r_ex-w_ring/2-d_coupling+w_port/2", "r_air/2"});
    model.component("comp1").geom("geom1").create("ls8", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls8").label("port 3");
    model.component("comp1").geom("geom1").feature("ls8").set("selresult", true);
    model.component("comp1").geom("geom1").feature("ls8").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls8")
         .set("coord1", new String[]{"- a/2-r_ex-w_ring/2-d_coupling-w_port/2", "-r_air/2"});
    model.component("comp1").geom("geom1").feature("ls8").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls8")
         .set("coord2", new String[]{"- a/2-r_ex-w_ring/2-d_coupling+w_port/2", "-r_air/2"});
    model.component("comp1").geom("geom1").create("ls9", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls9").label("port 1");
    model.component("comp1").geom("geom1").feature("ls9").set("selresult", true);
    model.component("comp1").geom("geom1").feature("ls9").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls9")
         .set("coord1", new String[]{"-(- a/2-r_ex-w_ring/2-d_coupling-w_port/2)", "-r_air/2"});
    model.component("comp1").geom("geom1").feature("ls9").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls9")
         .set("coord2", new String[]{"-(- a/2-r_ex-w_ring/2-d_coupling+w_port/2)", "-r_air/2"});
    model.component("comp1").geom("geom1").create("ls10", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls10").label("port 2");
    model.component("comp1").geom("geom1").feature("ls10").set("selresult", true);
    model.component("comp1").geom("geom1").feature("ls10").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls10")
         .set("coord1", new String[]{"-(- a/2-r_ex-w_ring/2-d_coupling-w_port/2)", "r_air/2"});
    model.component("comp1").geom("geom1").feature("ls10").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls10")
         .set("coord2", new String[]{"-(- a/2-r_ex-w_ring/2-d_coupling+w_port/2)", "r_air/2"});
    model.component("comp1").geom("geom1").create("c6", "Circle");
    model.component("comp1").geom("geom1").feature("c6").set("pos", new String[]{"- a/2", "a/2/sqrt(3)"});
    model.component("comp1").geom("geom1").feature("c6").set("r", "r_ex");
    model.component("comp1").geom("geom1").create("c7", "Circle");
    model.component("comp1").geom("geom1").feature("c7").set("pos", new String[]{"- a/2", "a/2/sqrt(3)"});
    model.component("comp1").geom("geom1").feature("c7").set("r", "r_in");
    model.component("comp1").geom("geom1").run();

    model.view().create("view2", 3);

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat1").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat2").selection()
         .set(2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 36, 37, 38, 39, 40, 41, 43, 44, 47);
    model.component("comp1").material("mat2").propertyGroup().create("RefractiveIndex", "Refractive index");

    model.component("comp1").cpl().create("intop1", "Integration");
    model.component("comp1").cpl().create("intop2", "Integration");
    model.component("comp1").cpl("intop2").selection().all();

    model.component("comp1").physics().create("ewfd_TE", "ElectromagneticWavesFrequencyDomain", "geom1");
    model.component("comp1").physics("ewfd_TE").identifier("ewfd_TE");
    model.component("comp1").physics("ewfd_TE").create("sctr1", "Scattering", 1);
    model.component("comp1").physics("ewfd_TE").feature("sctr1").selection().all();
    model.component("comp1").physics().create("ewfd_TM", "ElectromagneticWavesFrequencyDomain", "geom1");
    model.component("comp1").physics("ewfd_TM").identifier("ewfd_TM");
    model.component("comp1").physics("ewfd_TM").create("sctr1", "Scattering", 1);
    model.component("comp1").physics("ewfd_TM").feature("sctr1").selection()
         .set(1, 2, 3, 6, 9, 12, 13, 15, 18, 21, 22, 25, 30, 33, 36, 37, 41, 44, 45, 48, 55, 57, 60, 64, 65, 68, 72, 73, 74, 79, 80, 82, 88, 91, 94, 95, 96, 97, 98, 101, 102, 105, 106, 107, 108, 113, 114, 115, 116, 117, 118, 119, 120, 121, 123, 125, 126, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144);
    model.component("comp1").physics("ewfd_TM").create("imp1", "Impedance", 1);
    model.component("comp1").physics("ewfd_TM").feature("imp1").selection()
         .set(1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 86, 87, 89, 90, 92, 93, 94, 95, 96);
    model.component("comp1").physics("ewfd_TM").create("port1", "Port", 1);
    model.component("comp1").physics("ewfd_TM").feature("port1").selection().named("geom1_ls9_bnd");
    model.component("comp1").physics("ewfd_TM").create("port2", "Port", 1);
    model.component("comp1").physics("ewfd_TM").feature("port2").selection().named("geom1_ls10_bnd");
    model.component("comp1").physics("ewfd_TM").create("port3", "Port", 1);
    model.component("comp1").physics("ewfd_TM").feature("port3").selection().named("geom1_ls8_bnd");
    model.component("comp1").physics("ewfd_TM").create("port4", "Port", 1);
    model.component("comp1").physics("ewfd_TM").feature("port4").selection().named("geom1_ls7_bnd");

    model.component("comp1").mesh("mesh1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").create("size2", "Size");
    model.component("comp1").mesh("mesh1").create("ftri3", "FreeTri");
    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").feature("size1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("size1").selection()
         .set(2, 4, 5, 6, 7, 10, 12, 13, 14, 17, 18, 19, 21, 24, 25, 26, 27, 29, 31, 32, 33, 36, 37, 38, 39, 40, 41, 43, 47);
    model.component("comp1").mesh("mesh1").feature("size2").selection().geom("geom1", 1);
    model.component("comp1").mesh("mesh1").feature("size2").selection()
         .set(4, 5, 7, 8, 10, 11, 86, 87, 89, 90, 92, 93);
    model.component("comp1").mesh("mesh1").feature("ftri3").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftri3").selection().set(3);

    model.result().table("tbl1").comments("\u5168\u5c40\u8ba1\u7b97 1");

    model.component("comp1").view("view1").set("showmaterial", true);
    model.component("comp1").view("view1").axis().set("xmin", -1.0428158020658884E-5);
    model.component("comp1").view("view1").axis().set("xmax", 1.0428158020658884E-5);
    model.component("comp1").view("view1").axis().set("ymin", -4.399999852466863E-6);
    model.component("comp1").view("view1").axis().set("ymax", 4.399999852466863E-6);

    model.component("comp1").material("mat1").set("family", "air");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex")
         .set("n", new String[]{"n_air", "0", "0", "0", "n_air", "0", "0", "0", "n_air"});
    model.component("comp1").material("mat1").propertyGroup("RefractiveIndex")
         .set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat2").set("family", "gold");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex")
         .set("n", new String[]{"n_wg", "0", "0", "0", "n_wg", "0", "0", "0", "n_wg"});
    model.component("comp1").material("mat2").propertyGroup("RefractiveIndex")
         .set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});

    model.component("comp1").cpl("intop1").set("opname", "intwg");
    model.component("comp1").cpl("intop2").set("opname", "intall");

    model.component("comp1").physics("ewfd_TM").label("TM modes");
    model.component("comp1").physics("ewfd_TM").prop("components").set("components", "outofplane");
    model.component("comp1").physics("ewfd_TM").feature("imp1").active(false);
    model.component("comp1").physics("ewfd_TM").feature("port1").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port1").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port1")
         .set("LambdaPort", "-13.545027907991829-2.7391341646235574E7i");
    model.component("comp1").physics("ewfd_TM").feature("port2").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port2").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port2")
         .set("LambdaPort", "-13.545027909466-2.7391341646235634E7i");
    model.component("comp1").physics("ewfd_TM").feature("port3").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port3").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port3")
         .set("LambdaPort", "-13.5450279138884-2.7391341646235753E7i");
    model.component("comp1").physics("ewfd_TM").feature("port4").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port4").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port4")
         .set("LambdaPort", "-13.54502791020241-2.7391341646235783E7i");

    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "L_unit/2/6");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "3.7500E-9");
    model.component("comp1").mesh("mesh1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size1").set("hmax", "L_unit/2/n_wg/6");
    model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hmin", "1.1030E-9");
    model.component("comp1").mesh("mesh1").feature("size1").set("hminactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hcurveactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh1").feature("size1").set("hgradactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size2").set("hmax", "w_ring/30");
    model.component("comp1").mesh("mesh1").feature("size2").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hmin", "1.8750E-9");
    model.component("comp1").mesh("mesh1").feature("size2").set("hminactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hcurveactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hnarrowactive", true);
    model.component("comp1").mesh("mesh1").feature("size2").set("hgradactive", true);
    model.component("comp1").mesh("mesh1").run();

    model.study().create("std2");
    model.study("std2").create("param1", "Parametric");
    model.study("std2").create("bma", "BoundaryModeAnalysis");
    model.study("std2").create("bma1", "BoundaryModeAnalysis");
    model.study("std2").create("bma2", "BoundaryModeAnalysis");
    model.study("std2").create("bma3", "BoundaryModeAnalysis");
    model.study("std2").create("freq", "Frequency");
    model.study("std2").feature("bma")
         .set("activate", new String[]{"ewfd_TE", "off", "ewfd_TM", "on", "frame:spatial1", "on", "frame:material1", "on"});
    model.study("std2").feature("bma1")
         .set("activate", new String[]{"ewfd_TE", "off", "ewfd_TM", "on", "frame:spatial1", "on", "frame:material1", "on"});
    model.study("std2").feature("bma2")
         .set("activate", new String[]{"ewfd_TE", "off", "ewfd_TM", "on", "frame:spatial1", "on", "frame:material1", "on"});
    model.study("std2").feature("bma3")
         .set("activate", new String[]{"ewfd_TE", "off", "ewfd_TM", "on", "frame:spatial1", "on", "frame:material1", "on"});
    model.study("std2").feature("freq")
         .set("activate", new String[]{"ewfd_TE", "off", "ewfd_TM", "on", "frame:spatial1", "on", "frame:material1", "on"});

    model.sol().create("sol1");
    model.sol("sol1").study("std2");
    model.sol("sol1").attach("std2");
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
    model.sol().create("sol6");
    model.sol("sol6").study("std2");
    model.sol("sol6").label("Parametric Solutions 1");

    model.batch().create("p1", "Parametric");
    model.batch().create("p2", "Parametric");
    model.batch("p1").create("so1", "Solutionseq");
    model.batch("p2").create("so1", "Solutionseq");
    model.batch("p1").study("std2");
    model.batch("p2").study("std2");

    model.result().create("pg1", "PlotGroup2D");
    model.result().create("pg2", "PlotGroup1D");
    model.result().create("pg3", "PlotGroup2D");
    model.result().create("pg4", "PlotGroup2D");
    model.result("pg1").set("data", "dset6");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").feature("surf1").set("expr", "ewfd_TM.Ez");
    model.result("pg2").set("data", "dset6");
    model.result("pg2").create("glob1", "Global");
    model.result("pg3").set("data", "dset6");
    model.result("pg3").create("surf1", "Surface");
    model.result("pg3").feature("surf1").set("expr", "ewfd_TE.n_iso");
    model.result("pg4").set("data", "dset6");
    model.result("pg4").create("surf1", "Surface");
    model.result("pg4").feature("surf1").set("expr", "imag(ewfd_TM.Ez)");
    model.result().export().create("tbl1", "Table");

    model.study("std2").label("Study_TM");
    model.study("std2").feature("param1").label("Parametric Sweep");
    model.study("std2").feature("param1").set("sweeptype", "filled");
    model.study("std2").feature("param1").set("pname", new String[]{"f_sweep"});
    model.study("std2").feature("param1").set("plistarr", new String[]{"nu0+range(-0.002,0.001,0.002)"});
    model.study("std2").feature("param1").set("punit", new String[]{"1"});
    model.study("std2").feature("bma").set("modeFreq", "f_sweep*c_const/r_ex");
    model.study("std2").feature("bma").set("shift", "n_wg");
    model.study("std2").feature("bma1").set("PortName", "2");
    model.study("std2").feature("bma1").set("modeFreq", "f_sweep*c_const/r_ex");
    model.study("std2").feature("bma1").set("shift", "n_wg");
    model.study("std2").feature("bma2").set("PortName", "3");
    model.study("std2").feature("bma2").set("modeFreq", "f_sweep*c_const/r_ex");
    model.study("std2").feature("bma2").set("shift", "n_wg");
    model.study("std2").feature("bma3").set("PortName", "4");
    model.study("std2").feature("bma3").set("modeFreq", "f_sweep*c_const/r_ex");
    model.study("std2").feature("bma3").set("shift", "n_wg");
    model.study("std2").feature("freq").set("plist", "f_sweep*c_const/r_ex");

    model.sol("sol1").attach("std2");
    model.sol("sol1").feature("e1").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e1").set("neigs", 1);
    model.sol("sol1").feature("e1").set("shift", "n_wg");
    model.sol("sol1").feature("e1").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e1").feature("d1").label("Suggested Direct Solver (ewfd_TM)");
    model.sol("sol1").feature("st2").set("studystep", "bma1");
    model.sol("sol1").feature("v2").set("initmethod", "sol");
    model.sol("sol1").feature("v2").set("initsol", "sol1");
    model.sol("sol1").feature("v2").set("solnum", "auto");
    model.sol("sol1").feature("v2").set("notsolmethod", "sol");
    model.sol("sol1").feature("v2").set("notsol", "sol1");
    model.sol("sol1").feature("v2").set("notsolnum", "auto");
    model.sol("sol1").feature("e2").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e2").set("neigs", 1);
    model.sol("sol1").feature("e2").set("shift", "n_wg");
    model.sol("sol1").feature("e2").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e2").feature("d1").label("Suggested Direct Solver (ewfd_TM)");
    model.sol("sol1").feature("st3").set("studystep", "bma2");
    model.sol("sol1").feature("v3").set("initmethod", "sol");
    model.sol("sol1").feature("v3").set("initsol", "sol1");
    model.sol("sol1").feature("v3").set("solnum", "auto");
    model.sol("sol1").feature("v3").set("notsolmethod", "sol");
    model.sol("sol1").feature("v3").set("notsol", "sol1");
    model.sol("sol1").feature("v3").set("notsolnum", "auto");
    model.sol("sol1").feature("e3").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e3").set("neigs", 1);
    model.sol("sol1").feature("e3").set("shift", "n_wg");
    model.sol("sol1").feature("e3").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e3").feature("d1").label("Suggested Direct Solver (ewfd_TM)");
    model.sol("sol1").feature("st4").set("studystep", "bma3");
    model.sol("sol1").feature("v4").set("initmethod", "sol");
    model.sol("sol1").feature("v4").set("initsol", "sol1");
    model.sol("sol1").feature("v4").set("solnum", "auto");
    model.sol("sol1").feature("v4").set("notsolmethod", "sol");
    model.sol("sol1").feature("v4").set("notsol", "sol1");
    model.sol("sol1").feature("v4").set("notsolnum", "auto");
    model.sol("sol1").feature("e4").set("transform", "effective_mode_index");
    model.sol("sol1").feature("e4").set("neigs", 1);
    model.sol("sol1").feature("e4").set("shift", "n_wg");
    model.sol("sol1").feature("e4").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("e4").feature("d1").label("Suggested Direct Solver (ewfd_TM)");
    model.sol("sol1").feature("st5").set("studystep", "freq");
    model.sol("sol1").feature("v5").set("initmethod", "sol");
    model.sol("sol1").feature("v5").set("initsol", "sol1");
    model.sol("sol1").feature("v5").set("solnum", "auto");
    model.sol("sol1").feature("v5").set("notsolmethod", "sol");
    model.sol("sol1").feature("v5").set("notsol", "sol1");
    model.sol("sol1").feature("v5").set("notsolnum", "auto");
    model.sol("sol1").feature("v5").set("clistctrl", new String[]{"p1"});
    model.sol("sol1").feature("v5").set("cname", new String[]{"freq"});
    model.sol("sol1").feature("v5").set("clist", new String[]{"f_sweep*c_const/r_ex"});
    model.sol("sol1").feature("s1").set("control", "freq");
    model.sol("sol1").feature("s1").set("stol", 0.01);
    model.sol("sol1").feature("s1").set("probesel", "none");
    model.sol("sol1").feature("s1").feature("aDef").set("complexfun", true);
    model.sol("sol1").feature("s1").feature("p1").set("control", "freq");
    model.sol("sol1").feature("s1").feature("p1").set("pname", new String[]{"freq"});
    model.sol("sol1").feature("s1").feature("p1").set("plistarr", new String[]{"f_sweep*c_const/r_ex"});
    model.sol("sol1").feature("s1").feature("p1").set("punit", new String[]{"THz"});
    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");
    model.sol("sol1").feature("s1").feature("p1").set("preusesol", "auto");
    model.sol("sol1").feature("s1").feature("d1").label("Suggested Direct Solver (ewfd_TM)");
    model.sol("sol1").runAll();

    model.batch("p1").set("control", "param1");
    model.batch("p1").set("sweeptype", "filled");
    model.batch("p1").set("pname", new String[]{"f_sweep"});
    model.batch("p1").set("plistarr", new String[]{"nu0+range(-0.002,0.001,0.002)"});
    model.batch("p1").set("punit", new String[]{"1"});
    model.batch("p1").set("err", true);
    model.batch("p1").feature("so1")
         .set("param", new String[]{"\"f_sweep\",\"1.551169\"", "\"f_sweep\",\"1.552169\"", "\"f_sweep\",\"1.553169\""});
    model.batch("p1").run();
    model.batch("p2").set("control", "param1");
    model.batch("p2").set("sweeptype", "filled");
    model.batch("p2").set("pname", new String[]{"f_sweep"});
    model.batch("p2").set("plistarr", new String[]{"nu0+range(-0.002,0.001,0.002)"});
    model.batch("p2").set("punit", new String[]{"1"});
    model.batch("p2").set("err", true);
    model.batch("p2").feature("so1").set("seq", "sol1");
    model.batch("p2").feature("so1").set("psol", "sol6");
    model.batch("p2").feature("so1")
         .set("param", new String[]{"\"f_sweep\",\"1.550169\"", "\"f_sweep\",\"1.551169\"", "\"f_sweep\",\"1.552169\"", "\"f_sweep\",\"1.553169\"", "\"f_sweep\",\"1.554169\""});
    model.batch("p2").attach("std2");
    model.batch("p2").run();

    model.result("pg1").label("Electric Field (ewfd_TM)");
    model.result("pg1").set("looplevel", new int[]{1, 3});
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").feature("surf1").set("rangecoloractive", true);
    model.result("pg1").feature("surf1").set("rangecolormin", -188278.78832014967);
    model.result("pg1").feature("surf1").set("rangecolormax", 188492.2640868896);
    model.result("pg1").feature("surf1").set("colortable", "Wave");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("resolution", "normal");
    model.result("pg2").label("Reflectance, Transmittance, and Absorptance (ewfd_TM)");
    model.result("pg2").set("ylabel", "Reflectance, transmittance, and absorptance (1)");

    return model;
  }

  public static Model run2(Model model) {
    model.result("pg2").set("ylabelactive", true);
    model.result("pg2").feature("glob1")
         .set("expr", new String[]{"ewfd_TM.Rport_1", "ewfd_TM.Tport_2", "ewfd_TM.Tport_3", "ewfd_TM.Tport_4", "ewfd_TM.Ttotal", "ewfd_TM.RTtotal", "ewfd_TM.Atotal"});
    model.result("pg2").feature("glob1").set("unit", new String[]{"1", "1", "1", "1", "1", "1", "1"});
    model.result("pg2").feature("glob1")
         .set("descr", new String[]{"Reflectance, port 1", "Transmittance, port 2", "Transmittance, port 3", "Transmittance, port 4", "Total transmittance", "Total reflectance and transmittance", "Absorptance"});
    model.result("pg2").feature("glob1").set("titletype", "none");
    model.result("pg2").feature("glob1").set("xdatasolnumtype", "all");
    model.result("pg2").feature("glob1").set("xdata", "expr");
    model.result("pg2").feature("glob1").set("xdataexpr", "f_sweep");
    model.result("pg2").feature("glob1").set("xdataunit", "1");
    model.result("pg2").feature("glob1").set("xdatadescr", "");
    model.result("pg3").set("looplevel", new int[]{1, 2});
    model.result("pg3").feature("surf1").set("resolution", "normal");
    model.result("pg4").label("Electric Field (ewfd_TM) 1");
    model.result("pg4").set("looplevel", new int[]{1, 3});
    model.result("pg4").set("frametype", "spatial");
    model.result("pg4").feature("surf1").set("rangecoloractive", true);
    model.result("pg4").feature("surf1").set("rangecolormin", -188278.78832014967);
    model.result("pg4").feature("surf1").set("rangecolormax", 188492.2640868896);
    model.result("pg4").feature("surf1").set("colortable", "Wave");
    model.result("pg4").feature("surf1").set("smooth", "internal");
    model.result("pg4").feature("surf1").set("resolution", "normal");
    model.result().export("tbl1")
         .set("filename", "CouplingTriangles_double_triangle_around_peaks_v5.mphS_params.csv");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    run2(model);
  }

}
