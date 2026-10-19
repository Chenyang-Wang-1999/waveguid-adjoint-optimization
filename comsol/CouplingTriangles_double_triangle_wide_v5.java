/*
 * CouplingTriangles_double_triangle_wide_v5.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 19 2026, 20:32 by COMSOL 5.5.0.359. */
public class CouplingTriangles_double_triangle_wide_v5 {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/home/wangchenyang/654/research-data/FEM/adjoint_optimization/comsol");

    model.label("CouplingTriangles_double_triangle_wide_v5.mph");

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

    model.component("comp1").geom("geom1").create("pol1", "Polygon");
    model.component("comp1").geom("geom1").feature("pol1").set("type", "open");
    model.component("comp1").geom("geom1").feature("pol1").set("source", "file");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("filename", "F:\\ResearchWorks\\MyWork\\CouplingRings\\adjoint_optimize\\L_opt_0.5\\ex_line_pts.csv");
    model.component("comp1").geom("geom1").create("pol2", "Polygon");
    model.component("comp1").geom("geom1").feature("pol2").set("type", "open");
    model.component("comp1").geom("geom1").feature("pol2").set("source", "file");
    model.component("comp1").geom("geom1").feature("pol2")
         .set("filename", "F:\\ResearchWorks\\MyWork\\CouplingRings\\adjoint_optimize\\L_opt_0.5\\in_line_pts.csv");
    model.component("comp1").geom("geom1").create("mir2", "Mirror");
    model.component("comp1").geom("geom1").feature("mir2").set("keep", true);
    model.component("comp1").geom("geom1").feature("mir2").selection("input").set("pol1", "pol2");
    model.component("comp1").geom("geom1").create("mir3", "Mirror");
    model.component("comp1").geom("geom1").feature("mir3").set("keep", true);
    model.component("comp1").geom("geom1").feature("mir3").set("axis", new int[]{0, 1});
    model.component("comp1").geom("geom1").feature("mir3").selection("input").set("mir2", "pol1", "pol2");
    model.component("comp1").geom("geom1").create("sca1", "Scale");
    model.component("comp1").geom("geom1").feature("sca1").set("factor", "L_unit");
    model.component("comp1").geom("geom1").feature("sca1").selection("input").set("mir2", "mir3", "pol1", "pol2");
    model.component("comp1").geom("geom1").create("rot4", "Rotate");
    model.component("comp1").geom("geom1").feature("rot4").setIndex("rot", "90", 0);
    model.component("comp1").geom("geom1").feature("rot4").selection("input").set("sca1");
    model.component("comp1").geom("geom1").create("c1", "Circle");
    model.component("comp1").geom("geom1").feature("c1").set("pos", new String[]{"- a/2", "- a/2/sqrt(3)"});
    model.component("comp1").geom("geom1").feature("c1").set("r", "r_ex");
    model.component("comp1").geom("geom1").create("c2", "Circle");
    model.component("comp1").geom("geom1").feature("c2").set("pos", new String[]{"- a/2", "- a/2/sqrt(3)"});
    model.component("comp1").geom("geom1").feature("c2").set("r", "r_in");
    model.component("comp1").geom("geom1").create("mov1", "Move");
    model.component("comp1").geom("geom1").feature("mov1")
         .set("displx", "(l_triangle-d_ring_tri)/tan(pi/3)+rho_wg_cen/sin(pi/3)");
    model.component("comp1").geom("geom1").feature("mov1").selection("input").set("rot4");
    model.component("comp1").geom("geom1").create("rot5", "Rotate");
    model.component("comp1").geom("geom1").feature("rot5").set("keep", true);
    model.component("comp1").geom("geom1").feature("rot5").set("rot", "range(60,60,300)");
    model.component("comp1").geom("geom1").feature("rot5").selection("input").set("mov1");
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
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex2").set("rot5(12)", 100);
    model.component("comp1").geom("geom1").create("ls2", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls2").selection("vertex1").set("dif2(1)", 2);
    model.component("comp1").geom("geom1").feature("ls2").selection("vertex2").set("rot5(17)", 100);
    model.component("comp1").geom("geom1").create("ls3", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls3").selection("vertex1").set("dif2(1)", 5);
    model.component("comp1").geom("geom1").feature("ls3").selection("vertex2").set("rot5(26)", 100);
    model.component("comp1").geom("geom1").create("ls4", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls4").selection("vertex1").set("dif2(1)", 6);
    model.component("comp1").geom("geom1").feature("ls4").selection("vertex2").set("rot5(21)", 100);
    model.component("comp1").geom("geom1").create("ls11", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls11").selection("vertex1").set("rot5(7)", 100);
    model.component("comp1").geom("geom1").feature("ls11").selection("vertex2").set("rot5(36)", 100);
    model.component("comp1").geom("geom1").create("ls12", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls12").selection("vertex1").set("rot5(2)", 100);
    model.component("comp1").geom("geom1").feature("ls12").selection("vertex2").set("rot5(31)", 100);
    model.component("comp1").geom("geom1").create("rot6", "Rotate");
    model.component("comp1").geom("geom1").feature("rot6").set("keep", true);
    model.component("comp1").geom("geom1").feature("rot6").set("rot", "range(60,60,300)");
    model.component("comp1").geom("geom1").feature("rot6").selection("input")
         .set("ls1", "ls11", "ls12", "ls2", "ls3", "ls4");
    model.component("comp1").geom("geom1").create("c5", "Circle");
    model.component("comp1").geom("geom1").feature("c5").active(false);
    model.component("comp1").geom("geom1").feature("c5").set("r", "r_air");
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
    model.component("comp1").material("mat2").selection().set(2, 4, 5, 6, 12, 13, 16, 17, 18, 20, 21, 24);
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
         .set(1, 2, 3, 6, 9, 12, 13, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 449, 450, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 778, 779, 780, 781, 782, 783, 786, 787, 790, 791, 792, 793, 794, 795, 798, 799, 802, 803, 804, 805, 806, 807, 810, 811, 812, 813, 816, 817, 818, 819, 822, 823, 824, 825, 828, 829, 832, 833, 836, 837, 838, 839, 840, 841, 844, 845, 848, 849, 850, 851, 854, 855, 856, 857, 860, 861, 862, 863, 864, 865, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1014, 1015, 1018, 1019, 1022, 1023, 1024, 1027, 1030, 1031, 1032, 1033, 1038, 1039, 1040, 1041, 1046, 1047, 1048, 1049, 1054, 1055, 1058, 1059, 1064, 1065, 1068, 1069, 1074, 1075, 1078, 1079, 1082, 1083, 1088, 1089, 1092, 1093, 1096, 1097, 1098, 1099, 1104, 1105, 1106, 1107, 1112, 1113, 1114, 1115, 1118, 1119, 1122, 1123, 1126, 1127, 1130, 1131, 1132, 1133, 1138, 1139, 1140, 1141, 1142, 1143, 1146, 1147, 1150, 1151, 1152, 1153, 1156, 1157, 1158, 1159, 1160, 1161, 1166, 1167, 1168, 1169, 1170, 1171, 1174, 1175, 1178, 1179, 1180, 1181, 1182, 1183, 1188, 1189, 1190, 1191, 1192, 1193, 1196, 1197, 1200, 1201, 1202, 1203, 1206, 1207, 1208, 1209, 1212, 1213, 1216, 1217, 1218, 1219, 1220, 1221, 1224, 1225, 1226, 1227, 1230, 1231, 1234, 1235, 1236, 1237, 1240, 1241, 1244, 1245, 1246, 1247, 1250, 1251, 1252, 1253, 1258, 1259, 1260, 1261, 1262, 1263, 1266, 1267, 1270, 1271, 1272, 1273, 1276, 1277, 1280, 1281, 1282, 1283, 1284, 1285, 1288, 1289, 1292, 1293, 1294, 1295, 1296, 1297, 1300, 1301, 1302, 1303, 1306, 1307, 1308, 1309, 1312, 1313, 1314, 1315, 1316, 1317, 1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1330, 1331, 1332, 1333, 1334, 1335, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1355, 1358, 1359, 1360, 1361, 1362, 1363, 1366, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375, 1376, 1377, 1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1412, 1413, 1414, 1415, 1416, 1417, 1420, 1421, 1422, 1423, 1426, 1427, 1428, 1429, 1430, 1431, 1434, 1435, 1438, 1439, 1440, 1441, 1444, 1445, 1446, 1447, 1456, 1457, 1460, 1461, 1472, 1473, 1476, 1477, 1480, 1481, 1484, 1485, 1486, 1487, 1490, 1491, 1494, 1495, 1496, 1497, 1500, 1501, 1502, 1503, 1506, 1507, 1508, 1509, 1512, 1513, 1514, 1515, 1516, 1517, 1521, 1522, 1523, 1524, 1525, 1526, 1527, 1528, 1529, 1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 1538, 1539, 1540, 1541, 1542, 1543, 1544, 1545, 1546, 1547, 1548, 1549, 1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558, 1559, 1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 1570, 1571, 1572, 1573, 1574, 1575, 1576, 1577, 1578, 1579, 1580, 1581, 1582, 1583, 1584, 1585, 1586, 1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595, 1596, 1597, 1598, 1599, 1600, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608, 1609, 1610, 1611, 1612, 1613, 1614, 1615, 1616, 1617, 1618, 1619, 1620, 1621, 1622, 1623, 1624, 1625, 1626, 1627, 1628, 1629, 1630, 1631, 1632, 1633, 1634, 1635, 1636, 1637, 1638, 1639, 1640, 1641, 1642, 1643, 1644, 1645, 1646, 1647, 1648, 1649, 1650, 1651, 1652, 1653, 1654, 1655, 1656, 1657, 1658, 1659, 1660, 1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668, 1669, 1670, 1671, 1672, 1673, 1674, 1675, 1676, 1677, 1678, 1679, 1680, 1681, 1682, 1683, 1684, 1685, 1686, 1687, 1688, 1689, 1690, 1691, 1692, 1693, 1694, 1695, 1696, 1697, 1698, 1699, 1700, 1701, 1702, 1703, 1704, 1705, 1706, 1708, 1709, 1710, 1711, 1712, 1713, 1714, 1715, 1716, 1717, 1720, 1721, 1722, 1723, 1724, 1725, 1726, 1727, 1730, 1731, 1732, 1733, 1734, 1735, 1736, 1737, 1740, 1741, 1742, 1743, 1744, 1745, 1746, 1747, 1748, 1749, 1752, 1753, 1754, 1755, 1756, 1757, 1758, 1759, 1760, 1761, 1764, 1765, 1766, 1767, 1768, 1769, 1770, 1771, 1772, 1773, 1774, 1775, 1776, 1777, 1778, 1779, 1780, 1781, 1782, 1783, 1784, 1785, 1786, 1787, 1788, 1789, 1792, 1793, 1794, 1795, 1796, 1797, 1798, 1799, 1802, 1803, 1804, 1805, 1806, 1807, 1808, 1809, 1810, 1811, 1812, 1813, 1816, 1817, 1818, 1819, 1820, 1821, 1822, 1823, 1824, 1825, 1828, 1829, 1830, 1831, 1832, 1833, 1834, 1835, 1836, 1837, 1840, 1841, 1842, 1843, 1844, 1845, 1846, 1847, 1850, 1851, 1852, 1853, 1854, 1855, 1856, 1857, 1860, 1861, 1862, 1863, 1864, 1865, 1866, 1867, 1870, 1871, 1872, 1873, 1874, 1875, 1878, 1879, 1880, 1881, 1882, 1883, 1886, 1887, 1888, 1889, 1890, 1891, 1892, 1893, 1896, 1897, 1898, 1899, 1900, 1901, 1904, 1905, 1906, 1907, 1908, 1909, 1912, 1913, 1914, 1915, 1916, 1917, 1920, 1921, 1922, 1923, 1924, 1925, 1926, 1927, 1930, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2053, 2056, 2057, 2058, 2059, 2060, 2061, 2062, 2063, 2066, 2067, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2078, 2079, 2080, 2081, 2082, 2083, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2093, 2096, 2097, 2098, 2099, 2100, 2101, 2102, 2103, 2106, 2107, 2110, 2111, 2112, 2113, 2114, 2115, 2118, 2119, 2120, 2121, 2124, 2125, 2126, 2127, 2128, 2129, 2132, 2133, 2134, 2135, 2138, 2139, 2140, 2141, 2144, 2145, 2148, 2149, 2150, 2151, 2152, 2153, 2156, 2157, 2158, 2159, 2160, 2161, 2164, 2165, 2166, 2167, 2168, 2169, 2170, 2171, 2174, 2175, 2176, 2177, 2178, 2179, 2180, 2181, 2184, 2185, 2186, 2187, 2188, 2189, 2192, 2193, 2194, 2195, 2198, 2199, 2200, 2201, 2202, 2203, 2206, 2207, 2208, 2209, 2210, 2213, 2214, 2215, 2216, 2217, 2220, 2221, 2222, 2223, 2224, 2225, 2228, 2229, 2232, 2233, 2236, 2237, 2240, 2241, 2244, 2245, 2248, 2249, 2252, 2253, 2258, 2259, 2262, 2263, 2266, 2267, 2270, 2271, 2274, 2275, 2278, 2279, 2282, 2283, 2286, 2287, 2290, 2291, 2294, 2295, 2300, 2301, 2304, 2305, 2308, 2309, 2312, 2313, 2316, 2317, 2320, 2321, 2324, 2325, 2328, 2329, 2332, 2333, 2338, 2339, 2342, 2343, 2346, 2347, 2350, 2351, 2354, 2355, 2358, 2359, 2362, 2363, 2366, 2367, 2370, 2371, 2376, 2377, 2380, 2381, 2384, 2385, 2388, 2389, 2392, 2393, 2396, 2397, 2400, 2401, 2404, 2405, 2408, 2409, 2413, 2414, 2417, 2418, 2422, 2423, 2426, 2427, 2430, 2431, 2434, 2435, 2438, 2439, 2442, 2443, 2446, 2447, 2450, 2451, 2454, 2455, 2460, 2461, 2464, 2465, 2468, 2469, 2472, 2473, 2476, 2477, 2480, 2481, 2484, 2485, 2488, 2489, 2492, 2493, 2498, 2499, 2502, 2503, 2506, 2507, 2510, 2511, 2514, 2515, 2518, 2519, 2522, 2523, 2526, 2527, 2530, 2531, 2536, 2537, 2540, 2541, 2544, 2545, 2548, 2549, 2552, 2553, 2556, 2557, 2560, 2561, 2564, 2565, 2568, 2569, 2572, 2573, 2578, 2579, 2580, 2582, 2586, 2587, 2590, 2591, 2594, 2595, 2598, 2599, 2602, 2603, 2606, 2607, 2610, 2611, 2612, 2613, 2614, 2615, 2618, 2619, 2620, 2621, 2622, 2625, 2626, 2627, 2628, 2629, 2632, 2633, 2634, 2635, 2636, 2637, 2640, 2641, 2642, 2643, 2646, 2647, 2648, 2649, 2650, 2651, 2654, 2655, 2656, 2657, 2658, 2659, 2660, 2661, 2664, 2665, 2666, 2667, 2668, 2669, 2670, 2671, 2674, 2675, 2676, 2677, 2678, 2679, 2682, 2683, 2684, 2685, 2686, 2687, 2688, 2689, 2690, 2691, 2694, 2695, 2696, 2697, 2700, 2701, 2702, 2703, 2706, 2707, 2708, 2709, 2710, 2711, 2714, 2715, 2716, 2717, 2720, 2721, 2722, 2723, 2724, 2725, 2728, 2729, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 2739, 2742, 2743, 2744, 2745, 2746, 2747, 2748, 2749, 2752, 2753, 2754, 2755, 2756, 2757, 2760, 2761, 2762, 2763, 2764, 2765, 2766, 2767, 2768, 2769, 2772, 2773, 2774, 2775, 2776, 2777, 2778, 2779, 2782, 2783, 2784, 2785, 2786, 2787, 2788, 2789, 2790, 2791, 2794, 2795, 2796, 2797, 2798, 2799, 2800, 2801, 2802, 2803, 2804, 2805, 2808, 2809, 2810, 2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818, 2819, 2820, 2821, 2822, 2823, 2826, 2827, 2828, 2829, 2830, 2831, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2839, 2840, 2841, 2842, 2843, 2846, 2847, 2848, 2849, 2850, 2851, 2852, 2853, 2854, 2855, 2856, 2857, 2860, 2861, 2862, 2863, 2864, 2865, 2866, 2867, 2868, 2869, 2872, 2873, 2874, 2875, 2876, 2877, 2878, 2879, 2880, 2881, 2882, 2883, 2886, 2887, 2888, 2889, 2890, 2891, 2892, 2893, 2896, 2897, 2898, 2899, 2900, 2901, 2902, 2903, 2904, 2905, 2908, 2909, 2910, 2911, 2912, 2913, 2914, 2915, 2918, 2919, 2920, 2921, 2922, 2923, 2926, 2927, 2928, 2929, 2930, 2931, 2934, 2935, 2936, 2937, 2938, 2939, 2942, 2943, 2944, 2945, 2946, 2947, 2948, 2949, 2952, 2953, 2954, 2955, 2956, 2957, 2960, 2961, 2962, 2963, 2964, 2965, 2968, 2969, 2970, 2971, 2972, 2973, 2974, 2975, 2978, 2979, 2980, 2981, 2982, 2983, 2984, 2985, 2988, 2989, 2990, 2991, 2992, 2993, 2994, 2995, 2998, 2999, 3000, 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3010, 3011, 3012, 3013, 3014, 3015, 3016, 3017, 3018, 3019, 3022, 3023, 3024, 3025, 3026, 3027, 3028, 3029, 3030, 3031, 3032, 3033, 3036, 3037, 3038, 3039, 3040, 3041, 3042, 3043, 3046, 3047, 3048, 3049, 3050, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3059, 3060, 3061, 3062, 3063, 3064, 3065, 3066, 3067, 3068, 3069, 3070, 3071, 3074, 3075, 3076, 3077, 3078, 3079, 3080, 3081, 3082, 3083, 3086, 3087, 3088, 3089, 3090, 3091, 3092, 3093, 3094, 3095, 3098, 3099, 3100, 3101, 3102, 3103, 3104, 3105, 3108, 3109, 3110, 3111, 3112, 3113, 3114, 3115, 3118, 3119, 3121, 3122, 3123, 3124, 3125, 3126, 3127, 3128, 3129, 3130, 3131, 3132, 3133, 3134, 3135, 3136, 3137, 3138, 3139, 3140, 3141, 3142, 3143, 3144, 3145, 3146, 3147, 3148, 3149, 3150, 3151, 3152, 3153, 3154, 3155, 3156, 3157, 3158, 3159, 3160, 3161, 3162, 3163, 3164, 3165, 3166, 3167, 3168, 3169, 3170, 3171, 3172, 3173, 3174, 3175, 3176, 3177, 3178, 3179, 3180, 3181, 3182, 3183, 3184, 3185, 3186, 3187, 3188, 3189, 3190, 3191, 3192, 3193, 3194, 3195, 3196, 3197, 3198, 3199, 3200, 3201, 3202, 3203, 3204, 3205, 3206, 3207, 3208, 3209, 3210, 3211, 3212, 3213, 3214, 3215, 3216, 3217, 3218, 3219, 3220, 3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229, 3230, 3231, 3232, 3233, 3234, 3235, 3236, 3237, 3238, 3239, 3240, 3241, 3242, 3243, 3244, 3245, 3246, 3247, 3248, 3249, 3250, 3251, 3252, 3253, 3254, 3255, 3256, 3257, 3258, 3259, 3260, 3261, 3262, 3263, 3264, 3265, 3266, 3267, 3268, 3269, 3270, 3271, 3272, 3273, 3274, 3275, 3276, 3277, 3278, 3279, 3280, 3281, 3282, 3283, 3284, 3285, 3286, 3287, 3288, 3289, 3290, 3291, 3292, 3293, 3294, 3295, 3296, 3297, 3298, 3299, 3300, 3301, 3302, 3303, 3304, 3305, 3309, 3310, 3311, 3312, 3313, 3318, 3319, 3320, 3321, 3322, 3323, 3326, 3327, 3328, 3329, 3332, 3333, 3334, 3335, 3338, 3339, 3340, 3341, 3344, 3345, 3348, 3349, 3350, 3351, 3354, 3355, 3358, 3359, 3362, 3363, 3374, 3375, 3378, 3379, 3388, 3389, 3390, 3391, 3394, 3395, 3396, 3397, 3400, 3401, 3404, 3405, 3406, 3407, 3408, 3409, 3412, 3413, 3414, 3415, 3418, 3419, 3420, 3421, 3422, 3423, 3426, 3427, 3428, 3429, 3430, 3431, 3432, 3433, 3436, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3446, 3447, 3448, 3449, 3450, 3451, 3452, 3453, 3454, 3455, 3458, 3459, 3460, 3461, 3462, 3463, 3464, 3465, 3466, 3467, 3468, 3469, 3472, 3473, 3474, 3475, 3476, 3477, 3480, 3481, 3482, 3483, 3484, 3485, 3486, 3487, 3490, 3491, 3492, 3493, 3494, 3495, 3496, 3497, 3500, 3501, 3502, 3503, 3504, 3505, 3508, 3509, 3510, 3511, 3512, 3513, 3514, 3515, 3518, 3519, 3520, 3521, 3522, 3523, 3526, 3527, 3528, 3529, 3532, 3533, 3534, 3535, 3538, 3539, 3540, 3541, 3542, 3543, 3546, 3547, 3550, 3551, 3552, 3553, 3554, 3555, 3558, 3559, 3562, 3563, 3564, 3565, 3568, 3569, 3572, 3573, 3574, 3575, 3576, 3577, 3582, 3583, 3584, 3585, 3588, 3589, 3590, 3591, 3594, 3595, 3598, 3599, 3600, 3601, 3604, 3605, 3608, 3609, 3610, 3611, 3614, 3615, 3616, 3617, 3618, 3619, 3622, 3623, 3626, 3627, 3628, 3629, 3632, 3633, 3634, 3635, 3638, 3639, 3642, 3643, 3644, 3645, 3646, 3647, 3652, 3653, 3654, 3655, 3656, 3657, 3660, 3661, 3664, 3665, 3666, 3667, 3668, 3669, 3674, 3675, 3676, 3677, 3678, 3679, 3682, 3683, 3684, 3685, 3688, 3689, 3692, 3693, 3694, 3695, 3696, 3697, 3702, 3703, 3704, 3705, 3708, 3709, 3712, 3713, 3716, 3717, 3720, 3721, 3722, 3723, 3728, 3729, 3730, 3731, 3736, 3737, 3738, 3739, 3742, 3743, 3746, 3747, 3752, 3753, 3756, 3757, 3760, 3761, 3766, 3767, 3770, 3771, 3776, 3777, 3780, 3781, 3786, 3787, 3788, 3789, 3794, 3795, 3796, 3797, 3802, 3803, 3804, 3805, 3808, 3811, 3812, 3813, 3816, 3817, 3820, 3821, 3824, 3825, 3826, 3827, 3828, 3829, 3830, 3831, 3832, 3833, 3834, 3835, 3836, 3837, 3838, 3839, 3840, 3841, 3842, 3843, 3844, 3845, 3846, 3847, 3848, 3849, 3850, 3851, 3852, 3853, 3854, 3855, 3856, 3857, 3858, 3859, 3860, 3861, 3862, 3863, 3864, 3865, 3866, 3867, 3868, 3869, 3870, 3871, 3872, 3873, 3874, 3875, 3876, 3877, 3878, 3879, 3880, 3881, 3882, 3883, 3884, 3885, 3886, 3887, 3888, 3889, 3890, 3891, 3892, 3893, 3894, 3895, 3896, 3897, 3898, 3899, 3900, 3901, 3902, 3903, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 3915, 3916, 3917, 3918, 3919, 3920, 3921, 3922, 3923, 3924, 3925, 3926, 3927, 3928, 3929, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 3939, 3940, 3941, 3942, 3943, 3944, 3945, 3946, 3947, 3948, 3949, 3950, 3951, 3952, 3953, 3954, 3955, 3956, 3957, 3958, 3959, 3960, 3961, 3962, 3963, 3964, 3965, 3966, 3967, 3970, 3971, 3972, 3973, 3974, 3975, 3978, 3979, 3980, 3981, 3984, 3985, 3986, 3987, 3990, 3991, 3994, 3995, 3996, 3997, 3998, 3999, 4002, 4003, 4006, 4007, 4010, 4011, 4012, 4013, 4016, 4017, 4018, 4019, 4022, 4023, 4024, 4025, 4028, 4029, 4030, 4031, 4032, 4033, 4036, 4037, 4040, 4041, 4042, 4043, 4044, 4045, 4048, 4049, 4052, 4053, 4054, 4055, 4056, 4057, 4060, 4061, 4062, 4063, 4065, 4066, 4070, 4071, 4072, 4073, 4074, 4075, 4078, 4079, 4080, 4081, 4084, 4085, 4088, 4089, 4090, 4091, 4094, 4095, 4096, 4097, 4100, 4101, 4104, 4105, 4108, 4109, 4110, 4111, 4114, 4115, 4118, 4119, 4120, 4121, 4124, 4125, 4128, 4129, 4132, 4133, 4134, 4135, 4138, 4139, 4142, 4143, 4144, 4145, 4146, 4147, 4152, 4153, 4156, 4157, 4158, 4159, 4160, 4161, 4164, 4165, 4168, 4169, 4170, 4171, 4172, 4173, 4176, 4177, 4180, 4181, 4182, 4183, 4184, 4185, 4188, 4189, 4190, 4191, 4194, 4195, 4196, 4197, 4198, 4199, 4202, 4203, 4204, 4205, 4208, 4209, 4210, 4211, 4212, 4213, 4216, 4217, 4218, 4219, 4220, 4221, 4224, 4225, 4226, 4227, 4230, 4231, 4232, 4233, 4236, 4237, 4238, 4239, 4240, 4241, 4244, 4245, 4248, 4249, 4250, 4251, 4252, 4253, 4256, 4257, 4258, 4259, 4262, 4263, 4266, 4267, 4272, 4273, 4278, 4279, 4282, 4283, 4284, 4285, 4286, 4287, 4288, 4289, 4292, 4293, 4294, 4295, 4296, 4297, 4300, 4301, 4310, 4311, 4318, 4319, 4322, 4323, 4334, 4335, 4340, 4341, 4344, 4345, 4346, 4347, 4348, 4349, 4352, 4353, 4354, 4355, 4356, 4357, 4358, 4359, 4360, 4361, 4366, 4367, 4454, 4455, 4458, 4459, 4464, 4465, 4466, 4467, 4468, 4469, 4470, 4471, 4472, 4473, 4474, 4475, 4476, 4477, 4480, 4481, 4482, 4483, 4484, 4485, 4488, 4489, 4498, 4499, 4510, 4511, 4514, 4515, 4522, 4523, 4528, 4529, 4532, 4533, 4534, 4535, 4536, 4537, 4540, 4541, 4542, 4543, 4546, 4547, 4552, 4553, 4558, 4559, 4562, 4563, 4568, 4569, 4572, 4573, 4576, 4577, 4580, 4581, 4584, 4585, 4588, 4589, 4592, 4593, 4594, 4595, 4598, 4599, 4602, 4603, 4604, 4605, 4608, 4609, 4612, 4613, 4614, 4615, 4618, 4619, 4622, 4623, 4626, 4627, 4630, 4631, 4634, 4635, 4640, 4641, 4644, 4645, 4650, 4651, 4656, 4657, 4660, 4661, 4666, 4667, 4672, 4673, 4676, 4677, 4682, 4683, 4686, 4687, 4691, 4692, 4694, 4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4765, 4766, 4767, 4768, 4769, 4770, 4771, 4772, 4773, 4774, 4775, 4776, 4777, 4778, 4779, 4780, 4781, 4784, 4785, 4786, 4787, 4788, 4789, 4790, 4791, 4792, 4793, 4794, 4795, 4796, 4797, 4798, 4799, 4800, 4801, 4802, 4803, 4804, 4805, 4806, 4807, 4808, 4809, 4810, 4811, 4812, 4813, 4816, 4819, 4822, 4823, 4824, 4825, 4826, 4829, 4830, 4833, 4834, 4835, 4836, 4841, 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4851, 4853, 4854, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872);

    return model;
  }

  public static Model run2(Model model) {
    model.component("comp1").physics("ewfd_TM").create("imp1", "Impedance", 1);
    model.component("comp1").physics("ewfd_TM").feature("imp1").selection()
         .set(1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 4814, 4815, 4817, 4818, 4820, 4821, 4822, 4823, 4824);
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
         .set(2, 4, 5, 6, 9, 12, 13, 16, 17, 18, 20, 24);
    model.component("comp1").mesh("mesh1").feature("size2").selection().geom("geom1", 1);
    model.component("comp1").mesh("mesh1").feature("size2").selection()
         .set(4, 5, 7, 8, 10, 11, 4814, 4815, 4817, 4818, 4820, 4821);
    model.component("comp1").mesh("mesh1").feature("ftri3").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftri3").selection().set(3);

    model.result().table("tbl1").comments("\u5168\u5c40\u8ba1\u7b97 1");

    model.component("comp1").view("view1").set("showmaterial", true);
    model.component("comp1").view("view1").axis().set("xmin", -4.798366717295721E-6);
    model.component("comp1").view("view1").axis().set("xmax", 4.76515015179757E-6);
    model.component("comp1").view("view1").axis().set("ymin", -4.200672265142202E-6);
    model.component("comp1").view("view1").axis().set("ymax", 4.197449925413821E-6);

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
         .set("LambdaPort", "-13.54502791094043-2.739134164623559E7i");
    model.component("comp1").physics("ewfd_TM").feature("port2").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port2").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port2")
         .set("LambdaPort", "-13.545027905780222-2.739134164623563E7i");
    model.component("comp1").physics("ewfd_TM").feature("port3").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port3").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port3")
         .set("LambdaPort", "-13.545027908728272-2.7391341646235757E7i");
    model.component("comp1").physics("ewfd_TM").feature("port4").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd_TM").feature("port4").set("PortFreq", "f_sweep*c_const/r_ex");
    model.component("comp1").physics("ewfd_TM").feature("port4")
         .set("LambdaPort", "-13.545027910939522-2.739134164623579E7i");

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
    model.result("pg1").set("data", "dset6");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").feature("surf1").set("expr", "ewfd_TM.normE");
    model.result("pg2").set("data", "dset6");
    model.result("pg2").create("glob1", "Global");
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
    model.result("pg1").set("titletype", "manual");
    model.result("pg1").set("title", "Electric field |E|, freq = eval(f_sweep)");
    model.result("pg1").set("allowevalintitle", true);
    model.result("pg1").set("edgecolor", "custom");
    model.result("pg1")
         .set("customedgecolor", new double[]{0.1568627506494522, 0.1568627506494522, 0.4627451002597809});
    model.result("pg1").feature("surf1").set("rangecoloractive", true);
    model.result("pg1").feature("surf1").set("rangecolormax", 50000);
    model.result("pg1").feature("surf1").set("colortable", "Thermal");
    model.result("pg1").feature("surf1").set("smooth", "internal");
    model.result("pg1").feature("surf1").set("resolution", "normal");
    model.result("pg2").label("Reflectance, Transmittance, and Absorptance (ewfd_TM)");
    model.result("pg2").set("ylabel", "Reflectance, transmittance, and absorptance (1)");
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
    model.result().export("tbl1")
         .set("filename", "CouplingTriangles_double_triangle_around_peaks_v5.mphS_params.csv");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    run2(model);
  }

}
