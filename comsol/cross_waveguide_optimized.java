/* ============================================================================
 * cross_waveguide_optimized.java
 *
 *   ⚠ 由 python/export_comsol.py 从 design.h5 拼字符串生成，不是 COMSOL 的
 *     “Save as Model File for Java” 导出件（那份带一行 Model exported on ...
 *     by COMSOL x.y 的戳，像 cross_waveguide_test.java）。语句照抄了参考文件
 *     能抄的部分；其余（CubicBezier / Box 选择 / named 选择 / study().run）
 *     是按 COMSOL 文档和 API 约定推的，没编译运行过。换设计请重新生成。
 *
 * 与 cross_waveguide_test.java 的关系：物理、网格、研究、结果那一套照搬
 * （ewfd + outofplane；4 个 Numeric 端口 + 4 步 Boundary Mode Analysis；
 * 非端口外边界用 Scattering 吸收；材料用 RefractiveIndex 填 n = sqrt(eps)），
 * 只把几何从“两条直矩形交叉”换成优化后的形状：
 *
 *   正六边形外框                    1 个 Polygon
 *   优化壁 in/ex x 4 象限 x 11 段  88 个 CubicBezier
 *   固定直壁段                      8 个 LineSegment
 *
 * 控制点就是求解器里 Bezier_interp 用的控制多边形，逐位相同，不做插值。
 *
 * 单位：几何数字的单位是模型长度单位（默认 m）。L_unit 是“原始无量纲坐标里
 * 1 对应多少米”，只出现在频率表达式里；要换物理尺度用
 *     python3 ../python/export_comsol.py --length-scale <s>
 * 重新生成（坐标乘 s、L_unit 同时变成 s[m]），这样几何和频率一起缩放。
 *
 * 选择：端口 / 散射边界 / 材料域 / 尖端全部按坐标的 Box 选择，不写死实体编号
 * （编号随几何变化，坐标不会）。哪个框没框住，点开对应的 Selection 就看得见。
 *
 * 端口编号 -> 我们的 S 参数（ewfd.S<k>1 = S_params[k-1]）：
 *   port1 input     中点 (+0.6000, +0.3464)  <-> S_params[0]
 *   port2 output1   中点 (-0.6000, -0.3464)  <-> S_params[1]
 *   port3 output21  中点 (+0.6000, -0.3464)  <-> S_params[2]
 *   port4 output22  中点 (-0.6000, +0.3464)  <-> S_params[3]
 * ============================================================================ */

import com.comsol.model.*;
import com.comsol.model.util.*;

public class cross_waveguide_optimized {

  public static Model run() {
    Model model = ModelUtil.create("Model");
    model.label("cross_waveguide_optimized.mph");

    // ---------------- 参数 ----------------
    model.param().set("L_unit", "1.0[m]");
    model.param().set("T_unit", "L_unit/c_const");
    model.param().set("P_unit", "1[W]");
    model.param().set("E_unit", "P_unit*T_unit");
    model.param().set("S_unit", "P_unit/L_unit^2");
    model.param().set("f", "1.55");
    model.param().set("epsilon_wg", "11.559999999999999");
    model.param().set("epsilon_amb", "1.0");

    // ---------------- 元件与几何 ----------------
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 2);

    // 外框：正六边形，顶点在半径 0.8 处、角度 60k 度
    model.component("comp1").geom("geom1").create("pol1", "Polygon");
    model.component("comp1").geom("geom1").feature("pol1").set("x", "0.8*cos(2*pi/6*range(0,1,6))");
    model.component("comp1").geom("geom1").feature("pol1").set("y", "0.8*sin(2*pi/6*range(0,1,6))");

    // 优化壁 in：4 个象限 x 11 段三次 Bézier（控制点 = 控制多边形）
    model.component("comp1").geom("geom1").create("cb1", "CubicBezier");   // in / input / 第 1 段
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.1932050807568877, 0, 0);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.06535898384862243, 1, 0);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.18948604127893776, 0, 1);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.06321179540490118, 1, 1);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.18540090743613008, 0, 2);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.06730831640982436, 1, 2);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.18072631051327112, 0, 3);
    model.component("comp1").geom("geom1").feature("cb1").setIndex("p", 0.06856138121147924, 1, 3);
    model.component("comp1").geom("geom1").create("cb2", "CubicBezier");   // in / input / 第 2 段
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.18072631051327112, 0, 0);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.06856138121147924, 1, 0);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.17605171359041216, 0, 1);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.06981444601313412, 1, 1);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.17078765358750186, 0, 2);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.0682240546115207, 1, 2);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.16629830156225683, 0, 3);
    model.component("comp1").geom("geom1").feature("cb2").setIndex("p", 0.06643557019288378, 1, 3);
    model.component("comp1").geom("geom1").create("cb3", "CubicBezier");   // in / input / 第 3 段
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.16629830156225683, 0, 0);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.06643557019288378, 1, 0);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.1618089495370118, 0, 1);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.06464708577424685, 1, 1);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.15809430548943199, 0, 2);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.06266050833858641, 1, 2);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.15443897163953926, 0, 3);
    model.component("comp1").geom("geom1").feature("cb3").setIndex("p", 0.06027673395740225, 1, 3);
    model.component("comp1").geom("geom1").create("cb4", "CubicBezier");   // in / input / 第 4 段
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.15443897163953926, 0, 0);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.06027673395740225, 1, 0);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.15078363778964654, 0, 1);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.05789295957621809, 1, 1);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.14718761413744094, 0, 2);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.05511198824951022, 1, 2);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.14409551518170757, 0, 3);
    model.component("comp1").geom("geom1").feature("cb4").setIndex("p", 0.052003453456223946, 1, 3);
    model.component("comp1").geom("geom1").create("cb5", "CubicBezier");   // in / input / 第 5 段
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.14409551518170757, 0, 0);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.052003453456223946, 1, 0);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.1410034162259742, 0, 1);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.04889491866293767, 1, 1);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.1384152419667131, 0, 2);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.04545882040307301, 1, 2);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.13596292843099772, 0, 3);
    model.component("comp1").geom("geom1").feature("cb5").setIndex("p", 0.042165909613575324, 1, 3);
    model.component("comp1").geom("geom1").create("cb6", "CubicBezier");   // in / input / 第 6 段
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.13596292843099772, 0, 0);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.042165909613575324, 1, 0);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.13351061489528235, 0, 1);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.03887299882407764, 1, 1);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.13119416208311274, 0, 2);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.03572327550494694, 1, 2);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.1287888229871091, 0, 3);
    model.component("comp1").geom("geom1").feature("cb6").setIndex("p", 0.03261530883032029, 1, 3);
    model.component("comp1").geom("geom1").create("cb7", "CubicBezier");   // in / input / 第 7 段
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.1287888229871091, 0, 0);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.03261530883032029, 1, 0);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.12638348389110543, 0, 1);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.029507342155693634, 1, 1);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.1238892585112677, 0, 2);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.026441132125571032, 1, 2);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.1211609012851205, 0, 3);
    model.component("comp1").geom("geom1").feature("cb7").setIndex("p", 0.02360373316033877, 1, 3);
    model.component("comp1").geom("geom1").create("cb8", "CubicBezier");   // in / input / 第 8 段
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.1211609012851205, 0, 0);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.02360373316033877, 1, 0);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.1184325440589733, 0, 1);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.020766334195106505, 1, 1);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.11547005498651665, 0, 2);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.018157746294764583, 1, 2);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.11208899594265007, 0, 3);
    model.component("comp1").geom("geom1").feature("cb8").setIndex("p", 0.015780857406382565, 1, 3);
    model.component("comp1").geom("geom1").create("cb9", "CubicBezier");   // in / input / 第 9 段
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.11208899594265007, 0, 0);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.015780857406382565, 1, 0);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.10870793689878348, 0, 1);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.013403968518000547, 1, 1);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.10490830788350697, 0, 2);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.011258778641578435, 1, 2);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.10125063709372815, 0, 3);
    model.component("comp1").geom("geom1").feature("cb9").setIndex("p", 0.010513186603964487, 1, 3);
    model.component("comp1").geom("geom1").create("cb10", "CubicBezier");   // in / input / 第 10 段
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.10125063709372815, 0, 0);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.010513186603964487, 1, 0);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.09759296630394934, 0, 1);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.009767594566350539, 1, 1);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.09407725373966822, 0, 2);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.010421600367544756, 1, 2);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.08936694442780331, 0, 3);
    model.component("comp1").geom("geom1").feature("cb10").setIndex("p", 0.009145941898970826, 1, 3);
    model.component("comp1").geom("geom1").create("cb11", "CubicBezier");   // in / input / 第 11 段
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.08936694442780331, 0, 0);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.009145941898970826, 1, 0);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.0846566351159384, 0, 1);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.007870283430396896, 1, 1);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.07875172905648971, 0, 2);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.004664960692054816, 1, 2);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.07877388587848679, 0, 3);
    model.component("comp1").geom("geom1").feature("cb11").setIndex("p", 0.0, 1, 3);
    model.component("comp1").geom("geom1").create("cb12", "CubicBezier");   // in / output22 / 第 1 段
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", -0.1932050807568877, 0, 0);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", 0.06535898384862243, 1, 0);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", -0.18948604127893776, 0, 1);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", 0.06321179540490118, 1, 1);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", -0.18540090743613008, 0, 2);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", 0.06730831640982436, 1, 2);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", -0.18072631051327112, 0, 3);
    model.component("comp1").geom("geom1").feature("cb12").setIndex("p", 0.06856138121147924, 1, 3);
    model.component("comp1").geom("geom1").create("cb13", "CubicBezier");   // in / output22 / 第 2 段
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", -0.18072631051327112, 0, 0);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", 0.06856138121147924, 1, 0);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", -0.17605171359041216, 0, 1);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", 0.06981444601313412, 1, 1);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", -0.17078765358750186, 0, 2);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", 0.0682240546115207, 1, 2);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", -0.16629830156225683, 0, 3);
    model.component("comp1").geom("geom1").feature("cb13").setIndex("p", 0.06643557019288378, 1, 3);
    model.component("comp1").geom("geom1").create("cb14", "CubicBezier");   // in / output22 / 第 3 段
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", -0.16629830156225683, 0, 0);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", 0.06643557019288378, 1, 0);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", -0.1618089495370118, 0, 1);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", 0.06464708577424685, 1, 1);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", -0.15809430548943199, 0, 2);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", 0.06266050833858641, 1, 2);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", -0.15443897163953926, 0, 3);
    model.component("comp1").geom("geom1").feature("cb14").setIndex("p", 0.06027673395740225, 1, 3);
    model.component("comp1").geom("geom1").create("cb15", "CubicBezier");   // in / output22 / 第 4 段
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", -0.15443897163953926, 0, 0);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", 0.06027673395740225, 1, 0);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", -0.15078363778964654, 0, 1);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", 0.05789295957621809, 1, 1);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", -0.14718761413744094, 0, 2);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", 0.05511198824951022, 1, 2);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", -0.14409551518170757, 0, 3);
    model.component("comp1").geom("geom1").feature("cb15").setIndex("p", 0.052003453456223946, 1, 3);
    model.component("comp1").geom("geom1").create("cb16", "CubicBezier");   // in / output22 / 第 5 段
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", -0.14409551518170757, 0, 0);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", 0.052003453456223946, 1, 0);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", -0.1410034162259742, 0, 1);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", 0.04889491866293767, 1, 1);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", -0.1384152419667131, 0, 2);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", 0.04545882040307301, 1, 2);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", -0.13596292843099772, 0, 3);
    model.component("comp1").geom("geom1").feature("cb16").setIndex("p", 0.042165909613575324, 1, 3);
    model.component("comp1").geom("geom1").create("cb17", "CubicBezier");   // in / output22 / 第 6 段
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", -0.13596292843099772, 0, 0);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", 0.042165909613575324, 1, 0);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", -0.13351061489528235, 0, 1);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", 0.03887299882407764, 1, 1);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", -0.13119416208311274, 0, 2);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", 0.03572327550494694, 1, 2);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", -0.1287888229871091, 0, 3);
    model.component("comp1").geom("geom1").feature("cb17").setIndex("p", 0.03261530883032029, 1, 3);
    model.component("comp1").geom("geom1").create("cb18", "CubicBezier");   // in / output22 / 第 7 段
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", -0.1287888229871091, 0, 0);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", 0.03261530883032029, 1, 0);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", -0.12638348389110543, 0, 1);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", 0.029507342155693634, 1, 1);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", -0.1238892585112677, 0, 2);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", 0.026441132125571032, 1, 2);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", -0.1211609012851205, 0, 3);
    model.component("comp1").geom("geom1").feature("cb18").setIndex("p", 0.02360373316033877, 1, 3);
    model.component("comp1").geom("geom1").create("cb19", "CubicBezier");   // in / output22 / 第 8 段
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", -0.1211609012851205, 0, 0);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", 0.02360373316033877, 1, 0);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", -0.1184325440589733, 0, 1);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", 0.020766334195106505, 1, 1);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", -0.11547005498651665, 0, 2);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", 0.018157746294764583, 1, 2);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", -0.11208899594265007, 0, 3);
    model.component("comp1").geom("geom1").feature("cb19").setIndex("p", 0.015780857406382565, 1, 3);
    model.component("comp1").geom("geom1").create("cb20", "CubicBezier");   // in / output22 / 第 9 段
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", -0.11208899594265007, 0, 0);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", 0.015780857406382565, 1, 0);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", -0.10870793689878348, 0, 1);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", 0.013403968518000547, 1, 1);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", -0.10490830788350697, 0, 2);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", 0.011258778641578435, 1, 2);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", -0.10125063709372815, 0, 3);
    model.component("comp1").geom("geom1").feature("cb20").setIndex("p", 0.010513186603964487, 1, 3);
    model.component("comp1").geom("geom1").create("cb21", "CubicBezier");   // in / output22 / 第 10 段
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", -0.10125063709372815, 0, 0);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", 0.010513186603964487, 1, 0);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", -0.09759296630394934, 0, 1);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", 0.009767594566350539, 1, 1);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", -0.09407725373966822, 0, 2);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", 0.010421600367544756, 1, 2);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", -0.08936694442780331, 0, 3);
    model.component("comp1").geom("geom1").feature("cb21").setIndex("p", 0.009145941898970826, 1, 3);
    model.component("comp1").geom("geom1").create("cb22", "CubicBezier");   // in / output22 / 第 11 段
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", -0.08936694442780331, 0, 0);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", 0.009145941898970826, 1, 0);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", -0.0846566351159384, 0, 1);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", 0.007870283430396896, 1, 1);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", -0.07875172905648971, 0, 2);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", 0.004664960692054816, 1, 2);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", -0.07877388587848679, 0, 3);
    model.component("comp1").geom("geom1").feature("cb22").setIndex("p", 0.0, 1, 3);
    model.component("comp1").geom("geom1").create("cb23", "CubicBezier");   // in / output1 / 第 1 段
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.1932050807568877, 0, 0);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.06535898384862243, 1, 0);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.18948604127893776, 0, 1);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.06321179540490118, 1, 1);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.18540090743613008, 0, 2);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.06730831640982436, 1, 2);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.18072631051327112, 0, 3);
    model.component("comp1").geom("geom1").feature("cb23").setIndex("p", -0.06856138121147924, 1, 3);
    model.component("comp1").geom("geom1").create("cb24", "CubicBezier");   // in / output1 / 第 2 段
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.18072631051327112, 0, 0);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.06856138121147924, 1, 0);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.17605171359041216, 0, 1);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.06981444601313412, 1, 1);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.17078765358750186, 0, 2);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.0682240546115207, 1, 2);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.16629830156225683, 0, 3);
    model.component("comp1").geom("geom1").feature("cb24").setIndex("p", -0.06643557019288378, 1, 3);
    model.component("comp1").geom("geom1").create("cb25", "CubicBezier");   // in / output1 / 第 3 段
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.16629830156225683, 0, 0);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.06643557019288378, 1, 0);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.1618089495370118, 0, 1);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.06464708577424685, 1, 1);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.15809430548943199, 0, 2);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.06266050833858641, 1, 2);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.15443897163953926, 0, 3);
    model.component("comp1").geom("geom1").feature("cb25").setIndex("p", -0.06027673395740225, 1, 3);
    model.component("comp1").geom("geom1").create("cb26", "CubicBezier");   // in / output1 / 第 4 段
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.15443897163953926, 0, 0);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.06027673395740225, 1, 0);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.15078363778964654, 0, 1);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.05789295957621809, 1, 1);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.14718761413744094, 0, 2);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.05511198824951022, 1, 2);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.14409551518170757, 0, 3);
    model.component("comp1").geom("geom1").feature("cb26").setIndex("p", -0.052003453456223946, 1, 3);
    model.component("comp1").geom("geom1").create("cb27", "CubicBezier");   // in / output1 / 第 5 段
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.14409551518170757, 0, 0);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.052003453456223946, 1, 0);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.1410034162259742, 0, 1);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.04889491866293767, 1, 1);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.1384152419667131, 0, 2);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.04545882040307301, 1, 2);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.13596292843099772, 0, 3);
    model.component("comp1").geom("geom1").feature("cb27").setIndex("p", -0.042165909613575324, 1, 3);
    model.component("comp1").geom("geom1").create("cb28", "CubicBezier");   // in / output1 / 第 6 段
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.13596292843099772, 0, 0);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.042165909613575324, 1, 0);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.13351061489528235, 0, 1);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.03887299882407764, 1, 1);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.13119416208311274, 0, 2);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.03572327550494694, 1, 2);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.1287888229871091, 0, 3);
    model.component("comp1").geom("geom1").feature("cb28").setIndex("p", -0.03261530883032029, 1, 3);
    model.component("comp1").geom("geom1").create("cb29", "CubicBezier");   // in / output1 / 第 7 段
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.1287888229871091, 0, 0);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.03261530883032029, 1, 0);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.12638348389110543, 0, 1);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.029507342155693634, 1, 1);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.1238892585112677, 0, 2);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.026441132125571032, 1, 2);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.1211609012851205, 0, 3);
    model.component("comp1").geom("geom1").feature("cb29").setIndex("p", -0.02360373316033877, 1, 3);
    model.component("comp1").geom("geom1").create("cb30", "CubicBezier");   // in / output1 / 第 8 段
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.1211609012851205, 0, 0);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.02360373316033877, 1, 0);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.1184325440589733, 0, 1);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.020766334195106505, 1, 1);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.11547005498651665, 0, 2);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.018157746294764583, 1, 2);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.11208899594265007, 0, 3);
    model.component("comp1").geom("geom1").feature("cb30").setIndex("p", -0.015780857406382565, 1, 3);
    model.component("comp1").geom("geom1").create("cb31", "CubicBezier");   // in / output1 / 第 9 段
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.11208899594265007, 0, 0);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.015780857406382565, 1, 0);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.10870793689878348, 0, 1);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.013403968518000547, 1, 1);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.10490830788350697, 0, 2);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.011258778641578435, 1, 2);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.10125063709372815, 0, 3);
    model.component("comp1").geom("geom1").feature("cb31").setIndex("p", -0.010513186603964487, 1, 3);
    model.component("comp1").geom("geom1").create("cb32", "CubicBezier");   // in / output1 / 第 10 段
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.10125063709372815, 0, 0);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.010513186603964487, 1, 0);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.09759296630394934, 0, 1);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.009767594566350539, 1, 1);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.09407725373966822, 0, 2);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.010421600367544756, 1, 2);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.08936694442780331, 0, 3);
    model.component("comp1").geom("geom1").feature("cb32").setIndex("p", -0.009145941898970826, 1, 3);
    model.component("comp1").geom("geom1").create("cb33", "CubicBezier");   // in / output1 / 第 11 段
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.08936694442780331, 0, 0);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.009145941898970826, 1, 0);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.0846566351159384, 0, 1);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.007870283430396896, 1, 1);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.07875172905648971, 0, 2);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.004664960692054816, 1, 2);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", -0.07877388587848679, 0, 3);
    model.component("comp1").geom("geom1").feature("cb33").setIndex("p", 0.0, 1, 3);
    model.component("comp1").geom("geom1").create("cb34", "CubicBezier");   // in / output21 / 第 1 段
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", 0.1932050807568877, 0, 0);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", -0.06535898384862243, 1, 0);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", 0.18948604127893776, 0, 1);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", -0.06321179540490118, 1, 1);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", 0.18540090743613008, 0, 2);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", -0.06730831640982436, 1, 2);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", 0.18072631051327112, 0, 3);
    model.component("comp1").geom("geom1").feature("cb34").setIndex("p", -0.06856138121147924, 1, 3);
    model.component("comp1").geom("geom1").create("cb35", "CubicBezier");   // in / output21 / 第 2 段
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", 0.18072631051327112, 0, 0);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", -0.06856138121147924, 1, 0);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", 0.17605171359041216, 0, 1);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", -0.06981444601313412, 1, 1);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", 0.17078765358750186, 0, 2);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", -0.0682240546115207, 1, 2);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", 0.16629830156225683, 0, 3);
    model.component("comp1").geom("geom1").feature("cb35").setIndex("p", -0.06643557019288378, 1, 3);
    model.component("comp1").geom("geom1").create("cb36", "CubicBezier");   // in / output21 / 第 3 段
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", 0.16629830156225683, 0, 0);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", -0.06643557019288378, 1, 0);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", 0.1618089495370118, 0, 1);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", -0.06464708577424685, 1, 1);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", 0.15809430548943199, 0, 2);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", -0.06266050833858641, 1, 2);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", 0.15443897163953926, 0, 3);
    model.component("comp1").geom("geom1").feature("cb36").setIndex("p", -0.06027673395740225, 1, 3);
    model.component("comp1").geom("geom1").create("cb37", "CubicBezier");   // in / output21 / 第 4 段
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", 0.15443897163953926, 0, 0);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", -0.06027673395740225, 1, 0);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", 0.15078363778964654, 0, 1);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", -0.05789295957621809, 1, 1);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", 0.14718761413744094, 0, 2);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", -0.05511198824951022, 1, 2);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", 0.14409551518170757, 0, 3);
    model.component("comp1").geom("geom1").feature("cb37").setIndex("p", -0.052003453456223946, 1, 3);
    model.component("comp1").geom("geom1").create("cb38", "CubicBezier");   // in / output21 / 第 5 段
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", 0.14409551518170757, 0, 0);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", -0.052003453456223946, 1, 0);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", 0.1410034162259742, 0, 1);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", -0.04889491866293767, 1, 1);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", 0.1384152419667131, 0, 2);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", -0.04545882040307301, 1, 2);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", 0.13596292843099772, 0, 3);
    model.component("comp1").geom("geom1").feature("cb38").setIndex("p", -0.042165909613575324, 1, 3);
    model.component("comp1").geom("geom1").create("cb39", "CubicBezier");   // in / output21 / 第 6 段
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", 0.13596292843099772, 0, 0);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", -0.042165909613575324, 1, 0);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", 0.13351061489528235, 0, 1);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", -0.03887299882407764, 1, 1);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", 0.13119416208311274, 0, 2);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", -0.03572327550494694, 1, 2);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", 0.1287888229871091, 0, 3);
    model.component("comp1").geom("geom1").feature("cb39").setIndex("p", -0.03261530883032029, 1, 3);
    model.component("comp1").geom("geom1").create("cb40", "CubicBezier");   // in / output21 / 第 7 段
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", 0.1287888229871091, 0, 0);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", -0.03261530883032029, 1, 0);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", 0.12638348389110543, 0, 1);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", -0.029507342155693634, 1, 1);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", 0.1238892585112677, 0, 2);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", -0.026441132125571032, 1, 2);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", 0.1211609012851205, 0, 3);
    model.component("comp1").geom("geom1").feature("cb40").setIndex("p", -0.02360373316033877, 1, 3);
    model.component("comp1").geom("geom1").create("cb41", "CubicBezier");   // in / output21 / 第 8 段
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", 0.1211609012851205, 0, 0);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", -0.02360373316033877, 1, 0);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", 0.1184325440589733, 0, 1);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", -0.020766334195106505, 1, 1);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", 0.11547005498651665, 0, 2);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", -0.018157746294764583, 1, 2);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", 0.11208899594265007, 0, 3);
    model.component("comp1").geom("geom1").feature("cb41").setIndex("p", -0.015780857406382565, 1, 3);
    model.component("comp1").geom("geom1").create("cb42", "CubicBezier");   // in / output21 / 第 9 段
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", 0.11208899594265007, 0, 0);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", -0.015780857406382565, 1, 0);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", 0.10870793689878348, 0, 1);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", -0.013403968518000547, 1, 1);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", 0.10490830788350697, 0, 2);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", -0.011258778641578435, 1, 2);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", 0.10125063709372815, 0, 3);
    model.component("comp1").geom("geom1").feature("cb42").setIndex("p", -0.010513186603964487, 1, 3);
    model.component("comp1").geom("geom1").create("cb43", "CubicBezier");   // in / output21 / 第 10 段
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", 0.10125063709372815, 0, 0);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", -0.010513186603964487, 1, 0);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", 0.09759296630394934, 0, 1);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", -0.009767594566350539, 1, 1);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", 0.09407725373966822, 0, 2);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", -0.010421600367544756, 1, 2);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", 0.08936694442780331, 0, 3);
    model.component("comp1").geom("geom1").feature("cb43").setIndex("p", -0.009145941898970826, 1, 3);
    model.component("comp1").geom("geom1").create("cb44", "CubicBezier");   // in / output21 / 第 11 段
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", 0.08936694442780331, 0, 0);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", -0.009145941898970826, 1, 0);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", 0.0846566351159384, 0, 1);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", -0.007870283430396896, 1, 1);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", 0.07875172905648971, 0, 2);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", -0.004664960692054816, 1, 2);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", 0.07877388587848679, 0, 3);
    model.component("comp1").geom("geom1").feature("cb44").setIndex("p", 0.0, 1, 3);

    // 优化壁 ex：4 个象限 x 11 段三次 Bézier（控制点 = 控制多边形）
    model.component("comp1").geom("geom1").create("cb45", "CubicBezier");   // ex / input / 第 1 段
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.15320508075688766, 0, 0);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.13464101615137752, 1, 0);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.14845453031738873, 0, 1);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.13189828457633393, 1, 1);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.14311222642213256, 0, 2);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.13048974228906687, 1, 2);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.13790467738852558, 0, 3);
    model.component("comp1").geom("geom1").feature("cb45").setIndex("p", 0.1285821044913435, 1, 3);
    model.component("comp1").geom("geom1").create("cb46", "CubicBezier");   // ex / input / 第 2 段
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.13790467738852558, 0, 0);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.1285821044913435, 1, 0);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.1326971283549186, 0, 1);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.12667446669362012, 1, 1);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.12762433418296082, 0, 2);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.12426773338544043, 1, 2);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.12326274126960114, 0, 3);
    model.component("comp1").geom("geom1").feature("cb46").setIndex("p", 0.12147426934712238, 1, 3);
    model.component("comp1").geom("geom1").create("cb47", "CubicBezier");   // ex / input / 第 3 段
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.12326274126960114, 0, 0);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.12147426934712238, 1, 0);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.11890114835624147, 0, 1);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.11868080530880433, 1, 1);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.11525075670147987, 0, 2);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.11550061054034792, 1, 2);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.11154356908661658, 0, 3);
    model.component("comp1").geom("geom1").feature("cb47").setIndex("p", 0.11218947889452197, 1, 3);
    model.component("comp1").geom("geom1").create("cb48", "CubicBezier");   // ex / input / 第 4 段
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.11154356908661658, 0, 0);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.11218947889452197, 1, 0);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.1078363814717533, 0, 1);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.10887834724869601, 1, 1);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.10407239789678834, 0, 2);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.10543627872550049, 1, 2);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.1006210834800041, 0, 3);
    model.component("comp1").geom("geom1").feature("cb48").setIndex("p", 0.1018595470588904, 1, 3);
    model.component("comp1").geom("geom1").create("cb49", "CubicBezier");   // ex / input / 第 5 段
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.1006210834800041, 0, 0);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.1018595470588904, 1, 0);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.09716976906321988, 0, 1);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.09828281539228031, 1, 1);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.0940311238046164, 0, 2);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.09457142058225565, 1, 2);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.08993117368509654, 0, 3);
    model.component("comp1").geom("geom1").feature("cb49").setIndex("p", 0.09023665122235175, 1, 3);
    model.component("comp1").geom("geom1").create("cb50", "CubicBezier");   // ex / input / 第 6 段
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.08993117368509654, 0, 0);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.09023665122235175, 1, 0);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.08583122356557668, 0, 1);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.08590188186244785, 1, 1);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.08076996858514048, 0, 2);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.0809437379526647, 1, 2);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.07485053558341831, 0, 3);
    model.component("comp1").geom("geom1").feature("cb50").setIndex("p", 0.07785964362901006, 1, 3);
    model.component("comp1").geom("geom1").create("cb51", "CubicBezier");   // ex / input / 第 7 段
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.07485053558341831, 0, 0);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.07785964362901006, 1, 0);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.06893110258169614, 0, 1);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.07477554930535542, 1, 1);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.062153491558688045, 0, 2);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.0735655045678293, 1, 2);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.055979918128120934, 0, 3);
    model.component("comp1").geom("geom1").feature("cb51").setIndex("p", 0.07337671856094037, 1, 3);
    model.component("comp1").geom("geom1").create("cb52", "CubicBezier");   // ex / input / 第 8 段
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.055979918128120934, 0, 0);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.07337671856094037, 1, 0);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.04980634469755382, 0, 1);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.07318793255405144, 1, 1);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.04423680885942772, 0, 2);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.07402040527779971, 1, 2);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.04051820879411936, 0, 3);
    model.component("comp1").geom("geom1").feature("cb52").setIndex("p", 0.07478040775334212, 1, 3);
    model.component("comp1").geom("geom1").create("cb53", "CubicBezier");   // ex / input / 第 9 段
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.04051820879411936, 0, 0);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.07478040775334212, 1, 0);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.036799608728811, 0, 1);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.07554041022888452, 1, 1);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.03493194443632041, 0, 2);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.07622794245622105, 1, 2);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.028555320090624104, 0, 3);
    model.component("comp1").geom("geom1").feature("cb53").setIndex("p", 0.07788559309172745, 1, 3);
    model.component("comp1").geom("geom1").create("cb54", "CubicBezier");   // ex / input / 第 10 段
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.028555320090624104, 0, 0);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.07788559309172745, 1, 0);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.0221786957449278, 0, 1);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.07954324372723384, 1, 1);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.011293111346025786, 0, 2);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.08217101277091009, 1, 2);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.0069107481079664825, 0, 3);
    model.component("comp1").geom("geom1").feature("cb54").setIndex("p", 0.07902935883592173, 1, 3);
    model.component("comp1").geom("geom1").create("cb55", "CubicBezier");   // ex / input / 第 11 段
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.0069107481079664825, 0, 0);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.07902935883592173, 1, 0);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.0025283848699071782, 0, 1);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.07588770490093337, 1, 1);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.004649242792690582, 0, 2);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.06697662798728041, 1, 2);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.0, 0, 3);
    model.component("comp1").geom("geom1").feature("cb55").setIndex("p", 0.06691252465224272, 1, 3);
    model.component("comp1").geom("geom1").create("cb56", "CubicBezier");   // ex / output22 / 第 1 段
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", -0.15320508075688766, 0, 0);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", 0.13464101615137752, 1, 0);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", -0.14845453031738873, 0, 1);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", 0.13189828457633393, 1, 1);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", -0.14311222642213256, 0, 2);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", 0.13048974228906687, 1, 2);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", -0.13790467738852558, 0, 3);
    model.component("comp1").geom("geom1").feature("cb56").setIndex("p", 0.1285821044913435, 1, 3);
    model.component("comp1").geom("geom1").create("cb57", "CubicBezier");   // ex / output22 / 第 2 段
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", -0.13790467738852558, 0, 0);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", 0.1285821044913435, 1, 0);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", -0.1326971283549186, 0, 1);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", 0.12667446669362012, 1, 1);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", -0.12762433418296082, 0, 2);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", 0.12426773338544043, 1, 2);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", -0.12326274126960114, 0, 3);
    model.component("comp1").geom("geom1").feature("cb57").setIndex("p", 0.12147426934712238, 1, 3);
    model.component("comp1").geom("geom1").create("cb58", "CubicBezier");   // ex / output22 / 第 3 段
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", -0.12326274126960114, 0, 0);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", 0.12147426934712238, 1, 0);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", -0.11890114835624147, 0, 1);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", 0.11868080530880433, 1, 1);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", -0.11525075670147987, 0, 2);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", 0.11550061054034792, 1, 2);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", -0.11154356908661658, 0, 3);
    model.component("comp1").geom("geom1").feature("cb58").setIndex("p", 0.11218947889452197, 1, 3);
    model.component("comp1").geom("geom1").create("cb59", "CubicBezier");   // ex / output22 / 第 4 段
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", -0.11154356908661658, 0, 0);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", 0.11218947889452197, 1, 0);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", -0.1078363814717533, 0, 1);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", 0.10887834724869601, 1, 1);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", -0.10407239789678834, 0, 2);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", 0.10543627872550049, 1, 2);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", -0.1006210834800041, 0, 3);
    model.component("comp1").geom("geom1").feature("cb59").setIndex("p", 0.1018595470588904, 1, 3);
    model.component("comp1").geom("geom1").create("cb60", "CubicBezier");   // ex / output22 / 第 5 段
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", -0.1006210834800041, 0, 0);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", 0.1018595470588904, 1, 0);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", -0.09716976906321988, 0, 1);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", 0.09828281539228031, 1, 1);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", -0.0940311238046164, 0, 2);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", 0.09457142058225565, 1, 2);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", -0.08993117368509654, 0, 3);
    model.component("comp1").geom("geom1").feature("cb60").setIndex("p", 0.09023665122235175, 1, 3);
    model.component("comp1").geom("geom1").create("cb61", "CubicBezier");   // ex / output22 / 第 6 段
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", -0.08993117368509654, 0, 0);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", 0.09023665122235175, 1, 0);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", -0.08583122356557668, 0, 1);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", 0.08590188186244785, 1, 1);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", -0.08076996858514048, 0, 2);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", 0.0809437379526647, 1, 2);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", -0.07485053558341831, 0, 3);
    model.component("comp1").geom("geom1").feature("cb61").setIndex("p", 0.07785964362901006, 1, 3);
    model.component("comp1").geom("geom1").create("cb62", "CubicBezier");   // ex / output22 / 第 7 段
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", -0.07485053558341831, 0, 0);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", 0.07785964362901006, 1, 0);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", -0.06893110258169614, 0, 1);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", 0.07477554930535542, 1, 1);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", -0.062153491558688045, 0, 2);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", 0.0735655045678293, 1, 2);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", -0.055979918128120934, 0, 3);
    model.component("comp1").geom("geom1").feature("cb62").setIndex("p", 0.07337671856094037, 1, 3);
    model.component("comp1").geom("geom1").create("cb63", "CubicBezier");   // ex / output22 / 第 8 段
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", -0.055979918128120934, 0, 0);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", 0.07337671856094037, 1, 0);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", -0.04980634469755382, 0, 1);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", 0.07318793255405144, 1, 1);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", -0.04423680885942772, 0, 2);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", 0.07402040527779971, 1, 2);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", -0.04051820879411936, 0, 3);
    model.component("comp1").geom("geom1").feature("cb63").setIndex("p", 0.07478040775334212, 1, 3);
    model.component("comp1").geom("geom1").create("cb64", "CubicBezier");   // ex / output22 / 第 9 段
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", -0.04051820879411936, 0, 0);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", 0.07478040775334212, 1, 0);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", -0.036799608728811, 0, 1);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", 0.07554041022888452, 1, 1);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", -0.03493194443632041, 0, 2);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", 0.07622794245622105, 1, 2);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", -0.028555320090624104, 0, 3);
    model.component("comp1").geom("geom1").feature("cb64").setIndex("p", 0.07788559309172745, 1, 3);
    model.component("comp1").geom("geom1").create("cb65", "CubicBezier");   // ex / output22 / 第 10 段
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", -0.028555320090624104, 0, 0);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", 0.07788559309172745, 1, 0);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", -0.0221786957449278, 0, 1);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", 0.07954324372723384, 1, 1);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", -0.011293111346025786, 0, 2);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", 0.08217101277091009, 1, 2);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", -0.0069107481079664825, 0, 3);
    model.component("comp1").geom("geom1").feature("cb65").setIndex("p", 0.07902935883592173, 1, 3);
    model.component("comp1").geom("geom1").create("cb66", "CubicBezier");   // ex / output22 / 第 11 段
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", -0.0069107481079664825, 0, 0);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", 0.07902935883592173, 1, 0);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", -0.0025283848699071782, 0, 1);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", 0.07588770490093337, 1, 1);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", -0.004649242792690582, 0, 2);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", 0.06697662798728041, 1, 2);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", 0.0, 0, 3);
    model.component("comp1").geom("geom1").feature("cb66").setIndex("p", 0.06691252465224272, 1, 3);
    model.component("comp1").geom("geom1").create("cb67", "CubicBezier");   // ex / output1 / 第 1 段
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.15320508075688766, 0, 0);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.13464101615137752, 1, 0);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.14845453031738873, 0, 1);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.13189828457633393, 1, 1);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.14311222642213256, 0, 2);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.13048974228906687, 1, 2);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.13790467738852558, 0, 3);
    model.component("comp1").geom("geom1").feature("cb67").setIndex("p", -0.1285821044913435, 1, 3);
    model.component("comp1").geom("geom1").create("cb68", "CubicBezier");   // ex / output1 / 第 2 段
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.13790467738852558, 0, 0);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.1285821044913435, 1, 0);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.1326971283549186, 0, 1);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.12667446669362012, 1, 1);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.12762433418296082, 0, 2);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.12426773338544043, 1, 2);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.12326274126960114, 0, 3);
    model.component("comp1").geom("geom1").feature("cb68").setIndex("p", -0.12147426934712238, 1, 3);
    model.component("comp1").geom("geom1").create("cb69", "CubicBezier");   // ex / output1 / 第 3 段
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.12326274126960114, 0, 0);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.12147426934712238, 1, 0);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.11890114835624147, 0, 1);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.11868080530880433, 1, 1);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.11525075670147987, 0, 2);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.11550061054034792, 1, 2);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.11154356908661658, 0, 3);
    model.component("comp1").geom("geom1").feature("cb69").setIndex("p", -0.11218947889452197, 1, 3);
    model.component("comp1").geom("geom1").create("cb70", "CubicBezier");   // ex / output1 / 第 4 段
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.11154356908661658, 0, 0);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.11218947889452197, 1, 0);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.1078363814717533, 0, 1);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.10887834724869601, 1, 1);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.10407239789678834, 0, 2);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.10543627872550049, 1, 2);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.1006210834800041, 0, 3);
    model.component("comp1").geom("geom1").feature("cb70").setIndex("p", -0.1018595470588904, 1, 3);
    model.component("comp1").geom("geom1").create("cb71", "CubicBezier");   // ex / output1 / 第 5 段
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.1006210834800041, 0, 0);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.1018595470588904, 1, 0);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.09716976906321988, 0, 1);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.09828281539228031, 1, 1);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.0940311238046164, 0, 2);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.09457142058225565, 1, 2);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.08993117368509654, 0, 3);
    model.component("comp1").geom("geom1").feature("cb71").setIndex("p", -0.09023665122235175, 1, 3);
    model.component("comp1").geom("geom1").create("cb72", "CubicBezier");   // ex / output1 / 第 6 段
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.08993117368509654, 0, 0);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.09023665122235175, 1, 0);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.08583122356557668, 0, 1);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.08590188186244785, 1, 1);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.08076996858514048, 0, 2);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.0809437379526647, 1, 2);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.07485053558341831, 0, 3);
    model.component("comp1").geom("geom1").feature("cb72").setIndex("p", -0.07785964362901006, 1, 3);
    model.component("comp1").geom("geom1").create("cb73", "CubicBezier");   // ex / output1 / 第 7 段
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.07485053558341831, 0, 0);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.07785964362901006, 1, 0);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.06893110258169614, 0, 1);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.07477554930535542, 1, 1);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.062153491558688045, 0, 2);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.0735655045678293, 1, 2);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.055979918128120934, 0, 3);
    model.component("comp1").geom("geom1").feature("cb73").setIndex("p", -0.07337671856094037, 1, 3);
    model.component("comp1").geom("geom1").create("cb74", "CubicBezier");   // ex / output1 / 第 8 段
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.055979918128120934, 0, 0);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.07337671856094037, 1, 0);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.04980634469755382, 0, 1);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.07318793255405144, 1, 1);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.04423680885942772, 0, 2);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.07402040527779971, 1, 2);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.04051820879411936, 0, 3);
    model.component("comp1").geom("geom1").feature("cb74").setIndex("p", -0.07478040775334212, 1, 3);
    model.component("comp1").geom("geom1").create("cb75", "CubicBezier");   // ex / output1 / 第 9 段
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.04051820879411936, 0, 0);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.07478040775334212, 1, 0);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.036799608728811, 0, 1);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.07554041022888452, 1, 1);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.03493194443632041, 0, 2);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.07622794245622105, 1, 2);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.028555320090624104, 0, 3);
    model.component("comp1").geom("geom1").feature("cb75").setIndex("p", -0.07788559309172745, 1, 3);
    model.component("comp1").geom("geom1").create("cb76", "CubicBezier");   // ex / output1 / 第 10 段
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.028555320090624104, 0, 0);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.07788559309172745, 1, 0);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.0221786957449278, 0, 1);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.07954324372723384, 1, 1);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.011293111346025786, 0, 2);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.08217101277091009, 1, 2);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.0069107481079664825, 0, 3);
    model.component("comp1").geom("geom1").feature("cb76").setIndex("p", -0.07902935883592173, 1, 3);
    model.component("comp1").geom("geom1").create("cb77", "CubicBezier");   // ex / output1 / 第 11 段
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.0069107481079664825, 0, 0);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.07902935883592173, 1, 0);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.0025283848699071782, 0, 1);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.07588770490093337, 1, 1);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.004649242792690582, 0, 2);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.06697662798728041, 1, 2);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", 0.0, 0, 3);
    model.component("comp1").geom("geom1").feature("cb77").setIndex("p", -0.06691252465224272, 1, 3);
    model.component("comp1").geom("geom1").create("cb78", "CubicBezier");   // ex / output21 / 第 1 段
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", 0.15320508075688766, 0, 0);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", -0.13464101615137752, 1, 0);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", 0.14845453031738873, 0, 1);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", -0.13189828457633393, 1, 1);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", 0.14311222642213256, 0, 2);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", -0.13048974228906687, 1, 2);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", 0.13790467738852558, 0, 3);
    model.component("comp1").geom("geom1").feature("cb78").setIndex("p", -0.1285821044913435, 1, 3);
    model.component("comp1").geom("geom1").create("cb79", "CubicBezier");   // ex / output21 / 第 2 段
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", 0.13790467738852558, 0, 0);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", -0.1285821044913435, 1, 0);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", 0.1326971283549186, 0, 1);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", -0.12667446669362012, 1, 1);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", 0.12762433418296082, 0, 2);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", -0.12426773338544043, 1, 2);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", 0.12326274126960114, 0, 3);
    model.component("comp1").geom("geom1").feature("cb79").setIndex("p", -0.12147426934712238, 1, 3);
    model.component("comp1").geom("geom1").create("cb80", "CubicBezier");   // ex / output21 / 第 3 段
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", 0.12326274126960114, 0, 0);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", -0.12147426934712238, 1, 0);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", 0.11890114835624147, 0, 1);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", -0.11868080530880433, 1, 1);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", 0.11525075670147987, 0, 2);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", -0.11550061054034792, 1, 2);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", 0.11154356908661658, 0, 3);
    model.component("comp1").geom("geom1").feature("cb80").setIndex("p", -0.11218947889452197, 1, 3);
    model.component("comp1").geom("geom1").create("cb81", "CubicBezier");   // ex / output21 / 第 4 段
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", 0.11154356908661658, 0, 0);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", -0.11218947889452197, 1, 0);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", 0.1078363814717533, 0, 1);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", -0.10887834724869601, 1, 1);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", 0.10407239789678834, 0, 2);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", -0.10543627872550049, 1, 2);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", 0.1006210834800041, 0, 3);
    model.component("comp1").geom("geom1").feature("cb81").setIndex("p", -0.1018595470588904, 1, 3);
    model.component("comp1").geom("geom1").create("cb82", "CubicBezier");   // ex / output21 / 第 5 段
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", 0.1006210834800041, 0, 0);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", -0.1018595470588904, 1, 0);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", 0.09716976906321988, 0, 1);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", -0.09828281539228031, 1, 1);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", 0.0940311238046164, 0, 2);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", -0.09457142058225565, 1, 2);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", 0.08993117368509654, 0, 3);
    model.component("comp1").geom("geom1").feature("cb82").setIndex("p", -0.09023665122235175, 1, 3);
    model.component("comp1").geom("geom1").create("cb83", "CubicBezier");   // ex / output21 / 第 6 段
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", 0.08993117368509654, 0, 0);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", -0.09023665122235175, 1, 0);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", 0.08583122356557668, 0, 1);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", -0.08590188186244785, 1, 1);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", 0.08076996858514048, 0, 2);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", -0.0809437379526647, 1, 2);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", 0.07485053558341831, 0, 3);
    model.component("comp1").geom("geom1").feature("cb83").setIndex("p", -0.07785964362901006, 1, 3);
    model.component("comp1").geom("geom1").create("cb84", "CubicBezier");   // ex / output21 / 第 7 段
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", 0.07485053558341831, 0, 0);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", -0.07785964362901006, 1, 0);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", 0.06893110258169614, 0, 1);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", -0.07477554930535542, 1, 1);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", 0.062153491558688045, 0, 2);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", -0.0735655045678293, 1, 2);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", 0.055979918128120934, 0, 3);
    model.component("comp1").geom("geom1").feature("cb84").setIndex("p", -0.07337671856094037, 1, 3);
    model.component("comp1").geom("geom1").create("cb85", "CubicBezier");   // ex / output21 / 第 8 段
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", 0.055979918128120934, 0, 0);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", -0.07337671856094037, 1, 0);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", 0.04980634469755382, 0, 1);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", -0.07318793255405144, 1, 1);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", 0.04423680885942772, 0, 2);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", -0.07402040527779971, 1, 2);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", 0.04051820879411936, 0, 3);
    model.component("comp1").geom("geom1").feature("cb85").setIndex("p", -0.07478040775334212, 1, 3);
    model.component("comp1").geom("geom1").create("cb86", "CubicBezier");   // ex / output21 / 第 9 段
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", 0.04051820879411936, 0, 0);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", -0.07478040775334212, 1, 0);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", 0.036799608728811, 0, 1);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", -0.07554041022888452, 1, 1);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", 0.03493194443632041, 0, 2);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", -0.07622794245622105, 1, 2);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", 0.028555320090624104, 0, 3);
    model.component("comp1").geom("geom1").feature("cb86").setIndex("p", -0.07788559309172745, 1, 3);
    model.component("comp1").geom("geom1").create("cb87", "CubicBezier");   // ex / output21 / 第 10 段
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", 0.028555320090624104, 0, 0);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", -0.07788559309172745, 1, 0);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", 0.0221786957449278, 0, 1);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", -0.07954324372723384, 1, 1);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", 0.011293111346025786, 0, 2);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", -0.08217101277091009, 1, 2);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", 0.0069107481079664825, 0, 3);
    model.component("comp1").geom("geom1").feature("cb87").setIndex("p", -0.07902935883592173, 1, 3);
    model.component("comp1").geom("geom1").create("cb88", "CubicBezier");   // ex / output21 / 第 11 段
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", 0.0069107481079664825, 0, 0);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", -0.07902935883592173, 1, 0);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", 0.0025283848699071782, 0, 1);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", -0.07588770490093337, 1, 1);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", 0.004649242792690582, 0, 2);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", -0.06697662798728041, 1, 2);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", 0.0, 0, 3);
    model.component("comp1").geom("geom1").feature("cb88").setIndex("p", -0.06691252465224272, 1, 3);

    // 固定直壁段（长 0.492820323027551）：波导口 -> 壁的固定端
    model.component("comp1").geom("geom1").create("ls1", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls1").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls1").set("coord1", new String[]{"0.1932050807568877", "0.06535898384862243"});
    model.component("comp1").geom("geom1").feature("ls1").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls1").set("coord2", new String[]{"0.6200000000000001", "0.3117691453623979"});
    model.component("comp1").geom("geom1").create("ls2", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls2").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls2").set("coord1", new String[]{"0.15320508075688766", "0.13464101615137752"});
    model.component("comp1").geom("geom1").feature("ls2").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls2").set("coord2", new String[]{"0.5800000000000001", "0.381051177665153"});
    model.component("comp1").geom("geom1").create("ls3", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls3").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls3").set("coord1", new String[]{"-0.1932050807568877", "0.06535898384862243"});
    model.component("comp1").geom("geom1").feature("ls3").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls3").set("coord2", new String[]{"-0.6200000000000001", "0.3117691453623979"});
    model.component("comp1").geom("geom1").create("ls4", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls4").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls4").set("coord1", new String[]{"-0.15320508075688766", "0.13464101615137752"});
    model.component("comp1").geom("geom1").feature("ls4").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls4").set("coord2", new String[]{"-0.5800000000000001", "0.381051177665153"});
    model.component("comp1").geom("geom1").create("ls5", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls5").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls5").set("coord1", new String[]{"-0.1932050807568877", "-0.06535898384862243"});
    model.component("comp1").geom("geom1").feature("ls5").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls5").set("coord2", new String[]{"-0.6200000000000001", "-0.3117691453623979"});
    model.component("comp1").geom("geom1").create("ls6", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls6").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls6").set("coord1", new String[]{"-0.15320508075688766", "-0.13464101615137752"});
    model.component("comp1").geom("geom1").feature("ls6").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls6").set("coord2", new String[]{"-0.5800000000000001", "-0.381051177665153"});
    model.component("comp1").geom("geom1").create("ls7", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls7").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls7").set("coord1", new String[]{"0.1932050807568877", "-0.06535898384862243"});
    model.component("comp1").geom("geom1").feature("ls7").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls7").set("coord2", new String[]{"0.6200000000000001", "-0.3117691453623979"});
    model.component("comp1").geom("geom1").create("ls8", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls8").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls8").set("coord1", new String[]{"0.15320508075688766", "-0.13464101615137752"});
    model.component("comp1").geom("geom1").feature("ls8").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls8").set("coord2", new String[]{"0.5800000000000001", "-0.381051177665153"});

    // 端口边界的印记线：与六边形边共线，把每条边切成 5 段（只有中间 3 段是端口）。
    // 参考文件的 ls1 干的正是这件事。不加的话六边形边只会被波导口切成 3 段，
    // Box 就框不出“中间 6w”这个集合了。
    model.component("comp1").geom("geom1").create("ls9", "LineSegment");   // 端口印记线 input
    model.component("comp1").geom("geom1").feature("ls9").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls9").set("coord1", new String[]{"0.7200000000000001", "0.13856406460551018"});
    model.component("comp1").geom("geom1").feature("ls9").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls9").set("coord2", new String[]{"0.4800000000000001", "0.5542562584220407"});
    model.component("comp1").geom("geom1").create("ls10", "LineSegment");   // 端口印记线 output22
    model.component("comp1").geom("geom1").feature("ls10").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls10").set("coord1", new String[]{"-0.7200000000000001", "0.13856406460551018"});
    model.component("comp1").geom("geom1").feature("ls10").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls10").set("coord2", new String[]{"-0.4800000000000001", "0.5542562584220407"});
    model.component("comp1").geom("geom1").create("ls11", "LineSegment");   // 端口印记线 output1
    model.component("comp1").geom("geom1").feature("ls11").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls11").set("coord1", new String[]{"-0.7200000000000001", "-0.13856406460551018"});
    model.component("comp1").geom("geom1").feature("ls11").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls11").set("coord2", new String[]{"-0.4800000000000001", "-0.5542562584220407"});
    model.component("comp1").geom("geom1").create("ls12", "LineSegment");   // 端口印记线 output21
    model.component("comp1").geom("geom1").feature("ls12").set("specify1", "coord");
    model.component("comp1").geom("geom1").feature("ls12").set("coord1", new String[]{"0.7200000000000001", "-0.13856406460551018"});
    model.component("comp1").geom("geom1").feature("ls12").set("specify2", "coord");
    model.component("comp1").geom("geom1").feature("ls12").set("coord2", new String[]{"0.4800000000000001", "-0.5542562584220407"});

    // Form Union：这些曲线把六边形切成 5 个域（1 个波导 + 4 个空气楔）。
    // 原点附近那 4 条同材料分割线故意不放 —— 它们会把波导切成 4 块，对物理没影响。
    model.component("comp1").geom("geom1").run();

    // ---------------- 按坐标的 Box 选择 ----------------
    // 端口 input：3 条共线边，中点 (+0.6000, +0.3464)
    model.component("comp1").selection().create("sel_input", "Box");
    model.component("comp1").selection("sel_input").set("entitydim", 1);
    model.component("comp1").selection("sel_input").set("xmin", 0.4640000000000001);
    model.component("comp1").selection("sel_input").set("xmax", 0.7360000000000001);
    model.component("comp1").selection("sel_input").set("ymin", 0.12256406460551018);
    model.component("comp1").selection("sel_input").set("ymax", 0.5702562584220408);
    model.component("comp1").selection("sel_input").set("condition", "inside");
    // 端口 output1：3 条共线边，中点 (-0.6000, -0.3464)
    model.component("comp1").selection().create("sel_output1", "Box");
    model.component("comp1").selection("sel_output1").set("entitydim", 1);
    model.component("comp1").selection("sel_output1").set("xmin", -0.7360000000000001);
    model.component("comp1").selection("sel_output1").set("xmax", -0.4640000000000001);
    model.component("comp1").selection("sel_output1").set("ymin", -0.5702562584220408);
    model.component("comp1").selection("sel_output1").set("ymax", -0.12256406460551018);
    model.component("comp1").selection("sel_output1").set("condition", "inside");
    // 端口 output21：3 条共线边，中点 (+0.6000, -0.3464)
    model.component("comp1").selection().create("sel_output21", "Box");
    model.component("comp1").selection("sel_output21").set("entitydim", 1);
    model.component("comp1").selection("sel_output21").set("xmin", 0.4640000000000001);
    model.component("comp1").selection("sel_output21").set("xmax", 0.7360000000000001);
    model.component("comp1").selection("sel_output21").set("ymin", -0.5702562584220408);
    model.component("comp1").selection("sel_output21").set("ymax", -0.12256406460551018);
    model.component("comp1").selection("sel_output21").set("condition", "inside");
    // 端口 output22：3 条共线边，中点 (-0.6000, +0.3464)
    model.component("comp1").selection().create("sel_output22", "Box");
    model.component("comp1").selection("sel_output22").set("entitydim", 1);
    model.component("comp1").selection("sel_output22").set("xmin", -0.7360000000000001);
    model.component("comp1").selection("sel_output22").set("xmax", -0.4640000000000001);
    model.component("comp1").selection("sel_output22").set("ymin", 0.12256406460551018);
    model.component("comp1").selection("sel_output22").set("ymax", 0.5702562584220408);
    model.component("comp1").selection("sel_output22").set("condition", "inside");
    // 散射边界：y = +0.6928 那条整边 + 它两端相邻的 2 段端部
    model.component("comp1").selection().create("sel_scatter_top", "Box");
    model.component("comp1").selection("sel_scatter_top").set("entitydim", 1);
    model.component("comp1").selection("sel_scatter_top").set("xmin", -0.4960000000000001);
    model.component("comp1").selection("sel_scatter_top").set("xmax", 0.4960000000000001);
    model.component("comp1").selection("sel_scatter_top").set("ymin", 0.5382562584220407);
    model.component("comp1").selection("sel_scatter_top").set("ymax", 0.708820323027551);
    model.component("comp1").selection("sel_scatter_top").set("condition", "inside");
    // 散射边界：y = -0.6928 那条整边 + 它两端相邻的 2 段端部
    model.component("comp1").selection().create("sel_scatter_bottom", "Box");
    model.component("comp1").selection("sel_scatter_bottom").set("entitydim", 1);
    model.component("comp1").selection("sel_scatter_bottom").set("xmin", -0.4960000000000001);
    model.component("comp1").selection("sel_scatter_bottom").set("xmax", 0.4960000000000001);
    model.component("comp1").selection("sel_scatter_bottom").set("ymin", -0.7088203230275513);
    model.component("comp1").selection("sel_scatter_bottom").set("ymax", -0.5382562584220407);
    model.component("comp1").selection("sel_scatter_bottom").set("condition", "inside");
    // 波导域（框落在直段内部）
    model.component("comp1").selection().create("sel_wg", "Box");
    model.component("comp1").selection("sel_wg").set("entitydim", 2);
    model.component("comp1").selection("sel_wg").set("xmin", 0.15320508075688768);
    model.component("comp1").selection("sel_wg").set("xmax", 0.19320508075688766);
    model.component("comp1").selection("sel_wg").set("ymin", 0.07999999999999997);
    model.component("comp1").selection("sel_wg").set("ymax", 0.11999999999999998);
    model.component("comp1").selection("sel_wg").set("condition", "intersects");
    // 空气楔 0°（两臂之间，六边形顶点侧）
    model.component("comp1").selection().create("sel_air1", "Box");
    model.component("comp1").selection("sel_air1").set("entitydim", 2);
    model.component("comp1").selection("sel_air1").set("xmin", 0.45999999999999996);
    model.component("comp1").selection("sel_air1").set("xmax", 0.5);
    model.component("comp1").selection("sel_air1").set("ymin", -0.02);
    model.component("comp1").selection("sel_air1").set("ymax", 0.02);
    model.component("comp1").selection("sel_air1").set("condition", "intersects");
    // 空气楔 90°（贴上外边界）
    model.component("comp1").selection().create("sel_air2", "Box");
    model.component("comp1").selection("sel_air2").set("entitydim", 2);
    model.component("comp1").selection("sel_air2").set("xmin", -0.02);
    model.component("comp1").selection("sel_air2").set("xmax", 0.02);
    model.component("comp1").selection("sel_air2").set("ymin", 0.5342562584220407);
    model.component("comp1").selection("sel_air2").set("ymax", 0.5742562584220408);
    model.component("comp1").selection("sel_air2").set("condition", "intersects");
    // 空气楔 180°
    model.component("comp1").selection().create("sel_air3", "Box");
    model.component("comp1").selection("sel_air3").set("entitydim", 2);
    model.component("comp1").selection("sel_air3").set("xmin", -0.5);
    model.component("comp1").selection("sel_air3").set("xmax", -0.45999999999999996);
    model.component("comp1").selection("sel_air3").set("ymin", -0.02);
    model.component("comp1").selection("sel_air3").set("ymax", 0.02);
    model.component("comp1").selection("sel_air3").set("condition", "intersects");
    // 空气楔 270°（贴下外边界）
    model.component("comp1").selection().create("sel_air4", "Box");
    model.component("comp1").selection("sel_air4").set("entitydim", 2);
    model.component("comp1").selection("sel_air4").set("xmin", -0.02);
    model.component("comp1").selection("sel_air4").set("xmax", 0.02);
    model.component("comp1").selection("sel_air4").set("ymin", -0.5742562584220408);
    model.component("comp1").selection("sel_air4").set("ymax", -0.5342562584220407);
    model.component("comp1").selection("sel_air4").set("condition", "intersects");
    // 尖端 1：(+0.078774, +0.000000)
    model.component("comp1").selection().create("sel_tip1", "Box");
    model.component("comp1").selection("sel_tip1").set("entitydim", 0);
    model.component("comp1").selection("sel_tip1").set("xmin", 0.07542825964587464);
    model.component("comp1").selection("sel_tip1").set("xmax", 0.08211951211109893);
    model.component("comp1").selection("sel_tip1").set("ymin", -0.003345626232612136);
    model.component("comp1").selection("sel_tip1").set("ymax", 0.003345626232612136);
    model.component("comp1").selection("sel_tip1").set("condition", "inside");
    // 尖端 2：(-0.078774, +0.000000)
    model.component("comp1").selection().create("sel_tip2", "Box");
    model.component("comp1").selection("sel_tip2").set("entitydim", 0);
    model.component("comp1").selection("sel_tip2").set("xmin", -0.08211951211109893);
    model.component("comp1").selection("sel_tip2").set("xmax", -0.07542825964587464);
    model.component("comp1").selection("sel_tip2").set("ymin", -0.003345626232612136);
    model.component("comp1").selection("sel_tip2").set("ymax", 0.003345626232612136);
    model.component("comp1").selection("sel_tip2").set("condition", "inside");
    // 尖端 3：(+0.000000, +0.066913)
    model.component("comp1").selection().create("sel_tip3", "Box");
    model.component("comp1").selection("sel_tip3").set("entitydim", 0);
    model.component("comp1").selection("sel_tip3").set("xmin", -0.003345626232612136);
    model.component("comp1").selection("sel_tip3").set("xmax", 0.003345626232612136);
    model.component("comp1").selection("sel_tip3").set("ymin", 0.06356689841963058);
    model.component("comp1").selection("sel_tip3").set("ymax", 0.07025815088485486);
    model.component("comp1").selection("sel_tip3").set("condition", "inside");
    // 尖端 4：(+0.000000, -0.066913)
    model.component("comp1").selection().create("sel_tip4", "Box");
    model.component("comp1").selection("sel_tip4").set("entitydim", 0);
    model.component("comp1").selection("sel_tip4").set("xmin", -0.003345626232612136);
    model.component("comp1").selection("sel_tip4").set("xmax", 0.003345626232612136);
    model.component("comp1").selection("sel_tip4").set("ymin", -0.07025815088485486);
    model.component("comp1").selection("sel_tip4").set("ymax", -0.06356689841963058);
    model.component("comp1").selection("sel_tip4").set("condition", "inside");

    // ---------------- 材料（n = sqrt(eps)，与参考文件同一写法） ----------------
    model.component("comp1").material().create("mat_air1", "Common");
    model.component("comp1").material("mat_air1").selection().named("sel_air1");
    model.component("comp1").material("mat_air1").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat_air1").propertyGroup("RefractiveIndex").set("n", new String[]{"sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)"});
    model.component("comp1").material("mat_air1").propertyGroup("RefractiveIndex").set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material().create("mat_air2", "Common");
    model.component("comp1").material("mat_air2").selection().named("sel_air2");
    model.component("comp1").material("mat_air2").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat_air2").propertyGroup("RefractiveIndex").set("n", new String[]{"sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)"});
    model.component("comp1").material("mat_air2").propertyGroup("RefractiveIndex").set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material().create("mat_air3", "Common");
    model.component("comp1").material("mat_air3").selection().named("sel_air3");
    model.component("comp1").material("mat_air3").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat_air3").propertyGroup("RefractiveIndex").set("n", new String[]{"sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)"});
    model.component("comp1").material("mat_air3").propertyGroup("RefractiveIndex").set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material().create("mat_air4", "Common");
    model.component("comp1").material("mat_air4").selection().named("sel_air4");
    model.component("comp1").material("mat_air4").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat_air4").propertyGroup("RefractiveIndex").set("n", new String[]{"sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)", "0", "0", "0", "sqrt(epsilon_amb)"});
    model.component("comp1").material("mat_air4").propertyGroup("RefractiveIndex").set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material().create("mat_wg", "Common");
    model.component("comp1").material("mat_wg").selection().named("sel_wg");
    model.component("comp1").material("mat_wg").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat_wg").propertyGroup("RefractiveIndex").set("n", new String[]{"sqrt(epsilon_wg)", "0", "0", "0", "sqrt(epsilon_wg)", "0", "0", "0", "sqrt(epsilon_wg)"});
    model.component("comp1").material("mat_wg").propertyGroup("RefractiveIndex").set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});

    // ---------------- 物理场：ewfd + out-of-plane Ez ----------------
    model.component("comp1").physics().create("ewfd", "ElectromagneticWavesFrequencyDomain", "geom1");
    model.component("comp1").physics("ewfd").prop("components").set("components", "outofplane");
    model.component("comp1").physics("ewfd").create("port1", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port1").selection().named("sel_input");
    model.component("comp1").physics("ewfd").feature("port1").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port1").set("PortFreq", "f/T_unit");
    model.component("comp1").physics("ewfd").create("port2", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port2").selection().named("sel_output1");
    model.component("comp1").physics("ewfd").feature("port2").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port2").set("PortFreq", "f/T_unit");
    model.component("comp1").physics("ewfd").create("port3", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port3").selection().named("sel_output21");
    model.component("comp1").physics("ewfd").feature("port3").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port3").set("PortFreq", "f/T_unit");
    model.component("comp1").physics("ewfd").create("port4", "Port", 1);
    model.component("comp1").physics("ewfd").feature("port4").selection().named("sel_output22");
    model.component("comp1").physics("ewfd").feature("port4").set("PortType", "Numeric");
    model.component("comp1").physics("ewfd").feature("port4").set("PortFreq", "f/T_unit");
    // 两条非端口外边界：散射边界条件。绝不能留默认 —— ewfd 的默认外边界
    // 是 PEC（Ez = 0），那是完全不同的物理。
    model.component("comp1").physics("ewfd").create("sctr_top", "Scattering", 1);
    model.component("comp1").physics("ewfd").feature("sctr_top").selection().named("sel_scatter_top");
    model.component("comp1").physics("ewfd").create("sctr_bottom", "Scattering", 1);
    model.component("comp1").physics("ewfd").feature("sctr_bottom").selection().named("sel_scatter_bottom");

    // ---------------- 网格 ----------------
    // 全局 hmax = 空气里 λ/16，波导域再细一档，端口边按 h_port，4 个尖端单独压到 1e-3
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "0.04032258064516129");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "0.0004");
    // 波导域：hmax = 0.011859582542694497
    model.component("comp1").mesh("mesh1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("size1").selection().named("sel_wg");
    model.component("comp1").mesh("mesh1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size1").set("hmax", "0.011859582542694497");
    model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
    // 端口 input：hmax = 0.0037
    model.component("comp1").mesh("mesh1").create("size2", "Size");
    model.component("comp1").mesh("mesh1").feature("size2").selection().named("sel_input");
    model.component("comp1").mesh("mesh1").feature("size2").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size2").set("hmax", "0.0037");
    model.component("comp1").mesh("mesh1").feature("size2").set("hmaxactive", true);
    // 端口 output1：hmax = 0.0037
    model.component("comp1").mesh("mesh1").create("size3", "Size");
    model.component("comp1").mesh("mesh1").feature("size3").selection().named("sel_output1");
    model.component("comp1").mesh("mesh1").feature("size3").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size3").set("hmax", "0.0037");
    model.component("comp1").mesh("mesh1").feature("size3").set("hmaxactive", true);
    // 端口 output21：hmax = 0.0037
    model.component("comp1").mesh("mesh1").create("size4", "Size");
    model.component("comp1").mesh("mesh1").feature("size4").selection().named("sel_output21");
    model.component("comp1").mesh("mesh1").feature("size4").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size4").set("hmax", "0.0037");
    model.component("comp1").mesh("mesh1").feature("size4").set("hmaxactive", true);
    // 端口 output22：hmax = 0.0037
    model.component("comp1").mesh("mesh1").create("size5", "Size");
    model.component("comp1").mesh("mesh1").feature("size5").selection().named("sel_output22");
    model.component("comp1").mesh("mesh1").feature("size5").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size5").set("hmax", "0.0037");
    model.component("comp1").mesh("mesh1").feature("size5").set("hmaxactive", true);
    // 尖端 1：hmax = 0.001
    model.component("comp1").mesh("mesh1").create("size6", "Size");
    model.component("comp1").mesh("mesh1").feature("size6").selection().named("sel_tip1");
    model.component("comp1").mesh("mesh1").feature("size6").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size6").set("hmax", "0.001");
    model.component("comp1").mesh("mesh1").feature("size6").set("hmaxactive", true);
    // 尖端 2：hmax = 0.001
    model.component("comp1").mesh("mesh1").create("size7", "Size");
    model.component("comp1").mesh("mesh1").feature("size7").selection().named("sel_tip2");
    model.component("comp1").mesh("mesh1").feature("size7").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size7").set("hmax", "0.001");
    model.component("comp1").mesh("mesh1").feature("size7").set("hmaxactive", true);
    // 尖端 3：hmax = 0.001
    model.component("comp1").mesh("mesh1").create("size8", "Size");
    model.component("comp1").mesh("mesh1").feature("size8").selection().named("sel_tip3");
    model.component("comp1").mesh("mesh1").feature("size8").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size8").set("hmax", "0.001");
    model.component("comp1").mesh("mesh1").feature("size8").set("hmaxactive", true);
    // 尖端 4：hmax = 0.001
    model.component("comp1").mesh("mesh1").create("size9", "Size");
    model.component("comp1").mesh("mesh1").feature("size9").selection().named("sel_tip4");
    model.component("comp1").mesh("mesh1").feature("size9").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size9").set("hmax", "0.001");
    model.component("comp1").mesh("mesh1").feature("size9").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").run();

    // ---------------- 研究：4 步边值模分析 + 1 步频域 ----------------
    model.study().create("std1");
    model.study("std1").create("bma", "BoundaryModeAnalysis");
    model.study("std1").create("bma1", "BoundaryModeAnalysis");
    model.study("std1").create("bma2", "BoundaryModeAnalysis");
    model.study("std1").create("bma3", "BoundaryModeAnalysis");
    model.study("std1").create("freq", "Frequency");
    // 端口的数值模式：在 f 处求，shift 用解析板模的有效折射率当起点
    model.study("std1").feature("bma").set("modeFreq", "f/T_unit");
    model.study("std1").feature("bma").set("shift", "2.623");
    model.study("std1").feature("bma1").set("PortName", "2");
    model.study("std1").feature("bma1").set("modeFreq", "f/T_unit");
    model.study("std1").feature("bma1").set("shift", "2.623");
    model.study("std1").feature("bma2").set("PortName", "3");
    model.study("std1").feature("bma2").set("modeFreq", "f/T_unit");
    model.study("std1").feature("bma2").set("shift", "2.623");
    model.study("std1").feature("bma3").set("PortName", "4");
    model.study("std1").feature("bma3").set("modeFreq", "f/T_unit");
    model.study("std1").feature("bma3").set("shift", "2.623");
    model.study("std1").feature("freq").set("plist", "f/T_unit");
    model.study("std1").feature("freq").set("ngen", 1);
    model.study("std1").run();

    // ---------------- 结果 ----------------
    // gev1：R / T / A 功率守恒；gev2：S11..S41（复数值）
    model.result().table().create("tbl1", "Table");
    model.result().table().create("tbl2", "Table");
    model.result().numerical().create("gev1", "EvalGlobal");
    model.result().numerical().create("gev2", "EvalGlobal");
    model.result().numerical("gev1").set("probetag", "none");
    model.result().numerical("gev1").set("table", "tbl1");
    model.result().numerical("gev1").set("expr", new String[]{"ewfd.Rport_1", "ewfd.Rtotal", "ewfd.Tport_2", "ewfd.Tport_3", "ewfd.Tport_4", "ewfd.Ttotal", "ewfd.RTtotal", "ewfd.Atotal"});
    model.result().numerical("gev1").set("descr", new String[]{"R (port 1)", "R total", "T (port 2)", "T (port 3)", "T (port 4)", "T total", "R+T", "A (= 从散射边界漏掉的)"});
    model.result().numerical("gev2").set("probetag", "none");
    model.result().numerical("gev2").set("table", "tbl2");
    model.result().numerical("gev2").set("expr", new String[]{"ewfd.S11", "ewfd.S21", "ewfd.S31", "ewfd.S41"});
    model.result().numerical("gev2").set("descr", new String[]{"S11", "S21", "S31", "S41"});
    model.result().numerical("gev1").setResult();
    model.result().numerical("gev2").setResult();
    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").create("surf1", "Surface");
    model.result("pg1").feature("surf1").set("expr", "ewfd.normE");
    model.result("pg1").feature("surf1").set("colortable", "RainbowLight");

    return model;
  }

  public static void main(String[] args) {
    run();
  }
}
